"""
pileup 汇总，用于 后续肿瘤纯度评估、变异过滤（尤其是 FilterMutectCalls）中计算群体等位基因频率（allele frequency）等目的。
将每个.aligned.duplicates_marked.recalibrated.bam整合成对应的一个.pileups.table文件
##################################################################################################################################################
# GATK GetPileupSummaries 工作原理说明：
#
# 此步骤的目标是从肿瘤样本的 BAM 文件中，提取在一组已知人群常见变异位点处的等位基因支持情况，
# 为后续的污染度评估（CalculateContamination）提供输入数据。
#
# 处理逻辑如下：
# 1. 读取肿瘤样本的 BAM 文件（-I 参数），从中提取测序 reads。
# 2. 读取提供的 VCF 文件（-V 参数），该文件包含常见人群变异位点及其等位基因频率（如 ExAC 数据）。
# 3. 对 VCF 中每一个位点，在 BAM 文件中的BED文件限定区域进行 pileup 操作（即VCF和BED的交集）：
#    - 统计多少 reads 支持参考等位基因（REF）
#    - 多少 reads 支持变异等位基因（ALT）
#    - 多少 reads 支持其他非 REF/ALT 的碱基（other alt）
# 4.参考基因组序列（FASTA格式），用于确认每个位点的参考碱基，保证 BAM、VCF、BED 在坐标系统和染色体名称上的一致性；
# 5. 输出一个 TSV 格式的 pileup 汇总表（pileups.table），用于后续计算污染度：
#    - 包含位置、REF/ALT碱基、人群频率、支持数量等信息
#
# 输出的 pileups.table 是 CalculateContamination 的直接输入，
# 可以用于判断样本是否存在外源 DNA 混入（交叉污染）。
##################################################################################################################################################

🧬 pileup 文件包含的关键信息
| 列名               | 含义                              
| ----------------- | ---------------------------------------------------- 
| `contig`          | 染色体名（如 chr1、chrX）                          
| `position`        | 基因组上的位置（1-based）                           
| `refAllele`       | 参考基因组上的碱基                                  
| `altAllele`       | 目标变异（常见变异数据库中提供的 ALT）               
| `alleleFrequency` | ALT 等位基因在人群中的频率                  
| `refCount`        | 在该位点上观测到 REF 碱基的 reads 数量       
| `altCount`        | 在该位点上观测到 ALT 碱基的 reads 数量       
| `otherCount`      | 观测到除 REF/ALT 以外的其他碱基数量（可能是测序错误） 

📌 举个例子：
contig	 position	refAllele	altAllele	alleleFrequency	  refCount	   altCount	    otherCount
chr1	 123456	    A	        G	        0.015	           200	        15	         2
这个例子表示在 chr1 的第 123456 位点：
参考等位基因为 A，变异等位基因为 G；G 这个突变在人群中出现频率是 1.5%；在当前样本中有 200 个 reads 支持 A（ref），15 个支持 G（alt），还有 2 个支持其他碱基。
"""
#!/bin/bash
#SBATCH -J 07_Pileup                           # 作业名称
#SBATCH -N 1                                   # 所需节点数（通常为1个节点即可）
#SBATCH -n 16                                  # 请求使用的CPU核心数
#SBATCH -o 07_Pileup.o                         # 标准输出文件
#SBATCH -e 07_Pileup.e                         # 标准错误输出文件

############################ 环境激活 ############################

# 激活 Conda 环境，加载包含 GATK 的环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

############################ 控制并发进程数量 ############################

Nproc=4                                        # 同时允许最多4个任务并发执行（根据资源配置可调整）
Pfifo="/tmp/$$.fifo"                           # 创建临时命名管道文件，$$ 表示当前脚本进程ID
mkfifo $Pfifo                                  # 创建命名管道文件
exec 6<>$Pfifo                                 # 打开管道文件并绑定到文件描述符6
rm -f $Pfifo                                   # 删除文件名（文件描述符仍然可用）

# 在文件描述符6中写入Nproc个“令牌”，每个令牌代表一个可用的并发“槽位”
for ((i=1; i<=Nproc; i++)); do
    echo
done >&6

############################ 设置输入输出路径和参考文件 ############################

input_folder_path=/groups/g5840141/home/zengqianwen/WES/align                  # 输入：BQSR后的BAM文件目录
output_folder_path=/groups/g5840141/home/zengqianwen/WES/pileup               # 输出：pileup结果表格目录
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa

# 用于限制分析范围的目标区域BED文件（对于wes,主要用于限定外显子区域）
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed

# 公共种系变异VCF，用于获取 pileup（变异等位基因频率估计）
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz

# 确保输出文件夹存在
mkdir -p "$output_folder_path"

############################ 记录开始时间 ############################

echo "Start Time: $(date +"%Y-%m-%d %H:%M:%S")"

############################ 遍历BAM文件并执行Pileup ############################

# 遍历目录中所有BQSR后的BAM文件
for file in "$input_folder_path"/*.aligned.duplicates_marked.recalibrated.bam; do
    # 提取样本名（去除路径和后缀）
    filename=$(basename "$file")
    sample_name="${filename%.aligned.duplicates_marked.recalibrated.bam}"

    # 定义输出文件路径
    output_table="${output_folder_path}/${sample_name}.pileups.table"

    echo "Checking sample: $sample_name at $(date +"%Y-%m-%d %H:%M:%S")"

    # 如果该样本的pileup结果已存在，跳过处理
    if [ -e "$output_table" ]; then
        echo "[SKIP] Output exists: $output_table"
        continue
    fi

    # 从令牌池中获取一个令牌，控制并发任务数量
    read -u6
    {
        echo "[RUN] Running GetPileupSummaries for sample: $sample_name"

        # 执行GATK GetPileupSummaries
        gatk GetPileupSummaries \
            -R "$reference_path" \                           # 指定参考基因组FASTA
            -I "$file" \                                     # 输入BAM文件（BQSR后）
            -V "$germline_resource_path" \                   # 公共种系变异资源（VCF格式）
            -L "$germline_resource_path" \                   # 用于从VCF中提取目标区域的注释
            -L "$interval_path" \                            # 指定捕获区域BED文件
            --interval-set-rule INTERSECTION \               # 仅分析两个L参数交集部分
            -O "$output_table"                               # 输出pileup结果表格

        echo "[DONE] Completed: $sample_name"
        sleep 2

        echo >&6  # 归还令牌，释放槽位
    } &
done

############################ 等待所有任务完成 ############################

wait
exec 6>&-  # 关闭文件描述符6

############################ 输出完成信息 ############################

echo "End Time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "07_Pileup done!"
