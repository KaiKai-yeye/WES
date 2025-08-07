"""
pileup 汇总，用于 后续肿瘤纯度评估、变异过滤（尤其是 FilterMutectCalls）中计算群体等位基因频率（allele frequency）等目的。
将每个.aligned.duplicates_marked.recalibrated.bam整合成对应的一个.pileups.table文件
################################################输入文件################################################################################
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align  # 输入：BQSR校正后的BAM文件目录

output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/pileup  # 输出：pileup结果表格目录

# 参考基因组FASTA文件（GRCh38版本）
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa

# 目标区间BED文件（Agilent捕获芯片的外显子区域，限定分析范围）
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed

# 公共种系变异VCF（包含等位基因频率信息，用于计算pileup）
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz

##################################################################################################################################################
gatk GetPileupSummaries \
    -R "$reference_path" \  # -R: 指定参考基因组FASTA文件路径（必须与BAM文件的参考版本一致）
    -I "$file" \          # -I: 输入的BAM文件（经过BQSR校正的样本数据）
    -V "$germline_resource_path" \  # -V: 种系变异VCF文件（包含等位基因频率，用于计算pileup）
    -L "$interval_path" \  # -L: 目标区间BED文件（限定仅分析外显子等捕获区域，提高效率）
    -O "$output_table"; then  # -O: 输出的pileup结果表格路径（包含每个位点的等位基因计数等信息）
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
chr1	 123456	       A	          G	           0.015	           200	        15	         2
这个例子表示在 chr1 的第 123456 位点：
参考等位基因为 A，变异等位基因为 G；G 这个突变在人群中出现频率是 1.5%；在当前样本中有 200 个 reads 支持 A（ref），15 个支持 G（alt），还有 2 个支持其他碱基。
"""
#!/bin/bash
#SBATCH -J Pileup_0807                           # 作业名称
#SBATCH -N 1                                   # 所需节点数
#SBATCH -n 64                                 # 请求CPU核心数
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/7_Pileup_01.o  # 标准输出
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/7_Pileup_01.e  # 错误输出
############################ 环境激活 ############################
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/
############################ 控制并发进程数量 ############################
Nproc=16
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=Nproc; i++)); do
    echo
done >&6
############################ 设置路径 ############################
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/pileup
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
mkdir -p "$output_folder_path"
############################ 记录开始时间 ############################
echo "Start Time: $(date +"%Y-%m-%d %H:%M:%S")"
############################ 遍历BAM文件并执行Pileup ############################
for file in "$input_folder_path"/*.aligned.duplicates_marked.recalibrated.bam; do
    filename=$(basename "$file")
    sample_name="${filename%.aligned.duplicates_marked.recalibrated.bam}"
    output_table="${output_folder_path}/${sample_name}.pileups.table"
    echo "Checking sample: $sample_name at $(date +"%Y-%m-%d %H:%M:%S")"
    if [ -e "$output_table" ]; then
        echo "[SKIP] Output exists: $output_table"
        continue
    fi
    read -u6
    {
        echo "[RUN] Running GetPileupSummaries for sample: $sample_name"
        # 核心命令：添加错误捕获，仅针对gatk命令报错
        if ! gatk GetPileupSummaries \
            -R "$reference_path" \
            -I "$file" \
            -V "$germline_resource_path" \
            -L "$interval_path" \
            -O "$output_table"; then
            # 命令执行失败时输出详细错误信息
            echo "[ERROR] GetPileupSummaries failed for sample: $sample_name" >&2
            echo "[ERROR] 失败命令: gatk GetPileupSummaries -R $reference_path -I $file -V $germline_resource_path -L $interval_path -O $output_table" >&2
            # 若需要中断整个脚本，可取消下面一行的注释；否则仅标记该样本失败，继续执行其他样本
            # exit 1
        fi
        echo "[DONE] Completed: $sample_name (若上方有ERROR提示，则该样本处理失败)"
        sleep 2
        echo >&6
    } &
done
wait
exec 6>&-
echo "End Time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "07_Pileup done! (注意检查是否有样本处理失败的ERROR提示)"
