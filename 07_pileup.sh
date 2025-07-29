"""
pileup 汇总，用于 后续肿瘤纯度评估、变异过滤（尤其是 FilterMutectCalls）中计算群体等位基因频率（allele frequency）等目的。
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

# 用于限制分析范围的目标区域BED文件
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
