#!/bin/bash
#SBATCH -J BQSR_0805                        # SLURM作业名称
#SBATCH -N 1                                # 申请1个节点
#SBATCH -n 64                               # 申请64个核心（根据实际需求调整）
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/1_BQSR_0805.o  # 标准输出日志
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/1_BQSR_0805.e  # 错误输出日志

###################### Step 1: 加载环境 ######################
# 激活包含GATK的Conda环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

###################### Step 2: 初始化并发控制 ######################
Nproc=8  # 最大并行任务数（根据节点核心数调整）
Pfifo="/tmp/$$.fifo"  # 唯一命名管道，避免冲突
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo  # 删除文件实体，保留文件描述符

# 初始化并发令牌
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

###################### Step 3: 路径设定 ######################
input_folder_path="/groups/g5840141/home/zengqianwen/WES_2025/align"  # 输入BAM目录
output_folder_path="/groups/g5840141/home/zengqianwen/WES_2025/align"  # 输出目录
reference_path="/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38
.primary_assembly.genome.fa"  # 参考基因组
interval_path="/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padd
ed.bed"  # 目标区域BED

# 已知变异位点文件（BQSR核心输入）
known_sites_path_1="/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_s
apiens_assembly38.dbsnp138.vcf"
known_sites_path_2="/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Mills_
and_1000G_gold_standard.indels.hg38.vcf.gz"
known_sites_path_3="/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_s
apiens_assembly38.known_indels.vcf.gz"

echo "Start Time: $(date +"%Y-%m-%d %H:%M:%S")"
###################### Step 4: 批量执行BQSR ######################

# === 新增：计算总文件数并初始化计数器 ===
total_files=$(ls -1 "${input_folder_path}"/*.aligned.duplicate_marked.sorted.bam | wc -l)
processed_count=0
echo "Total samples to process: $total_files"
# === 新增结束 ===

for input_bam in "${input_folder_path}"/*.aligned.duplicate_marked.sorted.bam; do
    # 提取样本名（去除文件名后缀）
    filename=$(basename "$input_bam")
    sample_name="${filename%.aligned.duplicate_marked.sorted.bam}"

    # 输出文件路径
    recal_csv="${output_folder_path}/${sample_name}.recal_data.csv"  # 校准模型
    output_bam="${output_folder_path}/${sample_name}.aligned.duplicates_marked.recalibrated.bam"  # 校
准后BAM

    # === 修改：在循环开始时更新计数器 ===
    ((processed_count++))
    # === 修改结束 ===

    echo "Processing sample: $sample_name ($processed_count/$total_files) (start at: $(date +"%H:%M:%S"
))"

    # 跳过已完成的样本
    if [ -e "$output_bam" ] && [ -e "${output_bam}.bai" ]; then
        echo "[SKIP] $sample_name: Output files already exist."
        continue
    fi

    # 并发控制：获取一个令牌
    read -u6
    {
        ################## Step 4.1: 构建重校准模型 ##################
        echo "[RUN] $sample_name: Starting BaseRecalibrator"
        gatk BaseRecalibrator \
            -R "${reference_path}" \
            -I "${input_bam}" \
            -L "${interval_path}" \
            --use-original-qualities \
            --known-sites "${known_sites_path_1}" \
            --known-sites "${known_sites_path_2}" \
            --known-sites "${known_sites_path_3}" \
            -O "${recal_csv}"

        # 检查模型文件是否生成
        if [ ! -e "${recal_csv}" ]; then
            echo "[ERROR] $sample_name: BaseRecalibrator failed (no output CSV)"
            echo >&6  # 释放令牌
            continue
        fi

        ################## Step 4.2: 应用重校准模型 ##################
        echo "[RUN] $sample_name: Starting ApplyBQSR"
        gatk ApplyBQSR \
            -R "${reference_path}" \
            -I "${input_bam}" \
            -bqsr "${recal_csv}" \
            -O "${output_bam}" \
            --create-output-bam-md5 true \
            --use-original-qualities

        # 检查输出BAM是否生成
        if [ -e "$output_bam" ]; then
            # === 修改：计算并显示进度 ===
            progress=$(printf "%.2f%%" "$(echo "scale=4; $processed_count / $total_files * 100" | bc)")
            echo "[DONE] $sample_name: BQSR completed. Progress: $processed_count/$total_files ($progre
ss)"
            # === 修改结束 ===
        else
            echo "[ERROR] $sample_name: ApplyBQSR failed (no output BAM)"
        fi

        # 释放令牌
        echo >&6
    } &
done

###################### Step 5: 等待所有任务完成 ######################
wait
exec 6>&-  # 关闭文件描述符

echo "End Time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "BQSR pipeline finished!"

