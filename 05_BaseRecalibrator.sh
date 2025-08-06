"""
纠正测序仪报告的碱基质量分数中的系统性偏差,提升变异检测的准确性

输出:
文件	                                         描述
.recal_data.csv	                                 校准数据表，记录了每种错误类型的系统偏差
.aligned.duplicates_marked.recalibrated.bam      应用校准后的 BAM 文件，用于后续变异检测
.aligned.duplicates_marked.recalibrated.bai      索引
.aligned.duplicates_marked.recalibrated.bam.md5  是一个 MD5 校验文件，它的作用是记录对应 BAM 文件的校验和，用于确保该 BAM 文件在传输或存储过程中 没有被损坏或篡改。文件很小


"""
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
reference_path="/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa"  # 参考基因组
interval_path="/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed"  # 目标区域BED

# 已知变异位点文件（BQSR核心输入）
known_sites_path_1="/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf"
known_sites_path_2="/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
known_sites_path_3="/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_sapiens_assembly38.known_indels.vcf.gz"

echo "Start Time: $(date +"%Y-%m-%d %H:%M:%S")"

###################### Step 4: 批量执行BQSR ######################
for input_bam in "${input_folder_path}"/*.aligned.duplicate_marked.sorted.bam; do
    # 提取样本名（去除文件名后缀）
    filename=$(basename "$input_bam")
    sample_name="${filename%.aligned.duplicate_marked.sorted.bam}"
    
    # 输出文件路径
    recal_csv="${output_folder_path}/${sample_name}.recal_data.csv"  # 校准模型
    output_bam="${output_folder_path}/${sample_name}.aligned.duplicates_marked.recalibrated.bam"  # 校准后BAM

    echo "Processing sample: $sample_name (start at: $(date +"%H:%M:%S"))"

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
            echo "[DONE] $sample_name: BQSR completed"
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
