"""
纠正测序仪报告的碱基质量分数中的系统性偏差,提升变异检测的准确性

输出:
文件	                                           描述
.recal_data.csv	                                 校准数据表，记录了每种错误类型的系统偏差
.aligned.duplicates_marked.recalibrated.bam      应用校准后的 BAM 文件，用于后续变异检测
"""
#!/bin/bash
#SBATCH -J BQSR_0729                        # SLURM作业名称
#SBATCH -N 1                                # 申请1个节点
#SBATCH -n 16                               # 申请16个核心
#SBATCH -o BQSR.o                           # 标准输出日志文件
#SBATCH -e BQSR.e                           # 错误输出日志文件

###################### Step 1: 加载环境 ######################
# 激活 Conda 环境，其中包含 GATK 相关工具
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

###################### Step 2: 初始化并发控制 ######################
# 控制最大并行任务数（同时运行的样本数）
Nproc=2

# 创建唯一命名管道（FIFO）用于并发控制
Pfifo="/tmp/$$.fifo"         # "$$"表示当前脚本进程的PID，避免命名冲突
mkfifo $Pfifo
exec 6<>$Pfifo               # 将管道与文件描述符6关联（读写）
rm -f $Pfifo                 # 删除原始管道文件名，不影响FD使用

# 初始化并发令牌（每个空行代表一个任务槽位）
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

###################### Step 3: 路径设定 ######################
# 输入目录：包含已经排序并标记重复的 BAM 文件
input_folder_path=/groups/g5840141/home/zengqianwen/WES/align

# 输出目录：与输入一致，BQSR 输出 BAM 也存到这里
output_folder_path=/groups/g5840141/home/zengqianwen/WES/align

# 参考基因组 fasta（.fa）
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa

# 目标区域 BED 文件，限制校准只发生在指定捕获区域
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed

# GATK 推荐的已知变异位点文件（用于识别系统性测序误差）
known_sites_path_1=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf
known_sites_path_2=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_sites_path_3=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

echo "Start Time: $(date +"%Y-%m-%d %H:%M:%S")"

###################### Step 4: 遍历 BAM 文件并执行 BQSR ######################

# 查找所有符合条件的 BAM 文件（已排序并标记重复）
for input_bam in ${input_folder_path}/*.aligned.duplicate_marked.sorted.bam; do

    # 从文件名中提取样本名（例如：PHCC2417T1）
    filename=$(basename "$input_bam")
    sample_name="${filename%.aligned.duplicate_marked.sorted.bam}"#使用的是标记重复后的bam文件

    # 设定输出文件路径
    recal_csv=${output_folder_path}/${sample_name}.recal_data.csv  # 校准模型文件
    output_bam=${output_folder_path}/${sample_name}.aligned.duplicates_marked.recalibrated.bam  # 最终输出 BAM

    echo "Checking sample: $sample_name at $(date +"%Y-%m-%d %H:%M:%S")"

    # 如果已存在输出文件，则跳过该样本
    if [ -e "$output_bam" ]; then
        echo "[SKIP] $output_bam already exists."
    else
        echo "[RUN] Starting BQSR for $sample_name"

        read -u6  # 从FD6中读取一个令牌（代表获取一个可用并发槽位）
        {
            ################## Step 4.1: 构建重校准模型 ##################
            # 基于已知变异，GATK 将识别测序误差并建立校准模型
            gatk BaseRecalibrator \
                -R $reference_path \                      # 参考基因组
                -I $input_bam \                           # 输入 BAM
                -L $interval_path \                       # 限定区域（BED）
                --use-original-qualities \                # 使用原始质量值进行建模
                --known-sites $known_sites_path_1 \       # 常见 SNP 变异
                --known-sites $known_sites_path_2 \       # 高可信 InDel 变异
                --known-sites $known_sites_path_3 \       # 额外 InDel 数据
                -O $recal_csv                             # 输出：重校准模型

            ################## Step 4.2: 应用重校准模型 ##################
            # 使用上一步生成的模型对原 BAM 文件中的碱基质量值进行调整
            gatk ApplyBQSR \
                -R $reference_path \                      # 参考基因组
                -I $input_bam \                           # 原始 BAM
                -bqsr $recal_csv \                        # 重校准模型
                -O $output_bam \                          # 输出新的 BAM
                --add-output-sam-program-record \         # 记录处理软件信息
                --create-output-bam-md5 \                 # 创建 MD5 校验文件
                --use-original-qualities                  # 保留原始质量值

            echo "[DONE] BQSR completed for $sample_name"
            echo >&6  # 归还令牌
        } &
    fi
done

###################### Step 5: 等待所有任务完成 ######################
wait            # 等待所有后台任务执行完成
exec 6>&-       # 关闭并释放文件描述符6

echo "End Time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "BQSR step done!"
