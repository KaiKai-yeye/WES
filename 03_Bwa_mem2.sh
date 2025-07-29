#!/bin/bash
#SBATCH -J Bwa_mem2_0729
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o 04_Bwa_mem2.o
#SBATCH -e 04_Bwa_mem2.e

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

####################################
#### 并行控制：限制并发作业数 ####
####################################
Nproc=4                               # 最大并发数
Pfifo="/tmp/$$.fifo"                  # 命名管道
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

####################################
#### 路径配置 ####
####################################
input_folder_path=/groups/g5840141/home/zengqianwen/WES/trim
output_folder_path=/groups/g5840141/home/zengqianwen/WES/align
bwa_fa_reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/bwa-mem2/GRCh38.primary_assembly.genome.fa

# 开始标记
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "04_Bwa_mem2 start!"

####################################
#### 自动提取样本名并比对 ####
####################################
for file_1 in ${input_folder_path}/*_val_1.fq.gz; do
    # 获取样本名
    sample_name=$(basename ${file_1} _val_1.fq.gz)
    file_2=${input_folder_path}/${sample_name}_val_2.fq.gz
    bam_file=${output_folder_path}/${sample_name}.aligned.sorted.bam

    echo "[INFO] Processing sample: $sample_name"

    # 检查配对文件是否存在
    if [ ! -f "$file_2" ]; then
        echo "[WARNING] Paired file missing for $sample_name: $file_2 not found. Skipping."
        continue
    fi

    # 如果 BAM 文件已存在且非空，跳过
    if [ -s "$bam_file" ]; then
        echo "[INFO] $bam_file exists and is non-empty. Skipping."
        continue
    else
        echo "[INFO] $bam_file not found or is empty. Proceeding with alignment."
    fi

    # 提取平台信息用于 Read Group
    # -------------------------------------------
# 提取测序头部信息用于构建 Read Group 标签
#
# Illumina FASTQ 文件的第一行格式一般如下：
# @<instrument_id>:<run_id>:<flowcell_id>:<lane_number>:...
#
# 各字段解释：
# - instrument_id：测序仪编号（如 A00123）
# - run_id：第几次运行（如 98）
# - flowcell_id：芯片编号（如 H3WLVDSXX）
# - lane_number：芯片上的通道编号（如 2）
#
# 构建变量说明：
# - platform_unit：表示该组 reads 来源的测序单元（flowcell.run.lane）
# - read_group_id：给每个样本生成唯一的 ID（sample.flowcell.lane），供 GATK 等工具识别
#
# 这些信息会添加到比对命令中的 -R @RG 标签，用于下游分析如去重复、样本区分等
# -------------------------------------------

    lane_num=$(zless "$file_1" | head -n1 | cut -d ":" -f4)
    flow_id=$(zless "$file_1" | head -n1 | cut -d ":" -f3)
    run_id=$(zless "$file_1" | head -n1 | cut -d ":" -f2)
    instrument_id=$(zless "$file_1" | head -n1 | cut -d ":" -f1 | sed s/@//)

    platform_unit="${flow_id}.${run_id}.${lane_num}"
    read_group_id="${sample_name}.${flow_id}.${lane_num}"

    # 控制并发：领取令牌
    read -u6

    {
        echo "[INFO] Starting alignment for $sample_name"
        # 执行 BWA-MEM2 比对并排序输出 BAM 文件
        bwa-mem2 mem \
            -t 16 \                               # 使用 16 线程并行加速
            -K 100000000 \                        # 读取块大小为 100MB，优化性能
            -Y \                                  # 输出软剪接信息（为 GATK 兼容）
            -R "@RG\\tID:${read_group_id}\\tLB:WES\\tSM:${sample_name}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tPM:${instrument_id}" \  # 添加 Read Group 信息
            "$bwa_fa_reference_path" \           # 输入参考基因组
            "$file_1" "$file_2" | \                # 输入配对的 fastq 文件
        samtools sort \
            -@ 16 \                               # samtools 使用 16 个线程排序
            -m 1500M \                            # 每线程最大内存 1.5G
            --write-index \                       # 同时生成 BAM 索引（.bai 文件）
            -o "${bam_file}"                      # 输出 BAM 文件路径

        echo "[INFO] Alignment done for $sample_name"

        sleep 5
        echo >&6  # 归还令牌
    } &
done

wait             # 等待所有子进程
exec 6>&-        # 关闭管道

# 结束标记
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "Bwa_mem2 done!"
