----------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------fastqc--------------------------------------------------------------
#!/bin/bash
#SBATCH -J 02_FastQC                        # 任务名称
#SBATCH -N 1                                # 请求1个节点
#SBATCH -n 16                               # 请求16个CPU核心
#SBATCH -o 02_FastQC.o                      # 标准输出文件
#SBATCH -e 02_FastQC.e                      # 标准错误文件

# 加载并激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#######################
### 控制并发任务数量 ###
#######################
Nproc=12                                    # 最大并发数，根据核心数和每个任务线程数设置
Pfifo="/tmp/$$.fifo"                        # 使用当前进程号作为唯一的管道文件名
mkfifo $Pfifo                               # 创建命名管道
exec 6<>$Pfifo                              # 打开管道文件作为fd 6
rm -f $Pfifo                                # 删除原始管道文件（文件描述符仍然保留）

# 初始化并发令牌池：放入 $Nproc 个空行，相当于 N 个“令牌”
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

##############################
### 路径配置（请替换成所需路径） ###
##############################
# 输入文件夹路径，请替换成你的数据所在目录
input_folder_path=/path/to/your/input_folder

# FastQC 输出结果路径，请替换成你想保存结果的目录
output_folder_path=/path/to/your/fastqc_output_folder

# MultiQC 输出结果路径，请替换成你想保存整合结果的目录
multiqc_output_path=/path/to/your/multiqc_output_folder

# 打印任务开始时间
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "02_FastQC start!"

###########################
### 遍历所有样本并处理 ###
###########################
for file_1 in ${input_folder_path}/*_1.fq.gz; do
    sample_name=$(basename ${file_1} _1.fq.gz)
    file_2=$(echo ${file_1} | sed 's/_1.fq.gz/_2.fq.gz/')

    # 检查是否已经有对应 FastQC 报告文件，避免重复分析
    if ls ${output_folder_path}/${sample_name}*fastqc.html 1> /dev/null 2>&1; then
        echo "${sample_name} already run FastQC"
    else
        echo "${sample_name} does not run FastQC"

        read -u6  # 领取一个令牌（无令牌则阻塞等待）
        {
            # 执行 FastQC（为每个样本分配4线程）
            fastqc -t 4 \
                ${input_folder_path}/${sample_name}*.fq.gz \
                -o ${output_folder_path}

            sleep 5  # 可选：缓解IO压力

            echo >&6  # 归还令牌
        } &
    fi
done

# 等待所有后台任务完成
wait
exec 6>&-  # 关闭文件描述符

# 打印任务结束时间
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"

#####################
### 运行 MultiQC ###
#####################
echo "Running MultiQC..."

----------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------trim--------------------------------------------------------------
#!/bin/bash
#SBATCH -J 03_TrimGalore
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/TrimGalore_0727.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/TrimGalore_0727.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=8    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 6<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作
# 在fd6中放置$Nproc个空行作为令牌
for((i=1; i<=$Nproc; i++)); do
        echo
done >&6


# 使用Trim Galore去除低质量的reads和adaptor
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/RawData
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/trim

mkdir -p "$output_folder_path"

for file_1 in ${input_folder_path}/*_1.fq.gz; do
    sample_name=$(basename ${file_1} _1.fq.gz)          # 提取样本名，去掉后缀
    file_2=$(echo ${file_1} | sed s/_1.fq.gz/_2.fq.gz/) # 替换R1为R2的文件名

        current_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Current Time: $current_time"

        if [ -e ${output_folder_path}/${sample_name}*.fq.gz ]; then
                echo "${sample_name} already run TrimGalore"
        else
                echo "${sample_name} dose not run TrimGalore"
                
                read -u6 # 领取令牌, 控制进程数量
                {
                        trim_galore \
                                -j 4 \
                                --paired ${file_1} ${file_2} \
                                --output_dir ${output_folder_path}
        
                        sleep 5
                        echo >&6 # 归还令牌
                }&
        fi
done
wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
~

multiqc ${output_folder_path} -o ${multiqc_output_path}
echo "MultiQC done!"

----------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------Bwa_mem2（比对）-----------------------------------------------------
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

----------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------Markdup（标记重复）-----------------------------------------------------
#!/bin/bash
#SBATCH -J Markdup_0729
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o Markdup_0729.o  # STDOUT 输出日志
#SBATCH -e Markdup_0729.e  # STDERR 错误日志

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

# ----------------------------
# 并行控制：最大并发任务数
# ----------------------------
Nproc=4
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

# ----------------------------
# 路径设置
# ----------------------------
input_folder_path=/groups/g5840141/home/zengqianwen/WES/align
output_folder_path=/groups/g5840141/home/zengqianwen/WES/align
temp_dir_path=/groups/g5840141/home/zengqianwen/WES/tmp

# ----------------------------
# 当前时间记录
# ----------------------------
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "[START] Current Time: $current_time"

# ----------------------------
# 遍历所有 aligned.sorted.bam 文件
# ----------------------------
for file in ${input_folder_path}/*.aligned.sorted.bam; do
    sample_name=$(basename "$file" .aligned.sorted.bam)
    output_bam="${output_folder_path}/${sample_name}.aligned.duplicate_marked.sorted.bam"

    current_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "[INFO] Checking $sample_name - $current_time"

    if [ -s "$output_bam" ]; then
        echo "[SKIP] $output_bam already exists and is not empty."
        continue
    fi

    # 获取并发令牌
    read -u6
    {
        echo "[RUNNING] Marking duplicates for $sample_name"

        sambamba markdup \
            -t 16 \
            --tmpdir="${temp_dir_path}" \    # 指定临时文件存放目录，避免临时文件占用默认系统目录空间
            --hash-table-size=262144 \       # 设置内部哈希表大小，用于加速重复片段检测，单位是哈希桶数
            --overflow-list-size=67108864 \  # 设置溢出列表大小，单位是字节，用于处理哈希冲突时存储数据，增大可减少内存溢出错误
            "$file" \
            "$output_bam"

        echo "[DONE] $sample_name duplicate marking completed."
        sleep 5
        echo >&6  # 归还令牌
    } &
done

wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "[END] Current Time: $current_time"
echo "05_Markdup done!"

----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------Base Quality Score Recalibration（基因质量重校准）手动输入样本------------------------------------------
"""
根据参考和已知变异模型调整为每个碱基的质量值
假设某碱基原本被测序仪赋予质量值 Q=30（错误率 0.1%），但在过去经验中发现：
某特定测序仪，在某一循环、序列上下文中，会系统性地高估质量分数。
那么 BQSR 就会把它调整为，比如 Q=25（错误率 0.3%），更符合实际。
BQSR 让 BAM 文件里的碱基质量值更可靠，提升变异检测的准确率，但不会改变序列或对齐本身。
"""
#!/bin/bash
#SBATCH -J BQSR_0729                      # 作业名称
#SBATCH -N 1                              # 所需节点数（通常为1）
#SBATCH -n 16                             # 请求使用的核心数
#SBATCH -o BQSR.o                      # 标准输出文件
#SBATCH -e BQSR.e                      # 错误输出文件

# 激活环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

###################### 并行进程控制部分 ######################

Nproc=2                                   # 控制最多同时运行 2 个样本任务
Pfifo="/tmp/$$.fifo"                      # 使用当前进程号创建唯一命名管道
mkfifo $Pfifo
exec 6<>$Pfifo                            # 打开命名管道，绑定到文件描述符6
rm -f $Pfifo                              # 删除命名管道文件（fd依然可用）

# 写入令牌以控制并发
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

###################### 路径与资源设置 ######################

input_folder_path=/groups/g5840141/home/zengqianwen/WES/align
output_folder_path=/groups/g5840141/home/zengqianwen/WES/align

# 参考基因组和目标区域
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed

# 已知变异位点资源，用于 BQSR
known_sites_path_1=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_sapiens_assembly38.dbsnp138.vcf
known_sites_path_2=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_sites_path_3=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-resource-bundle/v0/Homo_sapiens_assembly38.known_indels.vcf.gz

# 样本列表（已手动指定）
Sample_list=(PHCC2417T1 PHCC2417T2 PHCC2417T3 PHCC2417T4 PHCC2417P RHCC4584T PHCC1011T PHCC1011L RHCC4619T RHCC4619P PHCC966T RHCC4664T RHCC4664P PHCC3603T RHCC4691T RHCC4691P PHCC972T PHCC972L RHCC4349T PHCC3369T PHCC3369P RHCC4121T RHCC4121L)

# 打印当前时间
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"

###################### 主循环处理每个样本 ######################

for sample_name in ${Sample_list[@]}; do
    file=${input_folder_path}/${sample_name}.aligned.duplicate_marked.sorted.bam

    current_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Current Time: $current_time"

    # 检查是否已处理完成，避免重复计算
    if [ -e ${output_folder_path}/${sample_name}.aligned.duplicates_marked.recalibrated.bam ]; then
        echo "${sample_name}.aligned.duplicates_marked.recalibrated.bam exists"
    else
        echo "${sample_name}.aligned.duplicates_marked.recalibrated.bam not exists"

        read -u6 # 领取令牌，控制并发任务数
        {
            ################## 第一步：生成重校准表 ##################

            gatk BaseRecalibrator \
                -R ${reference_path} \                                # 参考基因组
                -L ${interval_path} \                                  # 分析指定捕获区域
                -I ${file} \                                           # 输入BAM（已标记重复）
                --use-original-qualities \                             # 使用原始质量分值（OQ标签）
                -O ${output_folder_path}/${sample_name}.recal_data.csv \ # 输出：校准数据表
                --known-sites ${known_sites_path_1} \                  # dbSNP变异数据库
                --known-sites ${known_sites_path_2} \                  # Mills Indels数据库
                --known-sites ${known_sites_path_3}                    # 1000G Indels数据库

            ################## 第二步：应用重校准 ##################

            gatk ApplyBQSR \
                -R ${reference_path} \                                 # 参考基因组
                -I ${file} \                                           # 同一输入BAM
                -O ${output_folder_path}/${sample_name}.aligned.duplicates_marked.recalibrated.bam \  # 输出新BAM
                -bqsr ${output_folder_path}/${sample_name}.recal_data.csv \  # 使用上一步生成的校准数据
                --add-output-sam-program-record \                      # 在BAM中添加程序信息记录
                --create-output-bam-md5 \                              # 生成MD5文件，便于后续验证
                --use-original-qualities                                # 使用原始测序质量值

            echo "Successfully recalibrated: ${sample_name}"
            sleep 5
            echo >&6	# 归还令牌
        } &
    fi
done

###################### 等待所有后台任务完成 ######################

wait
exec 6>&-  # 关闭fd6

# 打印结束时间
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "06_BQSR done!"

----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------Base Quality Score Recalibration（基因质量重校准）自动提取上一步结果------------------------------------
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
    sample_name="${filename%.aligned.duplicate_marked.sorted.bam}"

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


----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------pileup自动提取上一步结果-------------------------------------------------------------
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


----------------------------------------------------------------------------------------------------------------------------------------------------
-----------------------------------------------------肿瘤污染度（Contamination）评估分析 -------------------------------------------------------------
