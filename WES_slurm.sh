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
--------------------------------------------------------------------------------Bwa_mem2--------------------------------------------------------------
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
