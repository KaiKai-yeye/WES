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
