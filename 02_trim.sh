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
