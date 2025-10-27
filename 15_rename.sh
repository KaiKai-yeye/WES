#!/bin/bash
#SBATCH -J 15_Rename
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -o 15_/groups/g5840141/home/zengqianwen/WES_2025/shell_script/Rename.o 
#SBATCH -e 15_/groups/g5840141/home/zengqianwen/WES_2025/shell_script/Rename.e

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=16    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 6<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作
# 在fd6中放置$Nproc个空行作为令牌
for((i=1; i<=$Nproc; i++)); do
	echo
done >&6


input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/annovar/data
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/annovar/Rename_data

while IFS=$'\t' read -r orignal_name final_name; do
    file=${input_folder_path}/${orignal_name}.hg38_multianno.txt
    if [ -f "${file}" ]; then  # 检查文件是否存在
        cp "${file}" "${output_folder_path}/${final_name}.hg38_multianno.txt"
    else
        echo "File ${file} not found"
    fi
done < hccgroup_rename.txt
