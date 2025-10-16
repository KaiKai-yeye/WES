#!/bin/bash
#SBATCH -J 13_Get_PASS
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/13_Get_PASS.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/13_Get_PASS.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=4    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 6<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作
# 在fd6中放置$Nproc个空行作为令牌
for((i=1; i<=$Nproc; i++)); do
        echo
done >&6


output_folder_path_filtered_2=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_2
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_2_passed #记得新建文件夹

mkdir -p $output_folder_path

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"


while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
        vcf_file=${output_folder_path_filtered_2}/${Tumor_sample_name}.filtered_2.vcf
        sample_name=$(basename ${vcf_file} .filtered_2.vcf)

    read -u6        # 领取令牌

                     {
                        awk '/^#/ {print $0; next} $7=="PASS" \
                                {print $0}' ${vcf_file} > ${output_folder_path}/${sample_name}.passed.vcf
                        sleep 5
                        echo >&6
                        sleep 5
                        echo >&6
                    }&

done < 08_mutect2_Sample_list.txt
wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "13_Get_PASS done!"
