#!/bin/bash
#SBATCH -J Calculate_Contamination_0918
#SBATCH -N 1
#SBATCH -n 36
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/07_Calculate_Contamination_2.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/07_Calculate_Contamination_2.e

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

##################### 并行控制设置 #####################
Nproc=6
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=Nproc; i++)); do echo; done >&6

##################### 目录路径设置 #####################
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/pileup
output_folder_path_tumor_segmentation=/groups/g5840141/home/zengqianwen/WES_2025/tumor-segmentation
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/contamination
mkdir -p "$output_folder_path_tumor_segmentation"
mkdir -p "$output_folder_path"

##################### 进度显示设置 #####################
# 统计总样本数
total_samples=0
existing_samples=0
while read -r Tumor_sample_name Normal_sample_name; do
    ((total_samples++))
    if [ -e ${output_folder_path}/${Tumor_sample_name}.contamination.table ]; then
        ((existing_samples++))
    fi
done < 07_Calculate_Contamination_Sample_list.txt

# 计算需要处理的样本数
samples_to_process=$((total_samples - existing_samples))

echo "========================================"
echo "样本统计信息："
echo "总样本数: $total_samples"
echo "已完成样本数: $existing_samples"
echo "需要处理样本数: $samples_to_process"
echo "========================================"

# 进度显示变量
processed_count=0
progress_lock="/tmp/progress_lock_$$"

# 进度显示函数
show_progress() {
    local current=$1
    local total=$2
    local percentage=$((current * 100 / total))
    local bar_length=50
    local filled_length=$((percentage * bar_length / 100))
    
    printf "\r进度: ["
    for ((i=0; i<filled_length; i++)); do printf "#"; done
    for ((i=filled_length; i<bar_length; i++)); do printf " "; done
    printf "] %d/%d (%d%%)" "$current" "$total" "$percentage"
}

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "开始时间: $current_time"

# 如果所有样本都已存在，直接退出
if [ $samples_to_process -eq 0 ]; then
    echo "所有样本的contamination.table文件都已存在，无需处理！"
    exit 0
fi

echo "开始处理样本..."

while read -r Tumor_sample_name Normal_sample_name; do
    tumor_file=${input_folder_path}/${Tumor_sample_name}.pileups.table
    matched_normal_file=${input_folder_path}/${Normal_sample_name}.pileups.table
    
    if [ -e ${output_folder_path}/${Tumor_sample_name}.contamination.table ]; then
        echo "${Tumor_sample_name}.contamination.table exists"
    else
        echo "${Tumor_sample_name}.contamination.table not exists"
        read -u6        # 领取令牌
        {
            # 记录单个样本开始时间
            sample_start_time=$(date +"%Y-%m-%d %H:%M:%S")
            echo "[$sample_start_time] 开始处理样本: ${Tumor_sample_name}"
            
            gatk CalculateContamination \
            -I ${tumor_file} \
            -matched ${matched_normal_file} \
            -tumor-segmentation ${output_folder_path_tumor_segmentation}/${Tumor_sample_name}.segments.table \
            -O ${output_folder_path}/${Tumor_sample_name}.contamination.table
            
            # 记录单个样本完成时间和更新进度
            sample_end_time=$(date +"%Y-%m-%d %H:%M:%S")
            echo "[$sample_end_time] 完成处理样本: ${Tumor_sample_name}"
            
            # 使用文件锁更新进度计数器
            (
                flock -x 200
                ((processed_count++))
                show_progress $processed_count $samples_to_process
                if [ $processed_count -eq $samples_to_process ]; then
                    echo ""  # 换行
                fi
            ) 200>"$progress_lock"
            
            sleep 5
            echo >&6
        }&
    fi
done < list_rerun.txt

wait
exec 6>&-

# 清理进度锁文件
rm -f "$progress_lock"

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "========================================"
echo "完成时间: $current_time"
echo "处理完成！共处理了 $samples_to_process 个样本"
echo "========================================"
echo "08_Calculate_Contamination done!"
