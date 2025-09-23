"""
在肿瘤突变检测（如 Mutect2）中，一个常见的问题是 “read orientation artifacts”（测序方向偏倚伪影）。

这些伪影通常来自 FFPE 样本（福尔马林固定-石蜡包埋），在建库或测序时产生。

典型特征是 某些碱基替换只在特定方向的测序 reads 中出现（比如正向 reads 上很多 C>A，但反向 reads 上没有）。

如果不处理，就会导致 假阳性突变。

📌 工具作用

gatk LearnReadOrientationModel 的功能是：

学习/建模 测序数据中的 read orientation bias。

输入是 CollectF1R2Counts 或 Mutect2 生成的 F1R2 统计文件（*.f1r2.tar.gz）。

输出是一个 read orientation model 文件（*.read-orientation-model.tar.gz），里面包含了该样本中可能出现的方向性伪影模式。
"""

#!/bin/bash
#SBATCH -J 10_LearnReadOrientationModel
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/LearnReadOrientationModel.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/LearnReadOrientationModel.e

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=16
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo

for((i=1; i<=$Nproc; i++)); do
    echo
done >&6

#### 路径和文件设置 ####
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/f1r2
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/LROM
sample_list="08_mutect2_Sample_list.txt"
progress_log="10_LearnReadOrientationModel_progress.log"

# 检查样本列表文件是否存在
if [ ! -f "$sample_list" ]; then
    echo "错误: 样本列表文件 $sample_list 不存在!"
    exit 1
fi

#### 初始化计数器和日志 ####
total_samples=$(wc -l < "$sample_list")
processed=0
skipped=0
completed=0
failed=0

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "[$current_time] 开始处理 LearnReadOrientationModel 任务" > $progress_log
echo "总样本数: $total_samples" | tee -a $progress_log
echo "并行进程数: $Nproc" | tee -a $progress_log
echo "=========================================" | tee -a $progress_log

#### 函数定义 (在使用前定义) ####

# 函数：显示进度条
show_progress() {
    local current=$1
    local total=$2
    local percent=$((current * 100 / total))
    local filled=$((percent / 2))
    local empty=$((50 - filled))
    
    printf "\r进度: ["
    printf "%${filled}s" | tr ' ' '='
    printf "%${empty}s" | tr ' ' '-'
    printf "] %d%% (%d/%d)" $percent $current $total
}

# 函数：更新统计信息
update_stats() {
    local status=$1
    case $status in
        "completed")
            ((completed++))
            ;;
        "skipped")
            ((skipped++))
            ;;
        "failed")
            ((failed++))
            ;;
    esac
    ((processed++))
    
    show_progress $processed $total_samples
    
    if [ $((processed % 10)) -eq 0 ] || [ $processed -eq $total_samples ]; then
        printf "\n当前统计: 已完成=%d, 已跳过=%d, 失败=%d, 总进度=%d/%d\n" \
               $completed $skipped $failed $processed $total_samples
        
        current_time=$(date +"%Y-%m-%d %H:%M:%S")
        echo "[$current_time] 进度更新: 已完成=$completed, 已跳过=$skipped, 失败=$failed, 总进度=$processed/$total_samples" >> $progress_log
    fi
}

#### 主循环 ####
while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    # 跳过空行和注释行
    [[ -z "$Tumor_sample_name" || "$Tumor_sample_name" =~ ^#.*$ ]] && continue
    
    output_file="${output_folder_path}/${Tumor_sample_name}.read-orientation-model.tar.gz"
    
    # 检查输出文件是否已存在
    if [ -e "$output_file" ]; then
        current_time=$(date +"%Y-%m-%d %H:%M:%S")
        echo "[$current_time] 跳过已存在文件: ${Tumor_sample_name}.read-orientation-model.tar.gz"
        echo "[$current_time] 跳过: $Tumor_sample_name (文件已存在)" >> $progress_log
        update_stats "skipped"
    else
        # 如果输出文件不存在，则检查输入文件
        tumor_file="${input_folder_path}/${Tumor_sample_name}.f1r2.tar.gz"
        
        if [ ! -e "$tumor_file" ]; then
            current_time=$(date +"%Y-%m-%d %H:%M:%S")
            echo "[$current_time] 警告: 输入文件不存在: $tumor_file"
            echo "[$current_time] 失败: $Tumor_sample_name (输入文件不存在)" >> $progress_log
            update_stats "failed"
            continue # 继续处理下一个样本
        fi
        
        # 如果输入文件存在，则开始处理
        current_time=$(date +"%Y-%m-%d %H:%M:%S")
        echo "[$current_time] 开始处理: ${Tumor_sample_name}"
        echo "[$current_time] 开始处理: $Tumor_sample_name" >> $progress_log
        
       read -u6 # 领取一个令牌
        {
            start_time=$(date +"%Y-%m-%d %H:%M:%S")
            echo "[$start_time] 正在运行 GATK LearnReadOrientationModel: $Tumor_sample_name" >> $progress_log
            
            # 运行 GATK 命令并根据退出状态判断成功或失败
            if gatk LearnReadOrientationModel \
                -I "${tumor_file}" \
                -O "${output_file}" >> "${progress_log}" 2>&1; then
                
                end_time=$(date +"%Y-%m-%d %H:%M:%S")
                echo "[$end_time] 成功完成: $Tumor_sample_name" >> $progress_log
                update_stats "completed"
            else
                end_time=$(date +"%Y-%m-%d %H:%M:%S")
                echo "[$end_time] 处理失败: $Tumor_sample_name" >> $progress_log
                update_stats "failed"
            fi
            
            sleep 2
            echo >&6 # 归还令牌
        } &
    fi # 这是与 `if [ -e "$output_file" ]` 对应的 `fi`
    
done < "$sample_list"

#### 等待与收尾 ####
echo -e "\n等待所有后台任务完成..."
wait

exec 6>&- # 关闭文件描述符

# 最终统计
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo ""
echo "========================================="
echo "[$current_time] LearnReadOrientationModel 任务完成!"
echo "最终统计:"
echo "  - 总样本数: $total_samples"
echo "  - 成功完成: $completed"
echo "  - 跳过文件: $skipped"
echo "  - 失败任务: $failed"

if [ $((total_samples - skipped)) -gt 0 ]; then
    success_rate=$((completed * 100 / (total_samples - skipped)))
    echo "  - 成功率: ${success_rate}%"
else
    echo "  - 成功率: N/A (所有文件都被跳过)"
fi

# 记录最终统计到日志
echo "=========================================" >> $progress_log
echo "[$current_time] 任务完成" >> $progress_log
echo "最终统计:" >> $progress_log
echo "  - 总样本数: $total_samples" >> $progress_log
echo "  - 成功完成: $completed" >> $progress_log
echo "  - 跳过文件: $skipped" >> $progress_log
echo "  - 失败任务: $failed" >> $progress_log
if [ $((total_samples - skipped)) -gt 0 ]; then
    success_rate=$((completed * 100 / (total_samples - skipped)))
    echo "  - 成功率: ${success_rate}%" >> $progress_log
else
    echo "  - 成功率: N/A (所有文件都被跳过)" >> $progress_log
fi

echo "详细日志已保存到: $progress_log"
