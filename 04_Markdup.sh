"""
标记但不移除PCR扩增过程中产生的或测序平台造成的重复序列
通常两个（或多个）read 被认为是重复的标准是：
                            具有相同的参考序列起始坐标
                            具有相同的方向（正链或负链）
                            同一个 read pair 的两端都符合上述条件（用于 PE 数据）
##########################################################################################################################
                                    sambamba markdup \
            -t 16 \
            --tmpdir="${temp_dir_path}" \    # 指定临时文件存放目录，避免临时文件占用默认系统目录空间
            --hash-table-size=262144 \       # 设置内部哈希表大小，用于加速重复片段检测，单位是哈希桶数
            --overflow-list-size=67108864 \  # 设置溢出列表大小，单位是字节，用于处理哈希冲突时存储数据，增大可减少内存溢出错误
            "$file" \
            "$output_bam"
###########################################################################################################################          

 
 sambamba markdup 执行后，其退出状态码保存在 $? 中：
0 表示成功，非 0 表示失败
if [[ $? -eq 0 ]] 检查成功与否

echo "[ERROR] ..." 方便你在日志中快速定位失败的样本
"""
#!/bin/bash
#SBATCH -J Markdup_0729
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/5_Markdup_0729.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/5_Markdup_0729.e

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

# ----------------------------
# 并行控制
# ----------------------------
Nproc=8
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
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align
temp_dir_path=/groups/g5840141/home/zengqianwen/WES_2025/tmp
mkdir -p "${temp_dir_path}" || { echo "[ERROR] 无法创建临时目录: ${temp_dir_path}"; exit 1; }

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
            --tmpdir="${temp_dir_path}" \
            --hash-table-size=262144 \
            --overflow-list-size=67108864 \
            "$file" \
            "$output_bam"

        if [[ $? -eq 0 ]]; then
            echo "[DONE] $sample_name duplicate marking completed."
        else
            echo "[ERROR] $sample_name duplicate marking failed with exit code $?"
        fi

        sleep 5
        echo >&6  # 归还令牌
    } &
done

# 等待所有后台任务完成
wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "[END] Current Time: $current_time"
echo "05_Markdup done!"
