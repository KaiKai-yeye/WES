#!/bin/bash
#SBATCH -J Pileup_0807                           # 作业名称
#SBATCH -N 1                                   # 所需节点数
#SBATCH -n 40                                 # 请求CPU核心数
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/7_Pileup_01.o  # 标准输出
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/7_Pileup_01.e  # 错误输出
############################ 环境激活 ############################
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/
############################ 控制并发进程数量 ############################
Nproc=10
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=Nproc; i++)); do
    echo
done >&6
############################ 设置路径 ############################
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/pileup
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
mkdir -p "$output_folder_path"

# 新增：定义错误标记目录，用于并发错误统计
error_flag_dir="${output_folder_path}/.error_flags"
# 清理上一次运行时可能残留的错误标记
rm -rf "$error_flag_dir"
mkdir -p "$error_flag_dir"

############################ 记录开始时间 ############################
echo "Start Time: $(date +"%Y-%m-%d %H:%M:%S")"

############################ 1. 预扫描，确定要处理的样本列表 ############################
echo "--- Pre-scan: Determining which samples to run ---"
declare -a files_to_process

for file in "$input_folder_path"/*.aligned.duplicates_marked.recalibrated.bam; do
    # 检查文件是否存在，防止通配符无匹配时报错
    [ -e "$file" ] || continue
    filename=$(basename "$file")
    sample_name="${filename%.aligned.duplicates_marked.recalibrated.bam}"
    output_table="${output_folder_path}/${sample_name}.pileups.table"

    if [ -e "$output_table" ]; then
        echo "[SKIP] Output exists for sample: $sample_name"
    else
        files_to_process+=("$file")
    fi
done

total_to_run=${#files_to_process[@]}
echo "----------------------------------------------------"
echo "Scan complete. Found $total_to_run samples to process in this run."
echo "----------------------------------------------------"

if [ $total_to_run -eq 0 ]; then
    echo "All samples have already been processed. Exiting."
    exit 0
fi

############################ 2. 遍历待处理列表并执行Pileup ############################
processed_count=0

for file in "${files_to_process[@]}"; do
    processed_count=$((processed_count + 1))
    
    filename=$(basename "$file")
    sample_name="${filename%.aligned.duplicates_marked.recalibrated.bam}"
    output_table="${output_folder_path}/${sample_name}.pileups.table"

    percentage=$((processed_count * 100 / total_to_run))
    printf "\rProgress: [%-40s] %3d%% (%d/%d) | Submitting: %s " \
    "$(printf '%.0s#' $(seq 1 $((percentage * 40 / 100))))" \
    "$percentage" \
    "$processed_count" \
    "$total_to_run" \
    "$sample_name"

    read -u6
    {
        if ! gatk GetPileupSummaries \
            -R "$reference_path" \
            -I "$file" \
            -V "$germline_resource_path" \
            -L "$interval_path" \
            -O "$output_table" > "${output_table}.log" 2>&1; then
            # 命令失败，在终端打印错误并创建一个错误标记文件
            echo -e "\n[ERROR] GetPileupSummaries failed for sample: $sample_name. Check log: ${output_table}.log" >&2
            touch "${error_flag_dir}/${sample_name}.error"
        fi
        echo >&6
    } &
done

wait
exec 6>&-

# 进度条结束后打印一个换行符
echo
echo "----------------------------------------------------"

############################ 3. 最终错误检查与总结 ############################
# 统计错误标记文件的数量
error_count=$(ls -1q "$error_flag_dir" | wc -l)

if [ "$error_count" -gt 0 ]; then
    echo "<d83d><dd25><d83d><dd25><d83d><dd25> WARNING: A total of $error_count sample(s) failed during processing. <d83d><dd25><d83d><dd25><d83d><dd25>"
    echo "Failed samples list:"
    # 打印出具体是哪些样本失败了
    for flag in "$error_flag_dir"/*.error; do
        echo "  - $(basename "$flag" .error)"
    done
    echo "Please check the corresponding .log files in the output directory for details."
else
    echo "✅✅✅ SUCCESS: All $total_to_run samples processed successfully without any errors. ✅✅✅"
fi

# 清理错误标记目录
rm -rf "$error_flag_dir"

echo "----------------------------------------------------"
echo "End Time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "07_Pileup script finished."
