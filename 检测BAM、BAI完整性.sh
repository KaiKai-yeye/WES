#!/bin/bash
#SBATCH -J check_bam
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/check_bam_%j.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/check_bam_%j.e

# 激活 conda 环境 (确保 samtools 已安装)
# 请确保这里的路径是正确的
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

# ==================== 配置 ====================
# 输入 BAM 文件目录
BAM_DIR="/groups/g5840141/home/zengqianwen/WES_2025/align/"
# 定义输出结果文件和错误日志文件
OUTPUT_FILE="${BAM_DIR}/bam_bai_check_results.tsv"
ERROR_LOG="${BAM_DIR}/bam_bai_check_errors.log"

# ==================== 主程序 ====================
echo "======================================"
echo " BAM/BAI check started."
echo " Start Time: $(date)"
echo "======================================"

# 初始化输出文件和错误日志
# 在表头中，我们将第一列改为 "FileName" 以便更清晰地反映其内容
echo -e "FileName\tSampleID\tBAM_Status\tBAI_Status" > "$OUTPUT_FILE"
# 清空之前的错误日志，以便本次运行生成新的日志
> "$ERROR_LOG"

# 增强：在开始循环前，先检查是否存在 .bam 文件
shopt -s nullglob
bam_files=(${BAM_DIR}/*.bam)
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "警告：在目录 ${BAM_DIR} 中没有找到任何 .bam 文件。"
    echo "脚本执行结束。"
    exit 0
fi
shopt -u nullglob # 恢复默认行为

# 遍历所有 bam 文件
for bam in ${BAM_DIR}/*.bam; do
    filename=$(basename "$bam")
    sample_id="${filename%%.*}"
    # 对应的索引文件名逻辑保持不变，因为它是正确的
    bai="${bam}.bai"

    bam_status="OK"
    bai_status="OK"

    # 1. 检查 BAM 文件完整性
    bam_error_msg=$(samtools quickcheck "$bam" 2>&1)
    if [ $? -ne 0 ]; then
        bam_status="BAD"
        echo -e "${filename}\tBAM_ERROR\t${bam_error_msg}" >> "$ERROR_LOG"
    fi

    # 2. 检查 BAI 是否存在且可用
    if [ ! -f "$bai" ]; then
        bai_status="MISSING"
    else
        bai_error_msg=$(samtools idxstats "$bam" >/dev/null 2>&1)
        if [ $? -ne 0 ]; then
            bai_status="BAD"
            echo -e "${filename}\tBAI_ERROR\t${bai_error_msg}" >> "$ERROR_LOG"
        fi
    fi

    # 将检查结果写入主输出文件
    # 同时记录完整文件名和提取出的SampleID，信息更全面
    echo -e "${filename}\t${sample_id}\t${bam_status}\t${bai_status}" >> "$OUTPUT_FILE"
done

echo "======================================"
echo " BAM/BAI check finished."
echo " 结果已保存至: $OUTPUT_FILE"

# 增强：如果错误日志文件不为空，则提示用户查看
if [ -s "$ERROR_LOG" ]; then
    echo " 检测到错误，详细信息请查看: $ERROR_LOG"
else
    rm "$ERROR_LOG"
fi

echo " End Time: $(date)"
echo "======================================"
```这个版本的脚本现在会输出一个包含四列的报告，例如：
`FileName                                                 SampleID      BAM_Status   BAI_Status`
`RHCC2439T1.aligned.duplicates_marked.recalibrated.bam    RHCC2439T1    OK           OK`
`RHCC2439T1.aligned.sorted.bam                            RHCC2439T1    OK           OK`


