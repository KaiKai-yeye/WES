#!/bin/bash
#SBATCH -J 11_FilterMutectCalls
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/11_FilterMutectCalls.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/11_FilterMutectCalls.e

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=16    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for((i=1; i<=$Nproc; i++)); do
    echo
done >&6

input_folder_path_vcf=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/unfiltered
input_folder_path_LROM=/groups/g5840141/home/zengqianwen/WES_2025/LROM
input_folder_path_contamination=/groups/g5840141/home/zengqianwen/WES_2025/contamination
input_folder_path_tumor_segmentation=/groups/g5840141/home/zengqianwen/WES_2025/tumor-segmentation
input_folder_path_bam=/groups/g5840141/home/zengqianwen/WES_2025/align
output_folder_path_filtered_1=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_1
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa

# 判断输出目录是否存在，不存在则创建
if [ ! -d "$output_folder_path_filtered_1" ]; then
    echo "Output directory not found, creating: $output_folder_path_filtered_1"
    mkdir -p "$output_folder_path_filtered_1"
fi

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"

#### 进度条相关 ####
total_samples=$(wc -l < 08_mutect2_Sample_list.txt)
completed=0

print_progress() {
    local width=50  # 进度条宽度
    local percent=$(( completed * 100 / total_samples ))
    local filled=$(( completed * width / total_samples ))
    local empty=$(( width - filled ))

    printf "\r[%-${width}s] %d%% (%d/%d)" \
        "$(printf '#%.0s' $(seq 1 $filled))" \
        "$percent" "$completed" "$total_samples"
}

#### 主循环 ####
while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    vcf_file=${input_folder_path_vcf}/${Tumor_sample_name}.unfiltered.vcf
    sample_name=$(basename ${vcf_file} .unfiltered.vcf)

    LROM_file=${input_folder_path_LROM}/${sample_name}.read-orientation-model.tar.gz
    contamination_file=${input_folder_path_contamination}/${sample_name}.contamination.table
    tumor_segmentation_file=${input_folder_path_tumor_segmentation}/${sample_name}.segments.table
    tumor_bam_file=${input_folder_path_bam}/${sample_name}.aligned.duplicates_marked.recalibrated.bam

    if [ -e ${output_folder_path_filtered_1}/${sample_name}.filtered_1.vcf ]; then
        ((completed++))
        print_progress
    else
        read -u6   # 占用一个令牌

        {
            gatk FilterMutectCalls \
                -R ${reference_path} \
                -V ${vcf_file} \
                -O ${output_folder_path_filtered_1}/${sample_name}.filtered_1.vcf \
                --ob-priors ${LROM_file} \
                --contamination-table ${contamination_file} \
                --tumor-segmentation ${tumor_segmentation_file}

            ((completed++))
            print_progress
            echo >&6   # 归还令牌
        }&
    fi
done < 08_mutect2_Sample_list.txt

wait
exec 6>&-
echo    # 换行，避免最后一行进度条和日志混在一起

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "11_FilterMutectCalls done!"
