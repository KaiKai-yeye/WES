#!/bin/bash
#SBATCH -J Bwa_mem2
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/Bwa_mem2_0730.o   # STDOUT 文件
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/Bwa_mem2_0730.e   # STDERR 文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=8    # 可同时运行的最大作业数（根据资源情况可调）
Pfifo="/tmp/$$.fifo"   # 用当前进程ID作为命名管道名，避免冲突
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=$Nproc; i++)); do
    echo
done >&6

input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/trim
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align
bwa_fa_reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/bwa-mem2/GRCh38.primary_assembly.genome.fa

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "04_Bwa_mem2 start!"

for file_1 in ${input_folder_path}/*_1_val_1.fq.gz; do
    # 获取样本名
    sample_name=$(basename "$file_1" _1_val_1.fq.gz)
    file_2="${input_folder_path}/${sample_name}_2_val_2.fq.gz"
    bam_file="${output_folder_path}/${sample_name}.aligned.sorted.bam"

    current_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Current Time: $current_time"

    if [ -e "$bam_file" ]; then
        echo "$sample_name: aligned.sorted.bam exists, skip."
        continue
    else
        echo "$sample_name: aligned.sorted.bam not exists, start processing."

        lane_num=$(zless "$file_1" | head -n1 | cut -d ":" -f4)
        flow_id=$(zless "$file_1" | head -n1 | cut -d ":" -f3)
        run_id=$(zless "$file_1" | head -n1 | cut -d ":" -f2)
        instrument_id=$(zless "$file_1" | head -n1 | cut -d ":" -f1 | sed s/@//)
        platform_unit=${flow_id}.${run_id}.${lane_num}
        id=${sample_name}.${flow_id}.${lane_num}

        read -u6
        {
            bwa-mem2 mem \
                -t 8 \
                -K 100000000 \
                -Y \
                -R "@RG\\tID:${id}\\tLB:WES\\tSM:${sample_name}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tPM:${instrument_id}" \
                "$bwa_fa_reference_path" \
                "$file_1" "$file_2" \
            | samtools sort -@ 8 \
                -m 1500M \
                --write-index \
                -o "$bam_file"

            sleep 5
            echo >&6
        } &
    fi
done

wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "04_Bwa_mem2 done!"
