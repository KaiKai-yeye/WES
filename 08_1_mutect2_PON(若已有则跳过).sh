#!/bin/bash
#SBATCH -J 08_Mutect2
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/HCCout/home/zengqianwen/WES_HCC_Normal/script/08_Mutect2.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e /groups/HCCout/home/zengqianwen/WES_HCC_Normal/script/08_Mutect2.e # 把报错结果STDERR保存在哪一个文件

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

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "09_Mutect2 done!"

############### PON_1_Mutect2 ###################
input_folder_path=/groups/HCCout/home/zengqianwen/WES_HCC_Normal/align
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed
output_folder_path=/groups/HCCout/home/zengqianwen/WES_HCC_Normal/pon_vcf

############### PON_1_Mutect2 ###################
mkdir -p ${output_folder_path}

# 获取所有目标文件列表
bam_files=(${input_folder_path}/*.aligned.duplicates_marked.recalibrated.bam)
total=${#bam_files[@]}
count=0

for file in "${bam_files[@]}"; do
    count=$((count+1))
    sample_name=$(basename "${file}" .aligned.duplicates_marked.recalibrated.bam)
    output_file=${output_folder_path}/${sample_name}.pon.vcf

    if [[ -f "${output_file}" ]]; then
        echo "[$count/$total] Skipping ${sample_name} → 已存在 ${output_file}"
        continue
    fi

    echo "[$count/$total] Running Mutect2 on ${sample_name} ..."
    read -u6
    {
        gatk Mutect2 \
            -R ${reference_path} \
            -I ${file} \
            -L ${interval_path} \
            -O "${output_file}" \
            --max-mnp-distance 0
            --native-pair-hmm-threads 4   # 可选：每个任务用4核

        echo "[$count/$total] Finished ${sample_name}"
        echo >&6
    } &
done
wait
exec 6>&-

echo "Current Time: $(date +"%Y-%m-%d %H:%M:%S")"
echo "09_Mutect2 done!"
