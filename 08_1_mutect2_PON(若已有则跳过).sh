"""
输出文件：
1547L.pon.vcf
    这是主要的输出文件，VCF 格式 (Variant Call Format)。
    它记录了在一组正常样本中检测到的可疑突变位点（通常是测序伪影或背景噪音）。
    在后续分析中，这个 PON 文件会作为输入，帮助在实际样本中剔除这些常见伪阳性突变，减少错误检出率。

1547L.pon.vcf.idx
    这是 vcf 文件的 索引文件，通常由 GATK 自动生成（基于 .vcf）。
    索引用于快速随机访问 VCF 文件的特定位置，而不需要每次都加载整个文件。
    在后续 Mutect2 或 FilterMutectCalls 等步骤中，程序会自动用到它。

1547L.pon.vcf.stats
    这是 Mutect2 生成的 统计信息文件。
    它包含该 VCF 文件的相关统计，比如候选突变的分布、不同过滤类别的数量、各种错误模型的参数等。
    通常用于质量评估和调试，而不是直接作为输入。
"""
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
