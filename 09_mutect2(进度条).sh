#!/bin/bash
#SBATCH -J 09_Mutect2
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/09_Mutect2.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/09_Mutect2.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=6    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 6<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作

# 在fd6中放置$Nproc个空行作为令牌
for((i=1; i<=$Nproc; i++)); do
        echo
done >&6

input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align              # BAM 文件目录（已经对齐、去重、重校准）
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/unfiltered # 得新建文件夹， .unfiltered.vcf 文件，保存所有未过滤的候选体细胞突变（原始变异）
output_folder_path_f1r2=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/f1r2  # 记得新建文件夹 .f1r2.tar.gz 文件，Mutect2 为 捕获链特异性错误信号生成的辅助数据，下一步 LearnReadOrientationModel 会使用这
些 
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa  # 人类参考基因组（GRCh38.p13，来自 GENCODE）
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed   # 外显子靶向捕获区域的 BED 文件（这里是 Agilent S07604514 panel）
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz  # VCF 文件，包含了 人群中已知的常见生殖系变异，帮助 Mutect2 区分体细胞变异与生殖系变异
pon_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz  # Panel of Normals（PON），从多个正常样本中构建的背景变异数据库，去除在多个正常样本中复现的
假阳性变异（比如文库构建造成的错配），需要根据具体实验进行具体构建


##################### 进度显示设置 #####################
# 统计总样本数
total_samples=0
existing_samples=0
while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    ((total_samples++))
    if [ -e ${output_folder_path}/${Tumor_sample_name}.unfiltered.vcf ]; then
        ((existing_samples++))
    fi
done < list_rerun.txt

# 计算需要处理的样本数
samples_to_process=$((total_samples - existing_samples))

echo "========================================"
echo "Mutect2 样本统计信息："
echo "总样本数: $total_samples"
echo "已完成样本数: $existing_samples"
echo "需要处理样本数: $samples_to_process"
echo "并行进程数: $Nproc"
echo "========================================"

# 进度显示变量
processed_count=0
progress_lock="/tmp/mutect2_progress_lock_$$"

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
    echo "所有样本的unfiltered.vcf文件都已存在，无需处理！"
    exit 0
fi

echo "开始Mutect2变异检测..."

while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    tumor_file=${input_folder_path}/${Tumor_sample_name}.aligned.duplicates_marked.recalibrated.bam
    matched_normal_file=${input_folder_path}/${Normal_sample_name}.aligned.duplicates_marked.recalibrated.bam
    
    if [ -e ${output_folder_path}/${Tumor_sample_name}.unfiltered.vcf ]; then
        echo "${Tumor_sample_name}.unfiltered.vcf exists"
    else
        echo "${Tumor_sample_name}.unfiltered.vcf not exists"  
        read -u6        # 领取令牌
        {
            # 记录单个样本开始时间
            sample_start_time=$(date +"%Y-%m-%d %H:%M:%S")
            echo "[$sample_start_time] 开始处理样本: ${Tumor_sample_name} vs ${Normal_sample_name}"
            
            gatk Mutect2 \
                --native-pair-hmm-threads 8 \
                -R ${reference_path} \
                --germline-resource ${germline_resource_path} \
                -pon ${pon_path} \
                -L ${interval_path} \
                -I ${tumor_file} \
                -I ${matched_normal_file} \
                -normal ${Normal_sample_name} \
                --f1r2-tar-gz ${output_folder_path_f1r2}/${Tumor_sample_name}.f1r2.tar.gz \
                -O ${output_folder_path}/${Tumor_sample_name}.unfiltered.vcf
            
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
echo "Mutect2处理完成！共处理了 $samples_to_process 个样本"
echo "输出文件保存在:"
echo "  - VCF文件: $output_folder_path"
echo "  - F1R2文件: $output_folder_path_f1r2"
echo "========================================"
echo "09_Mutect2 done!"
