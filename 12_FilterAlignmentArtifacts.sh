"""
1. 背景
      在体细胞突变检测（尤其是 WES/WGS）中，有一类假阳性来自 比对错误（alignment artifacts）：
      某些序列很难唯一比对（重复区域、低复杂度序列、同源基因家族）。
      比对算法可能错误地把 reads 放到错误的位置，导致“假突变”。
      特别常见于 Indel、靠近同源序列的 SNV 等区域。
      这些假阳性单靠 FilterMutectCalls 不一定能很好识别，所以 GATK 提供了一个额外的工具 —— FilterAlignmentArtifacts。

2. FilterAlignmentArtifacts 的作用
      该工具会利用 比对质量特征 和 重比对检查，识别可能由于比对错误造成的假突变。
      它主要使用 Panel of Normals (PON) 和 BAM 文件 来判断某些突变是否更可能是 alignment artifact。
      核心逻辑是：
            取出候选变异周围的 reads；
            使用一个 比对器（通常是 BWA 或 minimap2） 重新比对这些 reads 到参考基因组；
            比较新的比对和原始比对 → 如果发现变异支持信号消失/减弱 → 判定该突变可能是假阳性。

3. 输出结果
  输出：一个 VCF，区别是：在 FILTER 字段里追加标记，比如：artifact_in_normal，artifact_in_tumor，orientation_bias（某些场景下也可能在这里标记）
       如果突变被怀疑是假阳性，就会被标记。通过检查的突变仍然是 PASS。最后输出过滤后的VCF
"""
#!/bin/bash
#SBATCH -J 12_FilterAlignmentArtifacts
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/12_FilterAlignmentArtifacts.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/12_FilterAlignmentArtifacts.e

# 激活环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 调试信息 ####
echo "=== 脚本开始执行: FilterAlignmentArtifacts ==="
echo "当前工作目录: $(pwd)"
echo "GATK版本检查:"
gatk --version || echo "GATK 命令不可用"

#### 并行进程数量控制 ####
Nproc=12    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo

# 初始化令牌池
for((i=1; i<=$Nproc; i++)); do
    echo
done >&6

# --- 路径定义 (已根据上一步的输出进行更新) ---
input_folder_path_filtered_1=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_1
input_folder_path_bam=/groups/g5840141/home/zengqianwen/WES_2025/align
output_folder_path_filtered_2=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_2
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa
#  BWA-MEM2 工具为上面的 FASTA 文件生成的索引镜像文件 (Index Image)
bwa_img_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa.img
SAMPLE_LIST_FILE="08_mutect2_Sample_list.txt" # 使用与上一步一致的样本列表

echo "=== 路径与文件检查 ==="

# 1. 检查并创建输出目录
if [ ! -d "$output_folder_path_filtered_2" ]; then
    echo "输出目录不存在，正在创建: $output_folder_path_filtered_2"
    mkdir -p "$output_folder_path_filtered_2"
    if [ $? -ne 0 ]; then
        echo "错误：输出目录创建失败！"
        exit 1
    fi
else
    echo "输出目录已存在: $output_folder_path_filtered_2"
fi

# 2. 检查样本列表文件
if [ ! -f "$SAMPLE_LIST_FILE" ]; then
    echo "错误：样本列表文件不存在: $SAMPLE_LIST_FILE"
    exit 1
fi

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "开始时间: $current_time"

#### 主循环 ####
while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    # 从肿瘤样本名推导出基础样本名
    sample_name=${Tumor_sample_name}
    
    # 定义输入输出文件路径
    input_vcf=${input_folder_path_filtered_1}/${sample_name}.filtered_1.vcf
    tumor_bam_file=${input_folder_path_bam}/${sample_name}.aligned.duplicates_marked.recalibrated.bam
    output_vcf=${output_folder_path_filtered_2}/${sample_name}.filtered_2.vcf

    # 检查输出文件是否已存在，如果存在则跳过
    if [ -f "$output_vcf" ]; then
        echo "跳过 $sample_name - 输出文件已存在."
        continue
    fi
    
    # 检查必需的输入文件是否存在
    if [ ! -f "$input_vcf" ]; then
        echo "跳过 $sample_name - 缺少输入文件: $input_vcf"
        continue
    fi
    if [ ! -f "$tumor_bam_file" ]; then
        echo "跳过 $sample_name - 缺少BAM文件: $tumor_bam_file"
        continue
    fi

    read -u6   # 占用一个令牌，实现并行控制
    {
        echo "开始处理 $sample_name..."
        
        gatk FilterAlignmentArtifacts \
            -R ${reference_path} \
            --bwa-mem-index-image ${bwa_img_path} \
            -V ${input_vcf} \
            -I ${tumor_bam_file} \
            -O ${output_vcf}
        
        # 检查GATK命令的退出状态
        exit_code=$?
        if [ $exit_code -eq 0 ]; then
            echo "✓ 成功处理 ${sample_name}"
        else
            echo "✗ 处理 ${sample_name} 失败，错误码: $exit_code"
        fi

        echo >&6  # 归还令牌
    }&
done < "$SAMPLE_LIST_FILE"

wait
exec 6>&-

echo ""
echo "=== 所有样本处理完毕 ==="
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "结束时间: $current_time"
echo "12_FilterAlignmentArtifacts done!"
