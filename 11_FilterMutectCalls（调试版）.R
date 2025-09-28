#!/bin/bash
#SBATCH -J 11_FilterMutectCalls
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/11_FilterMutectCalls.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/11_FilterMutectCalls.e

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 调试信息 ####
echo "=== 脚本开始执行 ==="
echo "当前工作目录: $(pwd)"
echo "当前用户: $(whoami)"
echo "Python环境: $(which python)"
echo "GATK版本检查:"
gatk --version || echo "GATK 命令不可用"

#### 并行进程数量控制 ####
Nproc=16    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo

# 初始化令牌池
for((i=1; i<=$Nproc; i++)); do
    echo
done >&6

# 路径定义
input_folder_path_vcf=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/unfiltered
input_folder_path_LROM=/groups/g5840141/home/zengqianwen/WES_2025/LROM
input_folder_path_contamination=/groups/g5840141/home/zengqianwen/WES_2025/contamination
input_folder_path_tumor_segmentation=/groups/g5840141/home/zengqianwen/WES_2025/tumor-segmentation
input_folder_path_bam=/groups/g5840141/home/zengqianwen/WES_2025/align
output_folder_path_filtered_1=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_1
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa

echo "=== 路径检查 ==="
echo "输出目录: $output_folder_path_filtered_1"

# 判断输出目录是否存在，不存在则创建
if [ ! -d "$output_folder_path_filtered_1" ]; then
    echo "输出目录不存在，正在创建: $output_folder_path_filtered_1"
    mkdir -p "$output_folder_path_filtered_1"
    if [ $? -eq 0 ]; then
        echo "目录创建成功"
    else
        echo "目录创建失败，错误码: $?"
        exit 1
    fi
else
    echo "输出目录已存在"
fi

# 验证目录是否真的创建成功
if [ -d "$output_folder_path_filtered_1" ]; then
    echo "确认输出目录存在: $output_folder_path_filtered_1"
    ls -la "$output_folder_path_filtered_1"
else
    echo "错误：输出目录不存在！"
    exit 1
fi

echo "=== 输入文件检查 ==="
# 检查样本列表文件
if [ ! -f "08_mutect2_Sample_list.txt" ]; then
    echo "错误：样本列表文件不存在: 08_mutect2_Sample_list.txt"
    echo "当前目录内容:"
    ls -la
    exit 1
else
    echo "样本列表文件存在，内容:"
    cat 08_mutect2_Sample_list.txt
    echo "总样本数: $(wc -l < 08_mutect2_Sample_list.txt)"
fi

# 检查关键输入目录
echo "检查输入目录:"
for dir in "$input_folder_path_vcf" "$input_folder_path_LROM" "$input_folder_path_contamination" "$input_folder_path_tumor_segmentation"; do
    if [ -d "$dir" ]; then
        echo "✓ $dir 存在 ($(ls "$dir" | wc -l) 个文件)"
    else
        echo "✗ $dir 不存在"
    fi
done

# 检查参考基因组
if [ -f "$reference_path" ]; then
    echo "✓ 参考基因组存在: $reference_path"
else
    echo "✗ 参考基因组不存在: $reference_path"
fi

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "开始时间: $current_time"

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

echo "=== 开始处理样本 ==="
#### 主循环 ####
line_number=0
while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    ((line_number++))
    echo ""
    echo "--- 处理第 $line_number 行 ---"
    echo "Tumor样本: $Tumor_sample_name"
    echo "Normal样本: $Normal_sample_name"
    
    vcf_file=${input_folder_path_vcf}/${Tumor_sample_name}.unfiltered.vcf
    sample_name=$(basename ${vcf_file} .unfiltered.vcf)
    LROM_file=${input_folder_path_LROM}/${sample_name}.read-orientation-model.tar.gz
    contamination_file=${input_folder_path_contamination}/${sample_name}.contamination.table
    tumor_segmentation_file=${input_folder_path_tumor_segmentation}/${sample_name}.segments.table
    
    echo "样本名: $sample_name"
    echo "VCF文件: $vcf_file"
    echo "LROM文件: $LROM_file"
    echo "污染文件: $contamination_file"
    echo "分割文件: $tumor_segmentation_file"
    
    # 检查所有必需的输入文件
    missing_files=0
    for file in "$vcf_file" "$LROM_file" "$contamination_file" "$tumor_segmentation_file" "$reference_path"; do
        if [ ! -f "$file" ]; then
            echo "✗ 缺少文件: $file"
            ((missing_files++))
        else
            echo "✓ 文件存在: $file"
        fi
    done
    
    if [ $missing_files -gt 0 ]; then
        echo "跳过 $sample_name - 缺少 $missing_files 个必需文件"
        ((completed++))
        print_progress
        continue
    fi
    
    output_file=${output_folder_path_filtered_1}/${sample_name}.filtered_1.vcf
    echo "输出文件: $output_file"

  if [ -e "$output_file" ]; then
        echo "跳过 $sample_name - 输出文件已存在"
        ((completed++))
        print_progress
    else
        read -u6   # 占用一个令牌
        {
            echo "开始处理 $sample_name..."
            echo "执行命令:"
            echo "gatk FilterMutectCalls -R $reference_path -V $vcf_file -O $output_file --ob-priors $LROM_file --contamination-table $contamination_file --tumor-segmentation $tumor_segmentation_file"
            
            gatk FilterMutectCalls \
                -R ${reference_path} \
                -V ${vcf_file} \
                -O ${output_file} \
                --ob-priors ${LROM_file} \
                --contamination-table ${contamination_file} \
                --tumor-segmentation ${tumor_segmentation_file}
            
            exit_code=$?
            if [ $exit_code -eq 0 ]; then
                echo "✓ 成功处理 ${sample_name}"
                if [ -f "$output_file" ]; then
                    echo "✓ 输出文件已创建: $output_file"
                    ls -la "$output_file"
                else
                    echo "✗ 警告：命令执行成功但输出文件不存在"
                fi
            else
                echo "✗ 处理 ${sample_name} 失败，错误码: $exit_code"
            fi
            
            ((completed++))
            print_progress
            echo >&6   # 归还令牌
        }&
    fi
done < 08_mutect2_Sample_list.txt

wait
exec 6>&-

echo ""
echo "=== 最终检查 ==="
if [ -d "$output_folder_path_filtered_1" ]; then
    echo "输出目录内容:"
    ls -la "$output_folder_path_filtered_1"
    echo "输出文件总数: $(ls "$output_folder_path_filtered_1" | wc -l)"
else
    echo "输出目录不存在！"
fi

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "结束时间: $current_time"
echo "11_FilterMutectCalls done!"
