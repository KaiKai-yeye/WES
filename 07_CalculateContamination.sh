"""
# Step: CalculateContamination
#
# 🎯 目的：
#   评估肿瘤样本中混入正常细胞的污染比例（contamination fraction），
#   为后续 somatic mutation 的准确识别提供支持。
#
# 📥 输入文件：
#   1. Tumor pileup table (.pileups.table)
#      - 来源：GetPileupSummaries 步骤
#      - 包含肿瘤样本在公共种系变异位点上的等位基因频率统计信息
#
#   2. Matched normal pileup table (.pileups.table)
#      - 同样来自 GetPileupSummaries，来自匹配的正常样本（血液P或肝脏L）
#
#   3. Tumor segmentation file (.segments.table)
#      - 来源：GATK ModelSegments 步骤（在 CNV 拷贝数变异分析中生成）
#      - 表示肿瘤样本基因组中不同区域的 copy number 状态（如扩增或缺失），用于校正因肿瘤样本中拷贝数异常而造成的等位基因频率偏差。
#        若不考虑这些拷贝数改变区域，可能会错误地将其解释为样本污染，从而降低评估准确性。
#
#      📄 示例结构（TSV）：
#      CONTIG   START      END        NUM_POINTS_COPY_RATIO   MEAN_LOG2_COPY_RATIO（0 表示无变化，负值表示缺失，正值表示扩增）
#      chr1     10000      885000     245                     -0.142
#      chr1     885001     1260000    186                     0.027
#
#      ✅ 含义简述：
#      - 这些区段表示肿瘤样本中可能发生了拷贝数改变（如染色体缺失或扩增）的区域；
#      - CalculateContamination 利用这些信息排除由 copy number 改变导致的等位基因频率偏差，
#        从而更准确地区分真正的样本污染。
#
# 📤 输出文件：
#   - contamination.table：包含污染率估计值
#     主要字段：
#       * contamination：估算的污染比例（如 0.03 表示 3% 混入正常细胞）
#       * contamination_sd：标准差，反映置信度
#
# ⚙️ 原理流程简述：
#   - GATK 首先根据 tumor 和 normal 的 pileup 数据，在常见种系变异位点上
#     比较等位基因频率，识别可能的污染信号（如在应为纯合位点发现另一种等位基因）；
#   - 同时结合肿瘤样本的 segments.table 文件，识别并校正受拷贝数异常影响的区域；
#   - 最终估算肿瘤样本中的 contamination 值，用于后续过滤（如 FilterMutectCalls）。
#
# 📌 注意事项：
#   - pileup 表格必须来源于相同版本的参考基因组和变异数据库；
#   - tumor segmentation 文件必须由对应样本的拷贝数分析流程（ModelSegments）生成；
#   - 匹配的 normal 样本必须尽可能来自同一患者（如血液P、肝组织L），
#     否则污染率估算可能不准确。
"""
#!/bin/bash
#SBATCH -J Calculate_Contamination
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -o Calculate_Contamination.o
#SBATCH -e Calculate_Contamination.e

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

##################### 并行控制设置 #####################
Nproc=4
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=Nproc; i++)); do echo; done >&6

##################### 目录路径设置 #####################
input_folder_path=/groups/g5840141/home/zengqianwen/WES/pileup
input_folder_path_tumor_segmentation=/groups/g5840141/home/zengqianwen/WES/tumor-segmentation
output_folder_path=/groups/g5840141/home/zengqianwen/WES/contamination

mkdir -p "$output_folder_path"

##################### 记录开始时间 #####################
echo "Start Time: $(date '+%Y-%m-%d %H:%M:%S')"

##################### 样本对自动识别处理 #####################
for tumpr_file in ${input_folder_path}/*T*.pileups.table; do
    # 提取肿瘤样本名
    tumor_sample=$(basename "$tumpr_file" .pileups.table)

    # 优先匹配 L，再匹配 P 作为正常样本
    matched_normal_file=""
    for suffix in L P; do
        normal_candidate="${input_folder_path}/${tumor_sample%T*}${suffix}.pileups.table"
        if [[ -f "$normal_candidate" ]]; then
            matched_normal_file="$normal_candidate"
            normal_sample=$(basename "$matched_normal_file" .pileups.table)
            break
        fi
    done

    if [[ -z "$matched_normal_file" ]]; then
        echo "[WARN] No matched normal sample found for tumor sample: $tumor_sample"
        continue
    fi

    # 输出文件路径
    output_file="${output_folder_path}/${tumor_sample}.contamination.table"

    if [[ -e "$output_file" ]]; then
        echo "[SKIP] Contamination already exists for $tumor_sample"
        continue
    fi

    # 检查分段文件是否存在
    segment_file="${input_folder_path_tumor_segmentation}/${tumor_sample}.segments.table"
    if [[ ! -f "$segment_file" ]]; then
        echo "[SKIP] Tumor segmentation file not found for $tumor_sample"
        continue
    fi

    # 执行 contamination 计算
    read -u6
    {
        echo "[RUN] Processing $tumor_sample (normal: $normal_sample) at $(date '+%H:%M:%S')"
        gatk CalculateContamination \
            -I "$tumpr_file" \
            -matched "$matched_normal_file" \
            -tumor-segmentation "$segment_file" \
            -O "$output_file"

        sleep 2
        echo >&6
    } &
done

wait
exec 6>&-

##################### 记录结束时间 #####################
echo "End Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "08_Calculate_Contamination done!"
