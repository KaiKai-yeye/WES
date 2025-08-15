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
#      - 同样来自 GetPileupSummaries，来自匹配的正常样本（血液P或肝脏L），复发发没正常样本可以用同一个病人的原发正常
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
#   - tumor segmentation 文件 肿瘤分段文件
#
# ⚙️ 原理流程简述：
#   - GATK 首先根据 tumor 和 normal 的 pileup 数据，在常见种系变异位点上
#     比较等位基因频率，识别可能的污染信号（如在应为纯合位点发现另一种等位基因）；
#   - 最终估算肿瘤样本中的 contamination 值，用于后续过滤（如 FilterMutectCalls）。
#
# 📌 注意事项：
#   - pileup 表格必须来源于相同版本的参考基因组和变异数据库；
#   - 必须由对应样本的拷贝数分析流程（ModelSegments）生成；
#   - 匹配的 normal 样本必须尽可能来自同一患者（如血液P、肝组织L），
#     否则污染率估算可能不准确。
"""
#!/bin/bash
#SBATCH -J Calculate_Contamination_0814
#SBATCH -N 1
#SBATCH -n 36
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/07_Calculate_Contamination_1.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/07_Calculate_Contamination_1.e

# 激活 Conda 环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

##################### 并行控制设置 #####################
Nproc=12
Pfifo="/tmp/$$.fifo"
mkfifo $Pfifo
exec 6<>$Pfifo
rm -f $Pfifo
for ((i=1; i<=Nproc; i++)); do echo; done >&6

##################### 目录路径设置 #####################
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/pileup
output_folder_path_tumor_segmentation=/groups/g5840141/home/zengqianwen/WES_2025/tumor-segmentation
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/contamination
mkdir -p "$output_folder_path_tumor_segmentation"
mkdir -p "$output_folder_path"

# Sample_list=(PHCC2048T PHCC2048L RHCC2358T RHCC2358M PHCC3519T PHCC3519P RHCC4057Th RHCC4057L PHCC1889T PHCC1889L RHCC4055T RHCC4055P PHCC2992T PHCC2992L RHCC4050T RHCC4050P PHCC2071T PHCC2071P RHCC4173T1 RHCC4173T2 RHCC4173L)
#Sample_list=(PHCC2417T1 PHCC2417T2 PHCC2417T3 PHCC2417T4 PHCC2417P RHCC4584T PHCC1011T PHCC1011L RHCC4619T RHCC4619P PHCC966T RHCC4664T RHCC4664P PHCC3603T RHCC4691T RHCC4691P PHCC972T PHCC972L RHCC4349T RHCC3369T RHCC3369P RHCC4121T RHCC4121L)

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"

while read -r Tumor_sample_name Normal_sample_name; do
    tumor_file=${input_folder_path}/${Tumor_sample_name}.pileups.table
    matched_normal_file=${input_folder_path}/${Normal_sample_name}.pileups.table
    
    if [ -e ${output_folder_path}/${Tumor_sample_name}.contamination.table ]; then
        echo "${Tumor_sample_name}.contamination.table exists"
    else
        echo "${Tumor_sample_name}.contamination.table not exists"
        read -u6        # 领取令牌
        {
            gatk CalculateContamination \
            -I ${tumor_file} \
            -matched ${matched_normal_file} \
            -tumor-segmentation ${output_folder_path_tumor_segmentation}/${Tumor_sample_name}.segments.table \
            -O ${output_folder_path}/${Tumor_sample_name}.contamination.table
            sleep 5
            echo >&6
        }&
    fi
done < 07_Calculate_Contamination_Sample_list.txt

wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "08_Calculate_Contamination done!"
