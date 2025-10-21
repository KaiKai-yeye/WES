"""
通过Mutect2（一款体细胞突变检测工具）对体细胞突变VCF文件进行批量注释。

脚本的核心功能分为三个步骤：
    格式转换：将 .vcf 文件转换为 ANNOVAR 特定的 .avinput 格式。
    执行注释：使用 refGene 数据库对 .avinput 文件进行基因功能注释。
    结果合并：从每个样本的注释结果中提取前10列关键信息，并添加样本名，最后合并成一个总的 all_sample.csv 文件。
"""

#!/bin/bash
#SBATCH -J 14_ANNOVAR
#SBATCH -N 1
#SBATCH -n 48
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/14_ANNOVAR.o
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/14_ANNOVAR.e

# 激活Anaconda基础环境
source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

# 定义输入和输出路径及数据库路径
input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/mutect2/vcf/filtered_2_passed
output_folder_path1=/groups/g5840141/home/zengqianwen/WES_2025/annovar/mutect2
output_folder_path2=/groups/g5840141/home/zengqianwen/WES_2025/annovar/data
reference_path=/groups/g5840087/home/share/annovar/humandb/

# 打印开始时间
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "ANNOVAR start!"

## 1. 转换格式为avinput
for file in ${input_folder_path}/*.passed.vcf; do
    sample_name=$(basename ${file} .passed.vcf)
    
    perl /groups/g5840087/home/share/annovar/convert2annovar.pl \
        -format vcf4 \
        -includeinfo \
        -allsample \
        -withfreq \
        -filter pass \
        ${file} \
        -o ${output_folder_path1}/${sample_name}.passed_vcf.avinput
done

        #-format vcf4 \
        # 指定输入文件的格式为VCF v4.x
        #-includeinfo \
        # 在avinput文件中包含VCF文件的INFO字段信息
        #-allsample \
        # 对VCF文件中的所有样本（如果存在多个）进行处理
        #-withfreq \
        # 在avinput文件中包含等位基因频率信息
        #-filter pass \
        # 仅处理VCF文件中FILTER列标记为'PASS'的变异位点
        #${file} \
        # 指定输入的VCF文件
        #-o ${output_folder_path1}/${sample_name}.passed_vcf.avinput
        # 指定输出的avinput文件路径和名称

## 2. 注释
for file in ${output_folder_path1}/*.passed_vcf.avinput; do
    sample_name=$(basename ${file} .passed_vcf.avinput)
    
    perl /groups/g5840087/home/share/annovar/table_annovar.pl \
        ${file} \
        ${reference_path} \
        -buildver hg38 \
        -out ${output_folder_path2}/${sample_name} \
        -protocol refGene \
        -operation g \
        -nastring . \
        -remove
done
        #${file} \
        # 输入文件，即上一步生成的.avinput文件
        #${reference_path} \
        # ANNOVAR数据库的根目录
        #-buildver hg38 \
        # 指定基因组构建版本为hg38
        #-out ${output_folder_path2}/${sample_name} \
        # 指定输出文件的前缀名
        #-protocol refGene \
        # 指定要使用的注释协议/数据库（此处为refGene，用于基因注释）
        #-operation g \
        # 指定注释操作类型：'g' 代表 gene-based（基于基因的注释）
        #-nastring . \
        # 指定如何表示缺失值（Not Available），此处使用点号'.'
        #-remove
        # 在注释完成后，删除生成的临时文件

## 3. 合并注释文件
for file in ${output_folder_path2}/*.hg38_multianno.txt; do
    sample_name=$(basename ${file} .hg38_multianno.txt)
    cut -f '1-10' ${file} | sed '1d' | sed "s/$/,${sample_name}/" >> ${output_folder_path2}/all_sample.csv
done

# 打印结束时间
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "ANNOVAR done!"
