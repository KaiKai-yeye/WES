#!/bin/bash
#SBATCH -J 08_cretepon
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o /groups/HCCout/home/zengqianwen/WES_HCC_Normal/script/08_createpon.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e /groups/HCCout/home/zengqianwen/WES_HCC_Normal/script/08_createpon.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

input_folder_path=/groups/HCCout/home/zengqianwen/WES_HCC_Normal/pon_vcf
output_folder_path=/groups/HCCout/home/zengqianwen/WES_HCC_Normal/pon
output_pon_name=PON_HCC.vcf
pon_database_path=/groups/HCCout/home/zengqianwen/WES_HCC_Normal/pon/database
temp_dir_path=/groups/HCCout/home/zengqianwen/WES_HCC_Normal/pon/tmp
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz

# 创建必要的目录
mkdir -p ${output_folder_path}
mkdir -p ${temp_dir_path}

echo "Starting PON creation process..."
current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"

# 检查输入VCF文件数量
vcf_count=$(ls ${input_folder_path}/*.pon.vcf 2>/dev/null | wc -l)
echo "Found ${vcf_count} PON VCF files to process"

if [ ${vcf_count} -eq 0 ]; then
    echo "Error: No PON VCF files found in ${input_folder_path}"
    exit 1
fi

# 删除已存在的database目录（如果存在）
if [ -d "${pon_database_path}" ]; then
    echo "Removing existing database directory: ${pon_database_path}"
    rm -rf ${pon_database_path}
fi
echo "Step 1: Creating GenomicsDB workspace..."
gatk --java-options "-Xms32G -Xmx64G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport \
        $(for file in ${input_folder_path}/*.pon.vcf; do
    printf " \055V %s" ${file}
    done) \
        --genomicsdb-workspace-path ${pon_database_path} \
        --tmp-dir ${temp_dir_path} \
        --merge-input-intervals true \
        -L ${interval_path} \
        --max-num-intervals-to-import-in-parallel 32 \
        --batch-size 50

if [ $? -eq 0 ]; then
    echo "GenomicsDBImport completed successfully"
else
    echo "Error: GenomicsDBImport failed"
    exit 1
fi

echo "Step 2: Creating Somatic Panel of Normals..."
gatk CreateSomaticPanelOfNormals \
        -R ${reference_path} \
        --germline-resource ${germline_resource_path} \
        -V gendb://${pon_database_path} \
        -O ${output_folder_path}/${output_pon_name}

if [ $? -eq 0 ]; then
    echo "CreateSomaticPanelOfNormals completed successfully"
    echo "PON file created: ${output_folder_path}/${output_pon_name}"
else
    echo "Error: CreateSomaticPanelOfNormals failed"
    exit 1
fi

# 清理临时文件
echo "Cleaning up temporary files..."
rm -rf ${temp_dir_path}

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "PON creation process completed!"
~                                          
