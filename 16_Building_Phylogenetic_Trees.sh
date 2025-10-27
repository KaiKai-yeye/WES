#!/bin/bash
#SBATCH -J 16_Building_Phylogenetic_Trees
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o /groups/g5840141/home/zengqianwen/WES_2025/shell_script/16_Building_Phylogenetic_Trees.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e /groups/g5840141/home/zengqianwen/WES_2025/shell_script/16_Building_Phylogenetic_Trees.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/


input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/annovar/Rename_data
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/annovar/fa_out #记得新建文件夹
reference_path=/groups/g5840141/home/zengqianwen/WES/reference/gencode/release_40/GRCh38.p13/fa/hg38.fa
sample_list=(HCC001 HCC003 HCC004 HCC005 HCC006 HCC007 HCC008 HCC009 HCC010 HCC011 HCC012 HCC013 HCC014
             HCC015 HCC016 HCC017 HCC018 HCC019 HCC020 HCC021 HCC022 HCC023 HCC024 HCC025 HCC026 HCC027 
             HCC028 HCC029 HCC030 HCC031 HCC032 HCC033 HCC034 HCC035 HCC036 HCC037 HCC038 HCC039 HCC040 
             HCC041 HCC042 HCC043 HCC044 HCC045 HCC046 HCC047 HCC048 HCC049 HCC050 HCC051 HCC052 HCC053 
             HCC054 HCC055 HCC056 HCC057 HCC058 HCC059 HCC060 HCC061 HCC062 HCC063 HCC064 HCC065 HCC066 
             HCC067 HCC068 HCC069 HCC070 HCC071 HCC072 HCC073 HCC074 HCC075 HCC076 HCC077 HCC078 HCC079 
             HCC080 HCC081 HCC082 HCC083 HCC084 HCC085 HCC086 HCC087 HCC088 HCC089 HCC090 HCC091 HCC092 
             HCC093 HCC094 HCC095 HCC096)

for sample in ${sample_list[@]}; do
	perl /groups/g5840141/home/zengqianwen/WES/shell_script/consensus_seq.pl \
	${reference_path} \
	${input_folder_path} \
	${sample} >${output_folder_path}/${sample}.fa
done
