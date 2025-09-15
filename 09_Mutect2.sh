"""
🧬 Mutect2 是什么？
GATK Mutect2 是用于检测体细胞突变（somatic mutation）的核心工具，尤其适用于：
癌症数据分析
    肿瘤-正常（tumor-normal）配对分析（也支持单肿瘤样本）可用于 WES、WGS、panel 等数据
它基于 HaplotypeCaller 的重组装思想，能检测：
      SNVs（单核苷酸突变）、Indels（插入/缺失）

输入文件
1. BAM 文件（已处理的测序数据）

路径: /groups/g5840141/home/zengqianwen/WES/align/
文件格式: ${样本名}.aligned.duplicates_marked.recalibrated.bam
内容: 已完成比对、去重、重校准的测序数据
用途: 肿瘤样本和配对正常样本的输入数据

2. 参考基因组

路径: GRCh38.primary_assembly.genome.fa
用途: 人类参考基因组序列（GRCh38.p13版本）

3. 外显子捕获区域

路径: S07604514_Padded.bed
用途: 定义外显子测序的目标区域（Agilent S07604514 panel）

4. 生殖系变异资源

路径: small_exac_common_3.hg38.vcf.gz
用途: 包含人群中常见生殖系变异，帮助区分体细胞vs生殖系变异

5. Panel of Normals (PON)

路径: 1000g_pon.hg38.vcf.gz
用途: 正常样本背景变异数据库，过滤假阳性变异

6. 样本配对信息

文件: 08_Calculate_Contamination_Sample_list.txt
格式: 制表符分隔，包含肿瘤样本名和配对正常样本名

输出文件
1. 未过滤的VCF文件

路径: /groups/g5840141/home/zengqianwen/WES/mutect2/vcf/unfiltered/
文件格式: ${肿瘤样本名}.unfiltered.vcf
内容:

所有检测到的候选体细胞变异
包括可能的假阳性变异
需要后续过滤步骤



2. F1R2数据文件

路径: /groups/g5840141/home/zengqianwen/WES/mutect2/vcf/f1r2/
文件格式: ${肿瘤样本名}.f1r2.tar.gz
内容:

链特异性错误信号数据
用于后续的 LearnReadOrientationModel 分析
帮助识别和过滤测序伪影
      
"""
#!/bin/bash
#SBATCH -J 09_Mutect2
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -o 09_Mutect2.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e 09_Mutect2.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=2    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 6<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作
# 在fd6中放置$Nproc个空行作为令牌
for((i=1; i<=$Nproc; i++)); do
	echo
done >&6


input_folder_path=/groups/g5840141/home/zengqianwen/WES/align              # BAM 文件目录（已经对齐、去重、重校准）
output_folder_path=/groups/g5840141/home/zengqianwen/WES/mutect2/vcf/unfiltered # 得新建文件夹， .unfiltered.vcf 文件，保存所有未过滤的候选体细胞突变（原始变异）
output_folder_path_f1r2=/groups/g5840141/home/zengqianwen/WES/mutect2/vcf/f1r2  # 记得新建文件夹 .f1r2.tar.gz 文件，Mutect2 为 捕获链特异性错误信号生成的辅助数据，下一步 LearnReadOrientationModel 会使用这些 
reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/GRCh38.primary_assembly.genome.fa  # 人类参考基因组（GRCh38.p13，来自 GENCODE）
interval_path=/groups/g5840087/home/share/refGenome/reference/agilent/S07604514_hs_hg38/S07604514_Padded.bed   # 外显子靶向捕获区域的 BED 文件（这里是 Agilent S07604514 panel）
germline_resource_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz  # VCF 文件，包含了 人群中已知的常见生殖系变异，帮助 Mutect2 区分体细胞变异与生殖系变异
pon_path=/groups/g5840087/home/share/refGenome/reference/gatk/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz  # Panel of Normals（PON），从多个正常样本中构建的背景变异数据库，去除在多个正常样本中复现的假阳性变异（比如文库构建造成的错配），需要根据具体实验进行具体构建


# Sample_list=(PHCC2048T PHCC2048L RHCC2358T RHCC2358M PHCC3519T PHCC3519P RHCC4057Th RHCC4057L PHCC1889T PHCC1889L RHCC4055T RHCC4055P PHCC2992T PHCC2992L RHCC4050T RHCC4050P PHCC2071T PHCC2071P RHCC4173T1 RHCC4173T2 RHCC4173L)
Sample_list=(PHCC2417T1 PHCC2417T2 PHCC2417T3 PHCC2417T4 PHCC2417P RHCC4584T PHCC1011T PHCC1011L RHCC4619T RHCC4619P PHCC966T RHCC4664T RHCC4664P PHCC3603T RHCC4691T RHCC4691P PHCC972T PHCC972L RHCC4349T PHCC3369T PHCC3369P RHCC4121T RHCC4121L)

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"

while IFS=$'\t' read -r Tumor_sample_name Normal_sample_name; do
    tumor_file=${input_folder_path}/${Tumor_sample_name}.aligned.duplicates_marked.recalibrated.bam
    matched_normal_file=${input_folder_path}/${Normal_sample_name}.aligned.duplicates_marked.recalibrated.bam
    if [ -e ${output_folder_path}/${Tumor_sample_name}.unfiltered.vcf ]; then
                echo "${Tumor_sample_name}.unfiltered.vcf exists"
        else
                echo "${Tumor_sample_name}.unfiltered.vcf not exists"  
	        read -u6        # 领取令牌
	                {
	                    gatk Mutect2 \
                            --native-pair-hmm-threads 16 \
	                        -R ${reference_path} \
	                        --germline-resource ${germline_resource_path} \
	                        -pon ${pon_path} \
	                        -L ${interval_path} \
	                        -I ${tumor_file} \
	                        -I ${matched_normal_file} \
	                        -normal ${Normal_sample_name} \
	                        --f1r2-tar-gz ${output_folder_path_f1r2}/${Tumor_sample_name}.f1r2.tar.gz \
	                        -O ${output_folder_path}/${Tumor_sample_name}.unfiltered.vcf
	                        sleep 5
	                        echo >&6
	                }&
	fi
done < 08_Calculate_Contamination_Sample_list.txt
wait
exec 6>&-

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "09_Mutect2 done!"

