"""
完成比对，输出排序好的比对bam文件和bai索引
samtools 从版本 1.10 起默认生成的索引格式是 CSI（Coordinate Sorted Index），而非传统的 BAI。因此注意指定输出格式

| 特点                   | BAI索引            | CSI索引         
| --------------------- | ------------------ | ------------------------- 
| 最大单条染色体长度支持  | \~2GB              | 超过2GB，适合超大基因组 
| 兼容性                 | 更广，老软件支持    | 新软件支持更好       
| 结构复杂度             | 较简单              | 复杂，功能更强       

"""
#!/bin/bash
#SBATCH -J 04_Bwa_mem2
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -o 04_Bwa_mem2.o # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e 04_Bwa_mem2.e # 把报错结果STDERR保存在哪一个文件

source /opt/app/anaconda3/bin/activate
conda activate /home/zengqianwen/.conda/envs/gatk4-zqw/

#### 并行进程数量控制 ####
Nproc=8    # 可同时运行的最大作业数
Pfifo="/tmp/$$.fifo"   # 以PID为名, 防止创建命名管道时与已有文件重名，从而失败
mkfifo $Pfifo          # 创建命名管道
exec 6<>$Pfifo         # 以读写方式打开命名管道, 文件标识符fd为6
                       # fd可取除0, 1, 2,5外0-9中的任意数字
rm -f $Pfifo           # 删除文件, 也可不删除, 不影响后面操作
# 在fd6中放置$Nproc个空行作为令牌
for((i=1; i<=$Nproc; i++)); do
        echo
done >&6

input_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/trim
output_folder_path=/groups/g5840141/home/zengqianwen/WES_2025/align
bwa_fa_reference_path=/groups/g5840087/home/share/refGenome/reference/gencode/release_40/GRCh38.p13/fa/bwa-mem2/GRCh38.primary_assembly.genome.fa
# Sample_list=(PHCC2048T PHCC2048L RHCC2358T RHCC2358M PHCC3519T PHCC3519P RHCC4057Th RHCC4057L PHCC1889T PHCC1889L RHCC4055T RHCC4055P PHCC2992T PHCC2992L RHCC4050T RHCC4050P PHCC2071T PHCC2071P RHCC4173T1 RHCC4173T2 RHCC4173L)
# Sample_list=(PHCC2417T1 PHCC2417T2 PHCC2417T3 PHCC2417T4 PHCC2417P RHCC4584T PHCC1011T PHCC1011L RHCC4619T RHCC4619P PHCC966T RHCC4664T RHCC4664P PHCC3603T RHCC4691T RHCC4691P PHCC972T PHCC972L RHCC4349T PHCC3369T PHCC3369P RHCC4121T RHCC4121L)
# Sample_list=(PHCC1314T1 PHCC1314T2 PHCC1314L RHCC4319T RHCC4319P PHCC1333T PHCC1333L RHCC4366T RHCC4366P PHCC3635T PHCC3635P RHCC4376T RHCC4376L)

current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "04_Bwa_mem2 start!"


for file_1 in ${input_folder_path}/*_1_val_1.fq.gz; do
        file_2=$(echo ${file_1} | sed s/_1_val_1.fq.gz/_2_val_2.fq.gz/)

        sample_name=$(basename ${file_1} _1_val_1.fq.gz)

        current_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Current Time: $current_time"

        if [ -e ${output_folder_path}/${sample_name}.aligned.sorted.bam ]; then
                echo "${sample_name}.aligned.sorted.bam exists"
        else
                echo "${sample_name}.aligned.sorted.bam not exists"

                lane_num=$(zless ${file_1} | head -n1 | cut -d ":" -f4)
                flow_id=$(zless ${file_1} | head -n1 | cut -d ":" -f3)
                run_id=$(zless ${file_1} | head -n1 | cut -d ":" -f2)
                instrument_id=$(zless ${file_1} | head -n1 | cut -d ":" -f1 | sed s/@//)
                platform_unit=${flow_id}.${run_id}.${lane_num}
                id=${sample_name}.${flow_id}.${lane_num}

                read -u6 # 领取令牌, 控制进程数量
                {
                        bwa-mem2 mem \ 
                        -t 8 \ # 使用8个线程加速比对
                        -K 100000000 \ # 每次批处理100M的碱基，提高效率
                        -Y \ # 对大于一定长度的软剪切使用硬剪切（影响比对结果的CIGAR）
                        -R "@RG\\tID:${id}\\tLB:WES\\tSM:${sample_name}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tPM:${instrument_id}" \ # 添加Read Group信息，方便后续样本追踪
                        ${bwa_fa_reference_path} \ # 参考基因组索引路径
                        ${file_1} ${file_2} \ # 输入的成对FASTQ文件（PE测序）
                        | \
                        samtools sort -@ 8 \ #对 BWA 输出的比对结果进行排序和索引，排序是是让BAM文件中的reads在染色体上的顺序和参考基因组的顺序一致，索引是加速BAM文件的访问，，方便后续分析这两步方便了后续分析
                        -m 1500M \ # 每个线程最多使用1500M内存进行排序
                        -o ${output_folder_path}/${sample_name}.aligned.sorted.bam

                        # 单独生成BAI索引
                        samtools index -b -@ 8 ${output_folder_path}/${sample_name}.aligned.sorted.bam
                        
                        sleep 5
                        echo >&6 # 归还令牌
                }&
        fi
done
wait
exec 6>&-


current_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Current Time: $current_time"
echo "04_Bwa_mem2 done!"
