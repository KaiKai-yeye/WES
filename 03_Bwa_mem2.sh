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
                        -t 8 \
                        -K 100000000 \
                        -Y \
                        -R "@RG\\tID:${id}\\tLB:WES\\tSM:${sample_name}\\tPL:ILLUMINA\\tPU:${platform_unit}\\tPM:${instrument_id}" \
                        ${bwa_fa_reference_path} \
                        ${file_1} ${file_2} \
                        | \
                        samtools sort -@ 8 \
                        -m 1500M \
                        --write-index \
                        -o ${output_folder_path}/${sample_name}.aligned.sorted.bam##idx##${output_folder_path}/${sample_name}.aligned.sorted.bam.bai
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
