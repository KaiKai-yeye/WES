################################################################################
#
# 免疫编辑分析流程 (R 脚本)
#
# 目的: 
# 1. (TCGA部分) 利用TCGA泛癌数据建立一个“背景突变模型”。
#    - 该模型量化了每一种突变频谱（192种）产生同义突变、
#      非同义突变和新抗原的背景概率。
# 2. (自己队列部分) 将此背景模型应用到“自己队列”（如早发肝癌）上。
#    - 为每个病人计算一个“免疫编辑得分 (R)”，
#      R = (观测新抗原率) / (预期新抗原率)。
#    - R < 1 表明存在新抗原耗竭，即免疫编辑的证据。
#
# 依赖工具 (外部): 
# - netMHCpan (需要在R脚本之外单独运行)
#
################################################################################


################################################################################
# 0. 环境设置与库加载
################################################################################

# 加载所需 R 包
library(stringr) # 用于字符串操作
library(dplyr)   # 用于数据框操作 (filter, mutate, %>% 等)
library(tidyr)   # 用于数据整理 (pivot_longer)
library(readxl)  # 用于读取 Excel 文件 (本脚本中未明确使用，但可能为依赖)
library(BSgenome) # 用于操作基因组序列
library(BSgenome.Hsapiens.UCSC.hg38) # 人类 hg38 参考基因组
library(Biostrings) # 用于处理生物序列 (FASTA)
library(seqinr)     # 用于读写 FASTA 文件
library(doParallel) # 用于并行计算
library(ggplot2)    # 用于绘图
library(ggpubr)     # 用于在 ggplot2 中添加统计检验

# 清理工作环境
gc() # 垃圾回收，释放内存
rm(list = ls()) # 清除所有变量


################################################################################
#
# 第 1 部分: TCGA 泛癌背景模型计算
# 
# 目标: 统计192种突变频谱各自的 Nbackground (非同义/同义突变比率) 
#       和 Bbackground (新抗原/非同义突变比率)。
#
################################################################################

# --- 1.1 加载 TCGA 突变数据 (SNV) 和 HLA 分型数据 ---

print("第1部分: 开始加载 TCGA 数据...")

# 读入 TCGA 泛癌体细胞突变数据
# 注意: 请确保文件路径正确
SNVdata <- read.csv("D:/课题组/免疫编辑/data/GDC_TCGA_pancancer_somaticmut.csv") 
SNVdata <- SNVdata[,-1] # 移除第一列 (通常是行号)
# 提取12位TCGA病人ID
SNVdata$patient <- str_sub(SNVdata$Sample_ID, 1, 12) 
length(unique(SNVdata$patient)) # 10157

# 读入 TCGA HLA 分型数据
hla_data <- read.table("D:/课题组/免疫编辑/data/OptiTypeCallsHLA_20171207.txt", header = T, sep = '\t')
# 提取12位TCGA病人ID（和SNVdata统一）
hla_data$patient <- str_sub(hla_data$aliquot_id, 1, 12) 
length(unique(hla_data$patient)) # 8912

# --- 1.2 筛选共有样本 ---

# 保留同时具有 SNV 数据和 HLA 数据的病人
### 取交集
tmp <- intersect(hla_data$patient, SNVdata$patient) # 8368
hla_data1 <- subset(hla_data, hla_data$patient %in% tmp)
TCGA_SNV_data1 <- subset(SNVdata, SNVdata$patient %in% tmp)

print(paste("TCGA 共有样本数:", length(tmp)))

# --- 1.3 加工 HLA 数据 ---

# 确保每个病人ID唯一
any(duplicated(hla_data1$patient))
hla_data2 <- dplyr::distinct(hla_data1, patient, .keep_all = TRUE)

# 将宽数据转换为长数据 (A1, A2, B1, B2, C1, C2 合并为一列)
HLA_data_long <- hla_data2 %>% 
  pivot_longer(cols = A1:C2, names_to = 'Type', values_to = 'Allele')

# 将 HLA 分型转换为 netMHCpan 接受的输入格式 (例如: HLA-A*02:01)
HLA_data_long$Allele <- gsub('[*]', '', HLA_data_long$Allele)
HLA_data_long$Allele <- paste0("HLA-", HLA_data_long$Allele)

# --- 1.4 加工错义突变数据 (用于新抗原预测) ---

print("加工 TCGA 错义突变数据...")

# 提取所有“非同义突变”(nonsynonymous SNV)
table(TCGA_SNV_data1$ExonicFunc.refGene)
TCGA_missense_mutation <- TCGA_SNV_data1[TCGA_SNV_data1$ExonicFunc.refGene %in% 'nonsynonymous SNV',]
length(unique(TCGA_missense_mutation$patient)) # 8359

# 提取蛋白突变信息 (例如 'A243V')
TCGA_missense_mutation$object <- sapply(strsplit(unlist(TCGA_missense_mutation$Protein_change), 'p.'), '[', 2)
# 为每个突变创建唯一ID
TCGA_missense_mutation$mutation_no <- paste0("no", 1:nrow(TCGA_missense_mutation))

# --- 1.5 读入并合并蛋白序列 (FASTA) ---

# 1. 读入自己 cohort 的转录本编号对应的蛋白序列 (可能包含TCGA缺失的部分)
fastaFile1 <- readBStringSet("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/protein_missence_mut_v2.fasta", format="fasta")
protein_file1 <- as.data.frame(fastaFile1)
protein_file1$name <- fastaFile1@ranges@NAMES
if(F){ # 原始代码中的调试/检查块，保留原样
  any(duplicated(protein_file1$name))
  protein_file1 <- protein_file1 %>%
    distinct(name, .keep_all = T)
}

# 2. 读入TCGA其余的蛋白序列
fastaFile2 <- readBStringSet("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA_protein_missence_mut.fasta", format="fasta")
protein_file2 <- as.data.frame(fastaFile2)
protein_file2$name <- fastaFile2@ranges@NAMES

# 3. 合并蛋白序列库
protein_file <- rbind(protein_file1, protein_file2)
colnames(protein_file) <- c("x", "NM_number") # 'x' 列是蛋白序列
# 按转录本ID (NM_number) 去重
protein_file <- protein_file %>% distinct(NM_number, .keep_all = TRUE)

# 将蛋白序列 (x) 左连接到错义突变数据框
TCGA_missense_mutation <- dplyr::left_join(TCGA_missense_mutation, protein_file, by = "NM_number")

# --- 1.6 定义肽提取函数 (8/9/10-mer) ---

# 提取以突变位点为中心的8-mer肽序列
extract_8mer <- function(sequence, mutation) {
  pre_mutation <- substr(mutation, 1, 1) # 突变前氨基酸
  position <- as.numeric(gsub("[A-Za-z]", "", mutation)) # 突变位置
  post_mutation <- substr(mutation, nchar(mutation), nchar(mutation)) # 突变后氨基酸
  seq_len <- nchar(sequence)
  
  # 提取前7个氨基酸 (注意边界)
  start_pos <- max(1, position - 7)
  pre_sequence <- substr(sequence, start_pos, position - 1)
  
  # 提取后7个氨基酸 (注意边界)
  end_pos <- min(seq_len, position + 7)
  post_sequence <- substr(sequence, position + 1, end_pos)
  
  # 构建包含突变氨基酸的序列 (用于netMHCpan输入)
  result <- paste0(pre_sequence, post_mutation, post_sequence)
  return(result)
}

# 提取以突变位点为中心的9-mer肽序列
extract_9mer <- function(sequence, mutation) {
  pre_mutation <- substr(mutation, 1, 1)
  position <- as.numeric(gsub("[A-Za-z]", "", mutation))
  post_mutation <- substr(mutation, nchar(mutation), nchar(mutation))
  seq_len <- nchar(sequence)
  
  # 提取前8个氨基酸
  start_pos <- max(1, position - 8)
  pre_sequence <- substr(sequence, start_pos, position - 1)
  
  # 提取后8个氨基酸
  end_pos <- min(seq_len, position + 8)
  post_sequence <- substr(sequence, position + 1, end_pos)
  
  result <- paste0(pre_sequence, post_mutation, post_sequence)
  return(result)
}

# 提取以突变位点为中心的10-mer肽序列
extract_10mer <- function(sequence, mutation) {
  pre_mutation <- substr(mutation, 1, 1)
  position <- as.numeric(gsub("[A-Za-z]", "", mutation))
  post_mutation <- substr(mutation, nchar(mutation), nchar(mutation))
  seq_len <- nchar(sequence)
  
  # 提取前9个氨基酸
  start_pos <- max(1, position - 9)
  pre_sequence <- substr(sequence, start_pos, position - 1)
  
  # 提取后9个氨基酸
  end_pos <- min(seq_len, position + 9)
  post_sequence <- substr(sequence, position + 1, end_pos)
  
  result <- paste0(pre_sequence, post_mutation, post_sequence)
  return(result)
}

# --- 1.7 生成 netMHCpan 输入文件 (肽 FASTA) ---

print("生成 TCGA netMHCpan 肽 FASTA 输入文件...")
# 按病人遍历，为每个病人生成 8, 9, 10-mer 的 FASTA 文件
for (i in unique(TCGA_missense_mutation$patient)) {
  
  # 筛选该病人的突变，并计算 8/9/10-mer 序列
  tmp_df <- TCGA_missense_mutation %>% 
    dplyr::filter(patient == i) %>% 
    dplyr::mutate(result8 = mapply(extract_8mer, x, object),
                  result9 = mapply(extract_9mer, x, object),
                  result10 = mapply(extract_10mer, x, object))
  
  # 将序列写入 FASTA 文件
  for (n in 1:nrow(tmp_df)) {
    write.fasta(sequences = tmp_df[n,'result8'], names = paste0(tmp_df[n,'mutation_no'],"_",tmp_df[n,'patient']), open = 'a',
                file.out = paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_input_aa8/',i,'.fsa'))
    
    write.fasta(sequences = tmp_df[n,'result9'], names = paste0(tmp_df[n,'mutation_no'],"_",tmp_df[n,'patient']), open = 'a',
                file.out = paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_input_aa9/',i,'.fsa'))
    
    write.fasta(sequences = tmp_df[n,'result10'], names = paste0(tmp_df[n,'mutation_no'],"_",tmp_df[n,'patient']), open = 'a',
                file.out = paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_input_aa10/',i,'.fsa'))
  }
}

# --- 1.8 生成 netMHCpan 输入文件 (HLA) ---

print("生成 TCGA netMHCpan HLA 输入文件...")
# 筛选 HLA 数据，使其与错义突变数据的病人一致
HLA_data_long1 <- subset(HLA_data_long, HLA_data_long$patient %in% unique(TCGA_missense_mutation$patient))

# 为每个病人生成一个 HLA 文件，包含该病人所有的 HLA 等位基因 (逗号分隔)
for (i in HLA_data_long1$patient){
  test1 <- subset(HLA_data_long1, HLA_data_long1$patient == i)
  test_hla <- paste(test1$Allele, collapse = ",")
  write.table(test_hla, 
              file = paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/netMHCpan_input_hla/",i,"_hla.txt"),
              quote = F, col.names = F, row.names = F)
}


####################################################################
#
# !!! 外部工具运行步骤 !!!
#
# 此处需要暂停 R 脚本，转到 Linux/Shell 环境，
# 使用 netMHCpan 工具处理上一步生成的 FASTA 和 HLA 文件。
#
# (伪代码示例):
# for patient in $(ls .../netMHCpan_input_aa8/); do
#   netMHCpan -p .../netMHCpan_input_aa8/${patient} \
#             -a $(cat .../netMHCpan_input_hla/${patient}_hla.txt) \
#             -l 8 > .../netMHCpan_out_aa8/${patient}.xls
#   ... (对 9-mer 和 10-mer 重复此操作)
# done
#
####################################################################


# --- 1.9 解析 netMHCpan 输出结果 ---

print("解析 TCGA netMHCpan 输出结果...")

# 1. 找到 8, 9, 10-mer 均有输出结果的病人
tmp1 <- list.files(path = './Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_out_aa8', pattern = '*.xls')
ccname1 <- gsub('.xls', '', tmp1)
tmp2 <- list.files(path = './Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_out_aa9/', pattern = '*.xls')
ccname2 <- gsub('.xls', '', tmp2)
tmp3 <- list.files(path = './Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_out_aa10/', pattern = '*.xls')
ccname3 <- gsub('.xls', '', tmp3)

# 共同的病人列表
common_ccname <- Reduce(intersect, list(ccname1, ccname2, ccname3)) # 8342
print(paste("TCGA netMHCpan 完整结果样本数:", length(common_ccname)))

# 2. 初始化一个数据框，用于存储每个突变是否为新抗原 (1/0)
final_df <- data.frame(mut_no = character(), antigen_final = integer())

# 按病人遍历 netMHCpan 输出文件
for (patient in common_ccname) {
  
  # --- 处理 8-mer 结果 ---
  # 读入并清理 netMHCpan 的 .xls 输出 (它是制表符分隔的)
  tmp_file1 <- read.table(paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_out_aa8/", patient, ".xls"), header = F, fill = T, sep = '\t', col.names = paste0('V', seq_len(35)))
  tmp_file1 <- tmp_file1[-1,] # 移除第一行空行
  colnames(tmp_file1) <- tmp_file1[1,] # 设置表头
  tmp_file1 <- tmp_file1[-1,] # 移除表头行
  tmp_file1 <- tmp_file1[, str_detect(colnames(tmp_file1), 'ID|nM')] # 仅保留 ID 和亲和力 (nM) 列
  tmp_file1$mut_no <- str_split(tmp_file1$ID, "[_]", simplify = T)[, 1] # 提取突变ID
  tmp_file1 <- tmp_file1[, -1] # 移除原始ID列
  tmp_file1[, 1:(ncol(tmp_file1)-1)] <- lapply(tmp_file1[, 1:(ncol(tmp_file1)-1)], as.numeric) # 转换为数值
  
  # 检查每个突变ID (mut_no) 是否产生了新抗原
  unique_mut_no <- unique(tmp_file1$mut_no)
  result_df1 <- data.frame(mut_no = character(), antigen = integer(), stringsAsFactors = FALSE)
  for (mut in unique_mut_no) {
    subset_df <- tmp_file1[tmp_file1$mut_no == mut, ]
    # ** 关键定义: 如果任何一个肽-HLA组合的亲和力 < 500nM，则定义为新抗原 **
    if (any(subset_df[, 1:(ncol(tmp_file1)-1)] < 500, na.rm = TRUE)) {
      antigen_value <- 1 # 是新抗原
    } else {
      antigen_value <- 0 # 不是新抗原
    }
    result_df1 <- rbind(result_df1, data.frame(mut_no = mut, antigen = antigen_value))
  }
  
  # --- 处理 9-mer 结果 (逻辑同上) ---
  tmp_file2 <- read.table(paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_out_aa9/", patient, ".xls"), header = F, fill = T, sep = '\t', col.names = paste0('V', seq_len(35)))
  tmp_file2 <- tmp_file2[-1,]; colnames(tmp_file2) <- tmp_file2[1,]; tmp_file2 <- tmp_file2[-1,]
  tmp_file2 <- tmp_file2[, str_detect(colnames(tmp_file2), 'ID|nM')]
  tmp_file2$mut_no <- str_split(tmp_file2$ID, "[_]", simplify = T)[, 1]
  tmp_file2 <- tmp_file2[, -1]
  tmp_file2[, 1:(ncol(tmp_file2)-1)] <- lapply(tmp_file2[, 1:(ncol(tmp_file2)-1)], as.numeric)
  
  unique_mut_no <- unique(tmp_file2$mut_no)
  result_df2 <- data.frame(mut_no = character(), antigen = integer(), stringsAsFactors = FALSE)
  for (mut in unique_mut_no) {
    subset_df <- tmp_file2[tmp_file2$mut_no == mut, ]
    if (any(subset_df[, 1:(ncol(tmp_file2)-1)] < 500, na.rm = TRUE)) {
      antigen_value <- 1
    } else {
      antigen_value <- 0
    }
    result_df2 <- rbind(result_df2, data.frame(mut_no = mut, antigen = antigen_value))
  }
  
  # --- 处理 10-mer 结果 (逻辑同上) ---
  tmp_file3 <- read.table(paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/GDC_TCGA/netMHCpan_out_aa10/", patient, ".xls"), header = F, fill = T, sep = '\t', col.names = paste0('V', seq_len(35)))
  tmp_file3 <- tmp_file3[-1,]; colnames(tmp_file3) <- tmp_file3[1,]; tmp_file3 <- tmp_file3[-1,]
  tmp_file3 <- tmp_file3[, str_detect(colnames(tmp_file3), 'ID|nM')]
  tmp_file3$mut_no <- str_split(tmp_file3$ID, "[_]", simplify = T)[, 1]
  tmp_file3 <- tmp_file3[, -1]
  tmp_file3[, 1:(ncol(tmp_file3)-1)] <- lapply(tmp_file3[, 1:(ncol(tmp_file3)-1)], as.numeric)
  
  unique_mut_no <- unique(tmp_file3$mut_no)
  result_df3 <- data.frame(mut_no = character(), antigen = integer(), stringsAsFactors = FALSE)
  for (mut in unique_mut_no) {
    subset_df <- tmp_file3[tmp_file3$mut_no == mut, ]
    if (any(subset_df[, 1:(ncol(tmp_file3)-1)] < 500, na.rm = TRUE)) {
      antigen_value <- 1
    } else {
      antigen_value <- 0
    }
    result_df3 <- rbind(result_df3, data.frame(mut_no = mut, antigen = antigen_value))
  }
  
  # --- 合并 8, 9, 10-mer 结果 ---
  # 只要 8, 9, 10-mer 中任何一个被预测为新抗原，该突变 (mut_no) 就被视为新抗原
  
  union_mut_no <- Reduce(union, list(result_df1$mut_no, result_df2$mut_no, result_df3$mut_no))
  
  df1_union <- merge(data.frame(mut_no = union_mut_no), result_df1, by = "mut_no", all.x = TRUE)
  df2_union <- merge(df1_union, result_df2, by = "mut_no", all.x = TRUE)
  df3_union <- merge(df2_union, result_df3, by = "mut_no", all.x = TRUE)
  
  # 替换 NA (未在 8/9/10-mer 中检出的) 为 0 (非新抗原)
  df3_union$antigen.x[is.na(df3_union$antigen.x)] <- 0
  df3_union$antigen.y[is.na(df3_union$antigen.y)] <- 0
  df3_union$antigen[is.na(df3_union$antigen)] <- 0
  
  # 只要有一个为1，总和就>0
  df3_union$antigen_final <- ifelse(rowSums(df3_union[, c("antigen.x", "antigen.y", "antigen")]) > 0, 1, 0)
  
  # 将该病人的结果汇总到 final_df
  final_df <- rbind(final_df, df3_union[, c("mut_no", "antigen_final")])
}
colnames(final_df) <- c("mutation_no", "antigen")
# final_df <- final_df[-1,] # 移除初始化行 (如果使用初始化行)
# 确保 final_df 是唯一的
final_df <- final_df %>% distinct(mutation_no, .keep_all = TRUE) 

print("TCGA 新抗原预测结果解析完成。")

# --- 1.10 制作 TCGA 突变频谱 (192种) ---

# 定义一个辅助函数，用于获取突变频谱 (三核苷酸上下文)
get_spectra <- function(chr, pos, ref, alt, genome) {
  # 确保染色体名称在基因组对象中 (例如 "chr1" vs "1")
  chr_name <- chr
  if (!chr_name %in% seqnames(genome)) {
    chr_name <- paste0("chr", chr)
    if (!chr_name %in% seqnames(genome)) {
      return(NA) # 无法识别的染色体
    }
  }
  
  seq <- genome[[chr_name]]
  
  # 检查边界
  if (pos > 1 && pos < length(seq)) {
    # 提取前中后三个碱基
    triplet <- as.character(seq[(pos - 1):(pos + 1)])
    # 格式化为 A[C>T]G 形式
    spect <- paste0(triplet[1], "[", ref, ">", alt, "]", triplet[3])
    return(spect)
  } else {
    return(NA) # 边界情况
  }
}

# --- 1.10a 制作非同义突变 (Missense) 的频谱 ---
print("制作 TCGA 非同义突变频谱...")

# (为节约时间，原始代码中此处为从文件读取，这里注释掉循环，假设文件已生成)
# 这是一个非常耗时的循环，如果 TCGA_missense_mutation 很大
# TCGA_new <- data.frame(spectra = character())
# for (i in 1:nrow(TCGA_missense_mutation)) {
#   spect <- get_spectra(TCGA_missense_mutation$Chr[i],
#                        as.numeric(TCGA_missense_mutation$Start[i]),
#                        TCGA_missense_mutation$Ref[i],
#                        TCGA_missense_mutation$Alt[i],
#                        Hsapiens)
#   TCGA_new <- rbind(TCGA_new, data.frame(spectra = spect))
# }
# TCGA_missense_mutation <- cbind(TCGA_missense_mutation, TCGA_new)
# TCGA_missense_mutation <- na.omit(TCGA_missense_mutation) # 移除NA

# 直接读入已处理好的文件 (如原始代码所示)
TCGA_missense_mutation <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/TCGA_missense_mutation_addspectra.csv")


# --- 1.10b 制作同义突变 (Synonymous) 的频谱 ---
print("制作 TCGA 同义突变频谱...")

# 提取同义突变
TCGA_synonymous_data <- TCGA_SNV_data1[TCGA_SNV_data1$ExonicFunc.refGene == "synonymous SNV", c("Chr", "Start", "End", "Sample_ID", "Ref", "Alt", "patient")]
length(unique(TCGA_synonymous_data$patient)) # 8219

# (同上，这是一个耗时的循环，假设文件已生成)
# TCGA_new <- data.frame(spectra = character())
# for (i in 1:nrow(TCGA_synonymous_data)) {
#   spect <- get_spectra(TCGA_synonymous_data$Chr[i],
#                        as.numeric(TCGA_synonymous_data$Start[i]),
#                        TCGA_synonymous_data$Ref[i],
#                        TCGA_synonymous_data$Alt[i],
#                        Hsapiens)
#   TCGA_new <- rbind(TCGA_new, data.frame(spectra = spect))
# }
# TCGA_synonymous_data <- cbind(TCGA_synonymous_data, TCGA_new)
# TCGA_synonymous_data <- na.omit(TCGA_synonymous_data)

# 直接读入已处理好的文件 (如原始代码所示)
TCGA_synonymous_data <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/TCGA_synonymous_data_addspectra.csv")


# --- 1.11 统计每种频谱的突变总数 (建立背景模型) ---

print("开始统计 TCGA 背景模型概率...")

# 读入192种标准频谱
spectra <- read.table('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/spectra192.txt', header = FALSE, sep = '\t', col.names = "spectrum")

# --- 1.11a 统计同义突变 (syn_num) ---

# 筛选出有 netMHCpan 结果的病人的同义突变
TCGA_synonymous_data_8348 <- TCGA_synonymous_data[TCGA_synonymous_data$patient %in% common_ccname,]
length(unique(TCGA_synonymous_data_8348$patient)) ## 8193

# 按192种频谱计数
syn_count <- data.frame(spectrum = character(), syn_num = integer())
for (i in 1:nrow(spectra)){
  spec_value <- spectra[i, 1]
  spec1 <- subset(TCGA_synonymous_data_8348, TCGA_synonymous_data_8348$spectra == spec_value)
  syn1 <- nrow(spec1)
  a1 <- data.frame(spectrum = spec_value, syn_num = syn1)
  syn_count <- rbind(syn_count, a1)
}

# --- 1.11b 统计非同义突变 (nonsyn_num) ---

# 筛选出有 netMHCpan 结果的病人的非同义突变
TCGA_missense_mutation_8348 <- TCGA_missense_mutation[TCGA_missense_mutation$patient %in% common_ccname,]
length(unique(TCGA_missense_mutation_8348$patient)) # 8342

# 合并新抗原预测结果 (final_df)
TCGA_missense_mutation_8348 <- left_join(TCGA_missense_mutation_8348, final_df, by = "mutation_no")

# (原始代码中的调试块)
if(F){ 
  any(is.na(TCGA_missense_mutation_8348$antigen))
  a <- setdiff(TCGA_missense_mutation_8348$mutation_no, final_df$mutation_no)
  tmp <- TCGA_missense_mutation_8348[TCGA_missense_mutation_8348$mutation_no %in% a,]
  TCGA_missense_mutation_8348$antigen[is.na(TCGA_missense_mutation_8348$antigen)] <- 0
}

# 确保用于统计背景模型的两组病人完全一致
# (使用有同义突变数据的病人列表来筛选非同义突变数据)
common_patients_final <- unique(TCGA_synonymous_data_8348$patient)
TCGA_missense_mutation_8348_new <- TCGA_missense_mutation_8348[TCGA_missense_mutation_8348$patient %in% common_patients_final,]
length(unique(TCGA_missense_mutation_8348_new$patient)) # 确保病人集一致

# 按192种频谱计数
nonsyn_count <- data.frame(spectrum = character(), nonsyn_num = integer())
for (i in 1:nrow(spectra)){
  spec_value <- spectra[i, 1]
  spec1 <- subset(TCGA_missense_mutation_8348_new, TCGA_missense_mutation_8348_new$spectra == spec_value)
  syn1 <- nrow(spec1)
  a1 <- data.frame(spectrum = spec_value, nonsyn_num = syn1)
  nonsyn_count <- rbind(nonsyn_count, a1)
}

# --- 1.11c 统计新抗原 (antigen) ---

# 使用与上一步相同的、病人一致的非同义突变数据
antigen_count <- data.frame(spectrum = character(), antigen = integer())
for (i in 1:nrow(spectra)){
  spec_value <- spectra[i, 1]
  spec1 <- subset(TCGA_missense_mutation_8348_new, TCGA_missense_mutation_8348_new$spectra == spec_value)
  syn1 <- sum(spec1$antigen, na.rm = TRUE) # 计算新抗原总数 (1/0)
  a1 <- data.frame(spectrum = spec_value, antigen = syn1)
  antigen_count <- rbind(antigen_count, a1)
}

# --- 1.12 计算最终背景概率 (方法1 和 方法2) ---

print("计算 TCGA 背景模型...")

# --- 方法1 (Rooney et al. 2015 原始方法) ---
# Nbackground = 非同义 / 同义
# Bbackground = 新抗原 / 非同义
final_count <- left_join(antigen_count, nonsyn_count, by = "spectrum")
final_count <- left_join(final_count, syn_count, by = "spectrum")
final_count$Nbackground <- final_count$nonsyn_num / final_count$syn_num
final_count$Bbackground <- final_count$antigen / final_count$nonsyn_num
# 处理除以0的情况 (如果 syn_num 或 nonsyn_num 为0)
final_count$Nbackground[is.infinite(final_count$Nbackground)] <- 0
final_count$Bbackground[is.na(final_count$Bbackground)] <- 0

write.csv(final_count, "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/TCGA_final_count.csv", row.names = F)

# --- 方法2 (备用方法) ---
# Nbackground = 非同义 / (非同义 + 同义)
# Bbackground = 新抗原 / (非同义 + 同义)
final_count1 <- left_join(antigen_count, nonsyn_count, by = "spectrum")
final_count1 <- left_join(final_count1, syn_count, by = "spectrum")
final_count1$Nbackground <- final_count1$nonsyn_num / (final_count1$nonsyn_num + final_count1$syn_num)
final_count1$Bbackground <- final_count1$antigen / (final_count1$nonsyn_num + final_count1$syn_num)
# 处理除以0的情况
final_count1$Nbackground[is.na(final_count1$Nbackground)] <- 0
final_count1$Bbackground[is.na(final_count1$Bbackground)] <- 0

write.csv(final_count1, "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/TCGA_final_count1.csv", row.names = F)

print("第1部分: TCGA 背景模型构建完成。")


################################################################################
#
# 第 2 部分: 自己队列 (早发肝癌) 样本计算
# 
# 目标: 使用第1部分构建的 'final_count' (背景模型)，
#       为本队列的每个样本计算免疫编辑得分 (R)。
#
################################################################################

print("第2部分: 开始处理“自己队列”数据...")

# --- 2.1 加载并处理“自己队列”突变数据 ---

my_mutation <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/updated_mutation.csv")
my_mutation <- my_mutation[, c(1:29)]
# 解析 ANNOVAR 格式的 AAChange.refGene 列
my_mutation <- my_mutation %>%
  mutate(
    gene_info = str_split(AAChange.refGene, ",", simplify = TRUE)[, 1],
    Gene = str_extract(gene_info, "^[^:]+"),
    NM_number = str_extract(gene_info, "(?<=:)NM_\\d+"),
    Protein_change = str_extract(gene_info, "p\\.[A-Za-z]+\\d+[A-Za-z]+")
  )
# 筛选出同义和非同义 SNV
my_SNVdata <- my_mutation[my_mutation$ExonicFunc.refGene %in% c('nonsynonymous SNV', 'synonymous SNV'),]

# --- 2.2 加工“自己队列”错义突变数据 ---

my_missense_mutation <- my_SNVdata[my_SNVdata$ExonicFunc.refGene == "nonsynonymous SNV", 
                                   c("Chromosome", "i.Start_Position", "i.End_Position", 
                                     "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", 
                                     "Func.refGene", "Gene", "NM_number", "Protein_change")]
my_missense_mutation$object <- sapply(strsplit(unlist(my_missense_mutation$Protein_change), 'p.'), '[', 2)
my_missense_mutation$mutation_no <- paste0("no", 1:nrow(my_missense_mutation)) # 创建唯一突变ID

# 合并蛋白序列 (使用之前读入的 protein_file1)
colnames(protein_file1) <- c("x", "NM_number")
protein_file1 <- protein_file1 %>% distinct(NM_number, .keep_all = TRUE)
my_missense_mutation <- dplyr::left_join(my_missense_mutation, protein_file1, by = "NM_number")

# --- 2.3 生成 netMHCpan 输入文件 (肽 FASTA) - 自己队列 ---

print("生成“自己队列” netMHCpan 肽 FASTA 输入文件...")
# 逻辑同 1.7
for (i in unique(my_missense_mutation$Tumor_Sample_Barcode)) {
  
  tmp_df <- my_missense_mutation %>% 
    dplyr::filter(Tumor_Sample_Barcode == i) %>% 
    dplyr::mutate(result8 = mapply(extract_8mer, x, object),
                  result9 = mapply(extract_9mer, x, object),
                  result10 = mapply(extract_10mer, x, object))
  
  for (n in 1:nrow(tmp_df)) {
    write.fasta(sequences = tmp_df[n,'result8'], names = paste0(tmp_df[n,'mutation_no'],"_",tmp_df[n,'Tumor_Sample_Barcode']), open = 'a',
                file.out = paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_input_aa8/',i,'.fsa'))
    
    write.fasta(sequences = tmp_df[n,'result9'], names = paste0(tmp_df[n,'mutation_no'],"_",tmp_df[n,'Tumor_Sample_Barcode']), open = 'a',
                file.out = paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_input_aa9/',i,'.fsa'))
    
    write.fasta(sequences = tmp_df[n,'result10'], names = paste0(tmp_df[n,'mutation_no'],"_",tmp_df[n,'Tumor_Sample_Barcode']), open = 'a',
                file.out = paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_input_aa10/',i,'.fsa'))
  }
}

# --- 2.4 HLA 输入文件 - 自己队列 ---
# (注释: 假设“自己队列”的 HLA 输入文件已按 1.8 逻辑准备好)
print("“自己队列” HLA 输入文件已跳过 (假设已准备好)。")

####################################################################
#
# !!! 外部工具运行步骤 !!!
#
# 再次需要暂停 R 脚本，为“自己队列”运行 netMHCpan。
#
####################################################################

# --- 2.5 解析 netMHCpan 输出结果 - 自己队列 ---

print("解析“自己队列” netMHCpan 输出结果...")

# 1. 找到 8, 9, 10-mer 均有输出结果的病人
tmp1 <- list.files(path = './Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_out_aa8_new/', pattern = '*.xls')
ccname1 <- gsub('.xls', '', tmp1)
tmp2 <- list.files(path = './Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_out_aa9_new/', pattern = '*.xls')
ccname2 <- gsub('.xls', '', tmp2)
tmp3 <- list.files(path = './Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_out_aa10_new/', pattern = '*.xls')
ccname3 <- gsub('.xls', '', tmp3)

my_common_ccname <- Reduce(intersect, list(ccname1, ccname2, ccname3))
print(paste("“自己队列” netMHCpan 完整结果样本数:", length(my_common_ccname)))

# 2. 初始化数据框
my_final_df <- data.frame(mut_no = character(), antigen_final = integer())

# 遍历病人，逻辑与 1.9 完全相同
for (patient in my_common_ccname) {
  
  # --- 8-mer (mydata) ---
  tmp_file1 <- read.table(paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_out_aa8_new/", patient, ".xls"), header = F, fill = T, sep = '\t', col.names = paste0('V', seq_len(35)))
  tmp_file1 <- tmp_file1[-1,]; colnames(tmp_file1) <- tmp_file1[1,]; tmp_file1 <- tmp_file1[-1,]
  tmp_file1 <- tmp_file1[, str_detect(colnames(tmp_file1), 'ID|nM')]
  tmp_file1$mut_no <- str_split(tmp_file1$ID, "[_]", simplify = T)[, 1]
  tmp_file1 <- tmp_file1[, -1]
  tmp_file1[, 1:(ncol(tmp_file1)-1)] <- lapply(tmp_file1[, 1:(ncol(tmp_file1)-1)], as.numeric)
  
  unique_mut_no <- unique(tmp_file1$mut_no)
  result_df1 <- data.frame(mut_no = character(), antigen = integer(), stringsAsFactors = FALSE)
  for (mut in unique_mut_no) {
    subset_df <- tmp_file1[tmp_file1$mut_no == mut, ]
    if (any(subset_df[, 1:(ncol(tmp_file1)-1)] < 500, na.rm = TRUE)) {
      antigen_value <- 1
    } else {
      antigen_value <- 0
    }
    result_df1 <- rbind(result_df1, data.frame(mut_no = mut, antigen = antigen_value))
  }
  
  # --- 9-mer (mydata) ---
  tmp_file2 <- read.table(paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_out_aa9_new/", patient, ".xls"), header = F, fill = T, sep = '\t', col.names = paste0('V', seq_len(35)))
  tmp_file2 <- tmp_file2[-1,]; colnames(tmp_file2) <- tmp_file2[1,]; tmp_file2 <- tmp_file2[-1,]
  tmp_file2 <- tmp_file2[, str_detect(colnames(tmp_file2), 'ID|nM')]
  tmp_file2$mut_no <- str_split(tmp_file2$ID, "[_]", simplify = T)[, 1]
  tmp_file2 <- tmp_file2[, -1]
  tmp_file2[, 1:(ncol(tmp_file2)-1)] <- lapply(tmp_file2[, 1:(ncol(tmp_file2)-1)], as.numeric)
  
  unique_mut_no <- unique(tmp_file2$mut_no)
  result_df2 <- data.frame(mut_no = character(), antigen = integer(), stringsAsFactors = FALSE)
  for (mut in unique_mut_no) {
    subset_df <- tmp_file2[tmp_file2$mut_no == mut, ]
    if (any(subset_df[, 1:(ncol(tmp_file2)-1)] < 500, na.rm = TRUE)) {
      antigen_value <- 1
    } else {
      antigen_value <- 0
    }
    result_df2 <- rbind(result_df2, data.frame(mut_no = mut, antigen = antigen_value))
  }
  
  # --- 10-mer (mydata) ---
  tmp_file3 <- read.table(paste0("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/mydata/netMHCpan_out_aa10_new/", patient, ".xls"), header = F, fill = T, sep = '\t', col.names = paste0('V', seq_len(35)))
  tmp_file3 <- tmp_file3[-1,]; colnames(tmp_file3) <- tmp_file3[1,]; tmp_file3 <- tmp_file3[-1,]
  tmp_file3 <- tmp_file3[, str_detect(colnames(tmp_file3), 'ID|nM')]
  tmp_file3$mut_no <- str_split(tmp_file3$ID, "[_]", simplify = T)[, 1]
  tmp_file3 <- tmp_file3[, -1]
  tmp_file3[, 1:(ncol(tmp_file3)-1)] <- lapply(tmp_file3[, 1:(ncol(tmp_file3)-1)], as.numeric)
  
  unique_mut_no <- unique(tmp_file3$mut_no)
  result_df3 <- data.frame(mut_no = character(), antigen = integer(), stringsAsFactors = FALSE)
  for (mut in unique_mut_no) {
    subset_df <- tmp_file3[tmp_file3$mut_no == mut, ]
    if (any(subset_df[, 1:(ncol(tmp_file3)-1)] < 500, na.rm = TRUE)) {
      antigen_value <- 1
    } else {
      antigen_value <- 0
    }
    result_df3 <- rbind(result_df3, data.frame(mut_no = mut, antigen = antigen_value))
  }
  
  # --- 合并 (mydata) ---
  union_mut_no <- Reduce(union, list(result_df1$mut_no, result_df2$mut_no, result_df3$mut_no))
  df1_union <- merge(data.frame(mut_no = union_mut_no), result_df1, by = "mut_no", all.x = TRUE)
  df2_union <- merge(df1_union, result_df2, by = "mut_no", all.x = TRUE)
  df3_union <- merge(df2_union, result_df3, by = "mut_no", all.x = TRUE)
  
  df3_union$antigen.x[is.na(df3_union$antigen.x)] <- 0
  df3_union$antigen.y[is.na(df3_union$antigen.y)] <- 0
  df3_union$antigen[is.na(df3_union$antigen)] <- 0
  
  df3_union$antigen_final <- ifelse(rowSums(df3_union[, c("antigen.x", "antigen.y", "antigen")]) > 0, 1, 0)
  
  my_final_df <- rbind(my_final_df, df3_union[, c("mut_no", "antigen_final")])
}
colnames(my_final_df) <- c("mutation_no", "antigen")
my_final_df <- my_final_df %>% distinct(mutation_no, .keep_all = TRUE) # 确保唯一

print("“自己队列”新抗原预测结果解析完成。")

# --- 2.6 制作“自己队列”突变频谱 (192种) ---

# --- 2.6a 同义突变 (mydata) ---
print("制作“自己队列”同义突变频谱...")
my_synonymous_data <- my_SNVdata[my_SNVdata$ExonicFunc.refGene == "synonymous SNV", 
                                 c("Chromosome", "i.Start_Position", "i.End_Position", 
                                   "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2")]
# (耗时循环，同 1.10)
# new <- data.frame(spectra = character())
# for (i in 1:nrow(my_synonymous_data)) {
#   spect <- get_spectra(my_synonymous_data$Chromosome[i],
#                        as.numeric(my_synonymous_data$i.Start_Position[i]),
#                        my_synonymous_data$Reference_Allele[i],
#                        my_synonymous_data$Tumor_Seq_Allele2[i],
#                        Hsapiens)
#   new <- rbind(new, data.frame(spectra = spect))
# }
# my_synonymous_data <- cbind(my_synonymous_data, new)
# my_synonymous_data <- na.omit(my_synonymous_data)
# write.csv(my_synonymous_data, "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/my_synonymous_data.csv", row.names = F)

# 读入已处理好的文件
my_synonymous_data <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/my_synonymous_data.csv")

# 筛选有 netMHCpan 结果的病人
my_synonymous_data_413 <- my_synonymous_data[my_synonymous_data$Tumor_Sample_Barcode %in% my_common_ccname,]
length(unique(my_synonymous_data_413$Tumor_Sample_Barcode)) # 413

# --- 2.6b 非同义突变 (mydata) ---
print("制作“自己队列”非同义突变频谱...")
# (耗时循环，同 1.10)
# new <- data.frame(spectra = character())
# for (i in 1:nrow(my_missense_mutation)) {
#   spect <- get_spectra(my_missense_mutation$Chromosome[i],
#                        as.numeric(my_missense_mutation$i.Start_Position[i]),
#                        my_missense_mutation$Reference_Allele[i],
#                        my_missense_mutation$Tumor_Seq_Allele2[i],
#                        Hsapiens)
#   new <- rbind(new, data.frame(spectra = spect))
# }
# my_missense_mutation <- cbind(my_missense_mutation, new)
# my_missense_mutation <- na.omit(my_missense_mutation)
# write.csv(my_missense_mutation, "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/my_missense_mutation.csv", row.names = F)

# 读入已处理好的文件
my_missense_mutation <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/my_missense_mutation.csv")

# 筛选有 netMHCpan 结果的病人
my_missense_mutation_413 <- my_missense_mutation[my_missense_mutation$Tumor_Sample_Barcode %in% my_common_ccname,]
length(unique(my_missense_mutation_413$Tumor_Sample_Barcode)) # 413

# --- 2.7 汇总用于计算的最终数据 ---

# 确保同义和非同义突变数据病人一致
any(!unique(my_synonymous_data_413$Tumor_Sample_Barcode) %in% unique(my_missense_mutation_413$Tumor_Sample_Barcode)) # FALSE
any(!unique(my_missense_mutation_413$Tumor_Sample_Barcode) %in% unique(my_synonymous_data_413$Tumor_Sample_Barcode)) # FALSE

# 合并新抗原预测结果 (my_final_df) 到非同义突变数据
my_missense_mutation_413 <- left_join(my_missense_mutation_413, my_final_df, by = "mutation_no")

# (原始代码中的调试块)
if(F){
  any(is.na(my_missense_mutation_413$antigen))
  a <- setdiff(my_missense_mutation_413$mutation_no, my_final_df$mutation_no)
  tmp <- my_missense_mutation_413[my_missense_mutation_413$mutation_no %in% a,]
  my_missense_mutation_413$antigen[is.na(my_missense_mutation_413$antigen)] <- 0
}
# (可选) 写出中间文件
# write.csv(my_missense_mutation_413, "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/my_missense_mutation_add_antigen_413.csv", row.names = F)
# my_missense_mutation_413 <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/R_out_file/my_missense_mutation_add_antigen_413.csv")

# 格式化并合并同义和非同义突变数据
my_missense_mutation_413_new <- my_missense_mutation_413[, c("Chromosome", "i.Start_Position", "i.End_Position", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "before", "after", "spectra", "antigen")]
my_missense_mutation_413_new$Variant_Classification <- 'Missense_Mutation'

my_synonymous_data_413_new <- my_synonymous_data_413
my_synonymous_data_413_new$antigen <- 0 # 同义突变的新抗原定义为0
my_synonymous_data_413_new$Variant_Classification <- 'Silent'

# 合并
my_SNVdata_413 <- rbind(my_synonymous_data_413_new, my_missense_mutation_413_new)
colnames(my_SNVdata_413) <- c("Chromosome", "i.Start_Position", "i.End_Position", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "before", "after", "spectrums", "antigen", "Variant_Classification")

print("第2部分: “自己队列”数据准备完成。")


################################################################################
#
# 第 3 部分: 免疫编辑得分 (R) 计算与绘图
#
################################################################################

print("第3部分: 开始计算免疫编辑得分...")

# --- 3.1 加载病人分组信息 (例如年龄) ---
subgroup436 <- read.csv("./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/subgroup_436.csv", sep = ",")
subgroup436 <- subgroup436 %>%
  mutate(group = case_when(
    age <= 40 ~ "<=40y",
    age > 40 & age < 50 ~ "40-50y",
    age >= 50 & age <= 55 ~ "50-55y",
    age > 55 & age < 65 ~ "55-65y",
    age >= 65 ~ ">=65y"))

# --- 3.2 方法1: 免疫编辑计算 (Rooney et al. 2015 原始方法) ---

# R = (Bobs/Nobs) / (Bpred/Npred)
# Npred 和 Bpred 是基于 *同义 (Silent)* 突变的频谱推导的

immune_editing = function(sample, mutation_3base, TCGA_mut) {
  # 1. 筛选该病人的所有突变
  sample_mutation = mutation_3base[mutation_3base$Tumor_Sample_Barcode == sample,]
  
  # 2. 计算 Npred (预期非同义突变数)
  Npred = 0
  # 遍历该病人的 *同义* 突变频谱
  for (m in unique(sample_mutation[sample_mutation$Variant_Classification == 'Silent', 'spectrums'])) {
    # Npred = (该频谱的同义突变数) * (TCGA该频谱的 Nbackground)
    Npred = Npred + 
      nrow(sample_mutation[sample_mutation$Variant_Classification == 'Silent' & sample_mutation$spectrums == m,]) * as.numeric(TCGA_mut[TCGA_mut$spectrum == m, 'Nbackground'])
  }
  
  # 3. 计算 Bpred (预期新抗原数)
  Bpred = 0
  # 遍历该病人的 *同义* 突变频谱
  for (m in unique(sample_mutation[sample_mutation$Variant_Classification == 'Silent', 'spectrums'])) {
    # Bpred = (该频谱的同义突变数) * (Nbackground) * (Bbackground)
    Bpred = Bpred + 
      nrow(sample_mutation[sample_mutation$Variant_Classification == 'Silent' & sample_mutation$spectrums == m,]) * as.numeric(TCGA_mut[TCGA_mut$spectrum == m, 'Nbackground']) * as.numeric(TCGA_mut[TCGA_mut$spectrum == m, 'Bbackground'])
  }
  
  # 4. 计算 Bobs (观测新抗原数) 和 Nobs (观测非同义突变数)
  Bobs = nrow(sample_mutation[sample_mutation$antigen == '1',])
  Nobs = nrow(sample_mutation[sample_mutation$Variant_Classification == 'Missense_Mutation',])
  
  # 5. 计算 R (免疫编辑得分)
  R = (Bobs / Nobs) / (Bpred / Npred)
  
  out_df = data.frame(sample = sample, Bobs = Bobs, Nobs = Nobs, Bpred = Bpred, Npred = Npred, R = R)
  return(out_df)
}

# --- 3.3 并行计算 (方法1) ---
print("并行计算 (方法1)...")

mutation_3base <- my_SNVdata_413 # "自己队列"的完整突变数据
TCGA_mut <- final_count         # TCGA 背景模型 (方法1)

cl <- makeCluster(4) # 设置使用4个核心
registerDoParallel(cl)

# 按年龄分组并行计算
age_groups <- c("<=40y", "40-50y", "50-55y", "55-65y", ">=65y")
for (type in age_groups) {
  
  samplelist = subgroup436[subgroup436$group == type, "Tumor_Sample_Barcode"]
  
  # 使用 foreach 进行并行循环
  mymut = foreach(i = 1:length(samplelist), .combine = rbind) %dopar% {
    tryCatch(
      {
        result <- immune_editing(samplelist[i], mutation_3base, TCGA_mut)
        if (!is.data.frame(result) || nrow(result) == 0) {
          # 返回 NA 结构
          return(data.frame(sample = samplelist[i], Bobs = NA, Nobs = NA, Bpred = NA, Npred = NA, R = NA))
        }
        return(result)
      },
      error = function(e) {
        # 错误处理
        data.frame(sample = samplelist[i], Bobs = NA, Nobs = NA, Bpred = NA, Npred = NA, R = NA)
      }
    )
  }
  
  mymut$age_group = type
  write.csv(mymut, paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/result/方法1/mymut', type, '.csv'), row.names = F)
}
stopCluster(cl) # 关闭并行集群

# --- 3.4 绘图 (方法1) ---
print("绘图 (方法1)...")

mymut_all_list = list.files(path = "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/result/方法1/", full.names = T, pattern = "csv")
mymut_all_file = lapply(mymut_all_list, read.csv)
mymut = do.call(rbind, mymut_all_file)
mymut <- na.omit(mymut) # 移除计算失败的 NA 值

# 定义颜色和因子顺序
my_cols = c("#2C91E0", "#3ABF99", "#F0A73A", "#D2A8B2", "#886caf")
names(my_cols) = age_groups
mymut$age_group = factor(mymut$age_group, levels = age_groups)

# 绘图
ggplot(mymut, aes(x = age_group, y = R, color = age_group)) +
  stat_boxplot(geom = "errorbar", width = 0.2) + # 添加误差线
  geom_boxplot(outlier.fill = "white", outlier.color = "white") + # 箱线图
  scale_color_manual(values = my_cols) + # 自定义颜色
  geom_jitter(position = position_jitter(width = 0.25, height = 0), alpha = 0.5, size = 1) + # 散点
  stat_compare_means(comparisons = list(c("<=40y", "40-50y"), c("<=40y", "50-55y"), c("<=40y", "55-65y"), c("<=40y", ">=65y")),
                     method = "wilcox.test") + # Wilcoxon 检验
  theme_classic() +
  theme(axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold", colour = "black"),
        legend.position = "none") + # 隐藏图例 (NoLegend())
  labs(title = "", x = "Age Group", y = "Immune Editing Score (R)")


# --- 3.5 方法2: 免疫编辑计算 (备用方法) ---

# R = (Bobs/Nobs) / (Bpred/Npred)
# Npred 和 Bpred 是基于 *所有 (同义+非同义)* 突变的频谱推导的
# 注意: 这种方法在生物学上可能不如方法1合理，因为非同义突变频谱本身已受选择

immune_editing_v2 = function(sample, mutation_3base, TCGA_mut) {
  sample_mutation = mutation_3base[mutation_3base$Tumor_Sample_Barcode == sample,]
  
  # Npred: 基于 *所有* 突变的频谱计算
  Npred = 0
  for (m in unique(sample_mutation$spectrums)) {
    Npred = Npred + 
      nrow(sample_mutation[sample_mutation$spectrums == m,]) * as.numeric(TCGA_mut[TCGA_mut$spectrum == m, 'Nbackground'])
  }
  
  # Bpred: 基于 *所有* 突变的频谱计算
  Bpred = 0
  for (m in unique(sample_mutation$spectrums)) {
    Bpred = Bpred + 
      nrow(sample_mutation[sample_mutation$spectrums == m,]) * as.numeric(TCGA_mut[TCGA_mut$spectrum == m, 'Bbackground'])
  }
  
  Bobs = nrow(sample_mutation[sample_mutation$antigen == '1',])
  Nobs = nrow(sample_mutation[sample_mutation$Variant_Classification == 'Missense_Mutation',])
  
  R = (Bobs / Nobs) / (Bpred / Npred)
  
  out_df = data.frame(sample = sample, Bobs = Bobs, Nobs = Nobs, Bpred = Bpred, Npred = Npred, R = R)
  return(out_df)
}

# --- 3.6 并行计算 (方法2) ---
print("并行计算 (方法2)...")

mutation_3base <- my_SNVdata_413 # "自己队列"的完整突变数据
TCGA_mut <- final_count1        # TCGA 背景模型 (方法2)

cl <- makeCluster(4)
registerDoParallel(cl)

for (type in age_groups) {
  
  samplelist = subgroup436[subgroup436$group == type, "Tumor_Sample_Barcode"]
  
  mymut = foreach(i = 1:length(samplelist), .combine = rbind) %dopar% {
    tryCatch(
      {
        result <- immune_editing_v2(samplelist[i], mutation_3base, TCGA_mut)
        if (!is.data.frame(result) || nrow(result) == 0) {
          return(data.frame(sample = samplelist[i], Bobs = NA, Nobs = NA, Bpred = NA, Npred = NA, R = NA))
        }
        return(result)
      },
      error = function(e) {
        data.frame(sample = samplelist[i], Bobs = NA, Nobs = NA, Bpred = NA, Npred = NA, R = NA)
      }
    )
  }
  
  mymut$age_group = type
  write.csv(mymut, paste0('./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/result/方法2/mymut', type, '.csv'), row.names = F)
}
stopCluster(cl)

# --- 3.7 绘图 (方法2) ---
print("绘图 (方法2)...")

mymut_all_list = list.files(path = "./Desktop/早发肝癌/分析/WES分析/初步分析-436/免疫编辑/result/方法2/", full.names = T, pattern = "csv")
mymut_all_file = lapply(mymut_all_list, read.csv)
mymut = do.call(rbind, mymut_all_file)
mymut <- na.omit(mymut)

mymut$age_group = factor(mymut$age_group, levels = age_groups)

ggplot(mymut, aes(x = age_group, y = R, color = age_group)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(outlier.fill = "white", outlier.color = "white") +
  scale_color_manual(values = my_cols) +
  geom_jitter(position = position_jitter(width = 0.25, height = 0), alpha = 0.5, size = 1) +
  stat_compare_means(comparisons = list(c("<=40y", "40-50y"), c("<=40y", "50-55y"), c("<=40y", "55-65y"), c("<=40y", ">=65y")),
                     method = "wilcox.test") +
  theme_classic() +
  theme(axis.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 0.9),
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold", colour = "black"),
        legend.position = "none") +
  labs(title = "", x = "Age Group", y = "Immune Editing Score (R)")

print("脚本执行完毕。")
