#!/bin/tcsh

################################################################################
# 1. 环境配置 (请根据您的实际路径修改)
################################################################################

# !!! 请修改这里: 设置您的 netMHCpan-4.2b 程序的完整路径
# (这是您解压 netMHCpan-4.2bstatic.linux.tar.gz 后得到的目录)
set NETMHCpan_PATH = "/data01/home/pankaiyan/tools/netMHCpan-4.2"

# --- 您的项目主目录 (根据您的 R 脚本自动推断) ---
set PROJECT_DIR = "/data01/home/pankaiyan/ImmuneEditing"
set HLA_DIR = "$PROJECT_DIR/netMHCpan_input_hla"

# --- 定义要处理的肽段长度列表 ---
set PEPTIDE_LENGTHS = (8 9 10)

################################################################################
# 2. 检查配置
################################################################################

# 检查 netMHCpan 可执行文件是否存在
if (! -f "$NETMHCpan_PATH/netMHCpan") then
    echo "错误: 在以下路径找不到 netMHCpan 可执行文件:"
    echo "$NETMHCpan_PATH/netMHCpan"
    echo "请检查脚本第 7 行的 NETMHCpan_PATH 变量是否设置正确。"
    exit 1
endif

# 检查 HLA 目录是否存在
if (! -d "$HLA_DIR") then
    echo "错误: 找不到 HLA 目录: $HLA_DIR"
    echo "请检查 PROJECT_DIR 变量或确认您的 R 脚本已成功生成 HLA 文件。"
    exit 1
endif

echo "--- 配置检查通过 ---"
echo "netMHCpan 路径: $NETMHCpan_PATH"
echo "项目主目录: $PROJECT_DIR"
echo "将处理肽段长度: $PEPTIDE_LENGTHS"
echo "-------------------------"

################################################################################
# 3. 主循环 (自动处理指定长度的肽段)
################################################################################

# 记录总体开始时间
set TOTAL_START_TIME = `date +%s`

# 外层循环，遍历预定义的肽段长度列表
foreach PEPTIDE_LEN ($PEPTIDE_LENGTHS)

    # --- 动态设置当前长度的输入和输出目录 ---
    set PEPTIDE_DIR = "$PROJECT_DIR/netMHCpan_input_aa${PEPTIDE_LEN}"
    set OUT_DIR = "$PROJECT_DIR/netMHCpan_output_aa${PEPTIDE_LEN}"

    # 检查输入目录是否存在
    if (! -d "$PEPTIDE_DIR") then
        echo "警告: 找不到 ${PEPTIDE_LEN}-mer 的输入目录: $PEPTIDE_DIR"
        echo "跳过 ${PEPTIDE_LEN}-mer 的分析..."
        continue
    endif
    
    # 创建输出目录
    mkdir -p $OUT_DIR
    
    # 统计当前长度的文件总数
    set total_files = `ls -1 "$PEPTIDE_DIR"/*.fsa |& wc -l`
    
    if ($total_files == 0) then
        echo "警告: ${PEPTIDE_LEN}-mer 目录中没有找到 .fsa 文件，跳过..."
        continue
    endif
    
    echo "" # 换行
    echo "================================================="
    echo "=== 开始处理 ${PEPTIDE_LEN}-mer 预测 ==="
    echo "================================================="
    echo "输入肽段目录: $PEPTIDE_DIR"
    echo "输出目录: $OUT_DIR"
    echo "待处理文件数: $total_files"
    echo "-------------------------------------------------"

    # 记录当前长度的开始时间
    set LENGTH_START_TIME = `date +%s`
    
    # 初始化计数器
    set current_count = 0
    set success_count = 0
    set skip_count = 0

    # 内层循环: 处理 $PEPTIDE_DIR 目录中的每一个 .fsa 文件
    foreach fasta_file ($PEPTIDE_DIR/*.fsa)
    
        # 检查是否找到了文件
        if (! -f "$fasta_file") continue

        # 更新计数器
        @ current_count++

        # 从 FASTA 文件路径中提取病人ID
        set patient_id = `basename "$fasta_file" .fsa`

        # 构建对应的 HLA 文件路径
        set hla_file = "$HLA_DIR/${patient_id}_hla.txt"

        # 检查 HLA 文件是否存在
        if (! -f "$hla_file") then
            echo "[$current_count/$total_files] 警告: 找不到 $patient_id 的 HLA 文件，跳过..."
            @ skip_count++
            continue
        endif

        # (关键步骤) 从 HLA 文件中读取逗号分隔的等位基因列表
        set alleles_list = `cat "$hla_file"`

        # 构建输出文件的路径
        set output_file = "$OUT_DIR/${patient_id}_netMHCpan_out.xls"

        # 显示进度
        echo "[$current_count/$total_files] 正在处理 ${PEPTIDE_LEN}-mer: $patient_id ..."

        # 运行 netMHCpan 命令
        $NETMHCpan_PATH/netMHCpan -f "$fasta_file" \
                                 -a "$alleles_list" \
                                 -l $PEPTIDE_LEN \
                                 -s \
                                 -xls \
                                 -xlsfile "$output_file"
        # -f <fasta_file> 指定输入序列文件（FASTA 格式）。
        # -a <alleles_list> 指定要预测的 HLA/MHC 等位基因列表。
        # -l <PEPTIDE_LEN> 指定要预测的肽段长度，例如 8、9、10。
        # -s 表示使用 shell-friendly output（简化输出）。
        # -xls 将输出格式设置为 Excel（tab 分隔）风格。
        # -xlsfile <output_file> 指定输出文件的路径和文件名（通常是 .xls 或 .txt）。
  
        # 检查命令是否成功执行
        if ($status == 0) then
            @ success_count++
        endif
    end

    # 计算当前长度的处理时间
    set LENGTH_END_TIME = `date +%s`
    @ LENGTH_DURATION = $LENGTH_END_TIME - $LENGTH_START_TIME
    
    echo ""
    echo "--- ${PEPTIDE_LEN}-mer 处理完成 ---"
    echo "成功: $success_count | 跳过: $skip_count | 总计: $total_files"
    echo "耗时: ${LENGTH_DURATION} 秒"
    echo ""

end

# 计算总耗时
set TOTAL_END_TIME = `date +%s`
@ TOTAL_DURATION = $TOTAL_END_TIME - $TOTAL_START_TIME
@ TOTAL_MINUTES = $TOTAL_DURATION / 60

echo ""
echo "================================================="
echo "=== 所有分析已完成 ($PEPTIDE_LENGTHS) ==="
echo "================================================="
echo "总耗时: ${TOTAL_DURATION} 秒 (${TOTAL_MINUTES} 分钟)"
echo ""
