#!/bin/bash

#==============================================================================
# SCRIPT:       write_sra_info.sh
# DESCRIPTION:  此脚本用于自动化RNA测序数据分析的初始数据准备工作。
#               它会创建一个项目目录，复制原始的sra_runid文件，
#               并向该文件中添加项目ID和参考基因组信息。
#               在处理之前，脚本会执行数据完整性检查，核对原始.sra文件
#               数量与sra_runid.txt中列出的Run ID数量是否一致。
#
# USAGE:        ./prepare_rna_project.sh <CL_NAME> <REFERENCE_GENOME> <PROJECT_ID1> [PROJECT_ID2 ...]
#               示例: ./prepare_rna_project.sh PK-15 Pig_E_Sscrofa11.1 PRJNA876757 PRJNA876758
#
# PREREQUISITES:
#               - 原始数据目录结构应为:
#                 ../rna_rawdata/<PROJECT_ID_CL_NAME>/sra_runid.txt
#                 ../rna_rawdata/<PROJECT_ID_CL_NAME>/<sub_dirs_with_sra_files>/
#               - sra_runid.txt 文件每行包含一个Run ID。
#
# DATE:         2026-05-18 
#==============================================================================

# 检查命令行参数数量
if [ $# -lt 3 ]; then
    echo "错误: 参数不足!"
    echo "用法: $0 <CL_NAME> <REFERENCE_GENOME> <PROJECT_ID1> [PROJECT_ID2 ...]"
    echo "示例: $0 PK-15 Pig_E_Sscrofa11.1 PRJNA876757"
    echo "示例: $0 PK-15 Pig_E_Sscrofa11.1 PRJNA876757 PRJNA876758 PRJNA876759"
    exit 1
fi

# 从命令行参数获取变量
CL_NAME="$1"
REFERENCE_GENOME="$2"
shift 2  # 移除前两个参数，剩下的都是PROJECT_ID

# 获取所有PROJECT_ID
PROJECT_IDS=("$@")

echo "================================================"
echo "开始处理RNA测序项目"
echo "细胞系名称: ${CL_NAME}"
echo "参考基因组: ${REFERENCE_GENOME}"
echo "项目ID列表: ${PROJECT_IDS[@]}"
echo "================================================"

# 进入基础工作目录
cd /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/EQ_jobs || {
    echo "错误: 无法进入工作目录!"
    exit 1
}

# 循环处理每个PROJECT_ID
for PROJECT_ID in "${PROJECT_IDS[@]}"; do
    echo ""
    echo "--- 开始处理项目: ${PROJECT_ID} ---"
    
    # 构建项目名称
    PROJECT_NAME="${PROJECT_ID}_${CL_NAME}"
    
    # 创建目录并进入
    mkdir -p "${PROJECT_NAME}"
    cd "${PROJECT_NAME}" || {
        echo "错误: 无法进入项目目录 ${PROJECT_NAME}"
        continue
    }
    
    # 定义原始数据路径
    RAWDATA_BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_rawdata/${PROJECT_NAME}"
    ORIGINAL_SRA_RUNID_FILE="${RAWDATA_BASE_DIR}/sra_runid.txt"
    
    echo "--- 开始数据完整性检查 ---"
    
    # 检查原始数据目录是否存在
    if [ ! -d "${RAWDATA_BASE_DIR}" ]; then
        echo "警告: 原始数据目录 '${RAWDATA_BASE_DIR}' 不存在。跳过此项目。"
        cd ..  # 返回上级目录
        continue
    fi
    
    # 检查原始sra_runid.txt文件是否存在
    if [ ! -f "${ORIGINAL_SRA_RUNID_FILE}" ]; then
        echo "警告: 原始sra_runid文件 '${ORIGINAL_SRA_RUNID_FILE}' 不存在。跳过此项目。"
        cd ..  # 返回上级目录
        continue
    fi
    
    # 1. 统计原始数据目录及其子目录中的 .sra 文件数量
    echo "正在统计 '${RAWDATA_BASE_DIR}' 及其子目录中的 .sra 文件..."
    SRA_FILE_COUNT=$(find "${RAWDATA_BASE_DIR}" -type f -name "*.sra" 2>/dev/null | wc -l)
    
    # 2. 统计原始 sra_runid.txt 文件中的 Run ID 数量（即行数）
    echo "正在统计 '${ORIGINAL_SRA_RUNID_FILE}' 中的 Run ID 数量..."
    RUNID_COUNT=$(grep -c "[^[:space:]]" "${ORIGINAL_SRA_RUNID_FILE}" 2>/dev/null)
    
    echo "发现 .sra 文件数量: ${SRA_FILE_COUNT}"
    echo "发现 sra_runid.txt 中的 Run ID 数量: ${RUNID_COUNT}"
    
    # 比较数量是否一致
    if [ "${SRA_FILE_COUNT}" -eq "${RUNID_COUNT}" ]; then
        echo "数据完整性检查通过: .sra 文件数量与 Run ID 数量一致。"
    else
        echo "警告: 数据完整性检查失败! .sra 文件数量与 Run ID 数量不匹配。"
        echo "请检查原始数据目录 '${RAWDATA_BASE_DIR}' 和文件 '${ORIGINAL_SRA_RUNID_FILE}'。"
        echo "继续处理，但请注意数据可能不完整。"
    fi
    
    echo "--- 数据完整性检查结束 ---"
    
    # 复制文件
    cp "${ORIGINAL_SRA_RUNID_FILE}" sra_runid_prjid_ref.txt
    echo "当前sra_runid_prjid_ref.txt (同sra_runid.txt):"
    cat sra_runid_prjid_ref.txt
    
    # 添加项目ID列
    awk -v OFS='\t' '{print $0, "'"${PROJECT_NAME}"'"}' sra_runid_prjid_ref.txt > temp.txt && mv temp.txt sra_runid_prjid_ref.txt
    
    # 添加参考基因组名称列
    awk -v OFS='\t' '{print $0, "'"${REFERENCE_GENOME}"'"}' sra_runid_prjid_ref.txt > temp.txt && mv temp.txt sra_runid_prjid_ref.txt
    
    echo "增加项目ID列和参考基因组名称列后的sra_runid_prjid_ref.txt:"
    cat sra_runid_prjid_ref.txt
    
    # 统计sra_runid_prjid_ref.txt文件的行数
    echo "最终处理后的 sra_runid_prjid_ref.txt 文件行数:"
    wc -l sra_runid_prjid_ref.txt
    
    echo "--- 项目 ${PROJECT_ID} 处理完成 ---"
    
    # 返回上级目录，准备处理下一个项目
    cd ..
done

echo ""
echo "================================================"
echo "所有项目处理完成!"
echo "================================================"