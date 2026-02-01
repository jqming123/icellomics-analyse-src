#!/bin/bash

#==============================================================================
# SCRIPT:       prepare_rna_project.sh
# DESCRIPTION:  此脚本用于自动化RNA测序项目的初始数据准备工作。
#               它会创建一个项目目录，复制原始的sra_runid文件，
#               并向该文件中添加项目ID和参考基因组信息。
#               在处理之前，脚本会执行数据完整性检查，核对原始.sra文件
#               数量与sra_runid.txt中列出的Run ID数量是否一致。
#
# USAGE:        直接运行此脚本。在运行前，请根据需要修改脚本开头的
#               PROJECT_NAME 和 REFERENCE_GENOME 变量。
#
# PREREQUISITES:
#               - 原始数据目录结构应为:
#                 ../rna_rawdata/<PROJECT_NAME>/sra_runid.txt
#                 ../rna_rawdata/<PROJECT_NAME>/<sub_dirs_with_sra_files>/
#               - sra_runid.txt 文件每行包含一个Run ID。
#
# DATE:         2025-12-03 
#==============================================================================


# 定义变量
PROJECT_ID="PRJNA599947"
CL_NAME="CHO"
PROJECT_NAME="${PROJECT_ID}_${CL_NAME}"
REFERENCE_GENOME="CH_Ensemble"

cd /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/EQ_jobs
# 创建目录并进入
mkdir -p "${PROJECT_NAME}"
cd "${PROJECT_NAME}"

# 定义原始数据路径
RAWDATA_BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_rawdata/${PROJECT_NAME}"
ORIGINAL_SRA_RUNID_FILE="${RAWDATA_BASE_DIR}/sra_runid.txt"

echo "--- 开始数据完整性检查 ---"

# 检查原始数据目录是否存在
if [ ! -d "${RAWDATA_BASE_DIR}" ]; then
    echo "错误: 原始数据目录 '${RAWDATA_BASE_DIR}' 不存在。请检查路径或项目ID。"
    exit 1
fi

# 检查原始sra_runid.txt文件是否存在
if [ ! -f "${ORIGINAL_SRA_RUNID_FILE}" ]; then
    echo "错误: 原始sra_runid文件 '${ORIGINAL_SRA_RUNID_FILE}' 不存在。请检查路径或项目ID。"
    exit 1
fi

# 1. 统计原始数据目录及其子目录中的 .sra 文件数量
echo "正在统计 '${RAWDATA_BASE_DIR}' 及其子目录中的 .sra 文件..."
SRA_FILE_COUNT=$(find "${RAWDATA_BASE_DIR}" -type f -name "*.sra" | wc -l)

# 2. 统计原始 sra_runid.txt 文件中的 Run ID 数量（即行数）
echo "正在统计 '${ORIGINAL_SRA_RUNID_FILE}' 中的 Run ID 数量..."
RUNID_COUNT=$(wc -l < "${ORIGINAL_SRA_RUNID_FILE}") 

echo "发现 .sra 文件数量: ${SRA_FILE_COUNT}"
echo "发现 sra_runid.txt 中的 Run ID 数量: ${RUNID_COUNT}"

# 比较数量是否一致
if [ "${SRA_FILE_COUNT}" -eq "${RUNID_COUNT}" ]; then
    echo "数据完整性检查通过: .sra 文件数量与 Run ID 数量一致。"
else
    echo "错误: 数据完整性检查失败! .sra 文件数量与 Run ID 数量不匹配。"
    echo "请检查原始数据目录 '${RAWDATA_BASE_DIR}' 和文件 '${ORIGINAL_SRA_RUNID_FILE}'。"
    exit 1 # 如果检查失败，脚本将退出
fi

echo "--- 数据完整性检查结束 ---"
echo "" # 添加空行增加可读性


# 复制文件
cp "../../rna_rawdata/${PROJECT_NAME}/sra_runid.txt" sra_runid_prjid_ref.txt
echo "当前sra_runid_prjid_ref.txt (同sra_runid.txt) :"
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