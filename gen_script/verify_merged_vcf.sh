#!/bin/bash

# 脚本功能: 核对 VCF 合并前后变异位点的数量是否一致。
# 运行方式: 在 resources/src 目录下执行
#           bash verify_merged_vcf.sh <PROJECT_NAME> <type> <scope>
#
# 参数说明:
#   <PROJECT_NAME>: 你的项目名称
#   <type>:         要核对的变异类型, 'snp' 或 'indel'
#   <scope>:        要核对的范围, 'all' (所有) 或 'chrs_only' (仅主染色体)
#
# 示例:
#   bash verify_merged_vcf.sh PRJNA378939_CHO snp all
#   bash verify_merged_vcf.sh PRJNA378939_CHO indel chrs_only

# --- ANSI Color Codes ---
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# --- 参数检查 ---
if [ "$#" -ne 3 ]; then
    echo "用法: $0 <PROJECT_NAME> <type> <scope>" >&2
    echo "  <type>:  'snp' 或 'indel'" >&2
    echo "  <scope>: 'all' 或 'chrs_only'" >&2
    exit 1
fi

PROJECT_NAME="$1"
TYPE="$2"
SCOPE="$3"

echo "当前项目名称 (PROJECT_NAME): ${PROJECT_NAME}"

export REF_NAME="dont_need_ref"

# --- 引入项目配置文件 ---
if [ -f "./config.sh" ]; then
    source ./config.sh
elif [ -f "../config/config.sh" ]; then
    source ../config/config.sh
else
    echo -e "${RED}错误: 找不到配置文件 config.sh${NC}" >&2
    exit 1
fi

# --- 路径定义 ---
INPUT_DIR="${RESULTS_DIR}/06_vcf_filtered"
OUTPUT_DIR="${RESULTS_DIR}/07_vcf_merged"

# --- 根据参数确定要核对的文件 ---
VCF_SUFFIX=""
MERGED_FILE_NAME=""
SOURCE_FILES_PATTERN=""
DESCRIPTION=""

# 设置变异类型
case "$TYPE" in
    snp)
        VCF_SUFFIX="snps"
        ;;
    indel)
        VCF_SUFFIX="indels"
        ;;
    *)
        echo -e "${RED}错误: 无效的变异类型 '$TYPE'. 请使用 'snp' 或 'indel'。${NC}" >&2
        exit 1
        ;;
esac

# 设置范围和文件名
case "$SCOPE" in
#    all)
#        MERGED_FILE_NAME="all_${VCF_SUFFIX}_with_scaffolds.vcf.gz"
#        SOURCE_FILES_PATTERN="${INPUT_DIR}/*.filtered.${VCF_SUFFIX}.vcf.gz"
#        DESCRIPTION="所有 ${VCF_SUFFIX^^} (包括 unplaced scaffolds)"
#        ;;
    chrs_only)
        MERGED_FILE_NAME="main_chrs_${VCF_SUFFIX}.vcf.gz"
        SOURCE_FILES_PATTERN="${INPUT_DIR}/*.filtered.${VCF_SUFFIX}.vcf.gz"
        DESCRIPTION="仅主染色体 ${VCF_SUFFIX^^}"
        ;;
    *)
        echo -e "${RED}错误: 无效的范围 '$SCOPE'. 请使用 'all' 或 'chrs_only'。${NC}" >&2
        exit 1
        ;;
esac

MERGED_FILE_PATH="${OUTPUT_DIR}/${MERGED_FILE_NAME}"

# --- 开始验证 ---
echo -e "\n--- ${YELLOW}开始验证: ${DESCRIPTION}${NC} ---"

# 1. 检查合并后的文件是否存在
if [ ! -f "${MERGED_FILE_PATH}" ]; then
    echo -e "${RED}错误: 合并后的文件不存在: ${MERGED_FILE_PATH}${NC}" >&2
    exit 1
fi

# 2. 统计合并后的VCF文件中的位点数
echo "正在统计合并后文件中的位点数..."
# zgrep -v "^#" 用于排除header行; wc -l 用于计数
COUNT_AFTER=$(zgrep -vc "^#" "${MERGED_FILE_PATH}")
echo "合并后文件 (${MERGED_FILE_NAME}): ${COUNT_AFTER} 个位点"

# 3. 统计合并前所有源文件中的位点总数
echo -e "\n正在统计源文件中的位点总数..."
SOURCE_FILES=($(ls -1 ${SOURCE_FILES_PATTERN}))

if [ ${#SOURCE_FILES[@]} -eq 0 ]; then
    echo -e "${RED}错误: 找不到任何源文件，模式: ${SOURCE_FILES_PATTERN}${NC}" >&2
    exit 1
fi

TOTAL_COUNT_BEFORE=0
for file in "${SOURCE_FILES[@]}"; do
    # 使用 zgrep -c 可以更高效地计数
    count_in_file=$(zgrep -vc "^#" "$file")
    # 使用 ((...)) 进行算术运算
    (( TOTAL_COUNT_BEFORE += count_in_file ))
    # 打印每个文件的计数，方便调试
    # printf "  - %-50s : %d\n" "$(basename "$file")" "$count_in_file"
done
echo "所有源文件合计: ${TOTAL_COUNT_BEFORE} 个位点"


# 4. 比较并报告结果
echo -e "\n--- ${YELLOW}验证结果${NC} ---"
echo "合并前总计 (源文件): ${TOTAL_COUNT_BEFORE}"
echo "合并后总计 (目标文件): ${COUNT_AFTER}"

if [ "${TOTAL_COUNT_BEFORE}" -eq "${COUNT_AFTER}" ]; then
    echo -e "\n${GREEN}✅ 成功: 合并前后的变异位点数量一致！${NC}"
    exit 0
else
    echo -e "\n${RED}❌ 失败: 合并前后的变异位点数量不一致！${NC}"
    diff=$(( TOTAL_COUNT_BEFORE - COUNT_AFTER ))
    echo "差异: ${diff}"
    exit 1
fi