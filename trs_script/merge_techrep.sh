#!/bin/bash

# ==================== 使用说明 (Usage) ====================
#
# 运行此脚本时，请提供组名和细胞系名称作为参数。
#
# 示例:
# bash mergecount.sh PRJNA974014_a_stressed_Day5 CHO
#
# ==========================================================

# --- 检查输入参数 ---
# 现在需要检查是否提供了至少 2 个参数
if [ "$#" -ne 2 ]; then
    echo "错误: 请提供组名和细胞系名称作为参数。"
    echo "用法: $0 <group_name> <cl_name>"
    echo "示例: $0 PRJNA974014_a_stressed_Day5 CHO"
    exit 1
fi

# 从命令行获取参数
GROUP_NAME="$1"
CL_NAME="$2"


# ==================== 配置区 (Configuration) ====================
# 在这里修改你的基础路径和固定参数

# --- 核心路径 ---
# 脚本所在目录
SCRIPT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/trs_script"
# Conda 环境路径
CONDA_ENV="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/R_DESeq2"
# 包含所有组的根目录
GROUPS_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/groups"

# --- 打印确认信息 (可选) ---
echo "[INFO] Group Name: $GROUP_NAME"
echo "[INFO] Cell Line:  $CL_NAME"

# 输入目录 (RSEM结果所在位置)
INPUT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"

# --- mergeRSEM.py 脚本固定参数 ---
# 计数类型 (例如: expected_count, TPM, FPKM)
COUNT_TYPE="expected_count"
echo "[INFO] Count Type:  $COUNT_TYPE"


# ==================== 动态路径构建 (Dynamic Paths) ====================
# 根据输入的组名自动生成文件路径

# 工作目录
WORK_DIR="${GROUPS_DIR}/${CL_NAME}"

# 技术重复映射文件路径：无表头，两列，第一列run id，第二列合并后的新样本名
SAMPLE_MAP="${WORK_DIR}/${GROUP_NAME}/${GROUP_NAME}.techrep.tsv"
# 输出文件路径
OUTPUT_FILE="${WORK_DIR}/${GROUP_NAME}/${GROUP_NAME}.count.tsv"


# ==================== 执行区 (Execution) ====================

echo "================================================="
echo "组名 (Group Name): ${GROUP_NAME}"
echo "技术重复映射表 (Sample Map): ${SAMPLE_MAP}"
echo "输出文件 (Output File): ${OUTPUT_FILE}"
echo "================================================="

# --- 预执行检查 ---
if [ ! -f "${SAMPLE_MAP}" ]; then
    echo "错误: 找不到技术重复映射文件: ${SAMPLE_MAP}"
    exit 1
fi

echo "==> 切换到工作目录: ${WORK_DIR}"
cd "${WORK_DIR}" || exit 1

echo "==> 激活 Conda 环境..."
eval "$(mamba shell hook --shell bash)"
mamba activate "${CONDA_ENV}"

echo "==> 开始按技术重复映射表合并 RSEM counts..."
python "${SCRIPT_DIR}/mergeRSEM_techrep.py" \
  -i "${INPUT_DIR}" \
  -l "${SAMPLE_MAP}" \
  -o "${OUTPUT_FILE}" \
  -c "${COUNT_TYPE}"

echo "==> 任务完成，正在停用 Conda 环境..."
mamba deactivate

echo "==> 脚本执行完毕"