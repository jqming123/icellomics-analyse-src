#!/bin/bash

# ==============================================================================
# 脚本名称: generate_DEG_allsteps.sh
# 功能描述: 整合 mergecount.sh, DEG_and_basic_img.R, advance_img.R
#          自动生成 Slurm 提交脚本，并支持自定义节点和内存参数
# 使用方法:
#   1. 单个生成: bash generate_DEG_allsteps.sh <GROUP_NAME> <CL_NAME>
#   2. 批量生成: bash generate_DEG_allsteps.sh <config_file>
# ==============================================================================

# ==================== 1. 基础路径配置 ====================
# 原始脚本存放目录
SRC_SCRIPT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/trs_script"
# Conda 环境路径
CONDA_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/bin/activate"
CONDA_ENV="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/R_DESeq2"

# 作业与日志根目录
JOB_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/jobs"
LOG_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/DEG/logs"

# ==================== 2. Slurm 资源配置 (在此修改参数) ====================
PARTITION="corexd192"      # 队列/分区名称 (-p)
MEMORY="64G"             # 内存限制 (--mem)
CPU_CORES="8"            # CPU 核心数 (-c)
TIME_LIMIT="24:00:00"    # 时间限制 (可选)

# ==================== 3. 内部函数: 生成 Slurm 脚本 ====================
generate_job() {
    local G_NAME=$1
    local C_NAME=$2

    # 路径准备
    local JOB_DIR="${JOB_ROOT}/${C_NAME}"
    local LOG_DIR="${LOG_ROOT}/${C_NAME}"
    local JOB_FILE="${JOB_DIR}/${G_NAME}.sh"
    local LOG_FILE="${LOG_DIR}/${G_NAME}.log"

    mkdir -p "${JOB_DIR}"
    mkdir -p "${LOG_DIR}"

    # 写入 Slurm 脚本内容
    cat <<EOF > "${JOB_FILE}"
#!/bin/bash
#SBATCH -J DEG_${G_NAME}
#SBATCH -o ${LOG_FILE}
#SBATCH -p ${PARTITION}
#SBATCH --mem=${MEMORY}
#SBATCH -n 1
#SBATCH -c ${CPU_CORES}
#SBATCH -t ${TIME_LIMIT}

# --- 1. 环境准备 ---
echo "Task started at: \$(date)"
echo "Running on node: \$(hostname)"
source "${CONDA_PATH}" "${CONDA_ENV}"

# --- 2. 运行 mergecount.sh ---
echo "==> Step 1: Merging RSEM counts..."
bash "${SRC_SCRIPT_DIR}/mergecount.sh" "${G_NAME}" "${C_NAME}"

# --- 3. 运行 DEG_and_basic_img.R ---
echo "==> Step 2: Running DESeq2 and basic plotting..."
Rscript "${SRC_SCRIPT_DIR}/DEG_and_basic_img.R" "${G_NAME}" "${C_NAME}"

# --- 4. 运行 advance_img.R ---
echo "==> Step 3: Running advanced visualization..."
Rscript "${SRC_SCRIPT_DIR}/advance_img.R" "${G_NAME}" "${C_NAME}"

echo "Task finished at: \$(date)"
EOF

    chmod +x "${JOB_FILE}"
    echo "[SUCCESS] Job script generated: ${JOB_FILE} (Res: ${PARTITION}, ${MEMORY})"
}

# ==================== 4. 主逻辑 ====================

if [ "$#" -eq 1 ] && [ -f "$1" ]; then
    # 批量模式 (处理 TSV 文件)
    echo "[INFO] Entering batch mode using TSV file: $1"
    
    # 使用制表符作为分隔符读取文件
    while IFS=$'\t' read -r G_NAME C_NAME || [[ -n "$G_NAME" ]]; do
        # 清理变量：去除可能存在的 Windows 回车符及首尾空格
        G_NAME=$(echo "$G_NAME" | tr -d '\r' | xargs)
        C_NAME=$(echo "$C_NAME" | tr -d '\r' | xargs)

        # 忽略空行和以 # 开头的注释行
        [[ -z "${G_NAME}" || "${G_NAME}" =~ ^# ]] && continue
        
        if [[ -n "${G_NAME}" && -n "${C_NAME}" ]]; then
            generate_job "${G_NAME}" "${C_NAME}"
        else
            echo "[WARN] Skipping invalid line (missing columns): G_NAME='$G_NAME', C_NAME='$C_NAME'"
        fi
    done < "$1"

elif [ "$#" -eq 2 ]; then
    # 单个模式
    G_NAME=$1
    C_NAME=$2
    generate_job "${G_NAME}" "${C_NAME}"

else
    # 错误提示
    echo "错误: 参数输入不正确。"
    echo "用法:"
    echo "  单个任务: bash $0 <GROUP_NAME> <CL_NAME>"
    echo "  批量任务: bash $0 <config_file>"
    exit 1
fi