#!/bin/bash
#
# 脚本功能:
#   对一个已知有问题的 VCF 块 (例如 chunk_007.vcf.gz) 进行重叠分块。
#   为每个重叠块生成 VEP Slurm 作业脚本，以捕捉跨边界的上下文错误。
#
# 运行方式:
#   在 `resources/src/` 目录下执行 `bash debug_vep_overlapping_chunks.sh <PROJECT_NAME>`

set -eou pipefail

# --- 步骤 0: 获取项目名称并加载配置 ---
if [ -z "$1" ]; then
    echo "用法: $0 <PROJECT_NAME>" >&2
    echo "请提供项目名称作为第一个参数。" >&2
    exit 1
fi
PROJECT_NAME="$1"
export PROJECT_NAME

echo "当前项目名称 (PROJECT_NAME): ${PROJECT_NAME}"

# 引入项目配置文件
if [ -f "./config.sh" ]; then
    source ./config.sh
elif [ -f "../config.sh" ]; then
    source ../config.sh
else
    echo "错误: 找不到配置文件 config.sh" >&2
    exit 1
fi

# --- 步骤 1: 配置与路径定义 ---

# --- 用户可配置的参数 ---
# 输入的有问题的 VCF 块文件
# SOURCE_CHUNK_VCF="${RESULTS_DIR}/08_vep_annotated/vep_debug_chunks/vcf_chunks/chunk_007.vcf.gz"
SOURCE_CHUNK_VCF_LIST=(
    "${RESULTS_DIR}/08_vep_annotated/vep_debug_chunks/vcf_chunks/chunk_006.vcf.gz"
    "${RESULTS_DIR}/08_vep_annotated/vep_debug_chunks/vcf_chunks/chunk_008.vcf.gz"
)
# 每个块的大小 (行数)
CHUNK_SIZE=5000
# 每个块之间的重叠大小 (行数)
OVERLAP_SIZE=500

# --- 主处理循环: 遍历每个有问题的 VCF 文件 ---
for SOURCE_CHUNK_VCF in "${SOURCE_CHUNK_VCF_LIST[@]}"; do

    echo ""
    echo "======================================================================"
    echo "开始处理文件: ${SOURCE_CHUNK_VCF}"
    echo "======================================================================"

    # --- 自动生成的路径 (基于当前处理的文件) ---
    # 从输入文件名中提取基本名称 (例如, "chunk_007")
    CHUNK_BASENAME_ORIGINAL=$(basename "${SOURCE_CHUNK_VCF}" .vcf.gz)
    
    # 为当前文件创建唯一的调试目录，以防结果覆盖
    OVERLAP_DEBUG_DIR="${RESULTS_DIR}/08_vep_annotated/vep_debug_overlap_chunks/${CHUNK_BASENAME_ORIGINAL}"
    OVERLAP_VCF_DIR="${OVERLAP_DEBUG_DIR}/vcf_chunks"
    OVERLAP_JOB_DIR="${OVERLAP_DEBUG_DIR}/job_scripts"
    OVERLAP_LOG_DIR="${LOG_DIR}/04_vep_annotation_overlap_debug"

    echo "创建针对 ${CHUNK_BASENAME_ORIGINAL} 的调试目录..."
    mkdir -p "${OVERLAP_VCF_DIR}"
    mkdir -p "${OVERLAP_JOB_DIR}"
    mkdir -p "${OVERLAP_LOG_DIR}/out_file"
    mkdir -p "${OVERLAP_LOG_DIR}/err_file"

    # 检查源 VCF 文件是否存在
    if [ ! -f "${SOURCE_CHUNK_VCF}" ]; then
        echo "错误: 源 VCF 块文件不存在: ${SOURCE_CHUNK_VCF}" >&2
        echo "跳过此文件..."
        continue # 继续处理列表中的下一个文件
    fi

    # VEP 特定配置
    VEP_SPECIES="cricetulus_griseus_picr"

    # --- 步骤 2: 准备数据 ---
    # 提取头部和数据体
    HEADER_FILE="${OVERLAP_DEBUG_DIR}/vcf_header.txt"
    BODY_FILE="${OVERLAP_DEBUG_DIR}/vcf_body.txt"
    echo "正在提取 VCF 头部和数据体..."
    zcat "${SOURCE_CHUNK_VCF}" | grep '^#' > "${HEADER_FILE}"
    zcat "${SOURCE_CHUNK_VCF}" | grep -v '^#' > "${BODY_FILE}"
    TOTAL_LINES=$(wc -l < "${BODY_FILE}")
    echo "数据体总行数: ${TOTAL_LINES}"

    # --- 步骤 3: 创建重叠块并生成脚本 ---
    echo "正在创建重叠块并生成 Slurm 作业脚本..."
    STEP_SIZE=$((CHUNK_SIZE - OVERLAP_SIZE))
    CHUNK_INDEX=0

    for (( START_LINE=1; START_LINE<=TOTAL_LINES; START_LINE+=STEP_SIZE )); do
        # 生成更具描述性的块名称，包含原始文件名
        CHUNK_BASENAME=$(printf "${CHUNK_BASENAME_ORIGINAL}_overlap_%03d" ${CHUNK_INDEX})
        CHUNK_VCF_PATH="${OVERLAP_VCF_DIR}/${CHUNK_BASENAME}.vcf"
        CHUNK_VCF_GZ="${CHUNK_VCF_PATH}.gz"

        echo "--- 处理: ${CHUNK_BASENAME} (行 ${START_LINE} - $((START_LINE + CHUNK_SIZE - 1))) ---"

        # 3.1: 从数据体文件中提取当前块的内容
        sed -n "${START_LINE},$((START_LINE + CHUNK_SIZE - 1))p; $((START_LINE + CHUNK_SIZE))q" "${BODY_FILE}" > "${CHUNK_VCF_PATH}"
        
        if [ ! -s "${CHUNK_VCF_PATH}" ]; then
            rm "${CHUNK_VCF_PATH}"
            break
        fi

        # 3.2: 添加头部并压缩
        cat "${HEADER_FILE}" "${CHUNK_VCF_PATH}" | bgzip -c > "${CHUNK_VCF_GZ}"
        rm "${CHUNK_VCF_PATH}"

        # 3.3: 生成 VEP 作业脚本
        JOB_SCRIPT_PATH="${OVERLAP_JOB_DIR}/${CHUNK_BASENAME}_vep.sh"
        LOG_BASENAME="vep_debug_${CHUNK_BASENAME}"
        OUTPUT_VEP_VCF_FILE="${OVERLAP_VCF_DIR}/${CHUNK_BASENAME}.vep.vcf.gz"

        cat <<EOF > "${JOB_SCRIPT_PATH}"
#!/bin/bash
#SBATCH -p ${QUEUE_NAME}
#SBATCH -J VEP_${CHUNK_BASENAME}
#SBATCH -o ${OVERLAP_LOG_DIR}/out_file/${LOG_BASENAME}.out
#SBATCH -e ${OVERLAP_LOG_DIR}/err_file/${LOG_BASENAME}.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS_VEP}
#SBATCH --mem=40G
#SBATCH -t 24:00:00

set -eo pipefail

LOGFILE="${OVERLAP_LOG_DIR}/${LOG_BASENAME}.log"
touch "\${LOGFILE}"
exec > "\${LOGFILE}" 2>&1

echo "Job started on: \$(hostname) for chunk ${CHUNK_BASENAME}"
echo "Start time: " && date

source "${CONDA_PROFILE_PATH}"
conda activate "${VEP_ENV_NAME}"

vep \\
    --input_file "${CHUNK_VCF_GZ}" \\
    --output_file "${OUTPUT_VEP_VCF_FILE}" \\
    --offline --cache --merged \\
    --dir_cache "${VEP_CACHE_DIR}" \\
    --assembly "${GENOME_ASSEMBLY}" \\
    --species "${VEP_SPECIES}" \\
    --fasta "${REF_GENOME}" \\
    --format vcf \\
    --vcf \\
    --compress_output bgzip \\
    --force_overwrite \\
    --everything \\
    --pick \\
    --fork ${THREADS_VEP} \\
    --buffer_size 5000 \\
    --check_ref \\
    --warning_file "${OVERLAP_LOG_DIR}/${LOG_BASENAME}.warnings.txt" \\
    --stats_file "${OVERLAP_LOG_DIR}/${LOG_BASENAME}.stats.html"

tabix -p vcf "${OUTPUT_VEP_VCF_FILE}"

echo "End time: " && date
echo "Job finished for chunk ${CHUNK_BASENAME}."
EOF

        echo "  -> 已生成作业脚本: ${JOB_SCRIPT_PATH}"
        CHUNK_INDEX=$((CHUNK_INDEX + 1))
    done

    # --- 步骤 4: 清理当前文件的临时文件 ---
    rm "${BODY_FILE}"
    rm "${HEADER_FILE}"

    echo ""
    echo "--------------------------------------------------------"
    echo "文件 ${SOURCE_CHUNK_VCF} 处理完成。"
    echo "VCF 重叠块文件存放在: ${OVERLAP_VCF_DIR}"
    echo "对应的 Slurm 作业脚本存放在: ${OVERLAP_JOB_DIR}"
    echo "请提交为该文件生成的作业："
    echo "for script in ${OVERLAP_JOB_DIR}/*.sh; do sbatch \"\${script}\"; done"
    echo "--------------------------------------------------------"

done # 结束对 VCF 文件列表的循环

echo ""
echo "=========================================================="
echo "所有文件的重叠块分割与脚本生成完成！"
echo "请检查上面每个文件的输出，并分别提交生成的作业。"
echo "=========================================================="
echo ""