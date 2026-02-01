#!/bin/bash
#
# 运行方式:
#   在 `resources/src/` 目录下执行 `bash debug_vep_overlapping_chunks.sh <PROJECT_NAME>`
#
# 修改说明:
#   1. SOURCE_CHUNK_VCF 已改为数组 PROBLEM_VCF_CHUNKS，以适应多个源文件。
#   2. 所有引起报错的 VCF 文件的行将汇总到一个文件里。
#      - 主脚本会为每个源文件创建独立的输出目录。
#      - 主脚本会生成一个名为 `collect_problematic_vcf_lines.sh` 的辅助脚本。
#      - 你需要在所有 VEP 作业完成后，手动运行 `collect_problematic_vcf_lines.sh` 来汇总问题行。

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
# 定义最后一次调试的根目录
OVERLAP_DEBUG_ROOT_DIR="${RESULTS_DIR}/08_vep_annotated/vep_debug_final_round"
# 
OC_DIR="${RESULTS_DIR}/08_vep_annotated/vep_debug_overlap_chunks"
# 请在此处添加所有需要处理的 VCF 文件的路径。
declare -a PROBLEM_VCF_CHUNKS=(
    "${OC_DIR}/chunk_006/vcf_chunks/chunk_006_overlap_008.vcf.gz"
    "${OC_DIR}/chunk_008/vcf_chunks/chunk_008_overlap_009.vcf.gz"
)

echo "以下是所有包含 'Died in forked process' 错误的 VCF 块文件路径:"
for vcf_file in "${PROBLEM_VCF_CHUNKS[@]}"; do
    echo "- ${vcf_file}"
    # 您可以在这里添加对这些文件的进一步处理或检查命令
    # 例如：
    # zcat "${vcf_file}" | head -n 20 # 查看文件头部内容
    # zcat "${vcf_file}" | grep -v '^#' | wc -l # 查看数据行数
done

# 2. 将块大小减小到 500 行
CHUNK_SIZE=100
# 3. 将重叠大小减小到 10 行
OVERLAP_SIZE=50

# --- 自动生成的路径 ---
# 汇总所有可能导致 VEP 报错或警告的 VCF 行的最终文件
ALL_PROBLEM_LINES_FILE="${OVERLAP_DEBUG_ROOT_DIR}/all_problematic_vcf_lines.txt"

echo "创建重叠块调试根目录: ${OVERLAP_DEBUG_ROOT_DIR}..."
mkdir -p "${OVERLAP_DEBUG_ROOT_DIR}"
# 确保问题行汇总文件是空的或不存在，避免重复追加
> "${ALL_PROBLEM_LINES_FILE}"

# VEP 特定配置
VEP_SPECIES="cricetulus_griseus_picr"

# --- 循环处理每个源 VCF 文件 ---
# GLOBAL_CHUNK_INDEX 用于为所有生成的块提供一个全局唯一的索引，尽管文件路径已通过 SOURCE_ID 区分
GLOBAL_CHUNK_INDEX=0

for SOURCE_VCF_PATH in "${PROBLEM_VCF_CHUNKS[@]}"; do
    # 从文件名中提取一个安全的标识符，用于创建独立的目录和文件。
    # 例如：从 /path/to/my_file.vcf.gz 提取 my_file
    SOURCE_BASENAME=$(basename "${SOURCE_VCF_PATH}")
    SOURCE_ID=$(echo "${SOURCE_BASENAME}" | sed 's/\.vcf\.gz$//' | sed 's/[^a-zA-Z0-9_.-]/_/g') # 清理文件名，使其适合作为目录名

    echo "--- 开始处理源文件: ${SOURCE_VCF_PATH} (ID: ${SOURCE_ID}) ---"

    # 为当前源文件创建独立的输出目录结构
    OVERLAP_DEBUG_DIR="${OVERLAP_DEBUG_ROOT_DIR}/${SOURCE_ID}"
    OVERLAP_VCF_DIR="${OVERLAP_DEBUG_DIR}/vcf_chunks"
    OVERLAP_JOB_DIR="${OVERLAP_DEBUG_DIR}/job_scripts"
    # 日志目录也特定于每个源文件
    OVERLAP_LOG_DIR="${LOG_DIR}/04_vep_annotation_overlap_debug/${SOURCE_ID}"

    echo "创建当前源文件 (${SOURCE_ID}) 的重叠块调试目录..."
    mkdir -p "${OVERLAP_VCF_DIR}"
    mkdir -p "${OVERLAP_JOB_DIR}"
    mkdir -p "${OVERLAP_LOG_DIR}/out_file"
    mkdir -p "${OVERLAP_LOG_DIR}/err_file"

    # 检查源 VCF 文件是否存在
    if [ ! -f "${SOURCE_VCF_PATH}" ]; then
        echo "错误: 源 VCF 块文件不存在: ${SOURCE_VCF_PATH}" >&2
        continue # 跳过当前文件，处理下一个
    fi

    # --- 步骤 2: 准备数据 ---
    # 提取头部和数据体，文件名为当前源文件 ID 专用
    HEADER_FILE="${OVERLAP_DEBUG_DIR}/vcf_header_${SOURCE_ID}.txt"
    BODY_FILE="${OVERLAP_DEBUG_DIR}/vcf_body_${SOURCE_ID}.txt"
    echo "正在提取 VCF 头部和数据体到 ${OVERLAP_DEBUG_DIR}..."
    zcat "${SOURCE_VCF_PATH}" | grep '^#' > "${HEADER_FILE}"
    zcat "${SOURCE_VCF_PATH}" | grep -v '^#' > "${BODY_FILE}"
    TOTAL_LINES=$(wc -l < "${BODY_FILE}")
    echo "数据体总行数: ${TOTAL_LINES}"

    # --- 步骤 3: 创建重叠块并生成脚本 ---
    echo "正在创建重叠块并生成 Slurm 作业脚本..."
    STEP_SIZE=$((CHUNK_SIZE - OVERLAP_SIZE))
    CHUNK_INDEX=0 # 为每个源文件重置块索引

    for (( START_LINE=1; START_LINE<=TOTAL_LINES; START_LINE+=STEP_SIZE )); do
        # 块的基本名称包含源文件 ID 和当前块索引，确保唯一性
        CHUNK_BASENAME=$(printf "%s_%03d" "${SOURCE_ID}" ${CHUNK_INDEX})
        CHUNK_VCF_PATH="${OVERLAP_VCF_DIR}/${CHUNK_BASENAME}.vcf"
        CHUNK_VCF_GZ="${CHUNK_VCF_PATH}.gz"

        echo "--- 处理: ${CHUNK_BASENAME} (源文件: ${SOURCE_ID}, 行 ${START_LINE} - $((START_LINE + CHUNK_SIZE - 1))) ---"

        # 3.1: 从数据体文件中提取当前块的内容
        # 使用 sed 来精确提取行范围
        sed -n "${START_LINE},$((START_LINE + CHUNK_SIZE - 1))p; $((START_LINE + CHUNK_SIZE))q" "${BODY_FILE}" > "${CHUNK_VCF_PATH}"
        
        # 如果提取的文件为空 (到达末尾)，则停止
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
#SBATCH --mem=20G
#SBATCH -t 4:00:00

set -eo pipefail

LOGFILE="${OVERLAP_LOG_DIR}/${LOG_BASENAME}.log"
touch "\${LOGFILE}"
exec > "\${LOGFILE}" 2>&1

echo "Job started on: \$(hostname) for chunk ${CHUNK_BASENAME} (from source ${SOURCE_ID})"
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
    --buffer_size 1 \\
    --check_ref \\
    --warning_file "${OVERLAP_LOG_DIR}/${LOG_BASENAME}.warnings.txt" \\
    --stats_file "${OVERLAP_LOG_DIR}/${LOG_BASENAME}.stats.html"

# 检查 VEP 命令的退出状态。如果 VEP 失败，则 Slurm 作业也应标记为失败。
if [ \$? -ne 0 ]; then
    echo "错误: VEP 处理 chunk ${CHUNK_BASENAME} 失败！请检查日志文件: \${LOGFILE} 和 Slurm 错误文件。" >> "\${LOGFILE}"
    exit 1 # 确保 Slurm 标记作业失败
fi

tabix -p vcf "${OUTPUT_VEP_VCF_FILE}"

echo "End time: " && date
echo "Job finished for chunk ${CHUNK_BASENAME}."
EOF

        echo "  -> 已生成作业脚本: ${JOB_SCRIPT_PATH}"
        CHUNK_INDEX=$((CHUNK_INDEX + 1))
        GLOBAL_CHUNK_INDEX=$((GLOBAL_CHUNK_INDEX + 1)) # 递增全局索引
    done

    # --- 步骤 4: 清理当前源文件相关的临时文件 ---
    rm "${BODY_FILE}"
    rm "${HEADER_FILE}"
    echo "已清理源文件 ${SOURCE_ID} 的临时文件。"
    echo "--- 源文件 ${SOURCE_ID} 处理完成 ---"
    echo ""
done # 结束对 PROBLEM_VCF_CHUNKS 的循环

echo ""
echo "=========================================================="
echo "重叠块分割与脚本生成完成！"
echo "所有 VCF 重叠块文件和作业脚本按源文件ID存放在: ${OVERLAP_DEBUG_ROOT_DIR}/<SOURCE_ID>/"
echo "所有日志文件按源文件ID存放在: ${LOG_DIR}/04_vep_annotation_overlap_debug/<SOURCE_ID>/"
echo "=========================================================="
echo ""
echo "请提交这些新生成的作业来找到失败的重叠块："
echo "find ${OVERLAP_DEBUG_ROOT_DIR} -name '*_vep.sh' -print0 | xargs -0 -n 1 sbatch"
echo ""
