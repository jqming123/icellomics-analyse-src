#!/bin/bash
#
# 脚本功能:
#   1. 将一个大型的、已标准化的 VCF 文件分割成多个小块。
#   2. 为每个生成的小 VCF 块自动创建一个对应的 Slurm 作业脚本，用于运行 VEP 注释。
#   3. 此脚本旨在通过 "分而治之" 的策略，帮助定位导致 VEP 失败的具体 VCF 文件区域。
#
# 运行方式:
#   在 `resources/src/` 目录下执行 `bash debug_vep_by_chunking.sh <PROJECT_NAME>`

set -eou pipefail # 遇到错误或未定义变量立即退出

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
# 输入的、已经过 norm 处理的 VCF 文件 (这是要调试的目标文件)
SOURCE_VCF="${RESULTS_DIR}/08_vep_annotated/main_chrs_snps.normalized.vcf.gz"
# 每个块的行数
# LINES_PER_CHUNK=1000000
LINES_PER_CHUNK=100000

# --- 自动生成的路径 ---
# 用于存放分割后文件、脚本和日志的调试目录
DEBUG_DIR="${RESULTS_DIR}/08_vep_annotated/vep_debug_chunks"
# 存放 VCF 块的目录
SPLIT_DIR="${DEBUG_DIR}/vcf_chunks"
# 存放为每个块生成的 Slurm 作业脚本的目录
JOB_SCRIPT_DIR="${PROJECT_DIR}/02_jobs/vep_debug_jobs"
# 为调试任务创建专用的日志目录
DEBUG_LOG_DIR="${LOG_DIR}/04_vep_annotation_debug"

echo "创建调试目录..."
mkdir -p "${SPLIT_DIR}"
mkdir -p "${JOB_SCRIPT_DIR}"
mkdir -p "${DEBUG_LOG_DIR}/out_file"
mkdir -p "${DEBUG_LOG_DIR}/err_file"

# 检查源 VCF 文件是否存在
if [ ! -f "${SOURCE_VCF}" ]; then
    echo "错误: 源 VCF 文件不存在: ${SOURCE_VCF}" >&2
    exit 1
fi

# --- VEP 特定配置 (从 config.sh 继承) ---
VEP_SPECIES="cricetulus_griseus_picr" # 根据您的缓存名称确定

# --- 步骤 2: 提取 VCF 头部 ---
HEADER_FILE="${SPLIT_DIR}/vcf_header.txt"
echo "正在提取 VCF 头部到: ${HEADER_FILE}"
zcat "${SOURCE_VCF}" | grep '^#' > "${HEADER_FILE}"
if [ ! -s "${HEADER_FILE}" ]; then
    echo "错误: 未能从 ${SOURCE_VCF} 提取到头部信息，或头部为空。" >&2
    exit 1
fi

# --- 步骤 3: 分割 VCF 数据体 ---
echo "正在将 VCF 文件分割成每块 ${LINES_PER_CHUNK} 行..."
# -l: 按行数分割
# -d: 使用数字后缀 (00, 01, ...)
# -a 3: 后缀长度为3位 (000, 001, ...)
# --additional-suffix=.vcf: 给分割的文件加上 .vcf 后缀
zcat "${SOURCE_VCF}" | grep -v '^#' | split -l ${LINES_PER_CHUNK} -d -a 3 --additional-suffix=.vcf - "${SPLIT_DIR}/chunk_"

# --- 步骤 4: 为每个块添加头部、压缩并生成 VEP 作业脚本 ---
echo "正在为每个块添加头部、压缩并生成 Slurm 作业脚本..."
for chunk_file in ${SPLIT_DIR}/chunk_*.vcf; do
    CHUNK_BASENAME=$(basename "${chunk_file}" .vcf) # 例如: chunk_000
    CHUNK_VCF_GZ="${chunk_file}.gz"
    
    echo "--- 处理: ${CHUNK_BASENAME} ---"

    # 4.1: 将头部和数据块合并，然后用 bgzip 压缩
    echo "  -> 压缩文件: ${CHUNK_VCF_GZ}"
    cat "${HEADER_FILE}" "${chunk_file}" | bgzip -c > "${CHUNK_VCF_GZ}"
    rm "${chunk_file}" # 删除未压缩的临时文件

    # 4.2: 为这个块生成 VEP Slurm 作业脚本
    JOB_SCRIPT_PATH="${JOB_SCRIPT_DIR}/${CHUNK_BASENAME}_vep.sh"
    LOG_BASENAME="vep_debug_${CHUNK_BASENAME}"
    OUTPUT_VEP_VCF_FILE="${SPLIT_DIR}/${CHUNK_BASENAME}.vep.vcf.gz"

    # 使用 here-document 将多行内容写入 Slurm 脚本
    cat <<EOF > "${JOB_SCRIPT_PATH}"
#!/bin/bash
#SBATCH -p ${QUEUE_NAME}
#SBATCH -J VEP_${CHUNK_BASENAME}
#SBATCH -o ${DEBUG_LOG_DIR}/out_file/${LOG_BASENAME}.out
#SBATCH -e ${DEBUG_LOG_DIR}/err_file/${LOG_BASENAME}.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS_VEP}
#SBATCH --mem=${MEM_XLARGE}
#SBATCH -t 24:00:00 # 调试任务通常较快，可适当缩短时间

set -eo pipefail

LOGFILE="${DEBUG_LOG_DIR}/${LOG_BASENAME}.log"
touch "\${LOGFILE}"
exec > "\${LOGFILE}" 2>&1

echo "Job started on: \$(hostname)"
echo "Start time: " && date
echo "Annotating VCF chunk ${CHUNK_BASENAME} with VEP..."

source "${CONDA_PROFILE_PATH}"
conda activate "${VEP_ENV_NAME}"

# --- 运行 VEP 注释 ---
# 注意：由于输入文件已经是标准化的，我们直接运行 VEP，跳过 bcftools norm 步骤。
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
    --warning_file "${DEBUG_LOG_DIR}/${LOG_BASENAME}.warnings.txt" \\
    --stats_file "${DEBUG_LOG_DIR}/${LOG_BASENAME}.stats.html"

# --- 索引注释后的 VCF 文件 ---
echo "Indexing annotated VCF file..."
tabix -p vcf "${OUTPUT_VEP_VCF_FILE}"

echo "End time: " && date
echo "Job finished for chunk ${CHUNK_BASENAME}."
EOF

    echo "  -> 已生成作业脚本: ${JOB_SCRIPT_PATH}"

done

# --- 步骤 5: 完成 ---
echo ""
echo "=========================================================="
echo "分割与脚本生成完成！"
echo "所有 VCF 块文件存放在: ${SPLIT_DIR}"
echo "所有对应的 Slurm 作业脚本存放在: ${JOB_SCRIPT_DIR}"
echo "=========================================================="