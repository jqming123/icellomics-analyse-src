#!/bin/bash

# 脚本功能: 为主染色体的SNP和INDEL VCF文件生成Slurm作业脚本，用于VEP变异位点注释。
# 运行方式: 在 resources/src/ 目录下执行 `bash 05_generate_vep_annotation_slurm.sh <PROJECT_NAME> <REF_NAME>`
# 例如: bash 05_generate_vep_annotation_slurm.sh MyProject hg38_Ensembl

# --- 获取参数 ---
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "用法: $0 <PROJECT_NAME> <REF_NAME>" >&2
    echo "例子: $0 my_human_project hg38" >&2
    exit 1
fi

PROJECT_NAME="$1"
REF_NAME="$2"  # <-- 动态获取参考基因组名称
export PROJECT_NAME
export REF_NAME

echo "当前项目名称 (PROJECT_NAME): ${PROJECT_NAME}"
echo "当前参考基因组 (REF_NAME): ${REF_NAME}"

# --- 引入项目配置文件 ---
if [ -f "./config.sh" ]; then
    source ./config.sh
elif [ -f "../config/config.sh" ]; then
    source ../config/config.sh
else
    echo "错误: 找不到配置文件 config.sh" >&2
    exit 1
fi

echo "当前使用的队列 (QUEUE_NAME): ${QUEUE_NAME}"


# --- 路径定义 ---
# 输入: 07_vcf_merged 目录下的主染色体文件
INPUT_DIR="${RESULTS_DIR}/07_vcf_merged"
# 输出: 08_vep_annotated 目录
OUTPUT_DIR="${RESULTS_DIR}/08_vep_annotated"
# 脚本生成目录
JOB_SCRIPT_DIR="${PROJECT_DIR}/02_jobs/vep_annotation_jobs"

# 创建该阶段所需的输出和日志目录
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}/04_vep_annotation/err_file" "${LOG_DIR}/04_vep_annotation/out_file" "${JOB_SCRIPT_DIR}"

# --- 检查 VEP_SPECIES 是否已加载 ---
if [ -z "${VEP_SPECIES}" ]; then
    echo "错误: VEP_SPECIES 未在 config.sh 中定义。请检查您的配置文件。" >&2
    exit 1
fi

echo "正在为 ${PROJECT_NAME} (${VEP_SPECIES}) 生成 VEP 注释作业脚本..."

if [ "$REF_NAME" = "CriGri-PICRH-1.0" ] ; then
    cache_arg="--merged"
else
    cache_arg=""
fi

echo $cache_arg

# 定义要处理的文件类型
variant_types=("snps" "indels")

for type in "${variant_types[@]}"; do
    # 定义路径
    INPUT_VCF_FILE="${INPUT_DIR}/main_chrs_${type}.vcf.gz"
    NORMALIZED_VCF_FILE="${OUTPUT_DIR}/main_chrs_${type}.normalized.vcf.gz"
    OUTPUT_VEP_VCF_FILE="${OUTPUT_DIR}/main_chrs_${type}.vep.vcf.gz"
    
    JOB_SCRIPT_PATH="${JOB_SCRIPT_DIR}/vep_annotate_${type}.sh"
    LOG_BASENAME="vep_annotate_${type}"

    # 检查输入文件是否存在 (main_chrs_snps.vcf.gz)
    if [ ! -f "${INPUT_VCF_FILE}" ]; then
        echo "警告: 未找到输入文件 ${INPUT_VCF_FILE}，跳过 ${type}。" >&2
        continue
    fi
    
    # 生成 Slurm 脚本
    cat <<EOF > "${JOB_SCRIPT_PATH}"
#!/bin/bash
#SBATCH -p ${QUEUE_NAME}
#SBATCH -J VEP_${type}_${PROJECT_NAME}
#SBATCH -o ${LOG_DIR}/04_vep_annotation/${LOG_BASENAME}_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=${THREADS_VEP}
#SBATCH --mem=${MEM_XLARGE}
#SBATCH -t 240:00:00

set -eo pipefail

echo "Job: VEP Annotation - ${type}"
echo "Genome: ${REF_NAME} (${GENOME_ASSEMBLY})"
echo "Species: ${VEP_SPECIES}"
echo "Start time: \$(date)"

# 激活环境
source "${CONDA_PROFILE_PATH}"
conda activate "${VEP_ENV_NAME}"

# 变量定义
REF_GENOME="${REF_GENOME}"
INPUT_VCF="${INPUT_VCF_FILE}"
NORMALIZED_VCF="${NORMALIZED_VCF_FILE}"
OUTPUT_VEP_VCF="${OUTPUT_VEP_VCF_FILE}"

# Step 1: 标准化
echo "--- Step 1: Normalizing variants with bcftools norm ---"
bcftools norm \\
    -f "\${REF_GENOME}" \\
    -m -any \\
    -O z \\
    -o "\${NORMALIZED_VCF}" \\
    "\${INPUT_VCF}"

# Step 2: VEP 注释
echo "--- Step 2: Running VEP annotation ---"
vep \\
    --input_file "\${NORMALIZED_VCF}" \\
    --output_file "\${OUTPUT_VEP_VCF}" \\
    --offline --cache ${cache_arg}\\
    --dir_cache "${VEP_CACHE_DIR}" \\
    --assembly "${GENOME_ASSEMBLY}" \\
    --species "${VEP_SPECIES}" \\
    --fasta "\${REF_GENOME}" \\
    --format vcf \\
    --vcf \\
    --compress_output bgzip \\
    --force_overwrite \\
    --everything \\
    --pick \\
    --fork ${THREADS_VEP} \\
    --buffer_size 5000 \\
    --check_ref \\
    --warning_file "${LOG_DIR}/04_vep_annotation/${LOG_BASENAME}.warnings.txt" \\
    --stats_file "${LOG_DIR}/04_vep_annotation/${LOG_BASENAME}.stats.html" 

# Step 3: 建立索引
echo "--- Step 3: Indexing ---"
tabix -p vcf "\${OUTPUT_VEP_VCF}"

echo "End time: \$(date)"
EOF
    echo "  已成功生成: ${JOB_SCRIPT_PATH}"
done

echo "---- 完成 ----"