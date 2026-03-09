#!/bin/bash

# 脚本功能: 生成四个 Slurm 作业脚本，用于分别合并 SNP 和 INDEL 的 VCF 文件。
# 运行方式: 在 resources/src 目录下执行 `bash 04_generate_merge_vcf_slurm.sh <PROJECT_NAME> <REF_NAME>`
# 例如: bash 04_generate_merge_vcf_slurm.sh MyProject hg38

# --- 获取参数 ---
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "用法: $0 <PROJECT_NAME> <REF_NAME>" >&2
    echo "例子: $0 my_human_project hg38" >&2
    echo "例子: $0 my_cho_project CriGri-PICRH-1.0" >&2
    exit 1
fi

PROJECT_NAME="$1"
REF_NAME="$2"  # <-- 修改1：从命令行获取参考基因组名称
export PROJECT_NAME
export REF_NAME

echo "当前项目名称: ${PROJECT_NAME}"
echo "当前参考基因组: ${REF_NAME}"

# --- 引入项目配置文件 ---
# 注意：source config.sh 会根据我们上面设置的 REF_NAME 加载对应的 main_chrs 数组
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
INPUT_DIR="${RESULTS_DIR}/06_vcf_filtered"
OUTPUT_DIR="${RESULTS_DIR}/07_vcf_merged"
JOB_SCRIPT_DIR="${PROJECT_DIR}/02_jobs/merge_vcf_jobs"

# 创建目录
mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}/03_merge_vcf" "${JOB_SCRIPT_DIR}"

# --- 识别需要处理的染色体 ---
# <-- 修改2：不再使用 ls NC_* 匹配，直接使用 config.sh 里的 main_chrs 数组 -->
echo "正在检查主染色体 VCF 文件是否存在于 ${INPUT_DIR} ..."

MAIN_CHRS_AVAILABLE=()
for chr in "${main_chrs[@]}"; do
    FILE_PATH="${INPUT_DIR}/${chr}.filtered.snps.vcf.gz"
    if [ -f "$FILE_PATH" ]; then
        MAIN_CHRS_AVAILABLE+=("$chr")
    else
        echo "警告: 未找到染色体 $chr 的文件: $FILE_PATH"
    fi
done

if [ ${#MAIN_CHRS_AVAILABLE[@]} -eq 0 ]; then
    echo "错误: 在 ${INPUT_DIR} 中未找到任何主染色体的 VCF 文件。请检查路径或文件名。" >&2
    exit 1
fi
echo "成功识别到 ${#MAIN_CHRS_AVAILABLE[@]} 个主染色体文件。"

# --- 生成作业脚本 ---
echo "正在生成Slurm作业脚本..."

# Helper function (保持不变)
generate_merge_script() {
    local job_name="$1"
    local output_filename="$2"
    local input_files_str="$3"
    local script_path="$4"
    local log_basename="$5"

    cat <<EOF > "${script_path}"
#!/bin/bash
#SBATCH -p ${QUEUE_NAME}
#SBATCH -J ${job_name}_${PROJECT_NAME}
#SBATCH -o ${LOG_DIR}/03_merge_vcf/${log_basename}_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=${MEM_SMALL}
#SBATCH --time=240:00:00

echo "Job: ${job_name}"
echo "Start time: \$(date)"

source "${CONDA_PROFILE_PATH}"
conda activate "${GENOME_ENV_NAME}"

OUTPUT_DIR="${OUTPUT_DIR}"
FINAL_OUTPUT_VCF="\${OUTPUT_DIR}/${output_filename}"

echo "--- 开始合并VCF文件 ---"
bcftools concat -a -o "\${FINAL_OUTPUT_VCF}" -O z ${input_files_str}

echo "--- 建立索引 ---"
tabix -p vcf "\${FINAL_OUTPUT_VCF}"

echo "End time: \$(date)"
EOF
}


# 后面决定忽略unplaced scaffolds，所以直接跳过了生成unplaced scaffolds对应脚本的步骤
# 1 & 2. 合并所有 (含 Scaffolds) - 逻辑保持不变，因为 ls *.filtered 会匹配所有
# JOB_SCRIPT_SNP_ALL="${JOB_SCRIPT_DIR}/merge_snps_all.sh"
# INPUT_FILES_SNP_ALL="\$(ls -1 ${INPUT_DIR}/*.filtered.snps.vcf.gz | sort -V)"
# generate_merge_script "MergeSNP_All" "all_snps_with_scaffolds.vcf.gz" "${INPUT_FILES_SNP_ALL}" "${JOB_SCRIPT_SNP_ALL}" "merge_snps_all"

# JOB_SCRIPT_INDEL_ALL="${JOB_SCRIPT_DIR}/merge_indels_all.sh"
# INPUT_FILES_INDEL_ALL="\$(ls -1 ${INPUT_DIR}/*.filtered.indels.vcf.gz | sort -V)"
# generate_merge_script "MergeINDEL_All" "all_indels_with_scaffolds.vcf.gz" "${INPUT_FILES_INDEL_ALL}" "${JOB_SCRIPT_INDEL_ALL}" "merge_indels_all"

# 3. 仅合并主染色体 SNPs (使用上面筛选出的 MAIN_CHRS_AVAILABLE)
JOB_SCRIPT_SNP_CHRS="${JOB_SCRIPT_DIR}/merge_snps_chrs_only.sh"
SNP_CHRS_FILES=()
for chr in "${MAIN_CHRS_AVAILABLE[@]}"; do
    SNP_CHRS_FILES+=("${INPUT_DIR}/${chr}.filtered.snps.vcf.gz")
done
INPUT_FILES_SNP_CHRS=$(printf "%s " "${SNP_CHRS_FILES[@]}")
generate_merge_script "MergeSNP_Chrs" "main_chrs_snps.vcf.gz" "${INPUT_FILES_SNP_CHRS}" "${JOB_SCRIPT_SNP_CHRS}" "merge_snps_chrs"
echo "已生成 $JOB_SCRIPT_SNP_CHRS"

# 4. 仅合并主染色体 INDELs
JOB_SCRIPT_INDEL_CHRS="${JOB_SCRIPT_DIR}/merge_indels_chrs_only.sh"
INDEL_CHRS_FILES=()
for chr in "${MAIN_CHRS_AVAILABLE[@]}"; do
    INDEL_CHRS_FILES+=("${INPUT_DIR}/${chr}.filtered.indels.vcf.gz")
done
INPUT_FILES_INDEL_CHRS=$(printf "%s " "${INDEL_CHRS_FILES[@]}")
generate_merge_script "MergeINDEL_Chrs" "main_chrs_indels.vcf.gz" "${INPUT_FILES_INDEL_CHRS}" "${JOB_SCRIPT_INDEL_CHRS}" "merge_indels_chrs"
echo "已生成 $JOB_SCRIPT_INDEL_CHRS"

echo "---- 完成 ----"
echo "已生成针对 ${PROJECT_NAME} 的作业脚本。"