#!/bin/bash
# ==============================================================================
# 01_generate_mapping_jobs_slurm.sh
# 生成：1) 每个 run 的 mapping job 脚本（并行提交）；
#      2) 每个 biosample 的 merge + gVCF 调用脚本（合并 run.bam -> merged.bam -> 调用 run_single_sample_processing.sh）
# 假设此脚本放在项目的 resources/src 目录下运行或至少能访问同目录下的 config.sh
# 运行方式: 在 resources/src 目录下执行 `bash 01_generate_mapping_jobs_slurm.sh <PROJECT_NAME> <REF_NAME>`
# ==============================================================================
set -eo pipefail

# --- 获取参数 ---
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "用法: $0 <PROJECT_NAME> <REF_NAME>" >&2
    echo "例子: $0 MyProject hg38" >&2
    exit 1
fi

PROJECT_NAME="$1"
REF_NAME="$2"        # <-- 修改：从第二个参数获取
export PROJECT_NAME
export REF_NAME

echo "当前项目名称 (PROJECT_NAME): ${PROJECT_NAME}"
echo "当前参考基因组 (REF_NAME): ${REF_NAME}"

# 脚本绝对路径
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 引入项目配置（必须在设置好 REF_NAME 之后执行）
if [ -f "${SCRIPT_DIR}/config.sh" ]; then
    source "${SCRIPT_DIR}/config.sh"
else
    echo "错误: 找不到配置文件 ${SCRIPT_DIR}/config.sh" >&2
    exit 1
fi

echo "当前使用的队列 (QUEUE_NAME): ${QUEUE_NAME}"

# 输出目录定义（与你原来保持一致）
BAM_DIR="${RESULTS_DIR}/02_bam"
GVCF_DIR="${RESULTS_DIR}/03_gvcf"
MAP_LOG_DIR="${LOG_DIR}/01_mapping_gvcf"

# 脚本生成目录（分 run_jobs 与 merge_jobs）
JOB_SCRIPT_DIR="${PROJECT_DIR}/02_jobs/generated_mapping_jobs"
RUN_JOB_DIR="${JOB_SCRIPT_DIR}/run_jobs"
MERGE_JOB_DIR="${JOB_SCRIPT_DIR}/merge_jobs"

# 保证目录存在
mkdir -p "${BAM_DIR}" "${GVCF_DIR}" "${MAP_LOG_DIR}" "${RUN_JOB_DIR}" "${MERGE_JOB_DIR}"
mkdir -p "${PROJECT_DIR}/00_data/trimmed"

echo "Generating per-run job scripts in: ${RUN_JOB_DIR}"
echo "Generating per-sample merge scripts in: ${MERGE_JOB_DIR}"

# 读取 SAMPLE_LIST，文件格式假设：每行至少两列： <run_id> <sample_id> ...
# 生成每个 run 的 job 脚本（脚本名用 run_id.sh）
awk '{print $1, $2}' "${SAMPLE_LIST}" | while read -r RUN_ID SAMPLE_ID; do
    # skip blank lines
    [ -z "${RUN_ID}" ] && continue

    # ensure sample dir exists for BAMs
    SAMPLE_BAM_DIR="${BAM_DIR}/${SAMPLE_ID}"
    mkdir -p "${SAMPLE_BAM_DIR}"

    RUN_JOB="${RUN_JOB_DIR}/${RUN_ID}.sh"

    cat > "${RUN_JOB}" <<EOF
#!/bin/bash
#SBATCH -J map_${RUN_ID}                   
#SBATCH -p ${QUEUE_NAME}                   
#SBATCH -t 240:00:00                       
#SBATCH -N 1                               
#SBATCH --ntasks-per-node=1                
#SBATCH --cpus-per-task=${THREADS}         
#SBATCH --mem=${MEM_XLARGE}                 
#SBATCH -o ${LOG_DIR}/01_mapping_gvcf/${SAMPLE_ID}_${RUN_ID}_%j.log    

set -eo pipefail

# 确保在运行节点上加载项目配置与 conda 环境
export PROJECT_NAME="${PROJECT_NAME}" # 再次导出 PROJECT_NAME, 确保子作业环境也包含
export REF_NAME="${REF_NAME}"
SCRIPT_DIR="${SCRIPT_DIR}"
source "\${SCRIPT_DIR}/config.sh"
source "\${CONDA_PROFILE_PATH}"
conda activate "\${GENOME_ENV_NAME}"

RUN_ID="${RUN_ID}"
SAMPLE="${SAMPLE_ID}"

RAW_R1="\${RAW_DATA_DIR}/${RUN_ID}_r1.fq.gz"
RAW_R2="\${RAW_DATA_DIR}/${RUN_ID}_r2.fq.gz"
TRIM_DIR="\${PROJECT_DIR}/00_data/trimmed"
OUT_BAM="\${RESULTS_DIR}/02_bam/\${SAMPLE}/${RUN_ID}.bam"

mkdir -p "\${TRIM_DIR}"

# 如果 OUT_BAM 已经存在且带 index, 则跳过后续步骤
if [ -f "\${OUT_BAM}" ] && [ -f "\${OUT_BAM}.bai" ]; then
    echo "\$(date): \${OUT_BAM} exists and indexed - skipping."
    exit 0
fi

export _JAVA_OPTIONS="-Xmx16G"

# v3.1 fix: 自适应Trimmomatic参数, 检测reads长度避免短reads被100%丢弃
# 原因: HEADCROP:8 + MINLEN:50 要求reads至少58bp, 50bp的reads会被全部丢弃
DETECTED_LEN=\$(zcat "\${RAW_R1}" 2>/dev/null | head -4000 | awk 'NR%4==2{sum+=length(\$0); n++} END{if(n>0) printf "%d", sum/n}' || true)
if [ -z "\${DETECTED_LEN}" ] || [ "\${DETECTED_LEN}" -lt 1 ]; then
    DETECTED_LEN=100
fi
echo "\$(date): Detected average read length for ${RUN_ID}: \${DETECTED_LEN}bp"

if [ "\${DETECTED_LEN}" -le 58 ]; then
    # 短reads (<=58bp): 去掉HEADCROP, MINLEN设为reads长度的70%
    MIN_LEN=\$(( DETECTED_LEN * 7 / 10 ))
    [ "\${MIN_LEN}" -lt 36 ] && MIN_LEN=36
    TRIM_PARAMS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:\${MIN_LEN}"
    echo "\$(date): Using short-read trim params (no HEADCROP, MINLEN:\${MIN_LEN}) for \${DETECTED_LEN}bp reads"
else
    # 正常reads (>58bp): 保持原参数
    TRIM_PARAMS="LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 HEADCROP:8 MINLEN:50"
    echo "\$(date): Using standard trim params (HEADCROP:8, MINLEN:50) for \${DETECTED_LEN}bp reads"
fi

echo "\$(date): Trimming ${RUN_ID}"
trimmomatic PE -threads ${THREADS} \
  "\${RAW_R1}" "\${RAW_R2}" \
  "\${TRIM_DIR}/${RUN_ID}_clean_1.fastq" "\${TRIM_DIR}/${RUN_ID}_single_1.fastq" \
  "\${TRIM_DIR}/${RUN_ID}_clean_2.fastq" "\${TRIM_DIR}/${RUN_ID}_single_2.fastq" \
  \${TRIM_PARAMS}

echo "\$(date): Mapping ${RUN_ID}"
RG="@RG\\tID:${RUN_ID}\\tSM:${SAMPLE_ID}\\tLB:\\tPL:ILLUMINA\\tPU:${RUN_ID}"
bwa mem -M -R "\${RG}" "\${REF_GENOME}" \
    "\${TRIM_DIR}/${RUN_ID}_clean_1.fastq" "\${TRIM_DIR}/${RUN_ID}_clean_2.fastq" \
  | samtools view -b -q 20 -@ ${THREADS} - \
  | samtools sort -@ ${THREADS} -o "\${OUT_BAM}.tmp" -

mv "\${OUT_BAM}.tmp" "\${OUT_BAM}" # 将临时文件重命名为最终文件

# samtools addreplacerg -w \
#   -r "@RG\tID:${ID:-$RUN_ID}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:${PL:-ILLUMINA}\tPU:${PU:-$RUN_ID}" \
#   -O BAM \
#   -o "\${OUT_BAM}" \
#   "\${OUT_BAM}.tmp"

# rm "\${OUT_BAM}.tmp"


samtools index "\${OUT_BAM}"

echo "\$(date): Done ${RUN_ID}"
EOF

    chmod +x "${RUN_JOB}"
done

# 生成每个 sample 的 merge + 后处理脚本（合并该 sample 下所有 run 的 bam，然后调用 run_single_sample_processing.sh）
awk '{print $2}' "${SAMPLE_LIST}" | sort -u | while read -r SAMPLE_ID; do
    MERGE_JOB="${MERGE_JOB_DIR}/${SAMPLE_ID}_merge.sh"

    cat > "${MERGE_JOB}" <<EOF
#!/bin/bash
#SBATCH -J merge_${SAMPLE_ID}              
#SBATCH -p ${QUEUE_NAME}                   
#SBATCH -t 240:00:00                        
#SBATCH -N 1                               
#SBATCH --ntasks-per-node=1                
#SBATCH --cpus-per-task=${THREADS}         
#SBATCH --mem=${MEM_LARGE}                 
#SBATCH -o ${LOG_DIR}/01_mapping_gvcf/${SAMPLE_ID}_merge_%j.log

set -eo pipefail

export PROJECT_NAME="${PROJECT_NAME}" # 再次导出 PROJECT_NAME, 确保子作业环境也包含
export REF_NAME="${REF_NAME}"
SCRIPT_DIR="${SCRIPT_DIR}"
source "\${SCRIPT_DIR}/config.sh"
source "\${CONDA_PROFILE_PATH}"
conda activate "\${GENOME_ENV_NAME}"

export GATK_MARKDUP_XMX="${GATK_MARKDUP_XMX}"
export GATK_HC_XMX="${GATK_HC_XMX}"
export GATK_XMS="${GATK_XMS}"

echo "GATK_MARKDUP_XMX=\${GATK_MARKDUP_XMX}"
echo "GATK_HC_XMX=\${GATK_HC_XMX}"
echo "GATK_XMS=\${GATK_XMS}"

SAMPLE="${SAMPLE_ID}"
SAMPLE_BAM_DIR="\${RESULTS_DIR}/02_bam/\${SAMPLE}"
MERGED_BAM="\${SAMPLE_BAM_DIR}/\${SAMPLE}.merged.bam"
MERGED_BAM_INDEX="\${MERGED_BAM}.bai" # 定义索引文件路径

# 检查合并后的BAM文件及其索引是否存在, 若存在则跳过合并BAM文件及建索引步骤
if [ -f "\${MERGED_BAM}" ] && [ -f "\${MERGED_BAM_INDEX}" ]; then
    echo "\$(date): Merged BAM file and index already exist for \${SAMPLE}. Skipping merge and index steps."
else
    # 收集该样本下所有 run 的 bam(若无则退出)
    run_bams=(\$(ls "\${SAMPLE_BAM_DIR}"/*.bam 2>/dev/null || true))
    if [ \${#run_bams[@]} -eq 0 ]; then
        echo "\$(date): No run BAMs found for \${SAMPLE} in \${SAMPLE_BAM_DIR} - exiting"
        exit 1
    fi

    echo "\$(date): Merging \${#run_bams[@]} BAMs for \${SAMPLE} (all should have SM=\${SAMPLE})"
    samtools merge -@ ${THREADS} "\${MERGED_BAM}" "\${run_bams[@]}"
    samtools index "\${MERGED_BAM}"
    echo "\$(date): Merging and indexing completed for \${SAMPLE}"
fi

echo "\$(date): Calling run_single_sample_processing.sh for \${SAMPLE}"
bash "\${SCRIPT_DIR}/run_single_sample_processing.sh" \
    --sampleName "\${SAMPLE}" \
    --inputBam "\${MERGED_BAM}" \
    --ref "\${REF_GENOME}" \
    --outDir "\${RESULTS_DIR}/03_gvcf" \
    --logDir "${MAP_LOG_DIR}" \
    --projectName "${PROJECT_NAME}" # Add this line to pass PROJECT_NAME

echo "\$(date): Merge job finished for \${SAMPLE}"
EOF

    chmod +x "${MERGE_JOB}"
done

echo "Finished generating scripts."
echo " -> Run jobs at: ${RUN_JOB_DIR}"
echo " -> Merge jobs at: ${MERGE_JOB_DIR}"
