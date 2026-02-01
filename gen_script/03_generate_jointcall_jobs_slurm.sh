#!/bin/bash

# 脚本功能: 为每个染色体和未定位的scaffolds生成一个 Slurm 任务脚本，用于联合基因分型 (Joint Calling) 和变异硬过滤 (Hard Filtering)
# 运行方式: 在 resources/src/ 目录下执行 `bash 03_generate_jointcall_jobs_slurm.sh <PROJECT_NAME> <REF_NAME>`

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

# 引入项目配置文件
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "${SCRIPT_DIR}/config.sh" ]; then
    source "${SCRIPT_DIR}/config.sh"
else
    echo "错误: 找不到配置文件 ${SCRIPT_DIR}/config.sh" >&2
    exit 1
fi

# --- 路径定义  ---
# 输入: 上一步生成的 GenomicsDB 数据库目录
DB_DIR="${RESULTS_DIR}/04_genomicsdb"
# 输出: 定义存放原始 VCF 和过滤后 VCF 的目录
VCF_RAW_DIR="${RESULTS_DIR}/05_vcf_raw"
VCF_FILTERED_DIR="${RESULTS_DIR}/06_vcf_filtered"

# 脚本生成目录
JOB_SCRIPT_DIR="${PROJECT_DIR}/02_jobs/joint_calling_jobs"

# 创建该阶段所需的输出和日志目录
# 同时为日志文件创建存放目录
mkdir -p "${VCF_RAW_DIR}" "${VCF_FILTERED_DIR}" "${LOG_DIR}/02_joint_calling/err_file" "${LOG_DIR}/02_joint_calling/out_file" "${JOB_SCRIPT_DIR}"

# --- 【注意】要并行处理的染色体列表改为在config.sh中设置 ---
# 定义主要的染色体 (NCBI RefSeq accession)
# main_chrs=(
#     NC_048595.1
#     NC_048596.1
#     NC_048597.1
#     NC_048598.1
#     NC_048599.1
#     NC_048600.1
#     NC_048601.1
#     NC_048602.1
#     NC_048603.1
#     NC_048604.1
# )

# 添加未定位scaffolds标记
# regions_to_process=("${main_chrs[@]}" "unplaced_scaffolds")

echo "Generating VCF joint calling and filtering jobs for each region..."

for region in "${main_chrs[@]}"; do
    JOB_SCRIPT_PATH="${JOB_SCRIPT_DIR}/${region}_call.sh"
    
    # 使用 here-document 将多行内容写入 Slurm 脚本
    cat <<EOF > "${JOB_SCRIPT_PATH}"
#!/bin/bash
#SBATCH -p ${QUEUE_NAME}                  # 指定作业提交的分区 (队列)
#SBATCH -J CALL_${region}_${PROJECT_NAME}                 # 作业名称 (与上一个脚本命名风格一致)
#SBATCH -o ${LOG_DIR}/02_joint_calling/${region}_call_%j.log  # 标准输出重定向 (修正：.out文件通常用于标准输出)
#SBATCH --nodes=1                         # 作业申请 1 个节点
#SBATCH --ntasks-per-node=1               # 单节点启动 1 个任务
#SBATCH --cpus-per-task=${THREADS}        # 使用配置的线程数
#SBATCH --mem=${MEM_LARGE}                # 申请的内存大小
#SBATCH --time=240:00:00                # 任务运行最长时间 (与上一个脚本一致)

set -eo pipefail

echo "Job started on: \$(hostname)"
echo "Start time: " && date

# 激活 mamba 环境 (与上一个脚本一致)
source "${CONDA_PROFILE_PATH}"
conda activate "${GENOME_ENV_NAME}"

# Step 1: Joint Genotyping
gatk --java-options "-Xmx${MEM_LARGE_M4} -Xms${MEM_MEDIUM}" GenotypeGVCFs \\
    -R "${REF_GENOME}" \\
    -V "gendb://${DB_DIR}/${region}" \\
    -O "${VCF_RAW_DIR}/${region}.raw.vcf.gz"

# Step 2: Split SNPs and INDELs
gatk SelectVariants -R "${REF_GENOME}" -V "${VCF_RAW_DIR}/${region}.raw.vcf.gz" \\
    --select-type-to-include SNP -O "${VCF_RAW_DIR}/${region}.raw.snps.vcf.gz"

gatk SelectVariants -R "${REF_GENOME}" -V "${VCF_RAW_DIR}/${region}.raw.vcf.gz" \\
    --select-type-to-include INDEL -O "${VCF_RAW_DIR}/${region}.raw.indels.vcf.gz"

# ===================================================================
# 添加：统计缺失字段的变异数量
# ===================================================================
echo "=== 统计缺失过滤字段的变异数量 ==="

# 检查SNP文件是否存在
if [[ ! -f "${VCF_RAW_DIR}/${region}.raw.snps.vcf.gz" ]]; then
    echo "错误: SNP文件不存在 - ${VCF_RAW_DIR}/${region}.raw.snps.vcf.gz"
    exit 1
fi

# 统计SNP文件中缺失字段的变异
echo "统计SNP文件缺失字段情况:"
for field in QD MQ FS SOR; do
    # 修正: 转义$, 使得count和field变量在子脚本执行时被解析
    count=\$(zcat "${VCF_RAW_DIR}/${region}.raw.snps.vcf.gz" 2>/dev/null | 
            awk -v field="\$field" '
                !/^#/ && !match(\$0, field"=") {count++} 
                END {print count+0}'
            ) || count="ERROR"
    echo " - 缺少 \$field 字段的SNP数量: \$count"
done

# 检查INDEL文件是否存在
if [[ ! -f "${VCF_RAW_DIR}/${region}.raw.indels.vcf.gz" ]]; then
    echo "错误: INDEL文件不存在 - ${VCF_RAW_DIR}/${region}.raw.indels.vcf.gz"
    exit 1
fi

# 统计INDEL文件中缺失字段的变异
echo "统计INDEL文件缺失字段情况:"
for field in QD FS SOR; do
    # 修正: 转义$, 使得count和field变量在子脚本执行时被解析
    count=\$(zcat "${VCF_RAW_DIR}/${region}.raw.indels.vcf.gz" 2>/dev/null | 
            awk -v field="\$field" '
                !/^#/ && !match(\$0, field"=") {count++} 
                END {print count+0}'
            ) || count="ERROR"
    echo " - 缺少 \$field 字段的INDEL数量: \$count"
done

# 统计完全缺失所有关键字段的变异
echo "统计完全缺失所有关键过滤字段的变异:"
# SNP
# 修正: 转义$
count_snp=\$(zcat "${VCF_RAW_DIR}/${region}.raw.snps.vcf.gz" 2>/dev/null | 
            awk '!/^#/ && !/QD=/ && !/MQ=/ && !/FS=/ && !/SOR=/ {count++} 
                END {print count+0}'
            ) || count_snp="ERROR"
echo " - 完全缺失所有字段的SNP数量: \$count_snp"

# INDEL
# 修正: 转义$
count_indel=\$(zcat "${VCF_RAW_DIR}/${region}.raw.indels.vcf.gz" 2>/dev/null | 
              awk '!/^#/ && !/QD=/ && !/FS=/ && !/SOR=/ {count++} 
                  END {print count+0}'
              ) || count_indel="ERROR"
echo " - 完全缺失所有字段的INDEL数量: \$count_indel"
echo "======================================"


# Step 3: Filter SNPs
gatk VariantFiltration -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${region}.raw.snps.vcf.gz" \\
    -O "${VCF_RAW_DIR}/${region}.raw.snps.tagged.vcf.gz" \\
    --filter-expression "QD < 2.0" --filter-name "QD2" \\
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
    --filter-expression "FS > 60.0" --filter-name "FS60" \\
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \\
    --missing-values-evaluate-as-failing true

gatk SelectVariants -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${region}.raw.snps.tagged.vcf.gz" \\
    --exclude-filtered -O "${VCF_FILTERED_DIR}/${region}.filtered.snps.vcf.gz"

# Step 4: Filter INDELs 
gatk VariantFiltration -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${region}.raw.indels.vcf.gz" \\
    -O "${VCF_RAW_DIR}/${region}.raw.indels.tagged.vcf.gz" \\
    --filter-expression "QD < 2.0" --filter-name "QD2" \\
    --filter-expression "FS > 200.0" --filter-name "FS200" \\
    --filter-expression "SOR > 10.0" --filter-name "SOR10" \\
    --missing-values-evaluate-as-failing true

gatk SelectVariants -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${region}.raw.indels.tagged.vcf.gz" \\
    --exclude-filtered -O "${VCF_FILTERED_DIR}/${region}.filtered.indels.vcf.gz"

echo "End time: " && date
echo "Job finished."
EOF

done

echo "Done. Slurm scripts for joint calling are generated in ${JOB_SCRIPT_DIR}"