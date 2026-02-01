#!/bin/bash

# ==============================================================================
# 脚本功能: 为每个 GenomicsDB 分区生成一个 PBS 任务脚本，用于执行以下步骤：
#           1. 联合基因分型 (GenotypeGVCFs)
#           2. 变异硬过滤 (VariantFiltration)
# 前置条件: 必须已成功完成所有 `GenomicsDBImport` 任务。
# 运行方式: 在 01_scripts/ 目录下执行 `bash 03_generate_jointcall_jobs_by_scatter.sh`
# 脚本03_generate_jointcall_jobs_by_scatter.sh是脚本03_generate_jointcall_jobs.sh的修改版，主要是把 **“一个任务处理一个区域（染色体）”** 改成了 **“一个任务处理一批区域（多个scaffolds）”
# ==============================================================================

# 引入项目配置文件
source ./config.sh

# --- 路径定义 (自动推导) ---
# 输入: 上一步生成的 GenomicsDB 数据库目录
DB_DIR="${RESULTS_DIR}/04_genomicsdb/intervals_by_size"
# 输出: 定义存放原始 VCF 和过滤后 VCF 的目录
VCF_RAW_DIR="${RESULTS_DIR}/05_vcf_raw"
VCF_FILTERED_DIR="${RESULTS_DIR}/06_vcf_filtered"

# --- 主逻辑 ---

# 1. 创建该阶段所需的输出和日志目录
mkdir -p "$VCF_RAW_DIR" "$VCF_FILTERED_DIR" "$LOG_DIR/02_joint_calling"

# 2. 检查 GenomicsDB 数据库分区是否存在
#    我们通过查找名为 'scatter_*' 的目录来确定
db_partitions=(${DB_DIR}/scatter_*)
if [ ! -d "${db_partitions[0]}" ]; then
    echo "Error: Could not find any GenomicsDB partitions in '${DB_DIR}'."
    echo "Please ensure the previous GenomicsDBImport step has completed successfully."
    exit 1
fi

num_partitions=$(ls -d ${DB_DIR}/scatter_* | wc -l)
echo "Found ${num_partitions} GenomicsDB partitions. Generating a PBS job for each..."

# 3. 遍历所有 GenomicsDB 分区，生成对应的 PBS 任务脚本
for db_workspace in ${DB_DIR}/scatter_*; do
    # 从数据库路径中提取分区的ID，例如 'scatter_001'
    partition_id=$(basename "$db_workspace")
    
    JOB_SCRIPT_NAME="${partition_id}_call.sh"

    # 使用 here-document 将多行内容写入 PBS 脚本
    cat <<EOF > "${JOB_SCRIPT_NAME}"
#PBS -q core40
#PBS -l walltime=1001:00:00,nodes=1:ppn=${THREADS},mem=${MEM_LARGE}
#HSCHED -s hschedd
#PBS -o ${LOG_DIR}/02_joint_calling/${partition_id}_call.out
#PBS -e ${LOG_DIR}/02_joint_calling/${partition_id}_call.err
#PBS -N VCF_${partition_id}

LOGFILE="${LOG_DIR}/02_joint_calling/${partition_id}_call.log"
touch "\${LOGFILE}"
exec > "\${LOGFILE}" 2>&1

echo "Start time: " && date
echo "Running Joint Calling and Filtering for partition: ${partition_id}"

# 激活 mamba 环境
source /gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/etc/profile.d/conda.sh
conda activate ${MAMBA_ENV_NAME}

# Step 1: Joint Genotyping for the current partition
# 输入是当前分区的 GenomicsDB 路径
gatk --java-options "-Xmx${MEM_MEDIUM}" GenotypeGVCFs \\
    -R "${REF_GENOME}" \\
    -V "gendb://${db_workspace}" \\
    -O "${VCF_RAW_DIR}/${partition_id}.raw.vcf.gz"

# Step 2: Split SNPs and INDELs
gatk SelectVariants -R "${REF_GENOME}" -V "${VCF_RAW_DIR}/${partition_id}.raw.vcf.gz" \\
    --select-type-to-include SNP -O "${VCF_RAW_DIR}/${partition_id}.raw.snps.vcf.gz"

gatk SelectVariants -R "${REF_GENOME}" -V "${VCF_RAW_DIR}/${partition_id}.raw.vcf.gz" \\
    --select-type-to-include INDEL -O "${VCF_RAW_DIR}/${partition_id}.raw.indels.vcf.gz"

# Step 3: Filter SNPs
gatk VariantFiltration -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${partition_id}.raw.snps.vcf.gz" \\
    -O "${VCF_RAW_DIR}/${partition_id}.raw.snps.tagged.vcf.gz" \\
    --filter-expression "QD < 2.0" --filter-name "QD2" \\
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \\
    --filter-expression "FS > 60.0" --filter-name "FS60" \\
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

gatk SelectVariants -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${partition_id}.raw.snps.tagged.vcf.gz" \\
    --exclude-filtered -O "${VCF_FILTERED_DIR}/${partition_id}.filtered.snps.vcf.gz"

# Step 4: Filter INDELs
gatk VariantFiltration -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${partition_id}.raw.indels.vcf.gz" \\
    -O "${VCF_RAW_DIR}/${partition_id}.raw.indels.tagged.vcf.gz" \\
    --filter-expression "QD < 2.0" --filter-name "QD2" \\
    --filter-expression "FS > 200.0" --filter-name "FS200" \\
    --filter-expression "SOR > 10.0" --filter-name "SOR10"

gatk SelectVariants -R "${REF_GENOME}" \\
    -V "${VCF_RAW_DIR}/${partition_id}.raw.indels.tagged.vcf.gz" \\
    --exclude-filtered -O "${VCF_FILTERED_DIR}/${partition_id}.filtered.indels.vcf.gz"

echo "End time: " && date
EOF

done

echo "Done. ${num_partitions} PBS scripts for VCF calling and filtering are generated in the current directory."