#!/bin/bash

# 脚本功能: 为每个染色体生成一个 PBS 任务脚本，用于联合基因分型 (Joint Calling) 和变异硬过滤 (Hard Filtering)
# 运行方式: 在 01_scripts/ 目录下执行 `bash 03_generate_jointcall_jobs.sh`

# 引入项目配置文件
source ./config.sh


# --- 路径定义 (自动推导) ---
# 输入: 上一步生成的 GenomicsDB 数据库目录
DB_DIR="${RESULTS_DIR}/04_genomicsdb"
# 输出: 定义存放原始 VCF 和过滤后 VCF 的目录
VCF_RAW_DIR="${RESULTS_DIR}/05_vcf_raw"
VCF_FILTERED_DIR="${RESULTS_DIR}/06_vcf_filtered"
# 创建该阶段所需的输出和日志目录
mkdir -p "$VCF_RAW_DIR" "$VCF_FILTERED_DIR" "$LOG_DIR/02_joint_calling"

# 定义要并行处理的染色体列表
Run=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

echo "Generating VCF joint calling and filtering jobs for each chromosome..."

for chromosome in "${Run[@]}"; do
    # 为每条染色体生成一个独立的 PBS 脚本
    JOB_SCRIPT_NAME="${chromosome}_call.sh"

    # 使用 here-document 将多行内容写入 PBS 脚本
    cat <<EOF > "${JOB_SCRIPT_NAME}"
#PBS -q core40
#PBS -l walltime=1001:00:00,nodes=1:ppn=${THREADS},mem=${MEM_LARGE}
#HSCHED -s hschedd
#PBS -o ${LOG_DIR}/02_joint_calling/${chromosome}_call.out
#PBS -e ${LOG_DIR}/02_joint_calling/${chromosome}_call.err
#PBS -N VCF_${chromosome} # 任务命名

echo 'Start time: ' && date

# 激活 genome_env 环境
source /gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/etc/profile.d/conda.sh
mamba activate genome_env

# Step 1: Joint Genotyping
gatk --java-options "-Xmx${MEM_MEDIUM}" GenotypeGVCFs \
    -R "${REF_GENOME}" \
    -V "gendb://${DB_DIR}/${chromosome}" \
    -O "${VCF_RAW_DIR}/${chromosome}.raw.vcf.gz"

# Step 2: Split SNPs and INDELs
gatk SelectVariants -R "${REF_GENOME}" -V "${VCF_RAW_DIR}/${chromosome}.raw.vcf.gz" \
    --select-type-to-include SNP -O "${VCF_RAW_DIR}/${chromosome}.raw.snps.vcf.gz"

gatk SelectVariants -R "${REF_GENOME}" -V "${VCF_RAW_DIR}/${chromosome}.raw.vcf.gz" \
    --select-type-to-include INDEL -O "${VCF_RAW_DIR}/${chromosome}.raw.indels.vcf.gz"

# Step 3: Filter SNPs
gatk VariantFiltration -R "${REF_GENOME}" \
    -V "${VCF_RAW_DIR}/${chromosome}.raw.snps.vcf.gz" \
    -O "${VCF_RAW_DIR}/${chromosome}.raw.snps.tagged.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3"

gatk SelectVariants -R "${REF_GENOME}" \
    -V "${VCF_RAW_DIR}/${chromosome}.raw.snps.tagged.vcf.gz" \
    --exclude-filtered -O "${VCF_FILTERED_DIR}/${chromosome}.filtered.snps.vcf.gz"

# Step 4: Filter INDELs
gatk VariantFiltration -R "${REF_GENOME}" \
    -V "${VCF_RAW_DIR}/${chromosome}.raw.indels.vcf.gz" \
    -O "${VCF_RAW_DIR}/${chromosome}.raw.indels.tagged.vcf.gz" \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "SOR > 10.0" --filter-name "SOR10"

gatk SelectVariants -R "${REF_GENOME}" \
    -V "${VCF_RAW_DIR}/${chromosome}.raw.indels.tagged.vcf.gz" \
    --exclude-filtered -O "${VCF_FILTERED_DIR}/${chromosome}.filtered.indels.vcf.gz"

echo 'End time: ' && date
EOF

done

echo "Done. PBS scripts for VCF calling and filtering are generated in the current directory."

