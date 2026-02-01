#!/bin/bash
#PBS -N SRA2FASTQ
#PBS -q core40
#PBS -l walltime=100:00:00,nodes=1:ppn=8,mem=50gb
#HSCHED -s hschedd
#PBS -o /gpfs/zhaowm_group/gaoxiaojing/CellLine/genome_analysis/02_results/logs/sra_to_fastq.out
#PBS -e /gpfs/zhaowm_group/gaoxiaojing/CellLine/genome_analysis/02_results/logs/sra_to_fastq.err

echo "Start time:" && date

# 加载配置文件
source /gpfs/zhaowm_group/gaoxiaojing/CellLine/genome_analysis/01_scripts/config.sh

# 激活环境
source /gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/etc/profile.d/conda.sh
conda activate genome_env

# 输入输出目录
SRA_DIR=${PROJECT_DIR}/00_data/dna_sra
FASTQ_DIR=${PROJECT_DIR}/00_data/raw_fastq
mkdir -p ${FASTQ_DIR}

# 使用 fasterq-dump 转换 SRA -> FASTQ
# 遍历所有项目编号下的 .sra 文件
find ${SRA_DIR} -name "*.sra" | while read sra_file; do
    sample=$(basename "${sra_file}" .sra)
    echo "Processing ${sample} ..."

    # fasterq-dump 默认输出 fastq
    fasterq-dump \
        --split-files \
        --threads ${THREADS} \
        --outdir ${FASTQ_DIR} \
        "${sra_file}"

    # 压缩 fastq
    gzip -f ${FASTQ_DIR}/${sample}_1.fastq
    gzip -f ${FASTQ_DIR}/${sample}_2.fastq

    # 重命名为 .fq.gz 以保持流程一致
    mv ${FASTQ_DIR}/${sample}_1.fastq.gz ${FASTQ_DIR}/${sample}_r1.fq.gz
    mv ${FASTQ_DIR}/${sample}_2.fastq.gz ${FASTQ_DIR}/${sample}_r2.fq.gz
done

echo "All SRA files have been converted to FASTQ."
echo "End time:" && date

