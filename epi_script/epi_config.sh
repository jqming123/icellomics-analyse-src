#!/bin/bash

# ==================================
# Epigenetics Project Configuration
# ==================================

# --- Project Info ---
# PROJECT_NAME 和 REF_NAME 期望在 sourcing 此脚本之前作为环境变量被设置。
# 例如: export PROJECT_NAME="YourProject"; export REF_NAME="hg38_Ensemble"

# 1. 检查 PROJECT_NAME 是否已设置
if [ -z "${PROJECT_NAME}" ]; then
    echo "错误: PROJECT_NAME 未设置。请在 sourcing epi_config.sh 之前设置它。" >&2
    exit 1
fi

# 2. 检查 REF_NAME 是否已设置
if [ -z "${REF_NAME}" ]; then
    echo "错误: REF_NAME 未设置。请在 sourcing epi_config.sh 之前设置它 (例如: export REF_NAME=\"hg38_Ensemble\")。" >&2
    exit 1
fi

# --- Base Directories ---
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects"
RESOURCES_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources"

# --- Project-specific Paths ---
PROJECT_DIR="${BASE_DIR}/${PROJECT_NAME}"
export TMP_DIR="${PROJECT_DIR}/tmp" # 临时文件目录

# --- Reference Genome Configuration ---
# 根据 REF_NAME 分支设置具体路径
case "${REF_NAME}" in
    "CriGri-PICRH-1.0")
        ## CHO 仓鼠基因组配置 
        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0"
        # Bowtie2 索引的基础名
        export BOWTIE2_INDEX="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic"
        # 如果后续需要用到 FASTA 路径，可以参考下面的写法：
        REF_GENOME="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
        ;;

    "CH_Ensemble")
        ## CHO 仓鼠基因组配置 
        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0_Ensemble"
        # Bowtie2 索引的基础名
        export BOWTIE2_INDEX="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel"
        # 如果后续需要用到 FASTA 路径，可以参考下面的写法：
        REF_GENOME="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
        ;;

    "hg38_Ensemble")
        ## 人源细胞系 (Ensembl hg38) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/hg38_Ensemble"
        # Bowtie2 索引的基础名 (根据参考代码中的 fasta 文件名推断)
        export BOWTIE2_INDEX="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly"
        REF_GENOME="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        ;;

    *)
        # 如果输入的 REF_NAME 不在上述列表中，报错退出
        echo "错误: 未识别的基因组名称 '${REF_NAME}'。" >&2
        echo "当前支持的选项有: CriGri-PICRH-1.0, CH_Ensemble, hg38_Ensemble" >&2
        exit 1
        ;;
esac

# --- Software & Environment ---
# Conda 环境配置
CONDA_PROFILE_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
EPI_CONDA_ENV_NAME="ATAC_E4" # 表观遗传学分析专用环境

# --- Analysis Scripts ---
EPI_SCRIPT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script"

# --- SLURM Resource Configuration ---
THREADS=8
MEM_SUPERLARGE="100G"
MEM_LARGE="70G"
MEM_MEDIUM="40G"
QUEUE_NAME="corexd192"

