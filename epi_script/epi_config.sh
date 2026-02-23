#!/bin/bash

# ==================================
# Epigenetics Project Configuration
# ==================================

# --- Project Info ---
# PROJECT_NAME 和 REF_NAME 期望在 sourcing 此脚本之前作为环境变量被设置。
if [ -z "${PROJECT_NAME}" ]; then
    echo "错误: PROJECT_NAME 未设置。请在 sourcing epi_config.sh 之前设置它。" >&2
    exit 1
fi

if [ -z "${REF_NAME}" ]; then
    echo "错误: REF_NAME 未设置。请在 sourcing epi_config.sh 之前设置它 (例如: export REF_NAME=\"hg38_Ensemble\")。" >&2
    exit 1
fi

# --- Base Directories ---
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects"
RESOURCES_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources"

# --- Project-specific Paths ---
PROJECT_DIR="${BASE_DIR}/${PROJECT_NAME}"
export TMP_DIR="${PROJECT_DIR}/tmp" 

# --- Reference Genome Configuration ---
# 根据 REF_NAME 设置路径及基因组参数 (如 MACS2 用的 GSIZE)
case "${REF_NAME}" in
    "CriGri-PICRH-1.0") ## CHO 仓鼠基因组配置 (已弃用)
        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0"
        export BOWTIE2_INDEX="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic"
        export REF_GENOME="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
        # MACS2 genome size for CHO
        export GSIZE="2366634374"
        ;;

    "CH_Ensemble")  ## CHO 仓鼠基因组配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0_Ensemble"
        export BOWTIE2_INDEX="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
        # MACS2 genome size for CHO
        export GSIZE="2366634374"
        ;;

    "hg38_Ensemble")
        ## 人源细胞系 (Ensembl hg38) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/hg38_Ensemble"
        export BOWTIE2_INDEX="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly"
        export REF_GENOME="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        # MACS2 genome size: 'hs' is shortcut for 2.7e9 (human)
        export GSIZE="hs"
        ;;

    *)
        echo "错误: 未识别的基因组名称 '${REF_NAME}'。" >&2
        echo "当前支持的选项有: CriGri-PICRH-1.0, CH_Ensemble, hg38_Ensemble" >&2
        exit 1
        ;;
esac

# --- Software & Environment ---
CONDA_PROFILE_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
EPI_CONDA_ENV_NAME="ATAC_E4"

# --- Analysis Scripts ---
EPI_SCRIPT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script"

# --- SLURM Resource Configuration ---
THREADS=8
MEM_SUPERLARGE="100G"
MEM_LARGE="70G"
MEM_MEDIUM="40G"
QUEUE_NAME="corexd192"