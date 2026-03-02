#!/usr/bin/env bash

# ==================================
# scRNA Project Shared Configuration
# ==================================

# --- Base layout ---
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine"
PROJECT_ROOT="$BASE_DIR/scRNA_projects"
RESOURCES_ROOT="$BASE_DIR/resources"
GEN_TOOLKIT_DIR="$RESOURCES_ROOT/src/scRNA_script/GENtoolkit"

# --- GENtoolkit scripts ---
UP_PY="$GEN_TOOLKIT_DIR/GENToolkit_alter.py"
DOWN_PY="$GEN_TOOLKIT_DIR/GENToolkit_alter.py"

# 其他软件需要的路径
STAR_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/envs/scRNA_env/bin"


# --- Conda/Mamba ---
CONDA_SH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="scRNA_env"

# --- SLURM defaults ---
UP_CPUS="12"
UP_MEM="64G"
UP_TIME="48:00:00"
DOWN_CPUS="8"
DOWN_MEM="32G"
DOWN_TIME="24:00:00"

# 可选队列
# QUEUE_NAME="vmcore128"
QUEUE_NAME="corexd192"
# QUEUE_NAME="core56"



# --- Downstream defaults ---
WORK_PATH_DEFAULT_SUFFIX="3_expression_result/downstream"

# --- Reference genome mapping ---
# Expect REF_NAME to be set before sourcing this file.
# Example: export REF_NAME="hg38_Ensemble"
if [[ -z "${REF_NAME:-}" ]]; then
  echo "错误: REF_NAME 未设置。请在 sourcing config.sh 之前设置它 (例如: export REF_NAME=\"hg38_Ensemble\")." >&2
  exit 1
fi

REF_DIR="$RESOURCES_ROOT/ref_genome"

case "${REF_NAME}" in
  "CriGri-PICRH-1.0_Ensemble"|"CH_Ensemble")
    REF_GENOME="$REF_DIR/CriGri-PICRH-1.0_Ensemble/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
    REF_GTF="$REF_DIR/CriGri-PICRH-1.0_Ensemble/Cricetulus_griseus_picr.CriGri-PICRH-1.0.115.gtf"
    REF_BED="$REF_DIR/CriGri-PICRH-1.0_Ensemble/CriGri-PICRH-1.0_Ensemble.115.for_scRNA.ref.bed"
    HISAT2_INDEX="$REF_DIR/CriGri-PICRH-1.0_Ensemble/hisat2_index/genome"
    RSEM_INDEX="$REF_DIR/CriGri-PICRH-1.0_Ensemble/rsem.index/reference"
    STAR_INDEX="$REF_DIR/CriGri-PICRH-1.0_Ensemble/star.index"
    CELLRANGER_INDEX="$REF_DIR/CriGri-PICRH-1.0_Ensemble/cellranger_index"
    ;;

  "hg38_Ensemble")
    REF_GENOME="$REF_DIR/hg38_Ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
    REF_GTF="$REF_DIR/hg38_Ensemble/Homo_sapiens.GRCh38.115.gtf"
    REF_BED="$REF_DIR/hg38_Ensemble/Homo_sapiens.GRCh38.115.for_scRNA.ref.bed"
    HISAT2_INDEX="$REF_DIR/hg38_Ensemble/hisat2_index/genome"
    RSEM_INDEX="$REF_DIR/hg38_Ensemble/rsem.index/reference"
    STAR_INDEX="$REF_DIR/hg38_Ensemble/star.index"
    CELLRANGER_INDEX="$REF_DIR/hg38_Ensemble/cellranger_index"
    ;;

  *)
    echo "错误: 未识别的基因组名称 '${REF_NAME}'." >&2
    echo "当前支持: CriGri-PICRH-1.0, CriGri-PICRH-1.0_Ensemble, hg38_Ensemble" >&2
    exit 1
    ;;
esac

# Export commonly used reference variables
export REF_GENOME REF_GTF REF_BED HISAT2_INDEX RSEM_INDEX STAR_INDEX CELLRANGER_INDEX
