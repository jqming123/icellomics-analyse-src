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
# 根据 REF_NAME 设置路径及基因组参数 (如 MACS3 用的 GSIZE)
case "${REF_NAME}" in
#    "CriGri-PICRH-1.0") ## CHO 仓鼠基因组配置 (已弃用)
#        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0"
#        export BOWTIE2_INDEX="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic"
#        export REF_GENOME="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
#        # MACS3 effective genome size for CHO
#        export GSIZE="2366634374"
#        ;;

    "CH_Ensembl")  ## CHO 仓鼠基因组配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0_Ensembl"
        export BOWTIE2_INDEX="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
        # MACS3 effective genome size for CHO
        export GSIZE="2366634374"
        ;;

    "hg38_Ensembl")
        ## 人源细胞系 (Ensembl hg38) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/hg38_Ensembl"
        export BOWTIE2_INDEX="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly"
        export REF_GENOME="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        # MACS3 effective genome size: 'hs' is shortcut for 2,913,022,398 (GRCh38)
        export GSIZE="hs"
        ;;
    "Cattle_E_ARSUCD2")
        ## 牛 (Bos taurus, ARS-UCD2.0) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/Cattle_E_ARSUCD2"
        export BOWTIE2_INDEX="${REF_DIR}/Bos_taurus.ARS-UCD2.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa"
        # MACS3 effective genome size 
        export GSIZE="2770686120"
        ;;
    "Chicken_E_GRCg7b")
        ## 鸡 (Gallus gallus, GRCg7b) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/Chicken_E_GRCg7b"
        export BOWTIE2_INDEX="${REF_DIR}/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
        # MACS3 effective genome size (待填写)
        export GSIZE="1053332251"
        ;;
    "Dog_E_UUGSD")
        ## 狗 (Canis lupus familiaris, UU_Cfam_GSD_1.0) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/Dog_E_UUGSD"
        export BOWTIE2_INDEX="${REF_DIR}/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
        # MACS3 effective genome size (待填写)
        export GSIZE=""
        ;;
    "GreenMonkey_E_ChlSab1.1")
        ## 绿猴 (Chlorocebus sabaeus, ChlSab1.1) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/GreenMonkey_E_ChlSab1.1"
        export BOWTIE2_INDEX="${REF_DIR}/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa"
        # MACS3 effective genome size (待填写)
        export GSIZE=""
        ;;
    "Mouse_E_GRCm39")
        ## 小鼠 (Mus musculus, GRCm39) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/Mouse_E_GRCm39"
        export BOWTIE2_INDEX="${REF_DIR}/Mus_musculus.GRCm39.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Mus_musculus.GRCm39.dna.toplevel.fa"
        # MACS3 effective genome size (待填写，这个印象中有缩写可用)
        export GSIZE=""
        ;;
    "Pig_E_Sscrofa11.1")
        ## 猪 (Sus scrofa, Sscrofa11.1) 配置
        REF_DIR="${RESOURCES_DIR}/ref_genome/Pig_E_Sscrofa11.1"
        export BOWTIE2_INDEX="${REF_DIR}/Sus_scrofa.Sscrofa11.1.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
        # MACS3 effective genome size (待填写)
        export GSIZE=""
        ;;
    "dont_need_ref")
        # 用于不需要参考基因组的任务（如 SRA 转 FASTQ）
        export BOWTIE2_INDEX="NONE"
        export REF_GENOME="NONE"
        export GSIZE="0"
        ;;

    *)
        echo "错误: 未识别的基因组名称 '${REF_NAME}'。" >&2
        echo "当前支持的选项有: CriGri-PICRH-1.0, CH_Ensemble, hg38_Ensemble, Cattle_E_ARSUCD2, Chicken_E_GRCg7b, Dog_E_UUGSD, GreenMonkey_E_ChlSab1.1, Mouse_E_GRCm39, Pig_E_Sscrofa11.1" >&2
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
MEM_SUPERLARGE="120G"
MEM_LARGE="64G"
MEM_MEDIUM="32G"
MEM_SMALL="16G"

# QUEUE_NAME="corexd192"
QUEUE_NAME="core56"
# QUEUE_NAME="vmcore128"
