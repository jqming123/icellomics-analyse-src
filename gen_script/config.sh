#!/bin/bash

# ==================================
# Genome Project Configuration
# ==================================

# --- Project Info ---
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects"
REF_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
SRC_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/gen_script"

# PROJECT_NAME 现在期望在 sourcing 此脚本之前作为环境变量被设置。
# 例如，在主脚本中通过 `export PROJECT_NAME="YourProject"` 设置。
if [ -z "${PROJECT_NAME}" ]; then
    echo "错误: PROJECT_NAME 未设置。请在 sourcing config.sh 之前设置它。" >&2
    exit 1
fi

# --- Paths ---
PROJECT_DIR="${BASE_DIR}/${PROJECT_NAME}"
SAMPLE_LIST="${PROJECT_DIR}/00_data/sample_list.txt"
RAW_DATA_DIR="${PROJECT_DIR}/00_data/raw_fastq"
RESULTS_DIR="${PROJECT_DIR}/01_results"
LOG_DIR="${PROJECT_DIR}/03_logs"

# --- 参考基因组配置 ---
# 这里要改成，在生成脚本的文件里指定参考基因组名称
# REF_NAME="CriGri-PICRH-1.0"
# REF_NAME="hg38"

# 1. 检查 REF_NAME 是否已设置
if [ -z "${REF_NAME}" ]; then
    echo "错误: REF_NAME 未设置。请在 sourcing config.sh 之前设置它 (例如: export REF_NAME=\"hg38\")。" >&2
    exit 1
fi

# 2. 根据 REF_NAME 设置参考基因组相关文件的具体路径
case "${REF_NAME}" in
    "CriGri-PICRH-1.0")
        ## CHO NCBI 基因组配置 
        REF_GENOME="${REF_DIR}/CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.fai"
        ## 染色体名称 (NCBI RefSeq accession), 使用.fai文件中的名称
        main_chrs=(
            NC_048595.1
            NC_048596.1
            NC_048597.1
            NC_048598.1
            NC_048599.1
            NC_048600.1
            NC_048601.1
            NC_048602.1
            NC_048603.1
            NC_048604.1
        )
        VEP_SPECIES="cricetulus_griseus_picr"
        ## 基因组版本名称 (必须与VEP缓存中的文件夹名称匹配)
        export GENOME_ASSEMBLY="CriGri-PICRH-1.0"
        ;;

    "CH_Ensemble")
        ## CHO Ensemble基因组配置 
        REF_GENOME="${REF_DIR}/CriGri-PICRH-1.0_Ensemble/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/CriGri-PICRH-1.0_Ensemble/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa.fai"
        ## 染色体名称, 使用.fai文件中的名称
        main_chrs=(
            2
            3
            4
            5
            6
            7
            8
            9
            10
            X
        )
        VEP_SPECIES="cricetulus_griseus_picr"
        ## 基因组版本名称 (必须与VEP缓存中的文件夹名称匹配)
        export GENOME_ASSEMBLY="CriGri-PICRH-1.0"
        ;;

    "hg38_Ensemble")
        ## 人源细胞系 (Ensembl hg38) 配置
        REF_GENOME="${REF_DIR}/hg38_Ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        REF_FAI="${REF_DIR}//hg38_Ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.fai"
        main_chrs=(
            1
            2
            3
            4
            5
            6
            7
            8
            9
            10
            11
            12
            13
            14
            15
            16
            17
            18
            19
            20
            21
            22
            X
            Y
            MT
        )
        ## 这里填的名称后面要核实一下
        VEP_SPECIES="homo_sapiens"
        ## 基因组版本名称 (必须与VEP缓存中的文件夹名称匹配)
        export GENOME_ASSEMBLY="GRCh38"
        ;;
    *)
        # 3. 兜底处理：如果输入的 REF_NAME 不在上述列表中，报错退出
        echo "错误: 未识别的基因组名称 '${REF_NAME}'。" >&2
        echo "当前支持的选项有: CriGri-PICRH-1.0, hg38" >&2
        exit 1
        ;;
esac

# --- Software ---
CONDA_PROFILE_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
# 除VEP以外的软件全部安装在mamba环境中，运行各个脚本前都要先激活下面这个mamba环境
GENOME_ENV_NAME="genome_env"
# VEP单独用一个环境
VEP_ENV_NAME="vep_115"


# --- Resources ---
MEM_XLARGE="128G" # VEP注释任务申请的内存，建议100G或更高
MEM_LARGE="64G"
MEM_LARGE_M4="60G"
MEM_MEDIUM="40G"
MEM_SMALL="20G"

# VEP 是一个资源密集型工具，建议为其分配独立的、更大的资源
THREADS_VEP=16      # VEP注释可以使用的线程数
THREADS=8

# 可选队列
# QUEUE_NAME="vmcore128"
QUEUE_NAME="corexd192"
# QUEUE_NAME="core56"

# ==================================
# VEP Annotation Configuration
# ==================================

# --- VEP Settings ---
# VEP 缓存目录的绝对路径。这是VEP离线运行所必需的。
export VEP_CACHE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/vep_cache"

