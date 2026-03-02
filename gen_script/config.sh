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
#    "CriGri-PICRH-1.0") 
#        ## CHO NCBI 基因组配置 （已弃用）
#        REF_GENOME="${REF_DIR}/CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
#        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
#        REF_FAI="${REF_DIR}/CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna.fai"
#        ## 染色体名称 (NCBI RefSeq accession), 使用.fai文件中的名称
#        main_chrs=(
#            NC_048595.1
#            NC_048596.1
#            NC_048597.1
#            NC_048598.1
#            NC_048599.1
#            NC_048600.1
#            NC_048601.1
#            NC_048602.1
#            NC_048603.1
#            NC_048604.1
#        )
#        VEP_SPECIES="cricetulus_griseus_picr_merged"
#        ## 基因组版本名称 (必须与VEP缓存中的文件夹名称匹配)
#        export GENOME_ASSEMBLY="CriGri-PICRH-1.0"
#        ;;

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
        ## 基因组版本名称 (必须与VEP缓存中的文件夹名称匹配)
        VEP_SPECIES="cricetulus_griseus_picr"
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
        ## 基因组版本名称 (必须与VEP缓存中的文件夹名称匹配)
        VEP_SPECIES="homo_sapiens"
        export GENOME_ASSEMBLY="GRCh38"
        ;;
        
    "Cattle_ARS-UCD2.0_Ensemble")
        ## 家牛 (Bos taurus) ARS-UCD2.0 Ensembl 基因组配置
        REF_GENOME="${REF_DIR}/Cattle_E_ARSUCD2/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa"
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/Cattle_E_ARSUCD2/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa.fai"
        ## 染色体名称: 家牛有 29 条常染色体 (1-29)，以及 X, Y, MT
        main_chrs=({1..29} W Z MT)
        ## VEP 物种名称
        VEP_SPECIES="bos_taurus"
        ## 基因组版本名称 (对应 Ensembl ARS-UCD2.0)
        export GENOME_ASSEMBLY="ARS-UCD2.0"
        ;;
        
    "Chicken_E_GRCg7b")
        ## 鸡 (Ensembl GRCg7b) 配置
        REF_GENOME="${REF_DIR}/Chicken_E_GRCg7b/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
        
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/Chicken_E_GRCg7b/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.fai"
        
        ## 染色体名称, 使用.fai文件中的名称 (1-39, W, Z, MT)
        main_chrs=({1..39} W Z MT)
        
        ## 基因组版本名称 (用于VEP等工具)
        VEP_SPECIES="gallus_gallus"
        export GENOME_ASSEMBLY="bGalGal1.mat.broiler.GRCg7b"
        ;;
        
    "Dog_E_UUGSD")
        ## 家犬 (German Shepherd Dog - UU_Cfam_GSD_1.0) 配置
    
        REF_GENOME="${REF_DIR}/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
        
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/Dog_E_UUGSD/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa.fai"
        
        ## 染色体名称, 包含 1-38 号常染色体和 X 性染色体
        main_chrs=({1..38} X)
        
        ## 基因组版本名称 (对应 Ensembl 命名规范)
        VEP_SPECIES="canis_lupus_familiarisgsd"
        export GENOME_ASSEMBLY="UU_Cfam_GSD_1.0"
        ;;
        
    "GreenMonkey_E_ChlSab1.1")
        ## Green Monkey (Ensembl ChlSab1.1) 基因组配置
        REF_GENOME="${REF_DIR}/GreenMonkey_E_ChlSab1.1/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa"

        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/GreenMonkey_E_ChlSab1.1/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa.fai"

        ## 主染色体名称
        main_chrs=({1..29} X Y MT)

        ## VEP 物种名称
        VEP_SPECIES="chlorocebus_sabaeus"

        ## 基因组版本（必须匹配 VEP cache assembly 名称）
        export GENOME_ASSEMBLY="ChlSab1.1"
        ;;
        
    "Mouse_E_GRCm39")
        ## 小鼠 (Mus musculus - GRCm39 Ensembl) 配置
        REF_GENOME="${REF_DIR}/Mouse_E_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa"
        
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/Mouse_E_GRCm39/Mus_musculus.GRCm39.dna.toplevel.fa.fai"
        
        ## 染色体名称：包含 1-19 号常染色体，以及 X, Y, MT
        main_chrs=({1..19} X Y MT)
    
        ## 基因组版本名称
        VEP_SPECIES="mus_musculus"
        export GENOME_ASSEMBLY="GRCm39"
        ;;
        
    "Pig_E_Sscrofa11.1")
        ## 猪 (Sus scrofa - Sscrofa11.1 Ensembl) 配置
        REF_GENOME="${REF_DIR}/Pig_E_Sscrofa11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
        
        ## 参考基因组 FASTA 索引文件 (.fai) 的路径
        REF_FAI="${REF_DIR}/Pig_E_Sscrofa11.1/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.fai"
        
        ## 染色体名称：包含 1-18 号常染色体，以及 X, Y, MT
        main_chrs=({1..18} X Y MT)
        
        ## 基因组版本名称
        VEP_SPECIES="sus_scrofa"
        export GENOME_ASSEMBLY="Sscrofa11.1"
        ;;
        
        
    *)
        # 兜底处理：如果输入的 REF_NAME 不在上述列表中，报错退出
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

