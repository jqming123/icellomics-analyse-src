#!/bin/bash
# ============================================================
# 批量检查链特异性脚本
# 使用前请修改下方三个变量：
#   BIOPROJECT   : BioProject 编号，例如 PRJNA681243
#   CELL_LINE    : 细胞系名称，例如 MDCK
#   ANIMAL_REF   : ref_genome 下的基因组文件夹名，例如 Dog_E_UUGSD
# ============================================================

source /hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/bin/activate inferexp

BIOPROJECT="PRJNA979805"
CELL_LINE="HT1080"
ANIMAL_REF="hg38_Ensembl"

# ============================================================
# 以下路径根据上方变量自动生成，无需修改
# ============================================================

RUN_LIST="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_rawdata/${BIOPROJECT}_${CELL_LINE}/sra_runid.txt"

BED_FILE=$(ls /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/${ANIMAL_REF}/*.bed12 2>/dev/null | head -1)

BAM_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"

LOG_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/logs/${BIOPROJECT}_${CELL_LINE}"
LOG_FILE="${LOG_DIR}/strand_results.log"

# ============================================================
# 运行前检查
# ============================================================

if [[ ! -f "$RUN_LIST" ]]; then
    echo "ERROR: sra_runid.txt 不存在: $RUN_LIST"
    exit 1
fi

if [[ ! -f "$BED_FILE" ]]; then
    echo "ERROR: 未找到 .bed12 文件，请检查 ANIMAL_REF 是否正确: $ANIMAL_REF"
    exit 1
fi

mkdir -p "$LOG_DIR"

# Write run info to top of log file
{
    echo "BIOPROJECT=\"${BIOPROJECT}\""
    echo "CELL_LINE=\"${CELL_LINE}\""
    echo "ANIMAL_REF=\"${ANIMAL_REF}\""
    echo "BED_FILE=${BED_FILE}"
    echo "LOG_FILE=${LOG_FILE}"
    echo ""
} | tee "$LOG_FILE"

# ============================================================
# 批量运行 infer_experiment.py
# ============================================================

while read -r RUN_ID || [[ -n "$RUN_ID" ]]; do
    [[ -z "$RUN_ID" ]] && continue

    BAM="${BAM_DIR}/${RUN_ID}/star/${RUN_ID}.gsd.Aligned.sortedByCoord.out.bam"

    echo "=== ${RUN_ID} ===" | tee -a "$LOG_FILE"

    if [[ ! -f "$BAM" ]]; then
        echo "WARNING: BAM 文件不存在，跳过: $BAM" | tee -a "$LOG_FILE"
        continue
    fi

    for attempt in 1 2 3; do
        result=$(infer_experiment.py -r "$BED_FILE" -i "$BAM" 2>&1)
        if echo "$result" | grep -q "Fraction of reads"; then
            echo "$result" | tee -a "$LOG_FILE"
            break
        else
            echo "Attempt ${attempt} failed, retrying in 3s..." | tee -a "$LOG_FILE"
            sleep 3
        fi
    done
    echo "" | tee -a "$LOG_FILE"

done < "$RUN_LIST"

echo "完成，结果已保存到: $LOG_FILE"