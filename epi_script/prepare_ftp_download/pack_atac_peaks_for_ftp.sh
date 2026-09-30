#!/bin/bash
#SBATCH -p corexd192
#SBATCH -J pack_atac_peaks
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects/prepare_ftp_download/slurm_logs/pack_atac_peaks_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=12:00:00

set -euo pipefail

echo "============================================================"
echo "SLURM Job: pack_atac_peaks"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "============================================================"

BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine"
SRC_DIR="${BASE_DIR}/epigen_projects/prepare_ftp_download"
PACK_DIR="${SRC_DIR}/packs"
LOG_DIR="${SRC_DIR}/slurm_logs"

mkdir -p "$PACK_DIR"
mkdir -p "$LOG_DIR"

echo ""
echo "=== Packing ATAC peak files by cell line ==="
echo "Source: ${SRC_DIR}"
echo "Pack dir: ${PACK_DIR}"
echo ""

total_celllines=0

for cl_dir in "$SRC_DIR"/*/; do
    [ -d "$cl_dir" ] || continue
    cl_name=$(basename "$cl_dir")

    # 跳过非细胞系目录（packs, slurm_logs 等）
    case "$cl_name" in
        packs|slurm_logs) continue ;;
    esac

    zip_file="${PACK_DIR}/${cl_name}_atac_peaks.zip"

    echo "Packing: ${cl_name} -> ${zip_file}"

    # 进入源目录，用相对路径打包，解压后直接就是细胞系目录
    (
        cd "$SRC_DIR"
        zip -r "$zip_file" "$cl_name"/
    )

    zip_size=$(du -h "$zip_file" | cut -f1)
    total_celllines=$((total_celllines + 1))
    echo "  Done: ${cl_name} -- size: ${zip_size}"
done

echo ""
echo "============================================================"
echo "All ${total_celllines} cell lines packed successfully"
echo "Output directory: ${PACK_DIR}"
echo "Finish: $(date)"
echo "============================================================"
