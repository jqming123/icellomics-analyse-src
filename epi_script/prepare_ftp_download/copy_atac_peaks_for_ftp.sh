#!/bin/bash
set -euo pipefail

# ================================================================
# 复制 ATAC peak calling 结果文件到 prepare_ftp_download 目录
# 目标结构: prepare_ftp_download/<cell_line>/<bioproject>/<三种文件>
# 三种文件: *_peaks.narrowPeak, *_peaks.xls, *_summits.bed
# ================================================================

BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine"
EPIGEN_DIR="${BASE_DIR}/epigen_projects"
DEST_DIR="${EPIGEN_DIR}/prepare_ftp_download"

# 目标细胞系列表
CELLLINES=(
    CHO
    H9
    HEK293
    HT1080
    Vero
    Hela
    MDCK
    PK-15
    MRC-5
    WI-38
)

# 需要复制的三种文件后缀
SUFFIXES=("_peaks.narrowPeak" "_peaks.xls" "_summits.bed")

echo "============================================================"
echo "Copy ATAC peak calling results to prepare_ftp_download"
echo "Source: ${EPIGEN_DIR}"
echo "Destination: ${DEST_DIR}"
echo "Start: $(date)"
echo "============================================================"

total_copied=0
total_skipped=0

for CL in "${CELLLINES[@]}"; do
    echo ""
    echo "--- Processing cell line: ${CL} ---"

    cl_total=0

    # 遍历所有匹配该细胞系的 PRJ 项目目录
    for proj_dir in "${EPIGEN_DIR}"/PRJ*"${CL}"*/; do
        [ -d "$proj_dir" ] || continue
        proj_name=$(basename "$proj_dir")
        peak_dir="${proj_dir}/1_result/3_peak_calling"

        if [ ! -d "$peak_dir" ]; then
            echo "  SKIP: ${proj_name} -- no 3_peak_calling directory"
            total_skipped=$((total_skipped + 1))
            continue
        fi

        # 收集所有 biosample 的 prefix（通过 narrowPeak 文件推断）
        narrowpeak_files=("$peak_dir"/*_peaks.narrowPeak)

        if [ ! -e "${narrowpeak_files[0]}" ]; then
            echo "  SKIP: ${proj_name} -- no narrowPeak files found"
            total_skipped=$((total_skipped + 1))
            continue
        fi

        dest_proj_dir="${DEST_DIR}/${CL}/${proj_name}"
        mkdir -p "$dest_proj_dir"

        biosample_copied=0
        for np in "${narrowpeak_files[@]}"; do
            base_name=$(basename "$np" _peaks.narrowPeak)
            all_exist=true

            # 检查三种文件是否都存在
            for suf in "${SUFFIXES[@]}"; do
                src_file="${peak_dir}/${base_name}${suf}"
                if [ ! -f "$src_file" ]; then
                    echo "  WARN: ${proj_name}/${base_name}${suf} -- missing, skipping this biosample"
                    all_exist=false
                    break
                fi
            done

            if [ "$all_exist" = true ]; then
                for suf in "${SUFFIXES[@]}"; do
                    src_file="${peak_dir}/${base_name}${suf}"
                    cp "$src_file" "$dest_proj_dir/"
                done
                biosample_copied=$((biosample_copied + 1))
            fi
        done

        if [ "$biosample_copied" -gt 0 ]; then
            echo "  OK: ${proj_name} -- ${biosample_copied} biosamples copied"
            cl_total=$((cl_total + biosample_copied))
        else
            echo "  SKIP: ${proj_name} -- no complete biosamples"
            total_skipped=$((total_skipped + 1))
        fi
    done

    echo "  >>> ${CL}: ${cl_total} biosamples copied"
    total_copied=$((total_copied + cl_total))
done

echo ""
echo "============================================================"
echo "Summary: ${total_copied} biosamples copied, ${total_skipped} bioprojects skipped"
echo "Finish: $(date)"
echo "============================================================"
