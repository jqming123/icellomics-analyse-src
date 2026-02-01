#!/bin/bash
#
# 这个脚本用于在所有 VEP 作业完成后，汇总所有导致问题的 VCF 行。
#
# 问题被定义为：
#   VEP 作业的 .log 文件中不包含成功结束的标志 "Job finished"。
#
# 脚本会尝试从错误日志中解析出具体的出错VCF行号，并只提取该行。
# 如果无法解析出行号，则会回退到提取整个VCF块中的所有数据行。
#
# 运行方式:
#   在所有 Slurm 作业完成后，执行此脚本。
#

set -eou pipefail

# --- 配置区 ---
# 根日志目录，脚本将从此目录开始递归查找 .log 文件
PRJ_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB39258_CHO"
LOG_ROOT_DIR="${PRJ_DIR}/03_logs/04_vep_annotation_overlap_debug"
# 存放原始 VCF 输入块的根目录
VCF_CHUNKS_ROOT_DIR="${PRJ_DIR}/01_results/08_vep_annotated/vep_debug_final_round"

# 最终汇总所有问题 VCF 行的输出文件 (包含详细信息和上下文)
PROBLEM_LINES_OUTPUT_FILE="${VCF_CHUNKS_ROOT_DIR}/all_problematic_vcf_info.txt"
# 最终汇总所有问题 VCF 行的输出文件 (只包含 VCF 头和数据行，方便后续 VCF 工具处理)
PROBLEM_VCF_OUTPUT_FILE="${VCF_CHUNKS_ROOT_DIR}/all_problematic_vcf_lines.vcf"

# --- 配置区结束 ---

echo "正在清空并准备汇总文件: ${PROBLEM_LINES_OUTPUT_FILE}"
> "${PROBLEM_LINES_OUTPUT_FILE}" # 确保文件是空的

echo "正在清空并准备纯 VCF 汇总文件: ${PROBLEM_VCF_OUTPUT_FILE}"
> "${PROBLEM_VCF_OUTPUT_FILE}" # 确保文件是空的

echo "# 汇总所有导致 VEP 报错的 VCF 行" >> "${PROBLEM_LINES_OUTPUT_FILE}"
echo "# (通过查找不含 'Job finished' 的日志文件来确定)" >> "${PROBLEM_LINES_OUTPUT_FILE}"
echo "# 生成时间: $(date)" >> "${PROBLEM_LINES_OUTPUT_FILE}"

SUCCESS_MESSAGE="Job finished"
HEADER_WRITTEN=false # 标志，用于确保 VCF 头只写入一次

# 使用 grep -L (列出不包含匹配项的文件) 高效地找出所有失败任务的日志
# 使用 find ... -print0 | xargs -0 ... 来安全处理可能包含特殊字符的文件名
find "${LOG_ROOT_DIR}" -type f -name '*.log' -print0 | xargs -0 grep -L "${SUCCESS_MESSAGE}" | while read -r LOG_FILE; do
    echo "--- 发现问题日志: ${LOG_FILE} ---"

    # 从日志文件的父目录名中提取源块 ID (e.g., chunk_003_overlap_013)
    SOURCE_CHUNK_ID=$(basename "$(dirname "${LOG_FILE}")")

    # 从日志文件名中提取子块 ID (e.g., chunk_003_overlap_013_002)
    SUB_CHUNK_BASENAME=$(basename "${LOG_FILE}" .log)
    SUB_CHUNK_ID=$(echo "${SUB_CHUNK_BASENAME}" | sed 's/^vep_debug_//')

    # 构建原始输入 VCF 块的路径
    ORIGINAL_CHUNK_VCF_GZ="${VCF_CHUNKS_ROOT_DIR}/${SOURCE_CHUNK_ID}/vcf_chunks/${SUB_CHUNK_ID}.vcf.gz"

    if [ -f "${ORIGINAL_CHUNK_VCF_GZ}" ]; then
        # 如果 VCF 头尚未写入，则从第一个找到的 VCF 文件中提取并写入到纯 VCF 文件
        if ! ${HEADER_WRITTEN}; then
            echo "正在提取并写入 VCF 头信息到 ${PROBLEM_VCF_OUTPUT_FILE}..."
            # 提取所有以 '#' 开头的行作为 VCF 头
            zcat "${ORIGINAL_CHUNK_VCF_GZ}" | grep '^#' > "${PROBLEM_VCF_OUTPUT_FILE}"
            HEADER_WRITTEN=true
        fi

        # 尝试从日志文件中提取 VCF 的行号
        # 正则表达式匹配类似 "... <$fh> line 1234." 的模式，并捕获数字 1234
        ERROR_LINE_NUMBER=$(grep -oP '<\$fh> line \K[0-9]+' "${LOG_FILE}" | head -n 1)

        echo "" >> "${PROBLEM_LINES_OUTPUT_FILE}"
        echo "## 问题来源日志: ${LOG_FILE}" >> "${PROBLEM_LINES_OUTPUT_FILE}"
        echo "## 对应输入 VCF 块: ${ORIGINAL_CHUNK_VCF_GZ}" >> "${PROBLEM_LINES_OUTPUT_FILE}"

        if [[ -n "${ERROR_LINE_NUMBER}" ]]; then
            # 如果成功提取到行号，则只提取那一行
            echo "成功从日志中定位到错误行号: ${ERROR_LINE_NUMBER}。正在提取该行..."
            echo "## 精确提取的 VCF 数据行 (原始文件第 ${ERROR_LINE_NUMBER} 行):" >> "${PROBLEM_LINES_OUTPUT_FILE}"
            
            # 使用 zcat 和 sed 提取特定行
            PROBLEM_VCF_LINE=$(zcat "${ORIGINAL_CHUNK_VCF_GZ}" | sed -n "${ERROR_LINE_NUMBER}p")
            
            if [[ -n "${PROBLEM_VCF_LINE}" ]]; then
                echo "${PROBLEM_VCF_LINE}" >> "${PROBLEM_LINES_OUTPUT_FILE}"
                echo "${PROBLEM_VCF_LINE}" >> "${PROBLEM_VCF_OUTPUT_FILE}" # 写入到纯 VCF 文件
            else
                echo "# 警告: 未能在 VCF 文件中找到第 ${ERROR_LINE_NUMBER} 行。文件可能不完整。" >> "${PROBLEM_LINES_OUTPUT_FILE}"
            fi

        else
            # 如果没有提取到行号，则回退到提取整个文件 (只提取数据行)
            echo "警告: 未能从日志中解析出具体行号。将转储整个 VCF 块的数据行。"
            echo "# 警告: 未能从日志中解析出具体行号，以下是该块的全部 VCF 数据行:" >> "${PROBLEM_LINES_OUTPUT_FILE}"
            
            # 提取所有非注释行
            zcat "${ORIGINAL_CHUNK_VCF_GZ}" | grep -v '^#' >> "${PROBLEM_LINES_OUTPUT_FILE}"
            zcat "${ORIGINAL_CHUNK_VCF_GZ}" | grep -v '^#' >> "${PROBLEM_VCF_OUTPUT_FILE}" # 写入到纯 VCF 文件
        fi
        echo "" >> "${PROBLEM_LINES_OUTPUT_FILE}"
    else
        echo "警告: 找不到对应的输入 VCF 块: ${ORIGINAL_CHUNK_VCF_GZ} (来自日志: ${LOG_FILE})" >&2
    fi
done

echo ""
echo "=========================================================="
echo "问题 VCF 行汇总完成！"
echo "所有可能导致 VEP 报错的 VCF 行已汇总到以下文件:"
echo "  - 详细信息和上下文: ${PROBLEM_LINES_OUTPUT_FILE}"
echo "  - 纯 VCF 数据行 (包含 VCF 头): ${PROBLEM_VCF_OUTPUT_FILE}"
echo "请检查这些文件以进行进一步调试或后续处理。"
echo "=========================================================="