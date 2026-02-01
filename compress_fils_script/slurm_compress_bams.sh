#!/bin/bash
#SBATCH -p corexd192                  # 指定作业提交的分区 (队列)
#SBATCH -J slurm_compress_bams
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/trs_script/slurm_compress_bams_%j.log
#SBATCH --nodes=1                         # 作业申请 1 个节点
#SBATCH --ntasks-per-node=1               # 单节点启动 1 个任务
#SBATCH --cpus-per-task=8                 # 使用配置的线程数
#SBATCH --mem=60G                         # 申请的内存大小
#SBATCH --time=240:00:00                  # 任务运行最长时间

# 遇到错误立即退出(-e), 未定义变量时报错(-u), 管道中任一命令失败则返回失败(-o pipefail)
set -euo pipefail

# --- 函数定义: 获取当前时间戳 ---
get_timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

# --- 函数定义: 检查参考文件并创建索引 ---
check_reference() {
    local ref_file=$1
    local ref_name=$2

    echo "[$(get_timestamp)] 正在检查${ref_name}参考文件: '$ref_file'..."

    if [ ! -f "$ref_file" ]; then
        echo "[$(get_timestamp)] 错误：${ref_name}参考文件未找到于 '$ref_file'"
        return 1
    fi

    if [ ! -f "${ref_file}.fai" ]; then
        echo "[$(get_timestamp)] 创建${ref_name}参考索引..."
        if ! samtools faidx "$ref_file"; then
            echo "[$(get_timestamp)] 错误：创建${ref_name}索引失败"
            return 1
        fi
        echo "[$(get_timestamp)] ${ref_name}参考索引创建成功"
    else
        echo "[$(get_timestamp)] ${ref_name}参考索引已存在"
    fi

    return 0
}

# --- 用户需要修改的部分 ---
# 基因组参考序列（用于基因组比对文件）
GENOME_FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
# 转录本参考序列（用于转录组比对文件）
TRANSCRIPTOME_FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/Homo_sapiens.GRCh38.cdna.all.fa"
WORK_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"
# -------------------------

# 动态获取线程数, 若不在Slurm环境中则默认为8
THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[$(get_timestamp)] 脚本开始运行. 将使用 ${THREADS} 个线程."

# 检查参考文件
check_reference "$GENOME_FASTA" "基因组" || exit 1
check_reference "$TRANSCRIPTOME_FASTA" "转录本" || exit 1

echo "[$(get_timestamp)] 使用基因组参考: $GENOME_FASTA"
echo "[$(get_timestamp)] 使用转录本参考: $TRANSCRIPTOME_FASTA"
echo "[$(get_timestamp)] 正在切换工作目录至: $WORK_DIR"

# 激活 conda 环境
source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate "RNAseq_E4"

# 切换到工作目录，如果失败则退出
cd "${WORK_DIR}" || { echo "[$(get_timestamp)] 错误: 无法切换到工作目录 '$WORK_DIR'. 退出."; exit 1; }

echo "[$(get_timestamp)] 开始在当前目录的子文件夹中查找并压缩BAM文件..."
echo "======================================================="

# --- 进度条功能: 计算需要处理的总样本数 ---
mapfile -t VALID_SAMPLE_DIRS < <(find . -maxdepth 2 -type d -name "star" -printf '%h\n' | sort -u)
TOTAL_SAMPLES_TO_PROCESS=${#VALID_SAMPLE_DIRS[@]}

if [ "$TOTAL_SAMPLES_TO_PROCESS" -eq 0 ]; then
    echo "[$(get_timestamp)] 警告: 未检测到任何包含 'star' 子目录的样本目录需要处理. 脚本将退出."
    exit 0
fi

echo "[$(get_timestamp)] 检测到 ${TOTAL_SAMPLES_TO_PROCESS} 个样本目录包含 'star' 目录需要处理."
CURRENT_SAMPLE_COUNT=0

# 遍历预先识别的有效样本目录
for sample_dir in "${VALID_SAMPLE_DIRS[@]}"; do
    star_dir="${sample_dir}/star"
    CURRENT_SAMPLE_COUNT=$((CURRENT_SAMPLE_COUNT + 1))
    echo "[$(get_timestamp)] --- 正在处理样本: $sample_dir (${CURRENT_SAMPLE_COUNT}/${TOTAL_SAMPLES_TO_PROCESS}) ---"

    while IFS= read -r -d $'\0' bam_file; do
        cram_file="${bam_file%.bam}.cram"

        if [ -f "$cram_file" ]; then
            echo "[$(get_timestamp)] 提示: '$cram_file' 已存在，跳过。"
            continue
        fi

        bam_filename=$(basename "$bam_file")

        # <<< 核心修改: 根据BAM文件名选择不同的处理策略 >>>

        if [[ "$bam_filename" == *"Aligned.sortedByCoord.out.bam" ]]; then
            # --- 策略1: 对于已按基因组坐标排序的BAM文件，直接压缩 ---
            REFERENCE="$GENOME_FASTA"
            file_type="基因组"

            echo "[$(get_timestamp)] 正在压缩${file_type}比对文件: '$bam_file' -> '$cram_file'"
            echo "[$(get_timestamp)] 使用参考: $(basename "$REFERENCE")"

            if samtools view -@ "$THREADS" -C -T "$REFERENCE" -o "$cram_file" "$bam_file"; then
                echo "[$(get_timestamp)] 成功创建: '$cram_file'"
                echo "[$(get_timestamp)] 验证CRAM文件..."
                if samtools quickcheck "$cram_file"; then
                    echo "[$(get_timestamp)] CRAM文件验证通过"
                    rm "$bam_file"
                    echo "[$(get_timestamp)] 已删除原始BAM文件: '$bam_file'"
                else
                    echo "[$(get_timestamp)] 警告: CRAM文件验证失败，保留原始BAM文件"
                    rm -f "$cram_file"
                    echo "[$(get_timestamp)] 已删除可能损坏的CRAM文件: '$cram_file'"
                fi
            else
                echo "[$(get_timestamp)] 错误: 压缩 '$bam_file' 失败。"
                rm -f "$cram_file" # 清理可能创建的不完整文件
            fi

        elif [[ "$bam_filename" == *"Aligned.toTranscriptome.out.bam" ]]; then
            # --- 策略2: 对于按转录本比对(默认按读段名排序)的BAM文件，先排序再压缩 ---
            REFERENCE="$TRANSCRIPTOME_FASTA"
            file_type="转录组"
            sorted_bam_file="${bam_file%.bam}.coordsorted.tmp.bam"

            echo "[$(get_timestamp)] 提示: '${bam_filename}' 需要先按坐标排序..."
            
            # 步骤1: 使用 samtools sort 进行坐标排序
            if samtools sort -@ "$THREADS" -o "$sorted_bam_file" "$bam_file"; then
                echo "[$(get_timestamp)] 排序成功，生成临时文件: '$sorted_bam_file'"

                # 步骤2: 使用排序后的临时BAM文件进行压缩
                echo "[$(get_timestamp)] 正在压缩${file_type}比对文件: '$sorted_bam_file' -> '$cram_file'"
                echo "[$(get_timestamp)] 使用参考: $(basename "$REFERENCE")"
                
                if samtools view -@ "$THREADS" -C -T "$REFERENCE" -o "$cram_file" "$sorted_bam_file"; then
                    echo "[$(get_timestamp)] 成功创建: '$cram_file'"
                    echo "[$(get_timestamp)] 验证CRAM文件..."
                    
                    # 步骤3: 验证CRAM文件
                    if samtools quickcheck "$cram_file"; then
                        echo "[$(get_timestamp)] CRAM文件验证通过"
                        # 成功后删除原始BAM和临时排序BAM
                        rm "$bam_file" "$sorted_bam_file"
                        echo "[$(get_timestamp)] 已删除原始BAM及临时排序文件"
                    else
                        echo "[$(get_timestamp)] 警告: CRAM文件验证失败，保留原始BAM文件"
                        # 失败则删除损坏的CRAM和临时排序BAM
                        rm -f "$cram_file" "$sorted_bam_file"
                        echo "[$(get_timestamp)] 已删除损坏的CRAM及临时排序文件"
                    fi
                else
                    echo "[$(get_timestamp)] 错误: 压缩 '$sorted_bam_file' 失败"
                    rm -f "$sorted_bam_file" # 清理临时文件
                fi
            else
                echo "[$(get_timestamp)] 错误: 排序 '$bam_file' 失败，跳过此文件"
            fi
            
        else
            # --- 策略3: 对于未知类型的BAM文件，发出警告并跳过 ---
            echo "[$(get_timestamp)] 警告: 未知BAM文件类型 '$bam_filename'，无法确定参考序列和排序方式，跳过处理。"
        fi

    done < <(find "$star_dir" -maxdepth 1 -type f -name "*.bam" -print0)

    # 检查和报告该样本目录的文件状态
    bam_remains=$(find "$star_dir" -maxdepth 1 -type f -name "*.bam" | wc -l)
    if [ "$bam_remains" -eq 0 ]; then
        echo "[$(get_timestamp)] 提示: '${star_dir}' 中所有BAM文件已成功处理。"
    else
        echo "[$(get_timestamp)] 警告: '${star_dir}' 中仍有 ${bam_remains} 个BAM文件未处理或处理失败。"
    fi
    cram_count=$(find "$star_dir" -maxdepth 1 -type f -name "*.cram" | wc -l)
    echo "[$(get_timestamp)] 提示: '${star_dir}' 中现有 ${cram_count} 个CRAM文件。"
    echo ""

done

echo "======================================================="
echo "[$(get_timestamp)] 所有目录处理完毕。"
echo "[$(get_timestamp)] 脚本运行结束。"

echo "======================================================="
echo "[$(get_timestamp)] 汇总统计:"
echo "[$(get_timestamp)] 处理的样本总数: ${TOTAL_SAMPLES_TO_PROCESS}"
echo "[$(get_timestamp)] 当前工作目录下的CRAM文件总数: $(find -L . -type f -name "*.cram" | wc -l)"
echo "[$(get_timestamp)] 当前工作目录下的BAM文件总数: $(find -L . -type f -name "*.bam" | wc -l)"