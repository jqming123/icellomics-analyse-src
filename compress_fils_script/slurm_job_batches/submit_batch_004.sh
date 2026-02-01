#!/bin/bash
#SBATCH -p corexd192
#SBATCH -J slurm_compress_bams_batch_004 # 每个作业有唯一的名称
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/compress_fils_script/logs/slurm_compress_bams_batch_004_%j.log # 日志文件包含批次名
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=240:00:00

# 获取当前时间戳 (这里为了日志输出，再次定义，也可以依赖模板中的定义)
get_timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

echo "[2025-12-16 15:31:51] --- 开始运行SLURM批处理任务: batch_004 ---"
echo "[2025-12-16 15:31:51] 将切换到主工作目录: '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result'"
cd "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result" || { echo "[2025-12-16 15:31:51] 错误: 无法切换到主工作目录 '/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result'. 退出."; exit 1; }

# 定义此批次要处理的样本目录
VALID_SAMPLE_DIRS=(
    "./SRR25745311" "./SRR25745312" "./SRR25745315" "./SRR25745316" "./SRR25745317" "./SRR25745318" "./SRR25745319" "./SRR25745320" "./SRR25745321" "./SRR25745322" "./SRR25745323" "./SRR25745324" "./SRR25745326" "./SRR25745327" "./SRR25745328" "./SRR25745329" "./SRR25745330" "./SRR25745331" "./SRR25745332" "./SRR25745333" "./SRR25745334" "./SRR25745335" "./SRR25745336" "./SRR25745337" "./SRR25745338" "./SRR25745339" "./SRR25745340" "./SRR25745342" "./SRR25745343" "./SRR25745344" "./SRR25745345" "./SRR25745346" "./SRR25745347" "./SRR25745348" "./SRR25745349" "./SRR25745350" "./SRR25745351" "./SRR25745353" "./SRR25745354" "./SRR25745355" "./SRR25745363" "./SRR25745364" "./SRR25745365" "./SRR25745366" "./SRR25745367" "./SRR25745368" "./SRR25745370" "./SRR25745371" "./SRR25745372" "./SRR25745373" 
)

# --- 嵌入核心压缩逻辑 ---
# 文件名: compression_logic_template.sh
# 描述: 这是一个核心压缩逻辑模板，其内容将被复制到每个生成的Slurm脚本中。

# 遇到错误立即退出(-e), 未定义变量时报错(-u), 管道中任一命令失败则返回失败(-o pipefail)
set -eo pipefail

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

# --- 用户需要修改的部分（这些路径会直接嵌入到生成的脚本中） ---
# 基因组参考序列（用于基因组比对文件）
GENOME_FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
# 转录本参考序列（用于转录组比对文件）
TRANSCRIPTOME_FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/Homo_sapiens.GRCh38.cdna.all.fa"
# -------------------------

# 动态获取线程数, 若不在Slurm环境中则默认为8
THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "[$(get_timestamp)] 脚本开始运行. 将使用 ${THREADS} 个线程."

# --- 【注意】 VALID_SAMPLE_DIRS 数组和 MAIN_WORK_DIRECTORY 的切换将在生成的脚本中显式定义 ---
# 此处假定 VALID_SAMPLE_DIRS 已经包含要处理的样本目录列表 (例如: ./ERR10047670)
# 并且脚本已经切换到正确的主工作目录

TOTAL_SAMPLES_TO_PROCESS=${#VALID_SAMPLE_DIRS[@]}

if [ "$TOTAL_SAMPLES_TO_PROCESS" -eq 0 ]; then
    echo "[$(get_timestamp)] 警告: 未检测到任何样本目录需要处理. 脚本将退出."
    exit 0
fi

# 检查参考文件
check_reference "$GENOME_FASTA" "基因组" || exit 1
check_reference "$TRANSCRIPTOME_FASTA" "转录本" || exit 1

echo "[$(get_timestamp)] 使用基因组参考: $GENOME_FASTA"
echo "[$(get_timestamp)] 使用转录本参考: $TRANSCRIPTOME_FASTA"

# 激活 conda 环境
source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate "RNAseq_E4"

echo "[$(get_timestamp)] 开始处理当前批次的BAM文件..."
echo "======================================================="

echo "[$(get_timestamp)] 将处理 ${TOTAL_SAMPLES_TO_PROCESS} 个样本目录."
CURRENT_SAMPLE_COUNT=0

# 遍历预先识别的有效样本目录
for sample_dir in "${VALID_SAMPLE_DIRS[@]}"; do
    # sample_dir 是相对于 MAIN_WORK_DIRECTORY 的路径 (例如: ./ERR10047670)
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

        # <<< 根据BAM文件名选择不同的处理策略 >>>

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
# --- 核心压缩逻辑嵌入结束 ---

echo "[2025-12-16 15:31:51] --- SLURM批处理任务: batch_004 运行结束 ---"
