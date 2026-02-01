#!/bin/bash
# 文件名: generate_slurm_batches.sh
# 描述: 用于生成批量压缩BAM文件到CRAM文件的Slurm提交脚本。
# 这个脚本本身不使用Slurm，直接在终端运行即可。

set -euo pipefail

# --- 配置部分 (请根据您的实际路径修改) ---
SCRIPTS_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/compress_fils_script"
# 指向核心压缩逻辑模板的完整路径
# 这个文件不直接运行，它的内容会被复制到每个生成的SLURM脚本中。
COMPRESSION_LOGIC_TEMPLATE="${SCRIPTS_DIR}/compression_logic_template.sh"
# 存放生成的SLURM提交脚本的目录
GENERATED_SLURM_SCRIPTS_DIR="${SCRIPTS_DIR}/slurm_job_batches"
# 日志文件存放的基础目录
LOG_BASE_DIR="${SCRIPTS_DIR}/logs"

# 主要工作目录，所有样本目录都在这个目录之下
MAIN_WORK_DIRECTORY="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"
# 每个SLURM任务处理的样本目录数量
BATCH_SIZE=50
# ---------------------------------------------

# --- 函数定义: 获取当前时间戳 (只用于此生成脚本的日志) ---
get_timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

echo "[$(get_timestamp)] 脚本开始运行，用于生成批量SLURM压缩任务."
echo "[$(get_timestamp)] 核心压缩逻辑模板路径: $COMPRESSION_LOGIC_TEMPLATE"
echo "[$(get_timestamp)] 主工作目录: $MAIN_WORK_DIRECTORY"
echo "[$(get_timestamp)] 每批处理样本数: $BATCH_SIZE"
echo "[$(get_timestamp)] 生成的SLURM脚本将存放于: $GENERATED_SLURM_SCRIPTS_DIR"

# 检查模板文件是否存在
if [ ! -f "$COMPRESSION_LOGIC_TEMPLATE" ]; then
    echo "[$(get_timestamp)] 错误: 未找到核心压缩逻辑模板文件 '$COMPRESSION_LOGIC_TEMPLATE'."
    echo "[$(get_timestamp)] 请确保路径正确。"
    exit 1
fi

# 读取模板文件内容
COMPRESSION_CODE=$(cat "$COMPRESSION_LOGIC_TEMPLATE")

# 创建存放生成脚本的目录
echo "[$(get_timestamp)] 正在创建或确认生成脚本目录: '$GENERATED_SLURM_SCRIPTS_DIR'"
mkdir -p "$GENERATED_SLURM_SCRIPTS_DIR"
mkdir -p "$LOG_BASE_DIR" # 确保日志目录存在

# 切换到主工作目录以查找样本 (此脚本执行时切换，不是生成的脚本)
echo "[$(get_timestamp)] 切换到主工作目录: '$MAIN_WORK_DIRECTORY' 进行目录查找."
cd "$MAIN_WORK_DIRECTORY" || { echo "[$(get_timestamp)] 错误: 无法切换到主工作目录 '$MAIN_WORK_DIRECTORY'. 退出."; exit 1; }

# 查找所有需要处理的样本目录
echo "[$(get_timestamp)] 正在查找包含 'star' 子目录的样本目录..."
# find . -maxdepth 2 -type d -name "star" -printf '%h\n' 会找到例如 './ERR10047670' 这样的相对路径
mapfile -t ALL_SAMPLE_DIRS < <(find . -maxdepth 2 -type d -name "star" -printf '%h\n' | sort -u)

NUM_TOTAL_DIRS=${#ALL_SAMPLE_DIRS[@]}

if [ "$NUM_TOTAL_DIRS" -eq 0 ]; then
    echo "[$(get_timestamp)] 警告: 未检测到任何包含 'star' 子目录的样本目录需要处理. 脚本将退出."
    exit 0
fi

echo "[$(get_timestamp)] 成功检测到 ${NUM_TOTAL_DIRS} 个样本目录需要处理."

# 计算需要生成的批次数量 (向上取整)
NUM_BATCHES=$(( (NUM_TOTAL_DIRS + BATCH_SIZE - 1) / BATCH_SIZE ))
echo "[$(get_timestamp)] 将生成 ${NUM_BATCHES} 个SLURM提交脚本."

# 循环生成每个批次的SLURM脚本
for (( i=0; i<NUM_BATCHES; i++ )); do
    START_INDEX=$(( i * BATCH_SIZE ))
    
    # 获取当前批次的目录列表 (Bash 4+ 数组切片功能)
    BATCH_DIRS_RAW=("${ALL_SAMPLE_DIRS[@]:$START_INDEX:$BATCH_SIZE}")

    # 将目录列表格式化为 bash 数组定义字符串
    # 例如: '"./ERR10047670"' '"./ERR10047671"' ...
    BATCH_DIRS_FORMATTED=""
    for dir in "${BATCH_DIRS_RAW[@]}"; do
        BATCH_DIRS_FORMATTED+="\"$dir\" "
    done

    BATCH_NAME="batch_$(printf "%03d" $((i+1)))" # 例如: batch_001, batch_002
    GENERATED_SLURM_SCRIPT="${GENERATED_SLURM_SCRIPTS_DIR}/submit_${BATCH_NAME}.sh"
    
    echo "[$(get_timestamp)] 正在生成 SLURM 脚本: '$GENERATED_SLURM_SCRIPT' (包含 ${#BATCH_DIRS_RAW[@]} 个目录)..."

    # 使用heredoc写入SLURM脚本的内容
    # 注意：这里的 EOF 必须是裸的，不能有缩进
    cat > "$GENERATED_SLURM_SCRIPT" <<EOF
#!/bin/bash
#SBATCH -p corexd192
#SBATCH -J compress_bams_${BATCH_NAME} # 每个作业有唯一的名称
#SBATCH -o ${LOG_BASE_DIR}/slurm_compress_bams_${BATCH_NAME}_%j.log # 日志文件包含批次名
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=240:00:00

# 获取当前时间戳 (这里为了日志输出，再次定义，也可以依赖模板中的定义)
get_timestamp() {
    date +"%Y-%m-%d %H:%M:%S"
}

echo "[$(get_timestamp)] --- 开始运行SLURM批处理任务: ${BATCH_NAME} ---"
echo "[$(get_timestamp)] 将切换到主工作目录: '${MAIN_WORK_DIRECTORY}'"
cd "${MAIN_WORK_DIRECTORY}" || { echo "[$(get_timestamp)] 错误: 无法切换到主工作目录 '$MAIN_WORK_DIRECTORY'. 退出."; exit 1; }

# 定义此批次要处理的样本目录
VALID_SAMPLE_DIRS=(
    ${BATCH_DIRS_FORMATTED}
)

# --- 嵌入核心压缩逻辑 ---
${COMPRESSION_CODE}
# --- 核心压缩逻辑嵌入结束 ---

echo "[$(get_timestamp)] --- SLURM批处理任务: ${BATCH_NAME} 运行结束 ---"
EOF

    # 使生成的脚本可执行
    chmod +x "$GENERATED_SLURM_SCRIPT"
    
done

echo "[$(get_timestamp)] 所有 ${NUM_BATCHES} 个SLURM提交脚本已生成于 '$GENERATED_SLURM_SCRIPTS_DIR'."
echo "[$(get_timestamp)] 您可以通过遍历该目录并使用 'sbatch <脚本名>' 来提交它们。"
echo "[$(get_timestamp)] 例如: for script in ${GENERATED_SLURM_SCRIPTS_DIR}/submit_batch_*.sh; do sbatch \$script; done"
echo "[$(get_timestamp)] 脚本运行结束."