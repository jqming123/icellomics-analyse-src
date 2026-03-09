#!/bin/bash

# 脚本功能: 为每个染色体和未定位的scaffolds生成一个 Slurm 任务脚本，用于将所有样本的 gVCF 合并到 GenomicsDB 数据库。
# 运行方式: 在 resources/src 目录下执行 `bash 02_generate_genomicsdb_slurm.sh <PROJECT_NAME> <REF_NAME>`
# 例如：bash 02_generate_genomicsdb_slurm.sh PRJEB9185_CHO CriGri-PICRH-1.0

# --- 获取参数 ---
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "用法: $0 <PROJECT_NAME> <REF_NAME>" >&2
    echo "例子: $0 MyProject hg38" >&2
    exit 1
fi

PROJECT_NAME="$1"
REF_NAME="$2"        # <-- 修改：从第二个参数获取
export PROJECT_NAME
export REF_NAME

echo "当前项目名称 (PROJECT_NAME): ${PROJECT_NAME}"
echo "当前参考基因组 (REF_NAME): ${REF_NAME}"

# 引入项目配置文件
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [ -f "${SCRIPT_DIR}/config.sh" ]; then
    source "${SCRIPT_DIR}/config.sh"
else
    echo "错误: 找不到配置文件 ${SCRIPT_DIR}/config.sh" >&2
    exit 1
fi

echo "当前使用的队列 (QUEUE_NAME): ${QUEUE_NAME}"

# --- 路径定义 ---
# GenomicsDB 的主输出目录
DB_DIR="${RESULTS_DIR}/04_genomicsdb"
# gVCF 样本映射文件
# 格式：制表符分隔，共两列，第一列是biosample_id，第二列是该样本对应的gvcf文件
GVCF_LIST="${RESULTS_DIR}/03_gvcf/sample_map.txt"

# 参考基因组 FASTA 索引文件 (.fai) 的路径改为在config.sh中设置

# 脚本生成目录
JOB_SCRIPT_DIR="${PROJECT_DIR}/02_jobs/genomicsdb_jobs"

# 创建该阶段所需的输出和日志目录
mkdir -p "${DB_DIR}" "${DB_DIR}/tmp" "${LOG_DIR}/02_joint_calling" "${JOB_SCRIPT_DIR}"

# --- 【注意】要处理的染色体列表改为在config.sh中设置 ---
# 定义主要的染色体 (NCBI RefSeq accession)
# main_chrs=(
#     NC_048595.1
#     NC_048596.1
#     NC_048597.1
#     NC_048598.1
#     NC_048599.1
#     NC_048600.1
#     NC_048601.1
#     NC_048602.1
#     NC_048603.1
#     NC_048604.1
# )

# 【新增步骤】从参考基因组的 .fai 文件中提取所有未定位的 scaffolds (以 NW_ 开头的)
UNPLACED_SCAFFOLDS_LIST="${JOB_SCRIPT_DIR}/unplaced_scaffolds.list"

# 检查文件是否存在且格式正确（非空且包含NW_开头的行）
if [[ ! -f "${UNPLACED_SCAFFOLDS_LIST}" ]] || ! grep -q "^NW_" "${UNPLACED_SCAFFOLDS_LIST}" 2>/dev/null; then
    echo "Extracting unplaced scaffolds from ${REF_FAI}..."
    grep "^NW_" "${REF_FAI}" | cut -f1 > "${UNPLACED_SCAFFOLDS_LIST}"
    echo "Found $(wc -l < ${UNPLACED_SCAFFOLDS_LIST}) unplaced scaffolds. List saved to ${UNPLACED_SCAFFOLDS_LIST}"
else
    echo "Unplaced scaffolds list already exists: ${UNPLACED_SCAFFOLDS_LIST}"
    echo "Skipping extraction. Found $(wc -l < ${UNPLACED_SCAFFOLDS_LIST}) unplaced scaffolds."
fi

# 将主染色体列表和 "unplaced" 标记合并为一个处理列表
# 为每个主染色体创建一个job，并为所有 unplaced scaffolds 创建一个合并的 job
regions_to_process=("${main_chrs[@]}" "unplaced_scaffolds")


echo "Generating GenomicsDB import jobs for each main chromosome and all unplaced scaffolds (Slurm format)..."

for region in "${regions_to_process[@]}"; do
    
    # 根据 region 是染色体还是 unplaced scaffolds 集合来设定 GATK 的 -L 参数
    if [[ "${region}" == "unplaced_scaffolds" ]]; then
        # 如果是 unplaced scaffolds，则使用文件列表作为输入
        # 后面决定忽略unplaced scaffolds，所以直接跳过了生成unplaced scaffolds对应脚本的步骤
        # gatk_l_option="-L ${UNPLACED_SCAFFOLDS_LIST}"
        continue
    else
        # 否则，直接使用染色体名称
        gatk_l_option="-L ${region}"
    fi

    JOB_SCRIPT_PATH="${JOB_SCRIPT_DIR}/${region}_gendb.sh"
    cat <<EOF > "${JOB_SCRIPT_PATH}"
#!/bin/bash
#SBATCH -p ${QUEUE_NAME}                  # 指定作业提交的分区 (队列)
#SBATCH -J GDB_${region}_${PROJECT_NAME}              # 指定作业名称
#SBATCH -o ${LOG_DIR}/02_joint_calling/${region}_gendb_%j.log  
#SBATCH --nodes=1                         # 作业申请 1 个节点
#SBATCH --ntasks-per-node=1               # 单节点启动 1 个任务
#SBATCH --cpus-per-task=2                 # 单任务使用 2 个 CPU 核心 (对应 reader-threads)
#SBATCH --mem=${MEM_LARGE}                    # 申请的内存大小
#SBATCH --time=240:00:00                # 任务运行最长时间

# --- 作业执行内容 ---
echo "Job started on: \$(hostname)"
echo "Start time: " && date

# 激活 mamba 环境 (miniforge)
source "${CONDA_PROFILE_PATH}"
conda activate "${GENOME_ENV_NAME}"

gatk --java-options "-Xmx${MEM_LARGE_M4} -Xms${MEM_MEDIUM}" GenomicsDBImport \\
  --sample-name-map "${GVCF_LIST}" \\
  --genomicsdb-workspace-path "${DB_DIR}/${region}" \\
  --tmp-dir "${DB_DIR}/tmp" \\
  --reader-threads 2 \\
  --batch-size 50 \\
  ${gatk_l_option}

echo "End time: " && date
echo "Job finished."
EOF

done
echo "Done. Slurm scripts for GenomicsDB import are generated in ${JOB_SCRIPT_DIR}"