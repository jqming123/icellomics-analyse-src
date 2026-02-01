#!/bin/bash

# 脚本功能: 为每个染色体生成一个 PBS 任务脚本，用于将所有样本的 gVCF 合并到 GenomicsDB 数据库。
# 运行方式: 在 01_scripts/ 目录下执行 `bash 02_generate_genomicsdb_jobs.sh`

# 引入项目配置文件
source ./config.sh

# --- 路径定义 ---
# GenomicsDB 的主输出目录
DB_DIR="${RESULTS_DIR}/04_genomicsdb"
# gVCF 样本映射文件
GVCF_LIST="${RESULTS_DIR}/03_gvcf/sample_map.txt"
# **读取**由 `01a` 脚本生成的分区文件所在目录
INTERVAL_DIR="${DB_DIR}/intervals_by_size" 

# 创建该阶段所需的输出和日志目录
mkdir -p "$DB_DIR" "$DB_DIR/tmp" "$LOG_DIR/02_joint_calling"

# # 定义要并行处理的染色体列表
# chrs=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11)

# echo "Generating GenomicsDB import jobs for each chromosome..."

# for chromosome in "${chrs[@]}"; do
#     JOB_SCRIPT_NAME="${chromosome}_gendb.sh"
#     cat <<EOF > "${JOB_SCRIPT_NAME}"
# #PBS -q core40
# #PBS -l walltime=1001:00:00,nodes=1:ppn=2,mem=${MEM_LARGE}
# #HSCHED -s hschedd
# #PBS -o ${LOG_DIR}/02_joint_calling/${chromosome}_gendb.out
# #PBS -e ${LOG_DIR}/02_joint_calling/${chromosome}_gendb.err
# #PBS -N GDB_${chromosome}

# LOGFILE="${LOG_DIR}/02_joint_calling/${chromosome}_gendb.log"
# touch "\${LOGFILE}"
# exec > "\${LOGFILE}" 2>&1

# echo "Start time: " && date

# # 激活 mamba 环境 (miniforge)
# source /gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/etc/profile.d/conda.sh
# conda activate ${MAMBA_ENV_NAME}

# gatk --java-options "-Xmx${MEM_LARGE} -Xms${MEM_LARGE}" GenomicsDBImport \\
#   --sample-name-map "${GVCF_LIST}" \\
#   --genomicsdb-workspace-path "${DB_DIR}/${chromosome}" \\
#   --tmp-dir "${DB_DIR}/tmp" \\
#   --reader-threads 2 \\
#   --batch-size 50 \\
#   -L "${chromosome}"

# echo "End time: " && date
# EOF

# done
# echo "Done. PBS scripts for GenomicsDB import are generated in the current directory."

# 检查分区文件是否存在
if [ ! -d "$INTERVAL_DIR" ] || [ -z "$(ls -A "$INTERVAL_DIR")" ]; then
    echo "Error: Interval directory '${INTERVAL_DIR}' is empty or does not exist."
    echo "Please run the '01a_prepare_intervals_by_size.sh' script first."
    exit 1
fi

num_intervals=$(ls -1 "${INTERVAL_DIR}/scatter_"* | wc -l)
echo "Found ${num_intervals} interval files. Generating a PBS job for each..."

# 3. 遍历所有分区文件，生成对应的 PBS 任务脚本
for interval_file in ${INTERVAL_DIR}/scatter_*; do
    
    job_id=$(basename "$interval_file" | sed 's/scatter_//')
    JOB_SCRIPT_NAME="scatter_${job_id}_gendb.sh"
    cat <<EOF > "${JOB_SCRIPT_NAME}"
#PBS -q core40
#PBS -l walltime=1001:00:00,nodes=1:ppn=${THREADS},mem=${MEM_LARGE}
#HSCHED -s hschedd
#PBS -o ${LOG_DIR}/02_joint_calling/scatter_${job_id}_gendb.out
#PBS -e ${LOG_DIR}/02_joint_calling/scatter_${job_id}_gendb.err
#PBS -N GDB_scatter_${job_id}

LOGFILE="${LOG_DIR}/02_joint_calling/scatter_${job_id}_gendb.log"
touch "\${LOGFILE}"
exec > "\${LOGFILE}" 2>&1

echo "Start time: " && date
echo "Processing scaffolds from interval file: ${interval_file}"

source /gpfs/zhaowm_group/gaoxiaojing/software/miniforge3/etc/profile.d/conda.sh
conda activate ${MAMBA_ENV_NAME}

gatk --java-options "-Xmx${MEM_LARGE} -Xms${MEM_LARGE}" GenomicsDBImport \\
  --sample-name-map "${GVCF_LIST}" \\
  --genomicsdb-workspace-path "${DB_DIR}/scatter_${job_id}" \\
  --tmp-dir "${DB_DIR}/tmp" \\
  --reader-threads 2 \\
  --batch-size 50 \\
  -L "${interval_file}"

echo "End time: " && date
EOF

done

echo "Done. ${num_intervals} PBS scripts for GenomicsDB import are generated in the current directory."





