import os
import glob
import sys

# --- 用户可配置的参数 ---
# 项目根目录，所有相对路径都将基于此目录
PROJECT_ROOT = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects"
# PROJECT_NAME = "PRJNA378939_CHO"
if len(sys.argv) < 3:
    print("错误: 请在运行脚本时提供 PROJECT_NAME 参数。")
    print("用法: python 00_generate_sra2fastq_jobs_slurm.py <PROJECT_NAME> <REF_NAME>")
    print("用法: python 00_generate_sra2fastq_jobs_slurm.py PRJNA378939_CHO CriGri-PICRH-1.0")
    sys.exit(1) # 退出脚本，表示错误
PROJECT_NAME = sys.argv[1]
REF_NAME = sys.argv[2]

# 当前项目目录
PROJECT_DIR=os.path.join(PROJECT_ROOT, PROJECT_NAME)

# SRA文件所在的输入目录
SRA_INPUT_DIR = os.path.join(PROJECT_DIR, "00_data/dna_sra")

# 转换后的FASTQ文件将要存放的输出目录
FASTQ_OUTPUT_DIR = os.path.join(PROJECT_DIR, "00_data/raw_fastq")

# 临时文件存放目录
TMP_DIR=os.path.join(PROJECT_DIR, "00_data/tmp")

# SLURM作业日志文件 (.out 和 .err) 的存放目录
ALL_LOGS_DIR = os.path.join(PROJECT_DIR, "03_logs/s2q_log_file")

# 生成的SRA转FASTQ脚本的存放目录
# 建议创建一个新的目录来存放这些自动生成的脚本
GENERATED_SCRIPTS_DIR = os.path.join(PROJECT_DIR, "02_jobs/sra2fastq")

# 配置文件路径，其中应定义 THREADS 等变量
CONFIG_SH_PATH = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/gen_script/config.sh"

# Conda环境的profile路径和环境名称
CONDA_PROFILE_PATH = "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV_NAME = "genome_env"

# SLURM作业资源请求 (可以根据实际需求调整)
SLURM_PARTITION = "corexd192"      # 指定队列名称
SLURM_TIME = "240:00:00"        # 作业最长运行时间
SLURM_NODES = "1"               # 请求节点数 
SLURM_NTASKS_PER_NODE = "1"     # 每个节点启动的任务数 (通常为1，除非作业本身是多任务并行)
SLURM_CPUS_PER_TASK = "8"       # 每个任务的核心数 (与config.sh中的THREADS保持一致)
SLURM_MEM = "40gb"              # 请求内存

# --- 确保必要的目录存在 ---
os.makedirs(GENERATED_SCRIPTS_DIR, exist_ok=True)
os.makedirs(FASTQ_OUTPUT_DIR, exist_ok=True) # fasterq-dump 也会创建，但提前创建更稳妥
os.makedirs(TMP_DIR, exist_ok=True)

print(f"正在查找 {SRA_INPUT_DIR} 中的 .sra 文件...")

# --- 查找所有 .sra 文件 ---
# 使用 glob 模块查找所有匹配的文件
sra_files = glob.glob(os.path.join(SRA_INPUT_DIR, "**", "*.sra"), recursive=True)

if not sra_files:
    print(f"在目录 '{SRA_INPUT_DIR}' 中未找到任何 .sra 文件。请检查路径或文件是否存在。")
    exit()

print(f"共找到 {len(sra_files)} 个 .sra 文件。")

# --- 生成并保存每个SRA文件的转换脚本 ---
generated_script_paths = []
for sra_file_path in sra_files:
    # 从SRA文件路径中提取样本名称 (例如 SRR12345)
    sample_name = os.path.basename(sra_file_path).replace(".sra", "")

    # 为每个作业生成唯一的名称和日志文件
    job_name = f"S2F_{sample_name}"
    output_log = os.path.join(ALL_LOGS_DIR, f"{sample_name}.log")

    # 生成的脚本文件名
    script_filename = f"s2f_{sample_name}.sh"
    generated_script_full_path = os.path.join(GENERATED_SCRIPTS_DIR, script_filename)

    # 构建shell脚本内容
    # 注意：这里使用了 f-string 来方便地插入Python变量
    script_content = f"""#!/bin/bash
#SBATCH -J {job_name}                       # 作业名为 {job_name}
#SBATCH -p {SLURM_PARTITION}                # 作业提交的分区为 {SLURM_PARTITION}
#SBATCH -t {SLURM_TIME}                     # 任务运行的最长时间
#SBATCH -N {SLURM_NODES}                    # 作业申请 {SLURM_NODES} 个节点
#SBATCH --ntasks-per-node={SLURM_NTASKS_PER_NODE} # 单节点启动的任务数为 {SLURM_NTASKS_PER_NODE}
#SBATCH --cpus-per-task={SLURM_CPUS_PER_TASK} # 单任务使用的 CPU 核心数为 {SLURM_CPUS_PER_TASK}
#SBATCH --mem={SLURM_MEM}                   # 请求内存
#SBATCH -o {output_log}                     

echo "Start time:" && date

PROJECT_NAME="{PROJECT_NAME}"
REF_NAME="{REF_NAME}"
export PROJECT_NAME # 导出 PROJECT_NAME，以便 config.sh 和生成的子脚本能访问
export REF_NAME

# 加载配置文件
source {CONFIG_SH_PATH}

# 激活环境
source {CONDA_PROFILE_PATH}
conda init
conda activate {CONDA_ENV_NAME}

# 定义输出目录 (使用Python脚本提供的绝对路径，确保一致性)
FASTQ_DIR="{FASTQ_OUTPUT_DIR}"
mkdir -p ${{FASTQ_DIR}}
cd ${{FASTQ_DIR}}

TMP_DIR="{TMP_DIR}"
mkdir -p ${{TMP_DIR}}

# 当前处理的SRA文件和样本名
current_sra_file="{sra_file_path}"
sample="{sample_name}"

echo "Processing ${{sample}} from ${{current_sra_file}} ..."

# 使用 fasterq-dump 转换 SRA -> FASTQ
# fasterq-dump 默认输出 fastq
fasterq-dump \\
    --split-files \\
    --threads ${{THREADS}} \\
    --outdir ${{FASTQ_DIR}} \\
    -t ${{TMP_DIR}} \\
    "${{current_sra_file}}"

# 检查 fasterq-dump 是否成功
if [ $? -ne 0 ]; then
    echo "Error: fasterq-dump failed for ${{sample}}."
    exit 1
fi

# 压缩 fastq
echo "Compressing FASTQ files for ${{sample}}..."
gzip -f ${{FASTQ_DIR}}/${{sample}}_1.fastq
gzip -f ${{FASTQ_DIR}}/${{sample}}_2.fastq

# 重命名为 .fq.gz 以保持流程一致
echo "Renaming FASTQ files for ${{sample}}..."
mv ${{FASTQ_DIR}}/${{sample}}_1.fastq.gz ${{FASTQ_DIR}}/${{sample}}_r1.fq.gz
mv ${{FASTQ_DIR}}/${{sample}}_2.fastq.gz ${{FASTQ_DIR}}/${{sample}}_r2.fq.gz

echo "SRA file ${{sample}} has been successfully converted to FASTQ."
echo "End time:" && date
"""
    # 将脚本内容写入文件
    with open(generated_script_full_path, "w") as f:
        f.write(script_content)

    # 赋予脚本执行权限
    os.chmod(generated_script_full_path, 0o755) # 0o755 对应 rwxr-xr-x

    generated_script_paths.append(generated_script_full_path)
    print(f"已为 '{sample_name}' 生成脚本: {generated_script_full_path}")

print("\n--- 所有脚本已生成 ---")
print(f"生成的脚本位于: {GENERATED_SCRIPTS_DIR}")
print(f"FASTQ输出将存放在: {FASTQ_OUTPUT_DIR}")
