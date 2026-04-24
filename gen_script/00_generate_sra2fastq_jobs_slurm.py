import os
import glob
import sys

# --- 用户可配置的参数 ---
PROJECT_ROOT = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects"

# 仅接收 PROJECT_NAME 参数
if len(sys.argv) < 2:
    print("错误: 请在运行脚本时提供 PROJECT_NAME 参数。")
    print("用法: python 00_generate_sra2fastq_jobs_slurm.py <PROJECT_NAME>")
    sys.exit(1) 

PROJECT_NAME = sys.argv[1]
# 内部固定 REF_NAME，触发 config.sh 中不需要参考基因组的逻辑
REF_NAME = "dont_need_ref"

PROJECT_DIR = os.path.join(PROJECT_ROOT, PROJECT_NAME)
SRA_INPUT_DIR = os.path.join(PROJECT_DIR, "00_data/dna_sra")
FASTQ_OUTPUT_DIR = os.path.join(PROJECT_DIR, "00_data/raw_fastq")
TMP_DIR = os.path.join(PROJECT_DIR, "00_data/tmp")
ALL_LOGS_DIR = os.path.join(PROJECT_DIR, "03_logs/s2q_log_file")
GENERATED_SCRIPTS_DIR = os.path.join(PROJECT_DIR, "02_jobs/sra2fastq")

CONFIG_SH_PATH = "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/gen_script/config.sh"
CONDA_PROFILE_PATH = "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV_NAME = "genome_env"

# SLURM 资源
SLURM_PARTITION = "corexd192"
SLURM_TIME = "240:00:00"
SLURM_NODES = "1"
SLURM_NTASKS_PER_NODE = "1"
SLURM_CPUS_PER_TASK = "8"
SLURM_MEM = "64gb"

os.makedirs(GENERATED_SCRIPTS_DIR, exist_ok=True)
os.makedirs(FASTQ_OUTPUT_DIR, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(ALL_LOGS_DIR, exist_ok=True)

print(f"正在查找 {SRA_INPUT_DIR} 中的 .sra 文件...")
sra_files = glob.glob(os.path.join(SRA_INPUT_DIR, "**", "*.sra"), recursive=True)

if not sra_files:
    print(f"未找到 SRA 文件。")
    exit()

for sra_file_path in sra_files:
    sample_name = os.path.basename(sra_file_path).replace(".sra", "")
    job_name = f"S2F_{sample_name}"
    output_log = os.path.join(ALL_LOGS_DIR, f"{sample_name}.log")
    script_filename = f"s2f_{sample_name}.sh"
    generated_script_full_path = os.path.join(GENERATED_SCRIPTS_DIR, script_filename)

    script_content = f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -p {SLURM_PARTITION}
#SBATCH -t {SLURM_TIME}
#SBATCH -N {SLURM_NODES}
#SBATCH --ntasks-per-node={SLURM_NTASKS_PER_NODE}
#SBATCH --cpus-per-task={SLURM_CPUS_PER_TASK}
#SBATCH --mem={SLURM_MEM}
#SBATCH -o {output_log}

echo "Start time:" && date

PROJECT_NAME="{PROJECT_NAME}"
REF_NAME="{REF_NAME}"
export PROJECT_NAME
export REF_NAME

source {CONFIG_SH_PATH}
source {CONDA_PROFILE_PATH}
conda activate {CONDA_ENV_NAME}

FASTQ_DIR="{FASTQ_OUTPUT_DIR}"
mkdir -p ${{FASTQ_DIR}}
cd ${{FASTQ_DIR}}

TMP_DIR="{TMP_DIR}"
mkdir -p ${{TMP_DIR}}

current_sra_file="{sra_file_path}"
sample="{sample_name}"

echo "Processing ${{sample}} ..."

# 运行 fasterq-dump
fasterq-dump \\
    --split-files \\
    --threads {SLURM_CPUS_PER_TASK} \\
    --outdir ${{FASTQ_DIR}} \\
    -t ${{TMP_DIR}} \\
    "${{current_sra_file}}"

if [ $? -ne 0 ]; then
    echo "Error: fasterq-dump failed for ${{sample}}."
    exit 1
fi

# 处理输出文件：
# 情况 1: 双端测序 (会有 _1.fastq 和 _2.fastq)
if [[ -f "${{sample}}_1.fastq" && -f "${{sample}}_2.fastq" ]]; then
    echo "Detected Paired-end data for ${{sample}}."
    gzip -f "${{sample}}_1.fastq" && mv "${{sample}}_1.fastq.gz" "${{sample}}_r1.fq.gz"
    gzip -f "${{sample}}_2.fastq" && mv "${{sample}}_2.fastq.gz" "${{sample}}_r2.fq.gz"

# 情况 2: 单端测序 (只有 .fastq，没有序号)
elif [[ -f "${{sample}}.fastq" ]]; then
    echo "Detected Single-end data for ${{sample}}."
    gzip -f "${{sample}}.fastq" && mv "${{sample}}.fastq.gz" "${{sample}}.fq.gz"

else
    echo "Warning: Expected output files for ${{sample}} not found."
fi

echo "SRA file ${{sample}} conversion complete."
echo "End time:" && date
"""
    with open(generated_script_full_path, "w") as f:
        f.write(script_content)

    os.chmod(generated_script_full_path, 0o755)

print(f"\n脚本生成完毕，共 {len(sra_files)} 个。")