#!/bin/bash
#SBATCH -p vmcore128                  # 指定作业提交的分区 (队列)
#SBATCH -J VEP_test                    # 作业名称 (例如: VEP_snps, VEP_indels)
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB9185_CHO/03_logs/04_vep_annotation/out_file/vep_annotate_test.out  # 标准输出重定向
#SBATCH -e /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB9185_CHO/03_logs/04_vep_annotation/err_file/vep_annotate_test.err  # 错误输出重定向
#SBATCH --nodes=1                         # 作业申请 1 个节点
#SBATCH --ntasks-per-node=1               # 单节点启动 1 个任务
#SBATCH --cpus-per-task=8    # 使用为VEP配置的线程数
#SBATCH --mem=100G               # 申请的内存大小 (例如: 100G)
#SBATCH --time=120:00:00                 # 任务运行最长时间 (VEP注释可能耗时较长)

set -eo pipefail # 遇到错误立即退出，管道中的任何命令失败也退出
PRJ_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB9185_CHO"
# --- 作业执行内容 ---
LOGFILE="${PRJ_DIR}/01_results/08_vep_annotated/vep_debug_final_round/vcf_chunks/problem_line.log"
touch "${LOGFILE}"
exec > "${LOGFILE}" 2>&1 # 将所有标准输出和标准错误重定向到作业日志文件

echo "Job started on: $(hostname)"
echo "Start time: " && date
echo "Annotating snps with VEP..."

# 激活 mamba 环境 (VEP专用环境)
source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate "vep_115"

# 定义路径 
TEST_VCF="${PRJ_DIR}/01_results/08_vep_annotated/vep_debug_final_round/all_problematic_vcf_lines.vcf"
TEST_VEP_VCF="${PRJ_DIR}/01_results/08_vep_annotated/vep_debug_final_round/all_problematic_vcf_lines.vep.vcf.gz"

# 打印检测的行到日志
echo "当前检测的vcf位点: "
cat $TEST_VCF | grep -v '^#' 

# Run VEP annotation
echo "--- Running VEP annotation ---"
# --input_file: 输入VCF文件 (使用标准化后的VCF)
# --output_file: 输出注释后的VCF文件
# --cache: 使用本地缓存
# --dir_cache: 指定VEP缓存目录
# --assembly: 指定基因组装配版本
# --species: 指定物种拉丁名
# --fasta: 指定参考基因组FASTA文件 (VEP需要此文件来生成HGVS名称和检查参考序列)
# --format vcf: 输入文件格式为VCF
# --vcf: 输出文件格式为VCF
# --compress_output bgzip: 输出文件使用bgzip压缩
# --force_overwrite: 如果输出文件已存在则强制覆盖
# --everything: 启用所有常用注释选项 (包括Sift, PolyPhen, HGVS, Symbol, Gene biotype, AF等)
# --pick: 为每个变异选择一个“最严重”的后果 (通常是首选，简化输出)
# --fork: 启用多线程处理，提高速度
# --buffer_size: 内部缓冲区大小，可根据内存调整，默认5000通常足够
# --check_ref: 检查输入VCF中的参考等位基因是否与FASTA文件中的序列匹配
# --warning_file: 将警告信息写入单独的文件
# --stats_file --stats_html: 生成HTML格式的统计报告
vep \
    --input_file "${TEST_VCF}" \
    --output_file "${TEST_VEP_VCF}" \
    --offline --cache --merged\
    --dir_cache "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/vep_cache" \
    --assembly "CriGri-PICRH-1.0" \
    --species "cricetulus_griseus_picr" \
    --fasta "/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/CriGri-PICRH-1.0/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna" \
    --format vcf \
    --vcf \
    --compress_output bgzip \
    --force_overwrite \
    --everything \
    --pick \
    --fork 8 \
    --buffer_size 5000 \
    --check_ref \
    --warning_file "${PRJ_DIR}/03_logs/04_vep_annotation/vep_annotate_test.warnings.txt" \
    --stats_file "${PRJ_DIR}/03_logs/04_vep_annotation/vep_annotate_test.stats.html" \
    --stats_html


# Index the annotated VCF file
echo "--- Step 3: Indexing annotated VCF file ---"
tabix -p vcf "${TEST_VEP_VCF}"

echo "End time: " && date
echo "Job finished."
