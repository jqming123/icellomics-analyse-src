#!/bin/bash
#SBATCH -p vmcore128             
#SBATCH -J VEP_test              
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB39258_CHO/03_logs/out_file/filter_problematic_vcf.out  # 标准输出重定向
#SBATCH -e /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB39258_CHO/03_logs/err_file/filter_problematic_vcf.err  # 错误输出重定向
#SBATCH --nodes=1                
#SBATCH --ntasks-per-node=1      
#SBATCH --cpus-per-task=8    
#SBATCH --mem=100G               
#SBATCH --time=120:00:00        

PRJ_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/genome_projects/PRJEB39258_CHO"
LOGFILE="${PRJ_DIR}/01_results/08_vep_annotated/filter_problematic_vcf.log"
touch "${LOGFILE}"
exec > "${LOGFILE}" 2>&1 # 将所有标准输出和标准错误重定向到作业日志文件

# --- 配置 ---
# 原始的、有问题的 VCF 文件
SOURCE_VCF="${PRJ_DIR}/01_results/08_vep_annotated/main_chrs_snps.normalized.vcf.gz"

# 定义一个文件，其中包含要移除的完整 VCF 行。
# 脚本将从这些行中提取 CHROM 和 POS 信息。
# 示例:
# NC_048595.1 133527410 . G A ... (完整的VCF行)
# NC_048596.1 50608901 . T A ... (完整的VCF行)
PROBLEMATIC_VCF_FILE="${PRJ_DIR}/01_results/08_vep_annotated/vep_debug_final_round/all_problematic_vcf_lines.vcf"

# 定义一个新的、过滤后的输出文件路径
FILTERED_VCF="${PRJ_DIR}/01_results/08_vep_annotated/main_chrs_snps.normalized.filtered.vcf.gz"

# --- 检查先决条件 ---
if [ ! -f "${PROBLEMATIC_VCF_FILE}" ]; then
    echo "错误: 找不到有问题的 VCF 行文件: ${PROBLEMATIC_VCF_FILE}"
    echo "请创建此文件，并确保其包含要移除的完整 VCF 行。"
    exit 1
fi
if [ ! -f "${SOURCE_VCF}" ]; then
    echo "错误: 找不到源 VCF 文件: ${SOURCE_VCF}"
    exit 1
fi

# --- 执行命令 ---
echo "正在过滤 VCF 文件，移除有问题的行..."
echo "将使用文件 ${PROBLEMATIC_VCF_FILE} 中提取的 CHROM 和 POS 信息进行过滤。"

zcat "${SOURCE_VCF}" | \
    awk -v problematic_vcf_file="${PROBLEMATIC_VCF_FILE}" '
        BEGIN {
            # 在 awk 开始处理主 VCF 文件之前，读取 problematic_vcf_file
            # 将每一行的 CHROM_POS 组合作为键存储在 associative array (problematic) 中
            while ((getline < problematic_vcf_file) > 0) {
                # 忽略 problematic_vcf_file 中的头部行
                if ($0 ~ /^#/) {
                    continue
                }
                # 对于数据行，提取第1列 (CHROM) 和第2列 (POS)
                problematic[$1"_"$2] = 1
            }
            close(problematic_vcf_file) # 关闭 problematic_vcf_file
        }
        # 如果是主 VCF 文件的头部行 (以 # 开头)，直接打印并跳到下一行
        /^#/ {
            print
            next
        }
        # 对于主 VCF 文件的数据行，检查其 CHROM_POS 是否在 problematic 数组中
        {
            key = $1"_"$2 # 构建当前行的 CHROM_POS 键
            if (!(key in problematic)) {
                # 如果当前行的 CHROM_POS 不在 problematic 数组中，则打印该行
                print
            }
        }
    ' | \
    bgzip -c > "${FILTERED_VCF}"

# 检查上一个命令的退出状态
if [ $? -eq 0 ]; then
    echo "过滤完成！新的文件保存在: ${FILTERED_VCF}"
    echo "现在您可以用这个新文件重新运行 VEP。"
else
    echo "过滤过程中发生错误。请检查日志文件以获取详细信息。"
fi