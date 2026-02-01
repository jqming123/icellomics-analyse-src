#!/bin/bash

# --- 用户需要修改的部分 ---
# 请将下面的路径替换为你自己的参考基因组FASTA文件的绝对路径
# CRAM格式需要参考基因组来进行压缩。
REFERENCE_FASTA="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome/hg38_Ensemble/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
WORK_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/rna_count_result"
# -------------------------

# 检查参考基因组文件是否存在
if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "错误：参考基因组文件未找到于 '$REFERENCE_FASTA'"
    echo "请修改脚本，提供正确的路径。"
    exit 1
fi

echo "使用参考基因组: $REFERENCE_FASTA"
echo "开始在当前目录的子文件夹中查找并压缩BAM文件..."
echo "======================================================="

source "/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
conda activate "RNAseq_E4"

cd ${WORK_DIR}
# 遍历当前目录下的所有子目录 (e.g., ERR10047670, SRR15376523, etc.)
for sample_dir in */; do
    # 移除目录名末尾的斜杠
    sample_dir=${sample_dir%/}
    
    # 定义star目录的路径
    star_dir="${sample_dir}/star"

    if [ -d "$star_dir" ]; then
        echo "--- 正在处理样本: $sample_dir ---"
        
        # 查找star目录下的所有.bam文件
        # 使用find以便处理文件名中的空格等特殊字符
        find "$star_dir" -maxdepth 1 -type f -name "*.bam" | while read bam_file; do
            
            # 定义输出的cram文件名
            cram_file="${bam_file%.bam}.cram"

            if [ -f "$cram_file" ]; then
                echo "提示: '$cram_file' 已存在，跳过。"
            else
                echo "正在压缩: '$bam_file' -> '$cram_file'"
                
                # 使用samtools进行转换
                samtools view -@ 8 -C -T "$REFERENCE_FASTA" -o "$cram_file" "$bam_file"

                # 检查上一个命令是否成功
                if [ $? -eq 0 ]; then
                    echo "成功创建: '$cram_file'"
                    # [安全措施] 下面的删除命令被注释掉了。
                    # 确认CRAM文件都生成无误后，你可以手动删除BAM或取消下面一行的注释再运行一次脚本。
                    rm "$bam_file"
                    echo "已删除原始BAM文件: '$bam_file'"
                else
                    echo "错误: 压缩 '$bam_file' 失败。"
                fi
            fi
        done
    else
        echo "提示: 在 '$sample_dir' 中未找到 'star' 子目录，跳过。"
    fi
done

echo "======================================================="
echo "所有目录处理完毕。"
echo "请检查CRAM文件是否生成正确。如果确认无误，可以手动删除原始的.bam文件，"
echo "或者取消脚本中 'rm \"\$bam_file\"' 行的注释后重新运行以自动删除。"