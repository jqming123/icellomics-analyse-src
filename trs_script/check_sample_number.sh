#!/bin/bash

# 定义基础路径
# BASE_DIR="/gpfs/zhaowm_group/gaoxiaojing/CellLine/rna_shell/"
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/transcriptome_projects/EQ_jobs"


# 定义要检查的子目录列表
SUBDIRS=(
    "PRJNA378247"
    "PRJNA319417"
    "PRJNA316065"
    "PRJNA304606"
    "PRJNA255418"
    "PRJNA1123849"
    "PRJNA1088164"
    "PRJEB55916"
    "PRJEB48931"
    "PRJEB38542"
    "PRJEB37009"
    "PRJEB33024"
    "PRJEB30364"
)

echo "正在查看 ${BASE_DIR} 下指定子目录中的 .sh 文件数量..."
echo "-----------------------------------------------------"

# 遍历每个子目录
for dir_name in "${SUBDIRS[@]}"; do
    # 构建完整的目录路径
    FULL_PATH="${BASE_DIR}${dir_name}"

    # 检查目录是否存在
    if [ -d "$FULL_PATH" ]; then
        # 查找当前目录下的所有 .sh 文件并计数
        count=$(find "$FULL_PATH" -maxdepth 1 -type f -name "*.sh" | wc -l)
        echo "目录: ${dir_name} - .sh 文件数量: ${count}"
    else
        echo "目录: ${dir_name} - 不存在或不是一个目录。"
    fi
done

echo "-----------------------------------------------------"
echo "完成。"
