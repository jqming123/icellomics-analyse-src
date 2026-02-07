#!/usr/bin/env bash

REF_ROOT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
# 注意修改以下参数
GENOME_NAME=""
REF_GENOME_DIR="${REF_ROOT_PATH}/${GENOME_NAME}"
REF_GTF="${REF_GENOME_DIR}/Homo_sapiens.GRCh38.115.gtf"
REF_BED="${REF_GENOME_DIR}/Homo_sapiens.GRCh38.115.for_scRNA.ref.bed"

gtf2bed < ${REF_GTF} > ${REF_BED}