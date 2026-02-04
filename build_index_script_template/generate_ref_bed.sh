#!/usr/bin/env bash

REF_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/ref_genome"
REF_DIR="${REF_ROOT}/hg38_Ensemble"
REF_GTF="${REF_DIR}/Homo_sapiens.GRCh38.115.gtf"
REF_BED="${REF_DIR}/Homo_sapiens.GRCh38.115.for_scRNA.ref.bed"

gtf2bed < ${REF_GTF} > ${REF_BED}