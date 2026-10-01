#!/usr/bin/env bash
#SBATCH -J PRJNA593571_HEK293_cycle
#SBATCH -p core56
#SBATCH -c 12
#SBATCH --mem=128G
#SBATCH -t 48:00:00
#SBATCH -o /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/scRNA_projects/HEK293/PRJNA593571_HEK293/5_logs/PRJNA593571_HEK293_cycle_analysis.%j.log

set -eo pipefail

PROJECT_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/scRNA_projects/HEK293/PRJNA593571_HEK293"
SCRIPT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/scRNA_script/LTC/cycle_analysis"
CONDA_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3"

source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "single_cell_1"
echo "Successfully activated single_cell_1"

mkdir -p "${PROJECT_PATH}/annotation" "${PROJECT_PATH}/5_logs"

Rscript "${SCRIPT_DIR}/scRNA_downstream_cycle_analysis.R" \
  --project_path "${PROJECT_PATH}" \
  --genome_name "hg38_Ensembl" \
  --out_dir "annotation/preprocessing" \
  --annotation_script_dir "${SCRIPT_DIR}" \
  --annotation_output_dir "annotation/functional_state_results" \
  --annotation_ncores 12 \
  --run_functional_annotation TRUE \
  --species "Homo sapiens" \
  --db_species "HS" \
  --mt_pattern "^MT-" \
  --run_scrublet TRUE \
  --convert_to_h5ad TRUE \
  --resolution 0.8 \
  --nfeatures 2000 \
  --seed 1234 \
  --nFeature_min 200 \
  --nFeature_max 8000 \
  --nCount_min 500 \
  --nCount_max 50000 \
  --mt_threshold 15

Rscript "${SCRIPT_DIR}/compare_ucell_matched_rank.R" \
  --input-rds "${PROJECT_PATH}/annotation/functional_state_results/seurat/functional_state_scored_seurat.rds" \
  --output-dir "${PROJECT_PATH}/annotation/functional_state_results/ucell_rank_comparison" \
  --manifest "${PROJECT_PATH}/annotation/functional_state_results/tables/scoring_parameters.csv" \
  --cluster-col seurat_clusters \
  --ncores 12

echo "Completed: PRJNA593571_HEK293 QC, clustering and functional-state profiling"
echo "Results: ${PROJECT_PATH}/annotation"
