#!/bin/bash
# Generate a SLURM job script for running prepare_scrna_ds0004_import.py
# Usage:
#   bash generate_h5ad_parse_slurm_job.sh <project_cellline> [h5ad_path]
#
# The input h5ad is
#   <project>/annotation/preprocessing/h5ad/seurat_obj_functional_state_profiled.h5ad
# Projects that have not produced that file yet are SKIPPED (no dataset_id is
# allocated and no job script is generated; the script exits 0).
# Pass [h5ad_path] (or set H5AD_PATH_OVERRIDE) to force a specific file.
#
# Examples:
#   bash generate_h5ad_parse_slurm_job.sh PRJNA484547_HEK293
#   bash generate_h5ad_parse_slurm_job.sh PRJNA1008166_H9
#   bash generate_h5ad_parse_slurm_job.sh PRJNA781454_Hela
#
# The argument is parsed as <project_accession>_<cell_line>, e.g. PRJNA484547_HEK293
# -> project_accession=PRJNA484547, cell_line=HEK293
#
# singlecell_dataset_id is assigned per bioproject from bioproject_dataset_id_map.tsv
# Each bioproject gets a unique ID starting from 400001

set -euo pipefail

# ============================================================
# Path constants
# ============================================================
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine"
SCRIPTS_DIR="${BASE_DIR}/resources/src/scRNA_script"
IMPORT_SCRIPT="${SCRIPTS_DIR}/prepare_scrna_ds0004_import.py"
CELL_LINE_TSV="${BASE_DIR}/resources/cell_line_metadata.tsv"
DATASET_ID_MAP="${SCRIPTS_DIR}/bioproject_dataset_id_map.tsv"
H5AD_BASE="${BASE_DIR}/scRNA_projects"
OUTPUT_BASE="${BASE_DIR}/scRNA_projects/ds0004_import/parse_output"
JOB_OUTPUT_DIR="${BASE_DIR}/scRNA_projects/ds0004_import/slurm_jobs"
LOG_OUTPUT_DIR="${BASE_DIR}/scRNA_projects/ds0004_import/slurm_logs"

CONDA_INIT="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="single_cell_1"

# SLURM defaults
SLURM_PARTITION="corexd192"
SLURM_NODES="1"
SLURM_NTASKS="1"
SLURM_CPUS_PER_TASK="4"
SLURM_MEM="32G"
SLURM_TIME="24:00:00"

# ============================================================
# Parse input argument
# ============================================================
if [ $# -lt 1 ]; then
    echo "Usage: $0 <project_cellline> [h5ad_path]"
    echo "Example: $0 PRJNA484547_HEK293"
    echo "Example: $0 PRJNA593571_HEK293 /abs/path/to/seurat_obj_functional_state_profiled.h5ad"
    exit 1
fi

PROJECT_CELLLINE="$1"
# Optional override: 2nd positional argument, else H5AD_PATH_OVERRIDE env var
H5AD_PATH_OVERRIDE="${2:-${H5AD_PATH_OVERRIDE:-}}"

# Extract project_accession and cell_line_dir from the argument
# The format is <project_accession>_<cell_line_dir>
# e.g. PRJNA484547_HEK293 -> project_accession=PRJNA484547, cell_line_dir=HEK293
# e.g. PRJNA1008166_H9   -> project_accession=PRJNA1008166,   cell_line_dir=H9
# e.g. PRJNA781454_HeLa   -> project_accession=PRJNA781454,   cell_line_dir=HeLa
# e.g. PRJEB59449_H9      -> project_accession=PRJEB59449,    cell_line_dir=H9

# Some project accessions have underscores (rare), so match from the end
# Strategy: try PRJNA/PRJEB/SRP prefix, then the rest is cell_line
PROJECT_ACCESSION=$(echo "$PROJECT_CELLLINE" | grep -oP '^(PRJNA\d+|PRJEB\d+|SRP\d+)')

if [ -z "$PROJECT_ACCESSION" ]; then
    # Fallback: use everything before the last underscore as project_accession
    # But cell lines like HEK293, H9, Hela are single words without underscores
    CELL_LINE_DIR="${PROJECT_CELLLINE##*_}"
    PROJECT_ACCESSION="${PROJECT_CELLLINE%_*}"
else
    CELL_LINE_DIR="${PROJECT_CELLLINE#${PROJECT_ACCESSION}_}"
fi

# Map directory cell_line name to canonical name (for --cell-line parameter)
# The canonical names must match cell_line_name in resources/cell_line_metadata.tsv
case "$CELL_LINE_DIR" in
    HeLa|Hela)   CELL_LINE_CANONICAL="HeLa" ;;
    DF1|DF-1)    CELL_LINE_CANONICAL="DF-1" ;;
    WI38|WI-38)  CELL_LINE_CANONICAL="WI-38" ;;
    SP20|SP2/0)  CELL_LINE_CANONICAL="SP2/0" ;;
    MCR5|MRC-5)  CELL_LINE_CANONICAL="MRC-5" ;;
    *)           CELL_LINE_CANONICAL="$CELL_LINE_DIR" ;;
esac

# ============================================================
# Resolve input h5ad (functional-state-profiled result only)
#   Projects that have not produced
#     <project>/annotation/preprocessing/h5ad/seurat_obj_functional_state_profiled.h5ad
#   are SKIPPED (no dataset_id is allocated, no job script is generated).
#   An explicit path can be forced via the 2nd argument or H5AD_PATH_OVERRIDE.
# ============================================================
PROJECT_DIR="${H5AD_BASE}/${CELL_LINE_DIR}/${PROJECT_CELLLINE}"
FUNCTIONAL_H5AD="${PROJECT_DIR}/annotation/preprocessing/h5ad/seurat_obj_functional_state_profiled.h5ad"

if [ -n "${H5AD_PATH_OVERRIDE}" ]; then
    INPUT_H5AD="${H5AD_PATH_OVERRIDE}"
    if [ ! -f "${INPUT_H5AD}" ]; then
        echo "ERROR: override h5ad not found: ${INPUT_H5AD}"
        exit 1
    fi
elif [ -f "${FUNCTIONAL_H5AD}" ]; then
    INPUT_H5AD="${FUNCTIONAL_H5AD}"
else
    echo "SKIP: ${PROJECT_CELLLINE} has not produced the functional-state h5ad yet."
    echo "  Expected: ${FUNCTIONAL_H5AD}"
    echo "  No parsing job generated for this project."
    exit 0
fi

# ============================================================
# Assign unique singlecell_dataset_id per bioproject
# ============================================================
if [ ! -f "$DATASET_ID_MAP" ]; then
    # Create the mapping file with header if it doesn't exist
    echo "singlecell_dataset_id	project_accession	cell_line_name	status	created_at" > "$DATASET_ID_MAP"
fi

# Use flock for concurrent-safe access to the mapping file
LOCK_FILE="${SCRIPTS_DIR}/.bioproject_dataset_id_map.lock"
mkdir -p "${SCRIPTS_DIR}"

# Check if project_accession already has an ID in the mapping table

SINGLECELL_DATASET_ID=$(awk -F'\t' \
    -v proj="${PROJECT_ACCESSION}" \
    -v cell="${CELL_LINE_CANONICAL}" \
    'NR>1 && $2 == proj && $3 == cell {print $1; exit}' \
    "$DATASET_ID_MAP")

if [ -n "$SINGLECELL_DATASET_ID" ]; then
    echo "INFO: Found existing dataset_id=${SINGLECELL_DATASET_ID} for ${PROJECT_ACCESSION}_${CELL_LINE_CANONICAL}"
else
    # Assign new ID: find max ID and increment
    (
        flock -x 200
        MAX_ID=$(awk -F'\t' 'NR>1 && $1 ~ /^400[0-9]{3}$/ {if($1>max)max=$1} END{print max+0}' "$DATASET_ID_MAP")
        if [ "$MAX_ID" -lt 400001 ]; then
            MAX_ID=400000
        fi
        NEW_ID=$((MAX_ID + 1))
        CREATED_AT=$(date '+%Y-%m-%d')
        # Append new row to mapping file
        printf '%d\t%s\t%s\tpending\t%s\n' "$NEW_ID" "$PROJECT_ACCESSION" "$CELL_LINE_CANONICAL" "$CREATED_AT" >> "$DATASET_ID_MAP"
        echo "INFO: Assigned new dataset_id=${NEW_ID} for ${PROJECT_ACCESSION}_${CELL_LINE_CANONICAL}"
    ) 200>"$LOCK_FILE"

    # Read the newly assigned ID
    SINGLECELL_DATASET_ID=$(awk -F'\t' -v proj="${PROJECT_ACCESSION}" \
        'NR>1 && $2 == proj {print $1; exit}' "$DATASET_ID_MAP")
    # Read the newly assigned ID
    SINGLECELL_DATASET_ID=$(awk -F'\t' \
         -v proj="${PROJECT_ACCESSION}" \
         -v cell="${CELL_LINE_CANONICAL}" \
         'NR>1 && $2 == proj && $3 == cell {print $1; exit}' \
         "$DATASET_ID_MAP")

fi

OUTPUT_DIR="${OUTPUT_BASE}/${CELL_LINE_DIR}/${PROJECT_ACCESSION}"

# ============================================================
# Generate SLURM job script
# ============================================================
mkdir -p "$JOB_OUTPUT_DIR" "$LOG_OUTPUT_DIR" "$OUTPUT_DIR"

SLURM_SCRIPT="${JOB_OUTPUT_DIR}/run_parse_${PROJECT_CELLLINE}.sh"

cat > "$SLURM_SCRIPT" << SCRIPT_EOF
#!/bin/bash
#SBATCH -p ${SLURM_PARTITION}
#SBATCH -J parse_${PROJECT_CELLLINE}
#SBATCH -o ${LOG_OUTPUT_DIR}/parse_${PROJECT_CELLLINE}_%j.log
#SBATCH --nodes=${SLURM_NODES}
#SBATCH --ntasks-per-node=${SLURM_NTASKS}
#SBATCH --cpus-per-task=${SLURM_CPUS_PER_TASK}
#SBATCH --mem=${SLURM_MEM}
#SBATCH --time=${SLURM_TIME}

set -euo pipefail

echo "============================================================"
echo "SLURM Job: parse_${PROJECT_CELLLINE}"
echo "Job ID: \${SLURM_JOB_ID}"
echo "Node: \$(hostname)"
echo "Start: \$(date)"
echo "============================================================"

# Activate conda (set +u to work around conda deactivate script ZSH_VERSION issue)
set +u
source "${CONDA_INIT}"
conda activate "${CONDA_ENV}"
set -u

echo "Conda env: \$(conda info --envs | grep '*' || true)"
echo "Python: \$(which python)"

# Run import script
python "${IMPORT_SCRIPT}" \\
    --input-h5ad "${INPUT_H5AD}" \\
    --output-dir "${OUTPUT_DIR}" \\
    --singlecell-dataset-id ${SINGLECELL_DATASET_ID} \\
    --cell-line "${CELL_LINE_CANONICAL}" \\
    --project-accession "${PROJECT_ACCESSION}" \\
    --cell-line-tsv "${CELL_LINE_TSV}" \\
    --gene-symbol-source gene_name \\
    --raw-mode auto

exit_code=\$?

echo "============================================================"
echo "End: \$(date)"
echo "Exit code: \${exit_code}"
echo "Output dir: ${OUTPUT_DIR}"
echo "============================================================"

exit \${exit_code}
SCRIPT_EOF

chmod +x "$SLURM_SCRIPT"

echo "============================================================"
echo "Generated SLURM job script: ${SLURM_SCRIPT}"
echo "============================================================"
echo "  Project:          ${PROJECT_CELLLINE}"
echo "  Project Accession: ${PROJECT_ACCESSION}"
echo "  Cell Line (dir):   ${CELL_LINE_DIR}"
echo "  Cell Line (canonical): ${CELL_LINE_CANONICAL}"
echo "  Dataset ID:        ${SINGLECELL_DATASET_ID}"
echo "  Input h5ad:        ${INPUT_H5AD}"
echo "  Output dir:        ${OUTPUT_DIR}"
echo ""
echo "Submit with:  sbatch ${SLURM_SCRIPT}"
echo "============================================================"
