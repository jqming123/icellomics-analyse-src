#!/usr/bin/env bash
# Generate upstream/downstream SLURM scripts (GENtoolkit)
# Layout: /hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/scRNA_projects/<BioProjectID>/
#   1_raw/ 2_output/ 3_expression_result/ 4_jobs/ 5_logs/ project.conf
#   slurm scripts stored in 4_jobs/, logs stored in 5_logs/
#
# Project config (project.conf, PROJ_*):
#   PROJ_LIB_TYPE, PROJ_REF_NAME, PROJ_READ_TYPE, PROJ_SEQ_TYPE,
#   PROJ_INDEX_BUILD, PROJ_DESIGNATED_ALL, PROJ_SAMPLE_LIST,
#   PROJ_CELLRANGER_LOCALCORES, PROJ_CELLRANGER_LOCALMEM,
#   PROJ_DROP_REPORT, PROJ_EXPR_DATA, PROJ_META_DATA, PROJ_REF_PATH, PROJ_WORK_PATH
#
# Shared config (config.sh):
#   PROJECT_ROOT, RESOURCES_ROOT, GEN_TOOLKIT_DIR, CONDA_SH, CONDA_ENV,
#   UP_CPUS, UP_MEM, UP_TIME, DOWN_CPUS, DOWN_MEM, DOWN_TIME, QUEUE_NAME,
#   WORK_PATH_DEFAULT_SUFFIX, REF_* index paths
#
# Usage:
#   bash generate_scRNA_slurm.sh <project_dir_name>
#
# Notes:
#   - <project_dir_name> is the subdirectory name under $PROJECT_ROOT
#     (e.g., PRJDB6793_CHO).
#   - This script generates upstream/downstream SLURM scripts for that
#     single project only.

set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SHARED_CONFIG="$SCRIPT_DIR/config.sh"
PROJECT_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/scRNA_projects"
PRJ_CONFIG_NAME="project.conf"

if [[ $# -lt 1 ]]; then
  echo "[ERROR] Missing project directory name." >&2
  echo "Usage: bash $0 <project_dir_name>" >&2
  exit 1
fi

if [[ ! -d "$PROJECT_ROOT" ]]; then
  echo "[ERROR] PROJECT_ROOT not found: $PROJECT_ROOT" >&2
  exit 1
fi

proj_name="$1"
proj_dir="$PROJECT_ROOT/$proj_name"

if [[ ! -d "$proj_dir" ]]; then
  echo "[ERROR] Project directory not found: $proj_dir" >&2
  exit 1
fi

if [[ ! -f "$SHARED_CONFIG" ]]; then
  echo "[ERROR] Shared config not found: $SHARED_CONFIG" >&2
  exit 1
fi

cfg="$proj_dir/$PRJ_CONFIG_NAME"

if [[ ! -f "$cfg" ]]; then
  echo "[ERROR] Missing $PRJ_CONFIG_NAME in $proj_dir" >&2
  exit 1
fi


# load project config
set -a
# shellcheck disable=SC1090
source "$cfg"
# load shared config (requires REF_NAME)
REF_NAME="$PROJ_REF_NAME"
export REF_NAME
# shellcheck disable=SC1090
source "$SHARED_CONFIG"
set +a

if [[ -z "${PROJ_REF_NAME:-}" ]]; then
  echo "[ERROR] PROJ_REF_NAME not set in $PRJ_CONFIG_NAME" >&2
  exit 1
fi

if [[ ! -f "$UP_PY" ]]; then
  echo "[ERROR] GENtoolkit upstream script not found: $UP_PY" >&2
  exit 1
fi

if [[ ! -f "$DOWN_PY" ]]; then
  echo "[ERROR] GENtoolkit downstream script not found: $DOWN_PY" >&2
  exit 1
fi

if [[ -z "${PROJ_LIB_TYPE:-}" ]]; then
  echo "[ERROR] PROJ_LIB_TYPE not set in $PRJ_CONFIG_NAME" >&2
  exit 1
fi

raw_dir="$proj_dir/1_raw"
slurm_job_dir="$proj_dir/4_jobs"
log_dir="$proj_dir/5_logs"
mkdir -p "$log_dir"
mkdir -p "$slurm_job_dir"

up_slurm="$slurm_job_dir/${proj_name}_upstream.sh"
down_slurm="$slurm_job_dir/${proj_name}_downstream.sh"

# optional sample list
sample_args=""
if [[ "$PROJ_DESIGNATED_ALL" == "Designated_samples" && -n "$PROJ_SAMPLE_LIST" ]]; then
  sample_args="-da Designated_samples -sl $PROJ_SAMPLE_LIST"
else
  sample_args="-da $PROJ_DESIGNATED_ALL"
fi

# upstream command
case "$PROJ_LIB_TYPE" in
  "10X")
    upstream_cmd=(
      "python" "$UP_PY"
      "-blt" "10X"
      "-rgf" "$REF_GENOME"
      "-rgg" "$REF_GTF"
      "-ci" "$CELLRANGER_INDEX"
      "-rd" "$raw_dir"
      "-rt" "$PROJ_READ_TYPE"
      "-ib" "$PROJ_INDEX_BUILD"
      "-clc" "$PROJ_CELLRANGER_LOCALCORES"
      "-clm" "$PROJ_CELLRANGER_LOCALMEM"
      "$sample_args"
    )
    ;;
  "Smart-seq2"|"Bulk")
    upstream_cmd=(
      "python" "$UP_PY"
      "-blt" "$PROJ_LIB_TYPE"
      "-rgf" "$REF_GENOME"
      "-rgg" "$REF_GTF"
      "-bf" "$REF_BED"
      "-hi" "$HISAT2_INDEX"
      "-ri" "$RSEM_INDEX"
      "-rd" "$raw_dir"
      "-rt" "$PROJ_READ_TYPE"
      "-st" "$PROJ_SEQ_TYPE"
      "-ib" "$PROJ_INDEX_BUILD"
      "-sp" "$STAR_PATH"
      "$sample_args"
    )
    ;;
  "Drop-seq"|"inDrop_v1"|"inDrop_v2"|"inDrop_v3")
    # Determine config xml
    drop_cfg=""
    case "$PROJ_LIB_TYPE" in
        "Drop-seq") drop_cfg="$GEN_TOOLKIT_DIR/config/drop_seq.xml" ;;
        "inDrop_v1"|"inDrop_v2") drop_cfg="$GEN_TOOLKIT_DIR/config/indrop_v1_2.xml" ;;
        "inDrop_v3") drop_cfg="$GEN_TOOLKIT_DIR/config/indrop_v3.xml" ;;
    esac

    upstream_cmd=(
      "python" "$UP_PY"
      "-blt" "$PROJ_LIB_TYPE"
      "-rgf" "$REF_GENOME"
      "-rgg" "$REF_GTF"
      "-dc" "$drop_cfg"
      "-rd" "$raw_dir"
      "-rt" "$PROJ_READ_TYPE"
      "-ib" "$PROJ_INDEX_BUILD"
      "$sample_args"
    )
    if [[ -n "${PROJ_DROP_REPORT:-}" ]]; then
      upstream_cmd+=("-dm" "$PROJ_DROP_REPORT")
    fi
    ;;
  *)
    echo "[ERROR] Unsupported PROJ_LIB_TYPE=$PROJ_LIB_TYPE" >&2
    exit 1
    ;;
esac

# downstream defaults
local_work_suffix="$WORK_PATH_DEFAULT_SUFFIX"
if [[ -z "$local_work_suffix" || "$local_work_suffix" == "03.expression/downstream" ]]; then
  local_work_suffix="3_expression_result/downstream"
fi
WORK_PATH="${PROJ_WORK_PATH:-$proj_dir/$local_work_suffix}"
if [[ -z "${PROJ_EXPR_DATA:-}" ]]; then
  if [[ "$PROJ_LIB_TYPE" == "10X" ]]; then
    EXPR_DATA="$proj_dir/3_expression_result/outs/filtered_feature_bc_matrix"
  else
    EXPR_DATA="$proj_dir/3_expression_result/Project_GeneMat_rawCounts.txt"
  fi
else
  EXPR_DATA="$PROJ_EXPR_DATA"
fi

# downstream command (GENtoolkit)
downstream_cmd=(
  "python" "$DOWN_PY"
  "--stream" "down"
  "--BuildLibraryType" "$PROJ_LIB_TYPE"
  "--workpath" "$WORK_PATH"
  "--exprData" "$EXPR_DATA"
)

if [[ -n "${PROJ_META_DATA:-}" ]]; then
  downstream_cmd+=("--metaData" "$PROJ_META_DATA")
fi

if [[ -n "${PROJ_REF_PATH:-}" ]]; then
  downstream_cmd+=("--refpath" "$PROJ_REF_PATH")
fi

# Format commands with newlines for readability
up_cmd_formatted=""
for i in "${!upstream_cmd[@]}"; do
  val="${upstream_cmd[$i]}"
  if [[ $i -eq 0 ]]; then
    up_cmd_formatted="$val"
  elif [[ $i -eq 1 ]]; then
    up_cmd_formatted+=" $val"
  elif [[ "$val" == -* ]]; then
    up_cmd_formatted+=" \\"$'\n'"  $val"
  else
    up_cmd_formatted+=" $val"
  fi
done

down_cmd_formatted=""
for i in "${!downstream_cmd[@]}"; do
  val="${downstream_cmd[$i]}"
  if [[ $i -eq 0 ]]; then
    down_cmd_formatted="$val"
  elif [[ $i -eq 1 ]]; then
    down_cmd_formatted+=" $val"
  elif [[ "$val" == -* ]]; then
    down_cmd_formatted+=" \\"$'\n'"  $val"
  else
    down_cmd_formatted+=" $val"
  fi
done

cat > "$up_slurm" <<EOF
#!/usr/bin/env bash
#SBATCH -J ${proj_name}_up
#SBATCH -p ${QUEUE_NAME}
#SBATCH -c ${UP_CPUS}
#SBATCH --mem=${UP_MEM}
#SBATCH -t ${UP_TIME}
#SBATCH -o ${log_dir}/${proj_name}_up.%j.log

set -euo pipefail

if [[ -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
  conda activate "$CONDA_ENV"
  echo "Successfully activate $CONDA_ENV"
fi

cd "$GEN_TOOLKIT_DIR"

$up_cmd_formatted
EOF

cat > "$down_slurm" <<EOF
#!/usr/bin/env bash
#SBATCH -J ${proj_name}_down
#SBATCH -p ${QUEUE_NAME}
#SBATCH -c ${DOWN_CPUS}
#SBATCH --mem=${DOWN_MEM}
#SBATCH -t ${DOWN_TIME}
#SBATCH -o ${log_dir}/${proj_name}_down.%j.log

set -euo pipefail

if [[ -f "$CONDA_SH" ]]; then
  source "$CONDA_SH"
  conda activate "$CONDA_ENV"
  echo "Successfully activate $CONDA_ENV"
fi

cd "$GEN_TOOLKIT_DIR"

$down_cmd_formatted
EOF

echo "[OK] $proj_name -> $up_slurm, $down_slurm"
