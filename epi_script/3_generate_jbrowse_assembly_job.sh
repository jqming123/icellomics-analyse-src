#!/usr/bin/env bash

# @File       :3_generate_jbrowse_assembly_job.sh
# @Description:Generate one SLURM job script that prepares a JBrowse 2 assembly
#              bundle from an existing reference-genome directory.
# @Usage      :bash 3_generate_jbrowse_assembly_job.sh <REF_NAME>
# @Example    :bash 3_generate_jbrowse_assembly_job.sh hg38_Ensembl

if [ "$#" -ne 1 ]; then
    echo "Error: invalid argument count."
    echo "Usage: bash 3_generate_jbrowse_assembly_job.sh <REF_NAME>"
    exit 1
fi

export REF_NAME="$1"

CONFIG_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script/epi_config.sh"
if [ ! -f "${CONFIG_PATH}" ]; then
    echo "Error: config file not found at ${CONFIG_PATH}"
    exit 1
fi
source "${CONFIG_PATH}"

if [ -z "${REF_DIR:-}" ] || [ "${REF_DIR}" = "NONE" ]; then
    echo "Error: REF_DIR is not available for REF_NAME='${REF_NAME}' in ${CONFIG_PATH}"
    exit 1
fi

if [ -z "${JBROWSE_DATA_ROOT:-}" ]; then
    echo "Error: JBROWSE_DATA_ROOT is not defined in ${CONFIG_PATH}"
    exit 1
fi

if [ -z "${SPECIES_ID:-}" ] || [ -z "${ASSEMBLY_NAME:-}" ] || [ "${ASSEMBLY_NAME}" = "NONE" ]; then
    echo "Error: SPECIES_ID/ASSEMBLY_NAME are not available for REF_NAME='${REF_NAME}' in ${CONFIG_PATH}"
    exit 1
fi

if [ ! -d "${REF_DIR}" ]; then
    echo "Error: reference directory not found: ${REF_DIR}"
    exit 1
fi

JOB_DIR="${REF_DIR}/jbrowse2/2_jobs/assembly"
LOG_DIR="${REF_DIR}/jbrowse2/3_logs/assembly"
mkdir -p "${JOB_DIR}" "${LOG_DIR}"

JOB_SCRIPT_PATH="${JOB_DIR}/prepare_jbrowse_assembly_${REF_NAME}.sh"

cat > "${JOB_SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=jb2asm_${SPECIES_ID}
#SBATCH --partition=${QUEUE_NAME}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=${MEM_LARGE}
#SBATCH --time=1-00:00:00
#SBATCH --output=${LOG_DIR}/prepare_jbrowse_assembly_${REF_NAME}_%j.log

set -eo pipefail

echo "=========================================================="
echo "Job started on \$(date)"
echo "Reference name: ${REF_NAME}"
echo "Species ID: ${SPECIES_ID}"
echo "Assembly name: ${ASSEMBLY_NAME}"
echo "=========================================================="

export REF_NAME="${REF_NAME}"
source "${CONFIG_PATH}"

source "\${CONDA_PROFILE_PATH}"
conda activate "\${EPI_CONDA_ENV_NAME}"

REF_DIR="${REF_DIR}"
JBROWSE_DATA_ROOT="${JBROWSE_DATA_ROOT}"
ASSEMBLY_DIR="\${JBROWSE_DATA_ROOT}/assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}"

require_cmd() {
    command -v "\$1" >/dev/null 2>&1 || {
        echo "Error: required command not found: \$1" >&2
        exit 1
    }
}

# 智能识别基因组主 FASTA 文件（规避 gene, transcript, cdna 等干扰）
pick_genome_fasta() {
    local ext="\$1"
    local pattern file

    # 1. 优先尝试已知的标准基因组命名特征
    for pattern in "dna.toplevel" "dna_sm.primary_assembly" "dna.primary_assembly" "primary_assembly"; do
        for file in "\${REF_DIR}"/*"\${pattern}.\${ext}"; do
            if [ -f "\$file" ]; then
                printf '%s\n' "\$file"
                return 0
            fi
        done
    done

    # 2. 兜底匹配：查找任何 *.fa/*.fa.gz，但必须排除非基因组文件
    for file in "\${REF_DIR}"/*."\${ext}"; do
        if [ -f "\$file" ]; then
            case "\$(basename "\$file")" in
                *.gene.*|*.transcript.*|*.cdna.*)
                    continue
                    ;;
                *)
                    printf '%s\n' "\$file"
                    return 0
                    ;;
            esac
        fi
    done
    return 1
}

# 智能选择 GTF 文件（过滤 .bak 和临时文件）
pick_gtf_file() {
    local ext="\$1"
    local file
    for file in "\${REF_DIR}"/*."\${ext}"; do
        if [ -f "\$file" ]; then
            case "\$(basename "\$file")" in
                *.bak|*tmp*|*temp*)
                    continue
                    ;;
                *)
                    printf '%s\n' "\$file"
                    return 0
                    ;;
            esac
        fi
    done
    return 1
}

require_cmd samtools
require_cmd bgzip
require_cmd tabix
require_cmd gzip
require_cmd python
require_cmd sort

mkdir -p "\${ASSEMBLY_DIR}"
TMP_DIR=\$(mktemp -d)
trap 'rm -rf "\${TMP_DIR}"' EXIT

FASTA_GZ=\$(pick_genome_fasta "fa.gz") || true
FASTA_PLAIN=\$(pick_genome_fasta "fa") || true

if [ -n "\${FASTA_GZ:-}" ]; then
    if [ ! -f "\${FASTA_GZ}.fai" ] || [ ! -f "\${FASTA_GZ}.gzi" ]; then
        echo "Error: bgzip FASTA indexes are incomplete for \${FASTA_GZ}" >&2
        exit 1
    fi
    ln -sfn "\${FASTA_GZ}" "\${ASSEMBLY_DIR}/genome.fa.gz"
    ln -sfn "\${FASTA_GZ}.fai" "\${ASSEMBLY_DIR}/genome.fa.gz.fai"
    ln -sfn "\${FASTA_GZ}.gzi" "\${ASSEMBLY_DIR}/genome.fa.gz.gzi"
    cut -f1,2 "\${FASTA_GZ}.fai" > "\${ASSEMBLY_DIR}/chrom.sizes"
    SEQUENCE_ADAPTER="BgzipFastaAdapter"
    SEQUENCE_URI="assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genome.fa.gz"
    FAI_URI="assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genome.fa.gz.fai"
    GZI_URI="assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genome.fa.gz.gzi"
elif [ -n "\${FASTA_PLAIN:-}" ]; then
    if [ ! -f "\${FASTA_PLAIN}.fai" ]; then
        samtools faidx "\${FASTA_PLAIN}"
    fi
    ln -sfn "\${FASTA_PLAIN}" "\${ASSEMBLY_DIR}/genome.fa"
    ln -sfn "\${FASTA_PLAIN}.fai" "\${ASSEMBLY_DIR}/genome.fa.fai"
    cut -f1,2 "\${FASTA_PLAIN}.fai" > "\${ASSEMBLY_DIR}/chrom.sizes"
    SEQUENCE_ADAPTER="IndexedFastaAdapter"
    SEQUENCE_URI="assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genome.fa"
    FAI_URI="assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genome.fa.fai"
    GZI_URI=""
else
    echo "Error: no genome FASTA file was found under \${REF_DIR}" >&2
    exit 1
fi

GTF_GZ=\$(pick_gtf_file "gtf.gz") || true
GTF_PLAIN=\$(pick_gtf_file "gtf") || true

if [ -n "\${GTF_GZ:-}" ]; then
    gzip -dc "\${GTF_GZ}" > "\${TMP_DIR}/genes.input.gtf"
elif [ -n "\${GTF_PLAIN:-}" ]; then
    cp "\${GTF_PLAIN}" "\${TMP_DIR}/genes.input.gtf"
else
    echo "Error: no GTF file was found under \${REF_DIR}" >&2
    exit 1
fi

grep '^#' "\${TMP_DIR}/genes.input.gtf" > "\${TMP_DIR}/genes.header.gtf" || true
grep -v '^#' "\${TMP_DIR}/genes.input.gtf" | sort -t \$'\t' -k1,1 -k4,4n > "\${TMP_DIR}/genes.body.sorted.gtf"
cat "\${TMP_DIR}/genes.header.gtf" "\${TMP_DIR}/genes.body.sorted.gtf" | bgzip -f -c > "\${ASSEMBLY_DIR}/genes.gtf.gz"

tabix -f -p gff "\${ASSEMBLY_DIR}/genes.gtf.gz"

cat "\${TMP_DIR}/genes.header.gtf" "\${TMP_DIR}/genes.body.sorted.gtf" | python "\${EPI_SCRIPT_DIR}/gtf_to_gff3.py" - | bgzip -f -c > "\${ASSEMBLY_DIR}/genes.gff3.gz"

tabix -f -p gff "\${ASSEMBLY_DIR}/genes.gff3.gz"

cat > "\${ASSEMBLY_DIR}/assembly.meta.json" <<JSON
{
  "species_id": ${SPECIES_ID},
  "assembly_name": "${ASSEMBLY_NAME}",
  "ref_name": "${REF_NAME}",
  "sequence_adapter": "\${SEQUENCE_ADAPTER}",
  "annotation_adapter": "Gff3TabixAdapter",
  "sequence_uri": "\${SEQUENCE_URI}",
  "fai_uri": "\${FAI_URI}",
  "gzi_uri": "\${GZI_URI}",
  "gtf_uri": "assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genes.gtf.gz",
  "annotation_uri": "assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genes.gff3.gz",
  "annotation_index_uri": "assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/genes.gff3.gz.tbi",
  "chrom_sizes_uri": "assemblies/${SPECIES_ID}/${ASSEMBLY_NAME}/chrom.sizes"
}
JSON

set +u
conda deactivate

echo "=========================================================="
echo "Assembly bundle prepared at: \${ASSEMBLY_DIR}"
echo "Job finished on \$(date)"
echo "=========================================================="
EOF

chmod +x "${JOB_SCRIPT_PATH}"

echo "Done. Generated JBrowse 2 assembly job script:"
echo "  ${JOB_SCRIPT_PATH}"
echo ""
echo "Submit it with:"
echo "  sbatch ${JOB_SCRIPT_PATH}"
