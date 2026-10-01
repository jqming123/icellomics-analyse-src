#!/bin/bash

# ==================================
# Epigenetics Project Configuration
# ==================================

# --- Project Info ---
# PROJECT_NAME is optional. Project-scoped scripts should export it before
# sourcing this file. Assembly-only scripts may omit it.
if [ -n "${PROJECT_NAME:-}" ]; then
    if ! [[ "${PROJECT_NAME}" =~ ^PRJ[A-Z0-9]+_[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        echo "Error: PROJECT_NAME must match <bioproject_id>_<cellline>, e.g. PRJNA728969_HEK293. Current value: '${PROJECT_NAME}'" >&2
        exit 1
    fi
else
    echo "Notice: PROJECT_NAME is not set; current script does not use PROJECT_NAME." >&2
fi

if [ -z "${REF_NAME:-}" ]; then
    echo "Error: REF_NAME is not set. Export REF_NAME before sourcing epi_config.sh, e.g. export REF_NAME=\"hg38_Ensembl\"" >&2
    exit 1
fi

# --- Base Directories ---
BASE_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects"
RESOURCES_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources"
JBROWSE_DATA_ROOT="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects/jbrowse2_data"
export JBROWSE_DATA_ROOT

# --- Project-specific Paths ---
if [ -n "${PROJECT_NAME:-}" ]; then
    PROJECT_DIR="${BASE_DIR}/${PROJECT_NAME}"
    export PROJECT_DIR
    export TMP_DIR="${PROJECT_DIR}/tmp"
fi

# --- Reference Genome Configuration ---
# Configure reference-specific paths and MACS3 genome size values.
unset QC_THRESHOLD_PROFILE
case "${REF_NAME}" in
#    "CriGri-PICRH-1.0")
#        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0"
#        export BOWTIE2_INDEX="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic"
#        export REF_GENOME="${REF_DIR}/GCF_003668045.3_CriGri-PICRH-1.0_genomic.fna"
#        export GSIZE="2366634374"
#        ;;

    "CH_Ensembl")
        REF_DIR="${RESOURCES_DIR}/ref_genome/CriGri-PICRH-1.0_Ensembl"
        export BOWTIE2_INDEX="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Cricetulus_griseus_picr.CriGri-PICRH-1.0.dna.toplevel.fa"
        export GSIZE="2312971620"
        export SPECIES_ID="1"
        export ASSEMBLY_NAME="CriGri-PICRH-1.0"
        export MITO_CONTIGS="MT"
        ;;

    "hg38_Ensembl")
        REF_DIR="${RESOURCES_DIR}/ref_genome/hg38_Ensembl"
        export BOWTIE2_INDEX="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly"
        export REF_GENOME="${REF_DIR}/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        export GSIZE="hs"
        export SPECIES_ID="3"
        export ASSEMBLY_NAME="GRCh38.p14"
        export MITO_CONTIGS="MT"
        # This project uses Ensembl TSS annotations, so ENCODE's RefSeq-calibrated
        # TSS threshold is deliberately not applied yet.
        export QC_THRESHOLD_PROFILE="human_grch38_ensembl_phase23"
        ;;

    "Cattle_E_ARSUCD2")
        REF_DIR="${RESOURCES_DIR}/ref_genome/Cattle_E_ARSUCD2"
        export BOWTIE2_INDEX="${REF_DIR}/Bos_taurus.ARS-UCD2.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Bos_taurus.ARS-UCD2.0.dna.toplevel.fa"
        export GSIZE="2669382403"
        export SPECIES_ID="4"
        export ASSEMBLY_NAME="ARS-UCD2.0"
        export MITO_CONTIGS="MT"
        ;;

    "Chicken_E_GRCg7b")
        REF_DIR="${RESOURCES_DIR}/ref_genome/Chicken_E_GRCg7b"
        export BOWTIE2_INDEX="${REF_DIR}/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa"
        export GSIZE="1025963651"
        export SPECIES_ID="6"
        export ASSEMBLY_NAME="bGalGal1.mat.broiler.GRCg7b"
        export MITO_CONTIGS="MT"
        ;;

    "Dog_E_UUGSD")
        REF_DIR="${RESOURCES_DIR}/ref_genome/Dog_E_UUGSD"
        export BOWTIE2_INDEX="${REF_DIR}/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Canis_lupus_familiarisgsd.UU_Cfam_GSD_1.0.dna.toplevel.fa"
        export GSIZE="2372031271"
        export SPECIES_ID="5"
        export ASSEMBLY_NAME="UU_Cfam_GSD_1.0"
        export MITO_CONTIGS="MT"
        ;;

    "GreenMonkey_E_ChlSab1.1")
        REF_DIR="${RESOURCES_DIR}/ref_genome/GreenMonkey_E_ChlSab1.1"
        export BOWTIE2_INDEX="${REF_DIR}/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Chlorocebus_sabaeus.ChlSab1.1.dna.toplevel.fa"
        export GSIZE="2767805853"
        export SPECIES_ID="2"
        export ASSEMBLY_NAME="ChlSab1.1"
        export MITO_CONTIGS="MT"
        ;;

    "Mouse_E_GRCm39")
        REF_DIR="${RESOURCES_DIR}/ref_genome/Mouse_E_GRCm39"
        export BOWTIE2_INDEX="${REF_DIR}/Mus_musculus.GRCm39.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Mus_musculus.GRCm39.dna.toplevel.fa"
        export GSIZE="2495461690"
        export SPECIES_ID="7"
        export ASSEMBLY_NAME="GRCm39"
        export MITO_CONTIGS="MT"
        ;;

    "Pig_E_Sscrofa11.1")
        REF_DIR="${RESOURCES_DIR}/ref_genome/Pig_E_Sscrofa11.1"
        export BOWTIE2_INDEX="${REF_DIR}/Sus_scrofa.Sscrofa11.1.dna.toplevel"
        export REF_GENOME="${REF_DIR}/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
        export GSIZE="2455392186"
        export SPECIES_ID="8"
        export ASSEMBLY_NAME="Sscrofa11.1"
        export MITO_CONTIGS="MT"
        ;;

    "dont_need_ref")
        REF_DIR="NONE"
        export BOWTIE2_INDEX="NONE"
        export REF_GENOME="NONE"
        export GSIZE="0"
        export SPECIES_ID="0"
        export ASSEMBLY_NAME="NONE"
        export MITO_CONTIGS=""
        ;;

    *)
        echo "Error: unsupported REF_NAME '${REF_NAME}'." >&2
        echo "Supported options: CH_Ensembl, hg38_Ensembl, Cattle_E_ARSUCD2, Chicken_E_GRCg7b, Dog_E_UUGSD, GreenMonkey_E_ChlSab1.1, Mouse_E_GRCm39, Pig_E_Sscrofa11.1, dont_need_ref" >&2
        exit 1
        ;;
esac

# Non-human assemblies use diagnostic phase-1 thresholds until enough
# species-specific projects have been collected for calibration.
if [ -z "${QC_THRESHOLD_PROFILE:-}" ]; then
    export QC_THRESHOLD_PROFILE="custom_nonhuman_phase23"
fi

# --- QC Reference Resources ---
if [ "${REF_DIR}" != "NONE" ]; then
    export TSS_BED="${REF_DIR}/tss_bed/${REF_NAME}.tss.1bp.bed"
else
    export TSS_BED=""
fi

# Blacklists/exclusion sets are assembly-specific. Do not lift over a list from
# another assembly or silently describe an unfiltered peak set as filtered.
# A configured blacklist is mandatory only for a release-level analysis.  Core
# QC and existing-result QC must still be able to report their recoverable
# metrics when the resource or release-only outputs are unavailable.
export BLACKLIST_POLICY="release-required"
export BLACKLIST_BED=""
export BLACKLIST_STATUS="UNAVAILABLE_FOR_ASSEMBLY"
export BLACKLIST_SOURCE=""
case "${REF_NAME}" in
    "hg38_Ensembl")
        export BLACKLIST_BED="${RESOURCES_DIR}/blacklist_bed/hg38.blacklist.ensembl.bed.gz"
        export BLACKLIST_STATUS="CONFIGURED"
        export BLACKLIST_SOURCE="ENCODE ENCFF356LFX GRCh38 unified excludable regions; Ensembl contig naming"
        ;;
    "Mouse_E_GRCm39")
        export BLACKLIST_BED="${RESOURCES_DIR}/blacklist_bed/mm39.excluderanges.ensembl.bed.gz"
        export BLACKLIST_STATUS="CONFIGURED"
        export BLACKLIST_SOURCE="Ogata et al. 2023 mm39.excluderanges; AnnotationHub AH107321"
        ;;
    "dont_need_ref")
        export BLACKLIST_STATUS="NOT_APPLICABLE"
        ;;
    *)
        # No authoritative assembly-matched pre-generated exclusion set was
        # identified for the current assembly. Peak calling is allowed, but
        # provenance must state that no blacklist was applied.
        ;;
esac

# --- Execution scope and cost controls ---
# These are defaults rather than hard-coded behaviour in individual jobs.  A
# caller may override them in its environment when an explicitly fuller
# diagnostic run is wanted.
export REPRODUCIBILITY_MODE="${REPRODUCIBILITY_MODE:-auto}"
export PEAK_CAP="${PEAK_CAP:-500000}"
export ATAC_TSS_ARTIFACTS="${ATAC_TSS_ARTIFACTS:-profile}"

# --- Software & Environment ---
CONDA_PROFILE_PATH="/hpcdisk1/zhaowm_group/gaoxiaojing/softwares/miniforge3/etc/profile.d/conda.sh"
EPI_CONDA_ENV_NAME="ATAC_E4"

# --- Analysis Scripts ---
EPI_SCRIPT_DIR="/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/resources/src/epi_script"

# --- SLURM Resource Configuration ---
THREADS=16
MEM_SUPERLARGE="120G"
MEM_LARGE="64G"
MEM_MEDIUM="32G"
MEM_SMALL="16G"
# QC jobs compute the TSS offset histogram in-process (atac_qc/tss_qc.py), so
# the TSS stage no longer streams the pooled tagAlign into `bedtools coverage`
# and the old ~0.32G-per-1e6-records scaling no longer applies.  Keep MEM_QC
# generous anyway: it also covers samtools flagstat on legacy raw BAMs and
# in-memory fragment histograms.
MEM_QC="${MEM_QC:-160G}"
# Walltime for generated QC jobs.  Per-project overrides are allowed, e.g. for
# the largest legacy tagAligns that need many hours per biosample:
#   MEM_QC=240G QC_WALLTIME=7-00:00:00 bash 6_generate_atac_qc_job.sh ...
QC_WALLTIME="${QC_WALLTIME:-1-00:00:00}"

QUEUE_NAME="corexd192"
# QUEUE_NAME="core56"
# QUEUE_NAME="vmcore128"
