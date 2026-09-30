"""Common constants and formatting helpers for ATAC QC summaries."""

from __future__ import annotations

from pathlib import Path
from typing import Optional


DEFAULT_BASE_DIR = Path("/hpcdisk1/zhaowm_group/gaoxiaojing/CellLine/epigen_projects")
DEFAULT_MITO_CONTIGS = ("MT",)

QC_COLUMNS = [
    "project_name",
    "biosample_id",
    "biological_replicates",
    "biological_replicate_count",
    "legacy_sample_map",
    "qc_mode",
    "threshold_profile",
    "run_ids",
    "run_count",
    "layout",
    "fastp_raw_reads",
    "fastp_clean_reads",
    "picard_examined_reads",
    "picard_duplicate_reads",
    "picard_optical_duplicate_reads",
    "optical_duplicate_fraction",
    "duplication_fraction",
    "library_total_fragments",
    "library_distinct_fragments",
    "library_one_read_fragments",
    "library_two_read_fragments",
    "nrf",
    "pbc1",
    "pbc2",
    "library_complexity_file",
    "gc_bias_status",
    "raw_bam_total_reads",
    "raw_bam_mapped_reads",
    "mapping_rate",
    "raw_bam_paired_reads",
    "raw_bam_proper_pair_reads",
    "proper_pair_rate",
    "mapq_filter_threshold",
    "mapq_filter_evidence",
    "mapq_filter_status",
    "mapq_filter_operator",
    "mapq_excluded_sam_flags",
    "mapq_required_sam_flags",
    "mapq_raw_alignment_records",
    "mapq_post_filter_alignment_records",
    "mapq_post_final_filter_alignment_records",
    "mapq_filter_retention",
    "raw_mito_mapped_reads",
    "raw_mito_reads",
    "raw_mito_fraction",
    "nodup_pre_mito_reads",
    "final_nodup_reads",
    "final_usable_fraction",
    "mitochondrial_contig_reads",
    "mitochondrial_filter_removed_reads",
    "mitochondrial_filter_fraction",
    "tagalign_total_reads",
    "reads_in_peaks",
    "frip_fraction",
    "frip_percent",
    "fragment_total",
    "fragments_in_peaks",
    "fragment_frip_fraction",
    "frip_input_unit",
    "frip_qc_status",
    "frip_contig_intersection",
    "peak_count",
    "peak_total_bp",
    "peak_median_width",
    "frip_method",
    "fragment_count",
    "fragment_mean_size",
    "fragment_median_size",
    "fragment_p10_size",
    "fragment_p90_size",
    "nfr_fraction",
    "mono_nucleosome_fraction",
    "di_nucleosome_fraction",
    "fragment_periodicity_score",
    "fragment_histogram_file",
    "fragment_input_status",
    "blacklist_reads",
    "blacklist_fraction",
    "blacklist_qc_status",
    "blacklist_qc_scope",
    "blacklist_contig_intersection",
    "peak_provenance",
    "blacklist_peak_filter_status",
    "blacklist_config_status",
    "blacklist_config_bed",
    "blacklist_config_sha256",
    "reference_config_fai",
    "reference_config_fai_sha256",
    "peak_pipeline_mode",
    "peak_blacklist_applied",
    "peak_blacklist_status",
    "peak_blacklist_bed",
    "peak_blacklist_runtime_bed",
    "peak_blacklist_sha256",
    "peak_blacklist_source",
    "peak_blacklist_validation_status",
    "peak_blacklist_validation_record",
    "peak_blacklist_validation_bed",
    "peak_blacklist_validation_sha256",
    "peak_blacklist_validation_reference_fai",
    "peak_blacklist_validation_reference_fai_sha256",
    "peak_blacklist_validation_interval_count",
    "tss_enrichment_score",
    "tss_method",
    "tss_tagalign_shift_status",
    "tss_artifacts",
    "tss_profile_file",
    "tss_matrix_file",
    "tss_heatmap_file",
    "self_consistency_ratio",
    "rescue_ratio",
    "reproducibility_status",
    "reproducibility_json",
    "fingerprint_status",
    "fingerprint_json",
    "core_qc_completeness",
    "core_qc_status",
    "core_qc_warnings",
    "release_qc_status",
    "release_qc_warnings",
    "qc_status",
    "qc_warnings",
    "missing_files",
    "tagalign_file",
    "peak_file",
]


def split_values(text: str) -> list[str]:
    values: list[str] = []
    for token in text.replace(",", " ").split():
        token = token.strip()
        if token:
            values.append(token)
    return values


def as_int(value: object) -> Optional[int]:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def add_optional(total: Optional[int], value: Optional[int]) -> Optional[int]:
    if value is None:
        return total
    return (total or 0) + value


def format_int(value: Optional[int]) -> str:
    return "" if value is None else str(value)


def format_float(value: Optional[float], digits: int = 6) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"
