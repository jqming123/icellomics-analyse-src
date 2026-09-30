"""Profile-driven ATAC-seq QC status evaluation.

The pipeline distinguishes two deliberately different questions:

* ``core`` QC asks whether the evidence available for an experiment supports
  interpretation of the principal ATAC metrics (mapping, fragment sizes, TSS
  enrichment and FRiP).
* ``release`` QC asks whether a complete, current pipeline release has all
  publication/release artefacts (for example IDR, GC-bias and fingerprint
  outputs).

Keeping these questions separate is important for historical datasets.  A
legacy dataset can have a valid TSS profile or fragment-size distribution even
when it predates the current IDR/fingerprint workflow.  Such a dataset must be
reported as partially evidenced rather than having its available QC hidden by a
single ``INCOMPLETE`` release status.
"""

from __future__ import annotations

import json
import operator
from functools import lru_cache
from pathlib import Path


DEFAULT_THRESHOLDS = Path(__file__).resolve().parents[1] / "atac_qc_thresholds.json"

QC_MODES = frozenset({"existing-results", "core", "release"})

# These are the reviewer-facing, experiment-level measurements that a core
# ATAC-QC assessment can actually evaluate.  ``fragment_frip_fraction`` has a
# special alternative: a legacy tagAlign-record FRiP is useful as a diagnostic,
# but never receives the fragment-FRiP threshold.
CORE_REQUIRED_METRICS = (
    "mapping_rate",
    "fragment_count",
    "peak_count",
    "tss_enrichment_score",
)


def _as_float(value: object) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


@lru_cache(maxsize=8)
def load_profiles(path: str = str(DEFAULT_THRESHOLDS)) -> dict[str, dict[str, object]]:
    with Path(path).open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    profiles = data.get("profiles")
    if not isinstance(profiles, dict) or not profiles:
        raise ValueError(f"No QC profiles found in {path}")
    return profiles


def _evaluate_rule(metric: str, value: float, rule: dict[str, object]) -> tuple[str | None, str | None]:
    checks = (
        ("fail_below", operator.lt, "FAIL", "lt"),
        ("fail_at_or_below", operator.le, "FAIL", "le"),
        ("fail_above", operator.gt, "FAIL", "gt"),
        ("fail_at_or_above", operator.ge, "FAIL", "ge"),
        ("warn_below", operator.lt, "WARN", "lt"),
        ("warn_at_or_below", operator.le, "WARN", "le"),
        ("warn_above", operator.gt, "WARN", "gt"),
        ("warn_at_or_above", operator.ge, "WARN", "ge"),
    )
    for key, predicate, severity, operator_label in checks:
        boundary = _as_float(rule.get(key))
        if boundary is not None and predicate(value, boundary):
            return severity, f"{metric}_{operator_label}_{boundary:g}"
    return None, None


def evaluate_qc_status(
    metrics: dict[str, object],
    existing_warnings: list[str],
    profile_name: str = "human_grch38_ensembl_phase1",
    thresholds_path: Path | None = None,
) -> dict[str, str]:
    profiles = load_profiles(str(thresholds_path or DEFAULT_THRESHOLDS))
    if profile_name not in profiles:
        raise ValueError(f"Unknown QC threshold profile: {profile_name}")
    profile = profiles[profile_name]
    required = profile.get("required_metrics", [])
    incomplete = [
        f"missing_required_metric:{name}"
        for name in required
        if metrics.get(str(name)) is None or str(metrics.get(str(name))).strip() == ""
    ]
    if str(metrics.get("missing_files", "")).strip():
        incomplete.append("missing_required_files")
    for status_metric in ("gc_bias_status", "fingerprint_status", "reproducibility_status"):
        value = str(metrics.get(status_metric, "")).upper()
        if value in {"NOT_RUN", "INCOMPLETE", "NOT_APPLICABLE", "INVALID_JSON"}:
            incomplete.append(f"{status_metric}:{value.lower()}")

    # Assembly configurations with a verified exclusion resource must not
    # publish a release unless the reproducibility sidecar proves that the
    # selected peak call actually used it.  The post-hoc tagAlign overlap is a
    # separate diagnostic and cannot satisfy this gate.
    if str(metrics.get("blacklist_config_status", "")).upper() == "CONFIGURED":
        config_bed = str(metrics.get("blacklist_config_bed", "")).strip()
        config_sha256 = str(metrics.get("blacklist_config_sha256", "")).strip()
        config_fai = str(metrics.get("reference_config_fai", "")).strip()
        config_fai_sha256 = str(metrics.get("reference_config_fai_sha256", "")).strip()
        peak_sha256 = str(metrics.get("peak_blacklist_sha256", "")).strip()
        validation_sha256 = str(metrics.get("peak_blacklist_validation_sha256", "")).strip()
        if not config_bed:
            incomplete.append("blacklist_config_path_missing")
        if not config_sha256:
            incomplete.append("blacklist_config_sha256_missing")
        if not config_fai:
            incomplete.append("reference_config_fai_missing")
        if not config_fai_sha256:
            incomplete.append("reference_config_fai_sha256_missing")
        if str(metrics.get("peak_pipeline_mode", "")).lower() != "release":
            incomplete.append("peak_pipeline_not_release")
        if str(metrics.get("peak_blacklist_applied", "")).lower() != "true":
            incomplete.append("blacklist_filter_not_verified")
        if str(metrics.get("peak_blacklist_status", "")).upper() != "APPLIED":
            incomplete.append("blacklist_filter_status_not_applied")
        if not peak_sha256:
            incomplete.append("blacklist_filter_sha256_missing")
        elif config_sha256 and peak_sha256 != config_sha256:
            incomplete.append("blacklist_filter_sha256_config_mismatch")
        if str(metrics.get("peak_blacklist_validation_status", "")).upper() != "VERIFIED":
            incomplete.append("blacklist_validation_not_verified")
        if not str(metrics.get("peak_blacklist_validation_record", "")).strip():
            incomplete.append("blacklist_validation_record_missing")
        if config_bed and str(metrics.get("peak_blacklist_validation_bed", "")).strip() != config_bed:
            incomplete.append("blacklist_validation_path_config_mismatch")
        if not validation_sha256:
            incomplete.append("blacklist_validation_sha256_missing")
        elif config_sha256 and validation_sha256 != config_sha256:
            incomplete.append("blacklist_validation_sha256_config_mismatch")
        elif peak_sha256 and validation_sha256 != peak_sha256:
            incomplete.append("blacklist_validation_peak_sha256_mismatch")
        if not str(metrics.get("peak_blacklist_validation_reference_fai", "")).strip():
            incomplete.append("blacklist_validation_reference_fai_missing")
        elif config_fai and (
            str(metrics.get("peak_blacklist_validation_reference_fai", "")).strip()
            != config_fai
        ):
            incomplete.append("blacklist_validation_reference_fai_config_mismatch")
        if not str(metrics.get("peak_blacklist_validation_reference_fai_sha256", "")).strip():
            incomplete.append("blacklist_validation_reference_fai_sha256_missing")
        elif config_fai_sha256 and (
            str(metrics.get("peak_blacklist_validation_reference_fai_sha256", "")).strip()
            != config_fai_sha256
        ):
            incomplete.append("blacklist_validation_reference_fai_sha256_config_mismatch")

    failures: list[str] = []
    rule_warnings: list[str] = []
    for status_metric in ("gc_bias_status", "fingerprint_status"):
        value = str(metrics.get(status_metric, "")).upper()
        if value == "FAIL":
            failures.append(f"{status_metric}:fail")
        elif value == "WARN":
            rule_warnings.append(f"{status_metric}:warn")
    rules = profile.get("rules", {})
    if isinstance(rules, dict):
        for metric, raw_rule in rules.items():
            value = _as_float(metrics.get(metric))
            if value is None or not isinstance(raw_rule, dict):
                continue
            severity, message = _evaluate_rule(metric, value, raw_rule)
            if severity == "FAIL" and message:
                failures.append(message)
            elif severity == "WARN" and message:
                rule_warnings.append(message)

    if incomplete:
        status = "INCOMPLETE"
    elif failures:
        status = "FAIL"
    elif rule_warnings:
        status = "WARN"
    else:
        status = "PASS"
    messages = incomplete + failures + rule_warnings + list(existing_warnings)
    return {"qc_status": status, "qc_warnings": ";".join(dict.fromkeys(messages))}


def _append_message(messages: list[str], message: str) -> None:
    if message and message not in messages:
        messages.append(message)


def _profile_rules(profile_name: str, thresholds_path: Path | None) -> dict[str, object]:
    profiles = load_profiles(str(thresholds_path or DEFAULT_THRESHOLDS))
    if profile_name not in profiles:
        raise ValueError(f"Unknown QC threshold profile: {profile_name}")
    profile = profiles[profile_name]
    rules = profile.get("rules", {})
    return rules if isinstance(rules, dict) else {}


def evaluate_core_qc_status(
    metrics: dict[str, object],
    existing_warnings: list[str],
    profile_name: str = "human_grch38_ensembl_phase1",
    thresholds_path: Path | None = None,
) -> dict[str, str]:
    """Evaluate only evidence that belongs to core ATAC QC.

    Missing core inputs result in ``PARTIAL`` rather than ``INCOMPLETE``; no
    observed core measurement results in ``UNAVAILABLE``.  A numeric failure
    still takes precedence over incompleteness so that, for example, a poor
    observed mapping rate is not concealed by a missing IDR file.  The separate
    ``core_qc_completeness`` field uses the stable values ``COMPLETE``,
    ``PARTIAL`` and ``UNAVAILABLE`` for dashboards that need availability
    rather than an evaluative status.  Release-only artefacts and
    ``missing_files`` are intentionally not considered here.
    """

    missing: list[str] = []
    available_core_measurements = 0
    for name in CORE_REQUIRED_METRICS:
        if metrics.get(name) is None or str(metrics.get(name)).strip() == "":
            _append_message(missing, f"missing_core_metric:{name}")
        else:
            available_core_measurements += 1

    fragment_frip = _as_float(metrics.get("fragment_frip_fraction"))
    legacy_frip = _as_float(metrics.get("frip_fraction"))
    if fragment_frip is None and legacy_frip is None:
        _append_message(missing, "missing_core_metric:frip")
    else:
        available_core_measurements += 1

    fragment_input_status = str(metrics.get("fragment_input_status", "")).upper()
    if fragment_input_status == "PARTIAL":
        _append_message(missing, "partial_core_input:fragment_bams")
    elif fragment_input_status == "UNAVAILABLE":
        _append_message(missing, "missing_core_input:fragment_bams")

    warnings = list(existing_warnings)
    if fragment_frip is None and legacy_frip is not None:
        _append_message(warnings, "legacy_record_frip_diagnostic_only")

    failures: list[str] = []
    rule_warnings: list[str] = []
    rules = _profile_rules(profile_name, thresholds_path)

    # The profile can define more metrics than core QC is responsible for.
    # Evaluate only the rules that have unambiguous core semantics.  Fragment
    # FRiP rules are skipped for legacy tagAlign-record FRiP on purpose.
    core_rule_metrics = {"mapping_rate", "peak_count", "tss_enrichment_score"}
    if fragment_frip is not None:
        core_rule_metrics.add("fragment_frip_fraction")
    for metric in core_rule_metrics:
        raw_rule = rules.get(metric)
        value = _as_float(metrics.get(metric))
        if value is None or not isinstance(raw_rule, dict):
            continue
        severity, message = _evaluate_rule(metric, value, raw_rule)
        if severity == "FAIL" and message:
            _append_message(failures, message)
        elif severity == "WARN" and message:
            _append_message(rule_warnings, message)

    # A zero-count fragment-size result means a usable PE BAM was inspected
    # but no valid fragments were available.  It is more informative than a
    # generic missing-value warning and should not be allowed to pass.
    fragment_count = _as_float(metrics.get("fragment_count"))
    if fragment_count is not None and fragment_count <= 0:
        _append_message(failures, "fragment_count_at_or_below_0")

    # These provenance caveats do not invalidate measurements but should make
    # a fully measured core result visible as WARN rather than PASS.
    provenance_warning_prefixes = (
        "legacy_tagalign_record_frip_not_fragment_frip",
        "legacy_record_frip_diagnostic_only",
        "legacy_sample_map_biological_replicates_unknown",
        "mapq_filter_not_available_legacy",
        "mapq_filter_not_declared",
        "mapq_filter_threshold_declared_without_evidence",
        "mapq_filter_evidence_declared_without_threshold",
        "mapq_filter_evidence_path_missing",
        "tss_tagalign_shift_unknown_diagnostic_only",
        "tss_tagalign_not_tn5_shifted_diagnostic_only",
    )
    provenance_warnings = [
        message
        for message in warnings
        if message.startswith(provenance_warning_prefixes)
    ]

    if available_core_measurements == 0:
        completeness = "UNAVAILABLE"
    elif missing:
        completeness = "PARTIAL"
    else:
        completeness = "COMPLETE"

    if failures:
        status = "FAIL"
    elif completeness == "UNAVAILABLE":
        status = "UNAVAILABLE"
    elif completeness == "PARTIAL":
        status = "PARTIAL"
    elif rule_warnings or provenance_warnings:
        status = "WARN"
    else:
        status = "PASS"

    messages = missing + failures + rule_warnings + warnings
    return {
        "core_qc_completeness": completeness,
        "core_qc_status": status,
        "core_qc_warnings": ";".join(dict.fromkeys(messages)),
    }


def _apply_reproducibility_status(
    result: dict[str, str], metrics: dict[str, object]
) -> dict[str, str]:
    """Apply release-only IDR, MAPQ and fragment-size evidence gates."""

    status = result["qc_status"]
    messages = [value for value in result["qc_warnings"].split(";") if value]
    release_incomplete: list[str] = []

    # A release may not silently pass when the reviewer-facing fragment-size
    # measurement is absent.  It is meaningful only for PE libraries; a
    # single-end library is explicitly marked incomplete for this particular
    # release profile rather than being treated as an unqualified PASS.
    layout = str(metrics.get("layout", "")).upper()
    if layout == "PE":
        if str(metrics.get("fragment_input_status", "")).upper() != "COMPLETE":
            release_incomplete.append("fragment_size_input_not_complete")
        for field in ("fragment_count", "fragment_median_size"):
            if _as_float(metrics.get(field)) is None:
                release_incomplete.append(f"missing_release_fragment_size_metric:{field}")
    elif layout == "SE":
        release_incomplete.append("fragment_size_not_applicable_single_end")
    else:
        release_incomplete.append("fragment_size_layout_not_resolved")

    # New release runs must retain the structured mapping-filter provenance
    # written by ATAC_align.sh.  Historical data remains assessable in the
    # core/existing-results modes, where this function is not applied.
    if str(metrics.get("mapq_filter_status", "")).upper() != "EVIDENCE_FILE":
        release_incomplete.append("mapq_filter_provenance_not_verified")
    for field in (
        "mapq_filter_threshold",
        "mapq_filter_operator",
        "mapq_excluded_sam_flags",
        "mapq_raw_alignment_records",
        "mapq_post_filter_alignment_records",
        "mapq_post_final_filter_alignment_records",
    ):
        if str(metrics.get(field, "")).strip() == "":
            release_incomplete.append(f"missing_release_mapq_field:{field}")
    if layout == "PE" and not str(metrics.get("mapq_required_sam_flags", "")).strip():
        release_incomplete.append("missing_release_mapq_field:mapq_required_sam_flags")
    if str(metrics.get("peak_pipeline_mode", "")).lower() != "release":
        release_incomplete.append("peak_pipeline_not_release")

    if release_incomplete:
        status = "INCOMPLETE"
        for message in release_incomplete:
            _append_message(messages, message)
    reproducibility = str(metrics.get("reproducibility_status", "")).upper()
    if reproducibility in {
        "NOT_RUN_POOLED_ONLY", "PAIRWISE_IDR_ONLY", "NOT_ASSESSED", "INVALID_JSON"
    }:
        # A pooled-only call and pairwise-only diagnostic are useful outputs,
        # but neither establishes the release-level experiment reproducibility
        # that the full profile represents.
        status = "INCOMPLETE"
        _append_message(messages, f"reproducibility_not_assessed:{reproducibility.lower()}")
    elif status != "INCOMPLETE" and reproducibility == "FAIL":
        status = "FAIL"
        _append_message(messages, "idr_reproducibility_failed")
    elif status == "PASS" and reproducibility == "UNREPLICATED":
        status = "WARN"
        _append_message(messages, "single_biological_replicate")
    return {"qc_status": status, "qc_warnings": ";".join(messages)}


def evaluate_qc_statuses(
    metrics: dict[str, object],
    existing_warnings: list[str],
    profile_name: str = "human_grch38_ensembl_phase1",
    mode: str = "release",
    thresholds_path: Path | None = None,
) -> dict[str, str]:
    """Return independent core and release statuses plus active compatibility fields.

    ``qc_status`` and ``qc_warnings`` are retained for callers that consumed
    the former one-status TSV.  They select the core result in
    ``existing-results``/``core`` modes and the release result in ``release``
    mode.  ``existing-results`` never asserts that historical peaks have been
    reprocessed by the release workflow.
    """

    if mode not in QC_MODES:
        raise ValueError(f"Unknown QC mode: {mode}; expected one of {sorted(QC_MODES)}")

    core = evaluate_core_qc_status(
        metrics,
        existing_warnings,
        profile_name=profile_name,
        thresholds_path=thresholds_path,
    )
    if mode == "release":
        release = _apply_reproducibility_status(
            evaluate_qc_status(
                metrics,
                existing_warnings,
                profile_name=profile_name,
                thresholds_path=thresholds_path,
            ),
            metrics,
        )
    else:
        release = {
            "qc_status": "NOT_REQUESTED",
            "qc_warnings": f"release_qc_not_requested:{mode}",
        }

    active = release if mode == "release" else {
        "qc_status": core["core_qc_status"],
        "qc_warnings": core["core_qc_warnings"],
    }
    return {
        "core_qc_completeness": core["core_qc_completeness"],
        "core_qc_status": core["core_qc_status"],
        "core_qc_warnings": core["core_qc_warnings"],
        "release_qc_status": release["qc_status"],
        "release_qc_warnings": release["qc_warnings"],
        "qc_status": active["qc_status"],
        "qc_warnings": active["qc_warnings"],
    }
