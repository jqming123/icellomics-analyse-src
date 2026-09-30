#!/usr/bin/env python3
"""Generate machine-readable ATAC QC sidecars and a Markdown report.

The report generator deliberately tolerates historical projects that do not
have a current ``sample_run_map.tsv``.  In that case a supplied run manifest
is used for run-level provenance and every summary row still receives an
experiment-level JSON sidecar and a Markdown entry.
"""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from atac_qc.alignment_qc import parse_flagstat_output, read_idxstats, read_picard_dup_details
from atac_qc.fastp_qc import read_fastp_metrics
from atac_qc.sample_map import group_replicates, read_replicate_records


QC_MODES = ("existing-results", "core", "release")


def read_tsv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_tsv_row(path: Path) -> dict[str, str]:
    rows = read_tsv_rows(path)
    return rows[0] if rows else {}


def read_picard_table(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    header: list[str] | None = None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("##"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
            else:
                return dict(zip(header, fields))
    return {}


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def peak_count(path: Path) -> int | None:
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def manifest_path(value: str, manifest: Path) -> str:
    """Resolve a manifest-relative path for a provenance payload only."""

    if not value:
        return ""
    path = Path(value).expanduser()
    return str((manifest.parent / path).resolve()) if not path.is_absolute() else str(path)


def legacy_run_payloads(manifest: Path, project: Path) -> list[tuple[str, dict[str, Any]]]:
    """Return basic read-only run JSON sidecars from an optional legacy manifest."""

    payloads: list[tuple[str, dict[str, Any]]] = []
    for index, row in enumerate(read_tsv_rows(manifest), start=1):
        biosample = row.get("biosample_id", "")
        if not biosample or biosample.startswith("#"):
            continue
        run_id = row.get("run_id", "") or f"manifest_run_{index}"
        path_fields = {
            name: manifest_path(row.get(name, ""), manifest)
            for name in (
                "fastp_json", "raw_bam", "raw_mito_qc", "dup_metrics",
                "pre_mito_idxstats", "final_idxstats", "final_bam",
            )
            if row.get(name, "")
        }
        mapq_evidence = row.get("mapq_filter_evidence", "")
        if mapq_evidence:
            path_fields["mapq_filter_evidence"] = (
                mapq_evidence[len("literal:"):].strip()
                if mapq_evidence.startswith("literal:")
                else manifest_path(mapq_evidence, manifest)
            )
        payloads.append((run_id, {
            "level": "run",
            "project_name": project.name,
            "source": "legacy_run_manifest",
            "biosample_id": biosample,
            "biological_replicate": row.get("biological_replicate", ""),
            "run_id": run_id,
            "layout": row.get("layout", "AUTO"),
            "mapq_filter_threshold": row.get("mapq_filter_threshold", ""),
            "mapq_filter_evidence": row.get("mapq_filter_evidence", ""),
            "declared_files": path_fields,
        }))
    return payloads


def write_current_run_json(
    project: Path,
    json_root: Path,
    artifact_root: Path | None = None,
) -> dict[str, object]:
    """Write current-pipeline run/replicate JSON and return grouped records."""

    records = read_replicate_records(project / "0_data" / "sample_run_map.tsv")
    grouped = group_replicates(records)
    for record in records:
        bowtie = project / "1_result" / "1_alignment" / record.run_id / "bowtie2"
        flagstat_path = bowtie / f"{record.run_id}.raw.bam.flagstat"
        flagstat = (
            parse_flagstat_output(flagstat_path.read_text(encoding="utf-8", errors="replace"))
            if flagstat_path.exists() else {}
        )
        mapping_filter = bowtie / f"{record.run_id}.mapping_filter.provenance.json"
        mapping_filter_payload: object = None
        if mapping_filter.exists():
            try:
                mapping_filter_payload = json.loads(mapping_filter.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                mapping_filter_payload = {"status": "INVALID_JSON", "path": str(mapping_filter)}
        write_json(json_root / "run" / f"{record.run_id}.qc.json", {
            "level": "run",
            "project_name": project.name,
            "biosample_id": record.biosample_id,
            "biological_replicate": record.biological_replicate,
            "run_id": record.run_id,
            "layout": record.layout,
            "fastp": read_fastp_metrics(
                project / "1_result" / "1_alignment" / record.run_id / "fastp" /
                f"{record.run_id}_fastp.json"
            ),
            "raw_alignment_flagstat": flagstat,
            "raw_mito": read_tsv_row(bowtie / f"{record.run_id}.raw.mito.qc.tsv"),
            "mapping_filter_provenance": mapping_filter_payload,
            "bowtie2_log": str(bowtie / f"{record.run_id}.bowtie2.log"),
        })

    for biosample, replicates in grouped.items():
        for replicate, rep_records in replicates.items():
            rep_dir = project / "1_result" / "1_alignment" / "replicates" / biosample / replicate
            prefix = f"{biosample}.{replicate}"
            generated_complexity = (
                artifact_root / "library_complexity" / biosample / f"{prefix}.lib_complexity.qc"
                if artifact_root is not None else None
            )
            # Generated QC sidecars are mode-specific.  The read-only fallback
            # keeps older reports interpretable without ever writing under
            # ``1_alignment``.
            complexity_path = (
                generated_complexity
                if generated_complexity is not None and generated_complexity.exists()
                else rep_dir / f"{prefix}.lib_complexity.qc"
            )
            complexity = read_tsv_row(complexity_path)
            duplication = read_picard_dup_details(rep_dir / f"{prefix}.dup.qc")
            final_total, _ = read_idxstats(rep_dir / f"{prefix}.nodup.bam.idxstats", ("MT",))
            write_json(json_root / "replicate" / biosample / f"{replicate}.qc.json", {
                "level": "biological_replicate",
                "project_name": project.name,
                "biosample_id": biosample,
                "biological_replicate": replicate,
                "run_ids": [record.run_id for record in rep_records],
                "layout": rep_records[0].layout,
                "duplication": duplication,
                "library_complexity": complexity,
                "library_complexity_file": str(complexity_path) if complexity else "",
                "gc_bias": read_picard_table(rep_dir / f"{prefix}.gc_bias.summary.txt"),
                "gc_bias_status": (
                    (rep_dir / f"{prefix}.gc_bias.status.txt").read_text(encoding="utf-8").strip()
                    if (rep_dir / f"{prefix}.gc_bias.status.txt").exists() else "INCOMPLETE"
                ),
                "final_nodup_reads": final_total,
                "peak_count": peak_count(
                    project / "1_result" / "3_peak_calling" / "replicates" / biosample /
                    f"{prefix}_peaks.narrowPeak"
                ),
            })
    return grouped


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True, type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--mode", choices=QC_MODES, default="release")
    parser.add_argument("--run-manifest", type=Path)
    parser.add_argument("--biosample-manifest", type=Path)
    args = parser.parse_args()

    if args.mode != "existing-results" and (args.run_manifest or args.biosample_manifest):
        parser.error(
            "run/biosample manifests are supported only in existing-results mode; "
            "core and release reports use canonical pipeline outputs"
        )

    project = args.project_dir
    qc_dir = project / "1_result" / "5_qc"
    summary_path = args.summary or qc_dir / "atac_qc_summary.tsv"
    # Preserve the historic default release-report filename for direct callers;
    # generated jobs pass an explicit mode-qualified output path so that a
    # legacy baseline cannot overwrite a later release report.
    default_report = (
        qc_dir / "atac_qc_report.md"
        if args.mode == "release"
        else qc_dir / f"atac_qc_report_{args.mode}.md"
    )
    report_path = args.output or default_report
    if not summary_path.exists():
        raise SystemExit(f"QC summary does not exist: {summary_path}")

    summary_rows = read_tsv_rows(summary_path)
    # Keep mode-specific JSON sidecars beside the matching mode-specific
    # summary/report.  In particular, a release QC run must not overwrite the
    # JSON provenance used to interpret an existing-results baseline.
    json_root = qc_dir / "json" / args.mode
    map_file = project / "0_data" / "sample_run_map.tsv"
    if args.run_manifest:
        for run_id, payload in legacy_run_payloads(args.run_manifest, project):
            write_json(json_root / "run" / f"{run_id}.qc.json", payload)
    elif args.mode != "existing-results" and map_file.exists():
        write_current_run_json(project, json_root, qc_dir / "artifacts" / args.mode)

    # An experiment-sidecar is always emitted, including for manifest-only
    # legacy datasets with no recoverable run or biological-replicate records.
    for row in summary_rows:
        biosample = row.get("biosample_id", "")
        repro_path = Path(row["reproducibility_json"]) if row.get("reproducibility_json") else None
        fingerprint_path = Path(row["fingerprint_json"]) if row.get("fingerprint_json") else None
        def optional_json(path: Path | None) -> object:
            if path is None or not path.exists():
                return None
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                return {"status": "INVALID_JSON", "path": str(path)}

        write_json(json_root / "experiment" / f"{biosample}.qc.json", {
            "level": "experiment",
            "project_name": project.name,
            "biosample_id": biosample,
            "qc_mode": args.mode,
            "qc_summary": row,
            "reproducibility": optional_json(repro_path),
            "fingerprint": optional_json(fingerprint_path),
            "input_manifests": {
                "run_manifest": str(args.run_manifest or ""),
                "biosample_manifest": str(args.biosample_manifest or ""),
            },
        })

    generated = datetime.now(timezone.utc).isoformat()
    lines = [
        f"# ATAC-seq QC report: {project.name}", "",
        f"- Generated: `{generated}`",
        f"- QC mode: `{args.mode}`",
        f"- Summary source: `{summary_path}`",
        "- Report format: Markdown (no HTML report generated)",
    ]
    if args.run_manifest:
        lines.append(f"- Run manifest: `{args.run_manifest}`")
    if args.biosample_manifest:
        lines.append(f"- Biosample manifest: `{args.biosample_manifest}`")

    lines.extend([
        "", "## Experiment summary", "",
        "| Biosample | Core completeness | Core QC | Release QC | Mapping | MAPQ evidence | MAPQ retention | Fragment size | Fragment input | FRiP | FRiP status | TSS | TSS shift provenance | Peak blacklist | Blacklist overlap | Reproducibility | Active status |",
        "| --- | --- | --- | --- | ---: | --- | ---: | ---: | --- | ---: | --- | ---: | --- | --- | ---: | --- | --- |",
    ])
    for row in summary_rows:
        fields = [
            row.get("biosample_id", ""),
            row.get("core_qc_completeness", ""),
            row.get("core_qc_status", ""),
            row.get("release_qc_status", ""),
            row.get("mapping_rate", ""),
            row.get("mapq_filter_status", ""),
            row.get("mapq_filter_retention", ""),
            row.get("fragment_median_size", ""),
            row.get("fragment_input_status", ""),
            row.get("fragment_frip_fraction", "") or row.get("frip_fraction", ""),
            row.get("frip_qc_status", ""),
            row.get("tss_enrichment_score", ""),
            row.get("tss_tagalign_shift_status", ""),
            row.get("peak_blacklist_status", ""),
            row.get("blacklist_fraction", ""),
            row.get("reproducibility_status", ""),
            row.get("qc_status", ""),
        ]
        lines.append("| " + " | ".join(fields) + " |")

    lines.extend(["", "## Interpretation", ""])
    for row in summary_rows:
        warnings = row.get("qc_warnings", "") or "none"
        core_warnings = row.get("core_qc_warnings", "") or "none"
        lines.append(
            f"- **{row.get('biosample_id', '')}** — core `{row.get('core_qc_status', '')}` "
            f"(completeness `{row.get('core_qc_completeness', '')}`), release "
            f"`{row.get('release_qc_status', '')}`; active `{row.get('qc_status', '')}`. "
            f"Core notes: {core_warnings}. Active-mode notes: {warnings}"
        )

    lines.extend([
        "", "## Methods and boundaries", "",
        "- Mapping/MAPQ is reported as confirmed only when a script, log, or mapping-filter provenance sidecar supplies evidence; the TSV also records MAPQ operator, SAM flags, raw/post-filter records and retention when a new sidecar supplies them. A final BAM alone is not used to reconstruct historical filtering.",
        "- Fragment FRiP is the reviewer-comparable value. A tagAlign-only result is separately labelled legacy record-level FRiP and is not evaluated against a fragment-FRiP threshold.",
        "- TSS enrichment uses Tn5-shifted cut sites. Historical manifests must declare `tagalign_tn5_shifted=true` for a method-comparable TSS result; otherwise the score is explicitly diagnostic.",
        "- Blacklist overlap is a post-hoc tagAlign measurement. It never proves that a historical peak set was blacklist-filtered; `Peak blacklist` comes only from the peak/reproducibility provenance sidecar.",
        "- `core_qc_status` evaluates recoverable reviewer-facing evidence. `release_qc_status` evaluates release-only artefacts such as blacklist application, IDR, GC bias, fingerprint, and track provenance.",
    ])
    if args.mode == "existing-results":
        lines.extend([
            "- This is a read-only historical baseline. `release_qc_status=NOT_REQUESTED` does not invalidate a measured core metric and does not imply that old peaks were reprocessed.",
        ])
    lines.extend([
        "", "## Machine-readable outputs", "",
        f"- Run JSON: `{json_root / 'run'}`",
        f"- Replicate JSON: `{json_root / 'replicate'}`",
        f"- Experiment JSON: `{json_root / 'experiment'}`",
    ])

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote Markdown QC report: {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
