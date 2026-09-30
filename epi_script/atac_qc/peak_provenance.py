"""Select the JBrowse peak file and report truthful blacklist provenance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def select_peak_provenance(
    pooled_peak: Path,
    reproducible_peak: Path,
    reproducibility_json: Path,
) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    metadata_recorded = reproducibility_json.exists()
    if metadata_recorded:
        payload = json.loads(reproducibility_json.read_text(encoding="utf-8"))

    applied_value = payload.get("blacklist_applied")
    if isinstance(applied_value, bool):
        blacklist_applied: bool | None = applied_value
    elif payload.get("blacklist_bed") and payload.get("blacklist_sha256"):
        # Backward-compatible inference for phase-2/3 JSON written before the
        # explicit blacklist_applied field was introduced.
        blacklist_applied = True
    elif metadata_recorded:
        blacklist_applied = False
    else:
        blacklist_applied = None

    use_reproducible = (
        reproducible_peak.exists()
        and reproducible_peak.stat().st_size > 0
        and payload.get("reproducibility_status") == "PASS"
    )
    active_peak = reproducible_peak if use_reproducible else pooled_peak
    prefix = "true_replicate_idr_0.05" if use_reproducible else "pooled_macs3"
    if blacklist_applied is True:
        peak_mode = f"blacklist_filtered_{prefix}"
    elif blacklist_applied is False:
        peak_mode = f"{prefix}_no_blacklist"
    else:
        peak_mode = f"{prefix}_blacklist_unknown"

    status = payload.get("blacklist_status")
    if not status:
        status = "APPLIED" if blacklist_applied is True else (
            "NOT_RECORDED" if blacklist_applied is None else "NOT_CONFIGURED"
        )
    return {
        "active_peak_file": str(active_peak),
        "peak_mode": peak_mode,
        "blacklist_applied": blacklist_applied,
        "blacklist_status": status,
        "blacklist_bed": str(payload.get("blacklist_bed") or ""),
        "blacklist_sha256": str(payload.get("blacklist_sha256") or ""),
        "blacklist_source": str(payload.get("blacklist_source") or ""),
        "peak_pipeline_mode": str(payload.get("peak_pipeline_mode") or ""),
        "blacklist_validation_status": str(payload.get("blacklist_validation_status") or ""),
        "blacklist_validation_sha256": str(payload.get("blacklist_validation_sha256") or ""),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pooled-peak", required=True, type=Path)
    parser.add_argument("--reproducible-peak", required=True, type=Path)
    parser.add_argument("--reproducibility-json", required=True, type=Path)
    parser.add_argument("--emit-lines", action="store_true")
    args = parser.parse_args()
    result = select_peak_provenance(
        args.pooled_peak, args.reproducible_peak, args.reproducibility_json
    )
    if args.emit_lines:
        for key in (
            "active_peak_file", "peak_mode", "blacklist_applied",
            "blacklist_status", "blacklist_bed", "blacklist_sha256",
            "blacklist_source", "peak_pipeline_mode",
            "blacklist_validation_status", "blacklist_validation_sha256",
        ):
            value = result[key]
            if value is None:
                value = "unknown"
            elif isinstance(value, bool):
                value = str(value).lower()
            print(value)
    else:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
