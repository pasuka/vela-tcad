#!/usr/bin/env python3
"""Freeze portable SimpleMOS M8-A comparison evidence and attribution metadata."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import shutil
from typing import Any


REPO = Path(__file__).resolve().parents[1]
SOURCE_ROOT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
)
DEFAULT_DESTINATION = (
    REPO / "reference_tcad/simplemos_sentaurus2022/model_ablation"
)
EVIDENCE_PATH = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8a_model_ablation_evidence.json"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def sha256(path: Path) -> str:
    import hashlib
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def freeze_report(source: Path, destination: Path) -> dict[str, Any]:
    report = read_json(source)
    comparisons = destination / "comparisons"
    comparisons.mkdir(parents=True, exist_ok=True)
    for case in report["cases"]:
        current = Path(case["comparison_csv"])
        if not current.is_absolute():
            current = REPO / current
        target = comparisons / current.name
        shutil.copy2(current, target)
        case["comparison_csv"] = portable(target)
        case["comparison_csv_sha256"] = sha256(target)
    frozen = comparisons / "comparison_report.json"
    write_json(frozen, report)
    return {"path": portable(frozen), "sha256": sha256(frozen),
            "report": report}


def write_first_round_summary(report: dict[str, Any], output: Path) -> None:
    fields = ["variant", "status", "maximum_absolute_log10_ratio",
              "maximum_relative_error", "maximum_change_vs_full_dex",
              "deep_off_vg_0p05_change_vs_full_dex",
              "weak_inversion_change_vs_full_dex"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            effect = case["effect_vs_full"]
            writer.writerow({
                "variant": case["variant"], "status": case["status"],
                "maximum_absolute_log10_ratio": case[
                    "maximum_absolute_log10_ratio_above_floor"],
                "maximum_relative_error": case[
                    "maximum_relative_error_above_floor"],
                "maximum_change_vs_full_dex": effect[
                    "maximum_log10_ratio_change"],
                "deep_off_vg_0p05_change_vs_full_dex": effect[
                    "deep_off_vg_0p05_log10_ratio_change"],
                "weak_inversion_change_vs_full_dex": effect[
                    "weak_inversion_max_log10_ratio_change"],
            })


def write_confirmation_summary(report: dict[str, Any], output: Path) -> None:
    fields = ["device", "drain_voltage_V", "variant", "status",
              "maximum_absolute_log10_ratio", "maximum_relative_error",
              "maximum_change_vs_full_dex",
              "deep_off_vg_0p05_change_vs_full_dex",
              "weak_inversion_change_vs_full_dex"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            effect = case["effect_vs_full"]
            writer.writerow({
                "device": case["device"],
                "drain_voltage_V": case["drain_voltage_V"],
                "variant": case["variant"], "status": case["status"],
                "maximum_absolute_log10_ratio": case[
                    "maximum_absolute_log10_ratio_above_floor"],
                "maximum_relative_error": case[
                    "maximum_relative_error_above_floor"],
                "maximum_change_vs_full_dex": effect[
                    "maximum_log10_ratio_change"],
                "deep_off_vg_0p05_change_vs_full_dex": effect[
                    "deep_off_vg_0p05_log10_ratio_change"],
                "weak_inversion_change_vs_full_dex": effect[
                    "weak_inversion_max_log10_ratio_change"],
            })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--evidence", type=Path, default=EVIDENCE_PATH)
    args = parser.parse_args()
    source_root = args.source_root.resolve()
    destination = args.destination.resolve()
    first = freeze_report(
        source_root / "m8a_model_ablation/comparisons/comparison_report.json",
        destination / "first_round")
    confirmation = freeze_report(
        source_root / "m8a_confirmation/comparisons/comparison_report.json",
        destination / "confirmation")
    first_summary = destination / "first_round/ablation_summary.csv"
    confirmation_summary = destination / "confirmation/confirmation_summary.csv"
    write_first_round_summary(first["report"], first_summary)
    write_confirmation_summary(confirmation["report"], confirmation_summary)

    first_by_variant = {case["variant"]: case for case in first["report"]["cases"]}
    confirmation_by_variant = {}
    for variant in ("no_hfs", "phumob_only"):
        confirmation_by_variant[variant] = [
            case["effect_vs_full"]["maximum_log10_ratio_change"]
            for case in confirmation["report"]["cases"]
            if case["variant"] == variant
        ]
    hfs_consistent = all(value < 0.0 for value in
                         confirmation_by_variant["no_hfs"])
    figures = [
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_paired_idvg.png",
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_signed_residual.png",
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_gap_change.png",
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_model_response.png",
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_confirmation_full_idvg.png",
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_confirmation_residual.png",
        REPO / "docs/validation/figures/simplemos_m8a/simplemos_m8a_confirmation_improvement.png",
    ]
    implementations = [
        REPO / "reference_tcad/simplemos_sentaurus2022/simplemos_m8a_model_ablation_contract_v1.json",
        REPO / "reference_tcad/simplemos_sentaurus2022/simplemos_m8a_confirmation_contract_v1.json",
        REPO / "scripts/run_simplemos_m8a_model_ablation.py",
        REPO / "scripts/run_simplemos_m8a_confirmation.py",
        REPO / "scripts/plot_simplemos_m8a_model_ablation.py",
        REPO / "scripts/plot_simplemos_m8a_confirmation.py",
        REPO / "scripts/freeze_simplemos_m8a_evidence.py",
        REPO / "scripts/build_simplemos_m8a_report_artifact.py",
        REPO / "tests/regression/test_simplemos_m8a_model_ablation.py",
    ]
    first_full = first_by_variant["full"]
    first_no_hfs = first_by_variant["no_hfs"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m8a_model_attribution_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only paired model ablation; no SProcess claims",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": {
            "first_round": {"device": "n23", "drain_voltage_V": 0.05,
                            "variants": 8, "curves": 8, "points": 408},
            "confirmation": {"devices": ["n17", "n21"],
                             "drain_voltages_V": [0.05, 1.0],
                             "variants": 3, "curves": 12, "points": 612},
            "total_direct_bias_points": 1020,
            "interpolation": "forbidden",
        },
        "findings": {
            "first_round_full_maximum_error_dex": first_full[
                "maximum_absolute_log10_ratio_above_floor"],
            "first_round_no_hfs_maximum_error_dex": first_no_hfs[
                "maximum_absolute_log10_ratio_above_floor"],
            "first_round_no_hfs_change_dex": first_no_hfs[
                "effect_vs_full"]["maximum_log10_ratio_change"],
            "no_old_slotboom_maximum_error_dex": first_by_variant["no_bgn"][
                "maximum_absolute_log10_ratio_above_floor"],
            "confirmation_no_hfs_changes_dex": confirmation_by_variant["no_hfs"],
            "confirmation_phumob_only_changes_dex": confirmation_by_variant[
                "phumob_only"],
            "no_hfs_improves_all_confirmation_conditions": hfs_consistent,
        },
        "attribution": {
            "verified": [
                "SRH removal and plain SRH both worsen the n23 deep-off mismatch",
                "OldSlotboom removal causes a large cross-solver divergence",
                "No-HFS is the best first-round mismatch-reducing ablation",
            ],
            "likely": ([
                "Differences in HighFieldSaturation driving-force, smoothing, parameters, or limiting behavior contribute to the remaining mismatch",
            ] if hfs_consistent else [
                "HighFieldSaturation contributes to the n23 residual but is not stable across all confirmation conditions",
            ]),
            "unresolved": [
                "The experiment cannot isolate proprietary Sentaurus internal smoothing from parameter or driving-force differences",
                "Ablation evidence is diagnostic and does not prove one internal formula is the sole cause",
            ],
        },
        "artifacts": {
            "first_round_report": {"path": first["path"],
                                   "sha256": first["sha256"]},
            "first_round_summary": {"path": portable(first_summary),
                                    "sha256": sha256(first_summary)},
            "confirmation_report": {"path": confirmation["path"],
                                    "sha256": confirmation["sha256"]},
            "confirmation_summary": {"path": portable(confirmation_summary),
                                     "sha256": sha256(confirmation_summary)},
            "figures": [{"path": portable(path), "sha256": sha256(path)}
                        for path in figures],
        },
        "implementation_sha256": {portable(path): sha256(path)
                                  for path in implementations},
    }
    write_json(args.evidence.resolve(), evidence)
    print(json.dumps({"status": "frozen", "hfs_consistent": hfs_consistent,
                      "evidence": portable(args.evidence.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
