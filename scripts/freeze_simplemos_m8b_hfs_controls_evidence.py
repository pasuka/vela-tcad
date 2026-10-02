#!/usr/bin/env python3
"""Freeze portable SimpleMOS M8-B HFS-control evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8b_hfs_controls/comparisons/comparison_report.json"
)
DEFAULT_DESTINATION = (
    REPO / "reference_tcad/simplemos_sentaurus2022/hfs_controls"
)
DEFAULT_EVIDENCE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8b_hfs_controls_evidence.json"
)
DEFAULT_CONVERGENCE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8b_eparallel_convergence/convergence_report.json"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def resolve(path: str) -> Path:
    value = Path(path)
    return value if value.is_absolute() else REPO / value


def freeze_report(source: Path, destination: Path) -> dict[str, Any]:
    report = read_json(source)
    comparisons = destination / "comparisons"
    comparisons.mkdir(parents=True, exist_ok=True)
    for case in report["cases"]:
        for key in ("comparison_csv", "response_csv"):
            current = resolve(case[key])
            target = comparisons / current.name
            shutil.copy2(current, target)
            case[key] = portable(target)
            case[f"{key}_sha256"] = sha256(target)
    frozen = comparisons / "comparison_report.json"
    write_json(frozen, report)
    return report


def write_summary(report: dict[str, Any], output: Path) -> None:
    fields = ["device", "drain_voltage_V", "variant", "dimension",
              "maximum_absolute_log10_ratio_above_floor",
              "median_absolute_log10_ratio_above_floor",
              "maximum_error_change_dex", "median_error_change_dex",
              "maximum_sentaurus_response_dex",
              "median_sentaurus_response_dex"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            writer.writerow({
                "device": case["device"],
                "drain_voltage_V": case["drain_voltage_V"],
                "variant": case["variant"], "dimension": case["dimension"],
                "maximum_absolute_log10_ratio_above_floor": case[
                    "maximum_absolute_log10_ratio_above_floor"],
                "median_absolute_log10_ratio_above_floor": case[
                    "median_absolute_log10_ratio_above_floor"],
                "maximum_error_change_dex": case["effect_vs_default"][
                    "maximum_log10_error_change"],
                "median_error_change_dex": case["effect_vs_default"][
                    "median_log10_error_change"],
                "maximum_sentaurus_response_dex": case[
                    "maximum_sentaurus_response_dex"],
                "median_sentaurus_response_dex": case[
                    "median_sentaurus_response_dex"],
            })


def freeze_convergence(source: Path, destination: Path) -> dict[str, Any]:
    report = read_json(source)
    portable_report = dict(report)
    portable_cases = []
    for case in report["cases"]:
        portable_case = dict(case)
        portable_case.pop("deck", None)
        portable_case.pop("console_log", None)
        portable_cases.append(portable_case)
    portable_report["cases"] = portable_cases
    target = destination / "convergence/convergence_report.json"
    write_json(target, portable_report)
    return portable_report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination", type=Path,
                        default=DEFAULT_DESTINATION)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--convergence", type=Path,
                        default=DEFAULT_CONVERGENCE)
    args = parser.parse_args()
    destination = args.destination.resolve()
    contract = read_json(
        REPO / "reference_tcad/simplemos_sentaurus2022"
        / "simplemos_m8b_hfs_controls_contract_v1.json")
    report = freeze_report(args.source.resolve(), destination)
    convergence = freeze_convergence(
        args.convergence.resolve(), destination)
    summary = destination / "hfs_controls_summary.csv"
    write_summary(report, summary)

    summaries = report["variant_summaries"]
    non_null = [item for item in summaries
                if item["variant"] != "explicit_gradqf"]
    mean_changes = {
        item["variant"]: sum(item["maximum_error_changes_dex"])
        / item["condition_count"] for item in non_null
    }
    best = min(mean_changes, key=mean_changes.get)
    figures = [
        REPO / "docs/validation/figures/simplemos_m8b/simplemos_m8b_sentaurus_response.png",
        REPO / "docs/validation/figures/simplemos_m8b/simplemos_m8b_error_change.png",
        REPO / "docs/validation/figures/simplemos_m8b/simplemos_m8b_residual.png",
        REPO / "docs/validation/figures/simplemos_m8b/simplemos_m8b_best_control_idvg.png",
    ]
    implementations = [
        REPO / "reference_tcad/simplemos_sentaurus2022/simplemos_m8b_hfs_controls_contract_v1.json",
        REPO / "reference_tcad/simplemos_sentaurus2022/simplemos_m8b_eparallel_convergence_contract_v1.json",
        REPO / "scripts/run_simplemos_m8b_hfs_controls.py",
        REPO / "scripts/run_simplemos_m8b_eparallel_convergence.py",
        REPO / "scripts/plot_simplemos_m8b_hfs_controls.py",
        REPO / "scripts/freeze_simplemos_m8b_hfs_controls_evidence.py",
        REPO / "scripts/build_simplemos_m8b_report_artifact.py",
        REPO / "tests/regression/test_simplemos_m8b_hfs_controls.py",
        REPO / "tests/regression/test_simplemos_m8b_eparallel_convergence.py",
    ]
    frozen_report = destination / "comparisons/comparison_report.json"
    evidence = {
        "schema": "vela.simplemos.sdevice.m8b_hfs_controls_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only Sentaurus HFS controls; no SProcess claims",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": {
            "devices": ["n17", "n21"],
            "drain_voltages_V": [0.05, 1.0],
            "controls": 5, "curves": 20,
            "points_per_curve": 51, "total_direct_bias_points": 1020,
            "interpolation": "forbidden",
            "vela_basis": "fixed M8-A full-physics Vela curves",
            "numerical_retry_cases": [item["case"] for item in contract[
                "numerical_retry_policy"]["cases"]],
            "convergence_probe_cases": len(convergence["cases"]),
            "convergence_probe_complete_curves": sum(
                item["success"] for item in convergence["cases"]),
        },
        "findings": {
            "null_control_pass": report["null_control_pass"],
            "best_mean_error_change_control": best,
            "best_mean_maximum_error_change_dex": mean_changes[best],
            "material_in_all_conditions": [item["variant"] for item in non_null
                                             if item["material_in_all_conditions"]],
            "improves_all_conditions": [item["variant"] for item in non_null
                                         if item["improves_all_conditions"]],
            "worsens_all_conditions": [item["variant"] for item in non_null
                                        if item["worsens_all_conditions"]],
            "variant_summaries": summaries,
            "selected_eparallel_line_search_damping": convergence[
                "selected_probe"]["line_search_damping"],
            "eparallel_damping_robustness_pairs": convergence[
                "robustness_pairs"],
        },
        "attribution": {
            "verified": [
                "Explicit GradQuasiFermi is treated as a numerical null control",
                "All controls reuse the same TDR, contacts, drain biases, and exact 51-point gate lattice",
                "Vela curves are fixed while only one Sentaurus HFS control changes",
                "LineSearchDamping=0.01 restores both Vd=1 V Eparallel cases with all 51 exact gate points",
            ],
            "unresolved": [
                "A control response identifies a contributing option family but not a proprietary internal formula",
                "The RefDens experiment tests one value (1e8 cm^-3) rather than a full parameter sweep",
            ],
        },
        "artifacts": {
            "comparison_report": {"path": portable(frozen_report),
                                  "sha256": sha256(frozen_report)},
            "convergence_report": {
                "path": portable(destination / "convergence/convergence_report.json"),
                "sha256": sha256(destination / "convergence/convergence_report.json")},
            "summary": {"path": portable(summary), "sha256": sha256(summary)},
            "figures": [{"path": portable(path), "sha256": sha256(path)}
                        for path in figures],
        },
        "implementation_sha256": {portable(path): sha256(path)
                                  for path in implementations},
    }
    write_json(args.evidence.resolve(), evidence)
    print(json.dumps({"status": "frozen", "evidence": portable(
        args.evidence.resolve()), "best_control": best,
        "null_control_pass": report["null_control_pass"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
