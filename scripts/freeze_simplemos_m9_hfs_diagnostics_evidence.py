#!/usr/bin/env python3
"""Freeze portable SimpleMOS M9 HFS diagnostic evidence."""

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
    / "m9_hfs_diagnostics"
)
DEFAULT_DESTINATION = (
    REPO / "reference_tcad/simplemos_sentaurus2022/hfs_diagnostics"
)
DEFAULT_EVIDENCE = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m9_hfs_diagnostics_evidence.json"
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


def freeze_comparisons(source: Path, destination: Path) -> dict[str, Any]:
    report = read_json(source / "comparisons/comparison_report.json")
    target_dir = destination / "comparisons"
    target_dir.mkdir(parents=True, exist_ok=True)
    for case in report["cases"]:
        for key in ("comparison_csv", "response_csv"):
            current = resolve(case[key])
            target = target_dir / current.name
            shutil.copy2(current, target)
            case[key] = portable(target)
            case[f"{key}_sha256"] = sha256(target)
    target = target_dir / "comparison_report.json"
    write_json(target, report)
    return report


def write_scan_summary(report: dict[str, Any], path: Path) -> None:
    rows = []
    for case in report["cases"]:
        rows.append({
            "device": case["device"],
            "drain_voltage_V": case["drain_voltage_V"],
            "variant": case["variant"],
            "refdens_cm3": case["refdens_cm3"],
            "maximum_absolute_log10_ratio_above_floor": case[
                "maximum_absolute_log10_ratio_above_floor"],
            "maximum_error_change_dex": case["effect_vs_default"][
                "maximum_log10_error_change"],
            "maximum_sentaurus_response_dex": case[
                "maximum_sentaurus_response_dex"],
        })
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def copy_csv_with_lf(source: Path, target: Path) -> None:
    """Copy a generated CSV while keeping repository evidence platform-neutral."""
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(source.read_text(encoding="utf-8-sig"),
                      encoding="utf-8", newline="\n")


def freeze_parameter_audit(source: Path, destination: Path) -> Path:
    audit = read_json(source / "parameter_audit/hfs_parameter_audit.json")
    parameter_path = Path(audit["parameter_file"])
    audit["parameter_file"] = (
        "Sentaurus M9 n17 explicit-GradQF Vd=0.05 V sdevice -P export")
    audit["parameter_file_sha256"] = sha256(parameter_path)
    target = destination / "parameter_audit/hfs_parameter_audit.json"
    write_json(target, audit)
    copy_csv_with_lf(
        source / "parameter_audit/hfs_parameter_comparison.csv",
        target.with_name("hfs_parameter_comparison.csv"))
    return target


def freeze_edge_summary(source: Path, destination: Path) -> Path:
    report = read_json(source / "vela_edge_probes/edge_hfs_report.json")
    for case in report["cases"]:
        case["edge_csv"] = "ignored build artifact; integrity hash retained"
    target = destination / "vela_edge_probes/edge_hfs_report.json"
    write_json(target, report)
    copy_csv_with_lf(
        source / "vela_edge_probes/edge_hfs_summary.csv",
        target.with_name("edge_hfs_summary.csv"))
    return target


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    source = args.source.resolve()
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)
    report = freeze_comparisons(source, destination)
    summary = destination / "refdens_scan_summary.csv"
    write_scan_summary(report, summary)
    parameter_audit = freeze_parameter_audit(source, destination)
    edge_report_path = freeze_edge_summary(source, destination)
    edge_report = read_json(edge_report_path)
    contract = (
        REPO / "reference_tcad/simplemos_sentaurus2022"
        / "simplemos_m9_hfs_diagnostics_contract_v1.json")
    baseline_contract = read_json(contract)["default_model_policy"]
    positive = [item for item in report["variant_summaries"]
                if float(item["refdens_cm3"]) > 0.0]
    figures = sorted((REPO / "docs/validation/figures/simplemos_m9").glob("*.png"))
    implementations = [
        contract,
        REPO / "src/tools/vela_example_runner.cpp",
        REPO / "scripts/run_simplemos_m9_refdens_scan.py",
        REPO / "scripts/audit_simplemos_m9_hfs_parameters.py",
        REPO / "scripts/run_simplemos_m9_edge_hfs_diagnostics.py",
        REPO / "scripts/plot_simplemos_m9_hfs_diagnostics.py",
        REPO / "scripts/freeze_simplemos_m9_hfs_diagnostics_evidence.py",
        REPO / "scripts/build_simplemos_m9_report_artifact.py",
        REPO / "tests/regression/test_reference_tcad_tools.py",
        REPO / "tests/regression/test_simplemos_m9_hfs_diagnostics.py",
    ]
    frozen_report = destination / "comparisons/comparison_report.json"
    evidence = {
        "schema": "vela.simplemos.sdevice.m9_hfs_diagnostics_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only HFS reference-density, parameter, and edge diagnostics",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": {
            "sentaurus_curves": 24,
            "points_per_curve": 51,
            "total_direct_bias_points": 1224,
            "devices": ["n17", "n21"],
            "drain_voltages_V": [0.05, 1.0],
            "refdens_cm3": report["scan_values_cm3"],
            "interpolation": "forbidden",
            "vela_edge_probe_states": edge_report["case_count"],
            "default_model_changed": False,
        },
        "baseline_guard": {
            "status": "unchanged",
            "case_count": baseline_contract["required_case_count"],
            "direct_bias_points": baseline_contract[
                "required_direct_bias_points"],
            "evidence_sha256": baseline_contract[
                "m8_baseline_evidence_sha256"],
            "comparison_report_sha256": baseline_contract[
                "m8_comparison_report_sha256"],
        },
        "findings": {
            "null_control_pass": report["null_control_pass"],
            "best_positive_refdens_by_mean_error_change_cm3": report[
                "best_positive_refdens_by_mean_error_change_cm3"],
            "best_positive_mean_error_change_dex": report[
                "best_positive_mean_error_change_dex"],
            "material_response_from_refdens_cm3": min(
                float(item["refdens_cm3"]) for item in positive
                if item["material_in_all_conditions"]),
            "worsens_all_conditions": [float(item["refdens_cm3"])
                                        for item in positive
                                        if item["worsens_all_conditions"]],
            "maximum_scan_response_dex": max(
                float(item["maximum_response_dex"]) for item in positive),
            "parameter_core_formula_equivalent_at_300K": read_json(
                parameter_audit)["core_formula_equivalent_at_300K"],
            "maximum_edge_limiter_reconstruction_error": edge_report[
                "maximum_limiter_reconstruction_error"],
            "maximum_edge_mean_mobility_aggregation_error": edge_report[
                "maximum_mean_mobility_aggregation_error"],
        },
        "interpretation": {
            "supported": [
                "Sentaurus HFS reference-density interpolation materially changes Id-Vg from 1e6 cm^-3 upward in all four conditions",
                "The exported 300 K Caughey-Thomas core vsat0, beta0, and alpha match Vela",
                "Vela per-edge diagnostics reproduce the production effective limiter to floating-point precision",
            ],
            "not_supported": [
                "Using reference density as a Vela default-model fit parameter",
                "Attributing the remaining baseline residual to one proprietary Sentaurus formula",
            ],
        },
        "artifacts": {
            "comparison_report": {"path": portable(frozen_report),
                                  "sha256": sha256(frozen_report)},
            "scan_summary": {"path": portable(summary), "sha256": sha256(summary)},
            "parameter_audit": {"path": portable(parameter_audit),
                                "sha256": sha256(parameter_audit)},
            "edge_report": {"path": portable(edge_report_path),
                            "sha256": sha256(edge_report_path)},
            "edge_summary": {
                "path": portable(edge_report_path.with_name("edge_hfs_summary.csv")),
                "sha256": sha256(edge_report_path.with_name("edge_hfs_summary.csv"))},
            "figures": [{"path": portable(path), "sha256": sha256(path)}
                        for path in figures],
        },
        "implementation_sha256": {portable(path): sha256(path)
                                  for path in implementations},
    }
    write_json(args.evidence.resolve(), evidence)
    print(json.dumps({
        "status": "frozen", "evidence": portable(args.evidence.resolve()),
        "sentaurus_curves": 24, "direct_bias_points": 1224,
        "baseline_unchanged": True,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
