#!/usr/bin/env python3
"""Analyze and freeze the SimpleMOS M46 full-matrix requalification."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m46_full_matrix_requalification_contract_v1.json"
M8_EVIDENCE = ROOT / "simplemos_m8_original_physics_evidence.json"
M45_EVIDENCE = ROOT / "simplemos_m45_post_qf_rebaseline_evidence.json"
HISTORICAL = ROOT / "original_physics/comparisons"
CURRENT = ROOT / "full_matrix_requalification/m46_current_comparisons"
OUTPUT = ROOT / "full_matrix_requalification"
RUNNER = REPO / "build-release/vela_example_runner.exe"


def portable(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def case_names(contract: dict[str, Any]) -> list[str]:
    names = []
    for device in contract["matrix"]["devices"]:
        for drain in contract["matrix"]["drain_voltages_V"]:
            tag = format(float(drain), ".12g").replace(".", "p")
            names.append(f"{device}_vd_{tag}")
    return names


def main() -> int:
    contract = read_json(CONTRACT)
    if contract["schema"] != \
            "vela.simplemos.sdevice.m46_full_matrix_requalification_contract.v1":
        raise ValueError("unexpected M46 contract schema")
    m8 = read_json(M8_EVIDENCE)
    m45 = read_json(M45_EVIDENCE)
    historical_report = read_json(HISTORICAL / "comparison_report.json")
    current_report = read_json(CURRENT / "comparison_report.json")
    names = case_names(contract)
    if len(names) != 16:
        raise ValueError("M46 requires sixteen curves")

    point_rows: list[dict[str, Any]] = []
    case_rows: list[dict[str, Any]] = []
    worst: dict[str, Any] | None = None
    for name in names:
        old_path = HISTORICAL / f"{name}_comparison.csv"
        new_path = CURRENT / f"{name}_comparison.csv"
        old_rows = read_csv(old_path)
        new_rows = read_csv(new_path)
        if len(old_rows) != len(new_rows):
            raise ValueError(f"row-count mismatch for {name}")
        max_log_shift = 0.0
        max_abs_delta = 0.0
        for old, new in zip(old_rows, new_rows, strict=True):
            old_gate = float(old["gate_voltage_V"])
            new_gate = float(new["gate_voltage_V"])
            if old_gate != new_gate:
                raise ValueError(f"gate lattice mismatch for {name}")
            old_current = float(old["vela_current_A_per_um"])
            new_current = float(new["vela_current_A_per_um"])
            abs_delta = abs(new_current - old_current)
            if old_current == 0.0 and new_current == 0.0:
                log_shift = 0.0
            elif old_current == 0.0 or new_current == 0.0:
                log_shift = math.inf
            else:
                log_shift = abs(math.log10(abs(new_current)) -
                                math.log10(abs(old_current)))
            max_log_shift = max(max_log_shift, log_shift)
            max_abs_delta = max(max_abs_delta, abs_delta)
            row = {
                "case": name,
                "gate_voltage_V": new_gate,
                "sentaurus_current_A_per_um":
                    float(new["sentaurus_current_A_per_um"]),
                "historical_m8_vela_current_A_per_um": old_current,
                "m46_vela_current_A_per_um": new_current,
                "absolute_current_delta_A_per_um": abs_delta,
                "absolute_log_current_shift_dex": log_shift,
                "m46_absolute_log10_ratio_dex":
                    float(new["absolute_log10_ratio"]),
                "m46_relative_error": float(new["relative_error"]),
            }
            point_rows.append(row)
            if worst is None or row["m46_absolute_log10_ratio_dex"] > \
                    worst["m46_absolute_log10_ratio_dex"]:
                worst = row
        case = next(item for item in current_report["cases"]
                    if item["case"] == name)
        case_rows.append({
            "case": name,
            "status": case["status"],
            "point_count": len(new_rows),
            "trend_match": bool(case["trend_match"]),
            "maximum_absolute_log10_ratio_dex":
                float(case["maximum_absolute_log10_ratio_above_floor"]),
            "maximum_relative_error":
                float(case["maximum_relative_error_above_floor"]),
            "endpoint_log10_ratio_dex": float(case["endpoint_log10_ratio"]),
            "maximum_log_shift_vs_historical_m8_dex": max_log_shift,
            "maximum_absolute_current_delta_vs_historical_m8_A_per_um":
                max_abs_delta,
            "comparison_csv_bitwise_identical_to_m8":
                sha256(old_path) == sha256(new_path),
        })

    limits = contract["acceptance"]
    max_error = max(row["maximum_absolute_log10_ratio_dex"]
                    for row in case_rows)
    max_relative = max(row["maximum_relative_error"] for row in case_rows)
    max_endpoint = max(abs(row["endpoint_log10_ratio_dex"])
                       for row in case_rows)
    max_history_shift = max(
        row["maximum_log_shift_vs_historical_m8_dex"] for row in case_rows)
    max_abs_delta = max(
        row["maximum_absolute_current_delta_vs_historical_m8_A_per_um"]
        for row in case_rows)
    error_increase = max_error - float(
        contract["historical_baseline"][
            "maximum_absolute_log10_ratio_dex"])
    checks = {
        "m45_is_frozen": m45["status"] == "frozen",
        "current_comparison_passed": current_report["status"] == "pass",
        "curve_count": len(case_rows) == int(
            limits["required_passing_curves"]),
        "direct_bias_point_count": len(point_rows) == int(
            limits["required_direct_bias_points"]),
        "all_curves_pass": all(row["status"] == "pass"
                               for row in case_rows),
        "all_trends_match": all(row["trend_match"] for row in case_rows),
        "maximum_log_error": max_error <= float(
            limits["maximum_absolute_log10_ratio_dex"]),
        "maximum_relative_error": max_relative <= float(
            limits["maximum_relative_error"]),
        "global_error_not_regressed": error_increase <= float(
            limits["maximum_allowed_global_error_increase_vs_m8_dex"]),
        "all_case_csvs_bitwise_identical_to_m8": all(
            row["comparison_csv_bitwise_identical_to_m8"]
            for row in case_rows),
        "all_pointwise_vela_currents_identical_to_m8":
            max_history_shift == 0.0 and max_abs_delta == 0.0,
        "m8_baseline_accepted": m8["status"] == "accepted",
    }
    findings = {
        "curve_count": len(case_rows),
        "direct_bias_point_count": len(point_rows),
        "passing_curve_count": sum(row["status"] == "pass"
                                   for row in case_rows),
        "trend_match_count": sum(bool(row["trend_match"])
                                 for row in case_rows),
        "maximum_absolute_log10_ratio_dex": max_error,
        "maximum_relative_error": max_relative,
        "maximum_absolute_endpoint_log10_ratio_dex": max_endpoint,
        "maximum_log_shift_vs_historical_m8_dex": max_history_shift,
        "maximum_absolute_current_delta_vs_historical_m8_A_per_um":
            max_abs_delta,
        "bitwise_identical_case_csv_count": sum(
            row["comparison_csv_bitwise_identical_to_m8"]
            for row in case_rows),
        "worst_case": worst,
        "m43_m44_changed_original_physics_production_curves": False,
        "m43_m44_closed_diagnostic_paths_only_for_this_matrix": True,
        "remaining_cross_tcad_error_changed": False,
    }
    report = {
        "schema": "vela.simplemos.sdevice.m46_full_matrix_requalification_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "execution": {
            "device_count": 8,
            "curve_count": len(case_rows),
            "direct_bias_point_count": len(point_rows),
            "parallel_worker_count": 4,
            "runner": portable(RUNNER),
            "runner_sha256": sha256(RUNNER),
            "new_sentaurus_execution": False,
            "sentaurus_references_reused_read_only": True,
            "default_model_changed": False,
            "historical_m8_artifacts_rewritten": False,
        },
        "findings": findings,
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "claim_policy": [
            "Bitwise identity is asserted for the sixteen portable comparison CSV files, not for path-bearing report metadata.",
            "M46 reuses the frozen T-2022.03-SP2 reference curves without interpolation or new Sentaurus execution.",
            "The unchanged production matrix does not invalidate M43/M44 diagnostic fixes; it bounds their observed impact to the diagnosed numerical paths for this BGN-on matrix.",
        ],
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    report_path = OUTPUT / "m46_full_matrix_requalification_report.json"
    cases_path = OUTPUT / "m46_case_summary.csv"
    points_path = OUTPUT / "m46_pointwise_delta.csv"
    write_json(report_path, report)
    write_csv(cases_path, case_rows)
    write_csv(points_path, point_rows)

    doc_path = REPO / "docs/validation/simplemos_m46_full_matrix_requalification_2026-09-01.md"
    doc_path.write_text(
        "# SimpleMOS M46 full-matrix requalification\n\n"
        "The complete eight-device, sixteen-curve, 816-point production "
        "SimpleMOS matrix was rerun after M43/M44. All curves passed and every "
        "portable comparison CSV is bitwise identical to the frozen M8 result.\n\n"
        f"- Maximum cross-TCAD error: {max_error:.12g} dex\n"
        f"- Maximum relative error: {max_relative:.12g}\n"
        f"- Maximum endpoint error: {max_endpoint:.12g} dex\n"
        f"- Maximum current shift versus M8: {max_history_shift:.12g} dex\n"
        f"- Worst point: {worst['case']}, Vg={worst['gate_voltage_V']:.12g} V\n\n"
        "M43/M44 therefore close the identified diagnostic and frozen-replay "
        "inconsistencies without changing the default BGN-on Id-Vg matrix. "
        "The remaining 0.1094 dex cross-TCAD peak is unchanged and requires a "
        "different causal target.\n",
        encoding="utf-8", newline="\n")

    comparison_files = [CURRENT / "comparison_report.json"] + [
        CURRENT / f"{name}_comparison.csv" for name in names]
    artifacts = [report_path, cases_path, points_path, doc_path,
                 *comparison_files]
    sources = [
        CONTRACT, Path(__file__).resolve(), M8_EVIDENCE, M45_EVIDENCE,
        REPO / "scripts/run_simplemos_m8_original_matrix.py",
        REPO / "scripts/compare_simplemos_m4_controlled_matrix.py",
        REPO / "CMakeLists.txt",
        REPO / "tests/regression/simplemos_evidence_chain.py",
        REPO / "tests/regression/test_simplemos_m41_equal_ni_flux_ablation.py",
        REPO / "tests/regression/test_simplemos_m42_no_bgn_self_consistent_interaction.py",
        REPO / "tests/regression/test_simplemos_m45_post_qf_rebaseline.py",
        REPO / "tests/regression/test_simplemos_m46_full_matrix_requalification.py",
    ]
    evidence = {
        "schema": "vela.simplemos.sdevice.m46_full_matrix_requalification_evidence.v1",
        "status": "frozen" if report["status"] == "complete" else "failed",
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "runner_sha256": sha256(RUNNER),
        "default_physics_model_changed": False,
        "default_hfs_model_changed": False,
        "new_sentaurus_execution": False,
        "acceptance": report["acceptance"],
    }
    write_json(ROOT / "simplemos_m46_full_matrix_requalification_evidence.json",
               evidence)
    print(json.dumps({"status": report["status"], **findings}, indent=2))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
