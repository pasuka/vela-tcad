#!/usr/bin/env python3
"""Freeze compact, portable evidence for the SimpleMOS M11 mobility factorial."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                  / "m11_mobility_factorial")
DEFAULT_DESTINATION = (REPO / "reference_tcad/simplemos_sentaurus2022"
                       / "mobility_factorial")
DEFAULT_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "simplemos_m11_mobility_factorial_evidence.json")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def copy_file(source: Path, destination: Path) -> dict[str, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return {"path": portable(destination), "sha256": sha256(destination)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    source = args.source.resolve()
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)

    self_report = read_json(source / "self_consistent_factorial_report.json")
    comparison = read_json(source / "comparisons/comparison_report.json")
    frozen_report = read_json(source / "frozen_vela/frozen_factorial_report.json")
    reference_manifest = read_json(source / "sentaurus_reference/reference_manifest.json")

    artifacts: dict[str, Any] = {}
    artifacts["self_consistent_report"] = copy_file(
        source / "self_consistent_factorial_report.json",
        destination / "self_consistent_factorial_report.json")
    artifacts["self_consistent_effects"] = copy_file(
        source / "self_consistent_factorial_effects.csv",
        destination / "self_consistent_factorial_effects.csv")
    artifacts["frozen_report"] = copy_file(
        source / "frozen_vela/frozen_factorial_report.json",
        destination / "frozen_factorial_report.json")
    artifacts["frozen_metrics"] = copy_file(
        source / "frozen_vela/frozen_factorial_metrics.csv",
        destination / "frozen_factorial_metrics.csv")
    artifacts["frozen_effects"] = copy_file(
        source / "frozen_vela/frozen_factorial_effects.csv",
        destination / "frozen_factorial_effects.csv")

    frozen_comparison = dict(comparison)
    frozen_cases = []
    comparison_artifacts = []
    for case in comparison["cases"]:
        item = dict(case)
        source_csv = source / "comparisons" / Path(case["comparison_csv"]).name
        target_csv = destination / "comparisons" / source_csv.name
        frozen = copy_file(source_csv, target_csv)
        item["comparison_csv"] = frozen["path"]
        item["comparison_csv_sha256"] = frozen["sha256"]
        # Raw curves are intentionally represented by the qualified manifest
        # and the aligned comparison CSV instead of duplicated here.
        item["reference"] = f"sentaurus-manifest:{case['case']}"
        item["candidate"] = f"vela-self-consistent:{case['case']}"
        frozen_cases.append(item)
        comparison_artifacts.append(frozen)
    frozen_comparison["cases"] = frozen_cases
    comparison_target = destination / "comparison_report.json"
    write_json(comparison_target, frozen_comparison)
    artifacts["comparison_report"] = {
        "path": portable(comparison_target), "sha256": sha256(comparison_target)}
    artifacts["comparison_curves"] = comparison_artifacts

    figures = sorted((REPO / "docs/validation/figures/simplemos_m11").glob("*.png"))
    artifacts["figures"] = [
        {"path": portable(path), "sha256": sha256(path)} for path in figures]
    report_dir = REPO / "docs/validation/reports/simplemos_m11"
    report_files = [report_dir / "artifact.json", report_dir / "CHART_MAP.md",
                    report_dir / "simplemos_m11_mobility_factorial_report.html"]
    artifacts["portable_report"] = [
        {"path": portable(path), "sha256": sha256(path)}
        for path in report_files if path.is_file()]

    cases = comparison["cases"]
    data_quality = {
        "curve_count": len(cases),
        "direct_bias_points_per_solver": sum(int(case["point_count"]) for case in cases),
        "points_above_floor": sum(int(case["points_above_floor"]) for case in cases),
        "all_cases_pass": all(case["status"] == "pass" for case in cases),
        "all_trends_match": all(bool(case["trend_match"]) for case in cases),
        "all_finite": all(not case["failures"] for case in cases),
        "frozen_state_count": frozen_report["state_count"],
        "frozen_variant_evaluations": frozen_report["variant_evaluations"],
    }

    contract = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m11_mobility_factorial_contract_v1.json")
    implementations = [
        contract,
        REPO / "CMakeLists.txt",
        REPO / "src/physics/MobilityModel.cpp",
        REPO / "src/equation/CoupledDDAssembler.cpp",
        REPO / "include/vela/equation/AssemblerUtils.h",
        REPO / "include/vela/equation/ElementEdgeGssLauxAD.inl",
        REPO / "tests/test_mobility.cpp",
        REPO / "scripts/run_simplemos_m11_mobility_factorial.py",
        REPO / "scripts/plot_simplemos_m11_mobility_factorial.py",
        REPO / "scripts/freeze_simplemos_m11_mobility_factorial_evidence.py",
        REPO / "scripts/build_simplemos_m11_report_artifact.py",
        REPO / "tests/regression/test_simplemos_m11_mobility_factorial.py",
    ]
    evidence = {
        "schema": "vela.simplemos.sdevice.m11_mobility_factorial_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only 2^3 PhuMob x Enormal x HFS factorial",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": {
            "devices": ["n17", "n21"],
            "drain_voltages_V": [0.05, 1.0],
            "gate_sweep_V": {"start": 0.0, "step": 0.05, "stop": 2.5},
            "variants": 8,
            "self_consistent_curves": 32,
            "interpolation": "forbidden",
            "default_model_changed": False,
        },
        "data_quality": data_quality,
        "findings": {
            "self_consistent_aggregate_effects_dex": self_report["aggregate_effects"],
            "self_consistent_variant_aggregate_dex": self_report["variant_aggregate"],
            "self_consistent_condition_effects_dex": self_report["conditions"],
            "frozen_vela_aggregate_effects": frozen_report["aggregate_effects"],
            "frozen_vela_variant_aggregate": frozen_report["variant_aggregate"],
        },
        "interpretation": {
            "supported": [
                "In paired self-consistent maximum Id error, Enormal improves agreement in all four device/drain conditions",
                "In paired self-consistent maximum Id error, HFS increases mismatch in all four conditions and has the largest positive aggregate main effect",
                "PhuMob strongly improves the n21 low-drain condition but is nearly neutral at high drain",
                "On the same full-physics Sentaurus state, enabling all three Vela mobility components reduces mobility residual, with PhuMob the dominant low-field contribution",
                "Opposite HFS signs in self-consistent Id error and frozen-state mobility error identify drive construction, limiting, discretization, and nonlinear state feedback as the remaining investigation boundary",
            ],
            "not_supported": [
                "Treating the frozen Vela matrix as a frozen Sentaurus model re-evaluation",
                "Attributing all remaining mismatch to a single private Sentaurus HFS formula",
                "Changing the default Vela mobility stack from this diagnostic matrix alone",
            ],
        },
        "artifacts": artifacts,
        "source_integrity": {
            "sentaurus_reference_manifest_sha256": sha256(
                source / "sentaurus_reference/reference_manifest.json"),
            "sentaurus_reference_curve_count": len(reference_manifest["artifacts"]),
            "sentaurus_reference_all_51_points": all(
                int(item["point_count"]) == 51 for item in reference_manifest["artifacts"]),
            "raw_self_consistent_report_sha256": sha256(
                source / "self_consistent_factorial_report.json"),
            "raw_frozen_report_sha256": sha256(
                source / "frozen_vela/frozen_factorial_report.json"),
        },
        "implementation_sha256": {
            portable(path): sha256(path) for path in implementations},
    }
    write_json(args.evidence.resolve(), evidence)
    print(json.dumps({"status": "frozen", "curves": len(cases),
                      "states": frozen_report["state_count"],
                      "evidence": portable(args.evidence.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
