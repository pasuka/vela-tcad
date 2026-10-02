#!/usr/bin/env python3
"""Freeze portable M35 signed-AverageBox assembly evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/signed_average_box_assembly"
REPORT = ROOT / "m35_signed_average_box_assembly_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m35_signed_average_box_assembly_evidence.json")
ARTIFACTS = [
    "m35_signed_average_box_assembly_report.json",
    "m35_factorial_response.csv",
    "m35_factorial_main_effects.csv",
    "m35_interface_state_response.csv",
    "m35_sentaurus_measure_oracle.csv",
    "m35_m33_barycentric_comparison.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m35_signed_average_box_assembly_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/sentaurus_interface_box_probe/m34_interface_node_measures.csv",
    "scripts/run_simplemos_m35_signed_average_box_assembly.py",
    "scripts/plot_simplemos_m35_signed_average_box_assembly.py",
    "scripts/freeze_simplemos_m35_signed_average_box_assembly_evidence.py",
    "tests/regression/test_simplemos_m35_signed_average_box_assembly.py",
    "tests/test_mos_mixed_material.cpp",
    "tests/test_newton_solver.cpp",
    "docs/validation/simplemos_m35_signed_average_box_assembly_2026-08-31.md",
    "docs/validation/figures/simplemos_m35/simplemos_m35_signed_average_box_assembly.png",
    "include/vela/equation/AssemblerUtils.h",
    "include/vela/equation/DDAssembler.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/solver/NewtonSolver.cpp",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M35 report is not complete and accepted")
    all_on = report["all_region_resolved_variant"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m35_signed_average_box_assembly_evidence.v1",
        "status": "frozen",
        "summary": {
            "variant_count": report["acceptance"]["variant_count"],
            "converged_count": report["acceptance"]["converged_count"],
            "interface_node_count": report["sentaurus_measure_oracle"][
                "sentaurus_interface_node_count"],
            "maximum_measure_relative_error": report["sentaurus_measure_oracle"][
                "maximum_sentaurus_measure_relative_error"],
            "baseline_error_dex": report["factorial_response"][0]["absolute_error_dex"],
            "all_region_resolved_error_dex": all_on["absolute_error_dex"],
            "all_region_resolved_gap_closure_fraction": report["gap_closure"][
                "all_region_resolved_fraction"],
            "best_factorial_variant": report["best_variant"]["variant"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/signed_average_box_assembly/{name}",
            "sha256": sha256(ROOT / name),
        } for name in ARTIFACTS],
        "source_hashes": {name: sha256(REPO / name) for name in SOURCES},
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": str(OUTPUT),
                      **evidence["summary"]}))


if __name__ == "__main__":
    main()
