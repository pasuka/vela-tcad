#!/usr/bin/env python3
"""Freeze portable evidence for the M33 interface-assembly audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "region_resolved_interface_assembly")
REPORT = ROOT / "m33_region_resolved_interface_assembly_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m33_region_resolved_interface_assembly_evidence.json")
ARTIFACTS = [
    "m33_region_resolved_interface_assembly_report.json",
    "m33_factorial_response.csv",
    "m33_factorial_main_effects.csv",
    "m33_interface_geometry.csv",
    "m33_interface_state_response.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m33_region_resolved_interface_assembly_contract_v1.json",
    "scripts/run_simplemos_m33_region_resolved_interface_assembly.py",
    "scripts/plot_simplemos_m33_region_resolved_interface_assembly.py",
    "scripts/freeze_simplemos_m33_region_resolved_interface_assembly_evidence.py",
    "tests/regression/test_simplemos_m33_region_resolved_interface_assembly.py",
    "tests/test_mos_mixed_material.cpp",
    "docs/validation/simplemos_m33_region_resolved_interface_assembly_2026-08-31.md",
    "docs/validation/figures/simplemos_m33/simplemos_m33_region_resolved_interface_assembly.png",
    "include/vela/equation/AssemblerUtils.h",
    "include/vela/equation/DDAssembler.h",
    "include/vela/solver/NewtonSolver.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/solver/NewtonSolver.cpp",
    "src/simulation/DCSweep.cpp",
    "CMakeLists.txt",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M33 report is not complete and accepted")
    all_on = report["all_region_resolved_variant"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m33_region_resolved_interface_assembly_evidence.v1",
        "status": "frozen",
        "summary": {
            "variant_count": report["acceptance"]["variant_count"],
            "interface_edge_count": report["geometry"][
                "si_sio2_interface_edge_count"],
            "baseline_error_dex": report["factorial_response"][0][
                "absolute_error_dex"],
            "all_region_resolved_error_dex": all_on["absolute_error_dex"],
            "all_region_resolved_gap_closure_fraction": report[
                "gap_closure"]["all_region_resolved_fraction"],
            "best_factorial_variant": report["best_variant"]["variant"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/region_resolved_interface_assembly/{name}",
            "sha256": sha(ROOT / name),
        } for name in ARTIFACTS],
        "source_hashes": {name: sha(REPO / name) for name in SOURCES},
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": str(OUTPUT),
                      **evidence["summary"]}))


if __name__ == "__main__":
    main()
