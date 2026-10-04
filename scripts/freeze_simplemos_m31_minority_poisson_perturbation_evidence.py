#!/usr/bin/env python3
"""Freeze portable evidence for the SimpleMOS M31 perturbation audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/minority_poisson_perturbation"
REPORT = ROOT / "m31_minority_poisson_perturbation_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m31_minority_poisson_perturbation_evidence.json")
ARTIFACTS = [
    "m31_minority_poisson_perturbation_report.json",
    "m31_state_variant_response.csv",
    "m31_poisson_node_audit.csv",
    "m31_cross_block_summary.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m31_minority_poisson_perturbation_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m30_double_off_causal_closure_evidence.json",
    "include/vela/solver/NewtonSolver.h",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "scripts/run_simplemos_m31_minority_poisson_perturbation.py",
    "scripts/plot_simplemos_m31_minority_poisson_perturbation.py",
    "scripts/freeze_simplemos_m31_minority_poisson_perturbation_evidence.py",
    "tests/regression/test_simplemos_m31_minority_poisson_perturbation.py",
    "docs/validation/simplemos_m31_minority_poisson_perturbation_2026-08-31.md",
    "docs/validation/figures/simplemos_m31/simplemos_m31_minority_poisson_perturbation.png",
    "CMakeLists.txt",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M31 report is not complete and accepted")
    evidence = {
        "schema": "vela.simplemos.sdevice.m31_minority_poisson_perturbation_evidence.v1",
        "status": "frozen",
        "summary": {
            "hole_full_current_shift_dex": report["minority_hole"][
                "full_current_shift_dex"],
            "poisson_top10_l2_share": report["poisson_floor"]["top10_l2_share"],
            "electron_full_current_shift_dex": next(
                row["absolute_log10_current_shift_dex"]
                for row in report["state_variants"]
                if row["variant"] == "electron_full_replace"),
            "condition_estimates_skipped_count": report["acceptance"][
                "condition_estimates_skipped_count"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/minority_poisson_perturbation/{name}",
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
