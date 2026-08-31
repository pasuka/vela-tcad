#!/usr/bin/env python3
"""Freeze portable evidence for the SimpleMOS M30 causal-closure audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/double_off_causal_closure"
REPORT = ROOT / "m30_double_off_causal_closure_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m30_double_off_causal_closure_evidence.json")
ARTIFACTS = [
    "m30_double_off_causal_closure_report.json",
    "m30_state_operator_matrix.csv",
    "m30_native_gap_decomposition.csv",
    "m30_nodal_interaction.csv",
    "m30_nodal_interaction_summary.csv",
    "m30_solver_control_summary.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m30_double_off_causal_closure_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m29_bgn_srh_factorial_evidence.json",
    "scripts/run_simplemos_m30_double_off_causal_closure.py",
    "scripts/plot_simplemos_m30_double_off_causal_closure.py",
    "scripts/freeze_simplemos_m30_double_off_causal_closure_evidence.py",
    "tests/regression/test_simplemos_m30_double_off_causal_closure.py",
    "docs/validation/simplemos_m30_double_off_causal_closure_2026-08-31.md",
    "docs/validation/figures/simplemos_m30/simplemos_m30_double_off_causal_closure.png",
    "CMakeLists.txt",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M30 report is not complete and accepted")
    double_off = report["double_off_decomposition"]
    controls = report["solver_controls"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m30_double_off_causal_closure_evidence.v1",
        "status": "frozen",
        "summary": {
            "native_gap_dex": double_off["native_gap_dex"],
            "state_contribution_dex": double_off[
                "common_operator_state_contribution_dex"],
            "operator_extraction_contribution_dex": double_off[
                "sentaurus_state_vela_operator_contribution_dex"],
            "decomposition_closure_dex": double_off[
                "decomposition_closure_dex"],
            "solver_control_count": len(controls),
            "converged_solver_control_count": sum(
                bool(row["converged"]) for row in controls),
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/double_off_causal_closure/{name}",
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
