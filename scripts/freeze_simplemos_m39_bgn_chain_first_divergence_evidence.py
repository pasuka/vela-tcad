#!/usr/bin/env python3
"""Freeze portable M39 BGN-chain first-divergence evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/bgn_chain_first_divergence"
REPORT = ROOT / "m39_bgn_chain_first_divergence_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m39_bgn_chain_first_divergence_evidence.json")
ARTIFACTS = [
    "m39_bgn_chain_first_divergence_report.json",
    "m39_model_selection_ledger.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m39_bgn_chain_first_divergence_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m27_cross_solver_response_audit_evidence.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m29_bgn_srh_factorial_evidence.json",
    "scripts/run_simplemos_m39_bgn_chain_first_divergence.py",
    "scripts/freeze_simplemos_m39_bgn_chain_first_divergence_evidence.py",
    "tests/regression/test_simplemos_m39_bgn_chain_first_divergence.py",
    "docs/validation/simplemos_m39_bgn_chain_first_divergence_2026-09-01.md",
    "CMakeLists.txt"
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M39 report is not complete and accepted")
    first = report["first_divergence"]
    if (first["stage"] != "C-1_model_selection_contract" or
            first["observed"] != "BennettWilson"):
        raise RuntimeError("M39 did not preserve the qualified first divergence")
    evidence = {
        "schema": "vela.simplemos.sdevice.m39_bgn_chain_first_divergence_evidence.v1",
        "status": "frozen",
        "summary": {
            "first_divergence_stage": first["stage"],
            "m29_declared_model": first["declared"],
            "m29_effective_model": first["observed"],
            "m27_c1_total_max_abs_difference_eV": report[
                "m27_qualified_oldslotboom_chain"][
                    "c1_total_impurity_max_abs_difference_eV"],
            "m27_c3_max_joint_fit_rms": report[
                "m27_qualified_oldslotboom_chain"]["c3_max_joint_fit_rms"],
            "default_model_changed": False,
            "sentaurus_rerun": False,
            "vela_solve_rerun": False
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/bgn_chain_first_divergence/{name}",
            "sha256": sha256(ROOT / name)
        } for name in ARTIFACTS],
        "source_hashes": {name: sha256(REPO / name) for name in SOURCES}
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": str(OUTPUT),
                      **evidence["summary"]}))


if __name__ == "__main__":
    main()
