#!/usr/bin/env python3
"""Freeze portable evidence for the SimpleMOS M28 path audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/exact_bias_path_audit"
REPORT = ROOT / "m28_exact_bias_path_audit_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m28_exact_bias_path_audit_evidence.json")
ARTIFACTS = [
    "m28_exact_bias_path_audit_report.json",
    "m28_terminal_route_ledger.csv",
    "m28_state_route_ledger.csv",
    "m28_state_route_summary.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m28_exact_bias_path_audit_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m27_cross_solver_response_audit_evidence.json",
    "scripts/run_simplemos_m28_exact_bias_path_audit.py",
    "scripts/plot_simplemos_m28_exact_bias_path_audit.py",
    "scripts/freeze_simplemos_m28_exact_bias_path_audit_evidence.py",
    "tests/regression/test_simplemos_m28_exact_bias_path_audit.py",
    "docs/validation/simplemos_m28_exact_bias_path_audit_2026-08-31.md",
    "docs/validation/figures/simplemos_m28/simplemos_m28_exact_bias_path_audit.png",
    "CMakeLists.txt",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M28 report is not complete and accepted")
    evidence = {
        "schema": "vela.simplemos.sdevice.m28_exact_bias_path_audit_evidence.v1",
        "status": "frozen",
        "summary": {
            "new_sentaurus_state_count": report["execution"]["new_sentaurus_state_count"],
            "maximum_hfs_response_relative_difference_from_m8": report["terminal"]["maximum_hfs_response_relative_difference_from_m8"],
            "maximum_common_mode_offset_A_per_um": report["terminal"]["maximum_common_mode_offset_A_per_um"],
            "maximum_differential_offset_A_per_um": report["terminal"]["maximum_differential_offset_A_per_um"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{"path": f"reference_tcad/simplemos_sentaurus2022/exact_bias_path_audit/{name}",
                       "sha256": sha(ROOT / name)} for name in ARTIFACTS],
        "source_hashes": {name: sha(REPO / name) for name in SOURCES},
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": str(OUTPUT), **evidence["summary"]}))


if __name__ == "__main__":
    main()
