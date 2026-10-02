#!/usr/bin/env python3
"""Freeze portable M36 boundary/contact measure evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/boundary_contact_measure_audit"
REPORT = ROOT / "m36_boundary_contact_measure_audit_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m36_boundary_contact_measure_audit_evidence.json")
ARTIFACTS = [
    "m36_boundary_contact_measure_audit_report.json",
    "m36_local_measure_ledger.csv",
    "m36_node_measure_ledger.csv",
    "m36_node_class_summary.csv",
    "m36_contact_support_summary.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m36_boundary_contact_measure_audit_contract_v1.json",
    "scripts/run_simplemos_m36_boundary_contact_measure_audit.py",
    "scripts/freeze_simplemos_m36_boundary_contact_measure_audit_evidence.py",
    "tests/regression/test_simplemos_m36_boundary_contact_measure_audit.py",
    "docs/validation/simplemos_m36_boundary_contact_measure_audit_2026-08-31.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M36 report is not complete and accepted")
    evidence = {
        "schema": "vela.simplemos.sdevice.m36_boundary_contact_measure_audit_evidence.v1",
        "status": "frozen",
        "summary": {
            "node_count": report["mesh"]["node_count"],
            "triangle_count": report["mesh"]["triangle_count"],
            "maximum_si_node_relative_error": report["maximum_si_node_relative_error"],
            "maximum_local_relative_error": report["maximum_local_relative_error"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/boundary_contact_measure_audit/{name}",
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
