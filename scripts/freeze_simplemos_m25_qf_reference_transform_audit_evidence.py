#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M25 QF-reference transform evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "qf_reference_transform_audit")
REPORT = ROOT / "m25_qf_reference_transform_audit_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m25_qf_reference_transform_audit_evidence.json")
ARTIFACTS = [
    "m25_qf_reference_transform_audit_report.json",
    "m25_edge_reference_transform_ledger.csv",
    "m25_variant_edge_delta_ledger.csv",
    "m25_class_summary.csv",
]
UPSTREAM = [
    "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off/m22_summary.csv",
    "reference_tcad/simplemos_sentaurus2022/contact_boundary_audit/m24_drain_cut_edge_ledger.csv",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m23_numerical_resolvability_evidence.json",
]
FIGURE = ("docs/validation/figures/simplemos_m25/"
          "simplemos_m25_qf_reference_transform_audit.png")
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m25_qf_reference_transform_audit_contract_v1.json",
    "scripts/run_simplemos_m25_qf_reference_transform_audit.py",
    "scripts/plot_simplemos_m25_qf_reference_transform_audit.py",
    "scripts/freeze_simplemos_m25_qf_reference_transform_audit_evidence.py",
    "tests/regression/test_simplemos_m25_qf_reference_transform_audit.py",
    "docs/validation/simplemos_m25_qf_reference_transform_audit_2026-08-30.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(path: Path) -> dict[str, str]:
    return {"path": path.resolve().relative_to(REPO.resolve()).as_posix(),
            "sha256": sha256(path)}


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    topology = report["topology"]
    closure = report["closure"]
    if topology["reference_transition_edge_count"] != 51:
        raise ValueError("M25 reference-transition topology changed")
    if topology["drain_cut_reference_transition_edge_count"] != 0:
        raise ValueError("M25 drain cut unexpectedly crosses a reference basin")
    if not closure["all_acceptance_checks_pass"]:
        raise ValueError("M25 numerical closure regressed")
    evidence = {
        "schema": "vela.simplemos.sdevice.m25_qf_reference_transform_audit_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only n23 QF-reference and SG-log audit",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": report["execution"],
        "topology": topology,
        "closure": closure,
        "findings": report["findings"],
        "conclusions": report["conclusions"],
        "claim_policy": report["claim_policy"],
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "upstream_inputs": [entry(REPO / name) for name in UPSTREAM],
        "figures": [entry(REPO / FIGURE)],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n",
                        encoding="utf-8", newline="\n")
    print(json.dumps({
        "status": "frozen",
        "output": entry(EVIDENCE)["path"],
        "selected_edges": topology["unique_selected_edges"],
        "reference_transition_edges": topology[
            "reference_transition_edge_count"],
    }))


if __name__ == "__main__":
    main()
