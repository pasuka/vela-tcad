#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M24 contact-boundary evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/contact_boundary_audit"
REPORT = ROOT / "m24_contact_boundary_audit_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m24_contact_boundary_audit_evidence.json")
ARTIFACTS = ["m24_contact_boundary_audit_report.json",
             "m24_contact_node_ledger.csv", "m24_drain_cut_edge_ledger.csv",
             "m24_first_layer_continuity_ledger.csv"]
UPSTREAM = [
    "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off/m22_summary.csv",
    "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off/m22_key_state_self_consistent_edges.csv",
]
FIGURE = "docs/validation/figures/simplemos_m24/simplemos_m24_contact_boundary_audit.png"
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m24_contact_boundary_audit_contract_v1.json",
    "scripts/run_simplemos_m24_contact_boundary_audit.py",
    "scripts/plot_simplemos_m24_contact_boundary_audit.py",
    "scripts/freeze_simplemos_m24_contact_boundary_audit_evidence.py",
    "tests/regression/test_simplemos_m24_contact_boundary_audit.py",
    "docs/validation/simplemos_m24_contact_boundary_audit_2026-08-30.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(path: Path) -> dict[str, str]:
    return {"path": path.resolve().relative_to(REPO.resolve()).as_posix(),
            "sha256": sha256(path)}


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    closure = report["closure"]
    findings = report["findings"]
    if report["execution"]["read_only_probe_count"] != 4:
        raise ValueError("M24 probe matrix is incomplete")
    if closure["maximum_contact_qf_bias_error_V"] >= 1.0e-12:
        raise ValueError("M24 contact QF target regressed")
    if closure["maximum_sg_cut_current_relative_disagreement"] >= 1.0e-3:
        raise ValueError("M24 SG/current extraction disagreement regressed")
    if findings["maximum_full_no_hfs_contact_physical_state_delta"] != 0.0:
        raise ValueError("M24 HFS variants no longer share the contact state")
    evidence = {
        "schema": "vela.simplemos.sdevice.m24_contact_boundary_audit_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only n23 drain-contact boundary audit",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": report["execution"],
        "closure": closure,
        "findings": findings,
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
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "contact_nodes": 3, "cut_edges": 7,
                      "first_layer_nodes": 5}))


if __name__ == "__main__":
    main()
