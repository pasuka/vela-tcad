#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M27 cross-solver response evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "cross_solver_response_audit")
REPORT = ROOT / "m27_cross_solver_response_audit_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m27_cross_solver_response_audit_evidence.json")
SOURCE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m27_cross_solver_response_audit")
ARTIFACTS = [
    "m27_cross_solver_response_audit_report.json",
    "m27_node_response_ledger.csv",
    "m27_edge_response_ledger.csv",
    "m27_zone_summary.csv",
    "m27_terminal_response.csv",
]
UPSTREAM = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m22_n23_hfs_deep_off_evidence.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m26_first_layer_feedback_audit_evidence.json",
]
FIGURE = ("docs/validation/figures/simplemos_m27/"
          "simplemos_m27_cross_solver_response_audit.png")
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m27_cross_solver_response_audit_contract_v1.json",
    "scripts/run_simplemos_m27_cross_solver_response_audit.py",
    "scripts/plot_simplemos_m27_cross_solver_response_audit.py",
    "scripts/freeze_simplemos_m27_cross_solver_response_audit_evidence.py",
    "tests/regression/test_simplemos_m27_cross_solver_response_audit.py",
    "docs/validation/simplemos_m27_cross_solver_response_audit_2026-08-31.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(path: Path) -> dict[str, str]:
    return {"path": path.resolve().relative_to(REPO.resolve()).as_posix(),
            "sha256": sha256(path)}


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    execution = report["execution"]
    acceptance = report["acceptance"]
    response = report["terminal_response"]
    if execution["sentaurus_state_count"] != 2:
        raise ValueError("M27 Sentaurus state count changed")
    if execution["common_silicon_node_count"] != 942:
        raise ValueError("M27 common node topology changed")
    if not acceptance["all_checks_pass"]:
        raise ValueError("M27 data-integrity checks regressed")
    source_integrity = {
        "sentaurus_manifest_sha256": sha256(
            SOURCE / "sentaurus_manifest.json"),
        "sentaurus_export_manifest_sha256": sha256(
            SOURCE / "sentaurus_export_manifest.json"),
        "sentaurus_banner_sha256": sha256(SOURCE / "sentaurus_banner.txt")
        if (SOURCE / "sentaurus_banner.txt").is_file() else None,
        "sentaurus_states": {
            item["variant"]: {
                "tdr_sha256": item["tdr_sha256"],
                "current_file_sha256": item["current_file_sha256"],
                "field_manifest_sha256": item["field_manifest_sha256"],
            }
            for item in json.loads((SOURCE / "sentaurus_export_manifest.json")
                                   .read_text(encoding="utf-8"))["states"]
        },
    }
    evidence = {
        "schema": "vela.simplemos.sdevice.m27_cross_solver_response_audit_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only n23 Sentaurus-Vela HFS self-consistent response audit",
        "sentaurus_release": execution["sentaurus_release"],
        "execution": execution,
        "acceptance": acceptance,
        "terminal_response": response,
        "cross_solver_response": report["cross_solver_response"],
        "zones": report["zones"],
        "conclusions": report["conclusions"],
        "claim_policy": report["claim_policy"],
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "upstream_inputs": [entry(REPO / name) for name in UPSTREAM],
        "figures": [entry(REPO / FIGURE)],
        "source_integrity": source_integrity,
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n",
                        encoding="utf-8", newline="\n")
    print(json.dumps({
        "status": "frozen",
        "output": entry(EVIDENCE)["path"],
        "sentaurus_states": execution["sentaurus_state_count"],
        "terminal_response_relative_disagreement": response[
            "response_relative_disagreement"],
    }))


if __name__ == "__main__":
    main()
