#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M26 first-layer feedback evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "first_layer_feedback_audit")
REPORT = ROOT / "m26_first_layer_feedback_audit_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m26_first_layer_feedback_audit_evidence.json")
ARTIFACTS = [
    "m26_first_layer_feedback_audit_report.json",
    "m26_operator_adjoint_ledger.csv",
    "m26_operator_zone_summary.csv",
    "m26_direction_summary.csv",
    "m26_first_layer_fd_ledger.csv",
    "m26_first_layer_jacobian_ledger.csv",
    "m26_first_layer_response_ledger.csv",
    "m26_perturbation_summary.csv",
]
UPSTREAM = [
    "reference_tcad/simplemos_sentaurus2022/contact_boundary_audit/m24_first_layer_continuity_ledger.csv",
    "reference_tcad/simplemos_sentaurus2022/contact_boundary_audit/m24_drain_cut_edge_ledger.csv",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m24_contact_boundary_audit_evidence.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m25_qf_reference_transform_audit_evidence.json",
]
FIGURE = ("docs/validation/figures/simplemos_m26/"
          "simplemos_m26_first_layer_feedback_audit.png")
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m26_first_layer_feedback_audit_contract_v1.json",
    "scripts/run_simplemos_m26_first_layer_feedback_audit.py",
    "scripts/plot_simplemos_m26_first_layer_feedback_audit.py",
    "scripts/freeze_simplemos_m26_first_layer_feedback_audit_evidence.py",
    "tests/regression/test_simplemos_m26_first_layer_feedback_audit.py",
    "tests/test_newton_solver.cpp",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "docs/validation/simplemos_m26_first_layer_feedback_audit_2026-08-31.md",
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
    closure = report["closure"]
    if execution["first_layer_nodes"] != [1089, 1090, 1092, 1099, 1100]:
        raise ValueError("M26 drain-first-layer topology changed")
    if execution["perturbation_column_count"] != 12:
        raise ValueError("M26 complete one-ring perturbation basis changed")
    if not closure["all_acceptance_checks_pass"]:
        raise ValueError("M26 numerical closure regressed")
    evidence = {
        "schema": "vela.simplemos.sdevice.m26_first_layer_feedback_audit_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only n23 drain-first-layer self-consistent feedback audit",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "operator_directions": report["operator_directions"],
        "perturbation_summaries": report["perturbation_summaries"],
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
        "first_layer_nodes": execution["first_layer_nodes"],
        "maximum_first_order_relative_error": closure[
            "maximum_first_order_relative_error"],
    }))


if __name__ == "__main__":
    main()
