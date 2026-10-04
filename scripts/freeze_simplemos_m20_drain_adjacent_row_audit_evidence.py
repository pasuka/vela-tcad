#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M20 drain-adjacent row evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/drain_adjacent_row_audit"
REPORT = ROOT / "m20_drain_adjacent_row_audit_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m20_drain_adjacent_row_audit_evidence.json")
ARTIFACTS = [
    "m20_drain_adjacent_row_audit_report.json",
    "m20_state_summary.csv",
    "m20_node_balance.csv",
    "m20_node_delta.csv",
    "m20_incident_edge_jacobian.csv",
    "m20_key_state_node_ledger.csv",
]
FIGURES = [
    "m20_key_state_flux_balance.png",
    "m20_key_state_flux_conditioning.png",
    "m20_key_state_jacobian_partition.png",
]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m20_drain_adjacent_row_audit_contract_v1.json",
    "include/vela/solver/NewtonSolver.h",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "scripts/run_simplemos_m20_drain_adjacent_row_audit.py",
    "scripts/plot_simplemos_m20_drain_adjacent_row_audit.py",
    "scripts/freeze_simplemos_m20_drain_adjacent_row_audit_evidence.py",
    "tests/regression/test_simplemos_m20_drain_adjacent_row_audit.py",
    "docs/validation/simplemos_m20_drain_adjacent_row_audit_2026-08-30.md",
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
    key = report["findings"]["key_state"]
    if execution["state_count"] != 16 or execution["probe_count"] != 96:
        raise ValueError("M20 state/probe matrix is incomplete")
    if closure["maximum_raw_edge_flux_closure"] >= 1.0e-10:
        raise ValueError("M20 incident-edge flux closure failed")
    if closure["maximum_raw_term_sum_closure"] >= 1.0e-10:
        raise ValueError("M20 raw continuity-term closure failed")
    if closure["maximum_exact_jacobian_partition_relative_error"] >= 1.0e-12:
        raise ValueError("M20 exact Jacobian partition closure failed")
    if execution["new_sentaurus_execution"] or execution["default_model_changed"]:
        raise ValueError("M20 violated its read-only physical-model contract")
    if key["absolute_srh_delta_share"] >= 1.0e-15:
        raise ValueError("M20 key-state SRH exclusion changed")
    if not 0.49 < key["absolute_contact_flux_delta_share"] < 0.51:
        raise ValueError("M20 key-state contact/internal balance changed")

    evidence = {
        "schema": "vela.simplemos.sdevice.m20_drain_adjacent_row_audit_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only drain-adjacent electron-continuity and production-Jacobian audit",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "findings": report["findings"],
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "The exact incident SG edge sum and raw electron-continuity terms close across all audited drain-adjacent rows.",
                "At n21, Vd=0.05 V, Vg=0.8 V, the imported-state residual is a contact-versus-internal electron-flux imbalance; local SRH is negligible.",
                "Nodes 991 and 992 amplify the small remaining imbalance through roughly 1275x and 990x state-difference cancellation.",
                "The exact production Jacobian rows are almost entirely electron-quasi-Fermi coupled, with about 44.42 percent contact-column stiffness at the directly connected key nodes.",
            ],
            "not_supported": [
                "Attributing the residual to a wrong drain Dirichlet value, local SRH, or direct psi/phip row coupling.",
                "Using the edge-projection transport-Jacobian probe to reconstruct a transport-cell-vector production row.",
                "Claiming that the audit exposes the private Sentaurus Jacobian or face interpolation.",
                "Changing any production mobility, recombination, or contact default.",
            ],
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m20" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "probe_count": 96}))


if __name__ == "__main__":
    main()
