#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M21 incident-edge evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/incident_edge_balance"
REPORT = ROOT / "m21_incident_edge_balance_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m21_incident_edge_balance_evidence.json")
ARTIFACTS = [
    "m21_incident_edge_balance_report.json",
    "m21_state_summary.csv",
    "m21_edge_input_comparison.csv",
    "m21_node_factor_ledger.csv",
    "m21_key_state_incident_edges.csv",
    "m21_key_state_node_factor_ledger.csv",
]
FIGURES = [
    "m21_key_edge_flux_balance.png",
    "m21_key_edge_input_changes.png",
    "m21_key_factor_shares.png",
]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m21_incident_edge_balance_contract_v1.json",
    "include/vela/equation/CoupledDDAssembler.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "scripts/run_simplemos_m21_incident_edge_balance.py",
    "scripts/plot_simplemos_m21_incident_edge_balance.py",
    "scripts/freeze_simplemos_m21_incident_edge_balance_evidence.py",
    "tests/regression/test_simplemos_m21_incident_edge_balance.py",
    "docs/validation/simplemos_m21_incident_edge_balance_2026-08-30.md",
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
    findings = report["findings"]
    if execution["state_count"] != 16 or execution["sg_probe_count"] != 32:
        raise ValueError("M21 state/probe matrix is incomplete")
    if not execution["all_states_use_transport_cell_vector"]:
        raise ValueError("M21 production drive policy changed")
    if closure["maximum_stable_sg_flux_absolute_closure"] >= 1.0e-12:
        raise ValueError("M21 stable SG closure failed")
    if closure["maximum_edge_factor_absolute_closure"] >= 1.0e-12:
        raise ValueError("M21 edge-factor closure failed")
    if closure["maximum_m20_node_flux_delta_absolute_closure"] >= 1.0e-12:
        raise ValueError("M21-to-M20 node closure failed")
    if closure["maximum_geometry_difference"] != 0.0:
        raise ValueError("M21 frozen-state geometry changed")
    if execution["new_sentaurus_execution"] or execution["default_model_changed"]:
        raise ValueError("M21 violated its read-only physical-model contract")
    if findings["key_factor_absolute_shares"]["qf_log_imbalance"] <= 0.9999:
        raise ValueError("M21 key-state QF decision changed")
    if findings["all_state_qf_absolute_share_range"][0] <= 0.997:
        raise ValueError("M21 cross-state QF dominance changed")

    evidence = {
        "schema": "vela.simplemos.sdevice.m21_incident_edge_balance_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only production-SG input and incident-edge factor audit",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "findings": findings,
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "All audited states use the production transport_cell_vector mobility drive.",
                "The production stable SG decomposition exactly reproduces every selected edge flux.",
                "At the key state, 99.9936 percent of the absolute incident-edge response is assigned to QF log imbalance.",
                "The leading contact/internal pairs have 12.6 to 13.9 percent QF-drop changes but only ppm-scale mobility, density, and Bernoulli changes.",
                "Finite-volume length and couple are identical between paired frozen states.",
            ],
            "not_supported": [
                "Attributing the key first-layer residual to HFS mobility response, Bernoulli weights, carrier density, or geometry.",
                "Claiming that QF state attribution by itself identifies the upstream coupled-solve cause.",
                "Treating Vela intermediate quantities as private Sentaurus face outputs.",
                "Changing any production mobility, contact, or solver default.",
            ],
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m21" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "probe_count": 32}))


if __name__ == "__main__":
    main()
