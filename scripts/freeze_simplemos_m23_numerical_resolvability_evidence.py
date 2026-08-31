#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M23 numerical-resolvability evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/numerical_resolvability"
REPORT = ROOT / "m23_numerical_resolvability_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m23_numerical_resolvability_evidence.json")
ARTIFACTS = ["m23_numerical_resolvability_report.json",
             "m23_edge_resolvability.csv", "m23_state_summary.csv",
             "m23_solver_path_reproducibility.csv"]
UPSTREAM = [
    "reference_tcad/simplemos_sentaurus2022/incident_edge_balance/m21_key_state_incident_edges.csv",
    "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off/m22_key_state_self_consistent_edges.csv",
]
FIGURE = "docs/validation/figures/simplemos_m23/simplemos_m23_numerical_resolvability.png"
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m23_numerical_resolvability_contract_v1.json",
    "scripts/run_simplemos_m23_numerical_resolvability.py",
    "scripts/plot_simplemos_m23_numerical_resolvability.py",
    "scripts/freeze_simplemos_m23_numerical_resolvability_evidence.py",
    "tests/regression/test_simplemos_m23_numerical_resolvability.py",
    "docs/validation/simplemos_m23_numerical_resolvability_2026-08-30.md",
    "include/vela/discretization/ScharfetterGummel.h",
    "include/vela/equation/CoupledDDAssembler.h",
    "src/discretization/ScharfetterGummel.cpp",
    "src/equation/CoupledDDAssembler.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(path: Path) -> dict[str, str]:
    return {"path": path.resolve().relative_to(REPO.resolve()).as_posix(),
            "sha256": sha256(path)}


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    findings = report["findings"]
    if report["execution"]["edge_evaluation_count"] != 36:
        raise ValueError("M23 edge audit is incomplete")
    if not findings["all_edges_resolved_in_production_representation"]:
        raise ValueError("M23 found an unresolved production SG edge")
    if findings["maximum_production_vs_decimal_relative_error"] >= 1.0e-12:
        raise ValueError("M23 production/Decimal closure regressed")
    if findings["solver_path_spread"]["full"]["spread_dex"] >= 1.0e-4:
        raise ValueError("M23 HFS-on solver-path spread regressed")

    evidence = {
        "schema": "vela.simplemos.sdevice.m23_numerical_resolvability_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only SG/QF numerical resolvability at M21/M22 key states",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": report["execution"],
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
                      "edge_evaluations": 36, "solver_controls": 8}))


if __name__ == "__main__":
    main()
