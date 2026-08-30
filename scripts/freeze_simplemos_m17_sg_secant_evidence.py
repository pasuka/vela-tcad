#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M17 SG secant-factor evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/sg_secant"
REPORT = ROOT / "m17_sg_secant_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m17_sg_secant_evidence.json")
ARTIFACTS = [
    "m17_sg_secant_report.json",
    "m17_state_summary.csv",
    "m17_factor_contributions.csv",
    "m17_region_factor_contributions.csv",
    "m17_key_state_node_contributions.csv",
]
FIGURES = [
    "m17_factor_response_n21_vg0p8.png",
    "m17_factor_fraction_16_states.png",
    "m17_cancellation_reduction_n21_vg0p8.png",
    "m17_region_response_n21_vg0p8.png",
    "m17_nodal_qf_imbalance_n21_vg0p8.png",
]
IMPLEMENTATION = [
    "include/vela/equation/CoupledDDAssembler.h",
    "include/vela/solver/NewtonSolver.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m17_sg_secant_contract_v1.json",
    "scripts/run_simplemos_m17_sg_secant.py",
    "scripts/plot_simplemos_m17_sg_secant.py",
    "scripts/freeze_simplemos_m17_sg_secant_evidence.py",
    "tests/regression/test_simplemos_m17_sg_secant.py",
    "docs/validation/simplemos_m17_sg_secant_2026-08-29.md",
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
    if execution["state_count"] != 16 or execution["factor_probe_count"] != 16:
        raise ValueError("M17 state/probe matrix is incomplete")
    if execution["probe_variant_count"] != 8:
        raise ValueError("M17 three-factor vertices are incomplete")
    if closure["maximum_absolute_A_per_um"] >= 1.0e-12:
        raise ValueError("M17 absolute closure failed")
    if closure["maximum_relative"] >= 1.0e-4:
        raise ValueError("M17 all-state relative closure failed")
    if closure["maximum_strong_state_relative"] >= 1.0e-8:
        raise ValueError("M17 strong-state relative closure failed")
    if closure["maximum_endpoint_absolute_A_per_um"] != 0.0:
        raise ValueError("M17 production endpoints did not close exactly")
    if execution["new_sentaurus_execution"] or execution["default_model_changed"]:
        raise ValueError("M17 violated the read-only physical-model contract")

    evidence = {
        "schema": "vela.simplemos.sdevice.m17_sg_secant_evidence.v1",
        "status": "complete_with_bias_regime_separation",
        "scope": "SDevice-only stable SG secant-factor attribution of the M15 electron transport response",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "key_state": report["key_state"],
        "strong_state_fraction_range": report["strong_state_fraction_range"],
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "The stable mobility-times-secant-conductance-times-QF-imbalance factorization closes all 16 production endpoint responses.",
                "At n21, Vd=0.05 V, Vg=0.8 V, QF log imbalance supplies 92.89 percent of absolute factor response.",
                "At the key state, mobility supplies 2.10 percent and SG secant conductance supplies 5.01 percent of absolute factor response.",
                "The key-state cancellation amplification falls from 461.1x in the ungrouped M16 split to 1.166x in M17.",
                "At Vg=2.5 V, SG secant conductance dominates all four paired states with 75.16 to 93.04 percent of absolute response.",
                "Standalone Enormal diagnostics initialize the same interface geometry as the production residual, independent of diagnostic call order.",
                "The n21 low-drain weak-current discrepancy should next be localized through the channel quasi-Fermi imbalance and continuity state."
            ],
            "not_supported": [
                "Claiming a single Sentaurus-internal formula causes the QF-imbalance response.",
                "Relabeling the coupled SG secant conductance as a single electrostatic or carrier-density model.",
                "Claiming HFS has zero effect or changing a production model default.",
                "Using sub-fA relative responses as headline root-cause evidence.",
                "Treating pre-fix call-order-dependent M17 factor allocations as frozen evidence."
            ]
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m17" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "variant_count": 8}))


if __name__ == "__main__":
    main()
