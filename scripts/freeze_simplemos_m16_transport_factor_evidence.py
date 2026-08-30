#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M16 transport-factor evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/transport_factor"
REPORT = ROOT / "m16_transport_factor_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m16_transport_factor_evidence.json")
ARTIFACTS = [
    "m16_transport_factor_report.json",
    "m16_state_summary.csv",
    "m16_factor_contributions.csv",
    "m16_region_factor_contributions.csv",
    "m16_grouped_factor_contributions.csv",
    "m16_grouped_region_contributions.csv",
    "m16_key_state_node_contributions.csv",
]
FIGURES = [
    "m16_grouped_factor_response_n21_vg0p8.png",
    "m16_four_factor_cancellation_n21_vg0p8.png",
    "m16_grouped_factor_fraction_16_states.png",
    "m16_region_response_n21_vg0p8.png",
    "m16_nodal_sg_kernel_response_n21_vg0p8.png",
]
IMPLEMENTATION = [
    "include/vela/equation/CoupledDDAssembler.h",
    "include/vela/solver/NewtonSolver.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m16_transport_factor_contract_v1.json",
    "scripts/run_simplemos_m16_transport_factor.py",
    "scripts/plot_simplemos_m16_transport_factor.py",
    "scripts/freeze_simplemos_m16_transport_factor_evidence.py",
    "tests/regression/test_simplemos_m16_transport_factor.py",
    "docs/validation/simplemos_m16_transport_factor_2026-08-29.md",
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
        raise ValueError("M16 state/probe matrix is incomplete")
    if execution["probe_variant_count"] != 16:
        raise ValueError("M16 four-factor vertices are incomplete")
    if closure["maximum_absolute_A_per_um"] >= 1.0e-12:
        raise ValueError("M16 absolute closure failed")
    if closure["maximum_relative"] >= 1.0e-4:
        raise ValueError("M16 all-state relative closure failed")
    if closure["maximum_strong_state_relative"] >= 1.0e-8:
        raise ValueError("M16 strong-state relative closure failed")
    if closure["maximum_endpoint_absolute_A_per_um"] >= 1.0e-15:
        raise ValueError("M16 production endpoints exceeded roundoff closure")
    if execution["new_sentaurus_execution"] or execution["default_model_changed"]:
        raise ValueError("M16 violated the read-only physical-model contract")

    evidence = {
        "schema": "vela.simplemos.sdevice.m16_transport_factor_evidence.v1",
        "status": "complete_with_grouped_kernel_interpretation",
        "scope": "SDevice-only Shapley decomposition of the M15 electron SG transport response",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "key_state": report["key_state"],
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "All 16 states and all 16 four-factor vertices per state are complete.",
                "Production endpoint terms close exactly and the adjoint-weighted factor sum closes to M15 below 1e-12 A/um.",
                "At n21, Vd=0.05 V, Vg=0.8 V, the grouped SG state kernel contributes +13.0101 percent of baseline current and 97.67 percent of absolute grouped response.",
                "At the key state, mobility state contributes -0.0521 percent and mobility drive contributes -0.2588 percent of baseline current.",
                "Across Vg at or above 0.8 V, the SG state kernel supplies 80.50 to 98.68 percent of absolute grouped response.",
                "Independent Bernoulli and population interventions are cancellation-dominated and must be interpreted as a coupled kernel."
            ],
            "not_supported": [
                "Ranking standalone Bernoulli and carrier-population Shapley values as independent physical error sources.",
                "Claiming HFS has no interaction through the mobility-state factor.",
                "Treating Vela factor terms as native Sentaurus residual terms.",
                "Changing a default physical model from this diagnostic.",
                "Treating pre-fix call-order-dependent Enormal allocations as frozen evidence."
            ]
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m16" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "variant_count": 16}))


if __name__ == "__main__":
    main()
