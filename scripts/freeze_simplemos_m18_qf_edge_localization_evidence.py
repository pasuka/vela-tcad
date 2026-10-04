#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M18 QF edge-localization evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/qf_edge_localization"
REPORT = ROOT / "m18_qf_edge_localization_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m18_qf_edge_localization_evidence.json")
ARTIFACTS = [
    "m18_qf_edge_localization_report.json",
    "m18_state_summary.csv",
    "m18_bucket_factor_contributions.csv",
    "m18_key_state_edge_contributions.csv",
    "m18_feedback_ledger.csv",
]
FIGURES = [
    "m18_qf_bucket_support_n21_vg0p8.png",
    "m18_qf_edge_map_n21_vg0p8.png",
    "m18_qf_edge_pareto_n21_vg0p8.png",
    "m18_feedback_ledger_n21_vg0p8.png",
    "m18_qf_contact_fraction_16_states.png",
]
IMPLEMENTATION = [
    "include/vela/equation/CoupledDDAssembler.h",
    "include/vela/solver/NewtonSolver.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/solver/NewtonSolver.cpp",
    "src/tools/vela_example_runner.cpp",
    "tests/test_newton_solver.cpp",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m18_qf_edge_localization_contract_v1.json",
    "scripts/run_simplemos_m18_qf_edge_localization.py",
    "scripts/plot_simplemos_m18_qf_edge_localization.py",
    "scripts/freeze_simplemos_m18_qf_edge_localization_evidence.py",
    "tests/regression/test_simplemos_m18_qf_edge_localization.py",
    "docs/validation/simplemos_m18_qf_edge_localization_2026-08-30.md",
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
    if execution["state_count"] != 16 or execution["edge_probe_count"] != 16:
        raise ValueError("M18 state/probe matrix is incomplete")
    if execution["probe_variant_count"] != 8:
        raise ValueError("M18 three-factor vertices are incomplete")
    if closure["maximum_edge_to_node_absolute_A_per_um"] >= 1.0e-12:
        raise ValueError("M18 edge-to-node absolute closure failed")
    if closure["maximum_edge_to_node_relative"] >= 1.0e-8:
        raise ValueError("M18 edge-to-node relative closure failed")
    if closure["maximum_feedback_ledger_absolute_A_per_um"] >= 1.0e-12:
        raise ValueError("M18 feedback ledger closure failed")
    if execution["new_sentaurus_execution"] or execution["default_model_changed"]:
        raise ValueError("M18 violated the read-only physical-model contract")

    key = report["key_state"]
    contact_fraction = key["combined_contact_cut_absolute_edge_support_fraction"]
    channel_fraction = key["internal_channel_absolute_edge_support_fraction"]
    if not (contact_fraction > 0.80 and channel_fraction < 0.20):
        raise ValueError("M18 key-state localization decision changed")
    evidence = {
        "schema": "vela.simplemos.sdevice.m18_qf_edge_localization_evidence.v1",
        "status": "complete_with_bias_regime_separation",
        "scope": "SDevice-only native-edge localization of M17 QF imbalance and complete M15 feedback ledger",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "key_state": key,
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "All 16 native-edge Shapley games close to the regenerated M17 node-factor results.",
                "At n21, Vd=0.05 V, Vg=0.8 V, contact cuts carry 80.6363 percent of absolute QF response and the internal channel carries 19.3637 percent.",
                "The drain cut supplies effectively all material contact-cut response; one drain-cut edge supplies 66.5945 percent and the leading two supply 79.4133 percent of absolute QF support.",
                "The complete key-state feedback ledger is dominated by QF log imbalance, with smaller Poisson feedback and opposing mobility and SG conductance responses.",
                "Direct SRH operator response is negligible at the key state.",
                "The diagnostic Enormal geometry path is now independent of whether residual evaluation occurred first."
            ],
            "not_supported": [
                "Claiming SRH had no historical effect on the converged state.",
                "Equating the localized Vela edge response with a native Sentaurus residual decomposition.",
                "Claiming a private Sentaurus formula is the unique root cause.",
                "Changing production HFS, Enormal, SRH, or other model defaults.",
                "Using sub-fA off-state relative localization as headline evidence."
            ]
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m18" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "variant_count": 8}))


if __name__ == "__main__":
    main()
