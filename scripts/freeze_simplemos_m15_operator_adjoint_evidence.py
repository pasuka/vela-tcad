#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M15 operator-adjoint evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/operator_adjoint"
REPORT = ROOT / "m15_operator_adjoint_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m15_operator_adjoint_evidence.json")
ARTIFACTS = [
    "m15_operator_adjoint_report.json",
    "m15_state_summary.csv",
    "m15_component_contributions.csv",
    "m15_region_contributions.csv",
    "m15_key_state_node_contributions.csv",
    "m15_m11_mobility_crosscheck.csv",
]
FIGURES = [
    "m15_component_response_n21_vg0p8.png",
    "m15_region_response_n21_vg0p8.png",
    "m15_component_fraction_16_states.png",
    "m15_m11_mobility_crosscheck_n21_vg0p8.png",
    "m15_nodal_response_n21_vg0p8.png",
]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m15_operator_adjoint_contract_v1.json",
    "scripts/run_simplemos_m15_operator_adjoint.py",
    "scripts/plot_simplemos_m15_operator_adjoint.py",
    "scripts/freeze_simplemos_m15_operator_adjoint_evidence.py",
    "tests/regression/test_simplemos_m15_operator_adjoint.py",
    "docs/validation/simplemos_m15_operator_adjoint_2026-08-29.md",
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
    if execution["state_count"] != 16:
        raise ValueError("M15 state matrix is incomplete")
    if execution["carrier_term_probe_count"] != 32:
        raise ValueError("M15 carrier-term probe matrix is incomplete")
    if closure["maximum_absolute_A_per_um"] >= 1.0e-12:
        raise ValueError("M15 absolute operator-adjoint closure failed")
    if closure["maximum_relative"] >= 1.0e-8:
        raise ValueError("M15 relative operator-adjoint closure failed")
    if execution["new_sentaurus_execution"] or execution["cpp_changed"]:
        raise ValueError("M15 violated the read-only execution contract")
    if execution["default_model_changed"]:
        raise ValueError("M15 changed a default physical model")

    evidence = {
        "schema": "vela.simplemos.sdevice.m15_operator_adjoint_evidence.v1",
        "status": "complete_with_declared_operator_boundary",
        "scope": "SDevice-only Vela operator-residual attribution using the M14 terminal-current adjoint",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "key_state": report["key_state"],
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "All 16 paired states and 32 solved-equation carrier-term probes are complete.",
                "The equation-term sum reconstructs the M14 full-field adjoint relaxation within the frozen absolute and relative closure ceilings.",
                "At n21, Vd=0.05 V, Vg=0.8 V, electron transport supplies 95.66 percent of absolute equation-component response and Poisson supplies 4.34 percent.",
                "At the key state, channel response is +13.5254 percent of baseline current and drain response is -0.250314 percent; SRH is negligible.",
                "Electron-transport dominance identifies a broad Vela transport-operator mismatch, not HFS uniquely."
            ],
            "not_supported": [
                "Interpreting Vela residual terms at a mapped Sentaurus state as native Sentaurus equation residuals.",
                "Claiming HFS, PhuMob, electrostatics, or SRH as the unique cross-code root cause.",
                "Using sub-fA relative contributions as headline quantitative attribution.",
                "Changing a default physical model from this read-only diagnostic."
            ]
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m15" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "carrier_term_probe_count": 32}))


if __name__ == "__main__":
    main()
