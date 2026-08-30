#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M19 drain-cut audit evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/drain_cut_audit"
REPORT = ROOT / "m19_drain_cut_audit_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m19_drain_cut_audit_evidence.json")
ARTIFACTS = [
    "m19_drain_cut_audit_report.json",
    "m19_state_summary.csv",
    "m19_drain_cut_edges.csv",
    "m19_frozen_factor_contributions.csv",
    "m19_key_state_drain_cut_edges.csv",
]
FIGURES = [
    "m19_vg0p8_gap_closure.png",
    "m19_key_state_frozen_factors.png",
    "m19_key_state_edge_concentration.png",
    "m19_contact_boundary_error.png",
]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m19_drain_cut_audit_contract_v1.json",
    "scripts/run_simplemos_m19_drain_cut_audit.py",
    "scripts/plot_simplemos_m19_drain_cut_audit.py",
    "scripts/freeze_simplemos_m19_drain_cut_audit_evidence.py",
    "tests/regression/test_simplemos_m19_drain_cut_audit.py",
    "docs/validation/simplemos_m19_drain_cut_audit_2026-08-30.md",
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
    if execution["state_count"] != 16 or execution["sg_probe_count"] != 144:
        raise ValueError("M19 state/probe matrix is incomplete")
    if execution["probe_variants_per_state"] != 9:
        raise ValueError("M19 frozen game/full-state control is incomplete")
    if closure["maximum_contact_phin_bias_error_V"] >= 1.0e-12:
        raise ValueError("M19 drain-contact QF boundary closure failed")
    if (closure[
            "maximum_strong_state_sg_cut_current_reconstruction_relative_error"]
            >= 1.0e-6):
        raise ValueError("M19 strong-state SG cut-current reconstruction failed")
    if closure["maximum_shapley_closure_A_per_um"] >= 1.0e-18:
        raise ValueError("M19 frozen factor Shapley closure failed")
    if execution["new_sentaurus_execution"] or execution["default_model_changed"]:
        raise ValueError("M19 violated the read-only physical-model contract")
    if not findings["contact_boundary_phin_matches_bias_all_states"]:
        raise ValueError("M19 drain contact boundary decision changed")
    if findings["strong_state_max_abs_sentaurus_state_vela_sg_error_dex"] >= 5.0e-4:
        raise ValueError("M19 imported-state SG replay decision changed")
    low, high = findings["vg_0p8_full_state_gap_closed_fraction_range"]
    if low <= 0.97 or high >= 1.01:
        raise ValueError("M19 Vg=0.8 V gap closure decision changed")
    key = findings["key_state"]
    if key["frozen_factor_absolute_shares"][
            "drain_adjacent_interior_phin"] <= 0.999:
        raise ValueError("M19 key-state adjacent-QF attribution changed")

    evidence = {
        "schema": "vela.simplemos.sdevice.m19_drain_cut_audit_evidence.v1",
        "status": "complete_with_bias_regime_separation",
        "scope": "SDevice-only drain boundary, endpoint-QF, SG cut flux, and self-consistent feedback audit",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": execution,
        "closure": closure,
        "findings": findings,
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": {
            "supported": [
                "Vela and imported Sentaurus drain-contact electron quasi-Fermi values equal the configured drain bias in all 16 states.",
                "At Vg=0.8 V, the unchanged Vela drain-cut SG operator replaying the imported Sentaurus state is within 0.000444 dex of Sentaurus terminal current.",
                "The imported state closes 97.54 to 99.43 percent of the four Vg=0.8 V terminal-current gaps.",
                "At n21, Vd=0.05 V, Vg=0.8 V, the drain-adjacent free-node electron quasi-Fermi factor carries more than 99.99 percent of the frozen cut-state correction.",
                "The exact frozen adjacent-QF correction and the M18 self-consistent drain-cut QF response agree in magnitude.",
            ],
            "not_supported": [
                "Claiming every private Sentaurus internal SG face formula is identical to Vela.",
                "Using the guarded sub-fA mapped-state current reconstruction as headline evidence.",
                "Attributing the remaining adjacent-node QF state difference to one model before auditing its complete continuity-row balance.",
                "Changing HFS, Enormal, SRH, reference density, or another production default.",
            ],
        },
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m19" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION},
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 16, "probe_count": 144}))


if __name__ == "__main__":
    main()
