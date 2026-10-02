#!/usr/bin/env python3
"""Freeze reproducible SimpleMOS M22 n23 deep-off HFS evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off"
REPORT = ROOT / "m22_n23_hfs_deep_off_report.json"
EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m22_n23_hfs_deep_off_evidence.json")
ARTIFACTS = [
    "m22_n23_hfs_deep_off_report.json",
    "m22_summary.csv",
    "m22_drain_cut_frozen_edges.csv",
    "m22_key_state_drain_cut_edges.csv",
    "m22_drain_cut_self_consistent_edges.csv",
    "m22_key_state_self_consistent_edges.csv",
]
FIGURES = ["simplemos_m22_n23_hfs_deep_off.png"]
IMPLEMENTATION = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m22_n23_hfs_deep_off_contract_v1.json",
    "scripts/run_simplemos_m22_n23_hfs_deep_off.py",
    "scripts/plot_simplemos_m22_n23_hfs_deep_off.py",
    "scripts/freeze_simplemos_m22_n23_hfs_deep_off_evidence.py",
    "tests/regression/test_simplemos_m22_n23_hfs_deep_off.py",
    "docs/validation/simplemos_m22_n23_hfs_deep_off_2026-08-30.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def entry(path: Path) -> dict[str, str]:
    return {"path": path.resolve().relative_to(REPO.resolve()).as_posix(),
            "sha256": sha256(path)}


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    key = report["key_state"]
    findings = report["findings"]
    if report["state_count"] != 8 or report["frozen_formula_evaluations"] != 16:
        raise ValueError("M22 state/formula matrix is incomplete")
    if key["rerun_full_current_A_per_um"] != key["vela_full_current_A_per_um"]:
        raise ValueError("M22 full-current baseline did not reproduce exactly")
    if key["rerun_no_hfs_current_A_per_um"] != key["vela_no_hfs_current_A_per_um"]:
        raise ValueError("M22 no-HFS baseline did not reproduce exactly")
    if abs(findings["key_full_state_frozen_hfs_effect_dex"]) >= 1.0e-5:
        raise ValueError("M22 frozen direct HFS response changed")
    if findings["key_max_per_edge_sg_cancellation_condition"] <= 1.0e14:
        raise ValueError("M22 expected deep-off SG conditioning was not observed")
    if findings["key_max_frozen_drain_cut_mobility_effect_dex"] != 0.0:
        raise ValueError("M22 drain-cut direct mobility response changed")

    evidence = {
        "schema": "vela.simplemos.sdevice.m22_n23_hfs_deep_off_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only n23 deep-off HFS state-feedback attribution",
        "sentaurus_release": "T-2022.03-SP2",
        "execution": {
            "device": "n23",
            "drain_voltage_V": 0.05,
            "gate_voltages_V": [0.0, 0.05, 0.1, 0.15],
            "state_count": report["state_count"],
            "frozen_formula_evaluations": report[
                "frozen_formula_evaluations"],
            "new_sentaurus_execution": False,
            "default_model_changed": False
        },
        "key_state": key,
        "findings": findings,
        "interpretation_limits": report["interpretation_limits"],
        "conclusions": report["conclusions"],
        "artifacts": {name: entry(ROOT / name) for name in ARTIFACTS},
        "figures": [entry(REPO / "docs/validation/figures/simplemos_m22" / name)
                    for name in FIGURES],
        "implementation_sha256": {name: sha256(REPO / name)
                                  for name in IMPLEMENTATION}
    }
    EVIDENCE.write_text(json.dumps(evidence, indent=2) + "\n",
                        encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": entry(EVIDENCE)["path"],
                      "state_count": 8, "formula_evaluations": 16}))


if __name__ == "__main__":
    main()
