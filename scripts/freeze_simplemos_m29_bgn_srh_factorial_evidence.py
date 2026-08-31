#!/usr/bin/env python3
"""Freeze portable evidence for the SimpleMOS M29 BGN x SRH audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/bgn_srh_factorial"
REPORT = ROOT / "m29_bgn_srh_factorial_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m29_bgn_srh_factorial_evidence.json")
ARTIFACTS = [
    "m29_bgn_srh_factorial_report.json",
    "m29_factorial_matrix.csv",
    "m29_factor_effects.csv",
    "m29_node_factor_effects.csv",
    "m29_spatial_effect_summary.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m29_bgn_srh_factorial_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m28_exact_bias_path_audit_evidence.json",
    "scripts/run_simplemos_m29_bgn_srh_factorial.py",
    "scripts/plot_simplemos_m29_bgn_srh_factorial.py",
    "scripts/freeze_simplemos_m29_bgn_srh_factorial_evidence.py",
    "tests/regression/test_simplemos_m29_bgn_srh_factorial.py",
    "docs/validation/simplemos_m29_bgn_srh_factorial_2026-08-31.md",
    "docs/validation/figures/simplemos_m29/simplemos_m29_bgn_srh_factorial.png",
    "CMakeLists.txt",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M29 report is not complete and accepted")
    sent = report["factor_effects"]["sentaurus"]
    vela = report["factor_effects"]["vela"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m29_bgn_srh_factorial_evidence.v1",
        "status": "frozen",
        "summary": {
            "new_sentaurus_state_count": report["execution"]["new_sentaurus_state_count"],
            "baseline_absolute_gap_dex": report["matrix"]["baseline_absolute_gap_dex"],
            "best_gap_cell": report["matrix"]["best_gap_cell"],
            "best_absolute_gap_dex": report["matrix"]["best_absolute_gap_dex"],
            "sentaurus_interaction_dex": sent["bgn_srh_interaction_dex"],
            "vela_interaction_dex": vela["bgn_srh_interaction_dex"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/bgn_srh_factorial/{name}",
            "sha256": sha(ROOT / name),
        } for name in ARTIFACTS],
        "source_hashes": {name: sha(REPO / name) for name in SOURCES},
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": str(OUTPUT),
                      **evidence["summary"]}))


if __name__ == "__main__":
    main()
