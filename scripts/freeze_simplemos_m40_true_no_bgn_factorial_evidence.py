#!/usr/bin/env python3
"""Freeze portable M40 explicit no-BGN factorial evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/true_no_bgn_factorial"
REPORT = ROOT / "m40_true_no_bgn_factorial_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m40_true_no_bgn_factorial_evidence.json")
ARTIFACTS = [
    "m40_true_no_bgn_factorial_report.json",
    "m40_factorial_matrix.csv",
    "m40_factor_effects.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m40_true_no_bgn_factorial_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m39_bgn_chain_first_divergence_evidence.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m29_bgn_srh_factorial_evidence.json",
    "scripts/run_simplemos_m40_true_no_bgn_factorial.py",
    "scripts/freeze_simplemos_m40_true_no_bgn_factorial_evidence.py",
    "tests/regression/test_simplemos_m40_true_no_bgn_factorial.py",
    "docs/validation/simplemos_m40_true_no_bgn_factorial_2026-09-01.md",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M40 report is not complete and accepted")
    result = report["causal_result"]
    if not 0.87 < result["bgn_related_fraction_of_production_gap"] < 0.89:
        raise RuntimeError("M40 conditional BGN attribution changed")
    evidence = {
        "schema": "vela.simplemos.sdevice.m40_true_no_bgn_factorial_evidence.v1",
        "status": "frozen",
        "summary": {
            "manual_keyword": report["manual_confirmation"]["confirmed_keyword"],
            "manual_sha256": report["manual_confirmation"]["manual_sha256"],
            "production_gap_dex": result["production_gap_dex"],
            "true_no_bgn_srh_on_gap_dex": result["true_no_bgn_srh_on_gap_dex"],
            "bgn_response_mismatch_srh_on_dex": result[
                "bgn_response_mismatch_srh_on_dex"],
            "bgn_related_fraction_of_production_gap": result[
                "bgn_related_fraction_of_production_gap"],
            "factorial_main_effect_qualified": result[
                "factorial_main_effect_qualified"],
            "sentaurus_state_count": report["execution"]["new_sentaurus_state_count"],
            "vela_state_count": report["execution"]["new_vela_state_count"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/true_no_bgn_factorial/{name}",
            "sha256": sha256(ROOT / name),
        } for name in ARTIFACTS],
        "source_hashes": {name: sha256(REPO / name) for name in SOURCES},
    }
    OUTPUT.write_text(json.dumps(evidence, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": "frozen", "output": str(OUTPUT),
                      **evidence["summary"]}))


if __name__ == "__main__":
    main()
