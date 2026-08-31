#!/usr/bin/env python3
"""Freeze portable evidence for the SimpleMOS M32 scaled-QF audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/electron_qf_scaled_perturbation"
REPORT = ROOT / "m32_electron_qf_scaled_perturbation_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m32_electron_qf_scaled_perturbation_evidence.json")
ARTIFACTS = [
    "m32_electron_qf_scaled_perturbation_report.json",
    "m32_fraction_response.csv",
    "m32_cross_block_scaling.csv",
    "m32_drain_cut_edge_response.csv",
    "m32_cross_block_node_localization.csv",
    "m32_top_transport_edges.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m32_electron_qf_scaled_perturbation_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/simplemos_m31_minority_poisson_perturbation_evidence.json",
    "scripts/run_simplemos_m32_electron_qf_scaled_perturbation.py",
    "scripts/plot_simplemos_m32_electron_qf_scaled_perturbation.py",
    "scripts/freeze_simplemos_m32_electron_qf_scaled_perturbation_evidence.py",
    "tests/regression/test_simplemos_m32_electron_qf_scaled_perturbation.py",
    "docs/validation/simplemos_m32_electron_qf_scaled_perturbation_2026-08-31.md",
    "docs/validation/figures/simplemos_m32/simplemos_m32_electron_qf_scaled_perturbation.png",
    "CMakeLists.txt",
]


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M32 report is not complete and accepted")
    full = next(row for row in report["fraction_response"]
                if row["fraction"] == 1.0)
    evidence = {
        "schema": "vela.simplemos.sdevice.m32_electron_qf_scaled_perturbation_evidence.v1",
        "status": "frozen",
        "summary": {
            "fraction_count": report["acceptance"]["fraction_count"],
            "probe_count": (report["acceptance"]["functional_probe_count"]
                            + report["acceptance"]["sg_probe_count"]
                            + report["acceptance"]["cross_block_probe_count"]),
            "largest_tested_valid_fraction": report["linearization"][
                "largest_tested_valid_fraction"],
            "full_fraction_current_shift_dex": full[
                "absolute_log10_current_shift_dex"],
            "full_fraction_drain_cut_condition": full["drain_cut_condition"],
            "maximum_schur_relative_closure": report["acceptance"][
                "maximum_schur_relative_closure"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/electron_qf_scaled_perturbation/{name}",
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
