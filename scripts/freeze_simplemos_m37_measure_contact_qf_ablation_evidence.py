#!/usr/bin/env python3
"""Freeze portable M37 measure/contact-SG/QF-feedback evidence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/measure_contact_qf_ablation"
REPORT = ROOT / "m37_measure_contact_qf_ablation_report.json"
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "simplemos_m37_measure_contact_qf_ablation_evidence.json")
ARTIFACTS = [
    "m37_measure_contact_qf_ablation_report.json",
    "m37_boundary_measure_response.csv",
]
SOURCES = [
    "reference_tcad/simplemos_sentaurus2022/simplemos_m37_measure_contact_qf_ablation_contract_v1.json",
    "reference_tcad/simplemos_sentaurus2022/signed_average_box_assembly/m35_signed_average_box_assembly_report.json",
    "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off/m22_n23_hfs_deep_off_report.json",
    "reference_tcad/simplemos_sentaurus2022/n23_hfs_deep_off/m22_key_state_drain_cut_edges.csv",
    "reference_tcad/simplemos_sentaurus2022/first_layer_feedback_audit/m26_first_layer_feedback_audit_report.json",
    "scripts/run_simplemos_m37_measure_contact_qf_ablation.py",
    "scripts/freeze_simplemos_m37_measure_contact_qf_ablation_evidence.py",
    "tests/regression/test_simplemos_m37_measure_contact_qf_ablation.py",
    "tests/test_mos_mixed_material.cpp",
    "tests/test_newton_solver.cpp",
    "docs/validation/simplemos_m37_measure_contact_qf_ablation_2026-08-31.md",
    "include/vela/equation/AssemblerUtils.h",
    "include/vela/equation/DDAssembler.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/solver/NewtonSolver.cpp",
    "CMakeLists.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    report = json.loads(REPORT.read_text(encoding="utf-8"))
    if report["status"] != "complete" or not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError("M37 report is not complete and accepted")
    measure = report["boundary_measure_ablation"]
    contact = report["contact_sg_ablation"]
    feedback = report["quasi_fermi_feedback_ablation"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m37_measure_contact_qf_ablation_evidence.v1",
        "status": "frozen",
        "summary": {
            "boundary_scope_fraction_of_all_signed_log_response": measure[
                "boundary_scope_fraction_of_all_signed_log_response"],
            "contact_sg_direct_electron_delta_A_per_um": contact[
                "signed_direct_electron_current_delta_A_per_um"],
            "qf_adjoint_relaxation_fraction_of_actual": feedback[
                "adjoint_relaxation_fraction_of_actual"],
            "qf_first_order_relative_error": feedback[
                "first_order_relative_error"],
            "default_model_changed": report["execution"]["default_model_changed"],
        },
        "artifacts": [{
            "path": f"reference_tcad/simplemos_sentaurus2022/measure_contact_qf_ablation/{name}",
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
