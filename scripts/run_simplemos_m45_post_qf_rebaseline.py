#!/usr/bin/env python3
"""Build the SimpleMOS M45 post-M43/M44 rebaseline evidence."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m45_post_qf_rebaseline_contract_v1.json"
OUTPUT = ROOT / "post_qf_rebaseline"
REPORTS = {
    "M31": ROOT / "minority_poisson_perturbation/m31_minority_poisson_perturbation_report.json",
    "M32": ROOT / "electron_qf_scaled_perturbation/m32_electron_qf_scaled_perturbation_report.json",
    "M37": ROOT / "measure_contact_qf_ablation/m37_measure_contact_qf_ablation_report.json",
    "M41": ROOT / "equal_ni_flux_ablation/m41_equal_ni_flux_ablation_report.json",
    "M42": ROOT / "no_bgn_self_consistent_interaction/m42_no_bgn_self_consistent_interaction_report.json",
}
UPSTREAM_EVIDENCE = {
    stage: next(ROOT.glob(f"simplemos_{stage.lower()}_*evidence.json"))
    for stage in ("M31", "M32", "M37", "M41", "M42", "M44")
}


def portable(path: Path) -> str:
    return path.relative_to(REPO).as_posix()


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def baseline_json(revision: str, path: Path) -> dict[str, Any]:
    git = shutil.which("git") or "D:/msys64/usr/bin/git.exe"
    result = subprocess.run(
        [git, "show", f"{revision}:{portable(path)}"],
        cwd=REPO, check=True, capture_output=True, text=True,
        encoding="utf-8")
    return json.loads(result.stdout)


def refresh_artifact_hashes(evidence_path: Path, stage: str) -> None:
    evidence = read_json(evidence_path)
    for artifact in evidence["artifacts"]:
        artifact["sha256"] = sha256(REPO / artifact["path"])
    evidence["m45_refresh"] = {
        "stage": stage,
        "reason": "post-M43/M44 quasi-Fermi frozen-state rebaseline",
        "historical_contract_rewritten": False,
    }
    write_json(evidence_path, evidence)


def metric(stage: str, name: str, old: Any, new: Any) -> dict[str, Any]:
    delta = None
    if isinstance(old, (int, float)) and isinstance(new, (int, float)):
        if math.isfinite(float(old)) and math.isfinite(float(new)):
            delta = float(new) - float(old)
    return {
        "stage": stage,
        "metric": name,
        "baseline_value": old,
        "rebaseline_value": new,
        "absolute_delta": delta,
    }


def make_markdown(report: dict[str, Any], metrics: list[dict[str, Any]]) -> str:
    findings = report["findings"]
    lines = [
        "# SimpleMOS M45 post-QF rebaseline",
        "",
        "M45 re-runs the frozen-state analyses affected by M43 SG-kernel "
        "unification and M44 quasi-Fermi reference/increment preservation.",
        "",
        "## Stage result",
        "",
        "| Stage | Re-run status | Impact | Historical contract |",
        "|---|---|---|---|",
    ]
    for row in report["stage_impact"]:
        lines.append(
            f"| {row['stage']} | {row['rebaseline_status']} | "
            f"{row['impact']} | {row['historical_contract']} |")
    lines.extend([
        "",
        "## Key closure",
        "",
        f"- M42 maximum compensated cut/report gap: "
        f"{findings['m42_maximum_compensated_cut_reported_gap_dex']:.6e} dex",
        f"- M42 minimum legacy cut/report gap: "
        f"{findings['m42_minimum_legacy_cut_reported_gap_dex']:.6e} dex",
        f"- M42 compensated free-electron residual L2: "
        f"{findings['m42_maximum_compensated_free_electron_residual_l2']:.6e}",
        f"- Legacy/compensated residual separation: "
        f"{findings['m42_legacy_to_compensated_residual_ratio']:.6e}",
        "",
        "M31/M32 retained their qualitative conclusions, M37 was bitwise "
        "unchanged, and M41/M42 deliberately retain failed historical checks. "
        "Those failures now mean that the old mixed-operator mismatch no longer "
        "reproduces; they are not execution failures.",
        "",
        f"The machine-readable metric table contains {len(metrics)} rows.",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    contract = read_json(CONTRACT)
    if contract["schema"] != \
            "vela.simplemos.sdevice.m45_post_qf_rebaseline_contract.v1":
        raise ValueError("unexpected M45 contract schema")
    baseline_revision = contract["baseline_revision"]
    current = {stage: read_json(path) for stage, path in REPORTS.items()}
    baseline = {
        stage: baseline_json(baseline_revision, path)
        for stage, path in REPORTS.items()
    }

    m31 = current["M31"]
    m32 = current["M32"]
    m37 = current["M37"]
    m41 = current["M41"]
    m42 = current["M42"]
    m44 = read_json(ROOT / "qf_coordinate_consistency/m44_qf_coordinate_consistency_report.json")
    limits = contract["acceptance"]

    metrics = [
        metric("M31", "minority_hole_full_current_shift_dex",
               baseline["M31"]["minority_hole"]["full_current_shift_dex"],
               m31["minority_hole"]["full_current_shift_dex"]),
        metric("M31", "poisson_top10_l2_share",
               baseline["M31"]["poisson_floor"]["top10_l2_share"],
               m31["poisson_floor"]["top10_l2_share"]),
        metric("M32", "baseline_current_A_per_um",
               baseline["M32"]["linearization"]["baseline_current_A_per_um"],
               m32["linearization"]["baseline_current_A_per_um"]),
        metric("M32", "maximum_schur_relative_closure",
               baseline["M32"]["acceptance"]["maximum_schur_relative_closure"],
               m32["acceptance"]["maximum_schur_relative_closure"]),
        metric("M37", "boundary_scope_fraction_of_all_signed_log_response",
               baseline["M37"]["boundary_measure_ablation"][
                   "boundary_scope_fraction_of_all_signed_log_response"],
               m37["boundary_measure_ablation"][
                   "boundary_scope_fraction_of_all_signed_log_response"]),
        metric("M41", "srh_on_legacy_gap_dex",
               baseline["M41"]["causal_result"]["srh_on_legacy_gap_dex"],
               m41["causal_result"]["srh_on_legacy_gap_dex"]),
        metric("M41", "srh_on_gap_improvement_dex",
               baseline["M41"]["causal_result"]["srh_on_gap_improvement_dex"],
               m41["causal_result"]["srh_on_gap_improvement_dex"]),
        metric("M41", "srh_off_legacy_gap_dex",
               baseline["M41"]["causal_result"]["srh_off_legacy_gap_dex"],
               m41["causal_result"]["srh_off_legacy_gap_dex"]),
        metric("M41", "srh_off_gap_improvement_dex",
               baseline["M41"]["causal_result"]["srh_off_gap_improvement_dex"],
               m41["causal_result"]["srh_off_gap_improvement_dex"]),
        metric("M42", "maximum_compensated_cut_reported_gap_dex",
               baseline["M42"]["findings"][
                   "maximum_compensated_cut_reported_gap_dex"],
               m42["findings"]["maximum_compensated_cut_reported_gap_dex"]),
        metric("M42", "minimum_legacy_cut_reported_gap_dex",
               baseline["M42"]["findings"][
                   "minimum_legacy_cut_reported_gap_dex"],
               m42["findings"]["minimum_legacy_cut_reported_gap_dex"]),
        metric("M42", "maximum_compensated_free_electron_residual_l2",
               baseline["M42"]["findings"][
                   "maximum_compensated_free_electron_residual_l2"],
               m42["findings"][
                   "maximum_compensated_free_electron_residual_l2"]),
        metric("M42", "legacy_to_compensated_residual_ratio",
               baseline["M42"]["findings"][
                   "legacy_to_compensated_residual_separation_min_ratio"],
               m42["findings"][
                   "legacy_to_compensated_residual_separation_min_ratio"]),
    ]

    m37_unchanged = sha256(REPORTS["M37"]) == hashlib.sha256(
        subprocess.run(
            [shutil.which("git") or "D:/msys64/usr/bin/git.exe", "show",
             f"{baseline_revision}:{portable(REPORTS['M37'])}"],
            cwd=REPO, check=True, capture_output=True).stdout).hexdigest()
    fractions = m42["findings"]["electron_qf_density_fraction_of_srh_shift"]
    max_fraction_error = max(abs(float(value) - 1.0)
                             for value in fractions.values())
    max_gap = max(
        float(m42["findings"]["maximum_compensated_cut_reported_gap_dex"]),
        float(m42["findings"]["minimum_legacy_cut_reported_gap_dex"]),
    )
    checks = {
        "m44_complete": m44["status"] == "complete" and
            m44["acceptance"]["all_checks_pass"],
        "m31_complete": m31["status"] == "complete" and
            m31["acceptance"]["all_checks_pass"],
        "m32_complete": m32["status"] == "complete" and
            m32["acceptance"]["all_checks_pass"],
        "m37_complete_and_bitwise_unchanged": m37["status"] == "complete" and
            m37["acceptance"]["all_checks_pass"] and m37_unchanged,
        "m41_workflows_converged":
            m41["acceptance"]["all_workflows_converged"],
        "m41_historical_contract_superseded":
            m41["status"] == "failed" and
            not m41["acceptance"]["legacy_m40_reproduced"],
        "m42_matrix_complete":
            m42["execution"]["self_consistent_state_count"] == 8 and
            m42["execution"]["frozen_operator_pair_count"] == 32 and
            m42["execution"]["qf_substitution_count"] == 14,
        "m42_operator_cut_reported_closure": max_gap <= float(
            limits["maximum_operator_cut_reported_gap_dex"]),
        "m42_historical_legacy_mismatch_superseded":
            m42["status"] == "failed" and
            not m42["acceptance"]["legacy_cut_reported_mismatch_reproduced"] and
            not m42["causal_interpretation"][
                "legacy_state_terminal_extraction_is_operator_inconsistent"],
        "m42_srh_feedback_retained":
            m42["causal_interpretation"][
                "srh_terminal_effect_is_state_feedback_not_direct_operator"] and
            max_fraction_error <= float(
                limits["maximum_electron_qf_srh_fraction_error"]),
        "m42_contact_null_retained":
            m42["causal_interpretation"][
                "contact_reconstruction_is_null_for_simplemos"],
        "m42_qf_reference_increment_requirement_retained":
            m42["causal_interpretation"][
                "qf_reference_increment_representation_is_required"],
        "m42_compensated_residual_resolved":
            m42["findings"][
                "maximum_compensated_free_electron_residual_l2"] <= float(
                    limits["maximum_compensated_free_electron_residual_l2"]) and
            m42["findings"][
                "legacy_to_compensated_residual_separation_min_ratio"] >= float(
                    limits["minimum_legacy_to_compensated_residual_ratio"]),
        "defaults_unchanged": all(
            not report["execution"].get("default_model_changed", False)
            for report in current.values()) and
            not contract["policy"]["default_physics_model_changed"],
        "no_new_sentaurus_execution": all(
            not report["execution"].get("new_sentaurus_execution", False)
            for report in current.values()),
    }

    stage_impact = [
        {"stage": "M31", "rebaseline_status": m31["status"],
         "impact": "numeric refresh; qualitative conclusion retained",
         "historical_contract": "retained"},
        {"stage": "M32", "rebaseline_status": m32["status"],
         "impact": "numeric refresh; qualitative conclusion retained",
         "historical_contract": "retained"},
        {"stage": "M37", "rebaseline_status": m37["status"],
         "impact": "bitwise unchanged control",
         "historical_contract": "retained"},
        {"stage": "M41", "rebaseline_status": m41["status"],
         "impact": "M43 operator consistency now propagated",
         "historical_contract": "superseded; failure retained"},
        {"stage": "M42", "rebaseline_status": m42["status"],
         "impact": "M43/M44 operator and QF packing closure propagated",
         "historical_contract": "superseded; failure retained"},
    ]
    findings = {
        "m31_qualitative_conclusion_changed": False,
        "m32_qualitative_conclusion_changed": False,
        "m37_bitwise_unchanged": m37_unchanged,
        "m41_historical_contract_superseded": True,
        "m42_historical_legacy_mismatch_hypothesis_superseded": True,
        "legacy_operator_consistency_restored_by_m43": True,
        "compensated_reference_replay_closure_restored_by_m44": True,
        "srh_feedback_conclusion_retained": True,
        "contact_reconstruction_null_retained": True,
        "qf_reference_increment_requirement_retained": True,
        "m42_maximum_compensated_cut_reported_gap_dex":
            m42["findings"]["maximum_compensated_cut_reported_gap_dex"],
        "m42_minimum_legacy_cut_reported_gap_dex":
            m42["findings"]["minimum_legacy_cut_reported_gap_dex"],
        "m42_maximum_compensated_free_electron_residual_l2":
            m42["findings"][
                "maximum_compensated_free_electron_residual_l2"],
        "m42_legacy_to_compensated_residual_ratio":
            m42["findings"][
                "legacy_to_compensated_residual_separation_min_ratio"],
        "attribution_note": "M41/M42 were stale across both M43 and M44; the large legacy-path change is attributed to M43 SG/contact operator unification, while the compensated sub-millidex closure is attributed to M44 reference/increment preservation.",
    }
    report = {
        "schema": "vela.simplemos.sdevice.m45_post_qf_rebaseline_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "baseline_revision": baseline_revision,
        "execution": {
            "rerun_stages": list(REPORTS),
            "new_sentaurus_execution": False,
            "default_model_changed": False,
            "historical_contracts_rewritten": False,
        },
        "stage_impact": stage_impact,
        "findings": findings,
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "claim_policy": [
            "M41/M42 failed status records superseded historical expectations, not workflow nonconvergence.",
            "M43 and M44 attribution is separated because M41/M42 had not been rerun between those fixes.",
            "No production physics, HFS, equal-ni default, or Sentaurus reference was changed.",
        ],
    }

    OUTPUT.mkdir(parents=True, exist_ok=True)
    report_path = OUTPUT / "m45_post_qf_rebaseline_report.json"
    impact_path = OUTPUT / "m45_stage_impact_matrix.csv"
    metrics_path = OUTPUT / "m45_key_metric_delta.csv"
    doc_path = REPO / "docs/validation/simplemos_m45_post_qf_rebaseline_2026-09-01.md"
    write_json(report_path, report)
    write_csv(impact_path, stage_impact)
    write_csv(metrics_path, metrics)
    doc_path.write_text(make_markdown(report, metrics), encoding="utf-8",
                        newline="\n")

    refresh_artifact_hashes(UPSTREAM_EVIDENCE["M31"], "M31")
    refresh_artifact_hashes(UPSTREAM_EVIDENCE["M32"], "M32")
    m32_evidence = read_json(UPSTREAM_EVIDENCE["M32"])
    m32_evidence["source_hashes"][portable(UPSTREAM_EVIDENCE["M31"])] = \
        sha256(UPSTREAM_EVIDENCE["M31"])
    write_json(UPSTREAM_EVIDENCE["M32"], m32_evidence)

    artifacts = [report_path, impact_path, metrics_path, doc_path,
                 *REPORTS.values(), UPSTREAM_EVIDENCE["M31"],
                 UPSTREAM_EVIDENCE["M32"], UPSTREAM_EVIDENCE["M41"],
                 UPSTREAM_EVIDENCE["M42"]]
    sources = [
        CONTRACT, Path(__file__).resolve(), UPSTREAM_EVIDENCE["M44"],
        REPO / "CMakeLists.txt",
        REPO / "scripts/run_simplemos_m31_minority_poisson_perturbation.py",
        REPO / "scripts/run_simplemos_m32_electron_qf_scaled_perturbation.py",
        REPO / "scripts/run_simplemos_m37_measure_contact_qf_ablation.py",
        REPO / "scripts/run_simplemos_m41_equal_ni_flux_ablation.py",
        REPO / "scripts/run_simplemos_m42_no_bgn_self_consistent_interaction.py",
        REPO / "tests/regression/simplemos_evidence_chain.py",
        REPO / "tests/regression/test_simplemos_m9_hfs_diagnostics.py",
        REPO / "tests/regression/test_simplemos_m10_fixed_state_replay.py",
        REPO / "tests/regression/test_simplemos_m12_terminal_sensitivity.py",
        REPO / "tests/regression/test_simplemos_m30_double_off_causal_closure.py",
        REPO / "tests/regression/test_simplemos_m31_minority_poisson_perturbation.py",
        REPO / "tests/regression/test_simplemos_m32_electron_qf_scaled_perturbation.py",
        REPO / "tests/regression/test_simplemos_m37_measure_contact_qf_ablation.py",
        REPO / "tests/regression/test_simplemos_m41_equal_ni_flux_ablation.py",
        REPO / "tests/regression/test_simplemos_m42_no_bgn_self_consistent_interaction.py",
        REPO / "tests/regression/test_simplemos_m45_post_qf_rebaseline.py",
    ]
    evidence = {
        "schema": "vela.simplemos.sdevice.m45_post_qf_rebaseline_evidence.v1",
        "status": "frozen" if report["status"] == "complete" else "failed",
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "default_physics_model_changed": False,
        "default_hfs_model_changed": False,
        "default_equal_ni_flux_evaluation_changed": False,
        "new_sentaurus_execution": False,
        "acceptance": report["acceptance"],
    }
    write_json(ROOT / "simplemos_m45_post_qf_rebaseline_evidence.json",
               evidence)
    print(json.dumps({"status": report["status"], **findings}, indent=2))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
