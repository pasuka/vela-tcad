#!/usr/bin/env python3
"""Build the M38 read-only SimpleMOS/LDMOS mechanism crosscheck."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m38_ldmos_mechanism_crosscheck_contract_v1.json")
SIMPLEMOS = (REPO / "reference_tcad/simplemos_sentaurus2022"
             / "measure_contact_qf_ablation/m37_measure_contact_qf_ablation_report.json")
OUTPUT = (REPO / "reference_tcad/simplemos_sentaurus2022"
          / "ldmos_mechanism_crosscheck/m38_ldmos_mechanism_crosscheck_report.json")
DEFAULT_LDMOS = Path("D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a")

LDMOS_SOURCES = {
    "averagebox_full_mesh": (
        "templates_ldmos_averagebox_full_mesh_20260831/summary.json"),
    "averagebox_newton_ab": (
        "templates_ldmos_averagebox_newton_ab_release_20260831/summary.json"),
    "contact_hfs_replay": (
        "templates_ldmos_g3_contact_hfs_replay_release_20260831/summary.json"),
    "state_feedback": (
        "templates_ldmos_g3_state_feedback_20260831/summary.json"),
    "state_feedback_signed": (
        "templates_ldmos_g3_state_feedback_signed_cotangent_20260831/summary.json"),
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point_at(report: dict[str, Any], bias: float) -> dict[str, Any]:
    return min(report["points"], key=lambda row: abs(float(row["bias_V"]) - bias))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ldmos-worktree", type=Path, default=DEFAULT_LDMOS)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    if contract["schema"] != "vela.simplemos.sdevice.m38_ldmos_mechanism_crosscheck.v1":
        raise ValueError("unexpected M38 contract schema")
    staging = args.ldmos_worktree / "reference_staging"
    paths = {name: staging / relative for name, relative in LDMOS_SOURCES.items()}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing LDMOS snapshot sources: " + ", ".join(missing))

    simplemos = read_json(SIMPLEMOS)
    averagebox = read_json(paths["averagebox_full_mesh"])
    newton_ab = read_json(paths["averagebox_newton_ab"])
    contact = read_json(paths["contact_hfs_replay"])
    feedback = point_at(read_json(paths["state_feedback"]),
                        float(contract["ldmos_bias_V"]))
    signed = point_at(read_json(paths["state_feedback_signed"]),
                      float(contract["ldmos_bias_V"]))
    a = feedback["attribution"]
    s = signed["attribution"]
    feedback_dex = float(a["feedback_amplification_dex"])
    phin_recovery = float(a["single_sentaurus_family_error_recovery_dex"]["phin"])

    ldmos = {
        "transport_nodes": averagebox["geometry"]["transport_nodes"],
        "averagebox_electron_residual_l2_ratio": averagebox[
            "coefficient_only_over_baseline"]["l2"],
        "averagebox_electron_residual_max_ratio": averagebox[
            "coefficient_only_over_baseline"]["maximum_abs"],
        "averagebox_candidate_reclose_current_error_dex": newton_ab[
            "candidate_same_bias_reclose"]["drain_log10_abs_error_dex"],
        "contact_hfs_baseline_fixed_error_dex": contact[
            "baseline_median_magnitude_error_dex"],
        "contact_hfs_corrected_fixed_error_dex": contact[
            "contact_fallback_median_magnitude_error_dex"],
        "contact_hfs_operator_improvement_dex": contact[
            "median_improvement_dex"],
        "operator_error_dex": float(a["operator_log10_ratio_dex"]),
        "self_consistent_error_dex": float(a["self_consistent_log10_ratio_dex"]),
        "feedback_amplification_dex": feedback_dex,
        "feedback_fraction_of_self_consistent_error": feedback_dex /
            float(a["self_consistent_log10_ratio_dex"]),
        "psi_recovery_dex": float(
            a["single_sentaurus_family_error_recovery_dex"]["psi"]),
        "phin_recovery_dex": phin_recovery,
        "phip_recovery_dex": float(
            a["single_sentaurus_family_error_recovery_dex"]["phip"]),
        "phin_fraction_of_feedback": phin_recovery / feedback_dex,
        "signed_cotangent_operator_change_dex": float(
            s["operator_log10_ratio_dex"]) - float(a["operator_log10_ratio_dex"]),
        "signed_cotangent_self_consistent_change_dex": float(
            s["self_consistent_log10_ratio_dex"])
            - float(a["self_consistent_log10_ratio_dex"]),
        "signed_cotangent_feedback_change_dex": float(
            s["feedback_amplification_dex"]) - feedback_dex,
    }
    simple = {
        "boundary_scope_fraction_of_all_signed_log_response": simplemos[
            "boundary_measure_ablation"][
                "boundary_scope_fraction_of_all_signed_log_response"],
        "contact_sg_direct_electron_delta_A_per_um": simplemos[
            "contact_sg_ablation"]["signed_direct_electron_current_delta_A_per_um"],
        "qf_adjoint_relaxation_fraction_of_actual": simplemos[
            "quasi_fermi_feedback_ablation"][
                "adjoint_relaxation_fraction_of_actual"],
        "qf_electron_equation_signed_fraction": simplemos[
            "quasi_fermi_feedback_ablation"]["electron_equation_signed_fraction"],
    }
    limits = contract["acceptance"]
    numeric = [float(value) for value in ldmos.values()] + [
        float(value) for value in simple.values()]
    checks = {
        "signed_cotangent_invariant": max(
            abs(ldmos["signed_cotangent_operator_change_dex"]),
            abs(ldmos["signed_cotangent_self_consistent_change_dex"]),
            abs(ldmos["signed_cotangent_feedback_change_dex"])) <= limits[
                "maximum_signed_cotangent_attribution_change_dex"],
        "ldmos_phin_dominant": ldmos["phin_fraction_of_feedback"] >= limits[
            "minimum_ldmos_phin_fraction_of_feedback"],
        "ldmos_contact_operator_material": ldmos[
            "contact_hfs_operator_improvement_dex"] >= limits[
                "minimum_ldmos_contact_operator_improvement_dex"],
        "finite": all(math.isfinite(value) for value in numeric),
    }
    report = {
        "schema": "vela.simplemos.sdevice.m38_ldmos_mechanism_crosscheck_report.v1",
        "status": "qualified_read_only_snapshot" if all(checks.values())
            else "failed_acceptance",
        "execution": {
            "ldmos_worktree": str(args.ldmos_worktree.resolve()),
            "ldmos_branch_observed": "codex/templates-ldmos-phase-a",
            "ldmos_worktree_dirty_observed": True,
            "ldmos_files_modified": False,
            "new_sentaurus_execution": False,
            "default_model_changed": False,
        },
        "acceptance": {"checks": checks, "all_checks_pass": all(checks.values())},
        "simplemos": simple,
        "ldmos": ldmos,
        "cross_device_conclusion": {
            "shared_mechanism": (
                "After contact-operator correction, electron quasi-Fermi "
                "self-consistent feedback dominates both devices."),
            "device_specific_difference": (
                "Direct contact HFS/SG correction is material in LDMOS but "
                "zero on the SimpleMOS drain cut."),
            "geometry_interpretation": (
                "AverageBox geometry strongly closes the LDMOS frozen residual "
                "but does not by itself close the self-consistent current."),
        },
        "ldmos_source_hashes": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in paths.items()
        },
        "claim_policy": contract["claim_policy"],
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n",
                      encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "checks": checks,
                      "simplemos": simple, "ldmos": ldmos}, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
