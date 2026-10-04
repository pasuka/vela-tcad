#!/usr/bin/env python3
"""Run the SimpleMOS M44 quasi-Fermi coordinate-consistency audit."""

from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from typing import Any

import run_simplemos_m10_fixed_state_replay as m10
import run_simplemos_m40_true_no_bgn_factorial as m40


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "build-release/vela_example_runner.exe"
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m44_qf_coordinate_consistency_contract_v1.json")
M41_OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
              / "m41_equal_ni_flux_ablation/self_consistent")
M43_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m43_sg_kernel_consistency_evidence.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m44_qf_coordinate_consistency")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "qf_coordinate_consistency")
ELEMENTARY_CHARGE_C = 1.602176634e-19
MODES = {
    "legacy": "legacy_factor_difference",
    "compensated": "compensated_log_expm1",
}
SRH = {
    "srh_on": "bgn_off_srh_on",
    "srh_off": "bgn_off_srh_off",
}


def log_gap(left: float, right: float) -> float:
    left_abs = abs(left)
    right_abs = abs(right)
    if left_abs == 0.0 and right_abs == 0.0:
        return 0.0
    if left_abs <= 0.0 or right_abs <= 0.0:
        return math.inf
    return abs(math.log10(left_abs) - math.log10(right_abs))


def relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def validate_contract(contract: dict[str, Any]) -> None:
    if contract.get("schema") != \
            "vela.simplemos.sdevice.m44_qf_coordinate_consistency_contract.v1":
        raise ValueError("unexpected M44 contract schema")
    if int(contract["matrix"]["state_count"]) != 4:
        raise ValueError("M44 requires four frozen states")
    policy = contract["policy"]
    for key in ("default_physics_model_changed",
                "default_hfs_model_changed",
                "default_equal_ni_flux_evaluation_changed",
                "new_sentaurus_execution"):
        if policy[key]:
            raise ValueError(f"M44 policy forbids {key}")


def drain_nodes(config: dict[str, Any]) -> set[int]:
    mesh = m40.read_json(Path(config["mesh_file"]))
    for contact in mesh["contacts"]:
        if contact["name"].lower() == "drain":
            return {int(node) for node in contact["node_ids"]}
    raise ValueError("drain contact not found")


def run_case(mode: str, srh_name: str, force: bool) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    source = M41_OUTPUT / SRH[srh_name] / mode / "20_gate"
    source_config = source / "config.json"
    state = source / "state.csv"
    if not source_config.is_file() or not state.is_file():
        raise FileNotFoundError(f"missing M41 frozen state for {mode} {srh_name}")

    root = OUTPUT / mode / srh_name
    root.mkdir(parents=True, exist_ok=True)
    sg_csv = root / "sg_edges.csv"
    contact_csv = root / "contact_edges.csv"
    sg_config_path = root / "sg_edges.json"
    functional_config_path = root / "terminal_functional.json"
    functional_status_path = root / "terminal_functional_status.json"

    source_deck = m40.read_json(source_config)
    expected_mode = source_deck["solver"]["bandgap_narrowing"][
        "equal_ni_flux_evaluation"]
    if expected_mode != MODES[mode]:
        raise ValueError(f"M41 mode mismatch for {mode}: {expected_mode}")

    sg_config = m10.probe_config(
        source_config, state, sg_csv, 0.05, 0.05, "sg_edge_flux_probe")
    sg_config["simplemos_m44"] = {
        "mode": mode,
        "srh": srh_name,
        "read_only": True,
        "operator_role": "coupled_sg_edge_operator",
    }
    m40.write_json(sg_config_path, sg_config)
    if force or not sg_csv.is_file():
        m10.execute_runner(sg_config_path, RUNNER)

    functional_config = copy.deepcopy(sg_config)
    functional_config["simulation_type"] = "terminal_current_functional_probe"
    functional_config["contact"] = "drain"
    functional_config["contact_edge_output_csv"] = str(contact_csv)
    functional_config["simplemos_m44"]["operator_role"] = \
        "terminal_functional_and_contact_extractor"
    functional_config.pop("output_csv", None)
    m40.write_json(functional_config_path, functional_config)
    if force or not functional_status_path.is_file() or not contact_csv.is_file():
        status = m10.execute_runner(functional_config_path, RUNNER)
        m40.write_json(functional_status_path, status)
    else:
        status = m40.read_json(functional_status_path)

    sg_rows = m40.read_csv(sg_csv)
    sg_by_edge = {int(row["edge_id"]): row for row in sg_rows}
    edge_rows: list[dict[str, Any]] = []
    for contact_row in m40.read_csv(contact_csv):
        edge_id = int(contact_row["edge_id"])
        operator = sg_by_edge[edge_id]
        outward = float(contact_row["outward_sign"])
        operator_electron_current = (
            -ELEMENTARY_CHARGE_C * outward
            * float(operator["electron_particle_line_flux_per_m_s"]) * 1.0e-6)
        operator_hole_current = (
            -ELEMENTARY_CHARGE_C * outward
            * float(operator["hole_particle_line_flux_per_m_s"]) * 1.0e-6)
        operator_total_current = operator_electron_current - operator_hole_current
        contact_electron_current = float(contact_row["electron_current_A_per_um"])
        contact_hole_current = float(contact_row["hole_current_A_per_um"])
        contact_total_current = float(contact_row["total_current_A_per_um"])
        drive_operator = float(operator["electron_mobility_drive_V_m"])
        drive_contact = float(contact_row["electron_mobility_drive_V_m"])
        mobility_operator = float(operator["electron_mobility_m2_V_s"])
        mobility_contact = float(contact_row["electron_mobility_m2_V_s"])
        phin0_operator = float(operator["electron_sg_phin0_relative_V"])
        phin1_operator = float(operator["electron_sg_phin1_relative_V"])
        phin0_contact = float(contact_row["electron_sg_phin0_relative_V"])
        phin1_contact = float(contact_row["electron_sg_phin1_relative_V"])
        edge_rows.append({
            "mode": mode,
            "srh": srh_name,
            "edge_id": edge_id,
            "node0": int(contact_row["node0"]),
            "node1": int(contact_row["node1"]),
            "outward_sign": outward,
            "operator_electron_current_A_per_um": operator_electron_current,
            "contact_electron_current_A_per_um": contact_electron_current,
            "operator_hole_current_A_per_um": operator_hole_current,
            "contact_hole_current_A_per_um": contact_hole_current,
            "operator_total_current_A_per_um": operator_total_current,
            "contact_total_current_A_per_um": contact_total_current,
            "total_current_relative_difference": relative_difference(
                operator_total_current, contact_total_current),
            "operator_electron_mobility_drive_V_m": drive_operator,
            "contact_electron_mobility_drive_V_m": drive_contact,
            "mobility_drive_relative_difference": relative_difference(
                drive_operator, drive_contact),
            "operator_electron_mobility_m2_V_s": mobility_operator,
            "contact_electron_mobility_m2_V_s": mobility_contact,
            "mobility_relative_difference": relative_difference(
                mobility_operator, mobility_contact),
            "operator_phin0_relative_V": phin0_operator,
            "contact_phin0_relative_V": phin0_contact,
            "phin0_relative_difference_V": phin0_operator - phin0_contact,
            "operator_phin1_relative_V": phin1_operator,
            "contact_phin1_relative_V": phin1_contact,
            "phin1_relative_difference_V": phin1_operator - phin1_contact,
        })

    cut = m10.drain_cut_current(sg_rows, drain_nodes(source_deck))
    functional = float(status["current_A_per_um"])
    extractor = float(status["contact_current_extractor_A_per_um"])
    extractor_compensated = float(
        status["contact_current_extractor_compensated_A_per_um"])
    extractor_long_double = float(
        status["contact_current_extractor_long_double_A_per_um"])
    cut_total = float(cut["total_A_per_um"])
    summary = {
        "mode": mode,
        "equal_ni_flux_evaluation": MODES[mode],
        "srh": srh_name,
        "state": m40.portable(state),
        "drain_cut_A_per_um": cut_total,
        "terminal_functional_A_per_um": functional,
        "contact_extractor_A_per_um": extractor,
        "contact_extractor_compensated_A_per_um": extractor_compensated,
        "contact_extractor_long_double_A_per_um": extractor_long_double,
        "cut_vs_terminal_functional_gap_dex": log_gap(cut_total, functional),
        "terminal_functional_vs_contact_extractor_gap_dex": log_gap(
            functional, extractor),
        "terminal_functional_vs_contact_extractor_relative_difference":
            relative_difference(functional, extractor),
        "maximum_contact_edge_current_relative_difference": max(
            row["total_current_relative_difference"] for row in edge_rows),
        "maximum_qf_relative_coordinate_difference_V": max(
            max(abs(row["phin0_relative_difference_V"]),
                abs(row["phin1_relative_difference_V"]))
            for row in edge_rows),
        "maximum_mobility_drive_relative_difference": max(
            row["mobility_drive_relative_difference"] for row in edge_rows),
        "maximum_mobility_relative_difference": max(
            row["mobility_relative_difference"] for row in edge_rows),
        "positive_couple_contact_edge_count": len(edge_rows),
        "sg_crossing_edge_count": int(cut["crossing_edge_count"]),
        "sg_edges": m40.portable(sg_csv),
        "contact_edges": m40.portable(contact_csv),
        "terminal_status": m40.portable(functional_status_path),
    }
    return summary, edge_rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    contract = m40.read_json(CONTRACT)
    validate_contract(contract)
    if not M43_EVIDENCE.is_file():
        raise FileNotFoundError(M43_EVIDENCE)

    summaries: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    for mode in MODES:
        for srh in SRH:
            summary, case_edges = run_case(mode, srh, args.force)
            summaries.append(summary)
            edge_rows.extend(case_edges)

    PORTABLE.mkdir(parents=True, exist_ok=True)
    matrix_path = PORTABLE / "m44_qf_coordinate_consistency_matrix.csv"
    edge_matrix_path = PORTABLE / "m44_qf_coordinate_edge_matrix.csv"
    m10.write_csv(matrix_path, summaries)
    m10.write_csv(edge_matrix_path, edge_rows)

    acceptance_cfg = contract["acceptance"]
    max_cut_gap = max(row["cut_vs_terminal_functional_gap_dex"]
                      for row in summaries)
    max_extractor_gap = max(
        row["terminal_functional_vs_contact_extractor_gap_dex"]
        for row in summaries)
    max_edge_current = max(
        row["maximum_contact_edge_current_relative_difference"]
        for row in summaries)
    max_qf_coordinate = max(
        row["maximum_qf_relative_coordinate_difference_V"]
        for row in summaries)
    max_drive = max(
        row["maximum_mobility_drive_relative_difference"]
        for row in summaries)
    max_mobility = max(
        row["maximum_mobility_relative_difference"]
        for row in summaries)
    checks = {
        "state_count": len(summaries) == 4,
        "both_numerical_modes": {row["mode"] for row in summaries} == set(MODES),
        "cut_matches_terminal_functional": max_cut_gap <= float(
            acceptance_cfg["maximum_cut_vs_terminal_functional_gap_dex"]),
        "terminal_functional_matches_contact_extractor":
            max_extractor_gap <= float(
                acceptance_cfg[
                    "maximum_terminal_functional_vs_contact_extractor_gap_dex"]),
        "contact_edges_match": max_edge_current <= float(
            acceptance_cfg["maximum_contact_edge_current_relative_difference"]),
        "qf_relative_coordinates_match": max_qf_coordinate <= float(
            acceptance_cfg["maximum_qf_relative_coordinate_difference_V"]),
        "mobility_drives_match": max_drive <= float(
            acceptance_cfg["maximum_mobility_drive_relative_difference"]),
        "mobilities_match": max_mobility <= float(
            acceptance_cfg["maximum_mobility_relative_difference"]),
        "all_contact_cuts_nonempty": all(
            row["positive_couple_contact_edge_count"] > 0 for row in summaries),
    }
    baseline_gap = float(contract["baseline"][
        "m43_maximum_terminal_functional_vs_contact_extractor_gap_dex"])
    report = {
        "schema": "vela.simplemos.sdevice.m44_qf_coordinate_consistency_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "execution": {
            "frozen_state_replay": True,
            "state_count": len(summaries),
            "default_physics_model_changed": False,
            "default_hfs_model_changed": False,
            "default_equal_ni_flux_evaluation_changed": False,
            "new_sentaurus_execution": False,
        },
        "findings": {
            "maximum_cut_vs_terminal_functional_gap_dex": max_cut_gap,
            "maximum_terminal_functional_vs_contact_extractor_gap_dex":
                max_extractor_gap,
            "m43_baseline_gap_dex": baseline_gap,
            "minimum_gap_reduction_factor": baseline_gap /
                max(max_extractor_gap, 1.0e-300),
            "maximum_contact_edge_current_relative_difference":
                max_edge_current,
            "maximum_qf_relative_coordinate_difference_V": max_qf_coordinate,
            "maximum_mobility_drive_relative_difference": max_drive,
            "maximum_mobility_relative_difference": max_mobility,
            "mobility_reconstruction_was_root_cause": False,
            "reference_increment_packing_was_root_cause": True,
        },
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "artifacts": {
            "matrix": m40.portable(matrix_path),
            "edge_matrix": m40.portable(edge_matrix_path),
        },
        "policy": contract["policy"],
    }
    report_path = PORTABLE / "m44_qf_coordinate_consistency_report.json"
    m40.write_json(report_path, report)
    artifacts = [report_path, matrix_path, edge_matrix_path]
    sources = [
        CONTRACT, Path(__file__).resolve(), M43_EVIDENCE,
        REPO / "CMakeLists.txt",
        REPO / "src/equation/CoupledDDAssembler.cpp",
        REPO / "src/solver/NewtonSolver.cpp",
        REPO / "src/simulation/DCSweep.cpp",
        REPO / "include/vela/post/ContactCurrent.h",
        REPO / "src/post/ContactCurrent.cpp",
        REPO / "src/tools/vela_example_runner.cpp",
        REPO / "tests/test_newton_solver.cpp",
        REPO / "tests/regression/simplemos_evidence_chain.py",
        REPO / "tests/regression/test_simplemos_m31_minority_poisson_perturbation.py",
        REPO / "tests/regression/test_simplemos_m32_electron_qf_scaled_perturbation.py",
        REPO / "tests/regression/test_simplemos_m33_region_resolved_interface_assembly.py",
        REPO / "tests/regression/test_simplemos_m35_signed_average_box_assembly.py",
        REPO / "tests/regression/test_simplemos_m36_boundary_contact_measure_audit.py",
        REPO / "tests/regression/test_simplemos_m37_measure_contact_qf_ablation.py",
        REPO / "tests/regression/test_simplemos_m39_bgn_chain_first_divergence.py",
        REPO / "tests/regression/test_simplemos_m40_true_no_bgn_factorial.py",
    ]
    evidence = {
        "schema": "vela.simplemos.sdevice.m44_qf_coordinate_consistency_evidence.v1",
        "status": "frozen" if report["status"] == "complete" else "failed",
        "artifacts": [{"path": m40.portable(path), "sha256": m10.sha256(path)}
                      for path in artifacts],
        "source_hashes": {m40.portable(path): m10.sha256(path)
                          for path in sources},
        "default_physics_model_changed": False,
        "default_hfs_model_changed": False,
        "default_equal_ni_flux_evaluation_changed": False,
        "acceptance": report["acceptance"],
    }
    evidence_path = (REPO / "reference_tcad/simplemos_sentaurus2022"
                     / "simplemos_m44_qf_coordinate_consistency_evidence.json")
    m40.write_json(evidence_path, evidence)
    print(json.dumps(report["findings"], indent=2))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
