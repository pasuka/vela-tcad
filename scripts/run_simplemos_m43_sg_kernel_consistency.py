#!/usr/bin/env python3
"""Run the SimpleMOS M43 SG operator-consistency audit."""

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
            / "simplemos_m43_sg_kernel_consistency_contract_v1.json")
M41_OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
              / "m41_equal_ni_flux_ablation/self_consistent")
M42_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m42_no_bgn_self_consistent_interaction_evidence.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m43_sg_kernel_consistency")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "sg_kernel_consistency")
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


def validate_contract(contract: dict[str, Any]) -> None:
    if contract.get("schema") != \
            "vela.simplemos.sdevice.m43_sg_kernel_consistency_contract.v1":
        raise ValueError("unexpected M43 contract schema")
    if int(contract["matrix"]["state_count"]) != 4:
        raise ValueError("M43 requires four frozen states")
    policy = contract["policy"]
    if policy["default_physics_model_changed"]:
        raise ValueError("M43 may not change the default physics model")
    if policy["default_equal_ni_flux_evaluation_changed"]:
        raise ValueError("M43 may not change the default SG mode")


def drain_nodes(config: dict[str, Any]) -> set[int]:
    mesh_path = Path(config["mesh_file"])
    mesh = m40.read_json(mesh_path)
    for contact in mesh["contacts"]:
        if contact["name"].lower() == "drain":
            return {int(node) for node in contact["node_ids"]}
    raise ValueError("drain contact not found")


def run_case(mode: str, srh_name: str, force: bool) -> dict[str, Any]:
    source = M41_OUTPUT / SRH[srh_name] / mode / "20_gate"
    source_config = source / "config.json"
    state = source / "state.csv"
    if not source_config.is_file() or not state.is_file():
        raise FileNotFoundError(f"missing M41 frozen state for {mode} {srh_name}")

    root = OUTPUT / mode / srh_name
    root.mkdir(parents=True, exist_ok=True)
    sg_csv = root / "sg_edges.csv"
    sg_config_path = root / "sg_edges.json"
    functional_config_path = root / "terminal_functional.json"
    functional_status_path = root / "terminal_functional_status.json"

    source_deck = m40.read_json(source_config)
    expected_mode = source_deck["solver"]["bandgap_narrowing"][
        "equal_ni_flux_evaluation"]
    if expected_mode != MODES[mode]:
        raise ValueError(f"M41 mode mismatch for {mode}: {expected_mode}")

    sg_config = m10.probe_config(
        source_config, state, sg_csv, 0.05, 0.05,
        "sg_edge_flux_probe")
    sg_config["simplemos_m43"] = {
        "mode": mode, "srh": srh_name,
        "read_only": True, "operator_role": "edge_probe",
    }
    m40.write_json(sg_config_path, sg_config)
    if force or not sg_csv.is_file():
        m10.execute_runner(sg_config_path, RUNNER)

    functional_config = copy.deepcopy(sg_config)
    functional_config["simulation_type"] = "terminal_current_functional_probe"
    functional_config["contact"] = "drain"
    functional_config["simplemos_m43"]["operator_role"] = \
        "terminal_functional_and_contact_extractor"
    functional_config.pop("output_csv", None)
    m40.write_json(functional_config_path, functional_config)
    if force or not functional_status_path.is_file():
        status = m10.execute_runner(functional_config_path, RUNNER)
        m40.write_json(functional_status_path, status)
    else:
        status = m40.read_json(functional_status_path)

    cut = m10.drain_cut_current(
        m40.read_csv(sg_csv), drain_nodes(source_deck))
    functional = float(status["current_A_per_um"])
    extractor = float(status["contact_current_extractor_A_per_um"])
    cut_total = float(cut["total_A_per_um"])
    return {
        "mode": mode,
        "equal_ni_flux_evaluation": MODES[mode],
        "srh": srh_name,
        "state": m40.portable(state),
        "drain_cut_A_per_um": cut_total,
        "terminal_functional_A_per_um": functional,
        "contact_extractor_A_per_um": extractor,
        "cut_vs_terminal_functional_gap_dex": log_gap(cut_total, functional),
        "terminal_functional_vs_contact_extractor_gap_dex": log_gap(
            functional, extractor),
        "cut_minus_terminal_functional_A_per_um": cut_total - functional,
        "terminal_functional_minus_contact_extractor_A_per_um": (
            functional - extractor),
        "crossing_edge_count": int(cut["crossing_edge_count"]),
        "sg_edges": m40.portable(sg_csv),
        "terminal_status": m40.portable(functional_status_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    contract = m40.read_json(CONTRACT)
    validate_contract(contract)
    if not M42_EVIDENCE.is_file():
        raise FileNotFoundError(M42_EVIDENCE)

    rows = [run_case(mode, srh, args.force)
            for mode in MODES for srh in SRH]
    PORTABLE.mkdir(parents=True, exist_ok=True)
    matrix_path = PORTABLE / "m43_sg_kernel_consistency_matrix.csv"
    m10.write_csv(matrix_path, rows)
    acceptance_cfg = contract["acceptance"]
    max_cut_gap = max(row["cut_vs_terminal_functional_gap_dex"]
                      for row in rows)
    max_extractor_gap = max(
        row["terminal_functional_vs_contact_extractor_gap_dex"]
        for row in rows)
    checks = {
        "state_count": len(rows) == 4,
        "both_numerical_modes": {row["mode"] for row in rows} == set(MODES),
        "cut_matches_terminal_functional": max_cut_gap <= float(
            acceptance_cfg["maximum_cut_vs_terminal_functional_gap_dex"]),
        "terminal_functional_matches_contact_extractor":
            max_extractor_gap <= float(
                acceptance_cfg[
                    "maximum_terminal_functional_vs_contact_extractor_gap_dex"]),
        "all_contact_cuts_nonempty": all(
            row["crossing_edge_count"] > 0 for row in rows),
    }
    report = {
        "schema": "vela.simplemos.sdevice.m43_sg_kernel_consistency_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "execution": {
            "frozen_state_replay": True,
            "state_count": len(rows),
            "default_physics_model_changed": False,
            "default_equal_ni_flux_evaluation_changed": False,
            "new_sentaurus_execution": False,
        },
        "findings": {
            "maximum_cut_vs_terminal_functional_gap_dex": max_cut_gap,
            "maximum_terminal_functional_vs_contact_extractor_gap_dex":
                max_extractor_gap,
            "legacy_gap_closed": max(
                row["terminal_functional_vs_contact_extractor_gap_dex"]
                for row in rows if row["mode"] == "legacy") <= 0.005,
            "compensated_gap_not_regressed": max(
                row["terminal_functional_vs_contact_extractor_gap_dex"]
                for row in rows if row["mode"] == "compensated") <= 0.005,
        },
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "artifacts": {
            "matrix": m40.portable(matrix_path),
        },
        "policy": contract["policy"],
    }
    report_path = PORTABLE / "m43_sg_kernel_consistency_report.json"
    m40.write_json(report_path, report)
    artifacts = [report_path, matrix_path]
    sources = [
        CONTRACT, Path(__file__).resolve(), M42_EVIDENCE,
        REPO / "include/vela/discretization/ScharfetterGummel.h",
        REPO / "src/discretization/ScharfetterGummel.cpp",
        REPO / "src/equation/CoupledDDAssembler.cpp",
        REPO / "src/post/ContactCurrent.cpp",
        REPO / "src/tools/vela_example_runner.cpp",
    ]
    evidence = {
        "schema": "vela.simplemos.sdevice.m43_sg_kernel_consistency_evidence.v1",
        "status": "frozen" if report["status"] == "complete" else "failed",
        "artifacts": [{"path": m40.portable(path), "sha256": m10.sha256(path)}
                      for path in artifacts],
        "source_hashes": {m40.portable(path): m10.sha256(path)
                          for path in sources},
        "default_physics_model_changed": False,
        "default_equal_ni_flux_evaluation_changed": False,
        "acceptance": report["acceptance"],
    }
    evidence_path = (REPO / "reference_tcad/simplemos_sentaurus2022"
                     / "simplemos_m43_sg_kernel_consistency_evidence.json")
    m40.write_json(evidence_path, evidence)
    print(json.dumps(report["findings"], indent=2))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
