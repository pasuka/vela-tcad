#!/usr/bin/env python3
"""Run the M37 independent measure/contact-SG/QF-feedback ablation."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
from typing import Any


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m37_measure_contact_qf_ablation_contract_v1.json")
SOURCE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m8a_model_ablation/vela/full/workflow/vd_0p05")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m37_measure_contact_qf_ablation")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "measure_contact_qf_ablation")
M35 = (REPO / "reference_tcad/simplemos_sentaurus2022"
       / "signed_average_box_assembly/m35_signed_average_box_assembly_report.json")
M22 = (REPO / "reference_tcad/simplemos_sentaurus2022"
       / "n23_hfs_deep_off/m22_n23_hfs_deep_off_report.json")
M22_EDGES = (REPO / "reference_tcad/simplemos_sentaurus2022"
             / "n23_hfs_deep_off/m22_key_state_drain_cut_edges.csv")
M26 = (REPO / "reference_tcad/simplemos_sentaurus2022"
       / "first_layer_feedback_audit/m26_first_layer_feedback_audit_report.json")
RUNNER = REPO / "build-release/vela_example_runner.exe"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def stage_config(stage: str) -> Path:
    run_dir = OUTPUT / "external_boundary_contact_support_signed"
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = read_json(SOURCE / f"{stage}.json")
    cfg["solver"]["region_resolved_interface_assembly"] = {
        "enabled": False,
        "poisson_edge_coupling": True,
        "transport_edge_coupling": True,
        "transport_node_volume": False,
        "transport_signed_average_box_node_volume": True,
        "transport_signed_average_box_node_volume_scope":
            "external_boundary_contact_support",
    }
    cfg["output_csv"] = str((run_dir / f"{stage}.csv").resolve())
    cfg["log_file"] = str((run_dir / f"{stage}.log").resolve())
    cfg["sweep"]["write_state_file"] = str(
        (run_dir / f"{stage}_accepted_state.csv").resolve())
    if stage == "10_drain_ramp":
        cfg["sweep"]["initial_state_file"] = str(
            (run_dir / "00_equilibrium_accepted_state.csv").resolve())
    elif stage == "20_gate_sweep":
        cfg["sweep"]["initial_state_file"] = str(
            (run_dir / "10_drain_ramp_accepted_state.csv").resolve())
        cfg["sweep"].update({
            "start": 0.0,
            "stop": 0.05,
            "step": 0.05,
            "max_step": 0.05,
            "bias_points": [0.0, 0.05],
        })
    cfg["simplemos_m37"] = {
        "ablation": "boundary_measure",
        "scope": "external_boundary_contact_support",
        "default_model_changed": False,
    }
    path = run_dir / f"{stage}.json"
    write_json(path, cfg)
    return path


def execute(path: Path, runner: Path) -> None:
    result = subprocess.run(
        [str(runner), "--config", str(path)], cwd=REPO, text=True,
        capture_output=True, encoding="utf-8", errors="replace")
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {path}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")


def boundary_current(run: bool, runner: Path) -> tuple[float, bool]:
    configs = [stage_config(stage) for stage in (
        "00_equilibrium", "10_drain_ramp", "20_gate_sweep")]
    if run:
        for config in configs:
            execute(config, runner)
    curve = OUTPUT / "external_boundary_contact_support_signed/20_gate_sweep.csv"
    if not curve.exists():
        raise FileNotFoundError(f"missing {curve}; rerun with --run")
    row = min(read_csv(curve), key=lambda item: abs(float(item["bias_V"]) - 0.05))
    return abs(float(row["current_total_A_per_um"])), (
        row["converged"].lower() in ("1", "true"))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    if contract["schema"] != "vela.simplemos.sdevice.m37_measure_contact_qf_ablation.v1":
        raise ValueError("unexpected M37 contract schema")

    m35 = read_json(M35)
    m35_by_name = {row["variant"]: row for row in m35["factorial_response"]}
    boundary_value, boundary_converged = boundary_current(
        args.run, args.runner.resolve())
    baseline_value = float(m35_by_name["p1_t1_v0"]["current_A_per_um"])
    all_value = float(m35_by_name["p1_t1_v1"]["current_A_per_um"])
    sentaurus = float(m35["reference"]["sentaurus_current_A_per_um"])
    measure_rows = []
    for name, current, converged in (
        ("historical_node_volume", baseline_value, True),
        ("external_boundary_contact_support_signed", boundary_value,
         boundary_converged),
        ("all_transport_nodes_signed", all_value, True),
    ):
        measure_rows.append({
            "variant": name,
            "converged": converged,
            "current_A_per_um": current,
            "log_current_shift_vs_historical_dex": math.log10(
                current / baseline_value),
            "absolute_error_dex": abs(math.log10(current / sentaurus)),
        })
    boundary_shift = measure_rows[1]["log_current_shift_vs_historical_dex"]
    all_shift = measure_rows[2]["log_current_shift_vs_historical_dex"]

    m22 = read_json(M22)
    contact_edges = [row for row in read_csv(M22_EDGES)
                     if math.isclose(float(row["gate_voltage_V"]), 0.05)
                     and row["source_state_variant"] == "full"]
    contact_sg = {
        "drain_cut_edge_count": len(contact_edges),
        "maximum_mobility_log10_response_dex": max(
            abs(float(row["mobility_log10_ratio_on_over_off"]))
            for row in contact_edges),
        "signed_direct_electron_current_delta_A_per_um": sum(
            float(row["direct_electron_current_delta_A_per_um"])
            for row in contact_edges),
        "maximum_absolute_direct_electron_edge_delta_A_per_um": max(
            abs(float(row["direct_electron_current_delta_A_per_um"]))
            for row in contact_edges),
        "global_frozen_operator_delta_A_per_um": m22["key_state"][
            "full_state_direct_hfs_delta_A_per_um"],
        "global_frozen_operator_effect_dex": m22["key_state"][
            "full_state_frozen_hfs_effect_dex"],
    }

    m26 = read_json(M26)
    turn_on = next(row for row in m26["operator_directions"]
                   if row["direction"] == "turn_on")
    actual = float(turn_on["actual_self_consistent_current_change_A_per_um"])
    direct = float(turn_on["direct_frozen_operator_current_change_A_per_um"])
    relaxation = float(turn_on["adjoint_relaxation_current_change_A_per_um"])
    qf_feedback = {
        "actual_self_consistent_current_change_A_per_um": actual,
        "direct_frozen_operator_current_change_A_per_um": direct,
        "adjoint_relaxation_current_change_A_per_um": relaxation,
        "direct_fraction_of_actual": abs(direct) / max(abs(actual), 1e-300),
        "adjoint_relaxation_fraction_of_actual": abs(relaxation) /
            max(abs(actual), 1e-300),
        "first_order_relative_error": turn_on["first_order_relative_error"],
        "electron_equation_signed_fraction": m26["findings"][
            "turn_on_electron_equation_signed_fraction"],
    }

    expected = contract["acceptance"]
    baseline_replay = (
        abs(float(measure_rows[0]["current_A_per_um"]) - baseline_value) /
        baseline_value)
    all_replay = (
        abs(float(measure_rows[2]["current_A_per_um"]) - all_value) /
        all_value)
    numeric = [float(value) for row in measure_rows for key, value in row.items()
               if key not in ("variant", "converged")]
    numeric += [float(value) for value in contact_sg.values()]
    numeric += [float(value) for value in qf_feedback.values()]
    checks = {
        "boundary_variants": len(measure_rows) == expected["boundary_variant_count"],
        "boundary_converged": all(bool(row["converged"]) for row in measure_rows),
        "existing_endpoint_replay": max(baseline_replay, all_replay) <= expected[
            "maximum_existing_endpoint_replay_relative_error"],
        "drain_cut": len(contact_edges) >= expected["minimum_drain_cut_edges"],
        "contact_mobility": contact_sg[
            "maximum_mobility_log10_response_dex"] <= expected[
                "maximum_contact_cut_mobility_log10_response_dex"],
        "m26_first_order": qf_feedback["first_order_relative_error"] <= expected[
            "maximum_m26_first_order_relative_error"],
        "finite": all(math.isfinite(value) for value in numeric),
    }
    report = {
        "schema": "vela.simplemos.sdevice.m37_measure_contact_qf_ablation_report.v1",
        "status": "complete" if all(checks.values()) else "failed_acceptance",
        "execution": {
            "new_sentaurus_execution": False,
            "new_self_consistent_workflows": 1,
            "cpp_changed": True,
            "default_model_changed": False,
        },
        "acceptance": {"checks": checks, "all_checks_pass": all(checks.values())},
        "boundary_measure_ablation": {
            "variants": measure_rows,
            "boundary_scope_fraction_of_all_signed_log_response":
                boundary_shift / all_shift if all_shift else 0.0,
        },
        "contact_sg_ablation": contact_sg,
        "quasi_fermi_feedback_ablation": qf_feedback,
        "claim_policy": contract["claim_policy"],
        "artifacts": {
            "boundary_measure_response": portable(
                PORTABLE / "m37_boundary_measure_response.csv"),
        },
    }
    PORTABLE.mkdir(parents=True, exist_ok=True)
    write_csv(PORTABLE / "m37_boundary_measure_response.csv", measure_rows)
    write_json(PORTABLE / "m37_measure_contact_qf_ablation_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "checks": checks,
        "boundary_measure": report["boundary_measure_ablation"],
        "contact_sg": contact_sg,
        "qf_feedback": qf_feedback,
    }, indent=2))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
