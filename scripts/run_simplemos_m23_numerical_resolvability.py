#!/usr/bin/env python3
"""Audit SimpleMOS SG/QF numerical resolvability on the M21/M22 key states."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402


CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m23_numerical_resolvability_contract_v1.json")
M21_RAW = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
           / "m21_incident_edge_balance/n21_vd_0p05_vg_0p8")
M21_PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "incident_edge_balance/m21_key_state_incident_edges.csv")
M22_RAW = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
           / "m22_n23_hfs_deep_off")
M22_PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "n23_hfs_deep_off")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m23_numerical_resolvability")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "numerical_resolvability")
RUNNER = REPO / "build-release/vela_example_runner.exe"
KB = 1.380649e-23
Q = 1.602176634e-19
TEMPERATURE_K = 300.0


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
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def validate(contract: dict[str, Any]) -> None:
    if contract.get("schema") != "vela.simplemos.sdevice.m23_numerical_resolvability.v1":
        raise ValueError("unexpected M23 contract schema")
    if int(contract["precision"]["decimal_digits"]) < 80:
        raise ValueError("M23 requires at least 80 decimal digits")
    if contract["solver_path_controls"] != [
            "baseline_repeat", "tight_tolerance", "same_state_reclosure",
            "reverse_from_vg_0p1"]:
        raise ValueError("unexpected M23 solver-path controls")


def exact_float(value: str | float) -> Decimal:
    return Decimal.from_float(float(value))


def decimal_expm1(value: Decimal) -> Decimal:
    return value.exp() - Decimal(1)


def nextafter_step(value: float, step: int) -> float:
    result = value
    direction = math.inf if step > 0 else -math.inf
    for _ in range(abs(step)):
        result = math.nextafter(result, direction)
    return result


def state_map(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["node_id"]): row for row in read_csv(path)}


def edge_map(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["edge_id"]): row for row in read_csv(path)}


def qf_components(case_id: str, row: dict[str, str], states: dict[int, dict[str, str]]
                  ) -> tuple[Decimal, Decimal, float, float, str]:
    del states
    endpoint0 = float(row["electron_sg_phin0_relative_V"])
    endpoint1 = float(row["electron_sg_phin1_relative_V"])
    internal_delta = exact_float(endpoint1) - exact_float(endpoint0)
    representation = ("contact_basin_increment" if case_id == "m22_deep_off"
                      else "contact_basin_loaded_snapshot")
    absolute_delta = exact_float(row["phin1_V"]) - exact_float(row["phin0_V"])
    return internal_delta, absolute_delta, endpoint0, endpoint1, representation


def analyze_edge(case_id: str, state_variant: str, row: dict[str, str],
                 states: dict[int, dict[str, str]], precision: dict[str, Any]
                 ) -> dict[str, Any]:
    decimal_digits = int(precision["decimal_digits"])
    with localcontext() as context:
        context.prec = decimal_digits
        internal_delta, absolute_delta, endpoint0, endpoint1, representation = (
            qf_components(case_id, row, states))
        signed_difference = exact_float(row["electron_sg_signed_difference_m3"])
        reconstructed = exact_float(row["electron_sg_reconstructed_flux"])
        if signed_difference == 0:
            raise ValueError(f"edge {row['edge_id']} has zero subtractive SG difference")
        coefficient = reconstructed / signed_difference
        right_term = exact_float(row["electron_sg_right_term_m3"])
        left_term = exact_float(row["electron_sg_left_term_m3"])
        thermal_voltage = exact_float(KB * TEMPERATURE_K / Q)
        log_imbalance = exact_float(row["electron_sg_log_left_over_right"])
        right_factor_flux = exact_float(row["electron_sg_right_factor_flux"])
        decimal_stable = right_factor_flux * decimal_expm1(log_imbalance)
        decimal_direct = coefficient * (left_term - right_term)
        production_stable = exact_float(row["electron_sg_stable_flux"])
        production_direct = reconstructed
        scale = max(abs(decimal_stable), Decimal("1e-300"))
        production_relative_error = abs(production_stable - decimal_stable) / scale
        direct_relative_error = abs(production_direct - decimal_stable) / scale
        decimal_direct_relative_error = abs(decimal_direct - decimal_stable) / scale

        base_endpoint1 = exact_float(endpoint1)
        perturbed_fluxes: list[Decimal] = []
        for step0 in precision["endpoint_ulp_steps"]:
            for step1 in precision["endpoint_ulp_steps"]:
                perturbed0 = nextafter_step(endpoint0, int(step0))
                perturbed1 = nextafter_step(endpoint1, int(step1))
                delta0 = exact_float(perturbed0) - exact_float(endpoint0)
                delta1 = exact_float(perturbed1) - exact_float(endpoint1)
                perturbed_log = log_imbalance + (delta1 - delta0) / thermal_voltage
                perturbed_right_factor = right_factor_flux * (
                    -(exact_float(perturbed1) - base_endpoint1)
                    / thermal_voltage).exp()
                perturbed_fluxes.append(
                    perturbed_right_factor * decimal_expm1(perturbed_log))
        ulp_relative_envelope = max(
            abs(value - decimal_stable) / scale for value in perturbed_fluxes)
        endpoint_ulp = max(math.ulp(endpoint0), math.ulp(endpoint1))
        qf_drop_ulps = abs(internal_delta) / max(
            exact_float(endpoint_ulp), Decimal("1e-300"))
        absolute_loss = abs(absolute_delta - internal_delta)
        absolute_loss_relative = absolute_loss / max(
            abs(internal_delta), Decimal("1e-300"))
        collapsed = absolute_delta == 0 and internal_delta != 0
        resolved = (
            production_relative_error <= Decimal(str(
                precision["maximum_production_vs_decimal_relative_error"]))
            and ulp_relative_envelope <= Decimal(str(
                precision["maximum_ulp_flux_relative_envelope_for_resolved"]))
            and qf_drop_ulps >= Decimal(str(
                precision["minimum_qf_drop_ulps_for_resolved"])))
        return {
            "case": case_id,
            "state_variant": state_variant,
            "edge_id": int(row["edge_id"]),
            "node0": int(row["node0"]),
            "node1": int(row["node1"]),
            "qf_representation": representation,
            "internal_qf_drop_V": float(internal_delta),
            "absolute_export_qf_drop_V": float(absolute_delta),
            "absolute_export_collapsed": collapsed,
            "absolute_export_loss_relative": float(absolute_loss_relative),
            "qf_drop_ulps": float(qf_drop_ulps),
            "sg_cancellation_condition": float(
                row["electron_sg_cancellation_condition"]),
            "production_stable_flux": float(production_stable),
            "long_double_reference_flux": float(
                row["electron_sg_high_precision_reference_flux"]),
            "production_subtractive_flux": float(production_direct),
            "decimal100_stable_flux": float(decimal_stable),
            "decimal100_direct_from_rounded_terms": float(decimal_direct),
            "production_vs_decimal_relative_error": float(production_relative_error),
            "subtractive_vs_decimal_stable_relative_error": float(direct_relative_error),
            "decimal_rounded_terms_vs_stable_relative_error": float(
                decimal_direct_relative_error),
            "one_ulp_flux_relative_envelope": float(ulp_relative_envelope),
            "resolved_in_production_representation": resolved,
        }


def edge_audit(contract: dict[str, Any]) -> list[dict[str, Any]]:
    precision = contract["precision"]
    selected_m21 = {int(row["edge_id"]) for row in read_csv(M21_PORTABLE)}
    rows: list[dict[str, Any]] = []
    for variant, filename in (("baseline", "baseline_sg_edges.csv"),
                              ("sentaurus_import", "sentaurus_sg_edges.csv")):
        raw = edge_map(M21_RAW / filename)
        for edge_id in sorted(selected_m21):
            rows.append(analyze_edge(
                "m21_stable_weak_inversion", variant, raw[edge_id], {}, precision))

    selected_m22 = {int(row["edge_id"]) for row in read_csv(
        M22_PORTABLE / "m22_key_state_self_consistent_edges.csv")}
    for variant in ("full", "no_hfs"):
        raw = edge_map(M22_RAW / "frozen/vg_0p05"
                       / f"state_{variant}" / f"eval_{variant}" / "sg_edges.csv")
        states = state_map(M22_RAW / "self_consistent" / variant
                           / "vg_0p05/state.csv")
        for edge_id in sorted(selected_m22):
            rows.append(analyze_edge(
                "m22_deep_off", variant, raw[edge_id], states, precision))
    return rows


def solver_control_config(variant: str, control: str, output: Path,
                          contract: dict[str, Any]) -> tuple[Path, Path, Path]:
    source_dir = M22_RAW / "self_consistent" / variant / "vg_0p05"
    config = read_json(source_dir / "config.json")
    run_dir = output / "solver_paths" / variant / control
    run_dir.mkdir(parents=True, exist_ok=True)
    curve = run_dir / "curve.csv"
    state = run_dir / "state.csv"
    config["output_csv"] = str(curve.resolve())
    config["log_file"] = str((run_dir / "run.log").resolve())
    config["sweep"]["write_state_file"] = str(state.resolve())
    if control == "tight_tolerance":
        config["solver"].update(contract["tight_tolerance"])
    elif control == "same_state_reclosure":
        config["sweep"].update({
            "initial_state_file": str((source_dir / "state.csv").resolve()),
            "start": 0.05,
            "stop": 0.05,
            "bias_points": [0.05],
        })
    elif control == "reverse_from_vg_0p1":
        reverse_state = (M22_RAW / "self_consistent" / variant
                         / "vg_0p1/state.csv")
        config["sweep"].update({
            "initial_state_file": str(reverse_state.resolve()),
            "start": 0.1,
            "stop": 0.05,
            "step": -0.05,
            "bias_points": [0.1, 0.05],
        })
    elif control != "baseline_repeat":
        raise ValueError(f"unknown solver control: {control}")
    config["simplemos_m23"] = {
        "variant": variant,
        "control": control,
        "purpose": "numerical_resolvability",
    }
    path = run_dir / "config.json"
    write_json(path, config)
    return path, curve, state


def run_solver_control(task: tuple[str, str], output: Path, runner: Path,
                       contract: dict[str, Any]) -> dict[str, Any]:
    variant, control = task
    config, curve, state = solver_control_config(
        variant, control, output, contract)
    m10.execute_runner(config, runner)
    target = min(read_csv(curve), key=lambda row: abs(float(row["bias_V"]) - 0.05))
    if not math.isclose(float(target["bias_V"]), 0.05, abs_tol=1.0e-12):
        raise RuntimeError(f"{variant}/{control} did not reach Vg=0.05 V")
    return {
        "variant": variant,
        "control": control,
        "current_A_per_um": abs(float(target["current_total_A_per_um"])),
        "state": state,
        "config": config,
    }


def solver_path_audit(contract: dict[str, Any], output: Path, runner: Path
                      ) -> list[dict[str, Any]]:
    tasks = [(variant, control) for variant in ("full", "no_hfs")
             for control in contract["solver_path_controls"]]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(
            lambda task: run_solver_control(task, output, runner, contract), tasks))
    key_nodes = {int(row["node0"]) for row in read_csv(
        M22_PORTABLE / "m22_key_state_self_consistent_edges.csv")}
    key_nodes.update(int(row["node1"]) for row in read_csv(
        M22_PORTABLE / "m22_key_state_self_consistent_edges.csv"))
    rows: list[dict[str, Any]] = []
    for result in results:
        variant = result["variant"]
        baseline_curve = read_csv(M22_RAW / "self_consistent" / variant
                                  / "vg_0p05/curve.csv")
        baseline_target = min(
            baseline_curve, key=lambda row: abs(float(row["bias_V"]) - 0.05))
        baseline_current = abs(float(baseline_target["current_total_A_per_um"]))
        baseline_state = state_map(M22_RAW / "self_consistent" / variant
                                   / "vg_0p05/state.csv")
        candidate_state = state_map(Path(result["state"]))
        maximum_increment_delta = max(
            abs(float(candidate_state[node]["electron_qf_increment_V"])
                - float(baseline_state[node]["electron_qf_increment_V"]))
            for node in key_nodes)
        current = float(result["current_A_per_um"])
        rows.append({
            "variant": variant,
            "control": result["control"],
            "baseline_current_A_per_um": baseline_current,
            "control_current_A_per_um": current,
            "current_log10_change_dex": math.log10(current / baseline_current),
            "current_relative_change": current / baseline_current - 1.0,
            "maximum_key_node_qf_increment_change_V": maximum_increment_delta,
            "config": portable(Path(result["config"])),
            "state": portable(Path(result["state"])),
        })
    return rows


def summarize_edges(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    keys = sorted({(row["case"], row["state_variant"]) for row in rows})
    for case, variant in keys:
        group = [row for row in rows
                 if row["case"] == case and row["state_variant"] == variant]
        summaries.append({
            "case": case,
            "state_variant": variant,
            "edge_count": len(group),
            "maximum_sg_cancellation_condition": max(
                row["sg_cancellation_condition"] for row in group),
            "maximum_production_vs_decimal_relative_error": max(
                row["production_vs_decimal_relative_error"] for row in group),
            "maximum_subtractive_vs_decimal_stable_relative_error": max(
                row["subtractive_vs_decimal_stable_relative_error"] for row in group),
            "maximum_one_ulp_flux_relative_envelope": max(
                row["one_ulp_flux_relative_envelope"] for row in group),
            "minimum_qf_drop_ulps": min(row["qf_drop_ulps"] for row in group),
            "absolute_export_collapsed_edge_count": sum(
                int(row["absolute_export_collapsed"]) for row in group),
            "resolved_edge_count": sum(
                int(row["resolved_in_production_representation"]) for row in group),
        })
    return summaries


def analyze(contract: dict[str, Any], output: Path, portable_dir: Path,
            runner: Path) -> dict[str, Any]:
    edges = edge_audit(contract)
    solver_paths = solver_path_audit(contract, output, runner)
    summaries = summarize_edges(edges)
    write_csv(portable_dir / "m23_edge_resolvability.csv", edges)
    write_csv(portable_dir / "m23_state_summary.csv", summaries)
    write_csv(portable_dir / "m23_solver_path_reproducibility.csv", solver_paths)
    maximum_formula_error = max(
        row["production_vs_decimal_relative_error"] for row in edges)
    maximum_ulp_envelope = max(
        row["one_ulp_flux_relative_envelope"] for row in edges)
    maximum_condition = max(row["sg_cancellation_condition"] for row in edges)
    collapsed_edges = sum(int(row["absolute_export_collapsed"]) for row in edges)
    all_resolved = all(row["resolved_in_production_representation"] for row in edges)
    solver_spread = {}
    for variant in ("full", "no_hfs"):
        currents = [row["control_current_A_per_um"] for row in solver_paths
                    if row["variant"] == variant]
        solver_spread[variant] = {
            "minimum_current_A_per_um": min(currents),
            "maximum_current_A_per_um": max(currents),
            "spread_dex": math.log10(max(currents) / min(currents)),
            "maximum_absolute_relative_change_from_baseline": max(
                abs(row["current_relative_change"]) for row in solver_paths
                if row["variant"] == variant),
        }
    report = {
        "schema": "vela.simplemos.sdevice.m23_numerical_resolvability_report.v1",
        "status": "complete",
        "execution": {
            "edge_evaluation_count": len(edges),
            "solver_control_count": len(solver_paths),
            "new_sentaurus_execution": False,
            "default_model_changed": False,
            "decimal_digits": contract["precision"]["decimal_digits"],
        },
        "findings": {
            "maximum_raw_subtractive_sg_condition": maximum_condition,
            "maximum_production_vs_decimal_relative_error": maximum_formula_error,
            "maximum_one_ulp_flux_relative_envelope": maximum_ulp_envelope,
            "absolute_export_collapsed_edge_count": collapsed_edges,
            "all_edges_resolved_in_production_representation": all_resolved,
            "solver_path_spread": solver_spread,
        },
        "state_summaries": summaries,
        "conclusions": {
            "supported": [
                "The production contact-basin QF increments must be distinguished from rounded absolute QF exports.",
                "A large raw SG subtraction condition number does not imply equal loss in the factorized production flux.",
                "The Decimal-100 replay and endpoint-ULP envelope quantify formula evaluation and stored-state resolution separately.",
            ],
            "not_supported": [
                "Treating the raw subtractive SG condition number times machine epsilon as the production current error bound.",
                "Recovering missing state digits by applying multiprecision only after an absolute binary64 QF snapshot was written.",
                "Changing an HFS or solver default from this diagnostic.",
            ],
        },
        "artifacts": {
            "edge_csv": portable(portable_dir / "m23_edge_resolvability.csv"),
            "state_summary_csv": portable(portable_dir / "m23_state_summary.csv"),
            "solver_path_csv": portable(
                portable_dir / "m23_solver_path_reproducibility.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    write_json(portable_dir / "m23_numerical_resolvability_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--portable-dir", type=Path, default=PORTABLE)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = read_json(args.contract.resolve())
    validate(contract)
    report = analyze(contract, args.output_dir.resolve(),
                     args.portable_dir.resolve(), args.runner.resolve())
    print(json.dumps({
        "status": report["status"],
        "edge_evaluations": report["execution"]["edge_evaluation_count"],
        "solver_controls": report["execution"]["solver_control_count"],
        "all_edges_resolved": report["findings"][
            "all_edges_resolved_in_production_representation"],
        "maximum_formula_error": report["findings"][
            "maximum_production_vs_decimal_relative_error"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
