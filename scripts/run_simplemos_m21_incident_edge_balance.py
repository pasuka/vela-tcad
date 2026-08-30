#!/usr/bin/env python3
"""Run the SimpleMOS M21 production-SG incident-edge audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
M14 = ROOT / "m14_adjoint_attribution"
M18 = ROOT / "m18_qf_edge_localization"
M20 = REPO / "reference_tcad/simplemos_sentaurus2022/drain_adjacent_row_audit"
OUTPUT = ROOT / "m21_incident_edge_balance"
PORTABLE = REPO / "reference_tcad/simplemos_sentaurus2022/incident_edge_balance"
RUNNER = REPO / "build-release/vela_example_runner.exe"
KEY_STATE = "n21_vd_0p05_vg_0p8"
FACTORS = ("mobility", "sg_secant_conductance", "qf_log_imbalance")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def run_config(config: Path, runner: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run([str(runner), "--config", str(config)],
                               cwd=REPO, text=True, capture_output=True,
                               env=env, check=False)
    config.with_suffix(".stdout.txt").write_text(completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"M21 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def shapley_values(game: dict[int, float]) -> dict[str, float]:
    count = len(FACTORS)
    denominator = math.factorial(count)
    result: dict[str, float] = {}
    for index, factor in enumerate(FACTORS):
        others = [item for item in range(count) if item != index]
        value = 0.0
        for size in range(count):
            weight = (math.factorial(size) * math.factorial(count - size - 1)
                      / denominator)
            for subset in combinations(others, size):
                mask = sum(1 << item for item in subset)
                value += weight * (game[mask | (1 << index)] - game[mask])
        result[factor] = value
    return result


def f(row: dict[str, str], column: str) -> float:
    return float(row[column])


def signed_for_node(value: float, row: dict[str, str], node: int) -> float:
    if int(row["node0"]) == node:
        return value
    if int(row["node1"]) == node:
        return -value
    raise ValueError(f"node {node} not on edge {row['edge_id']}")


def log_ratio(candidate: float, baseline: float) -> float:
    return math.log10(max(abs(candidate), 1.0e-300)
                      / max(abs(baseline), 1.0e-300))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable", type=Path, default=PORTABLE)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    runner = args.runner.resolve()
    output = args.output.resolve()
    portable = args.portable.resolve()
    output.mkdir(parents=True, exist_ok=True)
    portable.mkdir(parents=True, exist_ok=True)

    m20_states = read_csv(M20 / "m20_state_summary.csv")
    m20_deltas = {(row["state"], int(row["node_id"])): row
                  for row in read_csv(M20 / "m20_node_delta.csv")}
    cases: list[dict[str, Any]] = []
    jobs: list[tuple[Path, str]] = []
    all_transport_cell_vector = True
    for state_row in m20_states:
        state = state_row["state"]
        target_nodes = {int(item) for item in state_row["adjacent_nodes"].split(";")}
        case_dir = output / state
        case_dir.mkdir(parents=True, exist_ok=True)
        base_config = read_json(M14 / state / "baseline_functional.json")
        discretization = base_config["solver"]["mobility"].get(
            "high_field_gradient_discretization", "edge_projection")
        all_transport_cell_vector &= discretization == "transport_cell_vector"
        mesh = read_json(Path(base_config["mesh_file"]))
        drain = next(item for item in mesh["contacts"]
                     if item["name"].lower() == "drain")
        contact_nodes = {int(node) for node in drain["node_ids"]}
        variants: dict[str, Path] = {}
        for variant, state_file in (
                ("baseline", M14 / state / "baseline_state.csv"),
                ("sentaurus", M14 / state / "field_all_state.csv")):
            csv_path = case_dir / f"{variant}_sg_edges.csv"
            config_path = case_dir / f"{variant}_sg_probe.json"
            deck = dict(base_config)
            deck.update({
                "simulation_type": "sg_edge_flux_probe",
                "state_file": str(state_file.resolve()),
                "output_csv": str(csv_path.resolve()),
                "simplemos_m21": {
                    "read_only": True,
                    "variant": variant,
                    "mobility_drive": discretization,
                },
            })
            deck.pop("contact", None)
            deck.pop("residual_output_csv", None)
            write_json(config_path, deck)
            variants[variant] = csv_path
            if (args.force or not csv_path.is_file()
                    or not config_path.with_suffix(".stdout.txt").is_file()):
                jobs.append((config_path, f"{state}:{variant}"))
        cases.append({
            "state": state,
            "state_row": state_row,
            "target_nodes": target_nodes,
            "contact_nodes": contact_nodes,
            "variants": variants,
            "discretization": discretization,
        })

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, runner): label
                   for config, label in jobs}
        for future in as_completed(futures):
            label = futures[future]
            status = future.result()
            if int(status["edge_count"]) <= 0:
                raise AssertionError(f"empty SG probe: {label}")
            print(f"completed {label}", flush=True)

    comparison_rows: list[dict[str, Any]] = []
    node_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    max_stable_closure = 0.0
    max_factor_closure = 0.0
    max_node_closure = 0.0
    max_geometry_difference = 0.0
    for case in cases:
        state = case["state"]
        baseline = {int(row["edge_id"]): row
                    for row in read_csv(case["variants"]["baseline"])}
        sentaurus = {int(row["edge_id"]): row
                     for row in read_csv(case["variants"]["sentaurus"])}
        selected = sorted(edge_id for edge_id, row in baseline.items()
                          if case["target_nodes"] & {
                              int(row["node0"]), int(row["node1"])})
        raw_factors = read_csv(M18 / state / "sg_secant_factor_edges.csv")
        factor_games: dict[int, dict[int, float]] = {}
        for row in raw_factors:
            factor_games.setdefault(int(row["edge_id"]), {})[
                int(row["mask"])] = f(row, "electron_flux")
        per_edge_factors: dict[int, dict[str, float]] = {}
        state_factor_abs = {factor: 0.0 for factor in FACTORS}
        state_factor_closure = 0.0
        for edge_id in selected:
            game = factor_games[edge_id]
            values = shapley_values(game)
            closure = sum(values.values()) - (game[7] - game[0])
            max_factor_closure = max(max_factor_closure, abs(closure))
            state_factor_closure = max(state_factor_closure, abs(closure))
            per_edge_factors[edge_id] = values
            for factor, value in values.items():
                state_factor_abs[factor] += abs(value)
        state_factor_abs_total = sum(state_factor_abs.values())

        state_stable_closure = 0.0
        for edge_id in selected:
            base = baseline[edge_id]
            sent = sentaurus[edge_id]
            if base["electron_sg_boltzmann_decomposition_available"] != "1" \
                    or sent["electron_sg_boltzmann_decomposition_available"] != "1":
                raise AssertionError(f"{state}:{edge_id}: Boltzmann decomposition absent")
            for row in (base, sent):
                closure = abs(f(row, "electron_sg_stable_flux")
                              - f(row, "electron_flux"))
                max_stable_closure = max(max_stable_closure, closure)
                state_stable_closure = max(state_stable_closure, closure)
            length_difference = f(sent, "length_m") - f(base, "length_m")
            couple_difference = f(sent, "couple_m") - f(base, "couple_m")
            max_geometry_difference = max(
                max_geometry_difference, abs(length_difference),
                abs(couple_difference))
            endpoints = {int(base["node0"]), int(base["node1"])}
            incident_targets = sorted(case["target_nodes"] & endpoints)
            edge_class = "drain_contact_cut" if (
                (int(base["node0"]) in case["contact_nodes"])
                != (int(base["node1"]) in case["contact_nodes"])) else "internal"
            factors = per_edge_factors[edge_id]
            factor_abs = sum(abs(value) for value in factors.values())
            row: dict[str, Any] = {
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "edge_id": edge_id,
                "node0": base["node0"],
                "node1": base["node1"],
                "incident_target_nodes": ";".join(map(str, incident_targets)),
                "edge_class": edge_class,
                "length_m": base["length_m"],
                "couple_m": base["couple_m"],
                "length_difference_m": length_difference,
                "couple_difference_m": couple_difference,
                "baseline_phin0_V": base["phin0_V"],
                "baseline_phin1_V": base["phin1_V"],
                "sentaurus_phin0_V": sent["phin0_V"],
                "sentaurus_phin1_V": sent["phin1_V"],
                "baseline_phin_drop_V": f(base, "phin1_V") - f(base, "phin0_V"),
                "sentaurus_phin_drop_V": f(sent, "phin1_V") - f(sent, "phin0_V"),
                "phin_drop_difference_V": ((f(sent, "phin1_V") - f(sent, "phin0_V"))
                                           - (f(base, "phin1_V") - f(base, "phin0_V"))),
                "baseline_electron_density0_m3": base["electron_density0_m3"],
                "baseline_electron_density1_m3": base["electron_density1_m3"],
                "sentaurus_electron_density0_m3": sent["electron_density0_m3"],
                "sentaurus_electron_density1_m3": sent["electron_density1_m3"],
                "maximum_endpoint_density_change_dex": max(
                    abs(log_ratio(f(sent, "electron_density0_m3"),
                                  f(base, "electron_density0_m3"))),
                    abs(log_ratio(f(sent, "electron_density1_m3"),
                                  f(base, "electron_density1_m3")))),
                "baseline_transport_cell_vector_drive_V_m":
                    base["electron_mobility_drive_V_m"],
                "sentaurus_transport_cell_vector_drive_V_m":
                    sent["electron_mobility_drive_V_m"],
                "transport_cell_vector_drive_change_dex": log_ratio(
                    f(sent, "electron_mobility_drive_V_m"),
                    f(base, "electron_mobility_drive_V_m")),
                "baseline_mobility_m2_V_s": base["electron_mobility_m2_V_s"],
                "sentaurus_mobility_m2_V_s": sent["electron_mobility_m2_V_s"],
                "mobility_change_dex": log_ratio(
                    f(sent, "electron_mobility_m2_V_s"),
                    f(base, "electron_mobility_m2_V_s")),
                "baseline_eta": base["electron_sg_eta"],
                "sentaurus_eta": sent["electron_sg_eta"],
                "eta_difference": f(sent, "electron_sg_eta") - f(base, "electron_sg_eta"),
                "baseline_bernoulli_minus_eta": base[
                    "electron_sg_bernoulli_minus_eta"],
                "sentaurus_bernoulli_minus_eta": sent[
                    "electron_sg_bernoulli_minus_eta"],
                "baseline_bernoulli_eta": base["electron_sg_bernoulli_eta"],
                "sentaurus_bernoulli_eta": sent["electron_sg_bernoulli_eta"],
                "maximum_bernoulli_weight_change_dex": max(
                    abs(log_ratio(f(sent, "electron_sg_bernoulli_minus_eta"),
                                  f(base, "electron_sg_bernoulli_minus_eta"))),
                    abs(log_ratio(f(sent, "electron_sg_bernoulli_eta"),
                                  f(base, "electron_sg_bernoulli_eta")))),
                "baseline_left_term_m3": base["electron_sg_left_term_m3"],
                "baseline_right_term_m3": base["electron_sg_right_term_m3"],
                "sentaurus_left_term_m3": sent["electron_sg_left_term_m3"],
                "sentaurus_right_term_m3": sent["electron_sg_right_term_m3"],
                "baseline_sg_cancellation_condition": base[
                    "electron_sg_cancellation_condition"],
                "sentaurus_sg_cancellation_condition": sent[
                    "electron_sg_cancellation_condition"],
                "baseline_electron_flux": base["electron_flux"],
                "sentaurus_electron_flux": sent["electron_flux"],
                "electron_flux_difference": f(sent, "electron_flux") - f(base, "electron_flux"),
                **{f"{factor}_flux_contribution": value
                   for factor, value in factors.items()},
                **{f"{factor}_absolute_share": abs(value) / max(factor_abs, 1.0e-300)
                   for factor, value in factors.items()},
                "factor_flux_closure": sum(factors.values())
                    - (f(sent, "electron_flux") - f(base, "electron_flux")),
            }
            comparison_rows.append(row)

        focus_nodes = sorted(case["target_nodes"])
        state_node_closure = 0.0
        for node in focus_nodes:
            incident = [row for row in comparison_rows
                        if row["state"] == state and
                        node in (int(row["node0"]), int(row["node1"]))]
            flux_delta = sum(signed_for_node(
                float(row["electron_flux_difference"]), row, node)
                for row in incident)
            factor_values = {
                factor: sum(signed_for_node(
                    float(row[f"{factor}_flux_contribution"]), row, node)
                    for row in incident)
                for factor in FACTORS
            }
            m20 = m20_deltas[(state, node)]
            m20_flux_delta = (f(m20, "delta_contact_flux")
                              + f(m20, "delta_internal_flux"))
            closure = flux_delta - m20_flux_delta
            max_node_closure = max(max_node_closure, abs(closure))
            state_node_closure = max(state_node_closure, abs(closure))
            node_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "node_id": node,
                "incident_edge_count": len(incident),
                "production_edge_flux_delta_sum": flux_delta,
                "m20_contact_plus_internal_flux_delta": m20_flux_delta,
                "m20_flux_delta_closure": closure,
                **{f"{factor}_signed_flux_contribution": value
                   for factor, value in factor_values.items()},
                "factor_sum_flux_contribution": sum(factor_values.values()),
            })

        state_rows.append({
            "state": state,
            "device": case["state_row"]["device"],
            "drain_voltage_V": case["state_row"]["drain_voltage_V"],
            "gate_voltage_V": case["state_row"]["gate_voltage_V"],
            "selected_edge_count": len(selected),
            "target_nodes": ";".join(map(str, focus_nodes)),
            "high_field_gradient_discretization": case["discretization"],
            "maximum_stable_sg_flux_absolute_closure": state_stable_closure,
            "maximum_edge_factor_absolute_closure": state_factor_closure,
            "maximum_m20_node_flux_delta_absolute_closure": state_node_closure,
            **{f"{factor}_absolute_flux_contribution": value
               for factor, value in state_factor_abs.items()},
            **{f"{factor}_absolute_share": value /
               max(state_factor_abs_total, 1.0e-300)
               for factor, value in state_factor_abs.items()},
        })

    write_csv(portable / "m21_edge_input_comparison.csv", comparison_rows)
    write_csv(portable / "m21_node_factor_ledger.csv", node_rows)
    write_csv(portable / "m21_state_summary.csv", state_rows)
    key_edges = [row for row in comparison_rows if row["state"] == KEY_STATE
                 and ({int(row["node0"]), int(row["node1"])} & {991, 992})]
    key_edges.sort(key=lambda row: abs(float(row["electron_flux_difference"])),
                   reverse=True)
    write_csv(portable / "m21_key_state_incident_edges.csv", key_edges)
    key_nodes = [row for row in node_rows if row["state"] == KEY_STATE
                 and int(row["node_id"]) in {991, 992}]
    write_csv(portable / "m21_key_state_node_factor_ledger.csv", key_nodes)

    key_factor_abs = {
        factor: sum(abs(float(row[f"{factor}_flux_contribution"]))
                    for row in key_edges)
        for factor in FACTORS
    }
    key_factor_total = sum(key_factor_abs.values())
    report = {
        "schema": "vela.simplemos.sdevice.m21_incident_edge_balance_report.v1",
        "execution": {
            "state_count": len(cases),
            "variant_count": 2,
            "sg_probe_count": len(cases) * 2,
            "selected_edge_rows": len(comparison_rows),
            "all_states_use_transport_cell_vector": all_transport_cell_vector,
            "new_sentaurus_execution": False,
            "diagnostic_cpp_changed": True,
            "default_model_changed": False,
        },
        "closure": {
            "maximum_stable_sg_flux_absolute_closure": max_stable_closure,
            "maximum_edge_factor_absolute_closure": max_factor_closure,
            "maximum_m20_node_flux_delta_absolute_closure": max_node_closure,
            "maximum_geometry_difference": max_geometry_difference,
        },
        "findings": {
            "key_state": KEY_STATE,
            "key_incident_edge_count": len(key_edges),
            "key_factor_absolute_sums": key_factor_abs,
            "key_factor_absolute_shares": {
                factor: value / max(key_factor_total, 1.0e-300)
                for factor, value in key_factor_abs.items()
            },
            "all_state_qf_absolute_share_range": [
                min(float(row["qf_log_imbalance_absolute_share"])
                    for row in state_rows),
                max(float(row["qf_log_imbalance_absolute_share"])
                    for row in state_rows),
            ],
            "key_nodes": key_nodes,
            "leading_key_edges": key_edges[:8],
        },
        "interpretation_limits": read_json(
            REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m21_incident_edge_balance_contract_v1.json"
        )["interpretation_limits"],
    }
    write_json(portable / "m21_incident_edge_balance_report.json", report)
    print(json.dumps({
        "status": "complete",
        "state_count": len(cases),
        "probe_count": len(cases) * 2,
        "selected_edge_rows": len(comparison_rows),
        "key_factor_absolute_shares": report["findings"][
            "key_factor_absolute_shares"],
        "maximum_stable_sg_flux_absolute_closure": max_stable_closure,
        "maximum_m20_node_flux_delta_absolute_closure": max_node_closure,
    }))


if __name__ == "__main__":
    main()
