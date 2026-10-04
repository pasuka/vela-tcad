#!/usr/bin/env python3
"""Run the SimpleMOS M19 drain-contact cut audit.

M19 keeps the converged Vela and exported Sentaurus states frozen.  It audits
the native drain contact cut, evaluates an exact three-factor endpoint-state
game with the production SG probe, and relates the frozen result to the M18
self-consistent adjoint ledger.
"""

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
BUILD_ROOT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
M10 = BUILD_ROOT / "m10_fixed_state_replay"
M14 = BUILD_ROOT / "m14_adjoint_attribution"
M18 = BUILD_ROOT / "m18_qf_edge_localization"
M18_PORTABLE = REPO / "reference_tcad/simplemos_sentaurus2022/qf_edge_localization"
M8_COMPARISON = REPO / "reference_tcad/simplemos_sentaurus2022/original_physics/comparisons"
OUTPUT = BUILD_ROOT / "m19_drain_cut_audit"
PORTABLE = REPO / "reference_tcad/simplemos_sentaurus2022/drain_cut_audit"
RUNNER = REPO / "build-release/vela_example_runner.exe"
FACTORS = (
    "drain_contact_endpoint_state",
    "drain_adjacent_interior_phin",
    "drain_adjacent_interior_other_state",
)
STATE_COLUMNS = ("psi", "phin", "phip", "electrons_m3", "holes_m3")
OTHER_STATE_COLUMNS = ("psi", "phip", "electrons_m3", "holes_m3")
KEY_STATE = "n21_vd_0p05_vg_0p8"
Q = 1.602176634e-19


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
        raise RuntimeError(f"M19 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def write_hybrid(path: Path, baseline: list[dict[str, str]],
                 replacement: dict[int, dict[str, str]],
                 contact_nodes: set[int], adjacent_nodes: set[int],
                 mask: int) -> None:
    rows: list[dict[str, Any]] = []
    columns = ("node_id",) + STATE_COLUMNS
    for source in baseline:
        node = int(source["node_id"])
        row = {column: source[column] for column in columns}
        if mask & 1 and node in contact_nodes:
            for field in STATE_COLUMNS:
                row[field] = replacement[node][field]
        if mask & 2 and node in adjacent_nodes:
            row["phin"] = replacement[node]["phin"]
        if mask & 4 and node in adjacent_nodes:
            for field in OTHER_STATE_COLUMNS:
                row[field] = replacement[node][field]
        rows.append(row)
    write_csv(path, rows)


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


def drain_support(mesh: dict[str, Any], state: str) -> tuple[set[int], set[tuple[int, int]]]:
    contacts = [item for item in mesh["contacts"] if item["name"].lower() == "drain"]
    if len(contacts) != 1:
        raise ValueError("mesh must contain exactly one drain contact")
    contact_nodes = {int(node) for node in contacts[0]["node_ids"]}
    cut: set[tuple[int, int]] = set()
    edge_path = M18 / state / "sg_secant_factor_edges.csv"
    for edge in read_csv(edge_path):
        if int(edge["mask"]) != 0:
            continue
        node0, node1 = int(edge["node0"]), int(edge["node1"])
        if (node0 in contact_nodes) != (node1 in contact_nodes):
            cut.add((min(node0, node1), max(node0, node1)))
    return contact_nodes, cut


def cut_fluxes(path: Path, contact_nodes: set[int],
               cut_pairs: set[tuple[int, int]]) -> tuple[dict[int, dict[str, Any]], dict[str, float]]:
    edges: dict[int, dict[str, Any]] = {}
    electron = 0.0
    hole = 0.0
    for row in read_csv(path):
        node0, node1 = int(row["node0"]), int(row["node1"])
        pair = (min(node0, node1), max(node0, node1))
        if pair not in cut_pairs:
            continue
        node0_contact = node0 in contact_nodes
        node1_contact = node1 in contact_nodes
        if node0_contact == node1_contact:
            raise AssertionError("drain cut edge classification drift")
        outward = 1.0 if node0_contact else -1.0
        electron_current = (-Q * outward
                            * float(row["electron_particle_line_flux_per_m_s"])
                            * 1.0e-6)
        hole_current = (-Q * outward
                        * float(row["hole_particle_line_flux_per_m_s"])
                        * 1.0e-6)
        electron += electron_current
        hole += hole_current
        edges[int(row["edge_id"])] = {
            **row,
            "outward_sign": outward,
            "electron_current_A_per_um": electron_current,
            "hole_current_A_per_um": hole_current,
            "total_current_A_per_um": electron_current - hole_current,
        }
    if len(edges) != len(cut_pairs):
        raise AssertionError(f"missing SG drain-cut edges in {path}")
    return edges, {
        "electron_A_per_um": electron,
        "hole_A_per_um": hole,
        "total_A_per_um": electron - hole,
    }


def curve_point(device: str, drain: str, gate: float) -> dict[str, str]:
    path = M8_COMPARISON / f"{device}_vd_{drain}_comparison.csv"
    rows = [row for row in read_csv(path)
            if math.isclose(float(row["gate_voltage_V"]), gate,
                            rel_tol=0.0, abs_tol=1.0e-12)]
    if len(rows) != 1:
        raise ValueError(f"missing direct M8 curve point: {device} {drain} {gate}")
    return rows[0]


def log_ratio(candidate: float, reference: float) -> float:
    return math.log10(max(abs(candidate), 1.0e-300)
                      / max(abs(reference), 1.0e-300))


def state_parts(state: str) -> tuple[str, str, float]:
    device, _, drain, _, gate = state.split("_", 4)
    return device, drain, float(gate.replace("p", "."))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable", type=Path, default=PORTABLE)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    output = args.output.resolve()
    portable = args.portable.resolve()
    runner = args.runner.resolve()
    output.mkdir(parents=True, exist_ok=True)
    portable.mkdir(parents=True, exist_ok=True)

    m10_report = read_json(M10 / "fixed_state_replay_report.json")
    m10_cases = {item["state"]: item for item in m10_report["cases"]}
    m18_states = read_csv(M18_PORTABLE / "m18_state_summary.csv")
    state_names = [row["state"] for row in m18_states]
    m18_buckets = {
        (row["state"], row["bucket"], row["factor"]): row
        for row in read_csv(M18_PORTABLE / "m18_bucket_factor_contributions.csv")
    }
    m18_ledger_rows = read_csv(M18_PORTABLE / "m18_feedback_ledger.csv")
    m18_ledger = {(row["state"], row["component"]): row
                  for row in m18_ledger_rows}

    cases: list[dict[str, Any]] = []
    jobs: list[tuple[Path, Path, str]] = []
    for state in state_names:
        case_dir = output / state
        case_dir.mkdir(parents=True, exist_ok=True)
        base_config = read_json(M14 / state / "baseline_functional.json")
        mesh_path = Path(base_config["mesh_file"])
        mesh = read_json(mesh_path)
        contact_nodes, cut_pairs = drain_support(mesh, state)
        adjacent_nodes = {
            node for pair in cut_pairs for node in pair if node not in contact_nodes}
        baseline_path = M14 / state / "baseline_state.csv"
        replacement_path = M10 / "replay" / state / "sentaurus_state_for_vela.csv"
        baseline = read_csv(baseline_path)
        replacement_rows = read_csv(replacement_path)
        replacement = {int(row["node_id"]): row for row in replacement_rows}
        variants: dict[int | str, Path] = {}
        for mask in range(8):
            state_path = case_dir / f"cut_mask_{mask}_state.csv"
            write_hybrid(state_path, baseline, replacement,
                         contact_nodes, adjacent_nodes, mask)
            variants[mask] = state_path
        variants["full_sentaurus"] = replacement_path
        probes: dict[int | str, Path] = {}
        for variant, state_path in variants.items():
            label = f"cut_mask_{variant}" if isinstance(variant, int) else str(variant)
            csv_path = case_dir / f"{label}_sg_edges.csv"
            config_path = case_dir / f"{label}_sg_probe.json"
            config = dict(base_config)
            config.update({
                "simulation_type": "sg_edge_flux_probe",
                "state_file": str(state_path.resolve()),
                "output_csv": str(csv_path.resolve()),
                "simplemos_m19": {"read_only": True, "variant": label},
            })
            config.pop("contact", None)
            config.pop("residual_output_csv", None)
            write_json(config_path, config)
            probes[variant] = csv_path
            if (args.force or not csv_path.is_file()
                    or not config_path.with_suffix(".stdout.txt").is_file()):
                jobs.append((config_path, csv_path, f"{state}:{label}"))
        cases.append({
            "state": state,
            "mesh": mesh,
            "mesh_path": mesh_path,
            "contact_nodes": contact_nodes,
            "cut_pairs": cut_pairs,
            "adjacent_nodes": adjacent_nodes,
            "baseline": {int(row["node_id"]): row for row in baseline},
            "replacement": replacement,
            "probes": probes,
            "drain_bias": next(float(item["bias"]) for item in base_config["contacts"]
                               if item["name"].lower() == "drain"),
        })

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, runner): label
                   for config, _expected, label in jobs}
        for future in as_completed(futures):
            label = futures[future]
            status = future.result()
            if int(status["edge_count"]) <= 0:
                raise AssertionError(f"empty M19 SG probe: {label}")
            print(f"completed {label}", flush=True)

    state_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    factor_rows: list[dict[str, Any]] = []
    maximum_reconstruction_relative = 0.0
    maximum_strong_reconstruction_relative = 0.0
    maximum_shapley_closure = 0.0
    maximum_contact_bias_error = 0.0
    key_edges: list[dict[str, Any]] = []
    for case in cases:
        state = case["state"]
        device, drain_tag, gate = state_parts(state)
        drain = case["drain_bias"]
        base = case["baseline"]
        sent = case["replacement"]
        contact_nodes = case["contact_nodes"]
        cut_pairs = case["cut_pairs"]
        cut_by_variant: dict[int | str, dict[int, dict[str, Any]]] = {}
        current_by_variant: dict[int | str, dict[str, float]] = {}
        for variant, path in case["probes"].items():
            cut_by_variant[variant], current_by_variant[variant] = cut_fluxes(
                path, contact_nodes, cut_pairs)

        curve = curve_point(device, drain_tag, gate)
        vela_curve_current = float(curve["vela_current_A_per_um"])
        sentaurus_current = float(curve["sentaurus_current_A_per_um"])
        baseline_current = current_by_variant[0]["total_A_per_um"]
        full_current = current_by_variant["full_sentaurus"]["total_A_per_um"]
        reconstruction_relative = abs(baseline_current - vela_curve_current) / max(
            abs(baseline_current), abs(vela_curve_current), 1.0e-300)
        maximum_reconstruction_relative = max(
            maximum_reconstruction_relative, reconstruction_relative)
        if gate >= 0.8:
            maximum_strong_reconstruction_relative = max(
                maximum_strong_reconstruction_relative, reconstruction_relative)

        game = {mask: current_by_variant[mask]["total_A_per_um"] for mask in range(8)}
        factor_values = shapley_values(game)
        closure = sum(factor_values.values()) - (game[7] - game[0])
        maximum_shapley_closure = max(maximum_shapley_closure, abs(closure))
        for factor, value in factor_values.items():
            factor_rows.append({
                "state": state,
                "device": device,
                "drain_voltage_V": drain,
                "gate_voltage_V": gate,
                "factor": factor,
                "direct_current_change_A_per_um": value,
                "relative_to_vela_current": value / baseline_current,
                "absolute_factor_share": abs(value) / max(
                    sum(abs(item) for item in factor_values.values()), 1.0e-300),
            })

        contact_bias_errors_vela = [abs(float(base[node]["phin"]) - drain)
                                    for node in contact_nodes]
        contact_bias_errors_sent = [abs(float(sent[node]["phin"]) - drain)
                                    for node in contact_nodes]
        contact_phin_differences = [abs(float(base[node]["phin"])
                                        - float(sent[node]["phin"]))
                                    for node in contact_nodes]
        maximum_contact_bias_error = max(
            maximum_contact_bias_error,
            *contact_bias_errors_vela, *contact_bias_errors_sent)

        comparison = {
            int(row["edge_id"]): row for row in read_csv(
                M10 / "replay" / state / "edge_replay_comparison.csv")}
        edge_factor_totals = {factor: 0.0 for factor in FACTORS}
        for edge_id in sorted(cut_by_variant[0]):
            edge_game = {mask: cut_by_variant[mask][edge_id]["total_current_A_per_um"]
                         for mask in range(8)}
            edge_factors = shapley_values(edge_game)
            for factor, value in edge_factors.items():
                edge_factor_totals[factor] += value
            base_edge = cut_by_variant[0][edge_id]
            sent_edge = cut_by_variant["full_sentaurus"][edge_id]
            node0, node1 = int(base_edge["node0"]), int(base_edge["node1"])
            length_m = float(base_edge["length_m"])
            vela_drop = float(base[node1]["phin"]) - float(base[node0]["phin"])
            sent_drop = float(sent[node1]["phin"]) - float(sent[node0]["phin"])
            replay = comparison[edge_id]
            row = {
                "state": state,
                "device": device,
                "drain_voltage_V": drain,
                "gate_voltage_V": gate,
                "edge_id": edge_id,
                "node0": node0,
                "node1": node1,
                "node0_is_drain_contact": int(node0 in contact_nodes),
                "node1_is_drain_contact": int(node1 in contact_nodes),
                "length_m": length_m,
                "midpoint_x_um": 0.5 * (float(base_edge["x0"]) + float(base_edge["x1"])) * 1.0e6,
                "midpoint_y_um": 0.5 * (float(base_edge["y0"]) + float(base_edge["y1"])) * 1.0e6,
                "vela_phin_drop_V": vela_drop,
                "sentaurus_phin_drop_V": sent_drop,
                "vela_minus_sentaurus_phin_drop_V": vela_drop - sent_drop,
                "vela_endpoint_secant_abs_V_cm": abs(vela_drop) / length_m / 100.0,
                "sentaurus_endpoint_secant_abs_V_cm": abs(sent_drop) / length_m / 100.0,
                "sentaurus_exported_gradqf_mean_V_cm": float(
                    replay["sentaurus_eGradQF_mean_magnitude_V_cm"]),
                "vela_transport_cell_vector_on_sentaurus_state_V_cm": float(
                    replay["vela_eGradQF_drive_V_cm"]),
                "vela_current_A_per_um": base_edge["total_current_A_per_um"],
                "sentaurus_state_vela_sg_current_A_per_um": sent_edge[
                    "total_current_A_per_um"],
                **{f"{factor}_direct_A_per_um": value
                   for factor, value in edge_factors.items()},
            }
            edge_rows.append(row)
            if state == KEY_STATE:
                key_edges.append(row)
        for factor in FACTORS:
            if not math.isclose(edge_factor_totals[factor], factor_values[factor],
                                rel_tol=1.0e-9, abs_tol=1.0e-24):
                raise AssertionError(f"edge/current Shapley mismatch: {state} {factor}")

        m18_qf_cut = float(m18_buckets[
            (state, "drain_contact_cut", "qf_log_imbalance")][
                "relative_to_baseline_current"])
        ledger = {component: float(m18_ledger[(state, component)][
            "relative_to_baseline_current"])
            for component in ("mobility", "sg_secant_conductance",
                              "qf_log_imbalance", "poisson", "srh",
                              "boundary_gauge", "other_carrier")}
        full_feedback = sum(ledger.values())
        denominator = baseline_current - sentaurus_current
        cut_gap_fraction = ((baseline_current - game[7]) / denominator
                            if abs(denominator) > 1.0e-300 else 0.0)
        full_gap_fraction = ((baseline_current - full_current) / denominator
                             if abs(denominator) > 1.0e-300 else 0.0)
        state_rows.append({
            "state": state,
            "device": device,
            "drain_voltage_V": drain,
            "gate_voltage_V": gate,
            "drain_contact_node_count": len(contact_nodes),
            "drain_cut_edge_count": len(cut_pairs),
            "drain_adjacent_interior_node_count": len(case["adjacent_nodes"]),
            "max_vela_contact_phin_bias_error_V": max(contact_bias_errors_vela),
            "max_sentaurus_contact_phin_bias_error_V": max(contact_bias_errors_sent),
            "max_contact_phin_difference_V": max(contact_phin_differences),
            "vela_self_consistent_current_A_per_um": baseline_current,
            "m8_vela_curve_current_A_per_um": vela_curve_current,
            "sentaurus_current_A_per_um": sentaurus_current,
            "terminal_signed_log10_error_dex": log_ratio(
                baseline_current, sentaurus_current),
            "sentaurus_state_vela_sg_current_A_per_um": full_current,
            "sentaurus_state_vela_sg_error_dex": log_ratio(
                full_current, sentaurus_current),
            "cut_endpoint_state_current_A_per_um": game[7],
            "cut_endpoint_state_gap_closed_fraction": cut_gap_fraction,
            "full_sentaurus_state_gap_closed_fraction": full_gap_fraction,
            "sg_cut_reconstruction_relative_error": reconstruction_relative,
            "cut_game_shapley_closure_A_per_um": closure,
            "m18_drain_cut_qf_response_relative_to_Id": m18_qf_cut,
            "m18_full_feedback_relative_to_Id": full_feedback,
            **{f"frozen_{factor}_relative_to_Id": value / baseline_current
               for factor, value in factor_values.items()},
        })

    write_csv(portable / "m19_state_summary.csv", state_rows)
    write_csv(portable / "m19_drain_cut_edges.csv", edge_rows)
    write_csv(portable / "m19_frozen_factor_contributions.csv", factor_rows)
    write_csv(portable / "m19_key_state_drain_cut_edges.csv", key_edges)

    key = next(row for row in state_rows if row["state"] == KEY_STATE)
    key_factors = {row["factor"]: row for row in factor_rows
                   if row["state"] == KEY_STATE}
    ranked_key_edges = sorted(
        key_edges,
        key=lambda row: abs(float(row["drain_adjacent_interior_phin_direct_A_per_um"])),
        reverse=True)
    strong_states = [row for row in state_rows if float(row["gate_voltage_V"]) >= 0.8]
    report = {
        "schema": "vela.simplemos.sdevice.m19_drain_cut_audit.report.v1",
        "status": "complete",
        "execution": {
            "state_count": len(state_rows),
            "sg_probe_count": len(state_rows) * 9,
            "probe_variants_per_state": 9,
            "probes_executed_this_run": len(jobs),
            "new_sentaurus_execution": False,
            "cpp_changed": False,
            "default_model_changed": False,
        },
        "closure": {
            "maximum_contact_phin_bias_error_V": maximum_contact_bias_error,
            "maximum_all_state_sg_cut_current_reconstruction_relative_error_guarded":
                maximum_reconstruction_relative,
            "maximum_strong_state_sg_cut_current_reconstruction_relative_error":
                maximum_strong_reconstruction_relative,
            "maximum_shapley_closure_A_per_um": maximum_shapley_closure,
        },
        "findings": {
            "contact_boundary_phin_matches_bias_all_states":
                maximum_contact_bias_error <= 1.0e-12,
            "strong_state_max_abs_sentaurus_state_vela_sg_error_dex": max(
                abs(float(row["sentaurus_state_vela_sg_error_dex"]))
                for row in strong_states),
            "vg_0p8_full_state_gap_closed_fraction_range": [
                min(float(row["full_sentaurus_state_gap_closed_fraction"])
                    for row in state_rows
                    if math.isclose(float(row["gate_voltage_V"]), 0.8)),
                max(float(row["full_sentaurus_state_gap_closed_fraction"])
                    for row in state_rows
                    if math.isclose(float(row["gate_voltage_V"]), 0.8)),
            ],
            "key_state": {
                **key,
                "frozen_factor_absolute_shares": {
                    factor: float(key_factors[factor]["absolute_factor_share"])
                    for factor in FACTORS
                },
                "top_adjacent_phin_edges": [
                    {
                        "edge_id": int(row["edge_id"]),
                        "direct_A_per_um": float(
                            row["drain_adjacent_interior_phin_direct_A_per_um"]),
                        "vela_phin_drop_V": float(row["vela_phin_drop_V"]),
                        "sentaurus_phin_drop_V": float(row["sentaurus_phin_drop_V"]),
                    }
                    for row in ranked_key_edges[:3]
                ],
            },
        },
        "interpretation_limits": [
            "The imported Sentaurus state uses exact original TDR node IDs without spatial interpolation.",
            "The exported Sentaurus GradQuasiFermi vector is diagnostic and is not treated as a disclosed private SG face formula.",
            "Frozen endpoint Shapley values and M18 adjoint responses are different mathematical layers and are reported separately.",
            "SG cut-current closure is accepted on Vg>=0.8 V states; sub-fA and weak-current mapped-state cancellation is not headline attribution evidence.",
        ],
    }
    write_json(portable / "m19_drain_cut_audit_report.json", report)
    print(json.dumps(report["closure"], indent=2), flush=True)


if __name__ == "__main__":
    main()
