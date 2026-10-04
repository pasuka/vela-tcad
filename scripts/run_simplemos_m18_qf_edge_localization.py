#!/usr/bin/env python3
"""Run M18 edge localization and feedback ledger for SimpleMOS QF response."""

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
M14_BUILD = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m14_adjoint_attribution")
M15 = REPO / "reference_tcad/simplemos_sentaurus2022/operator_adjoint"
M17 = REPO / "reference_tcad/simplemos_sentaurus2022/sg_secant"
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m18_qf_edge_localization")
PORTABLE = REPO / "reference_tcad/simplemos_sentaurus2022/qf_edge_localization"
RUNNER = REPO / "build-release/vela_example_runner.exe"
FACTORS = ("mobility", "sg_secant_conductance", "qf_log_imbalance")
BUCKETS = ("source_contact_cut", "drain_contact_cut",
           "substrate_contact_cut", "gate_contact_cut", "internal_source",
           "internal_channel", "internal_drain", "internal_body")
KEY_STATE = "n21_vd_0p05_vg_0p8"


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
    config.with_suffix(".stdout.txt").write_text(completed.stdout,
                                                   encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr,
                                                   encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"M18 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def shapley_values(game: dict[int, float]) -> dict[str, float]:
    count = len(FACTORS)
    denominator = math.factorial(count)
    values: dict[str, float] = {}
    for index, factor in enumerate(FACTORS):
        contribution = 0.0
        others = [item for item in range(count) if item != index]
        for size in range(count):
            weight = (math.factorial(size) * math.factorial(count - size - 1)
                      / denominator)
            for subset in combinations(others, size):
                mask = sum(1 << item for item in subset)
                contribution += weight * (game[mask | (1 << index)] - game[mask])
        values[factor] = contribution
    return values


def contact_names(value: str) -> set[str]:
    return {item.strip().lower() for item in value.split(";") if item.strip()}


def edge_bucket(row: dict[str, str]) -> str:
    endpoint0 = contact_names(row["node0_contacts"])
    endpoint1 = contact_names(row["node1_contacts"])
    for name in ("source", "drain", "substrate", "gate"):
        if (name in endpoint0) != (name in endpoint1):
            return f"{name}_contact_cut"
    x_um = float(row["midpoint_x"])
    y_um = float(row["midpoint_y"])
    if x_um > 0.25:
        return "internal_body"
    if y_um < -0.3:
        return "internal_source"
    if y_um > 0.3:
        return "internal_drain"
    return "internal_channel"


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

    m17_states = read_csv(M17 / "m17_state_summary.csv")
    m17_factors = read_csv(M17 / "m17_factor_contributions.csv")
    m17_factor_index = {(row["state"], row["factor"]): row
                        for row in m17_factors}
    m15_states = {row["state"]: row for row in read_csv(
        M15 / "m15_state_summary.csv")}
    m15_components = read_csv(M15 / "m15_component_contributions.csv")
    m15_component_index = {(row["state"], row["component"]): row
                           for row in m15_components}

    cases: list[dict[str, Any]] = []
    jobs: list[tuple[Path, Path, str]] = []
    for state_row in m17_states:
        state = state_row["state"]
        case_dir = output / state
        case_dir.mkdir(parents=True, exist_ok=True)
        source = read_json(M14_BUILD / state / "baseline_functional.json")
        node_csv = case_dir / "sg_secant_factor_nodes.csv"
        edge_csv = case_dir / "sg_secant_factor_edges.csv"
        config = case_dir / "qf_edge_probe.json"
        source.update({
            "simulation_type": "electron_transport_secant_factor_probe",
            "state_file": str((M14_BUILD / state / "baseline_state.csv").resolve()),
            "replacement_state_file": str(
                (M14_BUILD / state / "field_all_state.csv").resolve()),
            "output_csv": str(node_csv.resolve()),
            "edge_output_csv": str(edge_csv.resolve()),
            "simplemos_m18": {"read_only": True, "edge_shapley": True},
        })
        source.pop("residual_output_csv", None)
        write_json(config, source)
        jobs.append((config, edge_csv, state))
        cases.append({"state": state, "state_row": state_row,
                      "edge_csv": edge_csv})

    selected = [(config, expected, label) for config, expected, label in jobs
                if args.force or not expected.exists()
                or not config.with_suffix(".stdout.txt").exists()]
    edge_counts: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, runner): label
                   for config, _expected, label in selected}
        for future in as_completed(futures):
            label = futures[future]
            status = future.result()
            if status["variant_count"] != 8:
                raise ValueError(f"M18 variant matrix incomplete for {label}")
            edge_counts[label] = int(status["edge_row_count"]) // 8
            print(f"completed {label}", flush=True)

    summary_rows: list[dict[str, Any]] = []
    bucket_rows: list[dict[str, Any]] = []
    key_edge_rows: list[dict[str, Any]] = []
    ledger_rows: list[dict[str, Any]] = []
    maximum_edge_shapley_closure = 0.0
    for case in cases:
        state = case["state"]
        state_row = case["state_row"]
        raw_edges = read_csv(case["edge_csv"])
        by_mask: dict[int, dict[int, dict[str, str]]] = {}
        for row in raw_edges:
            by_mask.setdefault(int(row["mask"]), {})[int(row["edge_id"])] = row
        if set(by_mask) != set(range(8)):
            raise ValueError(f"M18 edge masks incomplete for {state}")
        edge_ids = sorted(by_mask[0])
        edge_counts.setdefault(state, len(edge_ids))
        adjoint = {int(row["node_id"]): row for row in read_csv(
            M14_BUILD / state / "adjoint.csv")}
        factor_totals = {factor: 0.0 for factor in FACTORS}
        factor_abs_edges = {factor: 0.0 for factor in FACTORS}
        buckets = {bucket: {factor: 0.0 for factor in FACTORS}
                   for bucket in BUCKETS}
        bucket_abs = {bucket: {factor: 0.0 for factor in FACTORS}
                      for bucket in BUCKETS}
        state_edge_rows: list[dict[str, Any]] = []
        for edge_id in edge_ids:
            game = {mask: float(by_mask[mask][edge_id]["electron_flux"])
                    for mask in range(8)}
            edge_shapley = shapley_values(game)
            maximum_edge_shapley_closure = max(
                maximum_edge_shapley_closure,
                abs(sum(edge_shapley.values()) - (game[7] - game[0])))
            base = by_mask[0][edge_id]
            node0 = int(base["node0"])
            node1 = int(base["node1"])
            lambda0 = 0.0 if int(base["node0_constrained"]) else float(
                adjoint[node0]["lambda_electron"])
            lambda1 = 0.0 if int(base["node1_constrained"]) else float(
                adjoint[node1]["lambda_electron"])
            weight = lambda1 - lambda0
            bucket = edge_bucket(base)
            weighted = {factor: weight * value
                        for factor, value in edge_shapley.items()}
            for factor, value in weighted.items():
                factor_totals[factor] += value
                factor_abs_edges[factor] += abs(value)
                buckets[bucket][factor] += value
                bucket_abs[bucket][factor] += abs(value)
            state_edge_rows.append({
                "state": state,
                "edge_id": edge_id,
                "node0": node0,
                "node1": node1,
                "bucket": bucket,
                "midpoint_x_um": base["midpoint_x"],
                "midpoint_y_um": base["midpoint_y"],
                "node0_contacts": base["node0_contacts"],
                "node1_contacts": base["node1_contacts"],
                "adjoint_edge_weight": weight,
                **{f"{factor}_A_per_um": value
                   for factor, value in weighted.items()},
                "total_A_per_um": sum(weighted.values()),
            })

        baseline_current = float(state_row["baseline_current_A_per_um"])
        closure_by_factor: dict[str, float] = {}
        for factor in FACTORS:
            reference = float(m17_factor_index[(state, factor)]
                              ["adjoint_weighted_current_A_per_um"])
            closure_by_factor[factor] = factor_totals[factor] - reference
        max_abs_closure = max(abs(value) for value in closure_by_factor.values())
        max_rel_closure = max(
            abs(closure_by_factor[factor]) / max(
                abs(factor_totals[factor]),
                abs(float(m17_factor_index[(state, factor)]
                          ["adjoint_weighted_current_A_per_um"])), 1.0e-300)
            for factor in FACTORS)
        for bucket in BUCKETS:
            for factor in FACTORS:
                value = buckets[bucket][factor]
                bucket_rows.append({
                    "state": state,
                    "device": state_row["device"],
                    "drain_voltage_V": state_row["drain_voltage_V"],
                    "gate_voltage_V": state_row["gate_voltage_V"],
                    "bucket": bucket,
                    "factor": factor,
                    "adjoint_weighted_current_A_per_um": value,
                    "relative_to_baseline_current": (
                        value / max(abs(baseline_current), 1.0e-300)),
                    "absolute_edge_support_A_per_um": bucket_abs[bucket][factor],
                    "absolute_edge_support_fraction": (
                        bucket_abs[bucket][factor]
                        / max(factor_abs_edges[factor], 1.0e-300)),
                })
        if state == KEY_STATE:
            qf_abs_total = factor_abs_edges["qf_log_imbalance"]
            state_edge_rows.sort(
                key=lambda row: abs(float(row["qf_log_imbalance_A_per_um"])),
                reverse=True)
            cumulative = 0.0
            for rank, row in enumerate(state_edge_rows, start=1):
                cumulative += abs(float(row["qf_log_imbalance_A_per_um"]))
                row["qf_absolute_rank"] = rank
                row["qf_cumulative_absolute_fraction"] = (
                    cumulative / max(qf_abs_total, 1.0e-300))
                key_edge_rows.append(row)

        ledger = {
            "mobility": factor_totals["mobility"],
            "sg_secant_conductance": factor_totals["sg_secant_conductance"],
            "qf_log_imbalance": factor_totals["qf_log_imbalance"],
            "poisson": float(m15_component_index[(state, "poisson")]
                             ["adjoint_weighted_current_A_per_um"]),
            "srh": sum(float(m15_component_index[(state, component)]
                            ["adjoint_weighted_current_A_per_um"])
                       for component in ("electron_srh", "hole_srh")),
            "boundary_gauge": sum(
                float(m15_component_index[(state, component)]
                      ["adjoint_weighted_current_A_per_um"])
                for component in ("electron_gauge", "electron_boundary",
                                  "hole_gauge", "hole_boundary")),
            "other_carrier": sum(
                float(m15_component_index[(state, component)]
                      ["adjoint_weighted_current_A_per_um"])
                for component in ("electron_impact", "hole_transport",
                                  "hole_impact")),
        }
        ledger_total = sum(ledger.values())
        ledger_reference = float(m15_states[state]["component_sum_A_per_um"])
        for component, value in ledger.items():
            ledger_rows.append({
                "state": state,
                "device": state_row["device"],
                "drain_voltage_V": state_row["drain_voltage_V"],
                "gate_voltage_V": state_row["gate_voltage_V"],
                "component": component,
                "adjoint_weighted_current_A_per_um": value,
                "relative_to_baseline_current": (
                    value / max(abs(baseline_current), 1.0e-300)),
            })
        summary_rows.append({
            "state": state,
            "device": state_row["device"],
            "drain_voltage_V": state_row["drain_voltage_V"],
            "gate_voltage_V": state_row["gate_voltage_V"],
            "edge_count": len(edge_ids),
            "maximum_factor_edge_to_node_closure_A_per_um": max_abs_closure,
            "maximum_factor_edge_to_node_relative_closure": max_rel_closure,
            "feedback_ledger_sum_A_per_um": ledger_total,
            "m15_component_sum_A_per_um": ledger_reference,
            "feedback_ledger_closure_A_per_um": ledger_total - ledger_reference,
        })

    write_csv(portable / "m18_state_summary.csv", summary_rows)
    write_csv(portable / "m18_bucket_factor_contributions.csv", bucket_rows)
    write_csv(portable / "m18_key_state_edge_contributions.csv", key_edge_rows)
    write_csv(portable / "m18_feedback_ledger.csv", ledger_rows)

    key_buckets = [row for row in bucket_rows
                   if row["state"] == KEY_STATE
                   and row["factor"] == "qf_log_imbalance"]
    key_ledger = {row["component"]: row for row in ledger_rows
                  if row["state"] == KEY_STATE}
    contact_buckets = {"source_contact_cut", "drain_contact_cut",
                       "substrate_contact_cut", "gate_contact_cut"}
    contact_abs_fraction = sum(
        float(row["absolute_edge_support_fraction"]) for row in key_buckets
        if row["bucket"] in contact_buckets)
    channel_abs_fraction = next(
        float(row["absolute_edge_support_fraction"]) for row in key_buckets
        if row["bucket"] == "internal_channel")
    top_edge_count_50 = next(
        int(row["qf_absolute_rank"]) for row in key_edge_rows
        if float(row["qf_cumulative_absolute_fraction"]) >= 0.5)
    report = {
        "schema": "vela.simplemos.sdevice.m18_qf_edge_localization.report.v1",
        "status": "complete",
        "execution": {
            "state_count": len(summary_rows),
            "edge_probe_count": len(jobs),
            "probe_variant_count": 8,
            "probes_executed_this_run": len(selected),
            "new_sentaurus_execution": False,
            "cpp_change_scope": "diagnostic-only edge output for the M17 secant-factor probe",
            "default_model_changed": False,
        },
        "closure": {
            "maximum_edge_to_node_absolute_A_per_um": max(
                float(row["maximum_factor_edge_to_node_closure_A_per_um"])
                for row in summary_rows),
            "maximum_edge_to_node_relative": max(
                float(row["maximum_factor_edge_to_node_relative_closure"])
                for row in summary_rows),
            "maximum_feedback_ledger_absolute_A_per_um": max(
                abs(float(row["feedback_ledger_closure_A_per_um"]))
                for row in summary_rows),
            "maximum_edge_shapley_flux_closure": maximum_edge_shapley_closure,
        },
        "key_state": {
            "state": KEY_STATE,
            "qf_bucket_relative_to_baseline_current": {
                row["bucket"]: float(row["relative_to_baseline_current"])
                for row in key_buckets
            },
            "qf_bucket_absolute_edge_support_fraction": {
                row["bucket"]: float(row["absolute_edge_support_fraction"])
                for row in key_buckets
            },
            "combined_contact_cut_absolute_edge_support_fraction":
                contact_abs_fraction,
            "internal_channel_absolute_edge_support_fraction":
                channel_abs_fraction,
            "top_edge_count_for_50_percent_absolute_qf_support": top_edge_count_50,
            "feedback_ledger_relative_to_baseline_current": {
                component: float(row["relative_to_baseline_current"])
                for component, row in key_ledger.items()
            },
        },
        "interpretation_limits": [
            "The adjoint edge response is first-order self-consistent sensitivity around the converged Vela state.",
            "A negligible direct SRH operator response does not prove SRH had no historical influence on the converged state.",
            "The edge partition localizes Vela operator response and is not a native Sentaurus residual decomposition.",
            "Sub-fA relative responses remain conditioning-sensitive and are not headline root-cause evidence.",
        ],
    }
    write_json(portable / "m18_qf_edge_localization_report.json", report)
    print(json.dumps({
        "status": "complete",
        "state_count": len(summary_rows),
        "closure": report["closure"],
        "key_state": report["key_state"],
    }, indent=2))


if __name__ == "__main__":
    main()
