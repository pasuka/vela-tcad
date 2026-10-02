#!/usr/bin/env python3
"""Run M16 four-factor adjoint decomposition of SimpleMOS electron SG transport."""

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
M15_PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "operator_adjoint")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m16_transport_factor_v2")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "transport_factor")
RUNNER = REPO / "build-release/vela_example_runner.exe"
FACTORS = ("mobility_state", "mobility_drive", "bernoulli_weights",
           "carrier_population")
GROUPED_FACTORS = ("mobility_state", "mobility_drive", "sg_state_kernel")
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
    config.with_suffix(".stdout.txt").write_text(completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"M16 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def region(x_um: float, y_um: float) -> str:
    if x_um > 0.25:
        return "body"
    if y_um < -0.3:
        return "source"
    if y_um > 0.3:
        return "drain"
    return "channel"


def shapley_values(game: dict[int, float],
                   factor_names: tuple[str, ...]) -> dict[str, float]:
    count = len(factor_names)
    denominator = math.factorial(count)
    values = {}
    for index, factor in enumerate(factor_names):
        terms: list[float] = []
        others = [item for item in range(count) if item != index]
        for size in range(count):
            weight = (math.factorial(size) * math.factorial(count - size - 1)
                      / denominator)
            for subset in combinations(others, size):
                mask = sum(1 << item for item in subset)
                terms.append(
                    weight * (game[mask | (1 << index)] - game[mask]))
        values[factor] = math.fsum(terms)
    return values


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

    states = read_csv(M15_PORTABLE / "m15_state_summary.csv")
    m15_components = read_csv(M15_PORTABLE / "m15_component_contributions.csv")
    m15_electron = {row["state"]: row for row in m15_components
                    if row["component"] == "electron_transport"}
    cases: list[dict[str, Any]] = []
    jobs: list[tuple[Path, Path, str]] = []
    for state_row in states:
        state = state_row["state"]
        case_dir = output / state
        case_dir.mkdir(parents=True, exist_ok=True)
        source = read_json(M14_BUILD / state / "baseline_functional.json")
        output_csv = case_dir / "transport_factor_nodes.csv"
        config = case_dir / "transport_factor_probe.json"
        source.update({
            "simulation_type": "electron_transport_factor_probe",
            "state_file": str((M14_BUILD / state / "baseline_state.csv").resolve()),
            "replacement_state_file": str(
                (M14_BUILD / state / "field_all_state.csv").resolve()),
            "output_csv": str(output_csv.resolve()),
            "simplemos_m16": {
                "read_only": True,
                "factor_count": 4,
                "allocation": "shapley",
            },
        })
        source.pop("residual_output_csv", None)
        write_json(config, source)
        jobs.append((config, output_csv, state))
        cases.append({"state": state, "state_row": state_row,
                      "output_csv": output_csv})

    selected = [(config, expected, label) for config, expected, label in jobs
                if args.force or not expected.exists()
                or not config.with_suffix(".stdout.txt").exists()]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, runner): label
                   for config, _expected, label in selected}
        for future in as_completed(futures):
            label = futures[future]
            status = future.result()
            if status["variant_count"] != 16:
                raise ValueError(f"M16 variant matrix incomplete for {label}")
            print(f"completed {label}", flush=True)

    summary_rows: list[dict[str, Any]] = []
    factor_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    grouped_factor_rows: list[dict[str, Any]] = []
    grouped_region_rows: list[dict[str, Any]] = []
    key_node_rows: list[dict[str, Any]] = []
    maximum_node_shapley_closure = 0.0
    maximum_grouped_node_shapley_closure = 0.0
    for case in cases:
        state = case["state"]
        probe_rows = read_csv(case["output_csv"])
        by_mask: dict[int, dict[int, dict[str, str]]] = {}
        for row in probe_rows:
            by_mask.setdefault(int(row["mask"]), {})[int(row["node_id"])] = row
        if set(by_mask) != set(range(16)):
            raise ValueError(f"M16 mask matrix incomplete for {state}")
        node_ids = sorted(by_mask[0])
        adjoint = {int(row["node_id"]): row for row in read_csv(
            M14_BUILD / state / "adjoint.csv")}
        factor_values = {factor: [] for factor in FACTORS}
        grouped_factor_values = {factor: [] for factor in GROUPED_FACTORS}
        by_region = {zone: {factor: 0.0 for factor in FACTORS}
                     for zone in ("source", "channel", "drain", "body")}
        grouped_by_region = {
            zone: {factor: 0.0 for factor in GROUPED_FACTORS}
            for zone in ("source", "channel", "drain", "body")}
        endpoint_weighted_values: list[float] = []
        for node in node_ids:
            game = {mask: float(by_mask[mask][node]["electron_flux"])
                    for mask in range(16)}
            node_shapley = shapley_values(game, FACTORS)
            grouped_game = {
                grouped_mask: game[
                    ((grouped_mask & 1) != 0)
                    + 2 * ((grouped_mask & 2) != 0)
                    + 12 * ((grouped_mask & 4) != 0)]
                for grouped_mask in range(8)
            }
            grouped_node_shapley = shapley_values(
                grouped_game, GROUPED_FACTORS)
            maximum_node_shapley_closure = max(
                maximum_node_shapley_closure,
                abs(math.fsum(node_shapley.values()) - (game[15] - game[0])))
            maximum_grouped_node_shapley_closure = max(
                maximum_grouped_node_shapley_closure,
                abs(math.fsum(grouped_node_shapley.values())
                    - (game[15] - game[0])))
            lambda_electron = float(adjoint[node]["lambda_electron"])
            x_um = float(by_mask[0][node]["x"])
            y_um = float(by_mask[0][node]["y"])
            zone = region(x_um, y_um)
            weighted_values = {}
            grouped_weighted_values = {}
            for factor, flux_value in node_shapley.items():
                value = -lambda_electron * flux_value
                factor_values[factor].append(value)
                by_region[zone][factor] += value
                weighted_values[factor] = value
            for factor, flux_value in grouped_node_shapley.items():
                value = -lambda_electron * flux_value
                grouped_factor_values[factor].append(value)
                grouped_by_region[zone][factor] += value
                grouped_weighted_values[factor] = value
            endpoint_weighted_values.append(
                -lambda_electron * (game[15] - game[0]))
            if state == KEY_STATE:
                key_node_rows.append({
                    "node_id": node,
                    "x_um": x_um,
                    "y_um": y_um,
                    "region": zone,
                    "lambda_electron": lambda_electron,
                    "baseline_electron_flux": game[0],
                    "replacement_electron_flux": game[15],
                    **{f"{factor}_A_per_um": value
                       for factor, value in weighted_values.items()},
                    **{f"grouped_{factor}_A_per_um": value
                       for factor, value in grouped_weighted_values.items()},
                    "total_A_per_um": math.fsum(weighted_values.values()),
                })

        factors = {factor: math.fsum(values)
                   for factor, values in factor_values.items()}
        grouped_factors = {factor: math.fsum(values)
                           for factor, values in grouped_factor_values.items()}
        endpoint_weighted = math.fsum(endpoint_weighted_values)
        total = math.fsum(factors.values())
        reference = float(m15_electron[state]["adjoint_weighted_current_A_per_um"])
        closure = total - reference
        relative_closure = abs(closure) / max(abs(total), abs(reference), 1.0e-300)
        endpoint_closure = endpoint_weighted - reference
        absolute_sum = math.fsum(abs(value) for value in factors.values())
        grouped_absolute_sum = math.fsum(
            abs(value) for value in grouped_factors.values())
        baseline_current = float(case["state_row"]["baseline_current_A_per_um"])
        for factor, value in factors.items():
            factor_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "factor": factor,
                "adjoint_weighted_current_A_per_um": value,
                "relative_to_baseline_current": value / max(abs(baseline_current), 1.0e-300),
                "absolute_contribution_fraction": abs(value) / max(absolute_sum, 1.0e-300),
            })
        for factor, value in grouped_factors.items():
            grouped_factor_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "factor": factor,
                "adjoint_weighted_current_A_per_um": value,
                "relative_to_baseline_current": value / max(abs(baseline_current), 1.0e-300),
                "absolute_contribution_fraction": (
                    abs(value) / max(grouped_absolute_sum, 1.0e-300)),
            })
        for zone, values in by_region.items():
            zone_total = math.fsum(values.values())
            region_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "region": zone,
                **{f"{factor}_A_per_um": value
                   for factor, value in values.items()},
                "total_A_per_um": zone_total,
                "relative_to_baseline_current": (
                    zone_total / max(abs(baseline_current), 1.0e-300)),
            })
        for zone, values in grouped_by_region.items():
            zone_total = math.fsum(values.values())
            grouped_region_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "region": zone,
                **{f"{factor}_A_per_um": value
                   for factor, value in values.items()},
                "total_A_per_um": zone_total,
                "relative_to_baseline_current": (
                    zone_total / max(abs(baseline_current), 1.0e-300)),
            })
        summary_rows.append({
            "state": state,
            "device": case["state_row"]["device"],
            "drain_voltage_V": case["state_row"]["drain_voltage_V"],
            "gate_voltage_V": case["state_row"]["gate_voltage_V"],
            "baseline_current_A_per_um": baseline_current,
            "factor_sum_A_per_um": total,
            "endpoint_weighted_A_per_um": endpoint_weighted,
            "m15_electron_transport_A_per_um": reference,
            "closure_A_per_um": closure,
            "relative_closure": relative_closure,
            "endpoint_closure_A_per_um": endpoint_closure,
            "dominant_absolute_factor": max(factors, key=lambda name: abs(factors[name])),
            "dominant_absolute_grouped_factor": max(
                grouped_factors, key=lambda name: abs(grouped_factors[name])),
        })

    write_csv(portable / "m16_state_summary.csv", summary_rows)
    write_csv(portable / "m16_factor_contributions.csv", factor_rows)
    write_csv(portable / "m16_region_factor_contributions.csv", region_rows)
    write_csv(portable / "m16_grouped_factor_contributions.csv", grouped_factor_rows)
    write_csv(portable / "m16_grouped_region_contributions.csv", grouped_region_rows)
    write_csv(portable / "m16_key_state_node_contributions.csv", key_node_rows)

    key_factors = {row["factor"]: row for row in factor_rows
                   if row["state"] == KEY_STATE}
    key_regions = {row["region"]: row for row in region_rows
                   if row["state"] == KEY_STATE}
    key_grouped_factors = {row["factor"]: row for row in grouped_factor_rows
                           if row["state"] == KEY_STATE}
    report = {
        "schema": "vela.simplemos.sdevice.m16_transport_factor.report.v1",
        "status": "complete",
        "execution": {
            "state_count": len(summary_rows),
            "factor_probe_count": len(jobs),
            "probe_variant_count": 16,
            "probes_executed_this_run": len(selected),
            "new_sentaurus_execution": False,
            "cpp_change_scope": "diagnostic-only production-equation factor probe",
            "default_model_changed": False,
        },
        "closure": {
            "maximum_absolute_A_per_um": max(abs(float(row["closure_A_per_um"]))
                                             for row in summary_rows),
            "maximum_relative": max(float(row["relative_closure"])
                                    for row in summary_rows),
            "maximum_strong_state_relative": max(
                float(row["relative_closure"]) for row in summary_rows
                if float(row["gate_voltage_V"]) >= 0.8),
            "maximum_endpoint_absolute_A_per_um": max(
                abs(float(row["endpoint_closure_A_per_um"]))
                for row in summary_rows),
            "maximum_node_shapley_flux_closure": maximum_node_shapley_closure,
            "maximum_grouped_node_shapley_flux_closure":
                maximum_grouped_node_shapley_closure,
            "reference": "M15 electron_transport adjoint response",
            "strong_state_definition": "Vg >= 0.8 V",
        },
        "key_state": {
            "state": KEY_STATE,
            "factor_relative_to_baseline_current": {
                factor: float(row["relative_to_baseline_current"])
                for factor, row in key_factors.items()
            },
            "factor_absolute_contribution_fraction": {
                factor: float(row["absolute_contribution_fraction"])
                for factor, row in key_factors.items()
            },
            "grouped_factor_relative_to_baseline_current": {
                factor: float(row["relative_to_baseline_current"])
                for factor, row in key_grouped_factors.items()
            },
            "grouped_factor_absolute_contribution_fraction": {
                factor: float(row["absolute_contribution_fraction"])
                for factor, row in key_grouped_factors.items()
            },
            "region_relative_to_baseline_current": {
                zone: float(row["relative_to_baseline_current"])
                for zone, row in key_regions.items()
            },
        },
        "interpretation_limits": [
            "The factors partition the Vela SG operator response; they are not native Sentaurus residual terms.",
            "Mobility-drive contribution covers Vela transport-cell-vector GradQuasiFermi construction and its HFS response, not an unpublished Sentaurus HFS formula.",
            "Mobility-state contribution retains low-field mobility, state, and HFS interactions.",
            "Sub-fA relative responses remain conditioning-sensitive and are not used for headline attribution.",
        ],
    }
    write_json(portable / "m16_transport_factor_report.json", report)
    print(json.dumps({
        "status": "complete",
        "state_count": len(summary_rows),
        "maximum_absolute_closure_A_per_um": report["closure"]["maximum_absolute_A_per_um"],
        "maximum_relative_closure": report["closure"]["maximum_relative"],
        "key_state": report["key_state"],
    }, indent=2))


if __name__ == "__main__":
    main()
