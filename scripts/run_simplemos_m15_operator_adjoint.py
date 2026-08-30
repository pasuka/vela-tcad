#!/usr/bin/env python3
"""Run M15 adjoint-weighted SimpleMOS operator-term attribution."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
M14_BUILD = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
             / "m14_adjoint_attribution")
M14_PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "adjoint_attribution")
M11_EFFECTS = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
               / "m11_mobility_factorial/frozen_vela/frozen_factorial_effects.csv")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m15_operator_adjoint")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "operator_adjoint")
RUNNER = REPO / "build-release/vela_example_runner.exe"
COMPONENT_COLUMNS = {
    "electron_transport": ("electron_flux", "lambda_electron"),
    "electron_srh": ("electron_recombination", "lambda_electron"),
    "electron_impact": ("electron_impact", "lambda_electron"),
    "electron_gauge": ("electron_gauge", "lambda_electron"),
    "electron_boundary": ("electron_boundary", "lambda_electron"),
    "hole_transport": ("hole_flux", "lambda_hole"),
    "hole_srh": ("hole_recombination", "lambda_hole"),
    "hole_impact": ("hole_impact", "lambda_hole"),
    "hole_gauge": ("hole_gauge", "lambda_hole"),
    "hole_boundary": ("hole_boundary", "lambda_hole"),
}


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


def read_status(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").strip().splitlines()
    return json.loads(lines[-1])


def run_config(config: Path, runner: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run([str(runner), "--config", str(config)],
                               cwd=REPO, text=True, capture_output=True,
                               env=env, check=False)
    config.with_suffix(".stdout.txt").write_text(completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"M15 probe failed for {config}: "
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
    output.mkdir(parents=True, exist_ok=True)
    portable.mkdir(parents=True, exist_ok=True)

    m14_states = read_csv(M14_PORTABLE / "m14_state_summary.csv")
    m14_substitutions = read_csv(M14_PORTABLE / "m14_substitution_summary.csv")
    m14_all = {row["state"]: row for row in m14_substitutions
               if row["variant"] == "field_all"}
    cases = []
    jobs: list[tuple[Path, Path, str]] = []
    for state_row in m14_states:
        state = state_row["state"]
        source = M14_BUILD / state / "baseline_functional.json"
        source_deck = read_json(source)
        case_dir = output / state
        case_dir.mkdir(parents=True, exist_ok=True)
        variants = {}
        for name, state_file in (
                ("baseline", M14_BUILD / state / "baseline_state.csv"),
                ("sentaurus", M14_BUILD / state / "field_all_state.csv")):
            config = case_dir / f"{name}_carrier_terms.json"
            csv_path = case_dir / f"{name}_carrier_terms.csv"
            deck = dict(source_deck)
            deck.update({
                "simulation_type": "newton_carrier_term_probe",
                "state_file": str(state_file.resolve()),
                "output_csv": str(csv_path.resolve()),
                "carrier_term_probe": {"solved_equation_terms": True},
                "m15": {"read_only": True, "state_role": name},
            })
            deck.pop("residual_output_csv", None)
            write_json(config, deck)
            jobs.append((config, csv_path, f"{state}:{name}"))
            variants[name] = {"config": config, "csv": csv_path}
        cases.append({"state": state, "state_row": state_row,
                      "variants": variants, "m14_all": m14_all[state]})

    selected = [(config, expected, label) for config, expected, label in jobs
                if args.force or not expected.exists()
                or not config.with_suffix(".stdout.txt").exists()]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, args.runner.resolve()): label
                   for config, _expected, label in selected}
        for future in as_completed(futures):
            label = futures[future]
            future.result()
            print(f"completed {label}", flush=True)

    component_rows: list[dict[str, Any]] = []
    region_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    key_node_rows: list[dict[str, Any]] = []
    for case in cases:
        state = case["state"]
        case_dir = M14_BUILD / state
        adjoint_rows = {int(row["node_id"]): row
                        for row in read_csv(case_dir / "adjoint.csv")}
        base_residual = {int(row["node_id"]): row
                         for row in read_csv(case_dir / "baseline_residual.csv")}
        sent_residual = {int(row["node_id"]): row
                         for row in read_csv(case_dir / "field_all_residual.csv")}
        base_terms = {int(row["node_id"]): row for row in read_csv(
            case["variants"]["baseline"]["csv"])}
        sent_terms = {int(row["node_id"]): row for row in read_csv(
            case["variants"]["sentaurus"]["csv"])}
        baseline_current = float(case["state_row"]["baseline_current_A_per_um"])
        contributions: dict[str, float] = {"poisson": 0.0}
        contributions.update({name: 0.0 for name in COMPONENT_COLUMNS})
        by_region = {name: {component: 0.0 for component in contributions}
                     for name in ("source", "channel", "drain", "body")}

        for node in sorted(adjoint_rows):
            adj = adjoint_rows[node]
            x_um = float(adj["x_m"]) * 1.0e6
            y_um = float(adj["y_m"]) * 1.0e6
            zone = region(x_um, y_um)
            delta_poisson = (float(sent_residual[node]["psi_residual"])
                             - float(base_residual[node]["psi_residual"]))
            poisson = -float(adj["lambda_poisson"]) * delta_poisson
            contributions["poisson"] += poisson
            by_region[zone]["poisson"] += poisson
            node_values = {"poisson": poisson}
            for component, (column, lambda_column) in COMPONENT_COLUMNS.items():
                delta = float(sent_terms[node][column]) - float(base_terms[node][column])
                value = -float(adj[lambda_column]) * delta
                contributions[component] += value
                by_region[zone][component] += value
                node_values[component] = value
            if state in {"n21_vd_0p05_vg_0", "n21_vd_0p05_vg_0p05",
                         "n21_vd_0p05_vg_0p8"}:
                key_node_rows.append({
                    "state": state, "node_id": node, "x_um": x_um, "y_um": y_um,
                    "region": zone, **node_values,
                    "total_A_per_um": sum(node_values.values()),
                })

        total = sum(contributions.values())
        reference = float(case["m14_all"]["adjoint_relaxation_current_change_A_per_um"])
        closure = total - reference
        relative_closure = abs(closure) / max(abs(reference), abs(total), 1.0e-300)
        absolute_sum = sum(abs(value) for value in contributions.values())
        for component, value in contributions.items():
            component_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "component": component,
                "adjoint_weighted_current_A_per_um": value,
                "relative_to_baseline_current": value / max(abs(baseline_current), 1.0e-300),
                "absolute_contribution_fraction": abs(value) / max(absolute_sum, 1.0e-300),
            })
        for zone, components in by_region.items():
            zone_total = sum(components.values())
            region_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "region": zone,
                "adjoint_weighted_current_A_per_um": zone_total,
                "relative_to_baseline_current": zone_total / max(abs(baseline_current), 1.0e-300),
                **{f"{name}_A_per_um": value for name, value in components.items()},
            })
        state_rows.append({
            "state": state,
            "device": case["state_row"]["device"],
            "drain_voltage_V": case["state_row"]["drain_voltage_V"],
            "gate_voltage_V": case["state_row"]["gate_voltage_V"],
            "terminal_signed_log10_error_dex": case["state_row"]["terminal_signed_log10_error_dex"],
            "baseline_current_A_per_um": baseline_current,
            "component_sum_A_per_um": total,
            "m14_reference_A_per_um": reference,
            "closure_A_per_um": closure,
            "relative_closure": relative_closure,
            "dominant_absolute_component": max(contributions, key=lambda name: abs(contributions[name])),
            "dominant_absolute_region": max(by_region,
                key=lambda name: abs(sum(by_region[name].values()))),
        })

    write_csv(portable / "m15_state_summary.csv", state_rows)
    write_csv(portable / "m15_component_contributions.csv", component_rows)
    write_csv(portable / "m15_region_contributions.csv", region_rows)
    write_csv(portable / "m15_key_state_node_contributions.csv", key_node_rows)

    factorial = [row for row in read_csv(M11_EFFECTS)
                 if row["response"] == "log10_terminal_current_A_per_um"
                 and row["term"] in {"phumob", "enormal", "hfs"}]
    write_csv(portable / "m15_m11_mobility_crosscheck.csv", factorial)
    key_state = "n21_vd_0p05_vg_0p8"
    key_components = {row["component"]: row for row in component_rows
                      if row["state"] == key_state}
    key_regions = {row["region"]: row for row in region_rows
                   if row["state"] == key_state}
    key_factorial = {row["term"]: float(row["effect"]) for row in factorial
                     if row["state"] == key_state}
    report = {
        "schema": "vela.simplemos.sdevice.m15_operator_adjoint.report.v1",
        "status": "complete",
        "execution": {
            "state_count": len(state_rows),
            "carrier_term_probe_count": 2 * len(state_rows),
            "new_sentaurus_execution": False,
            "cpp_changed": False,
            "default_model_changed": False,
        },
        "closure": {
            "maximum_absolute_A_per_um": max(abs(float(row["closure_A_per_um"]))
                                             for row in state_rows),
            "maximum_relative": max(float(row["relative_closure"])
                                    for row in state_rows),
            "reference": "M14 field_all adjoint relaxation",
        },
        "key_state": {
            "state": key_state,
            "terminal_signed_log10_error_dex": next(
                float(row["terminal_signed_log10_error_dex"]) for row in state_rows
                if row["state"] == key_state),
            "component_relative_to_baseline_current": {
                name: float(row["relative_to_baseline_current"])
                for name, row in key_components.items()},
            "region_relative_to_baseline_current": {
                name: float(row["relative_to_baseline_current"])
                for name, row in key_regions.items()},
            "m11_frozen_log10_current_factorial_effect_dex": key_factorial,
        },
        "interpretation_limits": [
            "This decomposes the Vela residual response to a mapped Sentaurus state; Sentaurus equation-term residuals are not available.",
            "Electron-transport dominance identifies the broad transport operator, not HFS uniquely.",
            "M11 factorial effects are cross-checks on frozen Vela formula response and are not Sentaurus operator terms.",
            "Sub-fA states remain conditioning-sensitive and are retained but not used for headline attribution.",
        ],
    }
    write_json(portable / "m15_operator_adjoint_report.json", report)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
