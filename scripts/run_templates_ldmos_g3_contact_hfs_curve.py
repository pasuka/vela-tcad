#!/usr/bin/env python3
"""Run the exact 31-point LDMOS G3 curve with contact HFS fallback."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from copy import deepcopy
from pathlib import Path


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    ucrt = r"D:\msys64\ucrt64\bin"
    env["PATH"] = ucrt + os.pathsep + env.get("PATH", "")
    return env


def expanded_bias_points(points: list[float], max_step: float) -> list[float]:
    if max_step <= 0.0:
        return points
    expanded = [points[0]]
    for start, stop in zip(points, points[1:]):
        segments = max(1, math.ceil(abs(stop - start) / max_step - 1.0e-12))
        expanded.extend(
            start + (stop - start) * part / segments
            for part in range(1, segments + 1))
    return expanded


def filter_exact_bias_rows(source: Path, destination: Path,
                           targets: list[float]) -> None:
    with source.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or "bias_V" not in reader.fieldnames:
            raise ValueError(f"missing bias_V column in {source}")
        rows = [row for row in reader if any(
            abs(float(row["bias_V"]) - target) <= 1.0e-10
            for target in targets)]
    if len(rows) != len(targets):
        raise ValueError(
            f"expected {len(targets)} exact rows in {source}, found {len(rows)}")
    with destination.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def prepare(baseline_path: Path, output_dir: Path,
            warmstart_max_step: float,
            max_iter: int | None) -> tuple[dict, list[float]]:
    config = deepcopy(json.loads(baseline_path.read_text(encoding="utf-8")))
    if config.get("simulation_type") != "dc_sweep":
        raise ValueError("contact-HFS curve requires a dc_sweep baseline")
    sweep = config["sweep"]
    if "predictor" in sweep or "continuation" in sweep:
        raise ValueError("contact-HFS qualification requires predictor disabled")
    target_points = [float(value) for value in sweep.get("bias_points", [])]
    if len(target_points) != 31:
        raise ValueError("contact-HFS qualification requires 31 exact bias points")
    mobility = config["solver"]["mobility"]
    if max_iter is not None:
        if max_iter <= 0:
            raise ValueError("max_iter must be positive")
        config["solver"]["max_iter"] = max_iter
    if mobility.get("high_field_driving_force") != "quasi_fermi_gradient":
        raise ValueError("contact-HFS qualification requires a GradQF baseline")
    if "ialmob" in json.dumps(mobility).lower():
        raise ValueError("contact-HFS qualification requires IALMob disabled")
    mobility.update({
        "contact_electric_field_fallback": True,
        "contact_electric_field_fallback_scope": "contact_node_cell",
        "contact_electric_field_fallback_mode": "cell_gradient_magnitude",
    })
    output_dir.mkdir(parents=True, exist_ok=True)
    config["_comment"] = (
        "Single-factor G3 contact-HFS qualification: interior GradQF plus "
        "contact-node-cell P1 ElectricField magnitude; IALMob and predictor off."
    )
    raw_suffix = "_raw" if warmstart_max_step > 0.0 else ""
    config["output_csv"] = str(
        (output_dir / f"g3_contact_hfs_idvg{raw_suffix}.csv").resolve())
    sweep["bias_points"] = expanded_bias_points(
        target_points, warmstart_max_step)
    sweep["write_state_file"] = str(
        (output_dir / "g3_contact_hfs_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str(
        (output_dir / "g3_contact_hfs_point").resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    if "terminal_balance" in diagnostics:
        diagnostics["terminal_balance"]["csv_file"] = str(
            (output_dir / f"terminal_balance{raw_suffix}.csv").resolve())
    if "srh_balance" in diagnostics:
        diagnostics["srh_balance"]["csv_file"] = str(
            (output_dir / f"srh_balance{raw_suffix}.csv").resolve())
    if "newton_history" in diagnostics:
        history = diagnostics["newton_history"]
        history["csv_file"] = str((output_dir / "newton_history.csv").resolve())
        history["attempts_csv_file"] = str(
            (output_dir / "newton_attempts.csv").resolve())
        history["iterations_csv_file"] = str(
            (output_dir / "newton_iterations.csv").resolve())
        history["rejected_state_directory"] = str(
            (output_dir / "rejected_states").resolve())
    return config, target_points


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--warmstart-max-step", type=float, default=0.0,
        help="Insert unscored exact warm-start points no farther apart than this value.")
    parser.add_argument(
        "--max-iter", type=int,
        help="Override only the Newton iteration budget; convergence gates are unchanged.")
    args = parser.parse_args()

    config, target_points = prepare(
        args.baseline_config, args.output_dir, args.warmstart_max_step,
        args.max_iter)
    config_path = args.output_dir / "g3_contact_hfs_idvg.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(args.runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((args.output_dir / "g3_contact_hfs_idvg.log").resolve())],
        text=True, capture_output=True, env=runner_environment(), check=False)
    (args.output_dir / "stdout.txt").write_text(
        completed.stdout, encoding="utf-8")
    (args.output_dir / "stderr.txt").write_text(
        completed.stderr, encoding="utf-8")
    curve_path = Path(config["output_csv"])
    balance_path = Path(
        config["sweep"]["diagnostics"]["srh_balance"]["csv_file"])
    if completed.returncode == 0 and args.warmstart_max_step > 0.0:
        exact_curve = args.output_dir / "g3_contact_hfs_idvg.csv"
        exact_balance = args.output_dir / "srh_balance.csv"
        filter_exact_bias_rows(curve_path, exact_curve, target_points)
        filter_exact_bias_rows(balance_path, exact_balance, target_points)
        curve_path = exact_curve
        balance_path = exact_balance
    summary = {
        "schema": "vela.templates_ldmos.g3_contact_hfs_curve_run.v1",
        "return_code": completed.returncode,
        "config": str(config_path.resolve()),
        "curve": str(curve_path.resolve()),
        "raw_curve": config["output_csv"],
        "balance": str(balance_path.resolve()),
        "raw_balance": config["sweep"]["diagnostics"]["srh_balance"]["csv_file"],
        "qualification_bias_points": len(target_points),
        "solver_bias_points": len(config["sweep"]["bias_points"]),
        "warmstart_max_step_V": args.warmstart_max_step,
        "max_iter": config["solver"]["max_iter"],
        "predictor": "disabled",
        "ialmob": "disabled",
    }
    (args.output_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
