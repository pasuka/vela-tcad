#!/usr/bin/env python3
"""Run the exact 31-point LDMOS G3 node-local plus AverageBox curve."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    ucrt = r"D:\msys64\ucrt64\bin"
    env["PATH"] = ucrt + os.pathsep + env.get("PATH", "")
    return env


def count_profile_edges(path: Path) -> int:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != ["node0", "node1", "couple_m"]:
            raise ValueError("AverageBox CSV header must be node0,node1,couple_m")
        return sum(1 for _ in reader)


def prepare(
    baseline_path: Path,
    profile_csv: Path,
    output_dir: Path,
    max_iter: int | None = None,
) -> tuple[dict[str, Any], int]:
    config = deepcopy(json.loads(baseline_path.read_text(encoding="utf-8")))
    if config.get("simulation_type") != "dc_sweep":
        raise ValueError("G3 AverageBox curve requires a dc_sweep baseline")
    sweep = config.get("sweep", {})
    if "predictor" in sweep or "continuation" in sweep:
        raise ValueError("G3 AverageBox qualification requires predictor disabled")
    points = [float(value) for value in sweep.get("bias_points", [])]
    if len(points) != 31:
        raise ValueError("G3 AverageBox qualification requires 31 exact points")

    solver = config.get("solver", {})
    mobility = solver.get("mobility", {})
    if "ialmob" in json.dumps(mobility).lower():
        raise ValueError("G3 AverageBox qualification requires IALMob disabled")
    if solver.get("impact_ionization", {}).get("model", "none") != "none":
        raise ValueError("G3 AverageBox qualification requires avalanche disabled")
    quantum = solver.get("quantum_potential", {})
    if quantum.get("enabled", False):
        raise ValueError("G3 AverageBox qualification requires QP disabled")
    if max_iter is not None:
        if max_iter <= 0:
            raise ValueError("max_iter must be positive")
        solver["max_iter"] = max_iter
    solver["contact_boundary_reconstruction"] = "legacy_node_local"

    profile_csv = profile_csv.resolve()
    expected_edges = count_profile_edges(profile_csv)
    geometry = config.setdefault("mesh_geometry", {})
    if geometry.get("node_volume_policy", "barycentric") != "barycentric":
        raise ValueError("G3 AverageBox qualification requires barycentric volumes")
    geometry["node_volume_policy"] = "barycentric"
    geometry["carrier_transport_couple_profile"] = (
        "templates_ldmos_external_averagebox"
    )
    geometry["external_averagebox_couples_file"] = str(profile_csv)
    geometry["external_averagebox_expected_edges"] = expected_edges

    output_dir.mkdir(parents=True, exist_ok=True)
    config["_comment"] = (
        "Templates/LDMOS G3 exact 31-point qualification: explicit node-local "
        "multipolarity source-short reconstruction plus carrier-only external "
        "AverageBox couples; predictor, IALMob, QP and avalanche off."
    )
    config["output_csv"] = str((output_dir / "g3_averagebox_idvg.csv").resolve())
    sweep["write_state_file"] = str(
        (output_dir / "g3_averagebox_state.csv").resolve()
    )
    sweep["write_state_every_point_prefix"] = str(
        (output_dir / "g3_averagebox_point").resolve()
    )
    diagnostics = sweep.setdefault("diagnostics", {})
    if "terminal_balance" in diagnostics:
        diagnostics["terminal_balance"]["csv_file"] = str(
            (output_dir / "terminal_balance.csv").resolve()
        )
    if "srh_balance" in diagnostics:
        diagnostics["srh_balance"]["csv_file"] = str(
            (output_dir / "srh_balance.csv").resolve()
        )
    if "newton_history" in diagnostics:
        history = diagnostics["newton_history"]
        history["csv_file"] = str((output_dir / "newton_history.csv").resolve())
        history["attempts_csv_file"] = str(
            (output_dir / "newton_attempts.csv").resolve()
        )
        history["iterations_csv_file"] = str(
            (output_dir / "newton_iterations.csv").resolve()
        )
        history["rejected_state_directory"] = str(
            (output_dir / "rejected_states").resolve()
        )
    return config, expected_edges


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--profile-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-iter", type=int)
    args = parser.parse_args()

    config, expected_edges = prepare(
        args.baseline_config, args.profile_csv, args.output_dir, args.max_iter
    )
    config_path = args.output_dir / "g3_averagebox_idvg.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(args.runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((args.output_dir / "g3_averagebox_idvg.log").resolve())],
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (args.output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (args.output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    summary = {
        "schema": "vela.templates_ldmos.g3_averagebox_curve_run.v1",
        "return_code": completed.returncode,
        "config": str(config_path.resolve()),
        "curve": config["output_csv"],
        "balance": config["sweep"]["diagnostics"]["srh_balance"]["csv_file"],
        "qualification_bias_points": len(config["sweep"]["bias_points"]),
        "solver_bias_points": len(config["sweep"]["bias_points"]),
        "max_iter": config["solver"]["max_iter"],
        "predictor": "disabled",
        "ialmob": "disabled",
        "contact_boundary_reconstruction": "legacy_node_local",
        "carrier_transport_couple_profile":
            "templates_ldmos_external_averagebox",
        "external_averagebox_edges": expected_edges,
        "external_averagebox_sha256": sha256(args.profile_csv),
    }
    (args.output_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
