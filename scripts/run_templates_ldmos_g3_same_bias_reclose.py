#!/usr/bin/env python3
"""Reclose one Templates/LDMOS G3 bias from a frozen Sentaurus state."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    if os.name == "nt":
        env["PATH"] = os.pathsep.join(
            [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin", env.get("PATH", "")]
        )
    env["VELA_LINEAR_SOLVER"] = "sparselu"
    return env


def prepare(
    baseline: dict[str, Any], state: Path, bias: float, output: Path,
    poisson_charge_volume_policy: str = "global",
) -> dict[str, Any]:
    config = deepcopy(baseline)
    mobility = config["solver"]["mobility"]
    if "ialmob" in json.dumps(mobility).lower():
        raise ValueError("same-bias G3 reclose requires IALMob disabled")
    if not mobility.get("contact_electric_field_fallback", False):
        raise ValueError("same-bias G3 reclose requires contact HFS fallback")
    config["simulation_type"] = "newton_solve_from_state"
    config["state_file"] = str(state.resolve())
    config["output_state_file"] = str((output / "reclosed_state.csv").resolve())
    config.pop("sweep", None)
    config.pop("output_csv", None)
    if poisson_charge_volume_policy not in {"global", "material_local"}:
        raise ValueError(
            "poisson_charge_volume_policy must be 'global' or 'material_local'")
    discretization = config.setdefault("discretization", {})
    discretization["poisson_charge_volume_policy"] = (
        poisson_charge_volume_policy
    )
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    config["_comment"] = (
        "Single-bias G3 reclose from a Sentaurus frozen state; IALMob and "
        "predictor disabled; all convergence ceilings retained; Poisson "
        f"charge volume={poisson_charge_volume_policy}."
    )
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--poisson-charge-volume-policy",
        choices=("global", "material_local"),
        default="global",
    )
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    config = prepare(
        baseline,
        args.state,
        args.bias,
        output,
        args.poisson_charge_volume_policy,
    )
    config_path = output / "reclose.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(args.runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((output / "reclose.log").resolve())],
        text=True, capture_output=True, env=runner_environment(), check=False,
    )
    (output / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    status: dict[str, Any] = {}
    for line in reversed(completed.stdout.splitlines()):
        try:
            status = json.loads(line)
            break
        except json.JSONDecodeError:
            continue
    summary = {
        "schema": "vela.templates_ldmos.g3_same_bias_reclose.v1",
        "return_code": completed.returncode,
        "bias_V": args.bias,
        "converged": bool(status.get("converged", False)),
        "iterations": status.get("iterations"),
        "initial_residual": status.get("initial_residual"),
        "final_residual": status.get("final_residual"),
        "final_block_residuals": status.get("final_block_residuals"),
        "drain_current_A_per_um": status.get(
            "contact_currents_A_per_um", {}
        ).get("drain"),
        "config": str(config_path.resolve()),
        "state": config["output_state_file"],
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
