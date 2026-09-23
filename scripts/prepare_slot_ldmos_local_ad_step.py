#!/usr/bin/env python3
"""Prepare an isolated coupled-load-line step using the legacy local-AD source."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument(
        "--base-config",
        default="simulation_ialmob_external_bv_ialmob_off.json",
    )
    parser.add_argument("--state-file", required=True)
    parser.add_argument("--inner-voltage-V", type=float, required=True)
    parser.add_argument("--outer-start-V", type=float, required=True)
    parser.add_argument("--outer-target-V", type=float, required=True)
    parser.add_argument(
        "--output-config",
        default="simulation_local_ad_coupled_step.json",
    )
    args = parser.parse_args()

    bundle = args.bundle.resolve()
    with (bundle / args.base_config).open(encoding="utf-8") as handle:
        document = json.load(handle)
    document["_comment"] = (
        "Isolated legacy-triangle local-AD and block-filter coupled step."
    )
    document["solver"]["impact_ionization"]["source_jacobian"] = "local_ad"
    sweep = document["sweep"]
    sweep["initial_state_file"] = args.state_file
    sweep["bias_points"] = [args.outer_start_V, args.outer_target_V]
    sweep["start"] = args.outer_start_V
    sweep["stop"] = args.outer_target_V
    sweep["output_csv"] = "diagnostics/local_ad_step/iv.csv"
    document["output_csv"] = "diagnostics/local_ad_step/iv.csv"
    external = sweep["external_circuit"]
    external["initial_inner_voltage_V"] = args.inner_voltage_V
    external["coupled_initial_outer_step_V"] = abs(
        args.outer_target_V - args.outer_start_V
    )
    external["coupled_max_outer_step_V"] = abs(
        args.outer_target_V - args.outer_start_V
    )
    sweep["boundary_control"]["resume"] = False
    sweep["boundary_control"]["evaluation_csv"] = (
        "diagnostics/local_ad_step/boundary_evaluations.csv"
    )
    sweep["boundary_control"]["checkpoint_directory"] = (
        "diagnostics/local_ad_step/checkpoints"
    )
    sweep["write_state_file"] = "diagnostics/local_ad_step/final_state.h5"
    sweep["write_state_every_point_prefix"] = (
        "diagnostics/local_ad_step/states/state"
    )
    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": "diagnostics/local_ad_step/newton_attempts.csv",
            "iterations_csv_file": "diagnostics/local_ad_step/newton_iterations.csv",
            "rejected_state_directory": "diagnostics/local_ad_step/rejected_states",
        }
    }
    (bundle / "diagnostics" / "local_ad_step").mkdir(
        parents=True, exist_ok=True
    )
    destination = bundle / args.output_config
    destination.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(destination)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
