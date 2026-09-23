#!/usr/bin/env python3
"""Prepare current-controlled threshold locators for strict IALMob A/B BVDS."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


OUTPUT_ROOT = "outputs/ialmob_ablation/direct_bordered_20260822_v5"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def prepare_case(
    bundle: Path,
    *,
    case: str,
    source_config: str,
    initial_state: str,
    switch_voltage: float,
    current_points: list[float],
) -> str:
    deck = read_json(bundle / source_config)
    relative_output = f"{OUTPUT_ROOT}/{case}_current_locator"
    (bundle / relative_output).mkdir(parents=True, exist_ok=True)
    deck["_comment"] = (
        "Current-controlled pre-location of the BVDS bracket. The device "
        "residual and local-AD physics are identical to the strict direct-"
        "bordered A/B deck; only the scalar continuation parameter changes."
    )
    deck["output_csv"] = f"{relative_output}/iv.csv"
    # Match the acceptance scale of the production coupled Newton solve.  The
    # direct-bordered checkpoints are converged to 1e-6 in the device block;
    # asking this locator bootstrap to re-polish them to 1e-9 caused an
    # unnecessary 80-iteration stall before current control even started.
    deck["solver"]["abstol"] = 1.0e-6
    sweep = deck["sweep"]
    sweep["external_circuit"]["enabled"] = False
    sweep["bias_points"] = [switch_voltage]
    sweep["start"] = switch_voltage
    sweep["stop"] = switch_voltage
    sweep["initial_state_file"] = initial_state
    sweep["write_state_file"] = f"{relative_output}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{relative_output}/states/state"
    sweep["voltage_to_current"] = {
        "enabled": True,
        "switch_voltage_V": switch_voltage,
        "current_direction": 1.0,
        "current_points_A_per_um": current_points,
        "current_tolerance_A_per_um": 1.0e-11,
        "voltage_tolerance_V": 1.0e-8,
        "max_inner_voltage_step_V": 0.5,
        "max_bracket_steps": 100,
        "max_iterations": 50,
    }
    sweep["boundary_control"].update(
        {
            "checkpoint_directory": f"{relative_output}/checkpoints",
            "evaluation_csv": f"{relative_output}/boundary_evaluations.csv",
            "resume": True,
            "adaptive_device_continuation": True,
            "predictor_max_step_factor": 4.0,
        }
    )
    sweep["diagnostics"]["newton_history"].update(
        {
            "attempts_csv_file": f"{relative_output}/newton_attempts.csv",
            "iterations_csv_file": f"{relative_output}/newton_iterations.csv",
            "rejected_state_directory": f"{relative_output}/rejected_states",
        }
    )
    output_name = f"simulation_direct_bordered_{case}_current_locator.json"
    write_json(bundle / output_name, deck)
    return output_name


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    outputs = {
        "ialmob_off": prepare_case(
            bundle,
            case="ialmob_off",
            source_config="simulation_direct_bordered_ialmob_off_bvds.json",
            initial_state=f"{OUTPUT_ROOT}/ialmob_off/states/state_bias_15p829142.h5",
            switch_voltage=15.829142186578611,
            current_points=[1.0e-8, 5.0e-8, 8.0e-8, 1.0e-7, 1.2e-7],
        ),
        "ialmob_on": prepare_case(
            bundle,
            case="ialmob_on",
            source_config="simulation_direct_bordered_ialmob_on_bvds.json",
            initial_state=f"{OUTPUT_ROOT}/ialmob_on_seed/final_state.h5",
            switch_voltage=0.8078552725248964,
            current_points=[1.0e-9, 1.0e-8, 5.0e-8, 8.0e-8, 1.0e-7, 1.2e-7],
        ),
    }
    print(json.dumps(outputs, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
