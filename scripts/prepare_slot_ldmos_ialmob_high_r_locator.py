#!/usr/bin/env python3
"""Prepare direct-bordered current locators for SLOT-LDMOS BVDS."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


OUTPUT_ROOT = "outputs/ialmob_ablation/direct_bordered_20260822_v5"
LOCATOR_R_OHM_UM = 1.0
PHYSICAL_R_OHM_UM = 1.0e12


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def remap_outer(inner: float, physical_outer: float) -> float:
    current = (physical_outer - inner) / PHYSICAL_R_OHM_UM
    return current


def prepare_corrected_support_seed(bundle: Path, *, branch_probe: bool = False) -> str:
    """Reclose the last off-state after triangle-GSS material-support fixes."""
    deck = read_json(bundle / "simulation_direct_bordered_ialmob_off_bvds.json")
    bias = 15.856737161516595
    output_name = (
        "ialmob_off_corrected_triangle_support_branch_probe"
        if branch_probe else "ialmob_off_corrected_triangle_support_seed"
    )
    output = f"{OUTPUT_ROOT}/{output_name}"
    (bundle / output).mkdir(parents=True, exist_ok=True)
    deck["_comment"] = (
        "Fixed-voltage recovery after excluding nontransport cells from the "
        "triangle-GSS avalanche residual. Old high-voltage checkpoints are "
        "initial guesses only and must be reclosed against the corrected PDE."
    )
    deck["output_csv"] = f"{output}/iv.csv"
    deck["solver"]["verbose"] = False
    deck["solver"]["impact_ionization"]["source_jacobian"] = "local_ad"
    sweep = deck["sweep"]
    sweep["bias_points"] = [bias, bias + 1.0e-3] if branch_probe else [bias]
    sweep["start"] = bias
    sweep["stop"] = sweep["bias_points"][-1]
    sweep["initial_state_file"] = (
        f"{OUTPUT_ROOT}/ialmob_off_corrected_triangle_support_seed/final_state.h5"
        if branch_probe else
        f"{OUTPUT_ROOT}/ialmob_off_direct_current_locator/"
        "states/state_bias_15p856737.h5"
    )
    sweep["write_state_file"] = f"{output}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["external_circuit"]["enabled"] = False
    sweep["continuation"] = {"arclength": {"enabled": False}}
    sweep["boundary_control"].update(
        {
            "checkpoint_directory": f"{output}/checkpoints",
            "evaluation_csv": f"{output}/boundary_evaluations.csv",
            "resume": False,
        }
    )
    sweep["diagnostics"]["newton_history"].update(
        {
            "attempts_csv_file": f"{output}/newton_attempts.csv",
            "iterations_csv_file": f"{output}/newton_iterations.csv",
            "rejected_state_directory": f"{output}/rejected_states",
        }
    )
    name = (
        "simulation_corrected_triangle_support_branch_probe_off.json"
        if branch_probe else "simulation_corrected_triangle_support_seed_off.json"
    )
    write_json(bundle / name, deck)
    return name


def prepare_case(
    bundle: Path,
    *,
    case: str,
    source_config: str,
    initial_state: str,
    initial_inner: float,
    initial_physical_outer: float,
    previous_state: str,
    previous_inner: float,
    previous_physical_outer: float,
    initial_step: float,
    single_step: bool,
    disable_secant: bool,
    linear_solver: str,
    output_variant: str = "",
) -> str:
    deck = read_json(bundle / source_config)
    smoke_variant = "_nosecant" if disable_secant else ""
    if linear_solver == "direct_bordered_qr":
        smoke_variant += "_qr"
    if single_step and output_variant:
        output_suffix = f"{case}{output_variant}_smoke{smoke_variant}"
    elif single_step:
        output_suffix = f"{case}_direct_current_locator_inexact_smoke{smoke_variant}"
    else:
        output_suffix = f"{case}_direct_current_locator{output_variant}"
    output = f"{OUTPUT_ROOT}/{output_suffix}"
    (bundle / output).mkdir(parents=True, exist_ok=True)
    (bundle / output / "states").mkdir(parents=True, exist_ok=True)
    initial_outer = initial_physical_outer
    active_initial_step = min(initial_step, 2.0e-11) if single_step else initial_step
    targets = (
        [initial_outer, initial_outer + active_initial_step]
        if single_step
        else [initial_outer, 1.0e-8, 5.0e-8, 8.0e-8, 1.0e-7, 1.2e-7]
    )
    deck["_comment"] = (
        "Direct-bordered current locator. The scalar row is I-Target=0, "
        "using coupled_voltage_coefficient=0 and R=1; endpoint states are "
        "subsequently reclosed with the physical R=1e12 load line."
    )
    deck["output_csv"] = f"{output}/iv.csv"
    # The production Newton solver is intentionally verbose by default.  A
    # locator retry can evaluate many rejected trial states, so keep the
    # evidence in the configured CSV diagnostics instead of flooding stdout
    # (which can cause an external runner to terminate before CSV buffers are
    # flushed).
    deck["solver"]["verbose"] = False
    sweep = deck["sweep"]
    sweep["bias_points"] = targets
    sweep["start"] = targets[0]
    sweep["stop"] = targets[-1]
    sweep["initial_state_file"] = initial_state
    sweep["write_state_file"] = f"{output}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep.pop("voltage_to_current", None)
    circuit = sweep["external_circuit"]
    circuit.update(
        {
            "enabled": True,
            "solver": "coupled_newton",
            "coupled_linear_solver": linear_solver,
            "current_directional_step": (
                1.0e-7 if output_variant else circuit.get(
                    "current_directional_step", 1.0e-5)
            ),
            "initial_inner_voltage_V": initial_inner,
            "resistance_ohm_um": LOCATOR_R_OHM_UM,
            "coupled_voltage_coefficient": 0.0,
            "residual_tolerance_V": (
                1.0e-16 if output_variant else 1.0e-12
            ),
            "voltage_tolerance_V": (
                1.0e-17 if output_variant else 1.0e-8
            ),
            "coupled_equation_tolerance": 1.0e-6,
            "coupled_initial_outer_step_V": active_initial_step,
            "coupled_min_outer_step_V": (
                1.0e-17 if output_variant else 1.0e-13
            ),
            "coupled_max_outer_step_V": 3.0e-8,
            "coupled_outer_growth_factor": 4.0,
            "coupled_outer_shrink_factor": 0.5,
            "coupled_max_step_retries": 20,
            # Locator-only globalization: permit a temporary PDE-block rise
            # while the load block moves off an already-converged device
            # state.  The physical-R endpoint reclosure restores the strict
            # production envelope of 2.
            "coupled_filter_envelope_factor": 2.0,
            "coupled_linearization_audit": single_step,
            "coupled_inexact_device_forcing": {
                "enabled": not bool(output_variant),
                "max_equation_tolerance": 5.0e-6,
                "load_activation_ratio": 100.0,
                "zero_converged_residual": not bool(output_variant),
            },
        }
    )
    if single_step:
        # A smoke test should produce a bounded accept/reject result.  Do not
        # allow the normal production retry budget (20 x 80 Newton steps) to
        # turn a developer regression check into a multi-hour continuation.
        circuit["max_iterations"] = 8
        circuit["coupled_max_step_retries"] = 0
    sweep["continuation"] = {"arclength": {"enabled": False}}
    if not disable_secant:
        sweep["continuation"]["arclength"].update(
            {
                "initial_secant_state_file": previous_state,
                "initial_secant_bias_V": previous_inner,
            }
        )
    sweep["boundary_control"].update(
        {
            "checkpoint_directory": f"{output}/checkpoints",
            "evaluation_csv": f"{output}/boundary_evaluations.csv",
            "resume": not bool(output_variant),
        }
    )
    sweep["diagnostics"]["newton_history"].update(
        {
            "attempts_csv_file": f"{output}/newton_attempts.csv",
            "iterations_csv_file": f"{output}/newton_iterations.csv",
            "rejected_state_directory": f"{output}/rejected_states",
        }
    )
    name = f"simulation_direct_bordered_{case}_high_r_locator.json"
    write_json(bundle / name, deck)
    return name


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--single-step", action="store_true")
    parser.add_argument("--disable-secant", action="store_true")
    parser.add_argument("--corrected-support-seed", action="store_true")
    parser.add_argument("--corrected-support-branch-probe", action="store_true")
    parser.add_argument("--use-corrected-support-seed", action="store_true")
    parser.add_argument(
        "--linear-solver",
        choices=("direct_bordered", "direct_bordered_qr"),
        default="direct_bordered",
    )
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    if args.corrected_support_seed:
        print(json.dumps({
            "ialmob_off_corrected_support_seed":
                prepare_corrected_support_seed(bundle)
        }, sort_keys=True))
        return 0
    if args.corrected_support_branch_probe:
        print(json.dumps({
            "ialmob_off_corrected_support_branch_probe":
                prepare_corrected_support_seed(bundle, branch_probe=True)
        }, sort_keys=True))
        return 0
    outputs = {
        "ialmob_off": prepare_case(
            bundle,
            case="ialmob_off",
            source_config="simulation_direct_bordered_ialmob_off_bvds.json",
            initial_state=(
                f"{OUTPUT_ROOT}/ialmob_off_corrected_triangle_support_branch_probe/final_state.h5"
                if args.use_corrected_support_seed else
                f"{OUTPUT_ROOT}/ialmob_off_direct_current_locator/states/state_bias_15p856737.h5"
            ),
            initial_inner=(
                15.857737161516594
                if args.use_corrected_support_seed else 15.856737161516595
            ),
            initial_physical_outer=(
                4.2444335526758775e-10
                if args.use_corrected_support_seed else 3.32361934414161e-9
            ),
            previous_state=(
                f"{OUTPUT_ROOT}/ialmob_off_corrected_triangle_support_seed/final_state.h5"
                if args.use_corrected_support_seed else
                f"{OUTPUT_ROOT}/ialmob_off_direct_current_locator/states/state_bias_15p880551.h5"
            ),
            previous_inner=(
                15.856737161516595
                if args.use_corrected_support_seed else 15.880550769659971
            ),
            previous_physical_outer=(
                4.244286026394933e-10
                if args.use_corrected_support_seed else 2.68361934414161e-9
            ),
            # The corrected branch has dI/dV roughly two orders of magnitude
            # smaller than the contaminated legacy checkpoint. Start with a
            # genuinely local tangent probe; production continuation can grow
            # this step after the first accepted point.
            initial_step=(2.0e-16 if args.use_corrected_support_seed else 2.56e-9),
            single_step=args.single_step,
            disable_secant=args.disable_secant,
            linear_solver=args.linear_solver,
            output_variant=(
                "_corrected_support" if args.use_corrected_support_seed else ""
            ),
        ),
        "ialmob_on": prepare_case(
            bundle,
            case="ialmob_on",
            source_config="simulation_direct_bordered_ialmob_on_bvds.json",
            initial_state=f"{OUTPUT_ROOT}/ialmob_on_seed/final_state.h5",
            initial_inner=0.8078552725248964,
            initial_physical_outer=5.9253738367672866e-11,
            previous_state="outputs/ialmob_ablation/probe_60v/ialmob_on/final_state.h5",
            previous_inner=0.8068552725248964,
            previous_physical_outer=64.9522443569,
            initial_step=2.0e-12,
            single_step=args.single_step,
            disable_secant=args.disable_secant,
            linear_solver=args.linear_solver,
        ),
    }
    print(json.dumps(outputs, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
