#!/usr/bin/env python3
"""Prepare strict direct-bordered IALMob off/on SLOT-LDMOS BVDS decks."""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path
from typing import Any


OFF_PREVIOUS_OUTER_V = 1188.82348632813
OFF_CURRENT_OUTER_V = 1191.94848632813
OFF_PREVIOUS_INNER_V = 15.721531914507386
OFF_CURRENT_INNER_V = 15.723336
ON_PREVIOUS_INNER_V = 0.8068552725248964
ON_SEED_INNER_V = ON_PREVIOUS_INNER_V + 1.0e-3
FINAL_OUTER_V = 110000.0
OUTPUT_ROOT = "outputs/ialmob_ablation/direct_bordered_20260822_v5"


def read_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def last_converged_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row.get("converged") == "1"]
    if not rows:
        raise RuntimeError(f"no converged row in {path}")
    return rows[-1]


def configure_case(
    template: dict[str, Any],
    *,
    case: str,
    output_directory: str,
    initial_state_file: str,
    initial_inner_voltage_v: float,
    first_outer_voltage_v: float,
    final_outer_voltage_v: float,
    initial_outer_step_v: float,
    circuit_solver: str = "coupled_newton",
    previous_state_file: str | None = None,
    previous_inner_voltage_v: float | None = None,
) -> dict[str, Any]:
    document = copy.deepcopy(template)
    document["_comment"] = (
        "Strict IALMob A/B BVDS continuation using the validated direct "
        "bordered device-circuit Newton solver."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "solver": "direct_bordered",
        "source_jacobian": "local_ad",
    }
    document["output_csv"] = f"{output_directory}/iv.csv"
    solver = document["solver"]
    solver["impact_ionization"]["source_jacobian"] = "local_ad"
    mobility = solver["mobility"]
    mobility.pop("surface", None)
    if case == "ialmob_off":
        mobility["model"] = "masetti_field"
    elif case == "ialmob_on":
        mobility["model"] = "masetti_field_lombardi"
        mobility["surface"] = {
            "surface_interface": ["Silicon_1", "Oxide_1"]
        }
    else:
        raise ValueError(f"unknown case {case!r}")

    sweep = document["sweep"]
    sweep["bias_points"] = [first_outer_voltage_v, final_outer_voltage_v]
    sweep["start"] = first_outer_voltage_v
    sweep["stop"] = final_outer_voltage_v
    sweep["initial_state_file"] = initial_state_file
    sweep["write_state_file"] = f"{output_directory}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{output_directory}/states/state"
    sweep["write_vtk"] = False
    sweep.pop("vtk_prefix", None)

    boundary = sweep["boundary_control"]
    boundary["checkpoint_directory"] = f"{output_directory}/checkpoints"
    boundary["evaluation_csv"] = f"{output_directory}/boundary_evaluations.csv"
    boundary["resume"] = True
    boundary["adaptive_device_continuation"] = True

    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": f"{output_directory}/newton_attempts.csv",
            "iterations_csv_file": f"{output_directory}/newton_iterations.csv",
            "rejected_state_directory": f"{output_directory}/rejected_states",
        }
    }

    circuit = sweep["external_circuit"]
    circuit.update(
        {
            "solver": circuit_solver,
            "coupled_linear_solver": "direct_bordered",
            "initial_inner_voltage_V": initial_inner_voltage_v,
            "resistance_ohm_um": 1.0e12,
            "coupled_apply_device_update_limit": True,
            "coupled_damping_factor": 0.5,
            "coupled_equation_tolerance": 1.0e-6,
            "coupled_line_search_mode": "residual_filter",
            "coupled_filter_gamma": 1.0e-4,
            "coupled_filter_envelope_factor": 2.0,
            "coupled_max_line_search_steps": 12,
            "coupled_initial_outer_step_V": initial_outer_step_v,
            "coupled_min_outer_step_V": 0.1,
            "coupled_max_outer_step_V": 30000.0,
            "coupled_outer_growth_factor": 4.0,
            "coupled_outer_shrink_factor": 0.5,
            "coupled_max_step_retries": 16,
            "max_inner_voltage_step_V": 1.0,
            "max_iterations": 80,
        }
    )

    arclength: dict[str, Any] = {"enabled": False}
    if previous_state_file is not None:
        if previous_inner_voltage_v is None:
            raise ValueError("previous inner voltage is required with a previous state")
        arclength.update(
            {
                "initial_secant_state_file": previous_state_file,
                "initial_secant_bias_V": previous_inner_voltage_v,
            }
        )
    sweep["continuation"] = {"arclength": arclength}
    return document


def prepare_output_directory(bundle: Path, relative_path: str) -> None:
    (bundle / relative_path).mkdir(parents=True, exist_ok=True)


def prepare(bundle: Path) -> dict[str, Any]:
    bundle = bundle.resolve()
    template = read_json(bundle / "simulation_local_ad_coupled_step.json")
    output_root = OUTPUT_ROOT

    off = configure_case(
        template,
        case="ialmob_off",
        output_directory=f"{output_root}/ialmob_off",
        initial_state_file=(
            f"{output_root}/ialmob_off/states/state_bias_15p723336.h5"
        ),
        initial_inner_voltage_v=OFF_CURRENT_INNER_V,
        first_outer_voltage_v=OFF_CURRENT_OUTER_V,
        final_outer_voltage_v=FINAL_OUTER_V,
        initial_outer_step_v=12.5,
        previous_state_file=(
            "diagnostics/direct_bordered_step/states/state_bias_15p721532.h5"
        ),
        previous_inner_voltage_v=OFF_PREVIOUS_INNER_V,
    )
    off_path = bundle / "simulation_direct_bordered_ialmob_off_bvds.json"
    write_json(off_path, off)
    prepare_output_directory(bundle, f"{output_root}/ialmob_off")

    on_seed = configure_case(
        template,
        case="ialmob_on",
        output_directory=f"{output_root}/ialmob_on_seed",
        initial_state_file="outputs/ialmob_ablation/probe_60v/ialmob_on/final_state.h5",
        initial_inner_voltage_v=ON_PREVIOUS_INNER_V,
        first_outer_voltage_v=ON_SEED_INNER_V,
        final_outer_voltage_v=ON_SEED_INNER_V,
        initial_outer_step_v=1.0e-3,
    )
    on_seed["sweep"]["external_circuit"]["enabled"] = False
    # This fixed-inner-voltage solve only supplies the second state required by
    # the direct-bordered secant predictor.  The accepted state still satisfies
    # the full local-AD residual; freezing the avalanche-source derivatives is
    # only a bootstrap Jacobian choice.  The strict A/B production sweep below
    # is regenerated with source_jacobian="local_ad".
    on_seed["solver"]["impact_ionization"]["source_jacobian"] = "frozen"
    on_seed["_ialmob_ablation"]["source_jacobian"] = "frozen_bootstrap_only"
    on_seed["solver"]["method"] = "newton"
    on_seed["sweep"]["bias_points"] = [ON_SEED_INNER_V]
    on_seed["sweep"]["start"] = ON_SEED_INNER_V
    on_seed["sweep"]["stop"] = ON_SEED_INNER_V
    on_seed_path = bundle / "simulation_direct_bordered_ialmob_on_seed.json"
    write_json(on_seed_path, on_seed)
    prepare_output_directory(bundle, f"{output_root}/ialmob_on_seed")
    prepare_output_directory(bundle, f"{output_root}/ialmob_on")

    on_main_path = bundle / "simulation_direct_bordered_ialmob_on_bvds.json"
    seed_iv = bundle / output_root / "ialmob_on_seed" / "iv.csv"
    on_main_written = False
    if seed_iv.exists() and seed_iv.stat().st_size > 0:
        terminal = last_converged_row(seed_iv)
        terminal_inner = float(terminal["bias_V"])
        if abs(terminal_inner - ON_SEED_INNER_V) > 1.0e-8:
            raise RuntimeError(
                f"on seed ended at {terminal_inner} V, expected "
                f"{ON_SEED_INNER_V} V"
            )
        terminal_current = float(terminal["current_total_A_per_um"])
        terminal_outer = terminal_inner + 1.0e12 * terminal_current
        on_main = configure_case(
            template,
            case="ialmob_on",
            output_directory=f"{output_root}/ialmob_on",
            initial_state_file=(
                f"{output_root}/ialmob_on_seed/final_state.h5"
            ),
            initial_inner_voltage_v=terminal_inner,
            first_outer_voltage_v=terminal_outer,
            final_outer_voltage_v=FINAL_OUTER_V,
            initial_outer_step_v=12.5,
            previous_state_file=(
                "outputs/ialmob_ablation/probe_60v/ialmob_on/final_state.h5"
            ),
            previous_inner_voltage_v=ON_PREVIOUS_INNER_V,
        )
        write_json(on_main_path, on_main)
        on_main_written = True

    manifest = {
        "schema": "vela.slot_ldmos.direct_bordered_ialmob_bvds.v1",
        "breakdown_criterion_A_per_um": 1.0e-7,
        "external_resistance_ohm_um": 1.0e12,
        "off_config": off_path.name,
        "on_seed_config": on_seed_path.name,
        "on_config": on_main_path.name if on_main_written else None,
        "on_main_ready": on_main_written,
        "output_root": output_root,
    }
    write_json(bundle / "direct_bordered_ialmob_bvds_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args.bundle), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
