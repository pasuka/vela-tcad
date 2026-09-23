#!/usr/bin/env python3
"""Prepare strict external-resistor IALMob on/off BVDS continuation decks."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import shutil
from pathlib import Path
from typing import Any

try:
    from prepare_slot_ldmos_vela_ialmob_ablation import (
        IALMobPreparationError,
        OXIDE_REGION,
        SILICON_REGION,
        normalized_physics,
        replace_output_prefix,
    )
except ModuleNotFoundError:  # Imported as scripts.* by regression tests.
    from scripts.prepare_slot_ldmos_vela_ialmob_ablation import (
        IALMobPreparationError,
        OXIDE_REGION,
        SILICON_REGION,
        normalized_physics,
        replace_output_prefix,
    )


BASE_CONFIG = "simulation_06_bvds_external_resistor_final.json"
SOURCE_OUTPUT = "outputs/stages/06_bvds_external_resistor_final"
CASES = ("ialmob_off", "ialmob_on")
MAX_INNER_STEP_V = 1.0
GUMMEL_MAX_ITER = 50


def last_converged_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle) if row.get("converged") == "1"]
    if not rows:
        raise IALMobPreparationError(f"no converged IV row in {path}")
    return rows[-1]


def build_case(
    base: dict[str, Any],
    case: str,
    initial_inner_voltage_V: float,
) -> dict[str, Any]:
    if case not in CASES:
        raise IALMobPreparationError(f"unknown case {case!r}")
    target_output = f"outputs/ialmob_ablation/external_bv/{case}"
    document = replace_output_prefix(copy.deepcopy(base), SOURCE_OUTPUT, target_output)
    document["_comment"] = (
        "Strict IALMob A/B BVDS continuation under the shared 1e12 ohm*um "
        "external resistor. Requested inner-voltage moves are capped at 1.0 V; "
        "the solver adaptively halves failed device steps."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "surface_interface": [SILICON_REGION, OXIDE_REGION],
    }

    mobility = document["solver"]["mobility"]
    document["solver"]["impact_ionization"]["source_jacobian"] = "local_ad"
    mobility.pop("surface", None)
    if case == "ialmob_off":
        mobility["model"] = "masetti_field"
    else:
        mobility["model"] = "masetti_field_lombardi"
        mobility["surface"] = {
            "surface_interface": [SILICON_REGION, OXIDE_REGION]
        }
    sweep = document["sweep"]
    document["solver"]["handoff"]["gummel_max_iter"] = GUMMEL_MAX_ITER
    sweep["initial_state_file"] = (
        f"outputs/ialmob_ablation/probe_60v/{case}/final_state.h5"
    )
    sweep["external_circuit"]["initial_inner_voltage_V"] = initial_inner_voltage_V
    sweep["external_circuit"]["solver"] = "coupled_newton"
    sweep["external_circuit"]["max_inner_voltage_step_V"] = MAX_INNER_STEP_V
    sweep["external_circuit"]["coupled_equation_tolerance"] = 1.0e-6
    sweep["external_circuit"]["current_directional_step"] = 1.0e-5
    sweep["external_circuit"]["coupled_damping_factor"] = 0.5
    sweep["external_circuit"]["coupled_max_line_search_steps"] = 12
    sweep["external_circuit"]["coupled_initial_outer_step_V"] = 25.0
    sweep["external_circuit"]["coupled_min_outer_step_V"] = 0.1
    sweep["external_circuit"]["coupled_max_outer_step_V"] = 5000.0
    sweep["external_circuit"]["coupled_outer_growth_factor"] = 1.5
    sweep["external_circuit"]["coupled_outer_shrink_factor"] = 0.5
    sweep["external_circuit"]["coupled_max_step_retries"] = 16
    sweep["external_circuit"]["coupled_apply_device_update_limit"] = False
    sweep["external_circuit"]["coupled_line_search_mode"] = "residual_filter"
    sweep["external_circuit"]["coupled_filter_gamma"] = 1.0e-4
    sweep["external_circuit"]["coupled_filter_envelope_factor"] = 1.25
    sweep["external_circuit"]["max_iterations"] = 80
    sweep["boundary_control"]["predictor_max_step_factor"] = 1.0
    sweep["boundary_control"]["adaptive_device_continuation"] = True
    sweep["boundary_control"]["resume"] = True
    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": f"{target_output}/newton_attempts.csv",
            "iterations_csv_file": f"{target_output}/newton_iterations.csv",
            "rejected_state_directory": f"{target_output}/rejected_states",
        }
    }
    sweep["write_vtk"] = False
    sweep.pop("vtk_prefix", None)
    return document


def normalized_external_physics(document: dict[str, Any]) -> dict[str, Any]:
    result = normalized_physics(document)
    result["sweep"]["external_circuit"]["initial_inner_voltage_V"] = "CASE_STATE"
    return result


def prepare(bundle: Path, seed_off_from_stage06: bool = False) -> dict[str, Any]:
    bundle = bundle.resolve()
    with (bundle / BASE_CONFIG).open(encoding="utf-8") as handle:
        base = json.load(handle)

    terminal = {
        case: last_converged_row(
            bundle / f"outputs/ialmob_ablation/probe_60v/{case}/iv.csv"
        )
        for case in CASES
    }
    documents = {
        case: build_case(base, case, float(terminal[case]["inner_voltage_V"]))
        for case in CASES
    }
    if normalized_external_physics(documents["ialmob_off"]) != normalized_external_physics(
        documents["ialmob_on"]
    ):
        raise IALMobPreparationError(
            "external-resistor A/B documents differ outside mobility and isolated paths"
        )

    cases: list[dict[str, Any]] = []
    for case, document in documents.items():
        output = bundle / "outputs" / "ialmob_ablation" / "external_bv" / case
        for directory in (
            output,
            output / "boundary_control_checkpoints",
            output / "states",
        ):
            directory.mkdir(parents=True, exist_ok=True)
        filename = f"simulation_ialmob_external_bv_{case}.json"
        (bundle / filename).write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        cases.append(
            {
                "case": case,
                "config": filename,
                "initial_inner_voltage_V": float(terminal[case]["inner_voltage_V"]),
            }
        )

    imported_off_evaluations = 0
    if seed_off_from_stage06:
        source = bundle / SOURCE_OUTPUT / "boundary_control_evaluations.csv"
        target = (
            bundle
            / "outputs/ialmob_ablation/external_bv/ialmob_off"
            / "boundary_control_evaluations.csv"
        )
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copy2(source, target)
        with target.open(newline="", encoding="utf-8") as handle:
            imported_off_evaluations = sum(
                1
                for row in csv.DictReader(handle)
                if row.get("device_converged") == "1" and row.get("state_file")
            )
        if imported_off_evaluations <= 0:
            raise IALMobPreparationError("no reusable off boundary evaluations found")

    manifest = {
        "schema": "vela.slot_ldmos.ialmob_external_bv.v1",
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "external_resistance_ohm_um": 1.0e12,
        "max_inner_voltage_step_V": MAX_INNER_STEP_V,
        "gummel_max_iter": GUMMEL_MAX_ITER,
        "predictor_max_step_factor": 1.0,
        "adaptive_device_continuation": True,
        "external_circuit_solver": "coupled_newton",
        "coupled_outer_voltage_continuation": {
            "initial_step_V": 25.0,
            "min_step_V": 0.1,
            "max_step_V": 5000.0,
            "growth_factor": 1.5,
            "shrink_factor": 0.5,
            "max_retries": 16,
        },
        "breakdown_criterion_A_per_um": 1.0e-7,
        "imported_off_evaluations": imported_off_evaluations,
        "cases": cases,
    }
    (bundle / "ialmob_external_bv_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--seed-off-from-stage06", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(args.bundle, args.seed_off_from_stage06),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
