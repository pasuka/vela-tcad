#!/usr/bin/env python3
"""Prepare paired intrinsic-voltage BVDS sweeps for the Vela IALMob A/B."""

from __future__ import annotations

import argparse
import copy
import json
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


def build_bvds_case(
    base: dict[str, Any], case: str, stop_voltage_V: int = 45
) -> dict[str, Any]:
    if case not in {"ialmob_off", "ialmob_on"}:
        raise IALMobPreparationError(f"unknown case {case!r}")
    if stop_voltage_V < 1:
        raise IALMobPreparationError("stop voltage must be at least 1 V")
    source_output = "outputs/stages/05_avalanche_on_60v"
    target_output = f"outputs/ialmob_ablation/intrinsic_bv/{case}"
    document = replace_output_prefix(copy.deepcopy(base), source_output, target_output)
    document["_comment"] = (
        "Strict intrinsic-drain-voltage IALMob A/B sweep. The pair shares "
        "voltage points, solver, avalanche, mesh, and numerical controls."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "surface_interface": [SILICON_REGION, OXIDE_REGION],
    }
    mobility = document["solver"]["mobility"]
    mobility.pop("surface", None)
    if case == "ialmob_off":
        mobility["model"] = "masetti_field"
    else:
        mobility["model"] = "masetti_field_lombardi"
        mobility["surface"] = {
            "surface_interface": [SILICON_REGION, OXIDE_REGION]
        }
    document["solver"]["handoff"]["gummel_max_iter"] = 0

    sweep = document["sweep"]
    sweep.pop("external_circuit", None)
    sweep.pop("boundary_control", None)
    activation_points = [0.85, 0.90, 0.95]
    half_volt_points = [
        0.5 * value for value in range(2, 2 * stop_voltage_V + 1)
    ]
    sweep["bias_points"] = activation_points + half_volt_points
    sweep["start"] = activation_points[0]
    sweep["stop"] = float(stop_voltage_V)
    sweep["initial_state_file"] = (
        f"outputs/ialmob_ablation/probe_60v/{case}/final_state.h5"
    )
    sweep.pop("diagnostics", None)
    sweep["write_vtk"] = False
    sweep.pop("vtk_prefix", None)
    return document


def prepare(bundle: Path, stop_voltage_V: int = 45) -> dict[str, Any]:
    bundle = bundle.resolve()
    with (bundle / "simulation_05_avalanche_on_60v.json").open(
        encoding="utf-8"
    ) as handle:
        base = json.load(handle)
    documents = {
        case: build_bvds_case(base, case, stop_voltage_V)
        for case in ("ialmob_off", "ialmob_on")
    }
    if normalized_physics(documents["ialmob_off"]) != normalized_physics(
        documents["ialmob_on"]
    ):
        raise IALMobPreparationError(
            "intrinsic BVDS A/B documents differ outside mobility and isolated paths"
        )
    cases: list[dict[str, Any]] = []
    for case, document in documents.items():
        output = bundle / "outputs" / "ialmob_ablation" / "intrinsic_bv" / case
        for directory in (output, output / "states"):
            directory.mkdir(parents=True, exist_ok=True)
        filename = f"simulation_ialmob_intrinsic_bv_{case}.json"
        (bundle / filename).write_text(
            json.dumps(document, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        cases.append({"case": case, "config": filename})
    manifest = {
        "schema": "vela.slot_ldmos.ialmob_intrinsic_bv.v1",
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "bias_points_V": documents["ialmob_off"]["sweep"]["bias_points"],
        "breakdown_criterion_A_per_um": 1.0e-7,
        "cases": cases,
    }
    (bundle / "ialmob_intrinsic_bv_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--stop-voltage", type=int, default=45)
    args = parser.parse_args()
    print(json.dumps(prepare(args.bundle, args.stop_voltage), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
