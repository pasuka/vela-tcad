#!/usr/bin/env python3
"""Prepare paired, avalanche-free SLOT-LDMOS parameter-parity ladders.

The two decks share the converged zero-bias Vela state and differ only in the
Enhanced Lombardi (IALMob) mobility term and their output paths.  Dense low-
voltage landmarks keep the depleted-region quasi-Fermi potentials on the Vela
device manifold instead of importing their numerically non-unique Sentaurus
values at 1 V.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


CASES = ("ialmob_off", "ialmob_on")
SILICON_REGION = "Silicon_1"
OXIDE_REGION = "Oxide_1"
OUTPUT_ROOT = "outputs/parameter_parity_avalanche_off_ladder"
DIRECT_RESTART_ROOT = "outputs/parameter_parity_avalanche_off_direct_restart"
DIRECT_POISSON_ROOT = "outputs/parameter_parity_direct_poisson_reclosure"
EQUILIBRIUM_STATE = "outputs/stages/00_equilibrium/final_state.csv"


class PreparationError(RuntimeError):
    """Raised when a strict parameter-parity ladder cannot be prepared."""


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def ladder_points(stop_voltage_V: float) -> list[float]:
    """Return the validated dense low-voltage device-manifold landmarks."""
    if stop_voltage_V <= 0.0 or stop_voltage_V > 60.0:
        raise PreparationError("stop voltage must be in (0, 60] V")
    landmarks = [
        1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5,
        1.0e-4, 2.0e-4, 5.0e-4,
        1.0e-3, 2.0e-3, 5.0e-3,
        0.01, 0.02, 0.03, 0.04, 0.05,
        0.06, 0.07, 0.08, 0.09, 0.10,
        0.12, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20,
        0.25, 0.30, 0.40, 0.50, 0.75, 1.0,
        1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.40, 1.50,
        1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0,
        15.0, 20.0, 30.0, 40.0, 50.0, 60.0,
    ]
    points = [point for point in landmarks if point < stop_voltage_V]
    if not points or points[-1] != stop_voltage_V:
        points.append(float(stop_voltage_V))
    return points


def _set_case_mobility(document: dict[str, Any], case: str) -> None:
    mobility = document["solver"]["mobility"]
    mobility.pop("surface", None)
    if case == "ialmob_off":
        mobility["model"] = "masetti_field"
    elif case == "ialmob_on":
        mobility["model"] = "masetti_field_lombardi"
        mobility["surface"] = {
            "surface_interface": [SILICON_REGION, OXIDE_REGION]
        }
    else:
        raise PreparationError(f"unknown case {case!r}")


def _replace_case_paths(value: Any, case: str, replacement: str) -> Any:
    needle = f"{OUTPUT_ROOT}/{case}"
    if isinstance(value, dict):
        return {
            key: _replace_case_paths(child, case, replacement)
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [_replace_case_paths(child, case, replacement) for child in value]
    if isinstance(value, str):
        return value.replace(needle, replacement)
    return value


def _normalized(document: dict[str, Any], case: str) -> dict[str, Any]:
    value = copy.deepcopy(document)
    value.pop("_comment", None)
    value.pop("_ialmob_ablation", None)
    mobility = value["solver"]["mobility"]
    mobility["model"] = "<controlled-mobility>"
    mobility.pop("surface", None)
    value = _replace_case_paths(value, case, "<output>")
    return _replace_literal_paths(
        value, f"{DIRECT_RESTART_ROOT}/{case}", "<output>")


def _replace_literal_paths(value: Any, needle: str, replacement: str) -> Any:
    if isinstance(value, dict):
        return {
            key: _replace_literal_paths(child, needle, replacement)
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [
            _replace_literal_paths(child, needle, replacement) for child in value
        ]
    if isinstance(value, str):
        return value.replace(needle, replacement)
    return value


def build_case(base: dict[str, Any], case: str, stop_voltage_V: float) -> dict[str, Any]:
    document = copy.deepcopy(base)
    output = f"{OUTPUT_ROOT}/{case}"
    document["_comment"] = (
        "Avalanche-free parameter-parity ladder from the common converged "
        "zero-bias Vela state. No state extrapolation is used."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "common_seed_state": EQUILIBRIUM_STATE,
        "impact_ionization": "disabled",
    }
    document["output_csv"] = f"{output}/iv.csv"

    solver = document["solver"]
    solver.pop("impact_ionization", None)
    solver["verbose"] = False
    solver["method"] = "newton"
    solver["max_iter"] = 40
    solver["line_search"] = True
    solver["line_search_mode"] = "merit"
    solver["max_update"] = 0.0
    solver["quasi_fermi_update_limit_mode"] = "componentwise_poisson_recorrect"
    solver["quasi_fermi_update_limit_V"] = 0.1
    solver["quasi_fermi_update_limit_minority_V"] = 0.05
    solver["carrier_regularization_scale"] = 1.0e-8
    solver["poisson_only_quasi_fermi_initialization"] = "preserve_initial"
    solver.setdefault("handoff", {})["gummel_max_iter"] = 0
    solver["handoff"]["newton_max_iter"] = 40
    solver["handoff"]["fallback"] = "none"
    solver["handoff"]["require_gummel_convergence"] = False
    mobility = solver["mobility"]
    mobility["high_field_driving_force"] = "quasi_fermi_gradient"
    mobility["high_field_gradient_discretization"] = "transport_cell_vector"
    mobility["jacobian_field_derivatives"] = True
    _set_case_mobility(document, case)

    points = ladder_points(stop_voltage_V)
    sweep = document["sweep"]
    sweep.pop("external_circuit", None)
    sweep["bias_points"] = points
    sweep["start"] = points[0]
    sweep["stop"] = points[-1]
    sweep["step"] = 1.0
    sweep["initial_step"] = points[0]
    sweep["min_step"] = 1.0e-10
    sweep["max_step"] = 0.25
    sweep["growth_factor"] = 1.6
    sweep["shrink_factor"] = 0.5
    sweep["max_retries"] = 20
    sweep["stop_on_failure"] = True
    sweep["initial_state_file"] = EQUILIBRIUM_STATE
    sweep["write_state_file"] = f"{output}/final_state.csv"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["write_vtk"] = False
    sweep["continuation"] = {
        "predictor": {"mode": "none"},
        "arclength": {"enabled": False},
    }
    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": f"{output}/newton_attempts.csv",
            "iterations_csv_file": f"{output}/newton_iterations.csv",
        }
    }
    return document


def prepare(bundle: Path, stop_voltage_V: float = 1.0) -> dict[str, Any]:
    bundle = bundle.resolve()
    base_path = bundle / "simulation_01_unit_resistor_1v.json"
    seed_path = bundle / EQUILIBRIUM_STATE
    if not base_path.is_file():
        raise PreparationError(f"base deck does not exist: {base_path}")
    if not seed_path.is_file():
        raise PreparationError(f"equilibrium seed does not exist: {seed_path}")
    base = _read_json(base_path)
    documents = {
        case: build_case(base, case, stop_voltage_V) for case in CASES
    }
    if _normalized(documents["ialmob_off"], "ialmob_off") != _normalized(
        documents["ialmob_on"], "ialmob_on"
    ):
        raise PreparationError("IALMob decks differ outside mobility and output paths")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / OUTPUT_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_parameter_parity_avalanche_off_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})

    manifest = {
        "schema": "vela.slot_ldmos.parameter_parity_avalanche_off_ladder.v1",
        "common_seed_state": EQUILIBRIUM_STATE,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "impact_ionization": "disabled",
        "bias_points_V": documents["ialmob_off"]["sweep"]["bias_points"],
        "cases": cases,
    }
    _write_json(bundle / "parameter_parity_avalanche_off_ladder_manifest.json", manifest)
    return manifest


def prepare_direct_restart(
    bundle: Path,
    seed_state: str,
    bias_voltage_V: float,
) -> dict[str, Any]:
    """Prepare a paired one-point Newton reclosure of an imported state."""
    if bias_voltage_V <= 0.0:
        raise PreparationError("direct-restart bias must be positive")
    bundle = bundle.resolve()
    base_path = bundle / "simulation_01_unit_resistor_1v.json"
    seed_path = bundle / seed_state
    if not base_path.is_file():
        raise PreparationError(f"base deck does not exist: {base_path}")
    if not seed_path.is_file():
        raise PreparationError(f"direct-restart seed does not exist: {seed_path}")
    base = _read_json(base_path)
    documents: dict[str, dict[str, Any]] = {}
    for case in CASES:
        document = build_case(base, case, bias_voltage_V)
        document = _replace_case_paths(
            document, case, f"{DIRECT_RESTART_ROOT}/{case}")
        document["_comment"] = (
            "One-point avalanche-free Newton reclosure of the imported "
            "Sentaurus high-voltage state."
        )
        document["_ialmob_ablation"]["common_seed_state"] = seed_state
        document["_ialmob_ablation"]["restart_bias_V"] = bias_voltage_V
        solver = document["solver"]
        solver["method"] = "gummel_newton"
        solver["max_iter"] = 80
        solver["handoff"]["gummel_max_iter"] = 100
        solver["handoff"]["newton_max_iter"] = 80
        sweep = document["sweep"]
        sweep["bias_points"] = [bias_voltage_V]
        sweep["start"] = bias_voltage_V
        sweep["stop"] = bias_voltage_V
        sweep["initial_step"] = 1.0e-4
        sweep["max_step"] = 1.0e-4
        sweep["max_retries"] = 0
        sweep["initial_state_file"] = seed_state
        documents[case] = document

    normalized_off = _normalized(documents["ialmob_off"], "ialmob_off")
    normalized_on = _normalized(documents["ialmob_on"], "ialmob_on")
    if normalized_off != normalized_on:
        raise PreparationError("direct-restart IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / DIRECT_RESTART_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_parameter_parity_direct_restart_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest = {
        "schema": "vela.slot_ldmos.parameter_parity_direct_restart.v1",
        "common_seed_state": seed_state,
        "bias_voltage_V": bias_voltage_V,
        "impact_ionization": "disabled",
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "cases": cases,
    }
    _write_json(bundle / "parameter_parity_direct_restart_manifest.json", manifest)
    return manifest


def prepare_direct_poisson_reclosure(
    bundle: Path,
    seed_state: str,
    bias_voltage_V: float,
) -> dict[str, Any]:
    """Prepare a Poisson-only high-voltage reclosure preserving imported QFs."""
    if bias_voltage_V <= 0.0:
        raise PreparationError("direct-Poisson bias must be positive")
    bundle = bundle.resolve()
    base_path = bundle / "simulation_01_unit_resistor_1v.json"
    seed_path = bundle / seed_state
    if not base_path.is_file():
        raise PreparationError(f"base deck does not exist: {base_path}")
    if not seed_path.is_file():
        raise PreparationError(f"direct-Poisson seed does not exist: {seed_path}")
    document = build_case(_read_json(base_path), "ialmob_off", bias_voltage_V)
    document = _replace_case_paths(
        document, "ialmob_off", DIRECT_POISSON_ROOT)
    document["_comment"] = (
        "Poisson-only reclosure of the imported high-voltage Sentaurus "
        "state while preserving electron and hole quasi-Fermi potentials."
    )
    document.pop("_ialmob_ablation", None)
    (bundle / DIRECT_POISSON_ROOT / "states").mkdir(
        parents=True, exist_ok=True)
    solver = document["solver"]
    solver["method"] = "poisson_only"
    solver["max_iter"] = 200
    solver["damping_psi"] = 0.2
    solver["poisson_only_quasi_fermi_initialization"] = "preserve_initial"
    solver.pop("handoff", None)
    sweep = document["sweep"]
    sweep["bias_points"] = [bias_voltage_V]
    sweep["start"] = bias_voltage_V
    sweep["stop"] = bias_voltage_V
    sweep["initial_state_file"] = seed_state
    sweep["max_retries"] = 0
    filename = "simulation_parameter_parity_direct_poisson_reclosure.json"
    _write_json(bundle / filename, document)
    manifest = {
        "schema": "vela.slot_ldmos.parameter_parity_direct_poisson.v1",
        "seed_state": seed_state,
        "bias_voltage_V": bias_voltage_V,
        "config": filename,
        "corrected_state": f"{DIRECT_POISSON_ROOT}/final_state.csv",
    }
    _write_json(bundle / "parameter_parity_direct_poisson_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--stop-voltage", type=float, default=1.0)
    parser.add_argument("--direct-restart-state")
    parser.add_argument("--direct-restart-bias", type=float)
    parser.add_argument("--direct-poisson-state")
    parser.add_argument("--direct-poisson-bias", type=float)
    args = parser.parse_args()
    if args.direct_poisson_state is not None:
        if args.direct_poisson_bias is None:
            parser.error("--direct-poisson-state requires --direct-poisson-bias")
        result = prepare_direct_poisson_reclosure(
            args.bundle, args.direct_poisson_state, args.direct_poisson_bias)
    elif args.direct_restart_state is not None:
        if args.direct_restart_bias is None:
            parser.error("--direct-restart-state requires --direct-restart-bias")
        result = prepare_direct_restart(
            args.bundle, args.direct_restart_state, args.direct_restart_bias)
    else:
        result = prepare(args.bundle, args.stop_voltage)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
