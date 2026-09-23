#!/usr/bin/env python3
"""Prepare fresh, paired SLOT-LDMOS IALMob voltage branches.

Both decks start from the same low-voltage device state.  The first requested
bias differs from the seed bias, so Vela must reclose the corrected nonlinear
system instead of recording the seed as an already-converged point.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


SILICON_REGION = "Silicon_1"
OXIDE_REGION = "Oxide_1"
SEED_INNER_V = 0.008374398259206585
SEED_STATE = "outputs/stages/04_avalanche_activation_1v/final_state.h5"
OUTPUT_ROOT = "outputs/ialmob_ablation/corrected_low_voltage_newton_20260823"
EXTENSION_ROOT = "outputs/ialmob_ablation/corrected_accelerated_20260823"
DENSE_EXTENSION_ROOT = "outputs/ialmob_ablation/corrected_dense_low_voltage_20260823"
POST_DENSE_ROOT = "outputs/ialmob_ablation/corrected_post_dense_20260823"
POINT_TWO_RECOVERY_ROOT = (
    "outputs/ialmob_ablation/corrected_point_two_recovery_20260823"
)
ONE_VOLT_EXTENSION_ROOT = (
    "outputs/ialmob_ablation/corrected_one_volt_extension_20260823"
)
HIGH_VOLTAGE_EXTENSION_ROOT = (
    "outputs/ialmob_ablation/corrected_high_voltage_manifold_20260823"
)
DEVICE_CORRECTOR_ROOT = "outputs/ialmob_ablation/device_manifold_corrector_20260823"


class PreparationError(RuntimeError):
    """Raised when a strict IALMob A/B deck cannot be constructed."""


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def bias_points(stop_voltage_V: float) -> list[float]:
    """Return sparse requested landmarks; DCSweep fills gaps adaptively."""
    if stop_voltage_V < 0.01:
        raise PreparationError("stop voltage must be at least 0.01 V")
    landmarks = [
        0.011, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 8.0,
        10.0, 12.0, 14.0, 15.0, 15.5, 15.75, 15.85, 16.0,
        18.0, 20.0, 25.0, 30.0, 35.0, 38.0, 40.0, 45.0,
    ]
    points = [value for value in landmarks if value < stop_voltage_V]
    if not points or points[-1] != stop_voltage_V:
        points.append(float(stop_voltage_V))
    return points


def _set_mobility(document: dict[str, Any], case: str) -> None:
    mobility = document["solver"]["mobility"]
    mobility.pop("surface", None)
    if case == "ialmob_off":
        mobility["model"] = "masetti_field"
    else:
        mobility["model"] = "masetti_field_lombardi"
        mobility["surface"] = {
            "surface_interface": [SILICON_REGION, OXIDE_REGION]
        }


def build_bootstrap_case(base: dict[str, Any], case: str) -> dict[str, Any]:
    """Build a low-voltage, avalanche-free mobility initialization solve."""
    if case not in {"ialmob_off", "ialmob_on"}:
        raise PreparationError(f"unknown case {case!r}")
    output = f"{OUTPUT_ROOT}/{case}_bootstrap"
    document = _replace_output_paths(
        copy.deepcopy(base),
        "outputs/stages/04_avalanche_activation_1v",
        output,
    )
    document["_comment"] = (
        "Low-voltage mobility bootstrap. Avalanche is omitted only for this "
        "10 mV initialization solve and restored in the production branch."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "phase": "avalanche_free_low_voltage_bootstrap",
        "common_seed_inner_voltage_V": SEED_INNER_V,
        "common_seed_state": SEED_STATE,
    }
    solver = document["solver"]
    solver["verbose"] = False
    solver.pop("impact_ionization", None)
    _set_mobility(document, case)
    sweep = document["sweep"]
    sweep["bias_points"] = [0.01]
    sweep["start"] = 0.01
    sweep["stop"] = 0.01
    sweep["initial_state_file"] = SEED_STATE
    sweep["write_state_file"] = f"{output}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["write_vtk"] = False
    sweep["continuation"] = {"arclength": {"enabled": False}}
    sweep.setdefault("external_circuit", {})["enabled"] = False
    sweep["boundary_control"] = {
        "checkpoint_directory": f"{output}/checkpoints",
        "evaluation_csv": f"{output}/boundary_evaluations.csv",
        "resume": False,
    }
    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": f"{output}/newton_attempts.csv",
            "iterations_csv_file": f"{output}/newton_iterations.csv",
            "rejected_state_directory": f"{output}/rejected_states",
        }
    }
    document["output_csv"] = f"{output}/iv.csv"
    return document


def _replace_output_paths(value: Any, old: str, new: str) -> Any:
    if isinstance(value, dict):
        return {
            key: _replace_output_paths(child, old, new)
            for key, child in value.items()
        }
    if isinstance(value, list):
        return [_replace_output_paths(child, old, new) for child in value]
    if isinstance(value, str):
        return value.replace(old, new)
    return value


def build_case(
    base: dict[str, Any], case: str, stop_voltage_V: float
) -> dict[str, Any]:
    if case not in {"ialmob_off", "ialmob_on"}:
        raise PreparationError(f"unknown case {case!r}")
    output = f"{OUTPUT_ROOT}/{case}"
    document = _replace_output_paths(
        copy.deepcopy(base),
        "outputs/stages/04_avalanche_activation_1v",
        output,
    )
    document["_comment"] = (
        "Fresh corrected intrinsic-voltage branch. Both IALMob cases share "
        "the 8.374 mV seed as an initial guess only; the 10 mV first point "
        "is solved again with local-AD avalanche and corrected material support."
    )
    document["_ialmob_ablation"] = {
        "case": case,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "common_seed_inner_voltage_V": SEED_INNER_V,
        "common_seed_state": SEED_STATE,
        "source_jacobian": "local_ad",
    }

    solver = document["solver"]
    solver["verbose"] = False
    solver["impact_ionization"]["source_jacobian"] = "local_ad"
    # The paired avalanche-free bootstraps are already converged at 10 mV.
    # Re-running the inherited 50 Gummel iterations makes the Lombardi branch
    # spend minutes assembling impact terms before the useful Newton handoff.
    solver.setdefault("handoff", {})["gummel_max_iter"] = 0
    _set_mobility(document, case)

    points = bias_points(stop_voltage_V)
    sweep = document["sweep"]
    sweep["bias_points"] = points
    sweep["start"] = points[0]
    sweep["stop"] = points[-1]
    sweep["step"] = 1.0
    sweep["initial_step"] = 0.01
    sweep["min_step"] = 1.0e-8
    sweep["max_step"] = 0.5
    sweep["max_retries"] = 26
    sweep["stop_on_failure"] = True
    sweep["initial_state_file"] = f"{OUTPUT_ROOT}/{case}_bootstrap/final_state.h5"
    sweep["write_state_file"] = f"{output}/final_state.h5"
    sweep["write_state_every_point_prefix"] = f"{output}/states/state"
    sweep["write_vtk"] = False
    sweep["continuation"] = {"arclength": {"enabled": False}}
    # Retain the circuit object as documented configuration evidence, but make
    # the production branch an intrinsic fixed-drain-voltage solve.
    sweep.setdefault("external_circuit", {})["enabled"] = False
    sweep["boundary_control"] = {
        "checkpoint_directory": f"{output}/checkpoints",
        "evaluation_csv": f"{output}/boundary_evaluations.csv",
        "resume": False,
    }
    sweep["diagnostics"] = {
        "newton_history": {
            "enabled": True,
            "attempts_csv_file": f"{output}/newton_attempts.csv",
            "iterations_csv_file": f"{output}/newton_iterations.csv",
            "rejected_state_directory": f"{output}/rejected_states",
        }
    }
    document["output_csv"] = f"{output}/iv.csv"
    return document


def _normalized_pair(document: dict[str, Any]) -> dict[str, Any]:
    value = copy.deepcopy(document)
    value.pop("_comment", None)
    value.pop("_ialmob_ablation", None)
    mobility = value["solver"]["mobility"]
    mobility["model"] = "<controlled-mobility>"
    mobility.pop("surface", None)
    value = _replace_output_paths(value, f"{OUTPUT_ROOT}/ialmob_off", "<output>")
    value = _replace_output_paths(value, f"{OUTPUT_ROOT}/ialmob_on", "<output>")
    value = _replace_output_paths(value, f"{EXTENSION_ROOT}/ialmob_off", "<output>")
    value = _replace_output_paths(value, f"{EXTENSION_ROOT}/ialmob_on", "<output>")
    value = _replace_output_paths(value, f"{DENSE_EXTENSION_ROOT}/ialmob_off", "<output>")
    value = _replace_output_paths(value, f"{DENSE_EXTENSION_ROOT}/ialmob_on", "<output>")
    value = _replace_output_paths(value, f"{POST_DENSE_ROOT}/ialmob_off", "<output>")
    value = _replace_output_paths(value, f"{POST_DENSE_ROOT}/ialmob_on", "<output>")
    value = _replace_output_paths(
        value, f"{POINT_TWO_RECOVERY_ROOT}/ialmob_off", "<output>"
    )
    value = _replace_output_paths(
        value, f"{POINT_TWO_RECOVERY_ROOT}/ialmob_on", "<output>"
    )
    value = _replace_output_paths(
        value, f"{ONE_VOLT_EXTENSION_ROOT}/ialmob_off", "<output>"
    )
    value = _replace_output_paths(
        value, f"{ONE_VOLT_EXTENSION_ROOT}/ialmob_on", "<output>"
    )
    value = _replace_output_paths(
        value, f"{HIGH_VOLTAGE_EXTENSION_ROOT}/ialmob_off", "<output>"
    )
    value = _replace_output_paths(
        value, f"{HIGH_VOLTAGE_EXTENSION_ROOT}/ialmob_on", "<output>"
    )
    return value


def prepare(bundle: Path, stop_voltage_V: float = 45.0) -> dict[str, Any]:
    bundle = bundle.resolve()
    base_path = bundle / "simulation_04_avalanche_activation_1v.json"
    base = _read_json(base_path)
    seed = bundle / SEED_STATE
    if not seed.is_file():
        raise PreparationError(f"common low-voltage seed does not exist: {seed}")

    bootstraps = {
        case: build_bootstrap_case(base, case)
        for case in ("ialmob_off", "ialmob_on")
    }
    documents = {
        case: build_case(base, case, stop_voltage_V)
        for case in ("ialmob_off", "ialmob_on")
    }
    if _normalized_pair(bootstraps["ialmob_off"]) != _normalized_pair(
        bootstraps["ialmob_on"]
    ):
        raise PreparationError("IALMob bootstraps differ outside mobility and output paths")
    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("IALMob pair differs outside the controlled mobility delta")

    cases: list[dict[str, str]] = []
    for case in ("ialmob_off", "ialmob_on"):
        output = bundle / OUTPUT_ROOT / case
        bootstrap_output = bundle / OUTPUT_ROOT / f"{case}_bootstrap"
        (output / "states").mkdir(parents=True, exist_ok=True)
        (bootstrap_output / "states").mkdir(parents=True, exist_ok=True)
        bootstrap_filename = (
            f"simulation_corrected_low_voltage_{case}_bootstrap.json"
        )
        _write_json(bundle / bootstrap_filename, bootstraps[case])
        filename = f"simulation_corrected_low_voltage_{case}.json"
        _write_json(bundle / filename, documents[case])
        cases.append({
            "case": case,
            "bootstrap_config": bootstrap_filename,
            "config": filename,
        })

    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_branches.v1",
        "common_seed_inner_voltage_V": SEED_INNER_V,
        "common_seed_state": SEED_STATE,
        "mobility_bootstrap_voltage_V": 0.01,
        "first_reclosed_voltage_V": documents["ialmob_off"]["sweep"]["bias_points"][0],
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "breakdown_criterion_A_per_um": 1.0e-7,
        "bias_points_V": documents["ialmob_off"]["sweep"]["bias_points"],
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_branches_manifest.json", manifest)
    return manifest


def prepare_accelerated_extension(
    bundle: Path,
    stop_voltage_V: float = 45.0,
    start_voltage_V: float = 0.05,
) -> dict[str, Any]:
    """Continue a verified paired branch with bounded, faster voltage steps."""
    if abs(start_voltage_V - 0.05) > 1.0e-12:
        raise PreparationError("the validated accelerated start is 0.05 V")
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    requested = [
        0.1, 0.5, 1.0, 2.0, 5.0, 8.0, 10.0, 12.0, 14.0,
        15.0, 15.5, 15.75, 15.85, 16.0, 18.0, 20.0, 25.0,
        30.0, 35.0, 38.0, 40.0, 45.0,
    ]
    points = [point for point in requested if point < stop_voltage_V]
    if not points or points[-1] != stop_voltage_V:
        points.append(float(stop_voltage_V))

    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        seed = (
            bundle / OUTPUT_ROOT / case / "states" / "state_bias_0p050000.h5"
        )
        if not seed.is_file():
            raise PreparationError(f"verified 50 mV extension seed is missing: {seed}")
        document = build_case(base, case, stop_voltage_V)
        document["solver"]["verbose"] = True
        document["solver"]["max_iter"] = 10
        document["solver"]["handoff"]["newton_max_iter"] = 10
        document = _replace_output_paths(
            document, f"{OUTPUT_ROOT}/{case}", f"{EXTENSION_ROOT}/{case}"
        )
        document["_comment"] = (
            "Accelerated continuation from the verified corrected 50 mV "
            "branch state. Secant prediction and <=1 V device steps reduce "
            "low-field work without relaxing nonlinear tolerances."
        )
        document["_ialmob_ablation"]["extension_seed_voltage_V"] = start_voltage_V
        sweep = document["sweep"]
        sweep["bias_points"] = points
        sweep["start"] = points[0]
        sweep["stop"] = points[-1]
        sweep["initial_state_file"] = (
            f"{OUTPUT_ROOT}/{case}/states/state_bias_0p050000.h5"
        )
        sweep["initial_step"] = 0.05
        sweep["growth_factor"] = 1.8
        sweep["max_step"] = 1.0
        sweep["continuation"] = {
            "predictor": {
                "mode": "secant",
                "fields": ["psi", "phin", "phip"],
                "max_extrapolation_ratio": 2.0,
            },
            "arclength": {"enabled": False},
        }
        documents[case] = document

    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("accelerated IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / EXTENSION_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_corrected_accelerated_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_accelerated.v1",
        "seed_voltage_V": start_voltage_V,
        "initial_step_V": 0.05,
        "growth_factor": 1.8,
        "maximum_step_V": 1.0,
        "predictor": "secant",
        "bias_points_V": points,
        "breakdown_criterion_A_per_um": 1.0e-7,
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_accelerated_manifest.json", manifest)
    return manifest


def prepare_dense_low_voltage_extension(
    bundle: Path,
    stop_voltage_V: float = 0.1,
    start_voltage_V: float = 0.05,
    step_voltage_V: float = 0.01,
) -> dict[str, Any]:
    """Follow the converged low-voltage device manifold with dense steps."""
    if abs(start_voltage_V - 0.05) > 1.0e-12:
        raise PreparationError("the validated dense-extension start is 0.05 V")
    if stop_voltage_V <= start_voltage_V:
        raise PreparationError("dense-extension stop must exceed 0.05 V")
    if step_voltage_V <= 0.0:
        raise PreparationError("dense-extension step must be positive")
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    count = int(round((stop_voltage_V - start_voltage_V) / step_voltage_V))
    points = [
        round(start_voltage_V + step_voltage_V * index, 12)
        for index in range(1, count + 1)
    ]
    if not points or abs(points[-1] - stop_voltage_V) > 1.0e-10:
        points.append(float(stop_voltage_V))

    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        seed = bundle / OUTPUT_ROOT / case / "states/state_bias_0p050000.h5"
        if not seed.is_file():
            raise PreparationError(f"verified 50 mV dense-extension seed is missing: {seed}")
        document = build_case(base, case, stop_voltage_V)
        document = _replace_output_paths(
            document, f"{OUTPUT_ROOT}/{case}", f"{DENSE_EXTENSION_ROOT}/{case}"
        )
        document["_comment"] = (
            "Dense device-manifold continuation from the paired converged "
            "50 mV states. No forcing or convergence tolerance is relaxed."
        )
        solver = document["solver"]
        solver["verbose"] = False
        solver["max_iter"] = 20
        solver["handoff"]["newton_max_iter"] = 20
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 4.0
        solver["quasi_fermi_update_limit_mode"] = "uniform_trust_region"
        solver["quasi_fermi_trust_region_growth_factor"] = 2.0
        solver["quasi_fermi_trust_region_max_multiplier"] = 4.0
        solver["quasi_fermi_trust_region_expansion_threshold"] = 0.75
        solver["quasi_fermi_trust_region_shrink_factor"] = 0.5
        solver["quasi_fermi_trust_region_min_multiplier"] = 0.125
        sweep = document["sweep"]
        sweep["bias_points"] = points
        sweep["start"] = points[0]
        sweep["stop"] = points[-1]
        sweep["initial_state_file"] = (
            f"{OUTPUT_ROOT}/{case}/states/state_bias_0p050000.h5"
        )
        sweep["initial_step"] = step_voltage_V
        sweep["growth_factor"] = 1.0
        sweep["max_step"] = step_voltage_V
        sweep["continuation"] = {
            "predictor": {
                "mode": "secant",
                "fields": ["psi", "phin", "phip"],
                "max_extrapolation_ratio": 1.0,
            },
            "arclength": {"enabled": False},
        }
        documents[case] = document

    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("dense IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / DENSE_EXTENSION_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        (output / "rejected_states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_corrected_dense_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_dense_extension.v1",
        "seed_voltage_V": start_voltage_V,
        "step_voltage_V": step_voltage_V,
        "bias_points_V": points,
        "forcing_relaxed": False,
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_dense_extension_manifest.json", manifest)
    return manifest


def prepare_post_dense_extension(
    bundle: Path,
    stop_voltage_V: float = 1.0,
) -> dict[str, Any]:
    """Extend paired, converged 0.1 V states with a gradual voltage grid."""
    if stop_voltage_V <= 0.1:
        raise PreparationError("post-dense stop must exceed 0.1 V")
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    requested = [0.12, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0]
    points = [point for point in requested if point < stop_voltage_V]
    if not points or abs(points[-1] - stop_voltage_V) > 1.0e-12:
        points.append(float(stop_voltage_V))

    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        seed = (
            bundle / DENSE_EXTENSION_ROOT / case /
            "states/state_bias_0p100000.h5"
        )
        if not seed.is_file():
            raise PreparationError(f"converged 0.1 V post-dense seed is missing: {seed}")
        document = build_case(base, case, stop_voltage_V)
        document = _replace_output_paths(
            document, f"{OUTPUT_ROOT}/{case}", f"{POST_DENSE_ROOT}/{case}"
        )
        document["_comment"] = (
            "Post-dense continuation from the paired converged 0.1 V states; "
            "the voltage grid grows gradually without relaxing forcing."
        )
        solver = document["solver"]
        solver["verbose"] = False
        solver["max_iter"] = 20
        solver["handoff"]["newton_max_iter"] = 20
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 4.0
        solver["quasi_fermi_update_limit_mode"] = "uniform_trust_region"
        solver["quasi_fermi_trust_region_growth_factor"] = 2.0
        solver["quasi_fermi_trust_region_max_multiplier"] = 4.0
        solver["quasi_fermi_trust_region_expansion_threshold"] = 0.75
        solver["quasi_fermi_trust_region_shrink_factor"] = 0.5
        solver["quasi_fermi_trust_region_min_multiplier"] = 0.125
        sweep = document["sweep"]
        sweep["bias_points"] = points
        sweep["start"] = points[0]
        sweep["stop"] = points[-1]
        sweep["initial_state_file"] = (
            f"{DENSE_EXTENSION_ROOT}/{case}/states/state_bias_0p100000.h5"
        )
        sweep["initial_step"] = 0.02
        sweep["growth_factor"] = 1.5
        sweep["max_step"] = 0.25
        sweep["continuation"] = {
            "predictor": {
                "mode": "secant",
                "fields": ["psi", "phin", "phip"],
                "max_extrapolation_ratio": 1.5,
            },
            "arclength": {"enabled": False},
        }
        documents[case] = document

    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("post-dense IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / POST_DENSE_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        (output / "rejected_states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_corrected_post_dense_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_post_dense.v1",
        "seed_voltage_V": 0.1,
        "bias_points_V": points,
        "forcing_relaxed": False,
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_post_dense_manifest.json", manifest)
    return manifest


def prepare_point_two_recovery(bundle: Path) -> dict[str, Any]:
    """Recover the paired 0.15--0.20 V branch without secant overshoot."""
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    points = [0.16, 0.17, 0.18, 0.19, 0.2]

    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        seed = (
            bundle / POST_DENSE_ROOT / case /
            "states/state_bias_0p150000.h5"
        )
        if not seed.is_file():
            raise PreparationError(
                f"converged 0.15 V point-two recovery seed is missing: {seed}"
            )
        document = build_case(base, case, points[-1])
        document = _replace_output_paths(
            document, f"{OUTPUT_ROOT}/{case}", f"{POINT_TWO_RECOVERY_ROOT}/{case}"
        )
        document["_comment"] = (
            "Device-manifold recovery from the paired converged 0.15 V states. "
            "Every 10 mV point starts from the preceding converged state; state "
            "extrapolation is disabled to avoid the observed 0.20 V secant overshoot."
        )
        solver = document["solver"]
        solver["verbose"] = False
        solver["max_iter"] = 20
        solver["handoff"]["newton_max_iter"] = 20
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 4.0
        solver["quasi_fermi_update_limit_mode"] = "uniform_trust_region"
        solver["quasi_fermi_trust_region_growth_factor"] = 2.0
        solver["quasi_fermi_trust_region_max_multiplier"] = 4.0
        solver["quasi_fermi_trust_region_expansion_threshold"] = 0.75
        solver["quasi_fermi_trust_region_shrink_factor"] = 0.5
        solver["quasi_fermi_trust_region_min_multiplier"] = 0.125
        sweep = document["sweep"]
        sweep["bias_points"] = points
        sweep["start"] = points[0]
        sweep["stop"] = points[-1]
        sweep["initial_state_file"] = (
            f"{POST_DENSE_ROOT}/{case}/states/state_bias_0p150000.h5"
        )
        sweep["initial_step"] = 0.01
        sweep["growth_factor"] = 1.0
        sweep["max_step"] = 0.01
        sweep["continuation"] = {
            "predictor": {"mode": "none"},
            "arclength": {"enabled": False},
        }
        documents[case] = document

    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("0.20 V recovery IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / POINT_TWO_RECOVERY_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        (output / "rejected_states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_corrected_point_two_recovery_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_point_two_recovery.v1",
        "seed_voltage_V": 0.15,
        "bias_points_V": points,
        "predictor": "none",
        "maximum_step_V": 0.01,
        "forcing_relaxed": False,
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_point_two_recovery_manifest.json", manifest)
    return manifest


def prepare_one_volt_extension(bundle: Path) -> dict[str, Any]:
    """Extend the paired 0.2 V anchors to 1 V without state extrapolation."""
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    points = [0.25, 0.3, 0.4, 0.5, 0.75, 1.0]
    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        seed = (
            bundle / POINT_TWO_RECOVERY_ROOT / case /
            "states/state_bias_0p200000.h5"
        )
        if not seed.is_file():
            raise PreparationError(
                f"converged 0.2 V one-volt extension seed is missing: {seed}"
            )
        document = build_case(base, case, points[-1])
        document = _replace_output_paths(
            document, f"{OUTPUT_ROOT}/{case}", f"{ONE_VOLT_EXTENSION_ROOT}/{case}"
        )
        document["_comment"] = (
            "Gradual 0.2--1 V device-manifold continuation. The preceding "
            "converged state is used directly at every landmark; secant state "
            "extrapolation remains disabled."
        )
        solver = document["solver"]
        solver["verbose"] = False
        solver["max_iter"] = 20
        solver["handoff"]["newton_max_iter"] = 20
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 4.0
        solver["quasi_fermi_update_limit_mode"] = "uniform_trust_region"
        solver["quasi_fermi_trust_region_growth_factor"] = 2.0
        solver["quasi_fermi_trust_region_max_multiplier"] = 4.0
        solver["quasi_fermi_trust_region_expansion_threshold"] = 0.75
        solver["quasi_fermi_trust_region_shrink_factor"] = 0.5
        solver["quasi_fermi_trust_region_min_multiplier"] = 0.125
        sweep = document["sweep"]
        sweep["bias_points"] = points
        sweep["start"] = points[0]
        sweep["stop"] = points[-1]
        sweep["initial_state_file"] = (
            f"{POINT_TWO_RECOVERY_ROOT}/{case}/states/state_bias_0p200000.h5"
        )
        sweep["initial_step"] = 0.05
        sweep["growth_factor"] = 1.0
        sweep["max_step"] = 0.25
        sweep["continuation"] = {
            "predictor": {"mode": "none"},
            "arclength": {"enabled": False},
        }
        documents[case] = document

    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("one-volt IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / ONE_VOLT_EXTENSION_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        (output / "rejected_states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_corrected_one_volt_extension_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_one_volt_extension.v1",
        "seed_voltage_V": 0.2,
        "bias_points_V": points,
        "predictor": "none",
        "forcing_relaxed": False,
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_one_volt_extension_manifest.json", manifest)
    return manifest


def prepare_high_voltage_extension(bundle: Path) -> dict[str, Any]:
    """Extend paired 1 V anchors to 12 V along the accepted device manifold."""
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    points = [
        1.05, 1.1, 1.2, 1.3, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0,
        6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
    ]
    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        seed = (
            bundle / ONE_VOLT_EXTENSION_ROOT / case /
            "states/state_bias_1p000000.h5"
        )
        if not seed.is_file():
            raise PreparationError(
                f"converged 1 V high-voltage extension seed is missing: {seed}"
            )
        document = build_case(base, case, points[-1])
        document = _replace_output_paths(
            document,
            f"{OUTPUT_ROOT}/{case}",
            f"{HIGH_VOLTAGE_EXTENSION_ROOT}/{case}",
        )
        document["_comment"] = (
            "High-voltage continuation from the paired converged 1 V states. "
            "Every nonlinear solve starts from the preceding accepted device "
            "state; state extrapolation and forcing relaxation are disabled."
        )
        solver = document["solver"]
        solver["verbose"] = False
        solver["max_iter"] = 20
        solver["handoff"]["newton_max_iter"] = 20
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 4.0
        solver["quasi_fermi_update_limit_mode"] = "uniform_trust_region"
        solver["quasi_fermi_trust_region_growth_factor"] = 2.0
        solver["quasi_fermi_trust_region_max_multiplier"] = 4.0
        solver["quasi_fermi_trust_region_expansion_threshold"] = 0.75
        solver["quasi_fermi_trust_region_shrink_factor"] = 0.5
        solver["quasi_fermi_trust_region_min_multiplier"] = 0.125
        sweep = document["sweep"]
        sweep["bias_points"] = points
        sweep["start"] = points[0]
        sweep["stop"] = points[-1]
        sweep["initial_state_file"] = (
            f"{ONE_VOLT_EXTENSION_ROOT}/{case}/states/state_bias_1p000000.h5"
        )
        sweep["initial_step"] = 0.05
        sweep["growth_factor"] = 1.25
        sweep["max_step"] = 0.25
        sweep["continuation"] = {
            "predictor": {"mode": "none"},
            "arclength": {"enabled": False},
        }
        documents[case] = document

    if _normalized_pair(documents["ialmob_off"]) != _normalized_pair(
        documents["ialmob_on"]
    ):
        raise PreparationError("high-voltage IALMob pair is not a strict A/B")

    cases: list[dict[str, str]] = []
    for case, document in documents.items():
        output = bundle / HIGH_VOLTAGE_EXTENSION_ROOT / case
        (output / "states").mkdir(parents=True, exist_ok=True)
        (output / "rejected_states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_corrected_high_voltage_extension_{case}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.corrected_ialmob_high_voltage_extension.v1",
        "seed_voltage_V": 1.0,
        "bias_points_V": points,
        "predictor": "none",
        "maximum_step_V": 0.25,
        "forcing_relaxed": False,
        "cases": cases,
    }
    _write_json(bundle / "corrected_ialmob_high_voltage_extension_manifest.json", manifest)
    return manifest


def prepare_device_corrector_chunk(
    bundle: Path,
    bias_voltage_V: float,
    chunk: int,
) -> dict[str, Any]:
    """Continue Newton on a positive-carrier rejected state at fixed bias."""
    if bias_voltage_V <= 0.0:
        raise PreparationError("corrector bias must be positive")
    if chunk < 1:
        raise PreparationError("corrector chunk must be at least one")
    bundle = bundle.resolve()
    base = _read_json(bundle / "simulation_04_avalanche_activation_1v.json")
    bias_tag = f"{bias_voltage_V:.6f}".replace("-", "m").replace(".", "p")
    cases: list[dict[str, str]] = []
    documents: dict[str, dict[str, Any]] = {}
    for case in ("ialmob_off", "ialmob_on"):
        if chunk == 1:
            source_state = (
                f"{EXTENSION_ROOT}/{case}/rejected_states/"
                f"attempt_1_bias_{bias_tag}_final.csv"
            )
        else:
            source_state = (
                f"{DEVICE_CORRECTOR_ROOT}/bias_{bias_tag}/{case}/"
                f"chunk_{chunk - 1:02d}/rejected_states/"
                f"attempt_1_bias_{bias_tag}_final.csv"
            )
        if not (bundle / source_state).is_file():
            raise PreparationError(f"corrector source state is missing: {source_state}")
        output = (
            f"{DEVICE_CORRECTOR_ROOT}/bias_{bias_tag}/{case}/chunk_{chunk:02d}"
        )
        document = build_case(base, case, bias_voltage_V)
        document = _replace_output_paths(
            document, f"{OUTPUT_ROOT}/{case}", output
        )
        document["_comment"] = (
            "Checkpointed device-manifold corrector. The input is the prior "
            "chunk's positive-carrier rejected final state at the same bias; "
            "nonlinear tolerances and forcing are unchanged."
        )
        document["_ialmob_ablation"].update({
            "corrector_bias_V": bias_voltage_V,
            "corrector_chunk": chunk,
            "corrector_source_state": source_state,
        })
        solver = document["solver"]
        solver["verbose"] = False
        solver["max_iter"] = 10
        solver["handoff"]["newton_max_iter"] = 10
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 4.0
        solver["quasi_fermi_update_limit_mode"] = "uniform_trust_region"
        solver["quasi_fermi_trust_region_growth_factor"] = 2.0
        solver["quasi_fermi_trust_region_max_multiplier"] = 4.0
        solver["quasi_fermi_trust_region_expansion_threshold"] = 0.75
        solver["quasi_fermi_trust_region_shrink_factor"] = 0.5
        solver["quasi_fermi_trust_region_min_multiplier"] = 0.125
        sweep = document["sweep"]
        sweep["bias_points"] = [bias_voltage_V]
        sweep["start"] = bias_voltage_V
        sweep["stop"] = bias_voltage_V
        sweep["initial_state_file"] = source_state
        sweep["continuation"] = {"arclength": {"enabled": False}}
        documents[case] = document

    # The state files differ by design; all configured physics still must be
    # the same controlled mobility delta.
    for case, document in documents.items():
        output = (
            bundle / DEVICE_CORRECTOR_ROOT / f"bias_{bias_tag}" / case /
            f"chunk_{chunk:02d}"
        )
        (output / "states").mkdir(parents=True, exist_ok=True)
        (output / "rejected_states").mkdir(parents=True, exist_ok=True)
        filename = f"simulation_device_corrector_{bias_tag}_{case}_chunk_{chunk:02d}.json"
        _write_json(bundle / filename, document)
        cases.append({"case": case, "config": filename})
    manifest: dict[str, Any] = {
        "schema": "vela.slot_ldmos.device_manifold_corrector.v1",
        "bias_voltage_V": bias_voltage_V,
        "chunk": chunk,
        "maximum_newton_iterations": 10,
        "acceptance": "normal production Newton convergence only",
        "forcing_relaxed": False,
        "quasi_fermi_update_limit_mode": "uniform_trust_region",
        "quasi_fermi_trust_region_growth_factor": 2.0,
        "quasi_fermi_trust_region_max_multiplier": 4.0,
        "quasi_fermi_trust_region_expansion_threshold": 0.75,
        "quasi_fermi_trust_region_shrink_factor": 0.5,
        "quasi_fermi_trust_region_min_multiplier": 0.125,
        "cases": cases,
    }
    _write_json(
        bundle / f"device_corrector_{bias_tag}_chunk_{chunk:02d}_manifest.json",
        manifest,
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--stop-voltage", type=float, default=45.0)
    parser.add_argument("--accelerated-extension", action="store_true")
    parser.add_argument("--dense-extension", action="store_true")
    parser.add_argument("--post-dense-extension", action="store_true")
    parser.add_argument("--point-two-recovery", action="store_true")
    parser.add_argument("--one-volt-extension", action="store_true")
    parser.add_argument("--high-voltage-extension", action="store_true")
    parser.add_argument("--corrector-bias", type=float)
    parser.add_argument("--corrector-chunk", type=int)
    args = parser.parse_args()
    if args.corrector_chunk is not None:
        if args.corrector_bias is None:
            parser.error("--corrector-chunk requires --corrector-bias")
        result = prepare_device_corrector_chunk(
            args.bundle, args.corrector_bias, args.corrector_chunk
        )
    elif args.high_voltage_extension:
        result = prepare_high_voltage_extension(args.bundle)
    elif args.one_volt_extension:
        result = prepare_one_volt_extension(args.bundle)
    elif args.point_two_recovery:
        result = prepare_point_two_recovery(args.bundle)
    elif args.post_dense_extension:
        result = prepare_post_dense_extension(args.bundle, args.stop_voltage)
    elif args.dense_extension:
        result = prepare_dense_low_voltage_extension(args.bundle, args.stop_voltage)
    elif args.accelerated_extension:
        result = prepare_accelerated_extension(args.bundle, args.stop_voltage)
    else:
        result = prepare(args.bundle, args.stop_voltage)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
