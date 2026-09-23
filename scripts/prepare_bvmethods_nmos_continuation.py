#!/usr/bin/env python3
"""Prepare a BVmethods NMOS pseudo-arclength run without changing physics.

The input is an already-qualified voltage or voltage-to-current Vela deck.
This tool retains its mesh, materials, contacts, and solver block verbatim and
replaces only the sweep controller with a restartable pseudo-arclength branch.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(text)


def prepare_config(
    base: dict[str, Any],
    initial_state: Path,
    initial_secant_state: Path,
    output_dir: Path,
    *,
    start_V: float = 6.056459,
    initial_secant_bias_V: float = 6.0,
    stop_V: float = 6.5,
    state_weight: float = 0.0,
    max_corrector_iterations: int = 20,
    corrector_tolerance: float = 1.0e-8,
    max_step_retries: int = 12,
    gummel_max_iter: int | None = None,
    arclength_initial_step: float = 0.01,
    arclength_min_step: float = 1.0e-6,
    arclength_max_step: float = 0.1,
    max_parameter_update: float = 0.02,
    source_jacobian: str | None = None,
    line_search_relative_increase_tolerance: float = 0.0,
) -> dict[str, Any]:
    if base.get("simulation_type") != "dc_sweep":
        raise ValueError("base config must use simulation_type=dc_sweep")
    solver = base.get("solver", {})
    if solver.get("method") not in {"newton", "gummel_newton"}:
        raise ValueError("pseudo-arclength requires a Newton-capable base config")
    impact = solver.get("impact_ionization", {})
    if impact.get("coupling_mode") != "self_consistent":
        raise ValueError("BVmethods continuation requires self-consistent avalanche")
    if not initial_state.is_file():
        raise FileNotFoundError(initial_state)
    if not initial_secant_state.is_file():
        raise FileNotFoundError(initial_secant_state)
    if initial_secant_bias_V >= start_V:
        raise ValueError("initial secant bias must be below the start voltage")
    if stop_V <= start_V:
        raise ValueError("stop voltage must exceed start voltage")
    if state_weight < 0.0:
        raise ValueError("state weight must be non-negative")
    if max_corrector_iterations <= 0:
        raise ValueError("max corrector iterations must be positive")
    if corrector_tolerance <= 0.0:
        raise ValueError("corrector tolerance must be positive")
    if max_step_retries < 0:
        raise ValueError("max step retries must be non-negative")
    if gummel_max_iter is not None and gummel_max_iter < 0:
        raise ValueError("Gummel max iterations must be non-negative")
    if not 0.0 < arclength_min_step <= arclength_initial_step <= arclength_max_step:
        raise ValueError("arclength steps must satisfy 0 < min <= initial <= max")
    if max_parameter_update <= 0.0:
        raise ValueError("maximum parameter update must be positive")
    if source_jacobian not in {None, "frozen", "finite_difference"}:
        raise ValueError("source Jacobian must be frozen or finite_difference")
    if line_search_relative_increase_tolerance < 0.0:
        raise ValueError("line-search relative increase tolerance must be non-negative")

    output_dir = output_dir.resolve()
    config = copy.deepcopy(base)
    if gummel_max_iter is not None:
        config["solver"].setdefault("handoff", {})["gummel_max_iter"] = (
            gummel_max_iter
        )
    config["output_csv"] = str(output_dir / "sweep.csv")
    config["sweep"] = {
        "mode": "bv_reverse",
        "start": start_V,
        "stop": stop_V,
        "step": 0.01,
        "initial_step": 0.01,
        "min_step": 1.0e-8,
        "max_step": 0.1,
        "growth_factor": 1.25,
        "shrink_factor": 0.5,
        "max_retries": 20,
        "stop_on_failure": True,
        "contact": "drain",
        "current_contact": "drain",
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str(output_dir / "last_state.h5"),
        "write_state_every_point_prefix": str(output_dir / "states" / "state"),
        "write_vtk": False,
        "breakdown": {
            "max_electric_field_V_per_m": 1.0e12,
            "current_jump_ratio": 1.0e12,
            "non_convergence": False,
        },
        "diagnostics": {
            "qf_bounds": {
                "enabled": True,
                "mode": "reject_and_recover",
                "margin_V": 1.0,
                "min_carrier_density_m3": 1.0e6,
                "csv_file": str(output_dir / "qf_bounds.csv"),
            },
            "newton_history": {
                "enabled": True,
                "csv_file": str(output_dir / "newton_history.csv"),
                "attempts_csv_file": str(output_dir / "newton_attempts.csv"),
                "iterations_csv_file": str(output_dir / "newton_iterations.csv"),
            },
        },
        "continuation": {
            "arclength": {
                "enabled": True,
                "predictor": "tangent",
                "initial_step": arclength_initial_step,
                "min_step": arclength_min_step,
                "max_step": arclength_max_step,
                "growth_factor": 1.2,
                "shrink_factor": 0.5,
                "max_corrector_iterations": max_corrector_iterations,
                "corrector_tolerance": corrector_tolerance,
                "max_step_retries": max_step_retries,
                "parameter_scale": 1.0,
                # Zero selects the continuation core's mesh-independent 1/N
                # state metric.  An explicit nonzero value remains available
                # for controlled metric-sensitivity studies.
                "state_weight": state_weight,
                "damping_factor": 0.5,
                "max_line_search_steps": 12,
                "line_search_relative_increase_tolerance":
                    line_search_relative_increase_tolerance,
                "max_parameter_update": max_parameter_update,
                "bias_finite_difference_step_V": 1.0e-4,
                "initial_secant_state_file":
                    str(initial_secant_state.resolve()),
                "initial_secant_bias_V": initial_secant_bias_V,
            }
        },
    }
    if source_jacobian is not None:
        config["sweep"]["continuation"]["arclength"]["source_jacobian"] = (
            source_jacobian
        )
    config["_validation_case"] = {
        "purpose": "BVmethods NMOS non-transient continuation closure",
        "source_config": "physics copied verbatim from supplied base config",
        "sentaurus_reference_BV_V": 6.383727168968036,
        "current_threshold_A_per_um": 1.0e-4,
        "acceptance_relative_error": 0.02,
        "physics_parameter_scaling": "none",
    }
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--initial-state", type=Path, required=True)
    parser.add_argument("--initial-secant-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-V", type=float, default=6.056459)
    parser.add_argument("--initial-secant-bias-V", type=float, default=6.0)
    parser.add_argument("--stop-V", type=float, default=6.5)
    parser.add_argument("--state-weight", type=float, default=0.0)
    parser.add_argument("--max-corrector-iterations", type=int, default=20)
    parser.add_argument("--corrector-tolerance", type=float, default=1.0e-8)
    parser.add_argument("--max-step-retries", type=int, default=12)
    parser.add_argument("--gummel-max-iter", type=int)
    parser.add_argument("--arclength-initial-step", type=float, default=0.01)
    parser.add_argument("--arclength-min-step", type=float, default=1.0e-6)
    parser.add_argument("--arclength-max-step", type=float, default=0.1)
    parser.add_argument("--max-parameter-update", type=float, default=0.02)
    parser.add_argument(
        "--source-jacobian", choices=("frozen", "finite_difference"))
    parser.add_argument(
        "--line-search-relative-increase-tolerance", type=float, default=0.0)
    args = parser.parse_args()

    base_config = args.base_config.resolve()
    base = json.loads(base_config.read_text(encoding="utf-8"))
    for key in ("mesh_file", "materials_file", "node_doping_file"):
        if key in base:
            path = Path(base[key])
            if not path.is_absolute():
                base[key] = str((base_config.parent / path).resolve())
    config = prepare_config(
        base, args.initial_state, args.initial_secant_state, args.output_dir,
        start_V=args.start_V,
        initial_secant_bias_V=args.initial_secant_bias_V,
        stop_V=args.stop_V,
        state_weight=args.state_weight,
        max_corrector_iterations=args.max_corrector_iterations,
        corrector_tolerance=args.corrector_tolerance,
        max_step_retries=args.max_step_retries,
        gummel_max_iter=args.gummel_max_iter,
        arclength_initial_step=args.arclength_initial_step,
        arclength_min_step=args.arclength_min_step,
        arclength_max_step=args.arclength_max_step,
        max_parameter_update=args.max_parameter_update,
        source_jacobian=args.source_jacobian,
        line_search_relative_increase_tolerance=
            args.line_search_relative_increase_tolerance)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "simulation.json"
    write_text_lf(output, json.dumps(config, indent=2) + "\n")
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
