#!/usr/bin/env python3
"""Run BVmethods NMOS external-resistor or voltage-to-current extraction.

The script freezes the validated IIC transport/avalanche definitions and only
changes avalanche coupling from postprocess-only to self-consistent, as required
by the two circuit/boundary methods. It never scales mobility or avalanche
parameters.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
RUN = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
VALIDATION = RUN / "vela_validation"
BASE = (
    VALIDATION
    / "btbt_e2_iic_qf_vector_branch_6p5_7p1_20260805/simulation.json"
)
INITIAL = (
    VALIDATION
    / "btbt_e2_adaptive_0_7_20260804/segment_0p0_0p5/"
      "states/accepted_state_bias_0p000000.csv"
)
INITIAL_INNER_V = 5.9
PREBIAS_POINTS_V = [
    0.0, 0.01, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0,
    *[index / 10.0 for index in range(11, 60)],
    5.25,
]
PREBIAS_POINTS_V.sort()
REFERENCES = {
    "external_resistor": 6.379791636301563,
    "voltage_to_current": 6.38318420057198,
}
CURRENT_THRESHOLD_A_PER_UM = 1.0e-4
VOLTAGE_TO_CURRENT_VOLTAGE_POINTS_V = [
    5.9, 5.925, 5.95, 5.9625, 5.975, 5.9875, 6.0
]


def absolute(path: Path) -> str:
    return str(path.resolve())


def configure_common(config: dict[str, Any], output: Path, initial_state: Path) -> None:
    config["output_csv"] = absolute(output / "sweep.csv")
    config["materials_file"] = absolute(
        REPO / "reference_tcad/bvmethods_sentaurus2018/vela/materials_sentaurus2018.json"
    )
    solver = config["solver"]
    solver["verbose"] = True
    impact = solver["impact_ionization"]
    impact.update(
        {
            "model": "van_overstraeten",
            "coupling_mode": "self_consistent",
            "driving_force": "eparallel",
            "generation": "current_density",
            "current_approximation": "nodal_vector_current_reconstructed",
            "current_magnitude_mode": "edge_scalar_abs",
            "source_volume_policy": "genius_truncated",
            "source_mapping_mode": "node_F_node_alpha_node_G",
        }
    )
    sweep = config["sweep"]
    sweep["initial_state_file"] = absolute(initial_state)
    sweep["write_vtk"] = False
    sweep["stop_on_failure"] = True
    sweep["write_state_file"] = absolute(output / "last_state.csv")
    sweep["write_state_every_point_prefix"] = absolute(
        output / "states" / "accepted_state"
    )
    sweep["boundary_control"] = {
        "evaluation_csv": absolute(output / "boundary_control_evaluations.csv"),
        "checkpoint_directory": absolute(output / "boundary_control_checkpoints"),
        "resume": True,
        # Cross-target scalar and state prediction may span up to three proven
        # 25 mV high-field continuation steps. Failed state predictions retain
        # the existing constant-warm-start fallback.
        "predictor_max_step_factor": 3.0,
        "preferred_max_evaluations": 3,
    }
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics.pop("sg_avalanche_edges", None)
    diagnostics.pop("terminal_current_method_compare", None)
    diagnostics.pop("continuity_balance", None)
    if "qf_bounds" in diagnostics:
        diagnostics["qf_bounds"]["csv_file"] = absolute(output / "qf_bounds.csv")
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": absolute(output / "newton_history.csv"),
        "attempts_csv_file": absolute(output / "newton_attempts.csv"),
        "iterations_csv_file": absolute(output / "newton_iterations.csv"),
    }
    carrier = solver.get("carrier_row_convergence", {})
    carrier["diagnostic_csv"] = absolute(output / "carrier_row_convergence.csv")
    carrier["trace_csv"] = absolute(output / "carrier_row_trace.csv")
    # Local row ratios remain an audit because near-zero high-field rows can
    # amplify round-off. Global electron/hole continuity closure is the hard
    # conservation gate for every accepted prebias and boundary-control state.
    solver["global_continuity_closure"] = {
        "mode": "enforce",
        "tolerance": 1.0e-2,
        "source_floor": 1.0e-14,
    }


def configure_external_resistor(config: dict[str, Any], max_points: int | None = None) -> None:
    sweep = config["sweep"]
    sweep["continuation"] = {
        "predictor": {
            "mode": "secant",
            "fields": ["psi", "phin", "phip"],
            "max_extrapolation_ratio": 4.0,
        }
    }
    # Continue forward from the 5.9 V self-consistent state. Earlier 67.5 V and
    # 206 V targets lie below that state's actual 10 Mohm*um load-line outer
    # voltage and would force an unnecessary backward branch traversal.
    # 406 V remains below the 0.1 mA/um extraction threshold while the later
    # targets bracket and pass it.
    outer_points = [406.0, 606.0, 806.0, 1006.0, 1206.0]
    if max_points is not None:
        outer_points = outer_points[:max_points]
    sweep.update(
        {
            "start": outer_points[0],
            "stop": outer_points[-1],
            "step": 200.0,
            "bias_points": outer_points,
            "external_circuit": {
                "mode": "series_resistor",
                "resistance_ohm_um": 1.0e7,
                "current_direction": 1.0,
                "initial_inner_voltage_V": INITIAL_INNER_V,
                # Use the same 25 mV high-field continuation step as the
                # current-boundary method.  The former 6.09959 V failure was
                # traced to a Fermi-Dirac diagnostic-source mismatch and an
                # insufficient continuous-Newton budget, not this step size.
                "max_inner_voltage_step_V": 0.025,
                "residual_tolerance_V": 1.0e-4,
                "voltage_tolerance_V": 1.0e-8,
                "max_bracket_steps": 120,
                "max_iterations": 40,
            },
        }
    )
    sweep.pop("voltage_to_current", None)


def configure_voltage_to_current(
    config: dict[str, Any], voltage_points: list[float] | None = None
) -> None:
    sweep = config["sweep"]
    if voltage_points is None:
        voltage_points = VOLTAGE_TO_CURRENT_VOLTAGE_POINTS_V
    # The voltage phase reaches 6.0 V at roughly 3e-5--4e-5 A/um.  Current
    # control must continue forward from that operating point; lower targets
    # would reverse the branch immediately after the boundary switch.
    current_points = [4.0e-5, 6.0e-5, 1.0e-4]
    sweep["continuation"] = {
        "predictor": {
            "mode": "secant",
            "fields": ["psi", "phin", "phip"],
            "max_extrapolation_ratio": 4.0,
        }
    }
    # Once two current-controlled states exist, allow the scalar secant to
    # propose up to three proven 12.5 mV increments at once.  The DD state is
    # extrapolated over the same interval and falls back to a constant warm
    # start if that prediction is rejected.  Trials at 50 mV remained
    # convergent but approached the unchanged 220-Newton budget above 6.15 V.
    sweep["boundary_control"]["predictor_max_step_factor"] = 3.0
    sweep.update(
        {
            "start": voltage_points[0],
            "stop": voltage_points[-1],
            "step": 0.03,
            "bias_points": voltage_points,
            "voltage_to_current": {
                "switch_voltage_V": 6.0,
                "current_direction": 1.0,
                "current_points_A_per_um": current_points,
                # The same high-field branch that needs 12.5 mV voltage-phase
                # continuation also rejects the first 6.0 -> 6.025 V current
                # bracket trial.  Keep current closure on the proven lattice.
                "max_inner_voltage_step_V": 0.0125,
                "current_tolerance_A_per_um": 1.0e-10,
                "voltage_tolerance_V": 1.0e-8,
                "max_bracket_steps": 120,
                "max_iterations": 40,
            },
        }
    )
    sweep.pop("external_circuit", None)


def prebias_state_path(output: Path, bias_V: float) -> Path:
    label = f"{bias_V:.6f}".replace("-", "m").replace(".", "p")
    return output / "states" / f"accepted_state_bias_{label}.csv"


def configure_prebias(config: dict[str, Any], points_V: list[float]) -> None:
    # This stage constructs only the low-bias handoff branch. Keep the full
    # self-consistent avalanche equations and record carrier-row quality, but
    # do not let the strict per-row gate stall an otherwise converged state.
    # The boundary-method stage is rebuilt from the pristine base and keeps
    # carrier-row mode="enforce".
    carrier = config["solver"].get("carrier_row_convergence", {})
    carrier["mode"] = "report"
    recovery = carrier.get("recovery", {})
    recovery["mode"] = "off"
    sweep = config["sweep"]
    sweep.update(
        {
            "start": points_V[0],
            "stop": points_V[-1],
            "step": 0.1,
            "bias_points": points_V,
            "continuation": {
                "predictor": {
                    "mode": "secant",
                    "fields": ["psi", "phin", "phip"],
                    "max_extrapolation_ratio": 2.0,
                }
            },
        }
    )
    sweep.pop("external_circuit", None)
    sweep.pop("voltage_to_current", None)


def configure_boundary_row_reporting(config: dict[str, Any]) -> None:
    """Keep full Newton equations while making the auxiliary row gate diagnostic."""
    solver = config["solver"]
    # The high-field 6.08709 -> 6.09959 V transition needs 191 continuous
    # Newton iterations to reach the unchanged residual and global-continuity
    # gates.  Restarting the same bias loses the monotone trajectory, so keep a
    # bounded margin in the boundary stage without changing IIC/prebias limits.
    solver["max_iter"] = max(int(solver.get("max_iter", 0)), 220)
    carrier = solver.get("carrier_row_convergence", {})
    carrier["mode"] = "report"
    recovery = carrier.get("recovery", {})
    recovery["mode"] = "off"


def run_runner(config: dict[str, Any], output: Path, runner: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    config_path = output / "simulation.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    command = [str(runner.resolve()), "--config", str(config_path.resolve())]
    with (output / "run.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            cwd=output,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
        return process.wait()


def extract_bv(sweep_csv: Path) -> float:
    with sweep_csv.open(newline="", encoding="utf-8-sig") as handle:
        rows = [row for row in csv.DictReader(handle) if row["converged"] == "1"]
    pairs = [
        (float(row["inner_voltage_V"]), abs(float(row["current_total_A_per_um"])))
        for row in rows
    ]
    for (v0, i0), (v1, i1) in zip(pairs, pairs[1:]):
        if i0 <= CURRENT_THRESHOLD_A_PER_UM <= i1:
            if i1 == i0:
                return v1
            fraction = (CURRENT_THRESHOLD_A_PER_UM - i0) / (i1 - i0)
            return v0 + fraction * (v1 - v0)
    raise RuntimeError("run did not bracket 1e-4 A/um")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--method",
        choices=("external_resistor", "voltage_to_current"),
        default="external_resistor",
    )
    parser.add_argument("--base", type=Path, default=BASE)
    parser.add_argument("--initial-state", type=Path, default=INITIAL)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--prebias-output",
        type=Path,
        help="reuse or create the common self-consistent prebias in this directory",
    )
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--max-points", type=int)
    parser.add_argument(
        "--probe-inner-voltage",
        type=float,
        help="run one self-consistent voltage point without boundary control",
    )
    parser.add_argument(
        "--trace-node",
        type=int,
        action="append",
        default=[],
        help="node id to include in carrier-row trace output (repeatable)",
    )
    parser.add_argument(
        "--probe-max-iter",
        type=int,
        help="override solver.max_iter for a single-point probe",
    )
    args = parser.parse_args()

    if args.output is None:
        args.output = VALIDATION / f"boundary_{args.method}_20260806"
    if args.runner is None:
        executable = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
        args.runner = REPO / "build-release" / executable

    base = json.loads(args.base.read_text(encoding="utf-8-sig"))
    args.output.mkdir(parents=True, exist_ok=True)
    if args.probe_inner_voltage is not None:
        config = copy.deepcopy(base)
        configure_common(config, args.output, args.initial_state)
        configure_boundary_row_reporting(config)
        if args.probe_max_iter is not None:
            config["solver"]["max_iter"] = args.probe_max_iter
        if args.trace_node:
            config["solver"]["carrier_row_convergence"]["trace_nodes"] = (
                args.trace_node
            )
        sweep = config["sweep"]
        sweep.update(
            {
                "start": args.probe_inner_voltage,
                "stop": args.probe_inner_voltage,
                "step": 0.025,
                "bias_points": [args.probe_inner_voltage],
            }
        )
        sweep.pop("external_circuit", None)
        sweep.pop("voltage_to_current", None)
        config["_validation_case"] = {
            "purpose": "BVmethods NMOS self-consistent boundary-method probe",
            "probe_inner_voltage_V": args.probe_inner_voltage,
            "physics_parameter_scaling": "none",
        }
        returncode = run_runner(config, args.output, args.runner)
        if returncode == 0:
            print(args.output / "sweep.csv")
        return returncode

    prebias_output = (
        args.prebias_output if args.prebias_output is not None
        else args.output / "prebias"
    )
    completed_prebias_points = [
        bias for bias in PREBIAS_POINTS_V
        if prebias_state_path(prebias_output, bias).exists()
    ]
    last_prebias_V = completed_prebias_points[-1] if completed_prebias_points else None
    prebias_ready = last_prebias_V == INITIAL_INNER_V
    if not prebias_ready:
        if last_prebias_V is None:
            prebias_points = PREBIAS_POINTS_V
            prebias_initial_state = args.initial_state
        else:
            start_index = PREBIAS_POINTS_V.index(last_prebias_V)
            prebias_points = PREBIAS_POINTS_V[start_index:]
            prebias_initial_state = prebias_state_path(prebias_output, last_prebias_V)
        prebias = copy.deepcopy(base)
        configure_common(prebias, prebias_output, prebias_initial_state)
        configure_prebias(prebias, prebias_points)
        prebias["_validation_case"] = {
            "purpose": "BVmethods NMOS low-bias self-consistent avalanche prebias",
            "terminal_bias_V": INITIAL_INNER_V,
            "physics_parameter_scaling": "none",
            "resumed_from_bias_V": last_prebias_V,
        }
        returncode = run_runner(prebias, prebias_output, args.runner)
        if returncode != 0:
            return returncode

    prebias_state = prebias_state_path(prebias_output, INITIAL_INNER_V)
    boundary_initial_state = prebias_state
    voltage_phase_points: list[float] | None = None
    resumed_voltage_V: float | None = None
    if args.method == "voltage_to_current":
        for voltage_V in reversed(VOLTAGE_TO_CURRENT_VOLTAGE_POINTS_V):
            candidate = prebias_state_path(args.output, voltage_V)
            if candidate.exists():
                boundary_initial_state = candidate
                resumed_voltage_V = voltage_V
                break
        voltage_phase_points = [
            voltage_V for voltage_V in VOLTAGE_TO_CURRENT_VOLTAGE_POINTS_V
            if resumed_voltage_V is None or voltage_V > resumed_voltage_V + 1.0e-12
        ]
        if not voltage_phase_points:
            voltage_phase_points = [VOLTAGE_TO_CURRENT_VOLTAGE_POINTS_V[-1]]
    config = copy.deepcopy(base)
    configure_common(config, args.output, boundary_initial_state)
    configure_boundary_row_reporting(config)
    if args.method == "external_resistor":
        configure_external_resistor(config, args.max_points)
    else:
        configure_voltage_to_current(config, voltage_phase_points)
    config["_validation_case"] = {
        "purpose": f"BVmethods NMOS {args.method} boundary extraction",
        "sentaurus_reference_BV_V": REFERENCES[args.method],
        "acceptance_relative_error": 0.03,
        "current_threshold_A_per_um": CURRENT_THRESHOLD_A_PER_UM,
        "physics_parameter_scaling": "none",
        "resumed_voltage_phase_from_V": resumed_voltage_V,
    }
    returncode = run_runner(config, args.output, args.runner)
    if returncode != 0:
        return returncode

    bv = extract_bv(args.output / "sweep.csv")
    reference = REFERENCES[args.method]
    relative_error = (bv - reference) / reference
    summary = {
        "method": args.method,
        "vela_bv_V": bv,
        "sentaurus_bv_V": reference,
        "delta_V": bv - reference,
        "relative_error": relative_error,
        "pass_3_percent": abs(relative_error) <= 0.03,
        "current_threshold_A_per_um": CURRENT_THRESHOLD_A_PER_UM,
        "mobility_or_avalanche_parameter_fit": False,
    }
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0 if summary["pass_3_percent"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
