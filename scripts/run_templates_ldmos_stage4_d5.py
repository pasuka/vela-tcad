#!/usr/bin/env python3
"""Run the exact-point Templates/LDMOS classical isothermal D5 Id-Vd curves."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any, Sequence


IDVD_MIN_INTERNAL_STEP_V = 1.0e-3
IDVD_MAX_RETRIES = 12
IDVD_INITIAL_INTERNAL_STEP_V = 2.5e-3
# A 2x-growth control repeatedly overshot the qualified 2.5--3.125 mV
# transfer size and spent more time shrinking than the fixed-step path.
IDVD_GROWTH_FACTOR = 1.0


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    return env


def read_points(path: Path) -> list[float]:
    with path.open(newline="", encoding="utf-8") as stream:
        raw_points = [float(row["bias_V"]) for row in csv.DictReader(stream)]
    points: list[float] = []
    for point in raw_points:
        if not points or point != points[-1]:
            points.append(point)
    duplicate_count = len(raw_points) - len(points)
    if duplicate_count not in {0, 1}:
        raise ValueError(f"unexpected duplicate CurrentPlot points: {duplicate_count}")
    if len(points) != 31 or points[0] != 0.0 or any(
            right <= left for left, right in zip(points, points[1:])):
        raise ValueError("Stage-4 D5 reference must contain 31 strictly increasing points from zero")
    return points


def redirect_diagnostics(config: dict[str, Any], output: Path, stem: str) -> None:
    sweep = config["sweep"]
    diagnostics = sweep.setdefault("diagnostics", {})
    if "terminal_balance" in diagnostics:
        diagnostics["terminal_balance"]["csv_file"] = str(
            (output / f"{stem}_terminal_balance.csv").resolve())
    if "srh_balance" in diagnostics:
        diagnostics["srh_balance"]["csv_file"] = str(
            (output / f"{stem}_srh_balance.csv").resolve())
    if "newton_history" in diagnostics:
        history = diagnostics["newton_history"]
        history["csv_file"] = str((output / f"{stem}_newton_history.csv").resolve())
        history["attempts_csv_file"] = str((output / f"{stem}_newton_attempts.csv").resolve())
        history["iterations_csv_file"] = str((output / f"{stem}_newton_iterations.csv").resolve())
        history["rejected_state_directory"] = str((output / f"{stem}_rejected_states").resolve())


def assert_d5_contract(config: dict[str, Any]) -> None:
    serialized = json.dumps(config.get("solver", {})).lower()
    for forbidden in ("ialmob", "quantum_potential", "impact_ionization\": {\"model\": \"okuto"):
        if forbidden in serialized:
            raise ValueError(f"D5 config contains forbidden physics: {forbidden}")
    if config["solver"].get("impact_ionization", {}).get("model") != "none":
        raise ValueError("D5 requires avalanche disabled")
    if config.get("discretization", {}).get("poisson_charge_volume_policy") != "material_local":
        raise ValueError("D5 requires qualified material-local Poisson charge volume")
    geometry = config.get("mesh_geometry", {})
    if geometry.get("carrier_transport_couple_profile") != "templates_ldmos_external_averagebox":
        raise ValueError("D5 requires qualified external AverageBox carrier couples")
    if "predictor" in config["sweep"] or "continuation" in config["sweep"]:
        raise ValueError("D5 qualification requires predictor and continuation disabled")


def enable_drain_branch_guard(config: dict[str, Any]) -> None:
    config["solver"]["contact_majority_qf_branch_guard_contacts"] = ["drain"]
    config["solver"]["contact_majority_qf_branch_drop_limit_V"] = 5.0e-11


def make_gate8_prebias(base: dict[str, Any], initial_state: Path,
                       output: Path) -> dict[str, Any]:
    config = deepcopy(base)
    for contact in config["contacts"]:
        contact["bias"] = 0.1 if contact["name"] == "drain" else (
            5.0 if contact["name"] == "gate" else 0.0)
    config["output_csv"] = str((output / "gate8_prebias.csv").resolve())
    config["sweep"].update({
        "contact": "gate", "current_contact": "drain", "start": 5.0,
        "stop": 8.0, "step": 0.5,
        "bias_points": [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0],
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str((output / "gate8_prebias_state.h5").resolve()),
    })
    config["sweep"].pop("write_state_every_point_prefix", None)
    redirect_diagnostics(config, output, "gate8_prebias")
    assert_d5_contract(config)
    return config


def make_gate_prebias_vd0(base: dict[str, Any], initial_state: Path,
                          output: Path) -> dict[str, Any]:
    config = deepcopy(base)
    for contact in config["contacts"]:
        contact["bias"] = 0.0
    points = [0.5 * index for index in range(17)]
    config["output_csv"] = str((output / "gate_prebias_vd0.csv").resolve())
    config["sweep"].update({
        "contact": "gate", "current_contact": "drain", "start": 0.0,
        "stop": 8.0, "step": 0.5, "bias_points": points,
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str((output / "gate_prebias_vd0_state.h5").resolve()),
        "write_state_every_point_prefix": str((output / "gate_prebias_vd0_point").resolve()),
    })
    redirect_diagnostics(config, output, "gate_prebias_vd0")
    assert_d5_contract(config)
    return config


def make_idvd(base: dict[str, Any], gate_V: float, initial_state: Path,
              points: list[float], output: Path,
              drain_branch_guard: bool = True) -> dict[str, Any]:
    config = deepcopy(base)
    # Rebase the internally stored QF increments at every accepted warm-start
    # state. The absolute QF values are unchanged; this preserves sub-ULP
    # Newton corrections once the drain sweep leaves equilibrium.
    config["solver"]["quasi_fermi_recenter_on_initial_state"] = True
    # Use a lagged HFS field in the nonlinear Jacobian. The converged residual
    # still evaluates the live HFS mobility; this only selects the qualified,
    # substantially cheaper quasi-Newton linearisation for production.
    config["solver"].setdefault("mobility", {})[
        "jacobian_field_derivatives"
    ] = False
    if drain_branch_guard:
        enable_drain_branch_guard(config)
    stem = f"d5_idvd_vg{int(gate_V)}"
    for contact in config["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = gate_V
        elif contact["name"] == "drain":
            contact["bias"] = 0.0
        else:
            contact["bias"] = 0.0
    config["output_csv"] = str((output / f"{stem}.csv").resolve())
    startup_bridge = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8]
    solver_points = sorted(set(points + startup_bridge))
    config["sweep"].update({
        "contact": "drain", "current_contact": "drain", "start": points[0],
        "stop": points[-1], "step": points[1] - points[0],
        "bias_points": solver_points,
        # bias_points defines the exact reporting grid.  Keep retry subdivision
        # independent of its coarse CurrentPlot spacing so a failed transfer
        # can qualify through the 10--100 mV startup region without weakening
        # the frozen nonlinear residual ceilings.
        "min_step": IDVD_MIN_INTERNAL_STEP_V,
        "max_retries": IDVD_MAX_RETRIES,
        "initial_step": IDVD_INITIAL_INTERNAL_STEP_V,
        "growth_factor": IDVD_GROWTH_FACTOR,
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str((output / f"{stem}_state.h5").resolve()),
        "write_state_every_point_prefix": str((output / f"{stem}_point").resolve()),
    })
    redirect_diagnostics(config, output, stem)
    assert_d5_contract(config)
    return config


def make_drain_zero_prebias(base: dict[str, Any], gate_V: float,
                            initial_state: Path, output: Path) -> dict[str, Any]:
    config = deepcopy(base)
    enable_drain_branch_guard(config)
    stem = f"drain_zero_vg{int(gate_V)}"
    for contact in config["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = gate_V
        elif contact["name"] == "drain":
            contact["bias"] = 0.1
        else:
            contact["bias"] = 0.0
    config["output_csv"] = str((output / f"{stem}.csv").resolve())
    config["sweep"].update({
        "contact": "drain", "current_contact": "drain", "start": 0.1,
        "stop": 0.0, "step": -0.025,
        "bias_points": [0.1, 0.075, 0.05, 0.025, 0.0],
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str((output / f"{stem}_state.h5").resolve()),
    })
    config["sweep"].pop("write_state_every_point_prefix", None)
    redirect_diagnostics(config, output, stem)
    assert_d5_contract(config)
    return config


def execute(runner: Path, config: Path, log: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve()), "--log", str(log.resolve())],
        text=True, capture_output=True, env=runner_environment(), check=False)
    config.with_suffix(".stdout.txt").write_text(completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return {"config": str(config.resolve()), "config_sha256": sha256(config),
            "return_code": completed.returncode,
            "status": "pass" if completed.returncode == 0 else "fail"}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--vg4-initial-state", type=Path, required=True)
    parser.add_argument("--vg5-initial-state", type=Path, required=True)
    parser.add_argument("--equilibrium-state", type=Path)
    parser.add_argument("--prebiased-dir", type=Path)
    parser.add_argument("--disable-drain-branch-guard", action="store_true")
    parser.add_argument("--max-iter", type=int)
    parser.add_argument("--vg4-reference", type=Path, required=True)
    parser.add_argument("--vg8-reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    base = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    if args.max_iter is not None:
        if args.max_iter <= 0:
            raise ValueError("--max-iter must be positive")
        base["solver"]["max_iter"] = args.max_iter
    assert_d5_contract(base)
    points4, points8 = read_points(args.vg4_reference), read_points(args.vg8_reference)
    if points4 != points8:
        raise ValueError("Vg=4 and Vg=8 D5 references must share exact drain points")
    if args.prebiased_dir:
        prebiased = args.prebiased_dir.resolve()
        required = [prebiased / f"gate_prebias_vd0_point_bias_{gate}p000000.h5"
                    for gate in (4, 8)]
        if not all(path.is_file() for path in required):
            raise FileNotFoundError("--prebiased-dir lacks Vg=4/8, Vd=0 checkpoint states")
        runs = [{"status": "pass", "return_code": 0,
                 "classification": "reused_hashed_prebias_checkpoint",
                 "state_sha256": {path.name: sha256(path) for path in required}}]
    elif args.equilibrium_state:
        gate_prebias = make_gate_prebias_vd0(base, args.equilibrium_state, output)
        gate_prebias_path = output / "gate_prebias_vd0.json"
        gate_prebias_path.write_text(
            json.dumps(gate_prebias, indent=2) + "\n", encoding="utf-8")
        runs = [execute(args.runner, gate_prebias_path,
                        output / "gate_prebias_vd0.log")]
    else:
        gate8 = make_gate8_prebias(base, args.vg5_initial_state, output)
        gate8_path = output / "gate8_prebias.json"
        gate8_path.write_text(json.dumps(gate8, indent=2) + "\n", encoding="utf-8")
        runs = [execute(args.runner, gate8_path, output / "gate8_prebias.log")]
    if runs[-1]["return_code"] != 0:
        (output / "run_summary.json").write_text(json.dumps({"runs": runs}, indent=2) + "\n")
        return runs[-1]["return_code"]
    zero_states: dict[float, Path] = {}
    if args.prebiased_dir and all(item["status"] == "pass" for item in runs):
        prebiased = args.prebiased_dir.resolve()
        zero_states = {
            4.0: prebiased / "gate_prebias_vd0_point_bias_4p000000.h5",
            8.0: prebiased / "gate_prebias_vd0_point_bias_8p000000.h5",
        }
    elif args.equilibrium_state and all(item["status"] == "pass" for item in runs):
        zero_states = {
            4.0: output / "gate_prebias_vd0_point_bias_4p000000.h5",
            8.0: output / "gate_prebias_vd0_point_bias_8p000000.h5",
        }
    elif all(item["status"] == "pass" for item in runs):
        for gate, initial in ((4.0, args.vg4_initial_state),
                              (8.0, output / "gate8_prebias_state.h5")):
            prebias = make_drain_zero_prebias(base, gate, initial, output)
            prebias_path = output / f"drain_zero_vg{int(gate)}.json"
            prebias_path.write_text(json.dumps(prebias, indent=2) + "\n", encoding="utf-8")
            runs.append(execute(args.runner, prebias_path,
                                output / f"drain_zero_vg{int(gate)}.log"))
            if runs[-1]["return_code"] != 0:
                break
            zero_states[gate] = output / f"drain_zero_vg{int(gate)}_state.h5"
    if all(item["status"] == "pass" for item in runs):
        for gate in (4.0, 8.0):
            config = make_idvd(
                base, gate, zero_states[gate], points4, output,
                drain_branch_guard=not args.disable_drain_branch_guard)
            path = output / f"d5_idvd_vg{int(gate)}.json"
            path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
            runs.append(execute(args.runner, path, output / f"d5_idvd_vg{int(gate)}.log"))
            if runs[-1]["return_code"] != 0:
                break
    summary = {
        "schema": "vela.templates_ldmos.stage4_d5_run.v1",
        "physics": "isothermal classical DD, HighFieldSaturation retained, IALMob/QP/avalanche off",
        "predictor": "disabled",
        "exact_bias_point_count": len(points4),
        "runs": runs,
    }
    (output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if all(item["status"] == "pass" for item in runs) else 2


if __name__ == "__main__":
    raise SystemExit(main())
