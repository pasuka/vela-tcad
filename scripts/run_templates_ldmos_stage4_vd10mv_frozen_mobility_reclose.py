#!/usr/bin/env python3
"""Run a controlled D5 Id-Vd reclose or full-curve solver experiment."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

from audit_templates_ldmos_stage4_vd10mv_newton import load_json, set_contact_bias, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gate-bias", type=float)
    parser.add_argument("--initial-state", type=Path)
    parser.add_argument("--full-curve", action="store_true")
    parser.add_argument(
        "--bias-points",
        type=float,
        nargs="+",
        help="Optional exact drain-bias sequence; overrides the 0/10 mV default.",
    )
    parser.add_argument(
        "--mobility-jacobian",
        choices=("frozen", "live"),
        default="frozen",
    )
    parser.add_argument("--recenter-qf", action="store_true")
    parser.add_argument(
        "--min-step",
        type=float,
        default=1.0e-3,
        help="Minimum internal continuation substep in volts for full-curve retries.",
    )
    parser.add_argument(
        "--max-retries",
        type=int,
        default=12,
        help="Maximum failed continuation retries before rejecting an exact bias point.",
    )
    parser.add_argument("--initial-step", type=float, default=2.5e-3)
    parser.add_argument("--growth-factor", type=float, default=1.0)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    config = load_json(args.base_config.resolve())
    config["solver"]["mobility"]["jacobian_field_derivatives"] = (
        args.mobility_jacobian == "live"
    )
    config["solver"]["quasi_fermi_recenter_on_initial_state"] = (
        args.recenter_qf
    )
    config["solver"]["max_iter"] = 80
    if args.gate_bias is not None:
        set_contact_bias(config, "gate", args.gate_bias)
    if args.initial_state is not None:
        initial_state = args.initial_state.resolve()
        if not initial_state.is_file():
            raise FileNotFoundError(initial_state)
        config["sweep"]["initial_state_file"] = str(initial_state)
    drain_start = args.bias_points[0] if args.bias_points else 0.0
    set_contact_bias(config, "drain", drain_start)
    if args.min_step <= 0.0:
        raise ValueError("--min-step must be positive")
    if args.max_retries < 0:
        raise ValueError("--max-retries must be non-negative")
    if args.initial_step < args.min_step:
        raise ValueError("--initial-step must be at least --min-step")
    if args.growth_factor < 1.0:
        raise ValueError("--growth-factor must be at least one")
    config["sweep"]["initial_step"] = args.initial_step
    config["sweep"]["growth_factor"] = args.growth_factor
    if args.bias_points:
        if len(args.bias_points) < 2:
            raise ValueError("--bias-points requires at least two values")
        if any(right <= left for left, right in zip(
                args.bias_points, args.bias_points[1:])):
            raise ValueError("--bias-points must be strictly increasing")
        config["sweep"]["bias_points"] = args.bias_points
        config["sweep"]["start"] = args.bias_points[0]
        config["sweep"]["stop"] = args.bias_points[-1]
        config["sweep"]["step"] = args.bias_points[1] - args.bias_points[0]
        config["sweep"]["min_step"] = args.min_step
        config["sweep"]["max_retries"] = args.max_retries
    elif not args.full_curve:
        config["sweep"]["bias_points"] = [0.0, 0.01]
        config["sweep"]["start"] = 0.0
        config["sweep"]["stop"] = 0.01
        config["sweep"]["step"] = 0.01
    else:
        # Explicit bias_points are the reporting grid, while failed transfers
        # may still need smaller internal continuation steps.  Do not inherit
        # min_step from the coarse 1.333 V nominal output spacing.
        config["sweep"]["min_step"] = args.min_step
        config["sweep"]["max_retries"] = args.max_retries
    config["output_csv"] = str(output / "curve.csv")
    config["sweep"]["write_state_file"] = str(output / "state.csv")
    config["sweep"]["write_state_every_point_prefix"] = str(output / "point")
    config["sweep"]["write_state_every_accepted_step_prefix"] = str(
        output / "accepted_step"
    )
    diagnostics = config["sweep"]["diagnostics"]
    diagnostics["terminal_balance"]["csv_file"] = str(output / "terminal_balance.csv")
    diagnostics["srh_balance"]["csv_file"] = str(output / "srh_balance.csv")
    diagnostics["newton_history"]["csv_file"] = str(output / "newton_history.csv")
    diagnostics["newton_history"]["attempts_csv_file"] = str(output / "newton_attempts.csv")
    diagnostics["newton_history"]["iterations_csv_file"] = str(output / "newton_iterations.csv")
    diagnostics["newton_history"]["rejected_state_directory"] = str(output / "rejected_states")
    config_path = output / "reclose.json"
    write_json(config_path, config)
    completed = subprocess.run(
        [str(args.runner.resolve()), "--config", str(config_path)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    (output / "reclose.log").write_text(completed.stdout, encoding="utf-8")
    summary = {
        "schema": "vela.templates_ldmos.stage4_idvd_reclose.v2",
        "return_code": completed.returncode,
        "config": str(config_path),
        "curve_exists": (output / "curve.csv").is_file(),
        "state_exists": (output / "state.csv").is_file(),
        "full_curve": args.full_curve,
        "mobility_jacobian": args.mobility_jacobian,
        "quasi_fermi_recenter_on_initial_state": args.recenter_qf,
        "requested_bias_points": args.bias_points,
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
