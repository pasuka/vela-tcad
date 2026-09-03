#!/usr/bin/env python3
"""Reclose only the Vg=4 V, Vd=0 -> 10 mV D5 transfer with frozen HFS Jacobian."""

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
    parser.add_argument("--full-curve", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    config = load_json(args.base_config.resolve())
    config["solver"]["mobility"]["jacobian_field_derivatives"] = False
    config["solver"]["max_iter"] = 80
    set_contact_bias(config, "drain", 0.0)
    if not args.full_curve:
        config["sweep"]["bias_points"] = [0.0, 0.01]
        config["sweep"]["start"] = 0.0
        config["sweep"]["stop"] = 0.01
        config["sweep"]["step"] = 0.01
    config["output_csv"] = str(output / "curve.csv")
    config["sweep"]["write_state_file"] = str(output / "state.csv")
    config["sweep"]["write_state_every_point_prefix"] = str(output / "point")
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
        "schema": "vela.templates_ldmos.stage4_vd10mv_frozen_mobility_reclose.v1",
        "return_code": completed.returncode,
        "config": str(config_path),
        "curve_exists": (output / "curve.csv").is_file(),
        "state_exists": (output / "state.csv").is_file(),
        "full_curve": args.full_curve,
    }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
