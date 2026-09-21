#!/usr/bin/env python3
"""Run frozen-state mobility-field Jacobian A/B at the D5 10 mV failure."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Any

from audit_templates_ldmos_stage4_vd10mv_newton import (
    STATE_LABELS, load_json, set_contact_bias, status_from_stdout, write_json,
)


def run(runner: Path, config_path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=False)
    config_path.with_suffix(".log").write_text(completed.stdout, encoding="utf-8")
    status = status_from_stdout(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(json.dumps(status, indent=2))
    return status


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--states-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    base = load_json(args.base_config.resolve())
    base.pop("sweep", None)
    set_contact_bias(base, "drain", 0.01)
    summary: dict[str, Any] = {
        "schema": "vela.templates_ldmos.stage4_vd10mv_mobility_jacobian_ab.v1",
        "states": {},
    }
    for label in STATE_LABELS:
        state = args.states_dir.resolve() / f"attempt_2_bias_0p010000_{label}.csv"
        summary["states"][label] = {}
        for variant, enabled in (("live", True), ("frozen", False)):
            config = json.loads(json.dumps(base))
            config["simulation_type"] = "newton_step_probe"
            config["state_file"] = str(state)
            config["solver"]["mobility"]["jacobian_field_derivatives"] = enabled
            config["output_csv"] = str(output / f"{label}_{variant}.csv")
            config_path = output / f"{label}_{variant}.json"
            write_json(config_path, config)
            summary["states"][label][variant] = run(args.runner.resolve(), config_path)
    for label in STATE_LABELS:
        live = summary["states"][label]["live"]
        frozen = summary["states"][label]["frozen"]
        summary["states"][label]["comparison"] = {
            "raw_step_norm_ratio_frozen_over_live":
                frozen["raw_step_norm"] / max(live["raw_step_norm"], 1.0e-300),
            "trial_phin_ratio_frozen_over_live":
                frozen["trial_block_residuals"]["phin"] /
                max(live["trial_block_residuals"]["phin"], 1.0e-300),
            "trial_psi_ratio_frozen_over_live":
                frozen["trial_block_residuals"]["psi"] /
                max(live["trial_block_residuals"]["psi"], 1.0e-300),
        }
    write_json(output / "summary.json", summary)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
