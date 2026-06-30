#!/usr/bin/env python3
"""Run the PN2D A2/B1.05 avalanche internal source current audit."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import run_vanoverstraeten_A2_B_scale_full_bv_matrix as bmatrix
import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
CASE_ID = "BV-A2-B1p05-internal-source-current-audit"
A_SCALE = 2.0
B_SCALE = 1.05
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]


def default_runner() -> Path:
    exe = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
    build = REPO / "build" / exe
    if build.exists():
        return build
    return REPO / "build-release" / exe


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=(
            REPO
            / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
            / "coarse_previous_full20/simulation_coarse_previous_full20_aligned.json"
        ),
    )
    parser.add_argument("--runner", type=Path, default=default_runner())
    parser.add_argument("--out-dir", type=Path, default=REPO / "build/diagnostics")
    parser.add_argument("--bias-points", default=",".join(f"{bias:g}" for bias in DEFAULT_BIASES))
    parser.add_argument("--case-id", default=CASE_ID)
    parser.add_argument("--a-scale", type=float, default=A_SCALE)
    parser.add_argument("--b-scale", type=float, default=B_SCALE)
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def build_audit_config(base_config: Path, out_dir: Path, biases: list[float], case_id: str, a_scale: float, b_scale: float) -> Path:
    case_root = out_dir / "avalanche_internal_source_current_audit_case"
    case_dir = case_root / case_id
    config_path = bmatrix.build_case_config(base_config, case_dir, case_id, b_scale, biases)
    config = bv.read_json(config_path)

    solver = config.setdefault("solver", {})
    impact = solver.setdefault("impact_ionization", {})
    if isinstance(impact, str):
        impact = {"model": impact}
        solver["impact_ionization"] = impact
    impact["model"] = "van_overstraeten"
    impact["parameter_set"] = "default"
    impact["driving_force"] = "quasi_fermi_gradient"
    impact["generation"] = "current_density"
    impact.setdefault("current_approximation", "grad_qf")
    impact["A_scale"] = a_scale
    impact["B_scale"] = b_scale

    sweep = config.setdefault("sweep", {})
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["avalanche_internal_source_current_audit"] = {
        "enabled": True,
        "csv_file": str((out_dir / "avalanche_internal_source_current_audit.csv").resolve()),
        "summary_file": str((out_dir / "avalanche_internal_source_current_audit_summary.md").resolve()),
    }
    config["_avalanche_internal_source_current_audit_metadata"] = {
        "case_id": case_id,
        "parameter_set": "default",
        "driving_force": "GradQuasiFermi",
        "A_scale": a_scale,
        "B_scale": b_scale,
        "csv_file": diagnostics["avalanche_internal_source_current_audit"]["csv_file"],
        "summary_file": diagnostics["avalanche_internal_source_current_audit"]["summary_file"],
    }

    bv.write_json(config_path, config)
    return config_path


def main() -> int:
    args = parse_args()
    biases = bv.parse_biases(args.bias_points)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    config_path = build_audit_config(args.base_config, args.out_dir, biases, args.case_id, args.a_scale, args.b_scale)

    run_return_code = 0
    if not args.skip_run:
        run_return_code = bv.run_case(
            args.runner,
            config_path,
            args.out_dir / "avalanche_internal_source_current_audit_case" / args.case_id,
        )
        if run_return_code != 0:
            return run_return_code

    print(args.out_dir / "avalanche_internal_source_current_audit.csv")
    print(args.out_dir / "avalanche_internal_source_current_audit_summary.md")
    print(config_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
