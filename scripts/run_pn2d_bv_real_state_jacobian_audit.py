#!/usr/bin/env python
"""Run PN2D BV Jacobian block audits on real sweep restart states."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple


REPO = Path(__file__).resolve().parents[1]
EXE = ".exe" if sys.platform.startswith("win") else ""
DEFAULT_BASE_CONFIG = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "vela"
    / "simulation_bv_minus20_avaljac.json"
)
DEFAULT_RUNNER = REPO / "build-release" / f"vela_example_runner{EXE}"
DEFAULT_OUT_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reports"
    / "bv_real_state_jacobian_audit"
)
DEFAULT_BIASES = [-11.5, -13.2, -18.0, -19.0, -20.0]


class BiasBlockReport(NamedTuple):
    bias: float
    block_csv: Path


def bias_token(bias: float) -> str:
    token = f"{abs(bias):.6f}".replace(".", "p")
    return ("m" if bias < 0 else "") + token


def state_path(out_dir: Path, bias: float) -> Path:
    return out_dir / "states" / f"bv_state_bias_{bias_token(bias)}.csv"


def block_path(out_dir: Path, bias: float) -> Path:
    return out_dir / "blocks" / f"jacobian_blocks_{bias_token(bias)}.csv"


def resolve_repo_path(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (REPO / path).resolve()


def write_combined_report(output: Path, reports: list[BiasBlockReport]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "bias_V",
                "block",
                "analytic_norm",
                "fd_norm",
                "diff_norm",
                "rel_diff",
            ],
        )
        writer.writeheader()
        for report in reports:
            with report.block_csv.open(newline="", encoding="utf-8") as block_handle:
                for row in csv.DictReader(block_handle):
                    merged = dict(row)
                    merged["bias_V"] = f"{report.bias:g}"
                    writer.writerow(merged)


def absolutize_base_paths(cfg: dict[str, object], base_config: Path) -> None:
    base_dir = base_config.resolve().parent
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        value = cfg.get(key)
        if isinstance(value, str) and value:
            path = Path(value)
            if not path.is_absolute():
                cfg[key] = str(base_dir / path)


def derive_snapshot_config(base_config: Path, out_dir: Path) -> Path:
    base_config = resolve_repo_path(base_config)
    out_dir = resolve_repo_path(out_dir)
    cfg = json.loads(base_config.read_text(encoding="utf-8"))
    absolutize_base_paths(cfg, base_config)
    sweep = cfg.setdefault("sweep", {})
    sweep.pop("bias_points", None)
    sweep["csv_file"] = str(out_dir / "snapshot_sweep.csv")
    sweep["write_state_every_point_prefix"] = str(out_dir / "states" / "bv_state")
    sweep["write_vtk"] = False
    cfg_path = out_dir / "snapshot_sweep.json"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return cfg_path


def _sweep_contact_name(cfg: dict[str, object]) -> str:
    sweep = cfg.get("sweep", {})
    if isinstance(sweep, dict) and isinstance(sweep.get("contact"), str):
        return str(sweep["contact"])
    return "Anode"


def _set_contact_bias(cfg: dict[str, object], contact_name: str, bias: float) -> None:
    contacts = cfg.get("contacts", [])
    if not isinstance(contacts, list):
        raise ValueError("config contacts must be a list")
    for contact in contacts:
        if isinstance(contact, dict) and contact.get("name") == contact_name:
            contact["bias"] = bias
            return
    raise ValueError(f"contact {contact_name!r} not found in config")


def derive_probe_config(base_config: Path, out_dir: Path, bias: float) -> Path:
    base_config = resolve_repo_path(base_config)
    out_dir = resolve_repo_path(out_dir)
    cfg = json.loads(base_config.read_text(encoding="utf-8"))
    absolutize_base_paths(cfg, base_config)
    cfg["simulation_type"] = "newton_jacobian_block_probe"
    cfg["state_file"] = str(state_path(out_dir, bias))
    cfg["output_csv"] = str(block_path(out_dir, bias))
    cfg["finite_difference_step"] = 1.0e-7
    _set_contact_bias(cfg, _sweep_contact_name(cfg), bias)
    cfg_path = out_dir / "configs" / f"jacobian_probe_{bias_token(bias)}.json"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return cfg_path


def run_runner(runner: Path, config: Path) -> None:
    subprocess.run([str(runner), "--config", str(config)], check=True, cwd=REPO)


def ensure_snapshot_states(runner: Path,
                           base_config: Path,
                           out_dir: Path,
                           biases: list[float],
                           force: bool) -> None:
    if not force and all(state_path(out_dir, bias).is_file() for bias in biases):
        return
    snapshot_config = derive_snapshot_config(base_config, out_dir)
    run_runner(runner, snapshot_config)
    missing = [bias for bias in biases if not state_path(out_dir, bias).is_file()]
    if missing:
        missing_text = ", ".join(f"{bias:g}" for bias in missing)
        raise RuntimeError(f"snapshot sweep did not write requested bias states: {missing_text}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", type=Path, default=DEFAULT_BASE_CONFIG)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--bias", action="append", type=float)
    parser.add_argument("--force-snapshot", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.base_config = resolve_repo_path(args.base_config)
    args.runner = resolve_repo_path(args.runner)
    args.out_dir = resolve_repo_path(args.out_dir)
    biases = args.bias if args.bias is not None else DEFAULT_BIASES

    args.out_dir.mkdir(parents=True, exist_ok=True)
    ensure_snapshot_states(
        args.runner,
        args.base_config,
        args.out_dir,
        biases,
        args.force_snapshot,
    )

    reports: list[BiasBlockReport] = []
    for bias in biases:
        probe_config = derive_probe_config(args.base_config, args.out_dir, bias)
        run_runner(args.runner, probe_config)
        reports.append(BiasBlockReport(bias, block_path(args.out_dir, bias)))
    write_combined_report(args.out_dir / "jacobian_blocks_real_state.csv", reports)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
