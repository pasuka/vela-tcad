#!/usr/bin/env python3
"""Prepare a focused PN2D BV restart sweep from an accepted Vela VTK state."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


RESTART_FIELDS = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--restart-vtk", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias-points", required=True)
    parser.add_argument("--restart-name", default="restart_from_vtk.csv")
    parser.add_argument("--config-name", default="simulation.json")
    parser.add_argument("--vtk-prefix-name", default="focused_restart")
    return parser.parse_args()


def parse_bias_points(raw: str) -> list[float]:
    points = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if len(points) < 2:
        raise SystemExit("--bias-points requires at least two values")
    if any(not math.isfinite(value) for value in points):
        raise SystemExit("--bias-points contains non-finite value")
    return points


def require_scalars(scalars: dict[str, list[float]]) -> None:
    missing = [
        name for name in (
            "Potential",
            "ElectronQuasiFermi",
            "HoleQuasiFermi",
            "Electrons",
            "Holes",
        )
        if name not in scalars
    ]
    if missing:
        raise RuntimeError(f"restart VTK is missing required scalars: {', '.join(missing)}")


def write_restart_csv(path: Path, vtk: Path) -> int:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    require_scalars(scalars)
    count = len(scalars["Potential"])
    for name in ("ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"):
        if len(scalars[name]) != count:
            raise RuntimeError(f"restart VTK scalar length mismatch for {name}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESTART_FIELDS)
        writer.writeheader()
        for node_id in range(count):
            writer.writerow({
                "node_id": node_id,
                "psi": scalars["Potential"][node_id],
                "phin": scalars["ElectronQuasiFermi"][node_id],
                "phip": scalars["HoleQuasiFermi"][node_id],
                "electrons_m3": scalars["Electrons"][node_id],
                "holes_m3": scalars["Holes"][node_id],
            })
    return count


def make_config(base_config: Path, out_dir: Path, restart_csv: Path, bias_points: list[float], vtk_prefix_name: str) -> dict[str, Any]:
    out_dir = out_dir.resolve()
    restart_csv = restart_csv.resolve()
    (out_dir / "vtk").mkdir(parents=True, exist_ok=True)
    config = json.loads(base_config.read_text(encoding="utf-8"))
    sweep = dict(config.get("sweep", {}))
    step = bias_points[1] - bias_points[0]
    sweep["start"] = bias_points[0]
    sweep["stop"] = bias_points[-1]
    sweep["step"] = step
    sweep["bias_points"] = bias_points
    sweep["csv_file"] = str(out_dir / "iv.csv")
    config["output_csv"] = str(out_dir / "iv.csv")
    sweep["initial_state_file"] = str(restart_csv)
    sweep["write_state_file"] = str(out_dir / "latest_state.csv")
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str(out_dir / "vtk" / vtk_prefix_name)
    diagnostics = dict(sweep.get("diagnostics", {}))
    newton_history = dict(diagnostics.get("newton_history", {}))
    if newton_history.get("enabled", False):
        newton_history["csv_file"] = str(out_dir / "newton_history.csv")
        diagnostics["newton_history"] = newton_history
    sg_edges = dict(diagnostics.get("sg_avalanche_edges", {}))
    if sg_edges.get("enabled", False):
        sg_edges["csv_file"] = str(out_dir / "sg_avalanche_edges.csv")
        diagnostics["sg_avalanche_edges"] = sg_edges
    if diagnostics:
        sweep["diagnostics"] = diagnostics
    config["sweep"] = sweep
    return config


def main() -> int:
    args = parse_args()
    bias_points = parse_bias_points(args.bias_points)
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    restart_csv = out_dir / args.restart_name
    node_count = write_restart_csv(restart_csv, args.restart_vtk)
    config = make_config(args.base_config, out_dir, restart_csv, bias_points, args.vtk_prefix_name)
    config_path = out_dir / args.config_name
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "config": str(config_path),
        "restart_csv": str(restart_csv),
        "node_count": node_count,
        "bias_points": bias_points,
        "vtk_prefix": config["sweep"]["vtk_prefix"],
    }
    (out_dir / "focused_restart_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
