#!/usr/bin/env python3
"""Prepare a focused PN2D BV restart sweep from an accepted Vela VTK state."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag
import state_archive


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--restart-vtk", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias-points", required=True)
    parser.add_argument("--restart-name", default="restart_from_vtk.h5")
    parser.add_argument("--config-name", default="simulation.json")
    parser.add_argument("--vtk-prefix-name", default="focused_restart")
    # On Python < 3.13 argparse treats values like "-13.0,-13.1" as option
    # flags because they start with "-" but are not single negative numbers.
    # Broaden the negative-number matcher so comma-separated negative lists
    # passed to --bias-points are accepted as values across Python versions.
    parser._negative_number_matcher = re.compile(
        r"^-\d+$|^-\d*\.\d+$|^-?\d[\d.eE+-]*(?:,-?\d[\d.eE+-]*)+$"
    )
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


def write_restart_hdf5(path: Path, vtk: Path, base_config: Path) -> int:
    scalars = sgdiag.parse_vtk_scalars(vtk)
    require_scalars(scalars)
    count = len(scalars["Potential"])
    for name in ("ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"):
        if len(scalars[name]) != count:
            raise RuntimeError(f"restart VTK scalar length mismatch for {name}")
    config = json.loads(base_config.read_text(encoding="utf-8"))
    mesh_path = Path(config["mesh_file"])
    if not mesh_path.is_absolute(): mesh_path = base_config.parent / mesh_path
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    if len(mesh["nodes"]) != count:
        raise ValueError("VTK state and configured mesh have different node counts")
    tokens = vtk.read_text().split()
    first = tokens.index("POINTS")
    if int(tokens[first+1]) != count:
        raise ValueError("VTK point count mismatch")
    coordinates = list(map(float, tokens[first+3:first+3+3*count]))
    for index, node in enumerate(mesh["nodes"]):
        if node["id"] != index or coordinates[3*index:3*index+3] != [node["x"], node["y"], 0.]:
            raise ValueError("VTK point order/coordinates do not match the configured mesh")
    fields = {target: scalars[name] for target, name in {
        "psi": "Potential", "phin": "ElectronQuasiFermi", "phip": "HoleQuasiFermi",
        "electrons_m3": "Electrons", "holes_m3": "Holes"}.items()}
    state_archive.write(path, fields, dict(mode="dd",
        mesh_sha256=state_archive.mesh_identity(mesh,
            1e-6 if config.get("scaling", {}).get("mode") == "unit_scaling" else 1.),
        potential_origin_V=config.get("potential_origin_V", 0.),
        source_vtk_sha256=hashlib.sha256(vtk.read_bytes()).hexdigest()))
    return count


def make_config(base_config: Path, out_dir: Path, restart_state: Path, bias_points: list[float], vtk_prefix_name: str) -> dict[str, Any]:
    out_dir = out_dir.resolve()
    restart_state = restart_state.resolve()
    (out_dir / "vtk").mkdir(parents=True, exist_ok=True)
    config = json.loads(base_config.read_text(encoding="utf-8"))
    sweep = dict(config.get("sweep", {}))
    for key in ("mesh_file", "materials_file", "node_doping_file"):
        if config.get(key) and not Path(config[key]).is_absolute():
            config[key] = str((base_config.parent / config[key]).resolve())
    config["state_format"] = "hdf5"
    step = bias_points[1] - bias_points[0]
    sweep["start"] = bias_points[0]
    sweep["stop"] = bias_points[-1]
    sweep["step"] = step
    sweep["bias_points"] = bias_points
    sweep["csv_file"] = str(out_dir / "iv.csv")
    config["output_csv"] = str(out_dir / "iv.csv")
    sweep["initial_state_file"] = str(restart_state)
    sweep["write_state_file"] = str(out_dir / "latest_state.h5")
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
    restart_state = out_dir / args.restart_name
    node_count = write_restart_hdf5(restart_state, args.restart_vtk, args.base_config)
    config = make_config(args.base_config, out_dir, restart_state, bias_points, args.vtk_prefix_name)
    config_path = out_dir / args.config_name
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "config": str(config_path),
        "restart_hdf5": str(restart_state),
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
