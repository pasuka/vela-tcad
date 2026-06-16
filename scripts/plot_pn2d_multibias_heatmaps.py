#!/usr/bin/env python3
"""Plot PN2D Sentaurus/Vela multibias state heatmaps."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import LogNorm, Normalize
from PIL import Image


QUANTITIES = {
    "potential": {
        "title": "Electrostatic potential",
        "unit": "V",
        "sentaurus_field": "ElectrostaticPotential",
        "vela_field": "Potential",
        "norm": "linear",
    },
    "electric_field": {
        "title": "Electric field magnitude",
        "unit": "V/cm",
        "sentaurus_field": "ElectricField",
        "vela_field": "__computed_electric_field__",
        "norm": "linear",
    },
    "electron_density": {
        "title": "Electron density",
        "unit": "cm^-3",
        "sentaurus_field": "eDensity",
        "vela_field": "Electrons",
        "vela_scale": 1.0e-6,
        "norm": "log",
    },
    "hole_density": {
        "title": "Hole density",
        "unit": "cm^-3",
        "sentaurus_field": "hDensity",
        "vela_field": "Holes",
        "vela_scale": 1.0e-6,
        "norm": "log",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--biases",
        default="0.5,1.0,1.5,2.0",
        help="Comma-separated biases to plot in volts.",
    )
    parser.add_argument(
        "--quantities",
        default=",".join(QUANTITIES),
        help="Comma-separated quantities to plot.",
    )
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def parse_quantities(raw: str) -> list[str]:
    quantities = [item.strip() for item in raw.split(",") if item.strip()]
    unknown = [item for item in quantities if item not in QUANTITIES]
    if unknown:
        raise ValueError(f"unknown quantities: {unknown}; expected one of {sorted(QUANTITIES)}")
    return quantities


def bias_token(bias: float) -> str:
    return f"{bias:g}".replace(".", "p")


def load_sentaurus_mesh(root: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes: dict[int, tuple[float, float]] = {}
    with (root / "nodes.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            nodes[int(row["id"])] = (float(row["x_um"]), float(row["y_um"]))
    triangles: list[list[int]] = []
    with (root / "elements.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            triangles.append([int(row["node0"]), int(row["node1"]), int(row["node2"])])
    node_ids = np.array(sorted(nodes), dtype=int)
    xy = np.array([nodes[int(node_id)] for node_id in node_ids], dtype=float)
    return xy[:, 0], xy[:, 1], np.array(triangles, dtype=int)


def load_sentaurus_scalar(root: Path, field: str) -> np.ndarray:
    values: dict[int, float] = {}
    path = root / "fields" / f"{field}_region0.csv"
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        component_cols = [name for name in reader.fieldnames or [] if name != "node_id"]
        for row in reader:
            comps = [float(row[name]) for name in component_cols if row.get(name, "") != ""]
            if len(comps) == 1:
                value = comps[0]
            else:
                value = math.sqrt(sum(component * component for component in comps))
            values[int(row["node_id"])] = value
    return np.array([values[node_id] for node_id in sorted(values)], dtype=float)


def parse_vtk(path: Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    points: list[tuple[float, float]] = []
    triangles: list[list[int]] = []
    scalars: dict[str, np.ndarray] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        if parts[0] == "POINTS":
            count = int(parts[1])
            for raw in lines[i + 1:i + 1 + count]:
                x, y, *_ = raw.split()
                points.append((float(x) * 1.0e6, float(y) * 1.0e6))
            i += 1 + count
            continue
        if parts[0] == "CELLS":
            count = int(parts[1])
            for raw in lines[i + 1:i + 1 + count]:
                cell = [int(item) for item in raw.split()]
                if cell and cell[0] == 3:
                    triangles.append(cell[1:4])
            i += 1 + count
            continue
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS",
                    "VECTORS",
                    "FIELD",
                    "CELL_DATA",
                    "POINT_DATA",
                }:
                    break
                values.extend(float(value) for value in next_parts)
                i += 1
            scalars[name] = np.array(values, dtype=float)
            continue
        i += 1
    return {
        "x": np.array([point[0] for point in points], dtype=float),
        "y": np.array([point[1] for point in points], dtype=float),
        "triangles": np.array(triangles, dtype=int),
        "scalars": scalars,
    }


def compute_node_electric_field_v_cm(
    x_um: np.ndarray,
    y_um: np.ndarray,
    triangles: np.ndarray,
    potential_v: np.ndarray,
) -> np.ndarray:
    x_m = x_um * 1.0e-6
    y_m = y_um * 1.0e-6
    accum = np.zeros_like(potential_v, dtype=float)
    counts = np.zeros_like(potential_v, dtype=float)
    for tri in triangles:
        ids = np.array(tri, dtype=int)
        a = np.column_stack([x_m[ids], y_m[ids], np.ones(3)])
        coeff = np.linalg.solve(a, potential_v[ids])
        magnitude_v_m = math.hypot(float(coeff[0]), float(coeff[1]))
        accum[ids] += magnitude_v_m * 1.0e-2
        counts[ids] += 1.0
    return np.divide(accum, counts, out=np.zeros_like(accum), where=counts > 0.0)


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in root.glob("*.vtk"):
        match = pattern.search(path.name)
        if match:
            result[round(float(match.group("bias")), 12)] = path
    return result


def finite_norm(values: list[np.ndarray], norm_kind: str) -> Normalize | LogNorm:
    merged = np.concatenate([np.asarray(value, dtype=float) for value in values])
    merged = merged[np.isfinite(merged)]
    if norm_kind == "log":
        merged = merged[merged > 0.0]
        return LogNorm(vmin=float(np.nanpercentile(merged, 1)), vmax=float(np.nanpercentile(merged, 99.5)))
    return Normalize(vmin=float(np.nanpercentile(merged, 1)), vmax=float(np.nanpercentile(merged, 99)))


def plot_quantity(
    quantity_key: str,
    spec: dict[str, Any],
    biases: list[float],
    sentaurus_exports: dict[float, Path],
    vela_vtks: dict[float, Path],
    out_dir: Path,
) -> dict[str, Any]:
    panels: list[dict[str, Any]] = []
    all_values: list[np.ndarray] = []
    for bias in biases:
        s_root = sentaurus_exports[bias]
        sx, sy, stri = load_sentaurus_mesh(s_root)
        svals = load_sentaurus_scalar(s_root, str(spec["sentaurus_field"]))
        vtk = parse_vtk(vela_vtks[round(bias, 12)])
        if spec["vela_field"] == "__computed_electric_field__":
            vvals = compute_node_electric_field_v_cm(
                vtk["x"], vtk["y"], vtk["triangles"], vtk["scalars"]["Potential"])
        else:
            vvals = vtk["scalars"][str(spec["vela_field"])] * float(spec.get("vela_scale", 1.0))
        panels.extend([
            {"bias": bias, "source": "Sentaurus", "x": sx, "y": sy, "tri": stri, "values": svals},
            {"bias": bias, "source": "Vela", "x": vtk["x"], "y": vtk["y"], "tri": vtk["triangles"], "values": vvals},
        ])
        all_values.extend([svals, vvals])

    norm = finite_norm(all_values, str(spec["norm"]))
    fig, axes = plt.subplots(len(biases), 2, figsize=(12, 3.1 * len(biases)), constrained_layout=True)
    cmap = "viridis" if spec["norm"] == "linear" else "magma"
    last_mesh = None
    for row, bias in enumerate(biases):
        for col, source in enumerate(["Sentaurus", "Vela"]):
            ax = axes[row, col] if len(biases) > 1 else axes[col]
            panel = next(item for item in panels if item["bias"] == bias and item["source"] == source)
            tri = mtri.Triangulation(panel["x"], panel["y"], panel["tri"])
            last_mesh = ax.tripcolor(tri, panel["values"], shading="gouraud", cmap=cmap, norm=norm)
            ax.triplot(tri, color="#1F2430", linewidth=0.16, alpha=0.34)
            ax.set_aspect("equal")
            ax.set_title(f"{source}  {bias:g} V", fontsize=10)
            ax.set_xlabel("x (um)")
            ax.set_ylabel("y (um)")
            ax.tick_params(labelsize=8)
    fig.suptitle(f"PN2D {spec['title']} comparison", fontsize=14, fontweight="semibold")
    cbar = fig.colorbar(last_mesh, ax=axes.ravel().tolist(), shrink=0.92, pad=0.015)
    cbar.set_label(f"{spec['title']} ({spec['unit']})")
    out_path = out_dir / f"pn2d_{quantity_key}_sentaurus_vela_heatmap.png"
    fig.savefig(out_path, dpi=220)
    plt.close(fig)
    with Image.open(out_path) as image:
        width, height = image.size
    return {
        "quantity": quantity_key,
        "path": str(out_path),
        "width_px": width,
        "height_px": height,
        "biases_V": biases,
    }


def main() -> int:
    args = parse_args()
    biases = parse_biases(args.biases)
    quantities = parse_quantities(args.quantities)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    sentaurus_exports = {
        bias: args.sentaurus_root / f"sentaurus_{bias_token(bias)}v"
        for bias in biases
    }
    missing = [str(path) for path in sentaurus_exports.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("missing Sentaurus exports: " + ", ".join(missing))
    vela_vtks = discover_vela_vtks(args.vela_vtk_root)
    missing_vtk = [bias for bias in biases if round(bias, 12) not in vela_vtks]
    if missing_vtk:
        raise FileNotFoundError(f"missing Vela VTK for biases: {missing_vtk}")

    manifest = {
        "schema": "pn2d.sentaurus_vela_heatmaps.v1",
        "outputs": [
            plot_quantity(key, spec, biases, sentaurus_exports, vela_vtks, args.out_dir)
            for key, spec in QUANTITIES.items()
            if key in quantities
        ],
    }
    manifest_path = args.out_dir / "pn2d_heatmap_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
