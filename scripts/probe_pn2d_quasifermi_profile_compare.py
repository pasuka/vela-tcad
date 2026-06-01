#!/usr/bin/env python3
"""Compare Sentaurus-exported potential/quasi-Fermi fields with a Vela VTK state.

This probe is node-based (same mesh node ids) and reports global + centerline-band
absolute-difference statistics. It is intended for pn2d diagnostics.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path


def read_sentaurus_scalar(fields_dir: Path, base_name: str, regions: list[int]) -> dict[int, float]:
    values: dict[int, float] = {}
    for region in regions:
        path = fields_dir / f"{base_name}_region{region}.csv"
        if not path.is_file():
            continue
        lines = path.read_text().splitlines()
        for raw in lines[1:]:
            if not raw.strip():
                continue
            node_id, component0, *_ = raw.split(",")
            values[int(node_id)] = float(component0)
    return values


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    npoints = 0
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "POINTS":
            npoints = int(parts[1])
            break

    fields: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < npoints:
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                if line.startswith("SCALARS") or line.startswith("CELL_DATA") or line.startswith("POINT_DATA"):
                    break
                values.extend(float(v) for v in line.split())
                i += 1
            fields[name] = values[:npoints]
            continue
        i += 1
    return fields


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    k = int(round((len(ordered) - 1) * p))
    return ordered[k]


def stats_for_pair(nodes: list[dict[str, float]],
                   sentaurus_values: dict[int, float],
                   vela_values: list[float],
                   centerline_halfwidth_um: float) -> dict[str, float | int | None]:
    common = sorted(set(sentaurus_values).intersection(range(len(vela_values))))
    if not common:
        return {"common_nodes": 0}

    diffs = [vela_values[i] - sentaurus_values[i] for i in common]
    abs_diffs = [abs(v) for v in diffs]

    ys = [nodes[i]["y"] for i in common]
    y_mid = 0.5 * (min(ys) + max(ys))
    band = [i for i in common if abs(nodes[i]["y"] - y_mid) <= centerline_halfwidth_um]
    band_abs = [abs(vela_values[i] - sentaurus_values[i]) for i in band]

    rmse = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    return {
        "common_nodes": len(common),
        "mean_abs_diff_V": sum(abs_diffs) / len(abs_diffs),
        "median_abs_diff_V": statistics.median(abs_diffs),
        "p95_abs_diff_V": quantile(abs_diffs, 0.95),
        "max_abs_diff_V": max(abs_diffs),
        "rmse_V": rmse,
        "centerline_band_y_um": y_mid,
        "centerline_band_halfwidth_um": centerline_halfwidth_um,
        "centerline_nodes": len(band),
        "centerline_mean_abs_diff_V": (sum(band_abs) / len(band_abs)) if band_abs else None,
        "centerline_max_abs_diff_V": max(band_abs) if band_abs else None,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vtk", type=Path, required=True)
    parser.add_argument("--sentaurus-fields-dir", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--centerline-halfwidth-um", type=float, default=0.01)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    mesh = json.loads(args.mesh.read_text())
    nodes = mesh["nodes"]
    vtk = parse_vtk_scalars(args.vtk)

    pairs = [
        ("ElectrostaticPotential", "Potential"),
        ("eQuasiFermiPotential", "ElectronQuasiFermi"),
        ("hQuasiFermiPotential", "HoleQuasiFermi"),
    ]

    report: dict[str, object] = {
        "note": (
            "Sentaurus IV TDR exported state is usually the final quasistationary point; "
            "if Vela VTK is from a different bias, treat results as trend-only."
        ),
        "mesh_nodes": len(nodes),
        "stats": {},
    }

    for sentaurus_name, vela_name in pairs:
        sentaurus_values = read_sentaurus_scalar(args.sentaurus_fields_dir, sentaurus_name, [0, 1])
        vela_values = vtk.get(vela_name, [])
        pair_stats = stats_for_pair(nodes, sentaurus_values, vela_values, args.centerline_halfwidth_um)
        pair_stats["vela_field"] = vela_name
        report["stats"][sentaurus_name] = pair_stats

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    print(args.output_json)
    print(json.dumps(report["stats"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
