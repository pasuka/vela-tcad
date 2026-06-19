#!/usr/bin/env python3
"""Estimate local ionization-integral path proxies from a Vela VTK state."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


FIELDS = [
    "bias_V",
    "path_id",
    "carrier",
    "start_x_um",
    "start_y_um",
    "end_x_um",
    "end_y_um",
    "max_field_V_per_m",
    "integral_alpha_dx",
    "dominant_material",
    "dominant_edge_or_cell_id",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtk", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--top", type=int, default=1)
    parser.add_argument(
        "--field-scale",
        type=float,
        default=100.0,
        help="Multiplier applied to VTK high-field scalars before alpha(E); Vela VTK stores these in V/cm scale.",
    )
    return parser.parse_args()


def read_mesh(path: Path) -> tuple[dict[int, dict[str, float]], list[dict[str, Any]], dict[int, str]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {"x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data["nodes"]
    }
    regions = {
        int(region["id"]): str(region.get("material") or region.get("name") or "")
        for region in data.get("regions", [])
    }
    triangles = []
    for triangle in data["triangles"]:
        region_id = int(triangle.get("region_id", 0))
        triangles.append({
            "id": int(triangle["id"]),
            "region_id": region_id,
            "node_ids": [int(node_id) for node_id in triangle["node_ids"]],
        })
    return nodes, triangles, regions


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a <= b else (b, a)


def edge_length_m(nodes: dict[int, dict[str, float]], a: int, b: int) -> float:
    ax = nodes[a]["x_um"] * 1.0e-6
    ay = nodes[a]["y_um"] * 1.0e-6
    bx = nodes[b]["x_um"] * 1.0e-6
    by = nodes[b]["y_um"] * 1.0e-6
    return math.hypot(bx - ax, by - ay)


def build_edges(nodes: dict[int, dict[str, float]], triangles: list[dict[str, Any]]) -> list[dict[str, Any]]:
    edge_map: dict[tuple[int, int], dict[str, Any]] = {}
    for triangle in triangles:
        ids = triangle["node_ids"]
        for i in range(3):
            a, b = edge_key(ids[i], ids[(i + 1) % 3])
            key = (a, b)
            if key not in edge_map:
                edge_map[key] = {
                    "edge_id": len(edge_map),
                    "node0": a,
                    "node1": b,
                    "cell_ids": [],
                    "region_ids": set(),
                    "length_m": edge_length_m(nodes, a, b),
                }
            edge_map[key]["cell_ids"].append(triangle["id"])
            edge_map[key]["region_ids"].add(triangle["region_id"])
    return list(edge_map.values())


def scalar_field(scalars: dict[str, list[float]], preferred: str, fallback: str) -> list[float]:
    if preferred in scalars:
        return scalars[preferred]
    if fallback in scalars:
        return scalars[fallback]
    raise KeyError(f"missing required scalar '{preferred}' or fallback '{fallback}'")


def best_edges(carrier: str,
               edges: list[dict[str, Any]],
               nodes: dict[int, dict[str, float]],
               regions: dict[int, str],
               field: list[float],
               top: int,
               field_scale: float) -> list[dict[str, Any]]:
    scored = []
    for edge in edges:
        n0 = int(edge["node0"])
        n1 = int(edge["node1"])
        avg_field = 0.5 * (abs(field[n0]) + abs(field[n1])) * field_scale
        if carrier == "electron":
            alpha = sgdiag.van_overstraeten_alpha(avg_field, 7.03e7, 7.03e7, 1.231e8, 1.231e8)
        else:
            alpha = sgdiag.van_overstraeten_alpha(avg_field, 1.582e8, 6.71e7, 2.036e8, 1.693e8)
        integral = alpha * float(edge["length_m"])
        region_id = sorted(edge["region_ids"])[0] if edge["region_ids"] else 0
        scored.append({
            "carrier": carrier,
            "edge": edge,
            "field": avg_field,
            "integral": integral,
            "material": regions.get(region_id, ""),
        })
    scored.sort(key=lambda item: item["integral"], reverse=True)
    return scored[:max(1, top)]


def row_for(path_id: int, bias: float, item: dict[str, Any], nodes: dict[int, dict[str, float]]) -> dict[str, Any]:
    edge = item["edge"]
    n0 = int(edge["node0"])
    n1 = int(edge["node1"])
    return {
        "bias_V": bias,
        "path_id": path_id,
        "carrier": item["carrier"],
        "start_x_um": nodes[n0]["x_um"],
        "start_y_um": nodes[n0]["y_um"],
        "end_x_um": nodes[n1]["x_um"],
        "end_y_um": nodes[n1]["y_um"],
        "max_field_V_per_m": item["field"],
        "integral_alpha_dx": item["integral"],
        "dominant_material": item["material"],
        "dominant_edge_or_cell_id": edge["edge_id"],
    }


def write_outputs(args: argparse.Namespace) -> None:
    nodes, triangles, regions = read_mesh(args.mesh)
    scalars = sgdiag.parse_vtk_scalars(args.vtk)
    edges = build_edges(nodes, triangles)
    electron_field = scalar_field(scalars, "ElectronHighFieldDrive", "ElectricField")
    hole_field = scalar_field(scalars, "HoleHighFieldDrive", "ElectricField")

    rows = []
    path_id = 0
    for carrier, field in [("electron", electron_field), ("hole", hole_field)]:
        for item in best_edges(carrier, edges, nodes, regions, field, args.top, args.field_scale):
            rows.append(row_for(path_id, args.bias, item, nodes))
            path_id += 1

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "ionization_integral_paths.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    by_carrier: dict[str, float] = defaultdict(float)
    for row in rows:
        by_carrier[str(row["carrier"])] = max(
            by_carrier[str(row["carrier"])],
            float(row["integral_alpha_dx"]),
        )
    summary = {
        "schema": "pn2d.ionization_integral_diagnostic.v1",
        "bias_V": args.bias,
        "vtk": str(args.vtk),
        "mesh": str(args.mesh),
        "path_count": len(rows),
        "path_model": "single dominant mesh-edge proxy",
        "field_scale": args.field_scale,
        "max_integral_alpha_dx": max((float(row["integral_alpha_dx"]) for row in rows), default=0.0),
        "max_integral_by_carrier": dict(by_carrier),
    }
    (args.out_dir / "ionization_integral_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    write_outputs(parse_args())


if __name__ == "__main__":
    main()
