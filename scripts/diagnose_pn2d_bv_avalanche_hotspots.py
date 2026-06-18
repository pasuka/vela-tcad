#!/usr/bin/env python3
"""Rank PN2D BV avalanche hotspots and nearby mesh geometry."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


HOTSPOT_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "avalanche_generation_m3_s",
    "electron_density_m3",
    "hole_density_m3",
    "electric_field_V_m",
    "electron_high_field_drive_V_m",
    "hole_high_field_drive_V_m",
    "node_control_volume_m2",
    "adjacent_element_count",
    "min_adjacent_angle_degrees",
    "max_adjacent_angle_degrees",
    "obtuse_adjacent_element_count",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtk", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--top", type=int, default=50)
    return parser.parse_args()


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    scalars: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
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
            scalars[name] = values
            continue
        i += 1
    return scalars


def triangle_area_m2(points_m: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = points_m
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) * 0.5


def angle_degrees(a: tuple[float, float],
                  b: tuple[float, float],
                  c: tuple[float, float]) -> float:
    bax = a[0] - b[0]
    bay = a[1] - b[1]
    bcx = c[0] - b[0]
    bcy = c[1] - b[1]
    norm_ba = math.hypot(bax, bay)
    norm_bc = math.hypot(bcx, bcy)
    if norm_ba == 0.0 or norm_bc == 0.0:
        return 0.0
    cos_theta = max(-1.0, min(1.0, (bax * bcx + bay * bcy) / (norm_ba * norm_bc)))
    return math.degrees(math.acos(cos_theta))


def triangle_angles(points_m: list[tuple[float, float]]) -> list[float]:
    return [
        angle_degrees(points_m[1], points_m[0], points_m[2]),
        angle_degrees(points_m[0], points_m[1], points_m[2]),
        angle_degrees(points_m[0], points_m[2], points_m[1]),
    ]


def read_mesh_geometry(path: Path) -> tuple[dict[int, dict[str, float]], dict[int, dict[str, Any]]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {"x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data["nodes"]
    }
    triangles = {
        int(triangle["id"]): {
            "node_ids": [int(node_id) for node_id in triangle["node_ids"]],
        }
        for triangle in data["triangles"]
    }
    return nodes, triangles


def geometry_by_node(nodes: dict[int, dict[str, float]],
                     triangles: dict[int, dict[str, Any]]) -> dict[int, dict[str, float]]:
    geometry = {
        node_id: {
            "node_control_volume_m2": 0.0,
            "adjacent_element_count": 0.0,
            "min_adjacent_angle_degrees": math.inf,
            "max_adjacent_angle_degrees": 0.0,
            "obtuse_adjacent_element_count": 0.0,
        }
        for node_id in nodes
    }
    for triangle in triangles.values():
        node_ids = triangle["node_ids"]
        points_m = [
            (nodes[node_id]["x_um"] * 1.0e-6, nodes[node_id]["y_um"] * 1.0e-6)
            for node_id in node_ids
        ]
        area = triangle_area_m2(points_m)
        angles = triangle_angles(points_m)
        obtuse = any(angle > 90.0 for angle in angles)
        for node_id in node_ids:
            item = geometry[node_id]
            item["node_control_volume_m2"] += area / 3.0
            item["adjacent_element_count"] += 1.0
            item["min_adjacent_angle_degrees"] = min(
                item["min_adjacent_angle_degrees"], min(angles))
            item["max_adjacent_angle_degrees"] = max(
                item["max_adjacent_angle_degrees"], max(angles))
            if obtuse:
                item["obtuse_adjacent_element_count"] += 1.0
    for item in geometry.values():
        if math.isinf(item["min_adjacent_angle_degrees"]):
            item["min_adjacent_angle_degrees"] = 0.0
    return geometry


def scalar_value(scalars: dict[str, list[float]], name: str, node_id: int) -> float | None:
    values = scalars.get(name)
    if values is None or node_id >= len(values):
        return None
    return values[node_id]


def vtk_field_cm_to_m(scalars: dict[str, list[float]], name: str, node_id: int) -> float | None:
    value = scalar_value(scalars, name, node_id)
    if value is None:
        return None
    return value * 100.0


def make_hotspot_rows(nodes: dict[int, dict[str, float]],
                      scalars: dict[str, list[float]],
                      geometry: dict[int, dict[str, float]],
                      top: int) -> list[dict[str, Any]]:
    avalanche = scalars.get("AvalancheGeneration")
    if avalanche is None:
        raise RuntimeError("VTK is missing required scalar AvalancheGeneration")
    ranked_ids = sorted(nodes, key=lambda node_id: abs(avalanche[node_id]), reverse=True)
    rows: list[dict[str, Any]] = []
    for node_id in ranked_ids[:max(0, top)]:
        node = nodes[node_id]
        geom = geometry[node_id]
        rows.append({
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "avalanche_generation_m3_s": avalanche[node_id],
            "electron_density_m3": scalar_value(scalars, "Electrons", node_id),
            "hole_density_m3": scalar_value(scalars, "Holes", node_id),
            "electric_field_V_m": vtk_field_cm_to_m(scalars, "ElectricField", node_id),
            "electron_high_field_drive_V_m": vtk_field_cm_to_m(
                scalars, "ElectronHighFieldDrive", node_id),
            "hole_high_field_drive_V_m": vtk_field_cm_to_m(
                scalars, "HoleHighFieldDrive", node_id),
            "node_control_volume_m2": geom["node_control_volume_m2"],
            "adjacent_element_count": int(geom["adjacent_element_count"]),
            "min_adjacent_angle_degrees": geom["min_adjacent_angle_degrees"],
            "max_adjacent_angle_degrees": geom["max_adjacent_angle_degrees"],
            "obtuse_adjacent_element_count": int(geom["obtuse_adjacent_element_count"]),
        })
    return rows


def concentration(values: list[float], n: int) -> float:
    total = sum(abs(value) for value in values)
    if total == 0.0:
        return 0.0
    return sum(abs(value) for value in values[:min(n, len(values))]) / total


def make_summary(rows: list[dict[str, Any]], node_count: int) -> dict[str, Any]:
    ordered_values = [float(row["avalanche_generation_m3_s"]) for row in rows]
    top = rows[0] if rows else {}
    return {
        "node_count": node_count,
        "reported_hotspot_count": len(rows),
        "top_node": {
            "node_id": int(top["node_id"]) if top else None,
            "avalanche_generation_m3_s": top.get("avalanche_generation_m3_s"),
            "x_um": top.get("x_um"),
            "y_um": top.get("y_um"),
        },
        "concentration": {
            "top1_fraction": concentration(ordered_values, 1),
            "top5_fraction": concentration(ordered_values, 5),
            "top20_fraction": concentration(ordered_values, 20),
        },
        "geometry": {
            "top1_node_control_volume_m2": top.get("node_control_volume_m2"),
            "top1_adjacent_element_count": top.get("adjacent_element_count"),
            "top1_min_adjacent_angle_degrees": top.get("min_adjacent_angle_degrees"),
            "top1_max_adjacent_angle_degrees": top.get("max_adjacent_angle_degrees"),
            "top1_obtuse_adjacent_element_count": top.get(
                "obtuse_adjacent_element_count"),
        },
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HOTSPOT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    if args.top <= 0:
        raise SystemExit("--top must be positive")
    nodes, triangles = read_mesh_geometry(args.mesh)
    scalars = parse_vtk_scalars(args.vtk)
    geometry = geometry_by_node(nodes, triangles)
    rows = make_hotspot_rows(nodes, scalars, geometry, args.top)
    summary = make_summary(rows, len(nodes))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "avalanche_hotspots.csv", rows)
    (args.out_dir / "avalanche_hotspots_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
