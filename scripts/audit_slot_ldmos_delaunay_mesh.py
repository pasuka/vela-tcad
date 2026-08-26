#!/usr/bin/env python3
"""Audit edge-level Delaunay eligibility on an imported SLOT-LDMOS mesh.

The audit is intentionally read-only.  It preserves Vela's existing
``require_non_obtuse`` policy and reports whether negative *local* cotangent
contributions become non-negative after the two triangles adjacent to an
interior edge are assembled.  That distinction is the evidence needed before
considering a separate Delaunay-qualified atomic avalanche profile.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


EDGE_FIELDS = [
    "node0", "node1", "mid_x_um", "mid_y_um", "length_um",
    "adjacent_cell_count", "cell_ids", "region_ids", "materials",
    "boundary_edge", "region_interface_edge", "material_interface_edge",
    "silicon_only", "local_cotangents", "local_negative_count",
    "aggregated_cotangent", "signed_dual_length_um",
    "signed_dual_over_primal", "delaunay_qualified",
    "repaired_by_edge_aggregation", "unresolved_negative_internal",
    "distance_to_hotspot_um",
]

NODE_FIELDS = [
    "node_id", "x_um", "y_um", "raw_signed_dual_area_um2",
    "mixed_voronoi_area_um2", "barycentric_area_um2",
    "raw_signed_dual_negative", "distance_to_hotspot_um",
]

CELL_FIELDS = [
    "cell_id", "region_id", "material", "node_ids", "centroid_x_um",
    "centroid_y_um", "area_um2", "min_angle_deg", "max_angle_deg",
    "edge_aspect_ratio", "obtuse", "distance_to_hotspot_um",
    "extreme_reasons",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--hotspot-x-um", type=float, default=0.0060533)
    parser.add_argument("--hotspot-y-um", type=float, default=4.42555)
    parser.add_argument("--cotangent-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--extreme-cell-limit", type=int, default=100)
    parser.add_argument("--sentaurus-non-delaunay-cells", type=int, default=2)
    return parser.parse_args()


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def cotangent(a: tuple[float, float], b: tuple[float, float],
              opposite: tuple[float, float]) -> float:
    ux, uy = a[0] - opposite[0], a[1] - opposite[1]
    vx, vy = b[0] - opposite[0], b[1] - opposite[1]
    cross = ux * vy - uy * vx
    if abs(cross) <= 1.0e-30:
        return 0.0
    return (ux * vx + uy * vy) / abs(cross)


def angle_degrees(center: tuple[float, float], a: tuple[float, float],
                  b: tuple[float, float]) -> float:
    ux, uy = a[0] - center[0], a[1] - center[1]
    vx, vy = b[0] - center[0], b[1] - center[1]
    scale = math.hypot(ux, uy) * math.hypot(vx, vy)
    if scale <= 1.0e-30:
        return 0.0
    cosine = max(-1.0, min(1.0, (ux * vx + uy * vy) / scale))
    return math.degrees(math.acos(cosine))


def triangle_area(points: list[tuple[float, float]]) -> float:
    a, b, c = points
    return 0.5 * abs(
        (b[0] - a[0]) * (c[1] - a[1])
        - (c[0] - a[0]) * (b[1] - a[1])
    )


def raw_circumcentric_shares(points: list[tuple[float, float]]) -> list[float]:
    shares: list[float] = []
    for index in range(3):
        left = (index + 1) % 3
        right = (index + 2) % 3
        cot_left = cotangent(points[index], points[right], points[left])
        cot_right = cotangent(points[index], points[left], points[right])
        shares.append(0.125 * (
            distance(points[index], points[right]) ** 2 * cot_left
            + distance(points[index], points[left]) ** 2 * cot_right
        ))
    return shares


def mixed_voronoi_shares(points: list[tuple[float, float]],
                         angles: list[float], area: float) -> list[float]:
    max_index = max(range(3), key=lambda index: angles[index])
    if angles[max_index] > 90.0:
        shares = [0.25 * area] * 3
        shares[max_index] = 0.5 * area
        return shares
    return [max(value, 0.0) for value in raw_circumcentric_shares(points)]


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def minimum_distance(rows: list[dict[str, Any]]) -> float | None:
    if not rows:
        return None
    return min(float(row["distance_to_hotspot_um"]) for row in rows)


def audit(mesh_path: Path, hotspot: tuple[float, float], tolerance: float,
          extreme_limit: int, sentaurus_non_delaunay_cells: int) -> tuple[
              dict[str, Any], list[dict[str, Any]], list[dict[str, Any]],
              list[dict[str, Any]]]:
    data = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = {
        int(item["id"]): (float(item["x"]), float(item["y"]))
        for item in data["nodes"]
    }
    regions = {int(item["id"]): item for item in data.get("regions", [])}
    triangles = data["triangles"]

    edge_records: dict[tuple[int, int], dict[str, Any]] = {}
    raw_node_area = defaultdict(float)
    mixed_node_area = defaultdict(float)
    barycentric_node_area = defaultdict(float)
    cell_rows_all: list[dict[str, Any]] = []

    for cell in triangles:
        cell_id = int(cell["id"])
        ids = [int(value) for value in cell["node_ids"][:3]]
        points = [nodes[node_id] for node_id in ids]
        area = triangle_area(points)
        angles = [
            angle_degrees(points[0], points[1], points[2]),
            angle_degrees(points[1], points[2], points[0]),
            angle_degrees(points[2], points[0], points[1]),
        ]
        lengths = [
            distance(points[0], points[1]),
            distance(points[1], points[2]),
            distance(points[2], points[0]),
        ]
        region_id = int(cell["region_id"])
        material = str(regions.get(region_id, {}).get("material", "unknown"))
        centroid = (sum(point[0] for point in points) / 3.0,
                    sum(point[1] for point in points) / 3.0)
        cell_rows_all.append({
            "cell_id": cell_id,
            "region_id": region_id,
            "material": material,
            "node_ids": ";".join(str(value) for value in ids),
            "centroid_x_um": centroid[0],
            "centroid_y_um": centroid[1],
            "area_um2": area,
            "min_angle_deg": min(angles),
            "max_angle_deg": max(angles),
            "edge_aspect_ratio": max(lengths) / min(lengths),
            "obtuse": max(angles) > 90.0,
            "distance_to_hotspot_um": distance(centroid, hotspot),
        })

        raw_shares = raw_circumcentric_shares(points)
        mixed_shares = mixed_voronoi_shares(points, angles, area)
        for index, node_id in enumerate(ids):
            raw_node_area[node_id] += raw_shares[index]
            mixed_node_area[node_id] += mixed_shares[index]
            barycentric_node_area[node_id] += area / 3.0

        for index in range(3):
            a, b = ids[index], ids[(index + 1) % 3]
            opposite = ids[(index + 2) % 3]
            key = edge_key(a, b)
            record = edge_records.setdefault(key, {
                "node0": key[0], "node1": key[1], "cells": [],
                "regions": [], "materials": [], "cotangents": [],
            })
            record["cells"].append(cell_id)
            record["regions"].append(region_id)
            record["materials"].append(material)
            record["cotangents"].append(
                cotangent(nodes[a], nodes[b], nodes[opposite]))

    edge_rows: list[dict[str, Any]] = []
    for key, record in sorted(edge_records.items()):
        point0, point1 = nodes[key[0]], nodes[key[1]]
        midpoint = ((point0[0] + point1[0]) / 2.0,
                    (point0[1] + point1[1]) / 2.0)
        length = distance(point0, point1)
        cotangents = record["cotangents"]
        aggregate = sum(cotangents)
        internal = len(record["cells"]) == 2
        regions_on_edge = sorted(set(record["regions"]))
        materials_on_edge = sorted(set(record["materials"]))
        # Match GeometryBuildReport semantics exactly: every negative local
        # contribution is counted, while the assembled edge qualification
        # retains a tolerance around a zero dual face.
        local_negative = sum(value < 0.0 for value in cotangents)
        qualified = internal and aggregate >= -tolerance
        edge_rows.append({
            "node0": key[0],
            "node1": key[1],
            "mid_x_um": midpoint[0],
            "mid_y_um": midpoint[1],
            "length_um": length,
            "adjacent_cell_count": len(record["cells"]),
            "cell_ids": ";".join(str(value) for value in record["cells"]),
            "region_ids": ";".join(str(value) for value in regions_on_edge),
            "materials": ";".join(materials_on_edge),
            "boundary_edge": len(record["cells"]) == 1,
            "region_interface_edge": len(regions_on_edge) > 1,
            "material_interface_edge": len(materials_on_edge) > 1,
            "silicon_only": materials_on_edge == ["Si"],
            "local_cotangents": ";".join(f"{value:.17g}" for value in cotangents),
            "local_negative_count": local_negative,
            "aggregated_cotangent": aggregate,
            "signed_dual_length_um": 0.5 * aggregate * length,
            "signed_dual_over_primal": 0.5 * aggregate,
            "delaunay_qualified": qualified,
            "repaired_by_edge_aggregation": internal and local_negative > 0 and qualified,
            "unresolved_negative_internal": internal and aggregate < -tolerance,
            "distance_to_hotspot_um": distance(midpoint, hotspot),
        })

    node_rows = [{
        "node_id": node_id,
        "x_um": point[0],
        "y_um": point[1],
        "raw_signed_dual_area_um2": raw_node_area[node_id],
        "mixed_voronoi_area_um2": mixed_node_area[node_id],
        "barycentric_area_um2": barycentric_node_area[node_id],
        "raw_signed_dual_negative": raw_node_area[node_id] < -tolerance,
        "distance_to_hotspot_um": distance(point, hotspot),
    } for node_id, point in sorted(nodes.items())]

    # Retain the union of the worst cells by each conditioning metric.
    extreme_ids: dict[int, set[str]] = defaultdict(set)
    metric_orders = {
        "largest_angle": sorted(cell_rows_all, key=lambda row: float(row["max_angle_deg"]), reverse=True),
        "smallest_angle": sorted(cell_rows_all, key=lambda row: float(row["min_angle_deg"])),
        "largest_aspect": sorted(cell_rows_all, key=lambda row: float(row["edge_aspect_ratio"]), reverse=True),
        "closest_hotspot": sorted(cell_rows_all, key=lambda row: float(row["distance_to_hotspot_um"])),
    }
    for reason, ordered in metric_orders.items():
        for row in ordered[:extreme_limit]:
            extreme_ids[int(row["cell_id"])].add(reason)
    cell_rows = []
    for row in cell_rows_all:
        cell_id = int(row["cell_id"])
        if cell_id in extreme_ids:
            row = dict(row)
            row["extreme_reasons"] = ";".join(sorted(extreme_ids[cell_id]))
            cell_rows.append(row)
    cell_rows.sort(key=lambda row: (float(row["distance_to_hotspot_um"]), int(row["cell_id"])))

    interior = [row for row in edge_rows if int(row["adjacent_cell_count"]) == 2]
    unresolved = [row for row in interior if row["unresolved_negative_internal"]]
    repaired = [row for row in interior if row["repaired_by_edge_aggregation"]]
    silicon_interior = [row for row in interior if row["silicon_only"]]
    silicon_unresolved = [row for row in silicon_interior if row["unresolved_negative_internal"]]
    interface = [row for row in edge_rows if row["material_interface_edge"]]
    interface_unresolved = [row for row in interface if row["unresolved_negative_internal"]]
    boundary = [row for row in edge_rows if row["boundary_edge"]]
    negative_nodes = [row for row in node_rows if row["raw_signed_dual_negative"]]
    obtuse_cells = [row for row in cell_rows_all if row["obtuse"]]
    summary = {
        "mesh": str(mesh_path.resolve()),
        "hotspot_um": {"x": hotspot[0], "y": hotspot[1]},
        "cotangent_tolerance": tolerance,
        "counts": {
            "nodes": len(nodes),
            "cells": len(cell_rows_all),
            "obtuse_cells": len(obtuse_cells),
            "edges": len(edge_rows),
            "interior_edges": len(interior),
            "boundary_edges": len(boundary),
            "region_interface_edges": sum(bool(row["region_interface_edge"]) for row in edge_rows),
            "material_interface_edges": len(interface),
            "edges_with_negative_local_cotangent": sum(int(row["local_negative_count"]) > 0 for row in edge_rows),
            "negative_local_cotangent_contributions": sum(int(row["local_negative_count"]) for row in edge_rows),
            "repaired_by_edge_aggregation": len(repaired),
            "unresolved_negative_interior_edges": len(unresolved),
            "unresolved_negative_silicon_interior_edges": len(silicon_unresolved),
            "unresolved_negative_material_interface_edges": len(interface_unresolved),
            "negative_raw_signed_dual_nodes": len(negative_nodes),
        },
        "fractions": {
            "interior_edges_delaunay_qualified": (
                (len(interior) - len(unresolved)) / len(interior) if interior else None),
            "negative_local_edges_repaired": (
                len(repaired) /
                sum(int(row["local_negative_count"]) > 0 for row in interior)
                if any(int(row["local_negative_count"]) > 0 for row in interior)
                else None),
        },
        "hotspot_distance_um": {
            "closest_obtuse_cell_centroid": minimum_distance(obtuse_cells),
            "closest_unresolved_negative_interior_edge": minimum_distance(unresolved),
            "closest_unresolved_negative_silicon_edge": minimum_distance(silicon_unresolved),
            "closest_material_interface_edge": minimum_distance(interface),
        },
        "sentaurus_crosscheck": {
            "reported_non_delaunay_cells": sentaurus_non_delaunay_cells,
            "vela_unresolved_negative_interior_edges": len(unresolved),
            "same_numeric_count": len(unresolved) == sentaurus_non_delaunay_cells,
            "note": (
                "Sentaurus reports cells whereas this audit reports assembled interior edges; "
                "matching counts are supportive but are not an identity proof."
            ),
        },
        "decision": {
            "new_profile_eligible_for_prototyping": len(silicon_unresolved) == 0,
            "reason": (
                "All Silicon interior edge cotangent sums are non-negative."
                if len(silicon_unresolved) == 0 else
                "Negative assembled Silicon dual faces require local fallback or remeshing."
            ),
            "runtime_behavior_changed": False,
        },
    }
    return summary, edge_rows, node_rows, cell_rows


def main() -> int:
    args = parse_args()
    if args.extreme_cell_limit <= 0:
        raise ValueError("--extreme-cell-limit must be positive")
    summary, edge_rows, node_rows, cell_rows = audit(
        args.mesh, (args.hotspot_x_um, args.hotspot_y_um),
        args.cotangent_tolerance, args.extreme_cell_limit,
        args.sentaurus_non_delaunay_cells)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "delaunay_edges.csv", edge_rows, EDGE_FIELDS)
    write_csv(args.out_dir / "delaunay_node_volumes.csv", node_rows, NODE_FIELDS)
    write_csv(args.out_dir / "delaunay_extreme_cells.csv", cell_rows, CELL_FIELDS)
    (args.out_dir / "delaunay_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
