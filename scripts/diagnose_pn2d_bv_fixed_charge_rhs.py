#!/usr/bin/env python3
"""Diagnose nodal versus control-volume-integrated fixed-charge RHS."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


Q = 1.602176634e-19
EPS = 1.0e-30

NODE_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "band",
    "contact_names",
    "control_volume_m2",
    "mixed_voronoi_volume_m2",
    "mixed_to_barycentric_volume_ratio",
    "nodal_net_doping_cm3",
    "integrated_net_doping_cm3",
    "mixed_volume_net_doping_cm3",
    "delta_net_doping_cm3",
    "mixed_volume_delta_net_doping_cm3",
    "integrated_to_nodal_ratio",
    "mixed_volume_to_nodal_ratio",
    "nodal_fixed_charge_C_per_m",
    "integrated_fixed_charge_C_per_m",
    "mixed_volume_fixed_charge_C_per_m",
    "delta_fixed_charge_C_per_m",
    "mixed_volume_delta_fixed_charge_C_per_m",
    "abs_delta_fixed_charge_C_per_m",
    "abs_mixed_volume_delta_fixed_charge_C_per_m",
]

BAND_FIELDS = [
    "band",
    "x_min_um",
    "x_max_um",
    "node_count",
    "median_nodal_net_doping_cm3",
    "median_integrated_net_doping_cm3",
    "median_mixed_volume_net_doping_cm3",
    "median_delta_net_doping_cm3",
    "median_mixed_volume_delta_net_doping_cm3",
    "median_integrated_to_nodal_ratio",
    "median_mixed_volume_to_nodal_ratio",
    "median_delta_fixed_charge_C_per_m",
    "median_mixed_volume_delta_fixed_charge_C_per_m",
    "median_abs_delta_fixed_charge_C_per_m",
    "median_abs_mixed_volume_delta_fixed_charge_C_per_m",
    "max_abs_delta_fixed_charge_C_per_m",
    "max_abs_delta_fixed_charge_node",
    "max_abs_mixed_volume_delta_fixed_charge_C_per_m",
    "max_abs_mixed_volume_delta_fixed_charge_node",
]

PROFILE_FIELDS = [
    "segment",
    "x_min_um",
    "x_max_um",
    "net_doping_cm3",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--bands",
        default=(
            "left_p:0.0:0.25,pre_junction_p:0.7:0.9,junction:0.95:1.05,"
            "post_junction_n:1.1:1.3,right_n:1.5:1.8"
        ),
        help="Comma-separated name:xmin:xmax bands in um.",
    )
    parser.add_argument(
        "--profile",
        default="auto",
        help=(
            "Doping profile source. Use 'auto' to infer x-piecewise constants "
            "from nodal net doping, or a comma-separated list of xmin:xmax:net_cm3."
        ),
    )
    return parser.parse_args()


def parse_bands(raw: str) -> list[tuple[str, float, float]]:
    result = []
    for item in raw.split(","):
        if not item.strip():
            continue
        name, x_min, x_max = item.split(":", 2)
        result.append((name.strip(), float(x_min), float(x_max)))
    return result


def read_mesh(path: Path) -> tuple[dict[int, dict[str, float]], list[dict[str, Any]], dict[int, set[str]]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {
            "x_um": float(node["x"]),
            "y_um": float(node["y"]),
        }
        for node in data["nodes"]
    }
    triangles = [
        {
            "id": int(triangle["id"]),
            "node_ids": [int(node_id) for node_id in triangle["node_ids"]],
        }
        for triangle in data["triangles"]
    ]
    contact_by_node: dict[int, set[str]] = defaultdict(set)
    for contact in data.get("contacts", []):
        name = str(contact.get("name", "contact"))
        for node_id in contact.get("node_ids", []):
            contact_by_node[int(node_id)].add(name)
    return nodes, triangles, contact_by_node


def read_doping(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    seen = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            if node_id < 0 or node_id >= node_count:
                raise RuntimeError(f"{path} references node {node_id}, expected 0..{node_count - 1}")
            donors[node_id] = float(row["donors_cm3"])
            acceptors[node_id] = float(row["acceptors_cm3"])
            seen.add(node_id)
    missing = sorted(set(range(node_count)) - seen)
    if missing:
        raise RuntimeError(f"{path} missing node ids: {missing[:8]}")
    return donors, acceptors


def triangle_area_m2(points_um: list[tuple[float, float]]) -> float:
    (ax, ay), (bx, by), (cx, cy) = points_um[:3]
    ax *= 1.0e-6
    ay *= 1.0e-6
    bx *= 1.0e-6
    by *= 1.0e-6
    cx *= 1.0e-6
    cy *= 1.0e-6
    return 0.5 * abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))


def angle_radians(a: tuple[float, float],
                  b: tuple[float, float],
                  c: tuple[float, float]) -> float:
    ax, ay = a[0] * 1.0e-6, a[1] * 1.0e-6
    bx, by = b[0] * 1.0e-6, b[1] * 1.0e-6
    cx, cy = c[0] * 1.0e-6, c[1] * 1.0e-6
    ux, uy = bx - ax, by - ay
    vx, vy = cx - ax, cy - ay
    ul = math.hypot(ux, uy)
    vl = math.hypot(vx, vy)
    if ul <= 0.0 or vl <= 0.0:
        return 0.0
    return math.acos(max(-1.0, min(1.0, (ux * vx + uy * vy) / (ul * vl))))


def cotangent_opposite(a: tuple[float, float],
                       b: tuple[float, float],
                       opp: tuple[float, float]) -> float:
    ax, ay = a[0] * 1.0e-6, a[1] * 1.0e-6
    bx, by = b[0] * 1.0e-6, b[1] * 1.0e-6
    ox, oy = opp[0] * 1.0e-6, opp[1] * 1.0e-6
    ux, uy = ax - ox, ay - oy
    vx, vy = bx - ox, by - oy
    cross = ux * vy - uy * vx
    if abs(cross) <= EPS:
        return 0.0
    return (ux * vx + uy * vy) / abs(cross)


def distance2_m2(a: tuple[float, float], b: tuple[float, float]) -> float:
    dx = (a[0] - b[0]) * 1.0e-6
    dy = (a[1] - b[1]) * 1.0e-6
    return dx * dx + dy * dy


def polygon_area_m2(poly: list[tuple[float, float]]) -> float:
    if len(poly) < 3:
        return 0.0
    twice = 0.0
    for idx, (x0, y0) in enumerate(poly):
        x1, y1 = poly[(idx + 1) % len(poly)]
        twice += (x0 * 1.0e-6) * (y1 * 1.0e-6) - (x1 * 1.0e-6) * (y0 * 1.0e-6)
    return 0.5 * abs(twice)


def clip_polygon_x_min(poly: list[tuple[float, float]], x_min: float) -> list[tuple[float, float]]:
    if not poly:
        return []
    result = []
    prev = poly[-1]
    prev_inside = prev[0] >= x_min - EPS
    for curr in poly:
        curr_inside = curr[0] >= x_min - EPS
        if curr_inside != prev_inside:
            denom = curr[0] - prev[0]
            if abs(denom) > EPS:
                t = (x_min - prev[0]) / denom
                result.append((x_min, prev[1] + t * (curr[1] - prev[1])))
        if curr_inside:
            result.append(curr)
        prev = curr
        prev_inside = curr_inside
    return result


def clip_polygon_x_max(poly: list[tuple[float, float]], x_max: float) -> list[tuple[float, float]]:
    if not poly:
        return []
    result = []
    prev = poly[-1]
    prev_inside = prev[0] <= x_max + EPS
    for curr in poly:
        curr_inside = curr[0] <= x_max + EPS
        if curr_inside != prev_inside:
            denom = curr[0] - prev[0]
            if abs(denom) > EPS:
                t = (x_max - prev[0]) / denom
                result.append((x_max, prev[1] + t * (curr[1] - prev[1])))
        if curr_inside:
            result.append(curr)
        prev = curr
        prev_inside = curr_inside
    return result


def clip_polygon_to_x_interval(
    poly: list[tuple[float, float]], x_min: float, x_max: float
) -> list[tuple[float, float]]:
    result = poly
    if math.isfinite(x_min):
        result = clip_polygon_x_min(result, x_min)
    if math.isfinite(x_max):
        result = clip_polygon_x_max(result, x_max)
    return result


def vertex_control_polygon(
    triangle_points: list[tuple[float, float]], vertex_index: int
) -> list[tuple[float, float]]:
    v = triangle_points[vertex_index]
    nxt = triangle_points[(vertex_index + 1) % 3]
    prv = triangle_points[(vertex_index + 2) % 3]
    centroid = (
        sum(point[0] for point in triangle_points) / 3.0,
        sum(point[1] for point in triangle_points) / 3.0,
    )
    return [
        v,
        ((v[0] + nxt[0]) / 2.0, (v[1] + nxt[1]) / 2.0),
        centroid,
        ((v[0] + prv[0]) / 2.0, (v[1] + prv[1]) / 2.0),
    ]


def infer_profile(
    nodes: dict[int, dict[str, float]],
    net_doping: list[float],
    raw_profile: str,
) -> list[tuple[float, float, float]]:
    if raw_profile.strip() != "auto":
        segments = []
        for item in raw_profile.split(","):
            if not item.strip():
                continue
            x_min, x_max, value = item.split(":", 2)
            segments.append((float(x_min), float(x_max), float(value)))
        if not segments:
            raise RuntimeError("--profile must contain at least one segment")
        return segments

    by_x: dict[float, list[float]] = defaultdict(list)
    for node_id, node in nodes.items():
        by_x[float(node["x_um"])].append(net_doping[node_id])
    columns = [
        (x_um, statistics.median(values))
        for x_um, values in sorted(by_x.items())
    ]
    if len(columns) == 1:
        return [(-math.inf, math.inf, columns[0][1])]

    boundaries = [
        (columns[idx][0] + columns[idx + 1][0]) / 2.0
        for idx in range(len(columns) - 1)
    ]
    raw_segments = []
    for idx, (_x_um, value) in enumerate(columns):
        x_min = -math.inf if idx == 0 else boundaries[idx - 1]
        x_max = math.inf if idx == len(columns) - 1 else boundaries[idx]
        raw_segments.append((x_min, x_max, value))

    collapsed: list[tuple[float, float, float]] = []
    for x_min, x_max, value in raw_segments:
        if collapsed and math.isclose(collapsed[-1][2], value, rel_tol=1.0e-10, abs_tol=1.0e6):
            collapsed[-1] = (collapsed[-1][0], x_max, collapsed[-1][2])
        else:
            collapsed.append((x_min, x_max, value))
    return collapsed


def band_for_x(x_um: float, bands: list[tuple[str, float, float]]) -> str:
    for name, x_min, x_max in bands:
        if x_min <= x_um <= x_max:
            return name
    return "outside"


def integrate_profile_over_polygon(
    poly: list[tuple[float, float]],
    profile: list[tuple[float, float, float]],
) -> tuple[float, float]:
    area = polygon_area_m2(poly)
    integral = 0.0
    for x_min, x_max, value_cm3 in profile:
        clipped = clip_polygon_to_x_interval(poly, x_min, x_max)
        integral += polygon_area_m2(clipped) * value_cm3
    return area, integral


def compute_integrated_net(
    nodes: dict[int, dict[str, float]],
    triangles: list[dict[str, Any]],
    profile: list[tuple[float, float, float]],
) -> tuple[list[float], list[float]]:
    node_count = max(nodes) + 1
    volumes = [0.0 for _ in range(node_count)]
    integrals = [0.0 for _ in range(node_count)]
    for triangle in triangles:
        node_ids = triangle["node_ids"]
        points = [(nodes[node_id]["x_um"], nodes[node_id]["y_um"]) for node_id in node_ids]
        for local_index, node_id in enumerate(node_ids):
            poly = vertex_control_polygon(points, local_index)
            area, integral = integrate_profile_over_polygon(poly, profile)
            volumes[node_id] += area
            integrals[node_id] += integral
    integrated_net = [
        integrals[node_id] / volumes[node_id] if volumes[node_id] > 0.0 else 0.0
        for node_id in range(node_count)
    ]
    return volumes, integrated_net


def compute_mixed_voronoi_volumes(
    nodes: dict[int, dict[str, float]],
    triangles: list[dict[str, Any]],
) -> list[float]:
    volumes = [0.0 for _ in range(max(nodes) + 1)]
    for triangle in triangles:
        node_ids = triangle["node_ids"]
        points = [(nodes[node_id]["x_um"], nodes[node_id]["y_um"]) for node_id in node_ids]
        area = triangle_area_m2(points)
        angles = [
            angle_radians(points[0], points[1], points[2]),
            angle_radians(points[1], points[2], points[0]),
            angle_radians(points[2], points[0], points[1]),
        ]
        obtuse = [angle > math.pi / 2.0 + 1.0e-12 for angle in angles]
        if any(obtuse):
            for local_index, node_id in enumerate(node_ids):
                volumes[node_id] += area / 2.0 if obtuse[local_index] else area / 4.0
            continue
        for local_index, node_id in enumerate(node_ids):
            next_index = (local_index + 1) % 3
            prev_index = (local_index + 2) % 3
            share = (
                distance2_m2(points[local_index], points[next_index])
                * cotangent_opposite(points[local_index], points[next_index], points[prev_index])
                + distance2_m2(points[local_index], points[prev_index])
                * cotangent_opposite(points[local_index], points[prev_index], points[next_index])
            ) / 8.0
            volumes[node_id] += max(share, 0.0)
    return volumes


def median(values: list[float]) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return statistics.median(finite) if finite else None


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    args = parse_args()
    nodes, triangles, contact_by_node = read_mesh(args.mesh)
    node_count = len(nodes)
    donors, acceptors = read_doping(args.doping_csv, node_count)
    nodal_net = [donor - acceptor for donor, acceptor in zip(donors, acceptors)]
    bands = parse_bands(args.bands)
    profile = infer_profile(nodes, nodal_net, args.profile)
    volumes, integrated_net = compute_integrated_net(nodes, triangles, profile)
    mixed_volumes = compute_mixed_voronoi_volumes(nodes, triangles)

    node_rows = []
    for node_id in sorted(nodes):
        node = nodes[node_id]
        nodal_charge = Q * (nodal_net[node_id] * 1.0e6) * volumes[node_id]
        integrated_charge = Q * (integrated_net[node_id] * 1.0e6) * volumes[node_id]
        mixed_charge = Q * (nodal_net[node_id] * 1.0e6) * mixed_volumes[node_id]
        mixed_ratio = (
            mixed_volumes[node_id] / volumes[node_id]
            if volumes[node_id] > 0.0
            else math.nan
        )
        mixed_net = nodal_net[node_id] * mixed_ratio if math.isfinite(mixed_ratio) else math.nan
        ratio = (
            integrated_net[node_id] / nodal_net[node_id]
            if abs(nodal_net[node_id]) > 0.0
            else math.nan
        )
        node_rows.append({
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "band": band_for_x(node["x_um"], bands),
            "contact_names": ";".join(sorted(contact_by_node.get(node_id, set()))),
            "control_volume_m2": volumes[node_id],
            "mixed_voronoi_volume_m2": mixed_volumes[node_id],
            "mixed_to_barycentric_volume_ratio": mixed_ratio,
            "nodal_net_doping_cm3": nodal_net[node_id],
            "integrated_net_doping_cm3": integrated_net[node_id],
            "mixed_volume_net_doping_cm3": mixed_net,
            "delta_net_doping_cm3": integrated_net[node_id] - nodal_net[node_id],
            "mixed_volume_delta_net_doping_cm3": mixed_net - nodal_net[node_id],
            "integrated_to_nodal_ratio": ratio,
            "mixed_volume_to_nodal_ratio": (
                mixed_net / nodal_net[node_id] if abs(nodal_net[node_id]) > 0.0 else math.nan
            ),
            "nodal_fixed_charge_C_per_m": nodal_charge,
            "integrated_fixed_charge_C_per_m": integrated_charge,
            "mixed_volume_fixed_charge_C_per_m": mixed_charge,
            "delta_fixed_charge_C_per_m": integrated_charge - nodal_charge,
            "mixed_volume_delta_fixed_charge_C_per_m": mixed_charge - nodal_charge,
            "abs_delta_fixed_charge_C_per_m": abs(integrated_charge - nodal_charge),
            "abs_mixed_volume_delta_fixed_charge_C_per_m": abs(mixed_charge - nodal_charge),
        })

    band_rows = []
    for name, x_min, x_max in bands:
        selected = [row for row in node_rows if row["band"] == name]
        if selected:
            max_row = max(selected, key=lambda row: float(row["abs_delta_fixed_charge_C_per_m"]))
            max_mixed_row = max(
                selected,
                key=lambda row: float(row["abs_mixed_volume_delta_fixed_charge_C_per_m"]),
            )
        else:
            max_row = {}
            max_mixed_row = {}
        band_rows.append({
            "band": name,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "node_count": len(selected),
            "median_nodal_net_doping_cm3": median([float(row["nodal_net_doping_cm3"]) for row in selected]),
            "median_integrated_net_doping_cm3": median([float(row["integrated_net_doping_cm3"]) for row in selected]),
            "median_mixed_volume_net_doping_cm3": median([
                float(row["mixed_volume_net_doping_cm3"]) for row in selected
            ]),
            "median_delta_net_doping_cm3": median([float(row["delta_net_doping_cm3"]) for row in selected]),
            "median_mixed_volume_delta_net_doping_cm3": median([
                float(row["mixed_volume_delta_net_doping_cm3"]) for row in selected
            ]),
            "median_integrated_to_nodal_ratio": median([
                float(row["integrated_to_nodal_ratio"]) for row in selected
            ]),
            "median_mixed_volume_to_nodal_ratio": median([
                float(row["mixed_volume_to_nodal_ratio"]) for row in selected
            ]),
            "median_delta_fixed_charge_C_per_m": median([
                float(row["delta_fixed_charge_C_per_m"]) for row in selected
            ]),
            "median_mixed_volume_delta_fixed_charge_C_per_m": median([
                float(row["mixed_volume_delta_fixed_charge_C_per_m"]) for row in selected
            ]),
            "median_abs_delta_fixed_charge_C_per_m": median([
                float(row["abs_delta_fixed_charge_C_per_m"]) for row in selected
            ]),
            "median_abs_mixed_volume_delta_fixed_charge_C_per_m": median([
                float(row["abs_mixed_volume_delta_fixed_charge_C_per_m"]) for row in selected
            ]),
            "max_abs_delta_fixed_charge_C_per_m": max_row.get("abs_delta_fixed_charge_C_per_m"),
            "max_abs_delta_fixed_charge_node": max_row.get("node_id"),
            "max_abs_mixed_volume_delta_fixed_charge_C_per_m": max_mixed_row.get(
                "abs_mixed_volume_delta_fixed_charge_C_per_m"),
            "max_abs_mixed_volume_delta_fixed_charge_node": max_mixed_row.get("node_id"),
        })

    profile_rows = [
        {
            "segment": idx,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "net_doping_cm3": value,
        }
        for idx, (x_min, x_max, value) in enumerate(profile)
    ]
    integrated_doping_rows = [
        {
            "node_id": row["node_id"],
            "donors_cm3": max(float(row["integrated_net_doping_cm3"]), 0.0),
            "acceptors_cm3": max(-float(row["integrated_net_doping_cm3"]), 0.0),
        }
        for row in node_rows
    ]
    mixed_volume_doping_rows = [
        {
            "node_id": row["node_id"],
            "donors_cm3": max(float(row["mixed_volume_net_doping_cm3"]), 0.0),
            "acceptors_cm3": max(-float(row["mixed_volume_net_doping_cm3"]), 0.0),
        }
        for row in node_rows
    ]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "fixed_charge_rhs_nodes.csv", node_rows, NODE_FIELDS)
    write_csv(args.out_dir / "fixed_charge_rhs_bands.csv", band_rows, BAND_FIELDS)
    write_csv(args.out_dir / "fixed_charge_rhs_profile.csv", profile_rows, PROFILE_FIELDS)
    write_csv(
        args.out_dir / "integrated_doping.csv",
        integrated_doping_rows,
        ["node_id", "donors_cm3", "acceptors_cm3"],
    )
    write_csv(
        args.out_dir / "mixed_volume_doping.csv",
        mixed_volume_doping_rows,
        ["node_id", "donors_cm3", "acceptors_cm3"],
    )
    top_rows = sorted(
        node_rows,
        key=lambda row: float(row["abs_delta_fixed_charge_C_per_m"]),
        reverse=True,
    )[:25]
    (args.out_dir / "fixed_charge_rhs_summary.json").write_text(
        json.dumps(clean_json({
            "mesh": str(args.mesh),
            "doping_csv": str(args.doping_csv),
            "node_rows": len(node_rows),
            "band_rows": len(band_rows),
            "profile_segments": profile_rows,
            "top_abs_delta_fixed_charge_nodes": top_rows,
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
