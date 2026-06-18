#!/usr/bin/env python3
"""Inspect PN2D junction node volumes and edge couplings near BV Poisson mismatch."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


NODE_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "selected_reason",
    "donors_cm3",
    "acceptors_cm3",
    "net_doping_cm3",
    "control_volume_m2",
    "control_volume_um2",
    "adjacent_cell_count",
    "incident_edge_count",
    "incident_couple_sum_m",
    "incident_couple_sum_um",
    "min_incident_edge_length_m",
    "max_incident_edge_length_m",
    "median_incident_edge_length_m",
    "contact_names",
]

EDGE_FIELDS = [
    "edge_id",
    "node0",
    "node1",
    "x0_um",
    "y0_um",
    "x1_um",
    "y1_um",
    "mid_x_um",
    "mid_y_um",
    "length_m",
    "couple_m",
    "couple_um",
    "couple_over_length",
    "adjacent_cell_count",
    "local_contribution_count",
    "negative_cotangent_count",
    "fallback_count",
    "local_couple_min_m",
    "local_couple_max_m",
    "local_couple_sum_m",
    "net_doping0_cm3",
    "net_doping1_cm3",
    "junction_crossing",
    "junction_touching",
    "contact_names",
    "selected_reason",
]

CELL_FIELDS = [
    "cell_id",
    "node_ids",
    "x_centroid_um",
    "y_centroid_um",
    "area_m2",
    "area_um2",
    "min_angle_deg",
    "max_angle_deg",
    "contains_focus_node",
    "contains_shoulder_x",
    "negative_cotangent_count",
    "fallback_count",
]

BAND_FIELDS = [
    "band",
    "x_min_um",
    "x_max_um",
    "node_count",
    "total_control_volume_m2",
    "median_control_volume_m2",
    "median_incident_couple_sum_m",
    "edge_count",
    "median_edge_couple_m",
    "median_couple_over_length",
    "negative_cotangent_edges",
    "fallback_edges",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path)
    parser.add_argument("--tdr-inventory-mesh", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--focus-nodes", default="955,1089")
    parser.add_argument("--x-window", default="0.82:1.18")
    parser.add_argument("--radius-um", type=float, default=0.08)
    parser.add_argument(
        "--bands",
        default="pre_shoulder:0.82:0.92,junction:0.95:1.05,post_shoulder:1.08:1.18",
    )
    return parser.parse_args()


def parse_nodes(raw: str) -> set[int]:
    return {int(item.strip()) for item in raw.split(",") if item.strip()}


def parse_window(raw: str) -> tuple[float, float]:
    left, right = raw.split(":", 1)
    return float(left), float(right)


def parse_bands(raw: str) -> list[tuple[str, float, float]]:
    result = []
    for item in raw.split(","):
        if not item.strip():
            continue
        name, x_min, x_max = item.split(":", 2)
        result.append((name.strip(), float(x_min), float(x_max)))
    return result


def read_mesh(path: Path) -> tuple[dict[int, dict[str, Any]], list[dict[str, Any]], dict[int, set[str]]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {
            "id": int(node["id"]),
            "x_um": float(node["x"]),
            "y_um": float(node["y"]),
        }
        for node in data["nodes"]
    }
    cells = [
        {
            "id": int(cell["id"]),
            "node_ids": [int(node_id) for node_id in cell["node_ids"]],
        }
        for cell in data["triangles"]
    ]
    contact_by_node: dict[int, set[str]] = defaultdict(set)
    for contact in data.get("contacts", []):
        name = str(contact.get("name", "contact"))
        for node_id in contact.get("node_ids", []):
            contact_by_node[int(node_id)].add(name)
    return nodes, cells, contact_by_node


def read_tdr_inventory_geometry(path: Path) -> tuple[list[tuple[float, float]], list[tuple[int, int, int]]]:
    data = json.loads(path.read_text())
    vertices = [(float(x), float(y)) for x, y in data["geometry"]["vertices"]]
    triangles: list[tuple[int, int, int]] = []
    for region in data["geometry"]["regions"]:
        for tri in region.get("triangles", []):
            triangles.append((int(tri[0]), int(tri[1]), int(tri[2])))
    return vertices, triangles


def load_doping(path: Path | None, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    if path is None:
        return donors, acceptors
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            donors[node_id] = float(row["donors_cm3"])
            acceptors[node_id] = float(row["acceptors_cm3"])
    return donors, acceptors


def triangle_area_m2(nodes: dict[int, dict[str, Any]], node_ids: list[int]) -> float:
    a, b, c = [nodes[node_id] for node_id in node_ids[:3]]
    ax, ay = a["x_um"] * 1.0e-6, a["y_um"] * 1.0e-6
    bx, by = b["x_um"] * 1.0e-6, b["y_um"] * 1.0e-6
    cx, cy = c["x_um"] * 1.0e-6, c["y_um"] * 1.0e-6
    return 0.5 * abs((bx - ax) * (cy - ay) - (cx - ax) * (by - ay))


def distance_m(nodes: dict[int, dict[str, Any]], node0: int, node1: int) -> float:
    n0 = nodes[node0]
    n1 = nodes[node1]
    dx = (n1["x_um"] - n0["x_um"]) * 1.0e-6
    dy = (n1["y_um"] - n0["y_um"]) * 1.0e-6
    return math.hypot(dx, dy)


def angle_deg_at(nodes: dict[int, dict[str, Any]], center: int, other0: int, other1: int) -> float:
    c = nodes[center]
    a = nodes[other0]
    b = nodes[other1]
    ux = (a["x_um"] - c["x_um"])
    uy = (a["y_um"] - c["y_um"])
    vx = (b["x_um"] - c["x_um"])
    vy = (b["y_um"] - c["y_um"])
    ul = math.hypot(ux, uy)
    vl = math.hypot(vx, vy)
    if ul <= 1.0e-30 or vl <= 1.0e-30:
        return 0.0
    cos_theta = max(-1.0, min(1.0, (ux * vx + uy * vy) / (ul * vl)))
    return math.degrees(math.acos(cos_theta))


def cotangent_at_opposite(nodes: dict[int, dict[str, Any]], a_id: int, b_id: int, opp_id: int) -> float:
    a = nodes[a_id]
    b = nodes[b_id]
    opp = nodes[opp_id]
    ux = (a["x_um"] - opp["x_um"]) * 1.0e-6
    uy = (a["y_um"] - opp["y_um"]) * 1.0e-6
    vx = (b["x_um"] - opp["x_um"]) * 1.0e-6
    vy = (b["y_um"] - opp["y_um"]) * 1.0e-6
    cross = ux * vy - uy * vx
    if abs(cross) <= 1.0e-30:
        return 0.0
    dot = ux * vx + uy * vy
    return dot / abs(cross)


def edge_key(node0: int, node1: int) -> tuple[int, int]:
    return (node0, node1) if node0 <= node1 else (node1, node0)


def build_geometry(nodes: dict[int, dict[str, Any]],
                   cells: list[dict[str, Any]]) -> tuple[list[float], list[dict[str, Any]], dict[int, list[int]]]:
    volumes = [0.0 for _ in range(max(nodes) + 1)]
    node_cells: dict[int, list[int]] = defaultdict(list)
    edge_index: dict[tuple[int, int], int] = {}
    edges: list[dict[str, Any]] = []
    for cell in cells:
        ids = cell["node_ids"]
        area = triangle_area_m2(nodes, ids)
        for node_id in ids:
            volumes[node_id] += area / 3.0
            node_cells[node_id].append(cell["id"])
        for k in range(3):
            a = ids[k]
            b = ids[(k + 1) % 3]
            key = edge_key(a, b)
            if key not in edge_index:
                edge_index[key] = len(edges)
                edges.append({
                    "edge_id": len(edges),
                    "node0": key[0],
                    "node1": key[1],
                    "length_m": distance_m(nodes, key[0], key[1]),
                    "couple_m": 0.0,
                    "cell_ids": [],
                    "local": [],
                })

    for cell in cells:
        ids = cell["node_ids"]
        area = triangle_area_m2(nodes, ids)
        for k in range(3):
            a = ids[k]
            b = ids[(k + 1) % 3]
            opp = ids[(k + 2) % 3]
            edge = edges[edge_index[edge_key(a, b)]]
            length = edge["length_m"]
            if length <= 1.0e-30:
                continue
            cot = cotangent_at_opposite(nodes, a, b, opp)
            local = 0.5 * cot * length
            fallback = False
            if cot < 0.0:
                local = area / (3.0 * length)
                fallback = True
            local = max(local, 0.0)
            edge["couple_m"] += local
            edge["cell_ids"].append(cell["id"])
            edge["local"].append({
                "cell_id": cell["id"],
                "opposite_node": opp,
                "cotangent": cot,
                "local_couple_m": local,
                "fallback": fallback,
            })
    return volumes, edges, node_cells


def selected_node_reasons(nodes: dict[int, dict[str, Any]],
                          focus: set[int],
                          x_window: tuple[float, float],
                          radius_um: float) -> dict[int, str]:
    reasons: dict[int, list[str]] = defaultdict(list)
    for node_id in focus:
        if node_id in nodes:
            reasons[node_id].append("focus")
    for node_id, node in nodes.items():
        x = node["x_um"]
        if x_window[0] <= x <= x_window[1]:
            reasons[node_id].append("x_window")
        for focus_id in focus:
            if focus_id not in nodes:
                continue
            focus_node = nodes[focus_id]
            distance = math.hypot(node["x_um"] - focus_node["x_um"], node["y_um"] - focus_node["y_um"])
            if distance <= radius_um:
                reasons[node_id].append(f"radius:{focus_id}")
    return {node_id: ";".join(values) for node_id, values in reasons.items()}


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def selected_edge_reason(edge: dict[str, Any], node_reasons: dict[int, str], x_window: tuple[float, float]) -> str:
    reasons = []
    if edge["node0"] in node_reasons or edge["node1"] in node_reasons:
        reasons.append("incident_selected_node")
    mid_x = 0.5 * (edge["x0_um"] + edge["x1_um"])
    if x_window[0] <= mid_x <= x_window[1]:
        reasons.append("mid_x_window")
    return ";".join(reasons)


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_node_rows(nodes: dict[int, dict[str, Any]],
                    volumes: list[float],
                    edges: list[dict[str, Any]],
                    node_cells: dict[int, list[int]],
                    contact_by_node: dict[int, set[str]],
                    node_reasons: dict[int, str],
                    donors: list[float],
                    acceptors: list[float]) -> list[dict[str, Any]]:
    incident: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for edge in edges:
        incident[edge["node0"]].append(edge)
        incident[edge["node1"]].append(edge)
    rows = []
    for node_id in sorted(node_reasons):
        item = nodes[node_id]
        incident_edges = incident[node_id]
        lengths = [edge["length_m"] for edge in incident_edges]
        couple_sum = sum(edge["couple_m"] for edge in incident_edges)
        rows.append({
            "node_id": node_id,
            "x_um": item["x_um"],
            "y_um": item["y_um"],
            "selected_reason": node_reasons[node_id],
            "donors_cm3": donors[node_id],
            "acceptors_cm3": acceptors[node_id],
            "net_doping_cm3": donors[node_id] - acceptors[node_id],
            "control_volume_m2": volumes[node_id],
            "control_volume_um2": volumes[node_id] * 1.0e12,
            "adjacent_cell_count": len(node_cells[node_id]),
            "incident_edge_count": len(incident_edges),
            "incident_couple_sum_m": couple_sum,
            "incident_couple_sum_um": couple_sum * 1.0e6,
            "min_incident_edge_length_m": min(lengths) if lengths else None,
            "max_incident_edge_length_m": max(lengths) if lengths else None,
            "median_incident_edge_length_m": median(lengths),
            "contact_names": ";".join(sorted(contact_by_node.get(node_id, set()))),
        })
    return rows


def build_edge_rows(nodes: dict[int, dict[str, Any]],
                    edges: list[dict[str, Any]],
                    contact_by_node: dict[int, set[str]],
                    node_reasons: dict[int, str],
                    x_window: tuple[float, float],
                    donors: list[float],
                    acceptors: list[float]) -> list[dict[str, Any]]:
    rows = []
    for edge in edges:
        n0 = nodes[edge["node0"]]
        n1 = nodes[edge["node1"]]
        edge["x0_um"] = n0["x_um"]
        edge["x1_um"] = n1["x_um"]
        edge["y0_um"] = n0["y_um"]
        edge["y1_um"] = n1["y_um"]
        reason = selected_edge_reason(edge, node_reasons, x_window)
        if not reason:
            continue
        local_values = [item["local_couple_m"] for item in edge["local"]]
        negative = sum(1 for item in edge["local"] if item["cotangent"] < 0.0)
        fallback = sum(1 for item in edge["local"] if item["fallback"])
        net0 = donors[edge["node0"]] - acceptors[edge["node0"]]
        net1 = donors[edge["node1"]] - acceptors[edge["node1"]]
        junction_crossing = (net0 < 0.0 < net1) or (net1 < 0.0 < net0)
        junction_touching = junction_crossing or (
            (net0 == 0.0 and net1 != 0.0) or (net1 == 0.0 and net0 != 0.0)
        )
        contacts = set(contact_by_node.get(edge["node0"], set())) | set(contact_by_node.get(edge["node1"], set()))
        rows.append({
            "edge_id": edge["edge_id"],
            "node0": edge["node0"],
            "node1": edge["node1"],
            "x0_um": n0["x_um"],
            "y0_um": n0["y_um"],
            "x1_um": n1["x_um"],
            "y1_um": n1["y_um"],
            "mid_x_um": 0.5 * (n0["x_um"] + n1["x_um"]),
            "mid_y_um": 0.5 * (n0["y_um"] + n1["y_um"]),
            "length_m": edge["length_m"],
            "couple_m": edge["couple_m"],
            "couple_um": edge["couple_m"] * 1.0e6,
            "couple_over_length": edge["couple_m"] / edge["length_m"] if edge["length_m"] > 0.0 else None,
            "adjacent_cell_count": len(edge["cell_ids"]),
            "local_contribution_count": len(edge["local"]),
            "negative_cotangent_count": negative,
            "fallback_count": fallback,
            "local_couple_min_m": min(local_values) if local_values else None,
            "local_couple_max_m": max(local_values) if local_values else None,
            "local_couple_sum_m": sum(local_values),
            "net_doping0_cm3": net0,
            "net_doping1_cm3": net1,
            "junction_crossing": junction_crossing,
            "junction_touching": junction_touching,
            "contact_names": ";".join(sorted(contacts)),
            "selected_reason": reason,
        })
    rows.sort(key=lambda row: (float(row["mid_x_um"]), float(row["mid_y_um"]), int(row["edge_id"])))
    return rows


def build_cell_rows(nodes: dict[int, dict[str, Any]],
                    cells: list[dict[str, Any]],
                    node_reasons: dict[int, str],
                    x_window: tuple[float, float]) -> list[dict[str, Any]]:
    rows = []
    for cell in cells:
        ids = cell["node_ids"]
        xs = [nodes[node_id]["x_um"] for node_id in ids]
        ys = [nodes[node_id]["y_um"] for node_id in ids]
        contains_focus = any("focus" in node_reasons.get(node_id, "") for node_id in ids)
        contains_window = any(x_window[0] <= x <= x_window[1] for x in xs)
        if not contains_focus and not contains_window:
            continue
        angles = [
            angle_deg_at(nodes, ids[0], ids[1], ids[2]),
            angle_deg_at(nodes, ids[1], ids[2], ids[0]),
            angle_deg_at(nodes, ids[2], ids[0], ids[1]),
        ]
        negative = 0
        fallback = 0
        for k in range(3):
            cot = cotangent_at_opposite(nodes, ids[k], ids[(k + 1) % 3], ids[(k + 2) % 3])
            if cot < 0.0:
                negative += 1
                fallback += 1
        area = triangle_area_m2(nodes, ids)
        rows.append({
            "cell_id": cell["id"],
            "node_ids": ";".join(str(node_id) for node_id in ids),
            "x_centroid_um": sum(xs) / 3.0,
            "y_centroid_um": sum(ys) / 3.0,
            "area_m2": area,
            "area_um2": area * 1.0e12,
            "min_angle_deg": min(angles),
            "max_angle_deg": max(angles),
            "contains_focus_node": contains_focus,
            "contains_shoulder_x": contains_window,
            "negative_cotangent_count": negative,
            "fallback_count": fallback,
        })
    rows.sort(key=lambda row: (float(row["x_centroid_um"]), float(row["y_centroid_um"]), int(row["cell_id"])))
    return rows


def build_band_rows(nodes: dict[int, dict[str, Any]],
                    node_rows: list[dict[str, Any]],
                    edge_rows: list[dict[str, Any]],
                    bands: list[tuple[str, float, float]]) -> list[dict[str, Any]]:
    del nodes
    rows = []
    for name, x_min, x_max in bands:
        selected_nodes = [
            row for row in node_rows
            if x_min <= float(row["x_um"]) <= x_max
        ]
        selected_edges = [
            row for row in edge_rows
            if x_min <= float(row["mid_x_um"]) <= x_max
        ]
        rows.append({
            "band": name,
            "x_min_um": x_min,
            "x_max_um": x_max,
            "node_count": len(selected_nodes),
            "total_control_volume_m2": sum(float(row["control_volume_m2"]) for row in selected_nodes),
            "median_control_volume_m2": median([float(row["control_volume_m2"]) for row in selected_nodes]),
            "median_incident_couple_sum_m": median([float(row["incident_couple_sum_m"]) for row in selected_nodes]),
            "edge_count": len(selected_edges),
            "median_edge_couple_m": median([float(row["couple_m"]) for row in selected_edges]),
            "median_couple_over_length": median([float(row["couple_over_length"]) for row in selected_edges]),
            "negative_cotangent_edges": sum(1 for row in selected_edges if int(row["negative_cotangent_count"]) > 0),
            "fallback_edges": sum(1 for row in selected_edges if int(row["fallback_count"]) > 0),
        })
    return rows


def compare_inventory(nodes: dict[int, dict[str, Any]],
                      cells: list[dict[str, Any]],
                      inventory_path: Path | None) -> dict[str, Any]:
    if inventory_path is None or not inventory_path.exists():
        return {"available": False}
    vertices, triangles = read_tdr_inventory_geometry(inventory_path)
    max_coord_delta = 0.0
    for node_id, node in nodes.items():
        if node_id >= len(vertices):
            continue
        x_um, y_um = vertices[node_id]
        max_coord_delta = max(
            max_coord_delta,
            abs(node["x_um"] - x_um),
            abs(node["y_um"] - y_um),
        )
    mesh_tris = {
        tuple(sorted(cell["node_ids"]))
        for cell in cells
    }
    inventory_tris = {tuple(sorted(tri)) for tri in triangles}
    return {
        "available": True,
        "vertex_count": len(vertices),
        "mesh_node_count": len(nodes),
        "triangle_count": len(triangles),
        "mesh_cell_count": len(cells),
        "max_coordinate_delta_um": max_coord_delta,
        "triangle_sets_match": mesh_tris == inventory_tris,
        "missing_in_mesh": len(inventory_tris - mesh_tris),
        "extra_in_mesh": len(mesh_tris - inventory_tris),
    }


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
    nodes, cells, contact_by_node = read_mesh(args.mesh)
    focus = parse_nodes(args.focus_nodes)
    x_window = parse_window(args.x_window)
    bands = parse_bands(args.bands)
    donors, acceptors = load_doping(args.doping_csv, len(nodes))
    volumes, edges, node_cells = build_geometry(nodes, cells)
    node_reasons = selected_node_reasons(nodes, focus, x_window, args.radius_um)
    node_rows = build_node_rows(
        nodes, volumes, edges, node_cells, contact_by_node, node_reasons, donors, acceptors)
    edge_rows = build_edge_rows(
        nodes, edges, contact_by_node, node_reasons, x_window, donors, acceptors)
    cell_rows = build_cell_rows(nodes, cells, node_reasons, x_window)
    band_rows = build_band_rows(nodes, node_rows, edge_rows, bands)
    inventory_comparison = compare_inventory(nodes, cells, args.tdr_inventory_mesh)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "junction_geometry_nodes.csv", node_rows, NODE_FIELDS)
    write_csv(args.out_dir / "junction_geometry_edges.csv", edge_rows, EDGE_FIELDS)
    write_csv(args.out_dir / "junction_geometry_cells.csv", cell_rows, CELL_FIELDS)
    write_csv(args.out_dir / "junction_geometry_bands.csv", band_rows, BAND_FIELDS)
    (args.out_dir / "junction_geometry_summary.json").write_text(
        json.dumps(clean_json({
            "mesh": str(args.mesh),
            "doping_csv": str(args.doping_csv) if args.doping_csv else None,
            "focus_nodes": sorted(focus),
            "x_window_um": {"min": x_window[0], "max": x_window[1]},
            "radius_um": args.radius_um,
            "selected_node_count": len(node_rows),
            "selected_edge_count": len(edge_rows),
            "selected_cell_count": len(cell_rows),
            "total_control_volume_m2_selected": sum(float(row["control_volume_m2"]) for row in node_rows),
            "negative_cotangent_edges_selected": sum(1 for row in edge_rows if int(row["negative_cotangent_count"]) > 0),
            "fallback_edges_selected": sum(1 for row in edge_rows if int(row["fallback_count"]) > 0),
            "inventory_comparison": inventory_comparison,
            "bands": band_rows,
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
