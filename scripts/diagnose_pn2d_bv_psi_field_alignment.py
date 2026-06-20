#!/usr/bin/env python3
"""Compare Vela and Sentaurus electrostatic potential around PN2D BV active nodes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from collections import deque
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "bias_V",
    "node_id",
    "node_class",
    "neighbor_depth",
    "x_um",
    "y_um",
    "contact_names",
    "donors_m3",
    "acceptors_m3",
    "net_doping_m3",
    "total_doping_m3",
    "vela_psi_V",
    "sentaurus_psi_V",
    "psi_delta_V",
    "psi_delta_residual_after_active_median_V",
    "vela_electron_density_m3",
    "sentaurus_electron_density_m3",
    "vela_hole_density_m3",
    "sentaurus_hole_density_m3",
    "ni_model_m3",
    "vela_inferred_ni_electron_m3",
    "vela_inferred_ni_hole_m3",
    "sentaurus_inferred_ni_electron_m3",
    "sentaurus_inferred_ni_hole_m3",
    "vela_inferred_ni_electron_over_model",
    "vela_inferred_ni_hole_over_model",
    "sentaurus_inferred_ni_electron_over_model",
    "sentaurus_inferred_ni_hole_over_model",
    "vela_net_charge_proxy_m3",
    "sentaurus_net_charge_proxy_m3",
    "net_charge_proxy_delta_m3",
]

CORRELATION_FIELDS = [
    "psi_delta_V",
    "psi_delta_residual_after_active_median_V",
    "net_doping_m3",
    "total_doping_m3",
    "net_charge_proxy_delta_m3",
    "vela_inferred_ni_hole_over_model",
    "sentaurus_inferred_ni_hole_over_model",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--edge-state-csv", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path)
    parser.add_argument("--vela-vtk-root", type=Path)
    parser.add_argument("--sentaurus-dir", type=Path)
    parser.add_argument("--sentaurus-root", type=Path)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--support-class", default="overlap")
    parser.add_argument("--neighbor-depth", type=int, default=1)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def read_csv(path: Path) -> list[dict[str, str]]:
    with open_path(path, "r", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(raw: Any, default: float = 0.0) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def matches_bias(row: dict[str, str], bias: float) -> bool:
    if row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row.get("bias_V")) - bias) <= 1.0e-9


def resolve_vela_vtk(args: argparse.Namespace) -> Path:
    if args.vela_vtk is not None:
        return args.vela_vtk
    if args.vela_vtk_root is None:
        raise SystemExit("provide --vela-vtk or --vela-vtk-root")
    vtk = fluxforms.discover_vela_vtks(args.vela_vtk_root).get(fluxforms.bias_key(args.bias))
    if vtk is None:
        raise SystemExit(f"missing Vela VTK for {args.bias:g} V under {args.vela_vtk_root}")
    return vtk


def resolve_sentaurus_dir(args: argparse.Namespace) -> Path:
    if args.sentaurus_dir is not None:
        return args.sentaurus_dir
    if args.sentaurus_root is None:
        raise SystemExit("provide --sentaurus-dir or --sentaurus-root")
    return fluxforms.sentaurus_dir(args.sentaurus_root, args.bias)


def adjacency_from_triangles(triangles: list[dict[str, Any]]) -> dict[int, set[int]]:
    adjacency: dict[int, set[int]] = {}
    for triangle in triangles:
        ids = [int(node_id) for node_id in triangle["node_ids"]]
        for index, node_id in enumerate(ids):
            adjacency.setdefault(node_id, set())
            for other in (ids[(index + 1) % 3], ids[(index + 2) % 3]):
                adjacency[node_id].add(int(other))
    return adjacency


def edge_nodes_by_id(nodes: dict[int, dict[str, float]], triangles: list[dict[str, Any]]) -> dict[int, tuple[int, int]]:
    result: dict[int, tuple[int, int]] = {}
    for edge in sgdiag.build_edges(nodes, triangles):
        result[int(edge["edge_id"])] = (int(edge["node0"]), int(edge["node1"]))
    return result


def active_seed_nodes(
    edge_state_csv: Path,
    nodes: dict[int, dict[str, float]],
    triangles: list[dict[str, Any]],
    variant: str,
    support_class: str,
    bias: float,
) -> tuple[set[int], set[int]]:
    support_nodes: set[int] = set()
    endpoint_nodes: set[int] = set()
    edge_lookup = edge_nodes_by_id(nodes, triangles)
    for row in read_csv(edge_state_csv):
        if row.get("variant") != variant:
            continue
        if row.get("support_class") != support_class:
            continue
        if not matches_bias(row, bias):
            continue
        support_raw = row.get("support_node_id", "")
        if support_raw != "":
            support_nodes.add(int(float(support_raw)))
        edge_raw = row.get("edge_id", "")
        if edge_raw != "":
            edge = edge_lookup.get(int(float(edge_raw)))
            if edge is not None:
                endpoint_nodes.update(edge)
    return support_nodes, endpoint_nodes


def selected_nodes_by_depth(adjacency: dict[int, set[int]], seeds: set[int], max_depth: int) -> dict[int, int]:
    depths = {node_id: 0 for node_id in seeds}
    queue: deque[int] = deque(sorted(seeds))
    while queue:
        node_id = queue.popleft()
        depth = depths[node_id]
        if depth >= max_depth:
            continue
        for neighbor in sorted(adjacency.get(node_id, set())):
            if neighbor in depths:
                continue
            depths[neighbor] = depth + 1
            queue.append(neighbor)
    return depths


def read_doping(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    seen = [False for _ in range(node_count)]
    for row in read_csv(path):
        node_id = int(row["node_id"])
        if "donors_m3" in row or "acceptors_m3" in row:
            donor = finite_float(row.get("donors_m3"))
            acceptor = finite_float(row.get("acceptors_m3"))
        else:
            donor = finite_float(row.get("donors_cm3")) * 1.0e6
            acceptor = finite_float(row.get("acceptors_cm3")) * 1.0e6
        donors[node_id] = donor
        acceptors[node_id] = acceptor
        seen[node_id] = True
    missing = [index for index, present in enumerate(seen) if not present]
    if missing:
        raise RuntimeError(f"doping CSV is missing node ids: {missing[:8]}")
    return donors, acceptors


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator in (None, 0.0):
        return None
    return numerator / denominator


def inferred_ni(state: dict[str, list[float]], node_id: int, vt: float, carrier: str) -> float:
    if carrier == "electron":
        density = state["n"][node_id]
        if density <= 0.0:
            return 0.0
        return density * sgdiag.limited_exp((state["phin"][node_id] - state["psi"][node_id]) / vt)
    if carrier == "hole":
        density = state["p"][node_id]
        if density <= 0.0:
            return 0.0
        return density * sgdiag.limited_exp((state["psi"][node_id] - state["phip"][node_id]) / vt)
    raise ValueError(carrier)


def node_class(node_id: int, support_nodes: set[int], endpoint_nodes: set[int], depth: int) -> str:
    if node_id in support_nodes:
        return "active_support"
    if node_id in endpoint_nodes:
        return "active_edge_endpoint"
    return f"neighbor_depth_{depth}"


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, contact_by_node = sgdiag.read_mesh(args.mesh)
    node_count = max(nodes) + 1
    support_nodes, endpoint_nodes = active_seed_nodes(
        args.edge_state_csv, nodes, triangles, args.variant, args.support_class, args.bias)
    seeds = set(support_nodes) | set(endpoint_nodes)
    if not seeds:
        raise SystemExit("no active seed nodes matched the requested filters")
    adjacency = adjacency_from_triangles(triangles)
    selected_depths = selected_nodes_by_depth(adjacency, seeds, args.neighbor_depth)
    donors, acceptors = read_doping(args.doping_csv, node_count)
    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    ni_model = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        node_count,
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    vela = fluxforms.load_vela_state(resolve_vela_vtk(args), node_count)
    sentaurus = fluxforms.load_sentaurus_state(resolve_sentaurus_dir(args), node_count)
    active_deltas = [
        sentaurus["psi"][node_id] - vela["psi"][node_id]
        for node_id in sorted(seeds)
    ]
    active_median = median(active_deltas) or 0.0
    rows: list[dict[str, Any]] = []
    for node_id, depth in sorted(selected_depths.items(), key=lambda item: (item[1], nodes[item[0]]["y_um"], nodes[item[0]]["x_um"], item[0])):
        donor = donors[node_id]
        acceptor = acceptors[node_id]
        vela_charge = vela["p"][node_id] - vela["n"][node_id] + donor - acceptor
        sentaurus_charge = sentaurus["p"][node_id] - sentaurus["n"][node_id] + donor - acceptor
        psi_delta = sentaurus["psi"][node_id] - vela["psi"][node_id]
        ni_node = ni_model[node_id]
        vela_ni_e = inferred_ni(vela, node_id, vt, "electron")
        vela_ni_h = inferred_ni(vela, node_id, vt, "hole")
        sentaurus_ni_e = inferred_ni(sentaurus, node_id, vt, "electron")
        sentaurus_ni_h = inferred_ni(sentaurus, node_id, vt, "hole")
        rows.append({
            "bias_V": args.bias,
            "node_id": node_id,
            "node_class": node_class(node_id, support_nodes, endpoint_nodes, depth),
            "neighbor_depth": depth,
            "x_um": nodes[node_id]["x_um"],
            "y_um": nodes[node_id]["y_um"],
            "contact_names": ",".join(sorted(contact_by_node.get(node_id, set()))),
            "donors_m3": donor,
            "acceptors_m3": acceptor,
            "net_doping_m3": donor - acceptor,
            "total_doping_m3": donor + acceptor,
            "vela_psi_V": vela["psi"][node_id],
            "sentaurus_psi_V": sentaurus["psi"][node_id],
            "psi_delta_V": psi_delta,
            "psi_delta_residual_after_active_median_V": psi_delta - active_median,
            "vela_electron_density_m3": vela["n"][node_id],
            "sentaurus_electron_density_m3": sentaurus["n"][node_id],
            "vela_hole_density_m3": vela["p"][node_id],
            "sentaurus_hole_density_m3": sentaurus["p"][node_id],
            "ni_model_m3": ni_node,
            "vela_inferred_ni_electron_m3": vela_ni_e,
            "vela_inferred_ni_hole_m3": vela_ni_h,
            "sentaurus_inferred_ni_electron_m3": sentaurus_ni_e,
            "sentaurus_inferred_ni_hole_m3": sentaurus_ni_h,
            "vela_inferred_ni_electron_over_model": safe_ratio(vela_ni_e, ni_node),
            "vela_inferred_ni_hole_over_model": safe_ratio(vela_ni_h, ni_node),
            "sentaurus_inferred_ni_electron_over_model": safe_ratio(sentaurus_ni_e, ni_node),
            "sentaurus_inferred_ni_hole_over_model": safe_ratio(sentaurus_ni_h, ni_node),
            "vela_net_charge_proxy_m3": vela_charge,
            "sentaurus_net_charge_proxy_m3": sentaurus_charge,
            "net_charge_proxy_delta_m3": sentaurus_charge - vela_charge,
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def pearson(xs: list[float], ys: list[float]) -> float | None:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return None
    mean_x = statistics.fmean(x for x, _ in pairs)
    mean_y = statistics.fmean(y for _, y in pairs)
    dx = [x - mean_x for x, _ in pairs]
    dy = [y - mean_y for _, y in pairs]
    denom = math.sqrt(sum(x * x for x in dx) * sum(y * y for y in dy))
    if denom == 0.0:
        return None
    return sum(x * y for x, y in zip(dx, dy)) / denom


def metric_summary(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"median": None, "min": None, "max": None}
    return {"median": statistics.median(values), "min": min(values), "max": max(values)}


def summarize_by_class(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    result: dict[str, dict[str, Any]] = {}
    for node_class_name in sorted({str(row.get("node_class", "")) for row in rows}):
        class_rows = [row for row in rows if row.get("node_class") == node_class_name]
        result[node_class_name] = {
            "count": len(class_rows),
            "psi_delta_V": metric_summary(finite_values(class_rows, "psi_delta_V")),
            "psi_delta_residual_after_active_median_V": metric_summary(
                finite_values(class_rows, "psi_delta_residual_after_active_median_V")
            ),
        }
    return result


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    y = finite_values(rows, "y_um")
    correlations: dict[str, dict[str, Any]] = {}
    for key in CORRELATION_FIELDS:
        values = finite_values(rows, key)
        correlations[key] = {"corr_y": pearson(y, values), **metric_summary(values)}
    active_deltas = [
        finite_float(row.get("psi_delta_V"))
        for row in rows
        if row.get("neighbor_depth") == 0
    ]
    return {
        "row_count": len(rows),
        "node_classes": summarize_by_class(rows),
        "active_seed_psi_delta_median_V": median(active_deltas),
        "correlations": correlations,
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
    rows = build_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "psi_field_alignment_nodes.csv", rows)
    with open_path(args.out_dir / "psi_field_alignment_summary.json", "w") as handle:
        json.dump(clean_json(summarize(rows)), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
