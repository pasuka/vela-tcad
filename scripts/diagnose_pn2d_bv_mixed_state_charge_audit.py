#!/usr/bin/env python3
"""Audit Poisson/space-charge impact of mixed high-density active support states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


Q = 1.602176634e-19

NODE_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "selected_support",
    "mixed_source",
    "contact_names",
    "distance_to_contact_um",
    "contact_bucket",
    "doping_sign",
    "node_volume_m2",
    "donors_m3",
    "acceptors_m3",
    "net_doping_m3",
    "vela_electrons_m3",
    "vela_holes_m3",
    "sentaurus_electrons_m3",
    "sentaurus_holes_m3",
    "mixed_electrons_m3",
    "mixed_holes_m3",
    "delta_electrons_m3",
    "delta_holes_m3",
    "baseline_mobile_charge_density_C_m3",
    "mixed_mobile_charge_density_C_m3",
    "delta_mobile_charge_density_C_m3",
    "baseline_net_charge_density_C_m3",
    "mixed_net_charge_density_C_m3",
    "delta_net_charge_density_C_m3",
    "baseline_mobile_charge_C_per_m",
    "mixed_mobile_charge_C_per_m",
    "delta_mobile_charge_C_per_m",
    "baseline_net_charge_C_per_m",
    "mixed_net_charge_C_per_m",
    "delta_net_charge_C_per_m",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--support-classes", default="false_negative")
    parser.add_argument("--near-contact-um", type=float, default=0.05)
    return parser.parse_args()


def parse_classes(raw: str) -> set[str]:
    classes = {item.strip() for item in raw.split(",") if item.strip()}
    if not classes:
        raise SystemExit("--support-classes must select at least one class")
    return classes


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_support(path: Path) -> dict[int, str]:
    result: dict[int, str] = {}
    for row in read_csv_rows(path):
        cls = row.get("support_class", "")
        if cls and cls != "inactive":
            result[int(row["node_id"])] = cls
    return result


def read_doping(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    donors = [0.0 for _ in range(node_count)]
    acceptors = [0.0 for _ in range(node_count)]
    seen = [False for _ in range(node_count)]
    for row in read_csv_rows(path):
        node_id = int(row["node_id"])
        if "donors_m3" in row or "acceptors_m3" in row:
            donor = float(row.get("donors_m3", 0.0) or 0.0)
            acceptor = float(row.get("acceptors_m3", 0.0) or 0.0)
        else:
            donor = float(row.get("donors_cm3", 0.0) or 0.0) * 1.0e6
            acceptor = float(row.get("acceptors_cm3", 0.0) or 0.0) * 1.0e6
        donors[node_id] = donor
        acceptors[node_id] = acceptor
        seen[node_id] = True
    missing = [idx for idx, present in enumerate(seen) if not present]
    if missing:
        raise RuntimeError(f"doping CSV is missing node ids: {missing[:8]}")
    return donors, acceptors


def read_scalar_csv_m3(path: Path, scale: float) -> list[float]:
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        value_col = next((name for name in reader.fieldnames or [] if name != "node_id"), None)
        if value_col is None:
            raise RuntimeError(f"{path} has no value column")
        for row in reader:
            values[int(row["node_id"])] = float(row[value_col]) * scale
    return [values[idx] for idx in range(max(values) + 1)]


def load_vela_densities(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    scalars = sgdiag.parse_vtk_scalars(path)
    for name in ("Electrons", "Holes"):
        if name not in scalars:
            raise RuntimeError(f"VTK is missing required scalar {name}")
        if len(scalars[name]) != node_count:
            raise RuntimeError(f"{name} has {len(scalars[name])} values, expected {node_count}")
    return list(scalars["Electrons"]), list(scalars["Holes"])


def load_sentaurus_densities(path: Path, node_count: int) -> tuple[list[float], list[float]]:
    electrons = read_scalar_csv_m3(path / "fields" / "eDensity_region0.csv", 1.0e6)
    holes = read_scalar_csv_m3(path / "fields" / "hDensity_region0.csv", 1.0e6)
    if len(electrons) != node_count or len(holes) != node_count:
        raise RuntimeError("Sentaurus density field length does not match mesh node count")
    return electrons, holes


def contact_names(contact_by_node: dict[int, set[str]], node_id: int) -> str:
    return ";".join(sorted(contact_by_node.get(node_id, set())))


def distance_um(nodes: dict[int, dict[str, float]], a: int, b: int) -> float:
    na = nodes[a]
    nb = nodes[b]
    return math.hypot(na["x_um"] - nb["x_um"], na["y_um"] - nb["y_um"])


def contact_distance_um(nodes: dict[int, dict[str, float]], contact_nodes: list[int], node_id: int) -> float | None:
    if not contact_nodes:
        return None
    return min(distance_um(nodes, node_id, contact_id) for contact_id in contact_nodes)


def contact_bucket(distance: float | None, is_contact: bool, near_contact_um: float) -> str:
    if is_contact:
        return "contact"
    if distance is not None and distance <= near_contact_um:
        return "near_contact"
    return "interior"


def doping_sign(net: float) -> str:
    if net > 0.0:
        return "n_type"
    if net < 0.0:
        return "p_type"
    return "compensated"


def charge_density(electrons: float, holes: float, donors: float, acceptors: float) -> tuple[float, float]:
    mobile = Q * (holes - electrons)
    net = Q * (holes - electrons + donors - acceptors)
    return mobile, net


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, contact_by_node = sgdiag.read_mesh(args.mesh)
    node_count = len(nodes)
    volumes = sgdiag.node_volumes(nodes, triangles)
    donors, acceptors = read_doping(args.doping_csv, node_count)
    vela_n, vela_p = load_vela_densities(args.vela_vtk, node_count)
    sent_n, sent_p = load_sentaurus_densities(args.sentaurus_dir, node_count)
    support = read_support(args.support_csv)
    selected_classes = parse_classes(args.support_classes)
    contact_nodes = sorted(contact_by_node)

    rows: list[dict[str, Any]] = []
    for node_id in sorted(support):
        cls = support[node_id]
        selected = cls in selected_classes
        mixed_n = sent_n[node_id] if selected else vela_n[node_id]
        mixed_p = sent_p[node_id] if selected else vela_p[node_id]
        base_mobile_density, base_net_density = charge_density(
            vela_n[node_id], vela_p[node_id], donors[node_id], acceptors[node_id])
        mixed_mobile_density, mixed_net_density = charge_density(
            mixed_n, mixed_p, donors[node_id], acceptors[node_id])
        volume = volumes[node_id]
        net_doping = donors[node_id] - acceptors[node_id]
        distance = contact_distance_um(nodes, contact_nodes, node_id)
        contacts = contact_names(contact_by_node, node_id)
        rows.append({
            "node_id": node_id,
            "x_um": nodes[node_id]["x_um"],
            "y_um": nodes[node_id]["y_um"],
            "support_class": cls,
            "selected_support": selected,
            "mixed_source": "sentaurus_density" if selected else "vela_density",
            "contact_names": contacts,
            "distance_to_contact_um": distance,
            "contact_bucket": contact_bucket(distance, bool(contacts), args.near_contact_um),
            "doping_sign": doping_sign(net_doping),
            "node_volume_m2": volume,
            "donors_m3": donors[node_id],
            "acceptors_m3": acceptors[node_id],
            "net_doping_m3": net_doping,
            "vela_electrons_m3": vela_n[node_id],
            "vela_holes_m3": vela_p[node_id],
            "sentaurus_electrons_m3": sent_n[node_id],
            "sentaurus_holes_m3": sent_p[node_id],
            "mixed_electrons_m3": mixed_n,
            "mixed_holes_m3": mixed_p,
            "delta_electrons_m3": mixed_n - vela_n[node_id],
            "delta_holes_m3": mixed_p - vela_p[node_id],
            "baseline_mobile_charge_density_C_m3": base_mobile_density,
            "mixed_mobile_charge_density_C_m3": mixed_mobile_density,
            "delta_mobile_charge_density_C_m3": mixed_mobile_density - base_mobile_density,
            "baseline_net_charge_density_C_m3": base_net_density,
            "mixed_net_charge_density_C_m3": mixed_net_density,
            "delta_net_charge_density_C_m3": mixed_net_density - base_net_density,
            "baseline_mobile_charge_C_per_m": base_mobile_density * volume,
            "mixed_mobile_charge_C_per_m": mixed_mobile_density * volume,
            "delta_mobile_charge_C_per_m": (mixed_mobile_density - base_mobile_density) * volume,
            "baseline_net_charge_C_per_m": base_net_density * volume,
            "mixed_net_charge_C_per_m": mixed_net_density * volume,
            "delta_net_charge_C_per_m": (mixed_net_density - base_net_density) * volume,
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = row.get(key)
        if isinstance(value, (int, float)) and math.isfinite(float(value)):
            values.append(float(value))
    return values


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    keys = [
        "delta_mobile_charge_C_per_m",
        "delta_net_charge_C_per_m",
        "baseline_net_charge_C_per_m",
        "mixed_net_charge_C_per_m",
        "node_volume_m2",
    ]
    result: dict[str, Any] = {"node_count": len(rows)}
    for key in keys:
        values = finite_values(rows, key)
        result[f"sum_{key}"] = sum(values)
        result[f"sum_abs_{key}"] = sum(abs(value) for value in values)
        result[f"median_{key}"] = statistics.median(values) if values else None
        result[f"max_abs_{key}"] = max((abs(value) for value in values), default=None)
    return result


def grouped_summary(rows: list[dict[str, Any]], key: str) -> dict[str, Any]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        groups[str(row[key])].append(row)
    return {group: summarize_rows(group_rows) for group, group_rows in sorted(groups.items())}


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    summary = {
        "mesh": str(args.mesh),
        "doping_csv": str(args.doping_csv),
        "support_csv": str(args.support_csv),
        "vela_vtk": str(args.vela_vtk),
        "sentaurus_dir": str(args.sentaurus_dir),
        "selected_support_classes": sorted(parse_classes(args.support_classes)),
        "all": summarize_rows(rows),
        "selected": summarize_rows([row for row in rows if row["selected_support"]]),
        "support_classes": grouped_summary(rows, "support_class"),
        "doping_sign": grouped_summary(rows, "doping_sign"),
        "contact_bucket": grouped_summary(rows, "contact_bucket"),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "mixed_state_charge_audit_nodes.csv", rows)
    (args.out_dir / "mixed_state_charge_audit_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())