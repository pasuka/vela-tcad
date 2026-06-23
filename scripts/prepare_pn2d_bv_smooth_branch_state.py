#!/usr/bin/env python3
"""Prepare a smooth PN2D BV branch-control initial state from a Vela VTK state.

This is a diagnostic helper. It shifts interior quasi-Fermi potentials with a
smooth support-distance weight while preserving contact Dirichlet nodes, then
writes a DCSweep-compatible initial_state_file.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


STATE_FIELDS = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
WEIGHT_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "weight",
    "distance_to_support_um",
    "selected_support",
    "contact_names",
    "electron_qf_shift_V",
    "hole_qf_shift_V",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--state-csv-name", default="smooth_branch_state.csv")
    parser.add_argument("--electron-qf-shift-v", type=float, required=True)
    parser.add_argument("--hole-qf-shift-v", type=float, required=True)
    parser.add_argument("--decay-length-um", type=float, required=True)
    parser.add_argument("--max-distance-um", type=float)
    parser.add_argument("--support-classes", default="false_negative")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def read_support(path: Path, selected_classes: set[str]) -> set[int]:
    selected: set[int] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cls = row.get("support_class", "")
            if cls in selected_classes:
                selected.add(int(row["node_id"]))
    return selected


def parse_classes(raw: str) -> set[str]:
    classes = {item.strip() for item in raw.split(",") if item.strip()}
    if not classes:
        raise SystemExit("--support-classes must select at least one class")
    return classes


def infer_ni(density: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if density <= 0.0:
        return 0.0
    if carrier == "electron":
        return density / sgdiag.limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return density / sgdiag.limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def shifted_density(ni: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if ni <= 0.0:
        return 0.0
    if carrier == "electron":
        return ni * sgdiag.limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return ni * sgdiag.limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def contact_names(contact_by_node: dict[int, set[str]], node_id: int) -> str:
    return ";".join(sorted(contact_by_node.get(node_id, set())))


def distance_um(nodes: dict[int, dict[str, float]], a: int, b: int) -> float:
    na = nodes[a]
    nb = nodes[b]
    return math.hypot(na["x_um"] - nb["x_um"], na["y_um"] - nb["y_um"])


def smooth_weight(distance: float, decay: float, max_distance: float | None) -> float:
    if max_distance is not None and distance > max_distance:
        return 0.0
    if distance == 0.0:
        return 1.0
    return math.exp(-((distance / decay) ** 2))


def load_vela_state(path: Path, node_count: int) -> dict[str, list[float]]:
    scalars = sgdiag.parse_vtk_scalars(path)
    mapping = {
        "psi": ("Potential", 1.0),
        "phin": ("ElectronQuasiFermi", 1.0),
        "phip": ("HoleQuasiFermi", 1.0),
        "n": ("Electrons", 1.0),
        "p": ("Holes", 1.0),
    }
    state: dict[str, list[float]] = {}
    for output, (source, scale) in mapping.items():
        if source not in scalars:
            raise RuntimeError(f"missing state field {source}")
        values = [value * scale for value in scalars[source]]
        if len(values) != node_count:
            raise RuntimeError(f"{source} has {len(values)} values, expected {node_count}")
        state[output] = values
    return state


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
    if args.decay_length_um <= 0.0 or not math.isfinite(args.decay_length_um):
        raise SystemExit("--decay-length-um must be positive and finite")
    max_distance = args.max_distance_um
    if max_distance is None:
        max_distance = 3.0 * args.decay_length_um
    if max_distance <= 0.0 or not math.isfinite(max_distance):
        raise SystemExit("--max-distance-um must be positive and finite")

    nodes, _, contact_by_node = sgdiag.read_mesh(args.mesh)
    selected_support = read_support(args.support_csv, parse_classes(args.support_classes))
    if not selected_support:
        raise SystemExit("selected support set is empty")
    missing = sorted(node_id for node_id in selected_support if node_id not in nodes)
    if missing:
        raise SystemExit(f"support nodes missing from mesh: {missing[:8]}")

    state = load_vela_state(args.vela_vtk, len(nodes))
    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    state_rows: list[dict[str, Any]] = []
    weight_rows: list[dict[str, Any]] = []
    contact_zero_count = 0
    nonzero_weight_count = 0
    max_weight = 0.0

    for node_id in range(len(nodes)):
        is_contact = bool(contact_by_node.get(node_id))
        nearest = min(distance_um(nodes, node_id, support_id) for support_id in selected_support)
        weight = 0.0 if is_contact else smooth_weight(nearest, args.decay_length_um, max_distance)
        if is_contact:
            contact_zero_count += 1
        if weight > 0.0:
            nonzero_weight_count += 1
        max_weight = max(max_weight, weight)

        psi = state["psi"][node_id]
        phin = state["phin"][node_id]
        phip = state["phip"][node_id]
        ni_e = infer_ni(state["n"][node_id], psi, phin, vt, "electron")
        ni_h = infer_ni(state["p"][node_id], psi, phip, vt, "hole")
        shifted_phin = phin + weight * args.electron_qf_shift_v
        shifted_phip = phip + weight * args.hole_qf_shift_v
        shifted_n = shifted_density(ni_e, psi, shifted_phin, vt, "electron")
        shifted_p = shifted_density(ni_h, psi, shifted_phip, vt, "hole")

        state_rows.append({
            "node_id": node_id,
            "psi": psi,
            "phin": shifted_phin,
            "phip": shifted_phip,
            "electrons_m3": shifted_n,
            "holes_m3": shifted_p,
        })
        node = nodes[node_id]
        weight_rows.append({
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "weight": weight,
            "distance_to_support_um": nearest,
            "selected_support": node_id in selected_support,
            "contact_names": contact_names(contact_by_node, node_id),
            "electron_qf_shift_V": weight * args.electron_qf_shift_v,
            "hole_qf_shift_V": weight * args.hole_qf_shift_v,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    state_csv = args.out_dir / args.state_csv_name
    write_csv(state_csv, state_rows, STATE_FIELDS)
    write_csv(args.out_dir / "smooth_branch_weights.csv", weight_rows, WEIGHT_FIELDS)
    summary = {
        "mesh": str(args.mesh),
        "vela_vtk": str(args.vela_vtk),
        "support_csv": str(args.support_csv),
        "state_csv": str(state_csv),
        "support_classes": sorted(parse_classes(args.support_classes)),
        "selected_support_node_count": len(selected_support),
        "contact_zero_weight_count": contact_zero_count,
        "nonzero_weight_node_count": nonzero_weight_count,
        "max_weight": max_weight,
        "decay_length_um": args.decay_length_um,
        "max_distance_um": max_distance,
        "electron_qf_shift_v": args.electron_qf_shift_v,
        "hole_qf_shift_v": args.hole_qf_shift_v,
    }
    (args.out_dir / "smooth_branch_state_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())