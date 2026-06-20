#!/usr/bin/env python3
"""Diagnose whether contact Dirichlet projection erases tiny QF drops."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


EDGE_FIELDS = [
    "bias_V",
    "contact",
    "edge_id",
    "node0",
    "node1",
    "contact_node",
    "interior_node",
    "outward_sign",
    "restart_contact_phip_V",
    "restart_interior_phip_V",
    "restart_edge_order_phip_drop_V",
    "restart_interior_minus_contact_phip_V",
    "projected_contact_phip_V",
    "projected_interior_phip_V",
    "projected_edge_order_phip_drop_V",
    "projected_interior_minus_contact_phip_V",
    "final_contact_phip_V",
    "final_interior_phip_V",
    "final_edge_order_phip_drop_V",
    "final_interior_minus_contact_phip_V",
    "reported_edge_order_phip_drop_V",
    "drop_source",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--restart-state", type=Path, required=True)
    parser.add_argument("--final-state", type=Path, required=True)
    parser.add_argument("--contact-edge-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument(
        "--zero-atol",
        type=float,
        default=1.0e-18,
        help="Absolute drop threshold treated as zero.",
    )
    return parser.parse_args()


def optional_float(raw: object) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value):
        return None
    return value


def optional_int(raw: object) -> int | None:
    value = optional_float(raw)
    if value is None:
        return None
    return int(value)


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_state(path: Path) -> dict[int, dict[str, float]]:
    state: dict[int, dict[str, float]] = {}
    for row in load_csv(path):
        node_id = optional_int(row.get("node_id"))
        if node_id is None:
            continue
        values: dict[str, float] = {}
        for key in ("psi", "phin", "phip"):
            value = optional_float(row.get(key))
            if value is not None:
                values[key] = value
        state[node_id] = values
    return state


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def near_bias(row: dict[str, str], bias: float) -> bool:
    row_bias = optional_float(row.get("bias_V"))
    return row_bias is not None and abs(row_bias - bias) <= 1.0e-9


def edge_nodes(row: dict[str, str]) -> tuple[int, int, int, int] | None:
    node0 = optional_int(row.get("node0"))
    node1 = optional_int(row.get("node1"))
    outward = optional_float(row.get("outward_sign"))
    if node0 is None or node1 is None or outward is None:
        return None
    if outward > 0.0:
        return node0, node1, node0, node1
    return node0, node1, node1, node0


def phip(state: dict[int, dict[str, float]], node: int) -> float | None:
    return state.get(node, {}).get("phip")


def edge_order_drop(
    node0: int,
    node1: int,
    phip0: float | None,
    phip1: float | None,
) -> float | None:
    del node0, node1
    if phip0 is None or phip1 is None:
        return None
    return phip1 - phip0


def classify_drop(
    restart_drop: float | None,
    projected_drop: float | None,
    final_drop: float | None,
    reported_drop: float | None,
    zero_atol: float,
) -> str:
    restart_nonzero = restart_drop is not None and abs(restart_drop) > zero_atol
    projected_zero = projected_drop is not None and abs(projected_drop) <= zero_atol
    final_zero = final_drop is not None and abs(final_drop) <= zero_atol
    projected_nonzero = projected_drop is not None and abs(projected_drop) > zero_atol
    final_nonzero = final_drop is not None and abs(final_drop) > zero_atol
    reported_disagrees = (
        reported_drop is not None and final_drop is not None
        and abs(reported_drop - final_drop) > zero_atol
    )
    if restart_nonzero and projected_zero and final_zero:
        return "contact_projection_erased_restart_drop"
    if restart_nonzero and projected_nonzero and final_zero:
        return "solver_convergence_erased_projected_drop"
    if final_nonzero:
        return "final_state_preserved_drop"
    if reported_disagrees:
        return "contact_edge_reporting_mismatch"
    if not restart_nonzero:
        return "restart_already_flat"
    return "undetermined"


def max_abs(values: list[float | None]) -> float:
    finite = [abs(value) for value in values if value is not None and math.isfinite(value)]
    return max(finite) if finite else 0.0


def main() -> None:
    args = parse_args()
    restart = load_state(args.restart_state)
    final = load_state(args.final_state)
    rows = [
        row for row in load_csv(args.contact_edge_csv)
        if row.get("current_contact") == args.contact and near_bias(row, args.bias)
    ]

    out_rows: list[dict[str, Any]] = []
    restart_drops: list[float | None] = []
    projected_drops: list[float | None] = []
    final_drops: list[float | None] = []
    drop_sources: Counter[str] = Counter()

    for row in rows:
        nodes = edge_nodes(row)
        if nodes is None:
            continue
        node0, node1, contact_node, interior_node = nodes
        outward = optional_float(row.get("outward_sign")) or 0.0

        restart0 = phip(restart, node0)
        restart1 = phip(restart, node1)
        final0 = phip(final, node0)
        final1 = phip(final, node1)
        restart_contact = phip(restart, contact_node)
        restart_interior = phip(restart, interior_node)
        final_contact = phip(final, contact_node)
        final_interior = phip(final, interior_node)

        projected_contact = args.bias
        projected_interior = restart_interior
        if outward > 0.0:
            projected0 = projected_contact
            projected1 = projected_interior
        else:
            projected0 = projected_interior
            projected1 = projected_contact

        restart_edge_drop = edge_order_drop(node0, node1, restart0, restart1)
        projected_edge_drop = edge_order_drop(node0, node1, projected0, projected1)
        final_edge_drop = edge_order_drop(node0, node1, final0, final1)
        reported_drop = edge_order_drop(
            node0,
            node1,
            optional_float(row.get("phip0")),
            optional_float(row.get("phip1")),
        )
        restart_ic_drop = (
            restart_interior - restart_contact
            if restart_interior is not None and restart_contact is not None
            else None
        )
        projected_ic_drop = (
            projected_interior - projected_contact
            if projected_interior is not None
            else None
        )
        final_ic_drop = (
            final_interior - final_contact
            if final_interior is not None and final_contact is not None
            else None
        )
        source = classify_drop(
            restart_ic_drop,
            projected_ic_drop,
            final_ic_drop,
            reported_drop,
            args.zero_atol,
        )
        drop_sources[source] += 1
        restart_drops.append(restart_ic_drop)
        projected_drops.append(projected_ic_drop)
        final_drops.append(final_ic_drop)
        out_rows.append({
            "bias_V": args.bias,
            "contact": args.contact,
            "edge_id": row.get("edge_id"),
            "node0": node0,
            "node1": node1,
            "contact_node": contact_node,
            "interior_node": interior_node,
            "outward_sign": outward,
            "restart_contact_phip_V": restart_contact,
            "restart_interior_phip_V": restart_interior,
            "restart_edge_order_phip_drop_V": restart_edge_drop,
            "restart_interior_minus_contact_phip_V": restart_ic_drop,
            "projected_contact_phip_V": projected_contact,
            "projected_interior_phip_V": projected_interior,
            "projected_edge_order_phip_drop_V": projected_edge_drop,
            "projected_interior_minus_contact_phip_V": projected_ic_drop,
            "final_contact_phip_V": final_contact,
            "final_interior_phip_V": final_interior,
            "final_edge_order_phip_drop_V": final_edge_drop,
            "final_interior_minus_contact_phip_V": final_ic_drop,
            "reported_edge_order_phip_drop_V": reported_drop,
            "drop_source": source,
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "contact_boundary_projection_edges.csv", EDGE_FIELDS, out_rows)
    dominant = drop_sources.most_common(1)[0][0] if drop_sources else "none"
    summary = {
        "bias_V": args.bias,
        "contact": args.contact,
        "edge_count": len(out_rows),
        "restart_nonzero_drop_edges": sum(
            1 for value in restart_drops if value is not None and abs(value) > args.zero_atol
        ),
        "projection_erased_restart_drop_edges": drop_sources.get(
            "contact_projection_erased_restart_drop", 0
        ),
        "solver_erased_projected_drop_edges": drop_sources.get(
            "solver_convergence_erased_projected_drop", 0
        ),
        "final_nonzero_drop_edges": sum(
            1 for value in final_drops if value is not None and abs(value) > args.zero_atol
        ),
        "max_abs_restart_interior_minus_contact_phip_V": max_abs(restart_drops),
        "max_abs_projected_interior_minus_contact_phip_V": max_abs(projected_drops),
        "max_abs_final_interior_minus_contact_phip_V": max_abs(final_drops),
        "drop_source_counts": dict(drop_sources),
        "dominant_drop_source": dominant,
    }
    (args.out_dir / "contact_boundary_projection_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
