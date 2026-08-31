#!/usr/bin/env python3
"""Audit a region-local Si/oxide interface-pair box coupling.

This is a read-only fixed-state probe.  It expands every shared material
interface node into a conceptual semiconductor master and oxide potential-only
slave, recomputes each edge's box contribution per region, and rescales only
the carrier-continuity edge flux by the semiconductor-side contribution.
Poisson continuity is not altered and no production solver option is created.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def edge_key(node0: int, node1: int) -> tuple[int, int]:
    return (node0, node1) if node0 < node1 else (node1, node0)


def triangle_area(
    a: tuple[float, float], b: tuple[float, float], c: tuple[float, float],
) -> float:
    return 0.5 * abs((b[0] - a[0]) * (c[1] - a[1])
                     - (c[0] - a[0]) * (b[1] - a[1]))


def cotangent_at_opposite(
    a: tuple[float, float], b: tuple[float, float], opposite: tuple[float, float],
) -> float:
    ux, uy = a[0] - opposite[0], a[1] - opposite[1]
    vx, vy = b[0] - opposite[0], b[1] - opposite[1]
    cross = ux * vy - uy * vx
    if abs(cross) < 1.0e-30:
        return 0.0
    return (ux * vx + uy * vy) / abs(cross)


def local_couple(
    a: tuple[float, float], b: tuple[float, float], opposite: tuple[float, float],
    *, fallback_negative: bool,
) -> tuple[float, bool]:
    length = math.hypot(b[0] - a[0], b[1] - a[1])
    if length < 1.0e-30:
        return 0.0, False
    cotangent = cotangent_at_opposite(a, b, opposite)
    if cotangent < 0.0:
        if not fallback_negative:
            return 0.0, True
        area = triangle_area(a, b, opposite)
        return area / (3.0 * length), True
    return 0.5 * cotangent * length, False


def region_local_geometry(
    mesh: dict[str, Any], transport_materials: set[str],
    *, fallback_negative: bool = True,
) -> tuple[dict[tuple[int, int], dict[str, Any]], list[dict[str, Any]]]:
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    region_by_id = {int(region["id"]): region for region in mesh["regions"]}
    edge_records: dict[tuple[int, int], dict[str, Any]] = defaultdict(
        lambda: {
            "all_couple": 0.0,
            "transport_couple": 0.0,
            "transport_cells": [],
            "insulator_cells": [],
            "negative_local_cells": [],
        }
    )
    region_nodes: dict[int, set[int]] = defaultdict(set)

    for cell in mesh["triangles"]:
        ids = [int(value) for value in cell["node_ids"]]
        region_id = int(cell["region_id"])
        region = region_by_id[region_id]
        material = str(region["material"]).lower()
        transport = material in transport_materials
        region_nodes[region_id].update(ids)
        for index in range(3):
            node0, node1 = ids[index], ids[(index + 1) % 3]
            opposite = ids[(index + 2) % 3]
            contribution, negative = local_couple(
                nodes[node0], nodes[node1], nodes[opposite],
                fallback_negative=fallback_negative,
            )
            record = edge_records[edge_key(node0, node1)]
            record["all_couple"] += contribution
            if transport:
                record["transport_couple"] += contribution
                record["transport_cells"].append(int(cell["id"]))
            else:
                record["insulator_cells"].append(int(cell["id"]))
            if negative:
                record["negative_local_cells"].append(int(cell["id"]))

    pair_records: list[dict[str, Any]] = []
    next_slave = len(nodes)
    for transport_region_id, transport_nodes in sorted(region_nodes.items()):
        transport_region = region_by_id[transport_region_id]
        if str(transport_region["material"]).lower() not in transport_materials:
            continue
        for insulator_region_id, insulator_nodes in sorted(region_nodes.items()):
            insulator_region = region_by_id[insulator_region_id]
            if str(insulator_region["material"]).lower() in transport_materials:
                continue
            for node in sorted(transport_nodes & insulator_nodes):
                pair_records.append({
                    "global_node_id": node,
                    "semiconductor_master_id": node,
                    "oxide_slave_id": next_slave,
                    "x": nodes[node][0],
                    "y": nodes[node][1],
                    "semiconductor_region": transport_region["name"],
                    "oxide_region": insulator_region["name"],
                    "master_dofs": "psi,phin,phip",
                    "slave_dofs": "psi",
                    "constraint": "psi_slave-psi_master=0",
                })
                next_slave += 1
    return dict(edge_records), pair_records


def norm_summary(values: Iterable[float]) -> dict[str, float]:
    data = list(values)
    return {
        "l1": sum(abs(value) for value in data),
        "l2": math.sqrt(sum(value * value for value in data)),
        "maximum_abs": max((abs(value) for value in data), default=0.0),
    }


def audit(
    mesh: dict[str, Any], sg_rows: list[dict[str, str]],
    carrier_rows: list[dict[str, str]], transport_materials: set[str],
    hotspot_nodes: set[int], *, coordinate_scale_m: float = 1.0e-6,
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]],
           list[dict[str, Any]]]:
    geometry, pairs = region_local_geometry(mesh, transport_materials)
    pair_nodes = {int(pair["global_node_id"]) for pair in pairs}
    original_flux = defaultdict(float)
    candidate_flux = defaultdict(float)
    original_abs = defaultdict(float)
    candidate_abs = defaultdict(float)
    edge_output: list[dict[str, Any]] = []
    global_couple_errors: list[float] = []

    for raw in sg_rows:
        edge_id = int(raw["edge_id"])
        node0, node1 = int(raw["node0"]), int(raw["node1"])
        record = geometry[edge_key(node0, node1)]
        all_couple = float(record["all_couple"])
        transport_couple = float(record["transport_couple"])
        ratio = transport_couple / all_couple if all_couple > 0.0 else 0.0
        flux = float(raw["electron_flux"])
        adjusted = flux * ratio
        original_flux[node0] += flux
        original_flux[node1] -= flux
        candidate_flux[node0] += adjusted
        candidate_flux[node1] -= adjusted
        original_abs[node0] += abs(flux)
        original_abs[node1] += abs(flux)
        candidate_abs[node0] += abs(adjusted)
        candidate_abs[node1] += abs(adjusted)
        reported_couple = float(raw["couple_m"])
        recomputed_m = all_couple * coordinate_scale_m
        scale = max(abs(reported_couple), abs(recomputed_m), 1.0e-300)
        global_couple_errors.append(abs(reported_couple - recomputed_m) / scale)
        if record["transport_cells"] and record["insulator_cells"]:
            edge_output.append({
                "edge_id": edge_id,
                "node0": node0,
                "node1": node1,
                "all_couple_m": recomputed_m,
                "transport_couple_m": transport_couple * coordinate_scale_m,
                "transport_fraction": ratio,
                "original_electron_flux": flux,
                "candidate_electron_flux": adjusted,
                "electron_flux_change": adjusted - flux,
                "transport_cells": ";".join(map(str, record["transport_cells"])),
                "insulator_cells": ";".join(map(str, record["insulator_cells"])),
                "negative_local_cells": ";".join(
                    map(str, record["negative_local_cells"])
                ),
            })

    node_output: list[dict[str, Any]] = []
    flux_reassembly_errors: list[float] = []
    flux_reassembly_abs_errors: list[float] = []
    for raw in carrier_rows:
        node = int(raw["node_id"])
        active_carrier_row = float(raw["electron_flux_abs_sum"]) > 0.0
        reported_flux = float(raw["electron_flux"])
        reassembled_flux = original_flux[node]
        if active_carrier_row:
            flux_scale = max(abs(reported_flux), abs(reassembled_flux), 1.0e-300)
            flux_reassembly_abs_errors.append(
                abs(reported_flux - reassembled_flux)
            )
            flux_reassembly_errors.append(
                abs(reported_flux - reassembled_flux) / flux_scale
            )
        nonflux = float(raw["electron_term_sum"]) - float(raw["electron_flux"])
        adjusted_residual = (
            candidate_flux[node] + nonflux
            if active_carrier_row else float(raw["electron_residual"])
        )
        node_output.append({
            "node_id": node,
            "x": float(raw["x"]),
            "y": float(raw["y"]),
            "is_interface_pair": node in pair_nodes,
            "is_frozen_hotspot": node in hotspot_nodes,
            "is_active_carrier_row": active_carrier_row,
            "reported_electron_flux": reported_flux,
            "reassembled_electron_flux": reassembled_flux,
            "candidate_electron_flux": candidate_flux[node],
            "reported_electron_residual": float(raw["electron_residual"]),
            "candidate_electron_residual": adjusted_residual,
            "original_flux_abs_sum": original_abs[node],
            "candidate_flux_abs_sum": candidate_abs[node],
        })

    baseline = norm_summary(
        row["reported_electron_residual"]
        for row in node_output if row["is_active_carrier_row"]
    )
    candidate = norm_summary(
        row["candidate_electron_residual"]
        for row in node_output if row["is_active_carrier_row"]
    )
    ratios = {
        key: candidate[key] / max(baseline[key], 1.0e-300)
        for key in baseline
    }
    ranked = sorted(
        node_output,
        key=lambda row: abs(row["candidate_electron_residual"]), reverse=True,
    )
    hotspot_detail = [
        row for row in node_output if row["node_id"] in hotspot_nodes
    ]
    hotspot_improved = all(
        abs(row["candidate_electron_residual"])
        <= 0.5 * max(abs(row["reported_electron_residual"]), 1.0e-300)
        for row in hotspot_detail
        if row["is_interface_pair"] and row["is_active_carrier_row"]
    )
    gate_pass = (
        ratios["l2"] <= 0.5
        and ratios["maximum_abs"] <= 0.5
        and hotspot_improved
    )
    summary = {
        "schema": "vela.templates_ldmos.interface_pair_box_audit.v1",
        "contract": {
            "mode": "read_only_fixed_state",
            "semiconductor_master_dofs": ["psi", "phin", "phip"],
            "oxide_slave_dofs": ["psi"],
            "potential_constraint": "psi_slave-psi_master=0",
            "carrier_coupling": "transport_region_local_box_contribution",
            "poisson_coupling": "unchanged_continuous_interface",
            "production_solver_modified": False,
        },
        "topology": {
            "global_nodes": len(mesh["nodes"]),
            "conceptual_interface_pairs": len(pairs),
            "conceptual_total_nodes": len(mesh["nodes"]) + len(pairs),
            "interface_transport_edges": len(edge_output),
        },
        "geometry_reproduction": {
            "maximum_relative_couple_error": max(global_couple_errors, default=0.0),
            "maximum_relative_node_flux_reassembly_error": max(
                flux_reassembly_errors, default=0.0
            ),
            "maximum_absolute_node_flux_reassembly_error": max(
                flux_reassembly_abs_errors, default=0.0
            ),
        },
        "baseline_electron_residual": baseline,
        "candidate_electron_residual": candidate,
        "candidate_over_baseline": ratios,
        "frozen_hotspots": hotspot_detail,
        "top_candidate_rows": ranked[:12],
        "same_bias_reclose_gate": {
            "requires_l2_ratio_at_most": 0.5,
            "requires_maximum_ratio_at_most": 0.5,
            "requires_each_interface_hotspot_ratio_at_most": 0.5,
            "passed": gate_pass,
        },
    }
    edge_output.sort(key=lambda row: abs(row["electron_flux_change"]), reverse=True)
    pairs.sort(key=lambda row: (row["global_node_id"], row["oxide_region"]))
    return summary, pairs, edge_output, node_output


def write_records(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("\n", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sg-edges", type=Path, required=True)
    parser.add_argument("--carrier-terms", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--transport-materials", default="Si,Silicon,PolySi")
    parser.add_argument("--hotspot-nodes", default="3747,10233,4492,4538")
    args = parser.parse_args()

    materials = {
        value.strip().lower()
        for value in args.transport_materials.split(",") if value.strip()
    }
    hotspots = {
        int(value) for value in args.hotspot_nodes.split(",") if value.strip()
    }
    summary, pairs, edges, nodes = audit(
        json.loads(args.mesh.read_text(encoding="utf-8")),
        read_csv(args.sg_edges), read_csv(args.carrier_terms), materials, hotspots,
    )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_records(output / "interface_pairs.csv", pairs)
    write_records(output / "interface_edges.csv", edges)
    write_records(output / "node_residuals.csv", nodes)
    summary["artifacts"] = {
        "interface_pairs": str((output / "interface_pairs.csv").resolve()),
        "interface_edges": str((output / "interface_edges.csv").resolve()),
        "node_residuals": str((output / "node_residuals.csv").resolve()),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str((output / "summary.json").resolve()),
        "candidate_over_baseline": summary["candidate_over_baseline"],
        "same_bias_reclose_gate": summary["same_bias_reclose_gate"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
