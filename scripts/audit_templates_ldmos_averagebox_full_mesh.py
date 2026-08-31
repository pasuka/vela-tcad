#!/usr/bin/env python3
"""Replay direct Sentaurus AverageBox data over every LDMOS silicon edge.

This is a read-only fixed-state audit.  It maps the complete grid-numbered
MeasureCoefficients.debug file to the exact imported topology, aggregates
elementwise coefficients and measures only over transport-material cells, and
reassembles every free electron-continuity row without changing Vela's mesh or
production solver policy.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.audit_templates_ldmos_averagebox_node4492 import (
    coefficient_edge,
    edge_key,
    norm,
    parse_debug_block,
    read_csv,
    write_csv,
)
from scripts.audit_templates_ldmos_interface_pair_box import region_local_geometry


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def triangle_area(points: list[tuple[float, float]]) -> float:
    a, b, c = points
    return 0.5 * abs(
        (b[0] - a[0]) * (c[1] - a[1])
        - (c[0] - a[0]) * (b[1] - a[1])
    )


def averagebox_geometry(
    mesh: dict[str, Any],
    measures: dict[int, list[float]],
    coefficients: dict[int, list[float]],
    transport_materials: set[str],
) -> tuple[dict[tuple[int, int], float], dict[int, float], dict[int, float],
           list[dict[str, Any]], dict[str, Any]]:
    points = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    regions = {int(region["id"]): region for region in mesh["regions"]}
    candidate_couple_um: dict[tuple[int, int], float] = defaultdict(float)
    candidate_measure_um2: dict[int, float] = defaultdict(float)
    vela_global_volume_um2: dict[int, float] = defaultdict(float)
    cell_rows: list[dict[str, Any]] = []
    transport_cells = 0
    missing: list[int] = []

    for cell in mesh["triangles"]:
        cell_id = int(cell["id"])
        nodes = [int(value) for value in cell["node_ids"]]
        local_points = [points[node] for node in nodes]
        area = triangle_area(local_points)
        for node in nodes:
            vela_global_volume_um2[node] += area / 3.0
        region = regions[int(cell["region_id"])]
        is_transport = str(region["material"]).lower() in transport_materials
        if not is_transport:
            continue
        transport_cells += 1
        if cell_id not in measures or cell_id not in coefficients:
            missing.append(cell_id)
            continue
        if len(measures[cell_id]) != 3 or len(coefficients[cell_id]) != 3:
            raise ValueError(f"transport cell {cell_id} is not a Tri3 debug record")
        for local_node, node in enumerate(nodes):
            candidate_measure_um2[node] += measures[cell_id][local_node]
        for coefficient_index in range(3):
            node0, node1, opposite = coefficient_edge(nodes, coefficient_index)
            length_um = math.dist(points[node0], points[node1])
            local_couple_um = coefficients[cell_id][coefficient_index] * length_um
            candidate_couple_um[edge_key(node0, node1)] += local_couple_um
            cell_rows.append({
                "cell_id": cell_id,
                "region_id": int(cell["region_id"]),
                "coefficient_index": coefficient_index,
                "node0": node0,
                "node1": node1,
                "opposite_local_vertex": opposite,
                "averagebox_coefficient": coefficients[cell_id][coefficient_index],
                "edge_length_um": length_um,
                "averagebox_local_couple_um": local_couple_um,
            })
    if missing:
        raise ValueError(f"debug oracle lacks transport cells: {missing[:20]}")
    metadata = {
        "transport_cells": transport_cells,
        "transport_edges": len(candidate_couple_um),
        "transport_nodes": len(candidate_measure_um2),
    }
    return (
        dict(candidate_couple_um), dict(candidate_measure_um2),
        dict(vela_global_volume_um2), cell_rows, metadata,
    )


def ratio(candidate: dict[str, float], baseline: dict[str, float]) -> dict[str, float]:
    return {
        key: candidate[key] / max(baseline[key], 1.0e-300)
        for key in baseline
    }


def audit(
    mesh: dict[str, Any],
    sg_rows: list[dict[str, str]],
    carrier_rows: list[dict[str, str]],
    measures: dict[int, list[float]],
    coefficients: dict[int, list[float]],
    transport_materials: set[str],
    frozen_hotspots: set[int],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]],
           list[dict[str, Any]]]:
    (
        candidate_couple_um, candidate_measure_um2, vela_volume_um2,
        cell_rows, geometry_metadata,
    ) = averagebox_geometry(
        mesh, measures, coefficients, transport_materials
    )
    production_geometry, interface_pairs = region_local_geometry(
        mesh, transport_materials
    )
    interface_nodes = {int(pair["global_node_id"]) for pair in interface_pairs}

    baseline_flux: dict[int, float] = defaultdict(float)
    candidate_flux: dict[int, float] = defaultdict(float)
    candidate_abs_sum: dict[int, float] = defaultdict(float)
    edge_rows: list[dict[str, Any]] = []
    production_couple_errors: list[float] = []
    for raw in sg_rows:
        node0, node1 = int(raw["node0"]), int(raw["node1"])
        key = edge_key(node0, node1)
        original = float(raw["electron_flux"])
        sent_couple_m = candidate_couple_um.get(key, 0.0) * 1.0e-6
        candidate = (
            float(raw["electron_scaled_flux_per_couple_m"]) * sent_couple_m
        )
        baseline_flux[node0] += original
        baseline_flux[node1] -= original
        candidate_flux[node0] += candidate
        candidate_flux[node1] -= candidate
        candidate_abs_sum[node0] += abs(candidate)
        candidate_abs_sum[node1] += abs(candidate)

        reported_couple = float(raw["couple_m"])
        reproduced = float(production_geometry[key]["all_couple"]) * 1.0e-6
        production_couple_errors.append(abs(reported_couple - reproduced))
        if original != 0.0 or candidate != 0.0:
            edge_rows.append({
                "edge_id": int(raw["edge_id"]),
                "node0": node0,
                "node1": node1,
                "is_interface_edge": bool(
                    production_geometry[key]["transport_cells"]
                    and production_geometry[key]["insulator_cells"]
                ),
                "production_couple_m": reported_couple,
                "averagebox_silicon_couple_m": sent_couple_m,
                "averagebox_over_production_couple": (
                    sent_couple_m / reported_couple
                    if reported_couple else (math.inf if sent_couple_m else 1.0)
                ),
                "production_electron_flux": original,
                "averagebox_electron_flux": candidate,
                "electron_flux_change": candidate - original,
            })

    node_rows: list[dict[str, Any]] = []
    flux_reassembly_errors: list[float] = []
    baseline_values: list[float] = []
    coefficient_values: list[float] = []
    measure_values: list[float] = []
    for raw in carrier_rows:
        node = int(raw["node_id"])
        active = float(raw["electron_flux_abs_sum"]) > 0.0
        reported_flux = float(raw["electron_flux"])
        reassembled_flux = baseline_flux[node]
        if active:
            flux_reassembly_errors.append(abs(reported_flux - reassembled_flux))

        recombination = float(raw["electron_recombination"])
        impact = float(raw["electron_impact"])
        gauge = float(raw["electron_gauge"])
        boundary = float(raw["electron_boundary"])
        nonflux = recombination + impact + gauge + boundary
        coefficient_residual = (
            candidate_flux[node] + nonflux
            if active else float(raw["electron_residual"])
        )
        volume_scale = (
            candidate_measure_um2.get(node, 0.0)
            / max(vela_volume_um2.get(node, 0.0), 1.0e-300)
        )
        measure_residual = (
            candidate_flux[node]
            + volume_scale * (recombination + impact)
            + gauge + boundary
            if active else float(raw["electron_residual"])
        )
        baseline = float(raw["electron_residual"])
        if active:
            baseline_values.append(baseline)
            coefficient_values.append(coefficient_residual)
            measure_values.append(measure_residual)
        node_rows.append({
            "node_id": node,
            "x": float(raw["x"]),
            "y": float(raw["y"]),
            "is_active_carrier_row": active,
            "is_interface_node": node in interface_nodes,
            "is_frozen_hotspot": node in frozen_hotspots,
            "production_volume_um2": vela_volume_um2.get(node, 0.0),
            "averagebox_silicon_measure_um2": candidate_measure_um2.get(node, 0.0),
            "measure_over_production_volume": volume_scale,
            "reported_electron_flux": reported_flux,
            "reassembled_production_electron_flux": reassembled_flux,
            "averagebox_electron_flux": candidate_flux[node],
            "averagebox_flux_abs_sum": candidate_abs_sum[node],
            "baseline_electron_residual": baseline,
            "coefficient_only_electron_residual": coefficient_residual,
            "coefficient_measure_electron_residual": measure_residual,
        })

    baseline_norm = norm(baseline_values)
    coefficient_norm = norm(coefficient_values)
    measure_norm = norm(measure_values)
    coefficient_ratio = ratio(coefficient_norm, baseline_norm)
    measure_ratio = ratio(measure_norm, baseline_norm)
    baseline_ranked = sorted(
        (row for row in node_rows if row["is_active_carrier_row"]),
        key=lambda row: abs(row["baseline_electron_residual"]), reverse=True,
    )
    candidate_ranked = sorted(
        (row for row in node_rows if row["is_active_carrier_row"]),
        key=lambda row: abs(row["coefficient_measure_electron_residual"]), reverse=True,
    )
    baseline_top_nodes = {int(row["node_id"]) for row in baseline_ranked[:12]}
    candidate_top_nodes = {int(row["node_id"]) for row in candidate_ranked[:12]}
    new_top_nodes = sorted(candidate_top_nodes - baseline_top_nodes)
    candidate_max_over_baseline_max = (
        measure_norm["maximum_abs"]
        / max(baseline_norm["maximum_abs"], 1.0e-300)
    )
    improved_rows = sum(
        abs(row["coefficient_measure_electron_residual"])
        < abs(row["baseline_electron_residual"])
        for row in node_rows if row["is_active_carrier_row"]
    )
    improved_twofold_rows = sum(
        abs(row["coefficient_measure_electron_residual"])
        <= 0.5 * abs(row["baseline_electron_residual"])
        for row in node_rows if row["is_active_carrier_row"]
    )
    frozen_detail = [
        row for row in node_rows
        if row["node_id"] in frozen_hotspots and row["is_active_carrier_row"]
    ]
    frozen_pass = all(
        abs(row["coefficient_measure_electron_residual"])
        <= 0.5 * max(abs(row["baseline_electron_residual"]), 1.0e-300)
        for row in frozen_detail
    )
    gate_pass = (
        measure_ratio["l2"] <= 0.5
        and measure_ratio["maximum_abs"] <= 0.5
        and frozen_pass
    )
    edge_rows.sort(key=lambda row: abs(row["electron_flux_change"]), reverse=True)
    summary = {
        "schema": "vela.templates_ldmos.averagebox_full_mesh_audit.v1",
        "contract": {
            "mode": "read_only_fixed_state",
            "transport_geometry": "direct_averagebox_transport_material_cells_only",
            "source_geometry": "direct_averagebox_transport_material_measure",
            "production_solver_modified": False,
        },
        "geometry": {
            **geometry_metadata,
            "interface_nodes": len(interface_nodes),
            "maximum_absolute_production_couple_reproduction_error_m": max(
                production_couple_errors, default=0.0
            ),
            "maximum_absolute_node_flux_reassembly_error": max(
                flux_reassembly_errors, default=0.0
            ),
        },
        "baseline_electron_residual": baseline_norm,
        "coefficient_only_electron_residual": coefficient_norm,
        "coefficient_only_over_baseline": coefficient_ratio,
        "coefficient_measure_electron_residual": measure_norm,
        "coefficient_measure_over_baseline": measure_ratio,
        "frozen_hotspots": frozen_detail,
        "baseline_top_rows": baseline_ranked[:12],
        "candidate_top_rows": candidate_ranked[:12],
        "row_outcomes": {
            "active_rows": len(baseline_values),
            "improved_rows": improved_rows,
            "improved_at_least_twofold_rows": improved_twofold_rows,
            "worsened_or_equal_rows": len(baseline_values) - improved_rows,
        },
        "hotspot_transfer": {
            "baseline_top12_nodes": sorted(baseline_top_nodes),
            "candidate_top12_nodes": sorted(candidate_top_nodes),
            "new_top12_nodes": new_top_nodes,
            "ranking_changed": bool(new_top_nodes),
            "candidate_maximum_over_baseline_maximum":
                candidate_max_over_baseline_max,
            "material_transfer_threshold": 0.5,
            "material_hotspot_transfer_detected":
                candidate_max_over_baseline_max > 0.5,
        },
        "external_profile_implementation_gate": {
            "requires_l2_ratio_at_most": 0.5,
            "requires_maximum_ratio_at_most": 0.5,
            "requires_each_frozen_hotspot_ratio_at_most": 0.5,
            "frozen_hotspots_passed": frozen_pass,
            "passed": gate_pass,
        },
    }
    return summary, cell_rows, edge_rows, node_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sg-edges", type=Path, required=True)
    parser.add_argument("--carrier-terms", type=Path, required=True)
    parser.add_argument("--measure-coefficients", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--transport-materials", default="Si,Silicon")
    parser.add_argument("--hotspot-nodes", default="3747,10233,4492,4538")
    args = parser.parse_args()

    debug = args.measure_coefficients.resolve()
    summary, cells, edges, nodes = audit(
        json.loads(args.mesh.read_text(encoding="utf-8")),
        read_csv(args.sg_edges),
        read_csv(args.carrier_terms),
        parse_debug_block(debug, "Measure"),
        parse_debug_block(debug, "Coefficients"),
        {value.strip().lower() for value in args.transport_materials.split(",")},
        {int(value) for value in args.hotspot_nodes.split(",") if value.strip()},
    )
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "cell_coefficients.csv", cells)
    write_csv(output / "edge_replay.csv", edges)
    write_csv(output / "node_residuals.csv", nodes)
    summary["oracle"] = {"path": str(debug), "sha256": sha256(debug)}
    summary["artifacts"] = {
        "cell_coefficients": str((output / "cell_coefficients.csv").resolve()),
        "edge_replay": str((output / "edge_replay.csv").resolve()),
        "node_residuals": str((output / "node_residuals.csv").resolve()),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str((output / "summary.json").resolve()),
        "coefficient_only_over_baseline": summary["coefficient_only_over_baseline"],
        "coefficient_measure_over_baseline": summary["coefficient_measure_over_baseline"],
        "hotspot_transfer": summary["hotspot_transfer"],
        "gate": summary["external_profile_implementation_gate"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
