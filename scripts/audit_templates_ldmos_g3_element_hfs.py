#!/usr/bin/env python3
"""Audit Sentaurus G3 element mobility against candidate HFS drives.

The diagnostic compares the saved element eMobility with Caughey--Thomas
predictions driven by element GradQuasiFermi and ElectricField.  It reports
the full silicon domain, the drain-node one-ring, and the exact drain-boundary
cell mask separately.  Saved TDR mobility remains diagnostic evidence rather
than a Vela assembly input.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Iterable


BIAS_POINTS_V = tuple(index / 6.0 for index in range(7))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(value for value in values if math.isfinite(value))
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def stats(values: Iterable[float]) -> dict[str, float | int]:
    samples = [value for value in values if math.isfinite(value)]
    return {
        "count": len(samples),
        "minimum": min(samples, default=math.nan),
        "median": statistics.median(samples) if samples else math.nan,
        "p95": percentile(samples, 0.95),
        "maximum": max(samples, default=math.nan),
    }


def scalar_cells(path: Path) -> dict[int, float]:
    return {
        int(row["cell_id"]): float(row["component0"])
        for row in read_csv(path)
    }


def vector_cells(path: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["cell_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_csv(path)
    }


def scalar_nodes(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(path)
    }


def magnitude(vector: tuple[float, float]) -> float:
    return math.hypot(*vector)


def caughey_thomas(
    low_field_mobility_cm2_V_s: float, drive_V_cm: float,
    saturation_velocity_cm_s: float, beta: float,
) -> float:
    normalized = low_field_mobility_cm2_V_s * abs(drive_V_cm) / saturation_velocity_cm_s
    return low_field_mobility_cm2_V_s / (
        1.0 + normalized ** beta
    ) ** (1.0 / beta)


def log_error(predicted: float, observed: float) -> float:
    return abs(math.log10(max(predicted, 1.0e-300) / max(observed, 1.0e-300)))


def cell_masks(mesh: dict) -> dict[str, set[int]]:
    silicon = {
        int(region["id"]) for region in mesh["regions"]
        if region["material"].lower() == "si"
    }
    drain = next(
        contact for contact in mesh["contacts"]
        if contact["name"].lower() == "drain"
    )
    drain_nodes = {int(node) for node in drain["node_ids"]}
    drain_edges = {
        frozenset(int(node) for node in edge)
        for edge in drain["edge_node_ids"]
    }
    silicon_cells = {
        int(cell["id"]) for cell in mesh["triangles"]
        if int(cell["region_id"]) in silicon
    }
    one_ring = {
        int(cell["id"]) for cell in mesh["triangles"]
        if int(cell["id"]) in silicon_cells
        and drain_nodes.intersection(int(node) for node in cell["node_ids"])
    }
    boundary = {
        int(cell["id"]) for cell in mesh["triangles"]
        if int(cell["id"]) in silicon_cells
        and any(
            edge.issubset(int(node) for node in cell["node_ids"])
            for edge in drain_edges
        )
    }
    return {
        "silicon_all": silicon_cells,
        "drain_node_one_ring": one_ring,
        "drain_boundary_cells": boundary,
    }


def audit_point(
    export_dir: Path, mesh: dict, bias: float,
    low_field: float, vsat: float, beta: float, output: Path,
) -> dict:
    fields = export_dir / "fields"
    mobility = scalar_cells(fields / "eMobility_region0_cells.csv")
    grad_qf = vector_cells(fields / "eGradQuasiFermi_region0_cells.csv")
    electric = vector_cells(fields / "ElectricField_region0_cells.csv")
    eparallel = scalar_nodes(fields / "eEparallel_region0.csv")
    masks = cell_masks(mesh)
    triangles = {int(cell["id"]): cell for cell in mesh["triangles"]}
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    common = sorted(set(mobility) & set(grad_qf) & set(electric))
    records = []
    for cell in common:
        ids = [int(node) for node in triangles[cell]["node_ids"]]
        qf = magnitude(grad_qf[cell])
        field = magnitude(electric[cell])
        nodal_parallel = statistics.mean(abs(eparallel[node]) for node in ids)
        qf_prediction = caughey_thomas(low_field, qf, vsat, beta)
        field_prediction = caughey_thomas(low_field, field, vsat, beta)
        parallel_prediction = caughey_thomas(
            low_field, nodal_parallel, vsat, beta
        )
        centroid_x = statistics.mean(nodes[node][0] for node in ids)
        centroid_y = statistics.mean(nodes[node][1] for node in ids)
        records.append({
            "bias_V": bias,
            "cell_id": cell,
            "centroid_x_m": centroid_x,
            "centroid_y_m": centroid_y,
            "drain_node_one_ring": int(cell in masks["drain_node_one_ring"]),
            "drain_boundary_cell": int(cell in masks["drain_boundary_cells"]),
            "observed_eMobility_cm2_V_s": mobility[cell],
            "element_grad_qf_V_cm": qf,
            "element_electric_field_V_cm": field,
            "mean_nodal_eEparallel_V_cm": nodal_parallel,
            "qf_prediction_cm2_V_s": qf_prediction,
            "electric_field_prediction_cm2_V_s": field_prediction,
            "nodal_eEparallel_prediction_cm2_V_s": parallel_prediction,
            "qf_error_dex": log_error(qf_prediction, mobility[cell]),
            "electric_field_error_dex": log_error(field_prediction, mobility[cell]),
            "nodal_eEparallel_error_dex": log_error(
                parallel_prediction, mobility[cell]
            ),
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)

    result: dict[str, dict] = {}
    by_cell = {int(row["cell_id"]): row for row in records}
    for name, cells in masks.items():
        selected = [by_cell[cell] for cell in sorted(cells & set(by_cell))]
        result[name] = {
            "cells": len(selected),
            "observed_eMobility_cm2_V_s": stats(
                float(row["observed_eMobility_cm2_V_s"]) for row in selected
            ),
            "element_grad_qf_V_cm": stats(
                float(row["element_grad_qf_V_cm"]) for row in selected
            ),
            "element_electric_field_V_cm": stats(
                float(row["element_electric_field_V_cm"]) for row in selected
            ),
            "mean_nodal_eEparallel_V_cm": stats(
                float(row["mean_nodal_eEparallel_V_cm"]) for row in selected
            ),
            "qf_prediction_error_dex": stats(
                float(row["qf_error_dex"]) for row in selected
            ),
            "electric_field_prediction_error_dex": stats(
                float(row["electric_field_error_dex"]) for row in selected
            ),
            "nodal_eEparallel_prediction_error_dex": stats(
                float(row["nodal_eEparallel_error_dex"]) for row in selected
            ),
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--export-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--low-field-mobility", type=float, default=1417.0)
    parser.add_argument("--saturation-velocity", type=float, default=1.07e7)
    parser.add_argument("--beta", type=float, default=1.109)
    args = parser.parse_args()

    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    summary = {
        "schema": "vela.templates_ldmos_g3_element_hfs_audit.v1",
        "formula": {
            "name": "Caughey-Thomas",
            "low_field_mobility_cm2_V_s": args.low_field_mobility,
            "saturation_velocity_cm_s": args.saturation_velocity,
            "beta": args.beta,
        },
        "plot_contract": {
            "eMobility_element": "available",
            "eGradQuasiFermi_element_vector": "available",
            "ElectricField_element_vector": "available",
            "eEparallel_element": "unavailable in T-2022.03; nodal mean diagnostic used",
        },
        "points": [],
    }
    for index in (1, 3, 5):
        bias = BIAS_POINTS_V[index]
        point = audit_point(
            args.export_root / f"export_{index:04d}", mesh, bias,
            args.low_field_mobility, args.saturation_velocity, args.beta,
            output / f"vg_{bias:.6f}_cells.csv".replace(".", "p", 1),
        )
        summary["points"].append({"bias_V": bias, "masks": point})
    contact_qf = [
        point["masks"]["drain_node_one_ring"]["qf_prediction_error_dex"]["median"]
        for point in summary["points"]
    ]
    contact_field = [
        point["masks"]["drain_node_one_ring"]
        ["electric_field_prediction_error_dex"]["median"]
        for point in summary["points"]
    ]
    summary["decision"] = {
        "contact_qf_median_error_dex": statistics.median(contact_qf),
        "contact_electric_field_median_error_dex": statistics.median(contact_field),
        "contact_hfs_support": (
            "electric_field_or_boundary_fallback_supported"
            if statistics.median(contact_field) < statistics.median(contact_qf)
            else "element_grad_qf_supported"
        ),
        "caveat": (
            "Saved mobility diagnoses Sentaurus support but is not imported as "
            "Vela assembly truth; G4 fixed-state replay supplies the assembly-real control."
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8",
    )
    print(json.dumps(summary["decision"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
