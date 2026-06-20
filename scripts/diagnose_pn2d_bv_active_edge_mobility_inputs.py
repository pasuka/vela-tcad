#!/usr/bin/env python3
"""Decompose active-edge mobility inputs for PN2D BV support nodes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "active_axis",
    "active_edge_count",
    "active_edge_ids",
    "electron_low_field_mobility_m2_V_s",
    "electron_mobility_limiter",
    "electron_limited_mobility_product_m2_V_s",
    "electron_final_mobility_m2_V_s",
    "sentaurus_electron_mobility_m2_V_s",
    "electron_product_over_final_mobility",
    "electron_final_over_low_field_mobility",
    "electron_final_over_sentaurus_mobility",
    "hole_low_field_mobility_m2_V_s",
    "hole_mobility_limiter",
    "hole_limited_mobility_product_m2_V_s",
    "hole_final_mobility_m2_V_s",
    "sentaurus_hole_mobility_m2_V_s",
    "hole_product_over_final_mobility",
    "hole_final_over_low_field_mobility",
    "hole_final_over_sentaurus_mobility",
    "vela_electric_field_scalar_avg_V_cm",
    "vela_electric_field_scalar_avg_V_m",
    "vela_electric_field_edge_abs_V_m",
    "sentaurus_electric_field_edge_abs_V_m",
    "electric_field_edge_abs_vela_over_sentaurus",
    "vela_electron_high_field_drive_avg_V_cm",
    "vela_electron_high_field_drive_avg_V_m",
    "sentaurus_electron_qf_field_edge_abs_V_m",
    "electron_drive_vela_over_sentaurus_qf",
    "vela_hole_high_field_drive_avg_V_cm",
    "vela_hole_high_field_drive_avg_V_m",
    "sentaurus_hole_qf_field_edge_abs_V_m",
    "hole_drive_vela_over_sentaurus_qf",
]


SUMMARY_METRICS = [
    "electron_low_field_mobility_m2_V_s",
    "electron_mobility_limiter",
    "electron_product_over_final_mobility",
    "electron_final_over_low_field_mobility",
    "electron_final_over_sentaurus_mobility",
    "hole_low_field_mobility_m2_V_s",
    "hole_mobility_limiter",
    "hole_product_over_final_mobility",
    "hole_final_over_low_field_mobility",
    "hole_final_over_sentaurus_mobility",
    "electric_field_edge_abs_vela_over_sentaurus",
    "electron_drive_vela_over_sentaurus_qf",
    "hole_drive_vela_over_sentaurus_qf",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--active-axis", choices=["x", "y"], default="x")
    parser.add_argument("--active-source-relative-threshold", type=float, default=1.0e-12)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def endpoint_avg(values: list[float], edge: dict[str, Any]) -> float:
    i = int(edge["node0"])
    j = int(edge["node1"])
    return 0.5 * (values[i] + values[j])


def edge_abs_field(values: list[float], edge: dict[str, Any]) -> float:
    h = float(edge["length_m"])
    if h <= 0.0:
        return 0.0
    i = int(edge["node0"])
    j = int(edge["node1"])
    return abs(values[j] - values[i]) / h


def active_edge_ids(active: list[dict[str, Any]]) -> str:
    return ";".join(str(int(item["edge_id"])) for item in active)


def area_weighted(active: list[dict[str, Any]], values: dict[int, float]) -> float | None:
    numerator = 0.0
    denominator = 0.0
    for item in active:
        edge_id = int(item["edge_id"])
        if edge_id not in values:
            continue
        area = float(item["endpoint_area_m2"])
        numerator += values[edge_id] * area
        denominator += area
    if denominator == 0.0:
        return None
    return numerator / denominator


def support_rows(path: Path) -> list[dict[str, str]]:
    return [
        row for row in read_csv_rows(path)
        if row.get("support_class") not in (None, "", "inactive")
    ]


def scalar_or_default(scalars: dict[str, list[float]], name: str, node_count: int, default: float = 0.0) -> list[float]:
    values = scalars.get(name)
    if values is None:
        return [default for _ in range(node_count)]
    if len(values) != node_count:
        raise RuntimeError(f"{name} has {len(values)} values, expected {node_count}")
    return values


def edge_maps(
    edges: dict[int, dict[str, Any]],
    edge_ids: list[int],
    vela_scalars: dict[str, list[float]],
    vela_state: dict[str, list[float]],
    sentaurus_state: dict[str, list[float]],
    node_count: int,
) -> dict[str, dict[int, float]]:
    electron_low = scalar_or_default(vela_scalars, "ElectronLowFieldMobility", node_count)
    hole_low = scalar_or_default(vela_scalars, "HoleLowFieldMobility", node_count)
    electron_limiter = scalar_or_default(vela_scalars, "ElectronMobilityLimiter", node_count, 1.0)
    hole_limiter = scalar_or_default(vela_scalars, "HoleMobilityLimiter", node_count, 1.0)
    electric_scalar = scalar_or_default(vela_scalars, "ElectricField", node_count)
    electron_drive = scalar_or_default(vela_scalars, "ElectronHighFieldDrive", node_count)
    hole_drive = scalar_or_default(vela_scalars, "HoleHighFieldDrive", node_count)
    maps: dict[str, dict[int, float]] = {key: {} for key in [
        "electron_low",
        "electron_limiter",
        "electron_product",
        "electron_final",
        "sentaurus_electron",
        "hole_low",
        "hole_limiter",
        "hole_product",
        "hole_final",
        "sentaurus_hole",
        "vela_electric_scalar",
        "vela_electric_scalar_si",
        "vela_electric_edge",
        "sentaurus_electric_edge",
        "electron_drive",
        "electron_drive_si",
        "sentaurus_electron_qf_edge",
        "hole_drive",
        "hole_drive_si",
        "sentaurus_hole_qf_edge",
    ]}
    for edge_id in edge_ids:
        edge = edges[edge_id]
        e_low = endpoint_avg(electron_low, edge)
        e_limiter = endpoint_avg(electron_limiter, edge)
        h_low = endpoint_avg(hole_low, edge)
        h_limiter = endpoint_avg(hole_limiter, edge)
        maps["electron_low"][edge_id] = e_low
        maps["electron_limiter"][edge_id] = e_limiter
        maps["electron_product"][edge_id] = e_low * e_limiter
        maps["electron_final"][edge_id] = endpoint_avg(vela_state["mun"], edge)
        maps["sentaurus_electron"][edge_id] = endpoint_avg(sentaurus_state["mun"], edge)
        maps["hole_low"][edge_id] = h_low
        maps["hole_limiter"][edge_id] = h_limiter
        maps["hole_product"][edge_id] = h_low * h_limiter
        maps["hole_final"][edge_id] = endpoint_avg(vela_state["mup"], edge)
        maps["sentaurus_hole"][edge_id] = endpoint_avg(sentaurus_state["mup"], edge)
        maps["vela_electric_scalar"][edge_id] = endpoint_avg(electric_scalar, edge)
        maps["vela_electric_scalar_si"][edge_id] = endpoint_avg(electric_scalar, edge) * 100.0
        maps["vela_electric_edge"][edge_id] = edge_abs_field(vela_state["psi"], edge)
        maps["sentaurus_electric_edge"][edge_id] = edge_abs_field(sentaurus_state["psi"], edge)
        maps["electron_drive"][edge_id] = endpoint_avg(electron_drive, edge)
        maps["electron_drive_si"][edge_id] = endpoint_avg(electron_drive, edge) * 100.0
        maps["sentaurus_electron_qf_edge"][edge_id] = edge_abs_field(sentaurus_state["phin"], edge)
        maps["hole_drive"][edge_id] = endpoint_avg(hole_drive, edge)
        maps["hole_drive_si"][edge_id] = endpoint_avg(hole_drive, edge) * 100.0
        maps["sentaurus_hole_qf_edge"][edge_id] = edge_abs_field(sentaurus_state["phip"], edge)
    return maps


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, _ = fluxforms.sgdiag.read_mesh(args.mesh)
    edges = {int(edge["edge_id"]): edge for edge in fluxforms.sgdiag.build_edges(nodes, triangles)}
    node_count = len(nodes)
    vtk = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    vela_scalars = fluxforms.sgdiag.parse_vtk_scalars(vtk)
    vela_state = fluxforms.load_vela_state(vtk, node_count)
    sentaurus_state = fluxforms.load_sentaurus_state(args.sentaurus_dir, node_count)
    active_by_node = fluxfac.load_active_cxx_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    rows_in = support_rows(args.support_csv)
    needed_edge_ids = sorted({
        int(item["edge_id"])
        for support in rows_in
        for item in active_by_node.get(int(support["node_id"]), [])
    })
    missing = [edge_id for edge_id in needed_edge_ids if edge_id not in edges]
    if missing:
        raise RuntimeError(f"active edge ids not found in mesh: {missing[:8]}")
    maps = edge_maps(edges, needed_edge_ids, vela_scalars, vela_state, sentaurus_state, node_count)

    rows: list[dict[str, Any]] = []
    for support in rows_in:
        node_id = int(support["node_id"])
        active = active_by_node.get(node_id, [])
        electron_low = area_weighted(active, maps["electron_low"])
        electron_limiter = area_weighted(active, maps["electron_limiter"])
        electron_product = area_weighted(active, maps["electron_product"])
        electron_final = area_weighted(active, maps["electron_final"])
        sentaurus_electron = area_weighted(active, maps["sentaurus_electron"])
        hole_low = area_weighted(active, maps["hole_low"])
        hole_limiter = area_weighted(active, maps["hole_limiter"])
        hole_product = area_weighted(active, maps["hole_product"])
        hole_final = area_weighted(active, maps["hole_final"])
        sentaurus_hole = area_weighted(active, maps["sentaurus_hole"])
        vela_electric_edge = area_weighted(active, maps["vela_electric_edge"])
        sentaurus_electric_edge = area_weighted(active, maps["sentaurus_electric_edge"])
        electron_drive = area_weighted(active, maps["electron_drive"])
        electron_drive_si = area_weighted(active, maps["electron_drive_si"])
        sentaurus_electron_qf = area_weighted(active, maps["sentaurus_electron_qf_edge"])
        hole_drive = area_weighted(active, maps["hole_drive"])
        hole_drive_si = area_weighted(active, maps["hole_drive_si"])
        sentaurus_hole_qf = area_weighted(active, maps["sentaurus_hole_qf_edge"])
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "active_axis": args.active_axis,
            "active_edge_count": len(active),
            "active_edge_ids": active_edge_ids(active),
            "electron_low_field_mobility_m2_V_s": electron_low,
            "electron_mobility_limiter": electron_limiter,
            "electron_limited_mobility_product_m2_V_s": electron_product,
            "electron_final_mobility_m2_V_s": electron_final,
            "sentaurus_electron_mobility_m2_V_s": sentaurus_electron,
            "electron_product_over_final_mobility": ratio(electron_product, electron_final),
            "electron_final_over_low_field_mobility": ratio(electron_final, electron_low),
            "electron_final_over_sentaurus_mobility": ratio(electron_final, sentaurus_electron),
            "hole_low_field_mobility_m2_V_s": hole_low,
            "hole_mobility_limiter": hole_limiter,
            "hole_limited_mobility_product_m2_V_s": hole_product,
            "hole_final_mobility_m2_V_s": hole_final,
            "sentaurus_hole_mobility_m2_V_s": sentaurus_hole,
            "hole_product_over_final_mobility": ratio(hole_product, hole_final),
            "hole_final_over_low_field_mobility": ratio(hole_final, hole_low),
            "hole_final_over_sentaurus_mobility": ratio(hole_final, sentaurus_hole),
            "vela_electric_field_scalar_avg_V_cm": area_weighted(active, maps["vela_electric_scalar"]),
            "vela_electric_field_scalar_avg_V_m": area_weighted(active, maps["vela_electric_scalar_si"]),
            "vela_electric_field_edge_abs_V_m": vela_electric_edge,
            "sentaurus_electric_field_edge_abs_V_m": sentaurus_electric_edge,
            "electric_field_edge_abs_vela_over_sentaurus": ratio(vela_electric_edge, sentaurus_electric_edge),
            "vela_electron_high_field_drive_avg_V_cm": electron_drive,
            "vela_electron_high_field_drive_avg_V_m": electron_drive_si,
            "sentaurus_electron_qf_field_edge_abs_V_m": sentaurus_electron_qf,
            "electron_drive_vela_over_sentaurus_qf": ratio(electron_drive_si, sentaurus_electron_qf),
            "vela_hole_high_field_drive_avg_V_cm": hole_drive,
            "vela_hole_high_field_drive_avg_V_m": hole_drive_si,
            "sentaurus_hole_qf_field_edge_abs_V_m": sentaurus_hole_qf,
            "hole_drive_vela_over_sentaurus_qf": ratio(hole_drive_si, sentaurus_hole_qf),
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {"count": len(rows)}
    for key in SUMMARY_METRICS:
        values = finite_values(rows, key)
        if values:
            result[f"{key}_median"] = statistics.median(values)
            result[f"{key}_min"] = min(values)
            result[f"{key}_max"] = max(values)
    return result


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        support_class: summarize([row for row in rows if row["support_class"] == support_class])
        for support_class in sorted({str(row["support_class"]) for row in rows})
    }


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
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    summary = {
        "row_count": len(rows),
        "bias_V": args.bias,
        "active_axis": args.active_axis,
        "support_classes": class_summary(rows),
        "all": summarize(rows),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_edge_mobility_inputs_nodes.csv", rows)
    (args.out_dir / "active_edge_mobility_inputs_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
