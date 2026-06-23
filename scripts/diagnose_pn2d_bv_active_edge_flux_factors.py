#!/usr/bin/env python3
"""Decompose active-edge SG flux factors for Vela and Sentaurus states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_flux_forms as fluxforms


ELEMENTARY_CHARGE_C = 1.602176634e-19

CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "active_axis",
    "active_edge_count",
    "active_edge_ids",
    "cxx_particle_flux_avg_m2_s",
    "vela_density_particle_flux_avg_m2_s",
    "sentaurus_density_particle_flux_avg_m2_s",
    "vela_qf_model_particle_flux_avg_m2_s",
    "sentaurus_qf_model_particle_flux_avg_m2_s",
    "sentaurus_current_particle_flux_m2_s",
    "cxx_particle_flux_over_vela_density_flux",
    "vela_density_particle_flux_over_sentaurus_state",
    "sentaurus_density_particle_flux_over_sentaurus_current",
    "vela_qf_model_particle_flux_over_sentaurus_state",
    "sentaurus_qf_model_particle_flux_over_sentaurus_current",
    "electron_density_geommean_vela_over_sentaurus",
    "hole_density_geommean_vela_over_sentaurus",
    "electron_mobility_vela_over_sentaurus",
    "hole_mobility_vela_over_sentaurus",
    "electron_qf_field_abs_vela_over_sentaurus",
    "hole_qf_field_abs_vela_over_sentaurus",
    "electric_field_abs_vela_over_sentaurus",
    "electron_density_flux_vela_over_sentaurus",
    "hole_density_flux_vela_over_sentaurus",
    "electron_qf_model_flux_vela_over_sentaurus",
    "hole_qf_model_flux_vela_over_sentaurus",
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
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
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


def edge_axis(row: dict[str, str]) -> str:
    dx = abs((optional_float(row.get("x1_um")) or 0.0) - (optional_float(row.get("x0_um")) or 0.0))
    dy = abs((optional_float(row.get("y1_um")) or 0.0) - (optional_float(row.get("y0_um")) or 0.0))
    return "x" if dx >= dy else "y"


def matches_bias(row: dict[str, str], bias: float) -> bool:
    row_bias = optional_float(row.get("bias_V"))
    if row_bias is None:
        return True
    return abs(row_bias - bias) <= 1.0e-9


def first_float(row: dict[str, str], *names: str) -> float | None:
    for name in names:
        value = optional_float(row.get(name))
        if value is not None:
            return value
    return None


def load_active_cxx_edges(edge_csv: Path, bias: float, active_axis: str, threshold_scale: float) -> dict[int, list[dict[str, Any]]]:
    by_node: dict[int, list[dict[str, Any]]] = {}
    for row in read_csv_rows(edge_csv):
        if not matches_bias(row, bias):
            continue
        electron_flux = abs(first_float(row, "electron_flux_proxy", "electron_flux_abs") or 0.0)
        hole_flux = abs(first_float(row, "hole_flux_proxy", "hole_flux_abs") or 0.0)
        electron_alpha = first_float(row, "electron_alpha_m_inv") or 0.0
        hole_alpha = first_float(row, "hole_alpha_m_inv") or 0.0
        alpha_flux = first_float(row, "source_integral", "alpha_flux")
        if alpha_flux is None:
            alpha_flux = electron_alpha * electron_flux + hole_alpha * hole_flux
        edge = {
            "edge_id": int(row["edge_id"]),
            "axis": edge_axis(row),
            "endpoint_area_m2": 0.5 * (first_float(row, "edge_area_proxy_m2", "edge_area_m2") or 0.0),
            "electron_flux_proxy": electron_flux,
            "hole_flux_proxy": hole_flux,
            "alpha_flux": alpha_flux,
        }
        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) not in (None, ""):
                by_node.setdefault(int(row[endpoint]), []).append(edge)
    active_by_node: dict[int, list[dict[str, Any]]] = {}
    for node_id, items in by_node.items():
        max_alpha_flux = max((abs(float(item["alpha_flux"])) for item in items), default=0.0)
        threshold = max_alpha_flux * threshold_scale
        active_by_node[node_id] = [
            item for item in items
            if item["axis"] == active_axis and abs(float(item["alpha_flux"])) > threshold
        ]
    return active_by_node


def load_current_fields(sentaurus_dir: Path) -> dict[int, dict[str, float]]:
    fields = {
        "electron": fluxforms.read_scalar_csv(sentaurus_dir / "fields" / "eCurrentDensity_region0.csv"),
        "hole": fluxforms.read_scalar_csv(sentaurus_dir / "fields" / "hCurrentDensity_region0.csv"),
    }
    return {
        node_id: {
            "electron_flux": abs(fields["electron"][node_id]) * 1.0e4 / ELEMENTARY_CHARGE_C,
            "hole_flux": abs(fields["hole"][node_id]) * 1.0e4 / ELEMENTARY_CHARGE_C,
        }
        for node_id in range(len(fields["electron"]))
    }


def discover_vtk(root: Path, bias: float) -> Path:
    vtks = fluxforms.discover_vela_vtks(root)
    key = fluxforms.bias_key(bias)
    if key not in vtks:
        raise FileNotFoundError(f"missing Vela VTK for {bias:g} V under {root}")
    return vtks[key]


def area_weighted(items: list[dict[str, Any]], values: dict[int, float]) -> float | None:
    numerator = 0.0
    denominator = 0.0
    for item in items:
        area = float(item["endpoint_area_m2"])
        numerator += values[int(item["edge_id"])] * area
        denominator += area
    if denominator == 0.0:
        return None
    return numerator / denominator


def cxx_area_weighted(items: list[dict[str, Any]], key: str) -> float | None:
    numerator = 0.0
    denominator = 0.0
    for item in items:
        area = float(item["endpoint_area_m2"])
        numerator += float(item[key]) * area
        denominator += area
    if denominator == 0.0:
        return None
    return numerator / denominator


def geommean(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.0
    return math.sqrt(a * b)


def edge_value_maps(rows: list[dict[str, Any]]) -> dict[str, dict[int, float]]:
    maps: dict[str, dict[int, float]] = {
        "electron_density_flux": {},
        "hole_density_flux": {},
        "particle_density_flux": {},
        "electron_qf_model_flux": {},
        "hole_qf_model_flux": {},
        "particle_qf_model_flux": {},
        "electron_density_geommean": {},
        "hole_density_geommean": {},
        "electron_mobility": {},
        "hole_mobility": {},
        "electron_qf_field_abs": {},
        "hole_qf_field_abs": {},
        "electric_field_abs": {},
    }
    for row in rows:
        edge_id = int(row["edge_id"])
        e_density = abs(float(row["electron_flux_density_abs"]))
        h_density = abs(float(row["hole_flux_density_abs"]))
        e_qf = abs(float(row["electron_flux_qf_model_abs"]))
        h_qf = abs(float(row["hole_flux_qf_model_abs"]))
        maps["electron_density_flux"][edge_id] = e_density
        maps["hole_density_flux"][edge_id] = h_density
        maps["particle_density_flux"][edge_id] = e_density + h_density
        maps["electron_qf_model_flux"][edge_id] = e_qf
        maps["hole_qf_model_flux"][edge_id] = h_qf
        maps["particle_qf_model_flux"][edge_id] = e_qf + h_qf
        maps["electron_density_geommean"][edge_id] = geommean(
            float(row["electron_density_0_m3"]),
            float(row["electron_density_1_m3"]),
        )
        maps["hole_density_geommean"][edge_id] = geommean(
            float(row["hole_density_0_m3"]),
            float(row["hole_density_1_m3"]),
        )
        maps["electron_mobility"][edge_id] = float(row["electron_mobility_m2_V_s"])
        maps["hole_mobility"][edge_id] = float(row["hole_mobility_m2_V_s"])
        maps["electron_qf_field_abs"][edge_id] = abs(float(row["electron_qf_field_V_m"]))
        maps["hole_qf_field_abs"][edge_id] = abs(float(row["hole_qf_field_V_m"]))
        maps["electric_field_abs"][edge_id] = abs(float(row["electric_field_V_m"]))
    return maps


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, contact_by_node = fluxforms.sgdiag.read_mesh(args.mesh)
    edges = {int(edge["edge_id"]): edge for edge in fluxforms.sgdiag.build_edges(nodes, triangles)}
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    ni_model = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    vela_state = fluxforms.load_vela_state(discover_vtk(args.vela_vtk_root, args.bias), len(nodes))
    sentaurus_state = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
    active_by_node = load_active_cxx_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    current_fields = load_current_fields(args.sentaurus_dir)

    support_rows = [
        row for row in read_csv_rows(args.support_csv)
        if row.get("support_class") not in (None, "", "inactive")
    ]
    needed_edge_ids = sorted({
        int(item["edge_id"])
        for support in support_rows
        for item in active_by_node.get(int(support["node_id"]), [])
    })
    missing = [edge_id for edge_id in needed_edge_ids if edge_id not in edges]
    if missing:
        raise RuntimeError(f"active edge ids not found in mesh: {missing[:8]}")

    vela_edge_rows = [
        fluxforms.edge_row(args.bias, "vela", edges[edge_id], nodes, contact_by_node, vela_state, ni_model, vt)
        for edge_id in needed_edge_ids
    ]
    sentaurus_edge_rows = [
        fluxforms.edge_row(args.bias, "sentaurus", edges[edge_id], nodes, contact_by_node, sentaurus_state, ni_model, vt)
        for edge_id in needed_edge_ids
    ]
    vela_maps = edge_value_maps(vela_edge_rows)
    sentaurus_maps = edge_value_maps(sentaurus_edge_rows)

    rows: list[dict[str, Any]] = []
    for support in support_rows:
        node_id = int(support["node_id"])
        active = active_by_node.get(node_id, [])
        cxx_electron = cxx_area_weighted(active, "electron_flux_proxy")
        cxx_hole = cxx_area_weighted(active, "hole_flux_proxy")
        cxx_particle = (
            cxx_electron + cxx_hole
            if cxx_electron is not None and cxx_hole is not None
            else None
        )
        vela_density_particle = area_weighted(active, vela_maps["particle_density_flux"])
        sentaurus_density_particle = area_weighted(active, sentaurus_maps["particle_density_flux"])
        vela_qf_particle = area_weighted(active, vela_maps["particle_qf_model_flux"])
        sentaurus_qf_particle = area_weighted(active, sentaurus_maps["particle_qf_model_flux"])
        sentaurus_current_particle = (
            current_fields[node_id]["electron_flux"] + current_fields[node_id]["hole_flux"]
            if node_id in current_fields
            else None
        )
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "active_axis": args.active_axis,
            "active_edge_count": len(active),
            "active_edge_ids": ";".join(str(int(item["edge_id"])) for item in active),
            "cxx_particle_flux_avg_m2_s": cxx_particle,
            "vela_density_particle_flux_avg_m2_s": vela_density_particle,
            "sentaurus_density_particle_flux_avg_m2_s": sentaurus_density_particle,
            "vela_qf_model_particle_flux_avg_m2_s": vela_qf_particle,
            "sentaurus_qf_model_particle_flux_avg_m2_s": sentaurus_qf_particle,
            "sentaurus_current_particle_flux_m2_s": sentaurus_current_particle,
            "cxx_particle_flux_over_vela_density_flux": ratio(cxx_particle, vela_density_particle),
            "vela_density_particle_flux_over_sentaurus_state": ratio(vela_density_particle, sentaurus_density_particle),
            "sentaurus_density_particle_flux_over_sentaurus_current": ratio(
                sentaurus_density_particle,
                sentaurus_current_particle,
            ),
            "vela_qf_model_particle_flux_over_sentaurus_state": ratio(vela_qf_particle, sentaurus_qf_particle),
            "sentaurus_qf_model_particle_flux_over_sentaurus_current": ratio(
                sentaurus_qf_particle,
                sentaurus_current_particle,
            ),
            "electron_density_geommean_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["electron_density_geommean"]),
                area_weighted(active, sentaurus_maps["electron_density_geommean"]),
            ),
            "hole_density_geommean_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["hole_density_geommean"]),
                area_weighted(active, sentaurus_maps["hole_density_geommean"]),
            ),
            "electron_mobility_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["electron_mobility"]),
                area_weighted(active, sentaurus_maps["electron_mobility"]),
            ),
            "hole_mobility_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["hole_mobility"]),
                area_weighted(active, sentaurus_maps["hole_mobility"]),
            ),
            "electron_qf_field_abs_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["electron_qf_field_abs"]),
                area_weighted(active, sentaurus_maps["electron_qf_field_abs"]),
            ),
            "hole_qf_field_abs_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["hole_qf_field_abs"]),
                area_weighted(active, sentaurus_maps["hole_qf_field_abs"]),
            ),
            "electric_field_abs_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["electric_field_abs"]),
                area_weighted(active, sentaurus_maps["electric_field_abs"]),
            ),
            "electron_density_flux_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["electron_density_flux"]),
                area_weighted(active, sentaurus_maps["electron_density_flux"]),
            ),
            "hole_density_flux_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["hole_density_flux"]),
                area_weighted(active, sentaurus_maps["hole_density_flux"]),
            ),
            "electron_qf_model_flux_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["electron_qf_model_flux"]),
                area_weighted(active, sentaurus_maps["electron_qf_model_flux"]),
            ),
            "hole_qf_model_flux_vela_over_sentaurus": ratio(
                area_weighted(active, vela_maps["hole_qf_model_flux"]),
                area_weighted(active, sentaurus_maps["hole_qf_model_flux"]),
            ),
        })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def median_or_none(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def mean_or_none(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    metric_keys = [
        "cxx_particle_flux_over_vela_density_flux",
        "vela_density_particle_flux_over_sentaurus_state",
        "sentaurus_density_particle_flux_over_sentaurus_current",
        "vela_qf_model_particle_flux_over_sentaurus_state",
        "sentaurus_qf_model_particle_flux_over_sentaurus_current",
        "electron_density_geommean_vela_over_sentaurus",
        "hole_density_geommean_vela_over_sentaurus",
        "electron_mobility_vela_over_sentaurus",
        "hole_mobility_vela_over_sentaurus",
        "electron_qf_field_abs_vela_over_sentaurus",
        "hole_qf_field_abs_vela_over_sentaurus",
        "electric_field_abs_vela_over_sentaurus",
        "electron_density_flux_vela_over_sentaurus",
        "hole_density_flux_vela_over_sentaurus",
        "electron_qf_model_flux_vela_over_sentaurus",
        "hole_qf_model_flux_vela_over_sentaurus",
    ]
    result: dict[str, Any] = {}
    for cls in sorted({str(row["support_class"]) for row in rows}):
        subset = [row for row in rows if row["support_class"] == cls]
        item: dict[str, Any] = {
            "count": len(subset),
            "active_edge_count_sum": sum(int(row["active_edge_count"]) for row in subset),
        }
        for key in metric_keys:
            values = finite_values(subset, key)
            item[f"{key}_median"] = median_or_none(values)
            item[f"{key}_mean"] = mean_or_none(values)
        result[cls] = item
    return result


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
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_edge_flux_factors_nodes.csv", rows)
    (args.out_dir / "active_edge_flux_factors_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
