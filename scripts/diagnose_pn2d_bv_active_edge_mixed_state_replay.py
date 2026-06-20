#!/usr/bin/env python3
"""Replay active-edge SG flux/source under mixed Vela/Sentaurus states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from copy import deepcopy
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "variant",
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "active_edge_count",
    "active_edge_ids",
    "electron_density_geommean",
    "hole_density_geommean",
    "electron_flux_avg_m2_s",
    "hole_flux_avg_m2_s",
    "particle_flux_avg_m2_s",
    "generation_avg_m3_s",
    "electron_density_geommean_over_sentaurus",
    "hole_density_geommean_over_sentaurus",
    "electron_flux_over_sentaurus",
    "hole_flux_over_sentaurus",
    "particle_flux_over_sentaurus",
    "generation_over_sentaurus",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--active-axis", choices=["x", "y"], default="x")
    parser.add_argument("--active-source-relative-threshold", type=float, default=1.0e-12)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--electron-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--hole-qf-shift-v", type=float, default=0.0)
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


def safe_geommean(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.0
    return math.sqrt(a * b)


def avg(values: list[float], i: int, j: int) -> float:
    return 0.5 * (values[i] + values[j])


def inferred_ni_e(state: dict[str, list[float]], node_id: int, vt: float) -> float:
    factor = fluxforms.sgdiag.limited_exp((state["psi"][node_id] - state["phin"][node_id]) / vt)
    return state["n"][node_id] / factor if factor > 0.0 else 0.0


def inferred_ni_h(state: dict[str, list[float]], node_id: int, vt: float) -> float:
    factor = fluxforms.sgdiag.limited_exp((state["phip"][node_id] - state["psi"][node_id]) / vt)
    return state["p"][node_id] / factor if factor > 0.0 else 0.0


def reconstruct_density_state(
    ni_source: dict[str, list[float]],
    psi: list[float],
    phin: list[float],
    phip: list[float],
    mobility_source: dict[str, list[float]],
    vt: float,
) -> dict[str, list[float]]:
    n = [
        ni_source["nie"][i] * fluxforms.sgdiag.limited_exp((psi[i] - phin[i]) / vt)
        for i in range(len(psi))
    ]
    p = [
        ni_source["nih"][i] * fluxforms.sgdiag.limited_exp((phip[i] - psi[i]) / vt)
        for i in range(len(psi))
    ]
    return {
        "psi": list(psi),
        "phin": list(phin),
        "phip": list(phip),
        "n": n,
        "p": p,
        "mun": list(mobility_source["mun"]),
        "mup": list(mobility_source["mup"]),
    }


def inferred_ni_state(state: dict[str, list[float]], vt: float) -> dict[str, list[float]]:
    return {
        "nie": [inferred_ni_e(state, i, vt) for i in range(len(state["psi"]))],
        "nih": [inferred_ni_h(state, i, vt) for i in range(len(state["psi"]))],
    }


def shifted(values: list[float], delta: float) -> list[float]:
    return [value + delta for value in values]


def make_variants(
    vela: dict[str, list[float]],
    sentaurus: dict[str, list[float]],
    vt: float,
    electron_shift: float,
    hole_shift: float,
) -> dict[str, dict[str, list[float]]]:
    vela_ni = inferred_ni_state(vela, vt)
    sentaurus_ni = inferred_ni_state(sentaurus, vt)
    return {
        "vela_baseline": deepcopy(vela),
        "sentaurus_baseline": deepcopy(sentaurus),
        "vela_psi_sentaurus_qf": reconstruct_density_state(
            vela_ni,
            vela["psi"],
            sentaurus["phin"],
            sentaurus["phip"],
            vela,
            vt,
        ),
        "sentaurus_psi_vela_qf": reconstruct_density_state(
            sentaurus_ni,
            sentaurus["psi"],
            vela["phin"],
            vela["phip"],
            vela,
            vt,
        ),
        "vela_qf_shift": reconstruct_density_state(
            vela_ni,
            vela["psi"],
            shifted(vela["phin"], electron_shift),
            shifted(vela["phip"], hole_shift),
            vela,
            vt,
        ),
        "vela_qf_shift_sentaurus_mobility": reconstruct_density_state(
            vela_ni,
            vela["psi"],
            shifted(vela["phin"], electron_shift),
            shifted(vela["phip"], hole_shift),
            sentaurus,
            vt,
        ),
    }


def active_items_for_node(
    active_by_node: dict[int, list[dict[str, Any]]],
    mesh_edges: dict[int, dict[str, Any]],
    node_id: int,
) -> list[dict[str, Any]]:
    items: list[dict[str, Any]] = []
    for item in active_by_node.get(node_id, []):
        edge = mesh_edges.get(int(item["edge_id"]))
        if edge is None:
            continue
        items.append({**item, "node0": int(edge["node0"]), "node1": int(edge["node1"]), "length_m": edge["length_m"]})
    return items


def edge_metrics(item: dict[str, Any], state: dict[str, list[float]], vt: float) -> dict[str, float]:
    i = int(item["node0"])
    j = int(item["node1"])
    h = float(item["length_m"])
    dpsi = state["psi"][j] - state["psi"][i]
    coef_n = avg(state["mun"], i, j) * vt / h if h > 0.0 else 0.0
    coef_p = avg(state["mup"], i, j) * vt / h if h > 0.0 else 0.0
    e_flux = abs(fluxforms.sgdiag.sg_electron_flux(state["n"][i], state["n"][j], dpsi, vt, coef_n))
    h_flux = abs(fluxforms.sgdiag.sg_hole_flux(state["p"][i], state["p"][j], dpsi, vt, coef_p))
    e_field = fluxforms.sgdiag.field_between(state["phin"], i, j, h)
    h_field = fluxforms.sgdiag.field_between(state["phip"], i, j, h)
    alpha_e = fluxforms.sgdiag.van_overstraeten_alpha(e_field, 7.03e7, 7.03e7, 1.231e8, 1.231e8)
    alpha_h = fluxforms.sgdiag.van_overstraeten_alpha(h_field, 1.582e8, 6.71e7, 2.036e8, 1.693e8)
    return {
        "electron_density_geommean": safe_geommean(state["n"][i], state["n"][j]),
        "hole_density_geommean": safe_geommean(state["p"][i], state["p"][j]),
        "electron_flux": e_flux,
        "hole_flux": h_flux,
        "generation": alpha_e * e_flux + alpha_h * h_flux,
    }


def area_weighted_metrics(items: list[dict[str, Any]], state: dict[str, list[float]], vt: float) -> dict[str, float | None]:
    sums = {
        "electron_density_geommean": 0.0,
        "hole_density_geommean": 0.0,
        "electron_flux": 0.0,
        "hole_flux": 0.0,
        "generation": 0.0,
    }
    area_sum = 0.0
    for item in items:
        area = float(item["endpoint_area_m2"])
        metrics = edge_metrics(item, state, vt)
        for key in sums:
            sums[key] += metrics[key] * area
        area_sum += area
    if area_sum == 0.0:
        return {key: None for key in sums}
    result: dict[str, float | None] = {key: value / area_sum for key, value in sums.items()}
    result["particle_flux"] = result["electron_flux"] + result["hole_flux"]  # type: ignore[operator]
    return result


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, _ = fluxforms.sgdiag.read_mesh(args.mesh)
    mesh_edges = {int(edge["edge_id"]): edge for edge in fluxforms.sgdiag.build_edges(nodes, triangles)}
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    vela = fluxforms.load_vela_state(fluxfac.discover_vtk(args.vela_vtk_root, args.bias), len(nodes))
    sentaurus = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
    variants = make_variants(vela, sentaurus, vt, args.electron_qf_shift_v, args.hole_qf_shift_v)
    active_by_node = fluxfac.load_active_cxx_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    support_rows = [
        row for row in read_csv_rows(args.support_csv)
        if row.get("support_class") not in (None, "", "inactive")
    ]
    rows: list[dict[str, Any]] = []
    for support in support_rows:
        node_id = int(support["node_id"])
        active = active_items_for_node(active_by_node, mesh_edges, node_id)
        reference = area_weighted_metrics(active, variants["sentaurus_baseline"], vt)
        for name, state in variants.items():
            metrics = area_weighted_metrics(active, state, vt)
            rows.append({
                "variant": name,
                "node_id": node_id,
                "x_um": support.get("x_um", ""),
                "y_um": support.get("y_um", ""),
                "support_class": support.get("support_class", ""),
                "active_edge_count": len(active),
                "active_edge_ids": ";".join(str(int(item["edge_id"])) for item in active),
                "electron_density_geommean": metrics["electron_density_geommean"],
                "hole_density_geommean": metrics["hole_density_geommean"],
                "electron_flux_avg_m2_s": metrics["electron_flux"],
                "hole_flux_avg_m2_s": metrics["hole_flux"],
                "particle_flux_avg_m2_s": metrics["particle_flux"],
                "generation_avg_m3_s": metrics["generation"],
                "electron_density_geommean_over_sentaurus": ratio(
                    metrics["electron_density_geommean"],
                    reference["electron_density_geommean"],
                ),
                "hole_density_geommean_over_sentaurus": ratio(
                    metrics["hole_density_geommean"],
                    reference["hole_density_geommean"],
                ),
                "electron_flux_over_sentaurus": ratio(metrics["electron_flux"], reference["electron_flux"]),
                "hole_flux_over_sentaurus": ratio(metrics["hole_flux"], reference["hole_flux"]),
                "particle_flux_over_sentaurus": ratio(metrics["particle_flux"], reference["particle_flux"]),
                "generation_over_sentaurus": ratio(metrics["generation"], reference["generation"]),
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
        "electron_density_geommean_over_sentaurus",
        "hole_density_geommean_over_sentaurus",
        "electron_flux_over_sentaurus",
        "hole_flux_over_sentaurus",
        "particle_flux_over_sentaurus",
        "generation_over_sentaurus",
    ]
    result: dict[str, Any] = {}
    for cls in sorted({str(row["support_class"]) for row in rows}):
        class_rows = [row for row in rows if row["support_class"] == cls]
        variants: dict[str, Any] = {}
        for variant in sorted({str(row["variant"]) for row in class_rows}):
            subset = [row for row in class_rows if row["variant"] == variant]
            item: dict[str, Any] = {
                "count": len(subset),
                "active_edge_count_sum": sum(int(row["active_edge_count"]) for row in subset),
            }
            for key in metric_keys:
                values = finite_values(subset, key)
                item[f"{key}_median"] = median_or_none(values)
                item[f"{key}_mean"] = mean_or_none(values)
            variants[variant] = item
        result[cls] = {"variants": variants}
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
        "electron_qf_shift_v": args.electron_qf_shift_v,
        "hole_qf_shift_v": args.hole_qf_shift_v,
        "support_classes": class_summary(rows),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_edge_mixed_state_replay_nodes.csv", rows)
    (args.out_dir / "active_edge_mixed_state_replay_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
