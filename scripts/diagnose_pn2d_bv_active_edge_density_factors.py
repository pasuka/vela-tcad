#!/usr/bin/env python3
"""Split active-edge density ratios into inferred-ni and Boltzmann factors."""

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
    "electron_density_geommean_vela_over_sentaurus",
    "electron_inferred_ni_geommean_vela_over_sentaurus",
    "electron_boltzmann_factor_geommean_vela_over_sentaurus",
    "electron_psi_minus_phin_delta_V",
    "hole_density_geommean_vela_over_sentaurus",
    "hole_inferred_ni_geommean_vela_over_sentaurus",
    "hole_boltzmann_factor_geommean_vela_over_sentaurus",
    "hole_phip_minus_psi_delta_V",
    "electron_density_left_vela_over_sentaurus",
    "electron_density_right_vela_over_sentaurus",
    "hole_density_left_vela_over_sentaurus",
    "hole_density_right_vela_over_sentaurus",
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


def electron_factor(state: dict[str, list[float]], node_id: int, vt: float) -> float:
    return fluxforms.sgdiag.limited_exp((state["psi"][node_id] - state["phin"][node_id]) / vt)


def hole_factor(state: dict[str, list[float]], node_id: int, vt: float) -> float:
    return fluxforms.sgdiag.limited_exp((state["phip"][node_id] - state["psi"][node_id]) / vt)


def inferred_ni(density: float, factor: float) -> float:
    if factor <= 0.0:
        return 0.0
    return density / factor


def edge_endpoint_metrics(
    state: dict[str, list[float]],
    item: dict[str, Any],
    vt: float,
) -> dict[str, float]:
    i = int(item["node0"])
    j = int(item["node1"])
    e_factor_i = electron_factor(state, i, vt)
    e_factor_j = electron_factor(state, j, vt)
    h_factor_i = hole_factor(state, i, vt)
    h_factor_j = hole_factor(state, j, vt)
    return {
        "electron_density_geommean": safe_geommean(state["n"][i], state["n"][j]),
        "hole_density_geommean": safe_geommean(state["p"][i], state["p"][j]),
        "electron_ni_geommean": safe_geommean(
            inferred_ni(state["n"][i], e_factor_i),
            inferred_ni(state["n"][j], e_factor_j),
        ),
        "hole_ni_geommean": safe_geommean(
            inferred_ni(state["p"][i], h_factor_i),
            inferred_ni(state["p"][j], h_factor_j),
        ),
        "electron_factor_geommean": safe_geommean(e_factor_i, e_factor_j),
        "hole_factor_geommean": safe_geommean(h_factor_i, h_factor_j),
        "electron_exponent_avg": 0.5 * (
            state["psi"][i] - state["phin"][i] + state["psi"][j] - state["phin"][j]
        ),
        "hole_exponent_avg": 0.5 * (
            state["phip"][i] - state["psi"][i] + state["phip"][j] - state["psi"][j]
        ),
        "electron_density_left": state["n"][i],
        "electron_density_right": state["n"][j],
        "hole_density_left": state["p"][i],
        "hole_density_right": state["p"][j],
    }


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


def state_metric_maps(
    active_items: list[dict[str, Any]],
    state: dict[str, list[float]],
    vt: float,
) -> dict[str, dict[int, float]]:
    result: dict[str, dict[int, float]] = {}
    for item in active_items:
        edge_id = int(item["edge_id"])
        metrics = edge_endpoint_metrics(state, item, vt)
        for key, value in metrics.items():
            result.setdefault(key, {})[edge_id] = value
    return result


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, _ = fluxforms.sgdiag.read_mesh(args.mesh)
    mesh_edges = {int(edge["edge_id"]): edge for edge in fluxforms.sgdiag.build_edges(nodes, triangles)}
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    vela_state = fluxforms.load_vela_state(fluxfac.discover_vtk(args.vela_vtk_root, args.bias), len(nodes))
    sentaurus_state = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
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
        active = []
        for item in active_by_node.get(node_id, []):
            edge = mesh_edges.get(int(item["edge_id"]))
            if edge is None:
                continue
            active.append({
                **item,
                "node0": int(edge["node0"]),
                "node1": int(edge["node1"]),
            })
        vela_maps = state_metric_maps(active, vela_state, vt)
        sentaurus_maps = state_metric_maps(active, sentaurus_state, vt)

        def weighted_ratio(metric: str) -> float | None:
            return ratio(
                area_weighted(active, vela_maps.get(metric, {})),
                area_weighted(active, sentaurus_maps.get(metric, {})),
            )

        e_exp_vela = area_weighted(active, vela_maps.get("electron_exponent_avg", {}))
        e_exp_sentaurus = area_weighted(active, sentaurus_maps.get("electron_exponent_avg", {}))
        h_exp_vela = area_weighted(active, vela_maps.get("hole_exponent_avg", {}))
        h_exp_sentaurus = area_weighted(active, sentaurus_maps.get("hole_exponent_avg", {}))

        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "active_axis": args.active_axis,
            "active_edge_count": len(active),
            "active_edge_ids": ";".join(str(int(item["edge_id"])) for item in active),
            "electron_density_geommean_vela_over_sentaurus": weighted_ratio("electron_density_geommean"),
            "electron_inferred_ni_geommean_vela_over_sentaurus": weighted_ratio("electron_ni_geommean"),
            "electron_boltzmann_factor_geommean_vela_over_sentaurus": weighted_ratio("electron_factor_geommean"),
            "electron_psi_minus_phin_delta_V": (
                e_exp_vela - e_exp_sentaurus
                if e_exp_vela is not None and e_exp_sentaurus is not None
                else None
            ),
            "hole_density_geommean_vela_over_sentaurus": weighted_ratio("hole_density_geommean"),
            "hole_inferred_ni_geommean_vela_over_sentaurus": weighted_ratio("hole_ni_geommean"),
            "hole_boltzmann_factor_geommean_vela_over_sentaurus": weighted_ratio("hole_factor_geommean"),
            "hole_phip_minus_psi_delta_V": (
                h_exp_vela - h_exp_sentaurus
                if h_exp_vela is not None and h_exp_sentaurus is not None
                else None
            ),
            "electron_density_left_vela_over_sentaurus": weighted_ratio("electron_density_left"),
            "electron_density_right_vela_over_sentaurus": weighted_ratio("electron_density_right"),
            "hole_density_left_vela_over_sentaurus": weighted_ratio("hole_density_left"),
            "hole_density_right_vela_over_sentaurus": weighted_ratio("hole_density_right"),
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
        "electron_density_geommean_vela_over_sentaurus",
        "electron_inferred_ni_geommean_vela_over_sentaurus",
        "electron_boltzmann_factor_geommean_vela_over_sentaurus",
        "electron_psi_minus_phin_delta_V",
        "hole_density_geommean_vela_over_sentaurus",
        "hole_inferred_ni_geommean_vela_over_sentaurus",
        "hole_boltzmann_factor_geommean_vela_over_sentaurus",
        "hole_phip_minus_psi_delta_V",
        "electron_density_left_vela_over_sentaurus",
        "electron_density_right_vela_over_sentaurus",
        "hole_density_left_vela_over_sentaurus",
        "hole_density_right_vela_over_sentaurus",
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
    write_csv(args.out_dir / "active_edge_density_factors_nodes.csv", rows)
    (args.out_dir / "active_edge_density_factors_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
