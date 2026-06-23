#!/usr/bin/env python3
"""Decompose SG edge-source to nodal source-density geometry factors."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "node_volume_m2",
    "incident_edge_count",
    "cxx_endpoint_area_sum_m2",
    "cxx_endpoint_area_over_node_volume",
    "cxx_active_endpoint_area_sum_m2",
    "cxx_active_endpoint_area_fraction",
    "cxx_alpha_flux_density_sum_m3_s",
    "cxx_node_source_integral",
    "cxx_effective_source_volume_m2",
    "cxx_effective_source_volume_over_node_volume",
    "cxx_source_density_proxy_m3_s",
    "cxx_source_density_proxy_over_sentaurus_generation",
    "cxx_active_source_density_proxy_m3_s",
    "cxx_active_source_density_proxy_over_sentaurus_generation",
    "sentaurus_impact_cm3_s",
    "sentaurus_generation_m3_s",
    "sentaurus_node_source_integral",
    "cxx_source_over_sentaurus_source",
    "source_ratio_factor_product",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=None)
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


def ratio(candidate: float, reference: float) -> float | None:
    if reference == 0.0:
        return None
    return candidate / reference


def point_m(nodes: dict[int, dict[str, float]], node_id: int) -> tuple[float, float]:
    node = nodes[node_id]
    return node["x_um"] * 1.0e-6, node["y_um"] * 1.0e-6


def triangle_area(points: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = points
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))


def load_node_volumes(mesh: Path) -> dict[int, float]:
    data = json.loads(mesh.read_text())
    nodes = {
        int(node["id"]): {"x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data.get("nodes", [])
    }
    volumes = {node_id: 0.0 for node_id in nodes}
    for triangle in data.get("triangles", []):
        ids = [int(node_id) for node_id in triangle["node_ids"]]
        area = triangle_area([point_m(nodes, node_id) for node_id in ids])
        for node_id in ids:
            volumes[node_id] += area / 3.0
    return volumes


def numeric_cell(row: dict[str, str], preferred: str = "component0") -> float:
    value = optional_float(row.get(preferred))
    if value is not None:
        return value
    for key, raw in row.items():
        if key == "node_id":
            continue
        value = optional_float(raw)
        if value is not None:
            return value
    return 0.0


def load_impact(sentaurus_dir: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted((sentaurus_dir / "fields").glob("ImpactIonization_region*.csv")):
        for row in read_csv_rows(path):
            result[int(row["node_id"])] = numeric_cell(row)
    return result


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None:
        return True
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


def load_edge_geometry(edge_csv: Path, bias: float | None, active_relative_threshold: float) -> dict[int, dict[str, float]]:
    contributions: dict[int, list[dict[str, float]]] = {}
    for row in read_csv_rows(edge_csv):
        if not matches_bias(row, bias):
            continue
        edge_area = first_float(row, "edge_area_proxy_m2", "edge_area_m2") or 0.0
        electron_flux = abs(first_float(row, "electron_flux_proxy", "electron_flux_abs") or 0.0)
        hole_flux = abs(first_float(row, "hole_flux_proxy", "hole_flux_abs") or 0.0)
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        alpha_flux_density = electron_alpha * electron_flux + hole_alpha * hole_flux
        edge_source = first_float(row, "edge_source_integral", "source_integral")
        if edge_source is None:
            electron_source = optional_float(row.get("electron_source_integral")) or 0.0
            hole_source = optional_float(row.get("hole_source_integral")) or 0.0
            edge_source = electron_source + hole_source
        endpoint_sources = {
            "node0": optional_float(row.get("node0_source_integral")),
            "node1": optional_float(row.get("node1_source_integral")),
        }

        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) in (None, ""):
                continue
            node_id = int(row[endpoint])
            source = endpoint_sources[endpoint]
            if source is None:
                source = 0.5 * edge_source
            contributions.setdefault(node_id, []).append({
                "endpoint_area_m2": 0.5 * edge_area,
                "alpha_flux_density_m3_s": alpha_flux_density,
                "source_integral": source,
            })
    result: dict[int, dict[str, float]] = {}
    for node_id, items in contributions.items():
        max_density = max((abs(item["alpha_flux_density_m3_s"]) for item in items), default=0.0)
        threshold = max_density * active_relative_threshold
        active_items = [
            item for item in items
            if abs(item["alpha_flux_density_m3_s"]) > threshold
        ]
        result[node_id] = {
            "incident_edge_count": float(len(items)),
            "cxx_endpoint_area_sum_m2": sum(item["endpoint_area_m2"] for item in items),
            "cxx_active_endpoint_area_sum_m2": sum(item["endpoint_area_m2"] for item in active_items),
            "cxx_alpha_flux_density_sum_m3_s": sum(item["alpha_flux_density_m3_s"] for item in items),
            "cxx_node_source_integral": sum(item["source_integral"] for item in items),
            "cxx_active_node_source_integral": sum(item["source_integral"] for item in active_items),
        }
    return result


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    volumes = load_node_volumes(args.mesh)
    impact = load_impact(args.sentaurus_dir)
    edge_geometry = load_edge_geometry(args.sg_edge_csv, args.bias, args.active_source_relative_threshold)
    rows: list[dict[str, Any]] = []
    for support in read_csv_rows(args.support_csv):
        if support.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(support["node_id"])
        node_volume = volumes.get(node_id, 0.0)
        cxx = edge_geometry.get(node_id, {})
        endpoint_area = float(cxx.get("cxx_endpoint_area_sum_m2", 0.0))
        alpha_flux_density = float(cxx.get("cxx_alpha_flux_density_sum_m3_s", 0.0))
        node_source = float(cxx.get("cxx_node_source_integral", 0.0))
        active_endpoint_area = float(cxx.get("cxx_active_endpoint_area_sum_m2", 0.0))
        active_node_source = float(cxx.get("cxx_active_node_source_integral", 0.0))
        source_density_proxy = ratio(node_source, endpoint_area)
        active_source_density_proxy = ratio(active_node_source, active_endpoint_area)
        effective_volume = ratio(node_source, alpha_flux_density)
        sentaurus_impact = impact.get(node_id)
        if sentaurus_impact is None:
            sentaurus_impact = optional_float(support.get("sentaurus_avalanche_cm3_s")) or 0.0
        sentaurus_generation = sentaurus_impact * 1.0e6
        sentaurus_source = sentaurus_generation * node_volume
        area_over_volume = ratio(endpoint_area, node_volume)
        density_ratio = (
            ratio(source_density_proxy, sentaurus_generation)
            if source_density_proxy is not None
            else None
        )
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "node_volume_m2": node_volume,
            "incident_edge_count": int(cxx.get("incident_edge_count", 0.0)),
            "cxx_endpoint_area_sum_m2": endpoint_area,
            "cxx_endpoint_area_over_node_volume": area_over_volume,
            "cxx_active_endpoint_area_sum_m2": active_endpoint_area,
            "cxx_active_endpoint_area_fraction": ratio(active_endpoint_area, endpoint_area),
            "cxx_alpha_flux_density_sum_m3_s": alpha_flux_density,
            "cxx_node_source_integral": node_source,
            "cxx_effective_source_volume_m2": effective_volume,
            "cxx_effective_source_volume_over_node_volume": (
                ratio(effective_volume, node_volume) if effective_volume is not None else None
            ),
            "cxx_source_density_proxy_m3_s": source_density_proxy,
            "cxx_source_density_proxy_over_sentaurus_generation": density_ratio,
            "cxx_active_source_density_proxy_m3_s": active_source_density_proxy,
            "cxx_active_source_density_proxy_over_sentaurus_generation": (
                ratio(active_source_density_proxy, sentaurus_generation)
                if active_source_density_proxy is not None
                else None
            ),
            "sentaurus_impact_cm3_s": sentaurus_impact,
            "sentaurus_generation_m3_s": sentaurus_generation,
            "sentaurus_node_source_integral": sentaurus_source,
            "cxx_source_over_sentaurus_source": ratio(node_source, sentaurus_source),
            "source_ratio_factor_product": (
                area_over_volume * density_ratio
                if area_over_volume is not None and density_ratio is not None
                else None
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
    result: dict[str, Any] = {}
    metric_keys = [
        "cxx_endpoint_area_over_node_volume",
        "cxx_active_endpoint_area_fraction",
        "cxx_source_density_proxy_over_sentaurus_generation",
        "cxx_active_source_density_proxy_over_sentaurus_generation",
        "cxx_effective_source_volume_over_node_volume",
        "cxx_source_over_sentaurus_source",
        "source_ratio_factor_product",
    ]
    for cls in sorted({str(row["support_class"]) for row in rows}):
        subset = [row for row in rows if row["support_class"] == cls]
        item: dict[str, Any] = {
            "count": len(subset),
            "incident_edge_count_sum": sum(int(row["incident_edge_count"]) for row in subset),
            "node_volume_m2_sum": sum(abs(float(row["node_volume_m2"] or 0.0)) for row in subset),
            "cxx_endpoint_area_sum_m2": sum(abs(float(row["cxx_endpoint_area_sum_m2"] or 0.0)) for row in subset),
            "cxx_active_endpoint_area_sum_m2": sum(
                abs(float(row["cxx_active_endpoint_area_sum_m2"] or 0.0)) for row in subset
            ),
            "cxx_node_source_integral": sum(abs(float(row["cxx_node_source_integral"] or 0.0)) for row in subset),
            "sentaurus_node_source_integral": sum(
                abs(float(row["sentaurus_node_source_integral"] or 0.0)) for row in subset
            ),
        }
        item["summed_cxx_source_over_sentaurus_source"] = ratio(
            item["cxx_node_source_integral"],
            item["sentaurus_node_source_integral"],
        )
        item["summed_endpoint_area_over_node_volume"] = ratio(
            item["cxx_endpoint_area_sum_m2"],
            item["node_volume_m2_sum"],
        )
        item["summed_active_endpoint_area_fraction"] = ratio(
            item["cxx_active_endpoint_area_sum_m2"],
            item["cxx_endpoint_area_sum_m2"],
        )
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
        "support_classes": class_summary(rows),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "source_geometry_nodes.csv", rows)
    (args.out_dir / "source_geometry_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
