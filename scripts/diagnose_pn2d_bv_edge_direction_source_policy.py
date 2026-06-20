#!/usr/bin/env python3
"""Classify SG source edges and compare nodal source-density policies."""

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
    "sentaurus_generation_m3_s",
    "incident_edge_count",
    "active_edge_count",
    "dominant_active_axis",
    "total_endpoint_area_m2",
    "active_endpoint_area_m2",
    "junction_normal_endpoint_area_m2",
    "junction_tangent_endpoint_area_m2",
    "junction_normal_active_area_m2",
    "junction_tangent_active_area_m2",
    "active_endpoint_area_fraction",
    "junction_normal_total_area_fraction",
    "junction_tangent_total_area_fraction",
    "junction_normal_active_area_fraction",
    "junction_tangent_active_area_fraction",
    "full_face_average_density_m3_s",
    "active_edge_average_density_m3_s",
    "active_edge_density_sum_m3_s",
    "dominant_axis_average_density_m3_s",
    "dominant_axis_density_sum_m3_s",
    "full_face_average_ratio_to_sentaurus",
    "active_edge_average_ratio_to_sentaurus",
    "active_edge_density_sum_ratio_to_sentaurus",
    "dominant_axis_average_ratio_to_sentaurus",
    "dominant_axis_density_sum_ratio_to_sentaurus",
]


POLICY_RATIO_KEYS = [
    "full_face_average_ratio_to_sentaurus",
    "active_edge_average_ratio_to_sentaurus",
    "active_edge_density_sum_ratio_to_sentaurus",
    "dominant_axis_average_ratio_to_sentaurus",
    "dominant_axis_density_sum_ratio_to_sentaurus",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--junction-normal-axis", choices=["x", "y"], default="x")
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


def ratio(candidate: float | None, reference: float) -> float | None:
    if candidate is None or reference == 0.0:
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


def edge_axis(row: dict[str, str]) -> str:
    dx = abs((optional_float(row.get("x1_um")) or 0.0) - (optional_float(row.get("x0_um")) or 0.0))
    dy = abs((optional_float(row.get("y1_um")) or 0.0) - (optional_float(row.get("y0_um")) or 0.0))
    return "x" if dx >= dy else "y"


def load_edge_contributions(edge_csv: Path, bias: float | None) -> dict[int, list[dict[str, float | str]]]:
    result: dict[int, list[dict[str, float | str]]] = {}
    for row in read_csv_rows(edge_csv):
        if not matches_bias(row, bias):
            continue
        area = optional_float(row.get("edge_area_proxy_m2")) or 0.0
        electron_flux = abs(optional_float(row.get("electron_flux_proxy")) or 0.0)
        hole_flux = abs(optional_float(row.get("hole_flux_proxy")) or 0.0)
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        density = electron_alpha * electron_flux + hole_alpha * hole_flux
        edge_source = optional_float(row.get("edge_source_integral"))
        if edge_source is None:
            edge_source = (optional_float(row.get("electron_source_integral")) or 0.0) + (
                optional_float(row.get("hole_source_integral")) or 0.0
            )
        endpoint_source = {
            "node0": optional_float(row.get("node0_source_integral")),
            "node1": optional_float(row.get("node1_source_integral")),
        }
        axis = edge_axis(row)
        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) in (None, ""):
                continue
            source = endpoint_source[endpoint]
            if source is None:
                source = 0.5 * edge_source
            result.setdefault(int(row[endpoint]), []).append({
                "axis": axis,
                "endpoint_area_m2": 0.5 * area,
                "source_integral": source,
                "source_density_m3_s": density,
            })
    return result


def area_sum(items: list[dict[str, float | str]]) -> float:
    return sum(float(item["endpoint_area_m2"]) for item in items)


def source_sum(items: list[dict[str, float | str]]) -> float:
    return sum(float(item["source_integral"]) for item in items)


def density_sum(items: list[dict[str, float | str]]) -> float:
    return sum(float(item["source_density_m3_s"]) for item in items)


def average_density(items: list[dict[str, float | str]]) -> float | None:
    area = area_sum(items)
    if area == 0.0:
        return None
    return source_sum(items) / area


def dominant_axis(active_items: list[dict[str, float | str]]) -> str:
    axis_density = {
        axis: density_sum([item for item in active_items if item["axis"] == axis])
        for axis in ["x", "y"]
    }
    if axis_density["x"] == axis_density["y"] == 0.0:
        return ""
    return "x" if axis_density["x"] >= axis_density["y"] else "y"


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    volumes = load_node_volumes(args.mesh)
    impacts = load_impact(args.sentaurus_dir)
    edge_contribs = load_edge_contributions(args.sg_edge_csv, args.bias)
    tangent_axis = "y" if args.junction_normal_axis == "x" else "x"
    rows: list[dict[str, Any]] = []
    for support in read_csv_rows(args.support_csv):
        if support.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(support["node_id"])
        items = edge_contribs.get(node_id, [])
        max_density = max((abs(float(item["source_density_m3_s"])) for item in items), default=0.0)
        threshold = max_density * args.active_source_relative_threshold
        active = [item for item in items if abs(float(item["source_density_m3_s"])) > threshold]
        normal = [item for item in items if item["axis"] == args.junction_normal_axis]
        tangent = [item for item in items if item["axis"] == tangent_axis]
        active_normal = [item for item in active if item["axis"] == args.junction_normal_axis]
        active_tangent = [item for item in active if item["axis"] == tangent_axis]
        dominant = dominant_axis(active)
        dominant_items = [item for item in active if item["axis"] == dominant]
        node_volume = volumes.get(node_id, 0.0)
        sentaurus_impact = impacts.get(node_id)
        if sentaurus_impact is None:
            sentaurus_impact = optional_float(support.get("sentaurus_avalanche_cm3_s")) or 0.0
        sentaurus_generation = sentaurus_impact * 1.0e6
        full_density = source_sum(items) / node_volume if node_volume != 0.0 else None
        active_average = average_density(active)
        active_density_sum = density_sum(active)
        dominant_average = average_density(dominant_items)
        dominant_density_sum = density_sum(dominant_items)
        total_area = area_sum(items)
        active_area = area_sum(active)
        normal_area = area_sum(normal)
        tangent_area = area_sum(tangent)
        active_normal_area = area_sum(active_normal)
        active_tangent_area = area_sum(active_tangent)
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "node_volume_m2": node_volume,
            "sentaurus_generation_m3_s": sentaurus_generation,
            "incident_edge_count": len(items),
            "active_edge_count": len(active),
            "dominant_active_axis": dominant,
            "total_endpoint_area_m2": total_area,
            "active_endpoint_area_m2": active_area,
            "junction_normal_endpoint_area_m2": normal_area,
            "junction_tangent_endpoint_area_m2": tangent_area,
            "junction_normal_active_area_m2": active_normal_area,
            "junction_tangent_active_area_m2": active_tangent_area,
            "active_endpoint_area_fraction": ratio(active_area, total_area),
            "junction_normal_total_area_fraction": ratio(normal_area, total_area),
            "junction_tangent_total_area_fraction": ratio(tangent_area, total_area),
            "junction_normal_active_area_fraction": ratio(active_normal_area, total_area),
            "junction_tangent_active_area_fraction": ratio(active_tangent_area, total_area),
            "full_face_average_density_m3_s": full_density,
            "active_edge_average_density_m3_s": active_average,
            "active_edge_density_sum_m3_s": active_density_sum,
            "dominant_axis_average_density_m3_s": dominant_average,
            "dominant_axis_density_sum_m3_s": dominant_density_sum,
            "full_face_average_ratio_to_sentaurus": ratio(full_density, sentaurus_generation),
            "active_edge_average_ratio_to_sentaurus": ratio(active_average, sentaurus_generation),
            "active_edge_density_sum_ratio_to_sentaurus": ratio(active_density_sum, sentaurus_generation),
            "dominant_axis_average_ratio_to_sentaurus": ratio(dominant_average, sentaurus_generation),
            "dominant_axis_density_sum_ratio_to_sentaurus": ratio(dominant_density_sum, sentaurus_generation),
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


def closest_policy(item: dict[str, Any]) -> str | None:
    best_name: str | None = None
    best_error: float | None = None
    for key in POLICY_RATIO_KEYS:
        value = optional_float(item.get(f"{key}_median"))
        if value is None:
            continue
        error = abs(math.log10(abs(value))) if value != 0.0 else math.inf
        if best_error is None or error < best_error:
            best_error = error
            best_name = key.removesuffix("_ratio_to_sentaurus")
    return best_name


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    metric_keys = [
        "active_endpoint_area_fraction",
        "junction_normal_total_area_fraction",
        "junction_tangent_total_area_fraction",
        "junction_normal_active_area_fraction",
        "junction_tangent_active_area_fraction",
        *POLICY_RATIO_KEYS,
    ]
    for cls in sorted({str(row["support_class"]) for row in rows}):
        subset = [row for row in rows if row["support_class"] == cls]
        item: dict[str, Any] = {
            "count": len(subset),
            "incident_edge_count_sum": sum(int(row["incident_edge_count"]) for row in subset),
            "active_edge_count_sum": sum(int(row["active_edge_count"]) for row in subset),
        }
        for key in metric_keys:
            values = finite_values(subset, key)
            item[f"{key}_median"] = median_or_none(values)
            item[f"{key}_mean"] = mean_or_none(values)
        item["closest_policy_by_median"] = closest_policy(item)
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
        "junction_normal_axis": args.junction_normal_axis,
        "support_classes": class_summary(rows),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "edge_direction_source_policy_nodes.csv", rows)
    (args.out_dir / "edge_direction_source_policy_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
