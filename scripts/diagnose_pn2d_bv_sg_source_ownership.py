#!/usr/bin/env python3
"""Backproject SG avalanche edge sources onto thresholded support nodes."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "sentaurus_active",
    "vela_active",
    "node_volume_m2",
    "sentaurus_node_source_integral",
    "vtk_node_source_integral",
    "reconstructed_node_source_integral",
    "reconstructed_electron_source_integral",
    "reconstructed_hole_source_integral",
    "reconstructed_over_vtk",
    "incident_edge_count",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--support-source-scale", type=float, default=1.0e6)
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


def get_numeric(row: dict[str, str], *names: str) -> float:
    for name in names:
        value = optional_float(row.get(name))
        if value is not None:
            return value
    return 0.0


def load_edge_backprojection(edge_csv: Path, bias: float | None) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    for row in read_csv_rows(edge_csv):
        if bias is not None:
            row_bias = optional_float(row.get("bias_V"))
            if row_bias is not None and abs(row_bias - bias) > 1.0e-9:
                continue
        source = get_numeric(row, "source_integral", "edge_source_integral")
        electron = get_numeric(row, "electron_source_integral")
        hole = get_numeric(row, "hole_source_integral")
        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) in (None, ""):
                continue
            node_id = int(row[endpoint])
            item = result.setdefault(node_id, {
                "reconstructed_node_source_integral": 0.0,
                "reconstructed_electron_source_integral": 0.0,
                "reconstructed_hole_source_integral": 0.0,
                "incident_edge_count": 0.0,
            })
            item["reconstructed_node_source_integral"] += 0.5 * source
            item["reconstructed_electron_source_integral"] += 0.5 * electron
            item["reconstructed_hole_source_integral"] += 0.5 * hole
            item["incident_edge_count"] += 1.0
    return result


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    volumes = load_node_volumes(args.mesh)
    backprojection = load_edge_backprojection(args.sg_edge_csv, args.bias)
    rows: list[dict[str, Any]] = []
    for support in read_csv_rows(args.support_csv):
        if support.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(support["node_id"])
        volume = volumes.get(node_id, 0.0)
        sentaurus_density = optional_float(support.get("sentaurus_avalanche_cm3_s")) or 0.0
        vela_density = optional_float(support.get("vela_avalanche_cm3_s")) or 0.0
        sentaurus_integral = sentaurus_density * args.support_source_scale * volume
        vtk_integral = vela_density * args.support_source_scale * volume
        item = backprojection.get(node_id, {})
        reconstructed = float(item.get("reconstructed_node_source_integral", 0.0))
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um"),
            "y_um": support.get("y_um"),
            "support_class": support.get("support_class"),
            "sentaurus_active": support.get("sentaurus_active", ""),
            "vela_active": support.get("vela_active", ""),
            "node_volume_m2": volume,
            "sentaurus_node_source_integral": sentaurus_integral,
            "vtk_node_source_integral": vtk_integral,
            "reconstructed_node_source_integral": reconstructed,
            "reconstructed_electron_source_integral": item.get(
                "reconstructed_electron_source_integral", 0.0),
            "reconstructed_hole_source_integral": item.get(
                "reconstructed_hole_source_integral", 0.0),
            "reconstructed_over_vtk": ratio(reconstructed, vtk_integral),
            "incident_edge_count": int(item.get("incident_edge_count", 0.0)),
        })
    return rows


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for cls in sorted({str(row["support_class"]) for row in rows}):
        subset = [row for row in rows if row["support_class"] == cls]
        sentaurus = sum(abs(float(row["sentaurus_node_source_integral"] or 0.0)) for row in subset)
        vtk = sum(abs(float(row["vtk_node_source_integral"] or 0.0)) for row in subset)
        reconstructed = sum(abs(float(row["reconstructed_node_source_integral"] or 0.0)) for row in subset)
        electron = sum(abs(float(row["reconstructed_electron_source_integral"] or 0.0)) for row in subset)
        hole = sum(abs(float(row["reconstructed_hole_source_integral"] or 0.0)) for row in subset)
        result[cls] = {
            "count": len(subset),
            "sentaurus_node_source_integral": sentaurus,
            "vtk_node_source_integral": vtk,
            "reconstructed_node_source_integral": reconstructed,
            "reconstructed_electron_source_integral": electron,
            "reconstructed_hole_source_integral": hole,
            "reconstructed_over_vtk": ratio(reconstructed, vtk),
            "vtk_over_sentaurus": ratio(vtk, sentaurus),
            "reconstructed_over_sentaurus": ratio(reconstructed, sentaurus),
        }
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
    write_csv(args.out_dir / "sg_source_ownership_nodes.csv", rows)
    (args.out_dir / "sg_source_ownership_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
