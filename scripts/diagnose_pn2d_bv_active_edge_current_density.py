#!/usr/bin/env python3
"""Compare junction-normal active SG edge fluxes with Sentaurus current densities."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


ELEMENTARY_CHARGE_C = 1.602176634e-19

CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "active_axis",
    "incident_edge_count",
    "active_edge_count",
    "active_endpoint_area_m2",
    "active_electron_flux_avg_m2_s",
    "active_hole_flux_avg_m2_s",
    "active_particle_flux_avg_m2_s",
    "active_generation_avg_m3_s",
    "active_weighted_alpha_m_inv",
    "sentaurus_impact_cm3_s",
    "sentaurus_generation_m3_s",
    "sentaurus_electron_current_density_A_cm2",
    "sentaurus_hole_current_density_A_cm2",
    "sentaurus_total_current_density_A_cm2",
    "sentaurus_electron_particle_flux_m2_s",
    "sentaurus_hole_particle_flux_m2_s",
    "sentaurus_eh_particle_flux_m2_s",
    "sentaurus_total_current_equiv_particle_flux_m2_s",
    "sentaurus_weighted_alpha_m_inv",
    "active_electron_flux_over_sentaurus_electron_flux",
    "active_hole_flux_over_sentaurus_hole_flux",
    "active_particle_flux_over_sentaurus_eh_flux",
    "active_particle_flux_over_sentaurus_total_current_equiv",
    "active_generation_over_sentaurus_generation",
    "active_weighted_alpha_over_sentaurus",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=None)
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


def ratio(candidate: float | None, reference: float) -> float | None:
    if candidate is None or reference == 0.0:
        return None
    return candidate / reference


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


def load_field(sentaurus_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted((sentaurus_dir / "fields").glob(f"{name}_region*.csv")):
        for row in read_csv_rows(path):
            result[int(row["node_id"])] = numeric_cell(row)
    return result


def load_sentaurus_fields(sentaurus_dir: Path) -> dict[int, dict[str, float]]:
    fields = {
        "impact_cm3_s": load_field(sentaurus_dir, "ImpactIonization"),
        "electron_current_A_cm2": load_field(sentaurus_dir, "eCurrentDensity"),
        "hole_current_A_cm2": load_field(sentaurus_dir, "hCurrentDensity"),
        "total_current_A_cm2": load_field(sentaurus_dir, "TotalCurrentDensity"),
    }
    node_ids: set[int] = set()
    for values in fields.values():
        node_ids.update(values.keys())
    return {
        node_id: {name: values.get(node_id, 0.0) for name, values in fields.items()}
        for node_id in sorted(node_ids)
    }


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
        electron_flux = abs(optional_float(row.get("electron_flux_proxy")) or 0.0)
        hole_flux = abs(optional_float(row.get("hole_flux_proxy")) or 0.0)
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        alpha_flux = electron_alpha * electron_flux + hole_alpha * hole_flux
        endpoint_area = 0.5 * (optional_float(row.get("edge_area_proxy_m2")) or 0.0)
        axis = edge_axis(row)
        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) in (None, ""):
                continue
            result.setdefault(int(row[endpoint]), []).append({
                "axis": axis,
                "endpoint_area_m2": endpoint_area,
                "electron_flux_m2_s": electron_flux,
                "hole_flux_m2_s": hole_flux,
                "alpha_flux_m3_s": alpha_flux,
            })
    return result


def area_weighted_average(items: list[dict[str, float | str]], key: str) -> float | None:
    area = sum(float(item["endpoint_area_m2"]) for item in items)
    if area == 0.0:
        return None
    return sum(float(item[key]) * float(item["endpoint_area_m2"]) for item in items) / area


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    sentaurus = load_sentaurus_fields(args.sentaurus_dir)
    edge_contribs = load_edge_contributions(args.sg_edge_csv, args.bias)
    rows: list[dict[str, Any]] = []
    for support in read_csv_rows(args.support_csv):
        if support.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(support["node_id"])
        items = edge_contribs.get(node_id, [])
        max_alpha_flux = max((abs(float(item["alpha_flux_m3_s"])) for item in items), default=0.0)
        threshold = max_alpha_flux * args.active_source_relative_threshold
        active = [
            item for item in items
            if item["axis"] == args.active_axis and abs(float(item["alpha_flux_m3_s"])) > threshold
        ]
        active_area = sum(float(item["endpoint_area_m2"]) for item in active)
        active_electron_flux = area_weighted_average(active, "electron_flux_m2_s")
        active_hole_flux = area_weighted_average(active, "hole_flux_m2_s")
        active_generation = area_weighted_average(active, "alpha_flux_m3_s")
        active_particle_flux = (
            active_electron_flux + active_hole_flux
            if active_electron_flux is not None and active_hole_flux is not None
            else None
        )
        active_alpha = (
            ratio(active_generation, active_particle_flux)
            if active_particle_flux is not None
            else None
        )
        sdev = sentaurus.get(node_id, {})
        impact_cm3_s = float(sdev.get("impact_cm3_s", 0.0))
        generation_m3_s = impact_cm3_s * 1.0e6
        electron_current = float(sdev.get("electron_current_A_cm2", 0.0))
        hole_current = float(sdev.get("hole_current_A_cm2", 0.0))
        total_current = float(sdev.get("total_current_A_cm2", 0.0))
        sentaurus_electron_flux = abs(electron_current) * 1.0e4 / ELEMENTARY_CHARGE_C
        sentaurus_hole_flux = abs(hole_current) * 1.0e4 / ELEMENTARY_CHARGE_C
        sentaurus_eh_flux = sentaurus_electron_flux + sentaurus_hole_flux
        sentaurus_total_equiv_flux = abs(total_current) * 1.0e4 / ELEMENTARY_CHARGE_C
        sentaurus_alpha = ratio(generation_m3_s, sentaurus_eh_flux)
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "active_axis": args.active_axis,
            "incident_edge_count": len(items),
            "active_edge_count": len(active),
            "active_endpoint_area_m2": active_area,
            "active_electron_flux_avg_m2_s": active_electron_flux,
            "active_hole_flux_avg_m2_s": active_hole_flux,
            "active_particle_flux_avg_m2_s": active_particle_flux,
            "active_generation_avg_m3_s": active_generation,
            "active_weighted_alpha_m_inv": active_alpha,
            "sentaurus_impact_cm3_s": impact_cm3_s,
            "sentaurus_generation_m3_s": generation_m3_s,
            "sentaurus_electron_current_density_A_cm2": electron_current,
            "sentaurus_hole_current_density_A_cm2": hole_current,
            "sentaurus_total_current_density_A_cm2": total_current,
            "sentaurus_electron_particle_flux_m2_s": sentaurus_electron_flux,
            "sentaurus_hole_particle_flux_m2_s": sentaurus_hole_flux,
            "sentaurus_eh_particle_flux_m2_s": sentaurus_eh_flux,
            "sentaurus_total_current_equiv_particle_flux_m2_s": sentaurus_total_equiv_flux,
            "sentaurus_weighted_alpha_m_inv": sentaurus_alpha,
            "active_electron_flux_over_sentaurus_electron_flux": ratio(active_electron_flux, sentaurus_electron_flux),
            "active_hole_flux_over_sentaurus_hole_flux": ratio(active_hole_flux, sentaurus_hole_flux),
            "active_particle_flux_over_sentaurus_eh_flux": ratio(active_particle_flux, sentaurus_eh_flux),
            "active_particle_flux_over_sentaurus_total_current_equiv": ratio(
                active_particle_flux,
                sentaurus_total_equiv_flux,
            ),
            "active_generation_over_sentaurus_generation": ratio(active_generation, generation_m3_s),
            "active_weighted_alpha_over_sentaurus": ratio(active_alpha, sentaurus_alpha) if sentaurus_alpha else None,
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
        "active_electron_flux_over_sentaurus_electron_flux",
        "active_hole_flux_over_sentaurus_hole_flux",
        "active_particle_flux_over_sentaurus_eh_flux",
        "active_particle_flux_over_sentaurus_total_current_equiv",
        "active_generation_over_sentaurus_generation",
        "active_weighted_alpha_over_sentaurus",
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
    write_csv(args.out_dir / "active_edge_current_density_nodes.csv", rows)
    (args.out_dir / "active_edge_current_density_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
