#!/usr/bin/env python3
"""Compare C++ SG edge flux/alpha diagnostics with Sentaurus current density."""

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
    "incident_edge_count",
    "cxx_electron_flux_sum",
    "cxx_hole_flux_sum",
    "cxx_particle_flux_sum",
    "cxx_electron_alpha_flux_sum",
    "cxx_hole_alpha_flux_sum",
    "cxx_weighted_alpha_m_inv",
    "cxx_electron_source_integral",
    "cxx_hole_source_integral",
    "cxx_edge_source_integral",
    "cxx_source_per_flux_proxy",
    "sentaurus_impact_cm3_s",
    "sentaurus_generation_m3_s",
    "sentaurus_electron_current_density_A_cm2",
    "sentaurus_hole_current_density_A_cm2",
    "sentaurus_electron_particle_flux_m2_s",
    "sentaurus_hole_particle_flux_m2_s",
    "sentaurus_particle_flux_sum_m2_s",
    "sentaurus_weighted_alpha_m_inv",
    "cxx_flux_over_sentaurus_flux",
    "cxx_weighted_alpha_over_sentaurus",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=None)
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


def load_edge_aggregates(edge_csv: Path, bias: float | None) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    for row in read_csv_rows(edge_csv):
        if not matches_bias(row, bias):
            continue
        electron_flux = abs(optional_float(row.get("electron_flux_proxy")) or 0.0)
        hole_flux = abs(optional_float(row.get("hole_flux_proxy")) or 0.0)
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        electron_source = optional_float(row.get("electron_source_integral")) or 0.0
        hole_source = optional_float(row.get("hole_source_integral")) or 0.0
        edge_source = optional_float(row.get("edge_source_integral")) or (electron_source + hole_source)
        node0_source = optional_float(row.get("node0_source_integral"))
        node1_source = optional_float(row.get("node1_source_integral"))

        for endpoint, endpoint_source in [("node0", node0_source), ("node1", node1_source)]:
            if row.get(endpoint) in (None, ""):
                continue
            node_id = int(row[endpoint])
            item = result.setdefault(node_id, {
                "incident_edge_count": 0.0,
                "cxx_electron_flux_sum": 0.0,
                "cxx_hole_flux_sum": 0.0,
                "cxx_electron_alpha_flux_sum": 0.0,
                "cxx_hole_alpha_flux_sum": 0.0,
                "cxx_electron_source_integral": 0.0,
                "cxx_hole_source_integral": 0.0,
                "cxx_edge_source_integral": 0.0,
            })
            item["incident_edge_count"] += 1.0
            item["cxx_electron_flux_sum"] += electron_flux
            item["cxx_hole_flux_sum"] += hole_flux
            item["cxx_electron_alpha_flux_sum"] += electron_alpha * electron_flux
            item["cxx_hole_alpha_flux_sum"] += hole_alpha * hole_flux
            item["cxx_electron_source_integral"] += 0.5 * electron_source
            item["cxx_hole_source_integral"] += 0.5 * hole_source
            item["cxx_edge_source_integral"] += endpoint_source if endpoint_source is not None else 0.5 * edge_source
    return result


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    sentaurus = load_sentaurus_fields(args.sentaurus_dir)
    edge_aggregates = load_edge_aggregates(args.sg_edge_csv, args.bias)
    rows: list[dict[str, Any]] = []
    for support in read_csv_rows(args.support_csv):
        if support.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(support["node_id"])
        cxx = edge_aggregates.get(node_id, {})
        cxx_electron_flux = float(cxx.get("cxx_electron_flux_sum", 0.0))
        cxx_hole_flux = float(cxx.get("cxx_hole_flux_sum", 0.0))
        cxx_particle_flux = cxx_electron_flux + cxx_hole_flux
        cxx_alpha_flux = (
            float(cxx.get("cxx_electron_alpha_flux_sum", 0.0))
            + float(cxx.get("cxx_hole_alpha_flux_sum", 0.0))
        )
        cxx_edge_source = float(cxx.get("cxx_edge_source_integral", 0.0))
        sdev = sentaurus.get(node_id, {})
        impact_cm3_s = float(sdev.get("impact_cm3_s", 0.0))
        generation_m3_s = impact_cm3_s * 1.0e6
        electron_current = float(sdev.get("electron_current_A_cm2", 0.0))
        hole_current = float(sdev.get("hole_current_A_cm2", 0.0))
        electron_flux_m2_s = abs(electron_current) * 1.0e4 / ELEMENTARY_CHARGE_C
        hole_flux_m2_s = abs(hole_current) * 1.0e4 / ELEMENTARY_CHARGE_C
        sentaurus_particle_flux = electron_flux_m2_s + hole_flux_m2_s
        cxx_weighted_alpha = ratio(cxx_alpha_flux, cxx_particle_flux)
        sentaurus_weighted_alpha = ratio(generation_m3_s, sentaurus_particle_flux)
        rows.append({
            "node_id": node_id,
            "x_um": support.get("x_um", ""),
            "y_um": support.get("y_um", ""),
            "support_class": support.get("support_class", ""),
            "incident_edge_count": int(cxx.get("incident_edge_count", 0.0)),
            "cxx_electron_flux_sum": cxx_electron_flux,
            "cxx_hole_flux_sum": cxx_hole_flux,
            "cxx_particle_flux_sum": cxx_particle_flux,
            "cxx_electron_alpha_flux_sum": cxx.get("cxx_electron_alpha_flux_sum", 0.0),
            "cxx_hole_alpha_flux_sum": cxx.get("cxx_hole_alpha_flux_sum", 0.0),
            "cxx_weighted_alpha_m_inv": cxx_weighted_alpha,
            "cxx_electron_source_integral": cxx.get("cxx_electron_source_integral", 0.0),
            "cxx_hole_source_integral": cxx.get("cxx_hole_source_integral", 0.0),
            "cxx_edge_source_integral": cxx_edge_source,
            "cxx_source_per_flux_proxy": ratio(cxx_edge_source, cxx_particle_flux),
            "sentaurus_impact_cm3_s": impact_cm3_s,
            "sentaurus_generation_m3_s": generation_m3_s,
            "sentaurus_electron_current_density_A_cm2": electron_current,
            "sentaurus_hole_current_density_A_cm2": hole_current,
            "sentaurus_electron_particle_flux_m2_s": electron_flux_m2_s,
            "sentaurus_hole_particle_flux_m2_s": hole_flux_m2_s,
            "sentaurus_particle_flux_sum_m2_s": sentaurus_particle_flux,
            "sentaurus_weighted_alpha_m_inv": sentaurus_weighted_alpha,
            "cxx_flux_over_sentaurus_flux": ratio(cxx_particle_flux, sentaurus_particle_flux),
            "cxx_weighted_alpha_over_sentaurus": (
                ratio(cxx_weighted_alpha, sentaurus_weighted_alpha)
                if cxx_weighted_alpha is not None and sentaurus_weighted_alpha is not None
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
        "cxx_flux_over_sentaurus_flux",
        "cxx_weighted_alpha_over_sentaurus",
        "cxx_weighted_alpha_m_inv",
        "sentaurus_weighted_alpha_m_inv",
        "cxx_source_per_flux_proxy",
    ]
    for cls in sorted({str(row["support_class"]) for row in rows}):
        subset = [row for row in rows if row["support_class"] == cls]
        item: dict[str, Any] = {
            "count": len(subset),
            "incident_edge_count_sum": sum(int(row["incident_edge_count"]) for row in subset),
            "cxx_particle_flux_sum": sum(abs(float(row["cxx_particle_flux_sum"] or 0.0)) for row in subset),
            "sentaurus_particle_flux_sum_m2_s": sum(
                abs(float(row["sentaurus_particle_flux_sum_m2_s"] or 0.0)) for row in subset
            ),
            "cxx_edge_source_integral": sum(abs(float(row["cxx_edge_source_integral"] or 0.0)) for row in subset),
            "sentaurus_generation_m3_s_sum": sum(
                abs(float(row["sentaurus_generation_m3_s"] or 0.0)) for row in subset
            ),
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
        "support_classes": class_summary(rows),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "cxx_edge_flux_current_nodes.csv", rows)
    (args.out_dir / "cxx_edge_flux_current_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
