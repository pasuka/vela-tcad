#!/usr/bin/env python3
"""Compute node-local source-component coefficients for carrier rows."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any


EDGE_FEATURE_FIELDS = [
    "dominant_active_axis",
    "incident_edge_count",
    "active_edge_count",
    "active_endpoint_area_fraction",
    "junction_normal_active_area_fraction",
    "junction_tangent_active_area_fraction",
    "active_edge_average_ratio_to_sentaurus",
    "active_edge_density_sum_ratio_to_sentaurus",
    "dominant_axis_average_ratio_to_sentaurus",
    "dominant_axis_density_sum_ratio_to_sentaurus",
]

NODE_FIELDS = [
    "variant",
    "bias_V",
    "node_id",
    "x",
    "y",
    "support_class",
    "impact_electron_source",
    "impact_hole_source",
    "impact_combined_source",
    "impact_electron_fraction",
    "impact_hole_fraction",
    "electron_required_source",
    "hole_required_source",
    "electron_required_over_combined",
    "hole_required_over_combined",
    "electron_required_over_electron_component",
    "hole_required_over_hole_component",
    "electron_min_norm_electron_coeff",
    "electron_min_norm_hole_coeff",
    "hole_min_norm_electron_coeff",
    "hole_min_norm_hole_coeff",
    *EDGE_FEATURE_FIELDS,
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-node-csv", type=Path, required=True)
    parser.add_argument("--carrier-term-csv", type=Path, default=None)
    parser.add_argument("--edge-feature-csv", type=Path, default=None)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def read_csv(path: Path) -> list[dict[str, str]]:
    with open_path(path, "r", newline="") as handle:
        return list(csv.DictReader(handle))


def finite_float(raw: Any, default: float = 0.0) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if denominator == 0.0:
        return None
    return numerator / denominator


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None:
        return True
    if row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def source_components(row: dict[str, str]) -> tuple[float, float, float]:
    electron = finite_float(row.get("impact_electron_source"))
    hole = finite_float(row.get("impact_hole_source"))
    combined = finite_float(row.get("impact_combined_source"))
    if combined == 0.0:
        impact = finite_float(row.get("electron_impact"))
        combined = abs(impact)
    if combined == 0.0 and (electron != 0.0 or hole != 0.0):
        combined = electron + hole
    return electron, hole, combined


def nonimpact_terms(row: dict[str, str], carrier: str) -> float:
    return (
        finite_float(row.get(f"{carrier}_flux"))
        + finite_float(row.get(f"{carrier}_recombination"))
        + finite_float(row.get(f"{carrier}_gauge"))
        + finite_float(row.get(f"{carrier}_boundary"))
    )


def min_norm_coefficients(required: float, electron_source: float, hole_source: float) -> tuple[float | None, float | None]:
    denom = electron_source * electron_source + hole_source * hole_source
    if denom == 0.0:
        return None, None
    return required * electron_source / denom, required * hole_source / denom


def load_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    rows = [
        row for row in read_csv(args.exact_node_csv)
        if row.get("variant") == args.variant and matches_bias(row, args.bias)
    ]
    if not rows:
        raise SystemExit("no exact rows matched the requested variant/bias")

    carrier_by_node: dict[str, dict[str, str]] = {}
    if args.carrier_term_csv is not None:
        carrier_by_node = {row.get("node_id", ""): row for row in read_csv(args.carrier_term_csv)}

    feature_by_node: dict[str, dict[str, str]] = {}
    if args.edge_feature_csv is not None:
        feature_by_node = {row.get("node_id", ""): row for row in read_csv(args.edge_feature_csv)}

    merged: list[dict[str, str]] = []
    component_fields = [
        "impact_electron_source",
        "impact_hole_source",
        "impact_combined_source",
    ]
    for row in rows:
        item = dict(row)
        carrier = carrier_by_node.get(item.get("node_id", ""))
        if carrier is not None:
            for field in component_fields:
                if item.get(field) in (None, "") and carrier.get(field) not in (None, ""):
                    item[field] = carrier[field]
        features = feature_by_node.get(item.get("node_id", ""))
        if features is not None:
            for field in EDGE_FEATURE_FIELDS:
                if features.get(field) not in (None, ""):
                    item[field] = features[field]
        merged.append(item)
    return merged


def coefficient_row(row: dict[str, str], variant: str, bias: float | None) -> dict[str, Any]:
    electron_source, hole_source, combined_source = source_components(row)
    electron_required = nonimpact_terms(row, "electron")
    hole_required = nonimpact_terms(row, "hole")
    electron_a, electron_b = min_norm_coefficients(electron_required, electron_source, hole_source)
    hole_a, hole_b = min_norm_coefficients(hole_required, electron_source, hole_source)
    item: dict[str, Any] = {
        "variant": variant,
        "bias_V": bias if bias is not None else row.get("bias_V"),
        "node_id": row.get("node_id"),
        "x": row.get("x", row.get("x_um", "")),
        "y": row.get("y", row.get("y_um", "")),
        "support_class": row.get("support_class", ""),
        "impact_electron_source": electron_source,
        "impact_hole_source": hole_source,
        "impact_combined_source": combined_source,
        "impact_electron_fraction": safe_ratio(electron_source, combined_source),
        "impact_hole_fraction": safe_ratio(hole_source, combined_source),
        "electron_required_source": electron_required,
        "hole_required_source": hole_required,
        "electron_required_over_combined": safe_ratio(electron_required, combined_source),
        "hole_required_over_combined": safe_ratio(hole_required, combined_source),
        "electron_required_over_electron_component": safe_ratio(electron_required, electron_source),
        "hole_required_over_hole_component": safe_ratio(hole_required, hole_source),
        "electron_min_norm_electron_coeff": electron_a,
        "electron_min_norm_hole_coeff": electron_b,
        "hole_min_norm_electron_coeff": hole_a,
        "hole_min_norm_hole_coeff": hole_b,
    }
    for field in EDGE_FEATURE_FIELDS:
        item[field] = row.get(field, "")
    return item


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    result: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            result.append(value)
    return result


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    metric_keys = [
        "impact_electron_fraction",
        "impact_hole_fraction",
        "electron_required_over_combined",
        "hole_required_over_combined",
        "electron_required_over_electron_component",
        "hole_required_over_hole_component",
        "electron_min_norm_electron_coeff",
        "electron_min_norm_hole_coeff",
        "hole_min_norm_electron_coeff",
        "hole_min_norm_hole_coeff",
        "active_endpoint_area_fraction",
        "junction_normal_active_area_fraction",
        "junction_tangent_active_area_fraction",
    ]
    summary: dict[str, Any] = {"count": len(rows)}
    for key in metric_keys:
        values = finite_values(rows, key)
        summary[f"{key}_median"] = statistics.median(values) if values else None
        summary[f"{key}_mean"] = statistics.fmean(values) if values else None
        summary[f"{key}_min"] = min(values) if values else None
        summary[f"{key}_max"] = max(values) if values else None
    return summary


def grouped_summary(rows: list[dict[str, Any]], key: str, prefix: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for value in sorted({str(row.get(key, "")) for row in rows if str(row.get(key, "")) != ""}):
        result[f"{prefix}:{value}"] = summarize([row for row in rows if str(row.get(key, "")) == value])
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
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    with open_path(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> int:
    args = parse_args()
    rows = [coefficient_row(row, args.variant, args.bias) for row in load_rows(args)]
    support_classes = grouped_summary(rows, "support_class", "support")
    feature_groups = grouped_summary(rows, "dominant_active_axis", "axis")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "source_component_coefficients_nodes.csv", rows)
    payload = clean_json({
        "variant": args.variant,
        "bias_V": args.bias,
        "row_count": len(rows),
        "support_classes": {
            key.removeprefix("support:"): value for key, value in support_classes.items()
        },
        "feature_groups": feature_groups,
    })
    write_text(
        args.out_dir / "source_component_coefficients_summary.json",
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
