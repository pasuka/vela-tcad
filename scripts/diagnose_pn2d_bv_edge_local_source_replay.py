#!/usr/bin/env python3
"""Replay incident-edge avalanche source assignments at support nodes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any


DIRECTIONS = ["left", "right", "down", "up"]
AXES = ["x", "y"]
CONSISTENCY_SOURCE_KEYS = [
    "current_cxx_endpoint_source",
    "sentaurus_support_generation_source",
    "sentaurus_edgeavg_generation_source",
    "sentaurus_support_current_scaled_source",
    "sentaurus_edgeavg_current_scaled_source",
    "sentaurus_edgeavg_electron_current_scaled_source",
    "sentaurus_edgeavg_hole_current_scaled_source",
    "sentaurus_edgeavg_abs_electron_current_scaled_source",
    "sentaurus_edgeavg_abs_hole_current_scaled_source",
]

BASE_FIELDS = [
    "variant",
    "bias_V",
    "node_id",
    "x",
    "y",
    "support_class",
    "incident_edge_count",
    "source_scale",
    "all_electron_source",
    "all_hole_source",
    "all_combined_source",
    "x_combined_source",
    "y_combined_source",
    "left_combined_source",
    "right_combined_source",
    "down_combined_source",
    "up_combined_source",
    "electron_required_source",
    "hole_required_source",
    "electron_current_residual",
    "hole_current_residual",
    "electron_double_all_residual",
    "hole_double_all_residual",
    "electron_combined_plus_electron_residual",
    "hole_combined_plus_hole_residual",
    "electron_x_double_residual",
    "hole_x_double_residual",
    "electron_required_over_all_combined",
    "hole_required_over_all_combined",
    "electron_double_all_residual_over_combined",
    "hole_combined_plus_hole_residual_over_combined",
]

CONSISTENCY_FIELDS = CONSISTENCY_SOURCE_KEYS + [
    f"{carrier}_required_over_{key}"
    for key in CONSISTENCY_SOURCE_KEYS
    for carrier in ["electron", "hole"]
]

HYBRID_FIELDS = [
    "hybrid_e2cxx_h_edgeavg_current_electron_residual",
    "hybrid_e2cxx_h_edgeavg_current_hole_residual",
    "hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined",
    "hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined",
    "hybrid_e2cxx_h_edgeavg_current_residual_l2",
]

CSV_FIELDS = BASE_FIELDS + CONSISTENCY_FIELDS + HYBRID_FIELDS + [
    f"{direction}_{component}_source"
    for direction in DIRECTIONS
    for component in ["electron", "hole"]
] + [
    f"{axis}_{component}_source"
    for axis in AXES
    for component in ["electron", "hole"]
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-node-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--edge-current-csv", type=Path, default=None)
    parser.add_argument("--carrier-term-csv", type=Path, default=None)
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


def nonimpact_terms(row: dict[str, str], carrier: str) -> float:
    return (
        finite_float(row.get(f"{carrier}_flux"))
        + finite_float(row.get(f"{carrier}_recombination"))
        + finite_float(row.get(f"{carrier}_gauge"))
        + finite_float(row.get(f"{carrier}_boundary"))
    )


def exact_rows(path: Path, variant: str, bias: float | None) -> list[dict[str, str]]:
    rows = [
        row for row in read_csv(path)
        if row.get("variant") == variant and matches_bias(row, bias)
    ]
    if not rows:
        raise SystemExit("no exact rows matched the requested variant/bias")
    return rows


def carrier_source_by_node(path: Path | None) -> dict[int, float]:
    if path is None:
        return {}
    result: dict[int, float] = {}
    for row in read_csv(path):
        node_id_raw = row.get("node_id")
        if node_id_raw in (None, ""):
            continue
        combined = finite_float(row.get("impact_combined_source"))
        if combined == 0.0:
            combined = abs(finite_float(row.get("electron_impact")))
        result[int(node_id_raw)] = combined
    return result


def edge_axis(dx: float, dy: float) -> str:
    return "x" if abs(dx) >= abs(dy) else "y"


def edge_direction(dx: float, dy: float) -> str:
    if abs(dx) >= abs(dy):
        return "right" if dx > 0.0 else "left"
    return "up" if dy > 0.0 else "down"


def endpoint_coordinates(row: dict[str, str], endpoint: str) -> tuple[float, float]:
    suffix = "0" if endpoint == "node0" else "1"
    return (
        finite_float(row.get(f"x{suffix}_um")),
        finite_float(row.get(f"y{suffix}_um")),
    )


def edge_source(row: dict[str, str]) -> tuple[float, float, float]:
    electron = finite_float(row.get("electron_source_integral"))
    hole = finite_float(row.get("hole_source_integral"))
    combined = finite_float(row.get("edge_source_integral"))
    if combined == 0.0:
        combined = electron + hole
    return electron, hole, combined


def build_edge_contributions(edge_csv: Path, bias: float | None) -> dict[int, list[dict[str, Any]]]:
    result: dict[int, list[dict[str, Any]]] = {}
    for row in read_csv(edge_csv):
        if not matches_bias(row, bias):
            continue
        electron, hole, combined = edge_source(row)
        if electron == 0.0 and hole == 0.0 and combined == 0.0:
            continue
        for endpoint, neighbor_endpoint in [("node0", "node1"), ("node1", "node0")]:
            if row.get(endpoint) in (None, "") or row.get(neighbor_endpoint) in (None, ""):
                continue
            node_id = int(row[endpoint])
            x, y = endpoint_coordinates(row, endpoint)
            nx, ny = endpoint_coordinates(row, neighbor_endpoint)
            dx = nx - x
            dy = ny - y
            result.setdefault(node_id, []).append({
                "edge_id": row.get("edge_id", ""),
                "axis": edge_axis(dx, dy),
                "direction": edge_direction(dx, dy),
                "electron_source": 0.5 * electron,
                "hole_source": 0.5 * hole,
                "combined_source": 0.5 * combined,
            })
    return result


def empty_consistency_sums() -> dict[str, float]:
    return {key: 0.0 for key in CONSISTENCY_SOURCE_KEYS}


def build_consistency_source_sums(path: Path | None, bias: float | None) -> dict[int, dict[str, float]]:
    if path is None:
        return {}
    result: dict[int, dict[str, float]] = {}
    for row in read_csv(path):
        if not matches_bias(row, bias):
            continue
        node_id_raw = row.get("support_node_id")
        if node_id_raw in (None, ""):
            continue
        node_id = int(node_id_raw)
        endpoint_area = finite_float(row.get("endpoint_area_m2"))
        cxx_endpoint_source = finite_float(row.get("cxx_endpoint_source_integral_s_inv"))
        if cxx_endpoint_source == 0.0:
            cxx_endpoint_source = finite_float(row.get("cxx_generation_density_m3_s")) * endpoint_area
        cxx_particle_flux = finite_float(row.get("cxx_particle_flux_m2_s"))
        support_current_ratio = safe_ratio(
            finite_float(row.get("sentaurus_support_particle_flux_m2_s")),
            cxx_particle_flux,
        )
        edgeavg_current_ratio = safe_ratio(
            finite_float(row.get("sentaurus_edgeavg_particle_flux_m2_s")),
            cxx_particle_flux,
        )
        edgeavg_electron_ratio = safe_ratio(
            finite_float(row.get("sentaurus_edgeavg_electron_current_flux_m2_s")),
            cxx_particle_flux,
        )
        edgeavg_hole_ratio = safe_ratio(
            finite_float(row.get("sentaurus_edgeavg_hole_current_flux_m2_s")),
            cxx_particle_flux,
        )
        sums = result.setdefault(node_id, empty_consistency_sums())
        sums["current_cxx_endpoint_source"] += cxx_endpoint_source
        sums["sentaurus_support_generation_source"] += (
            finite_float(row.get("sentaurus_support_generation_m3_s")) * endpoint_area
        )
        sums["sentaurus_edgeavg_generation_source"] += (
            finite_float(row.get("sentaurus_edgeavg_generation_m3_s")) * endpoint_area
        )
        if support_current_ratio is not None:
            sums["sentaurus_support_current_scaled_source"] += (
                cxx_endpoint_source * support_current_ratio
            )
        if edgeavg_current_ratio is not None:
            sums["sentaurus_edgeavg_current_scaled_source"] += (
                cxx_endpoint_source * edgeavg_current_ratio
            )
        if edgeavg_electron_ratio is not None:
            sums["sentaurus_edgeavg_electron_current_scaled_source"] += (
                cxx_endpoint_source * edgeavg_electron_ratio
            )
            sums["sentaurus_edgeavg_abs_electron_current_scaled_source"] += (
                cxx_endpoint_source * abs(edgeavg_electron_ratio)
            )
        if edgeavg_hole_ratio is not None:
            sums["sentaurus_edgeavg_hole_current_scaled_source"] += (
                cxx_endpoint_source * edgeavg_hole_ratio
            )
            sums["sentaurus_edgeavg_abs_hole_current_scaled_source"] += (
                cxx_endpoint_source * abs(edgeavg_hole_ratio)
            )
    return result


def empty_sums() -> dict[str, float]:
    result = {
        "all_electron_source": 0.0,
        "all_hole_source": 0.0,
        "all_combined_source": 0.0,
    }
    for direction in DIRECTIONS:
        result[f"{direction}_electron_source"] = 0.0
        result[f"{direction}_hole_source"] = 0.0
        result[f"{direction}_combined_source"] = 0.0
    for axis in AXES:
        result[f"{axis}_electron_source"] = 0.0
        result[f"{axis}_hole_source"] = 0.0
        result[f"{axis}_combined_source"] = 0.0
    return result


def aggregate_sources(items: list[dict[str, Any]]) -> dict[str, float]:
    sums = empty_sums()
    for item in items:
        electron = float(item["electron_source"])
        hole = float(item["hole_source"])
        combined = float(item["combined_source"])
        direction = str(item["direction"])
        axis = str(item["axis"])
        sums["all_electron_source"] += electron
        sums["all_hole_source"] += hole
        sums["all_combined_source"] += combined
        sums[f"{direction}_electron_source"] += electron
        sums[f"{direction}_hole_source"] += hole
        sums[f"{direction}_combined_source"] += combined
        sums[f"{axis}_electron_source"] += electron
        sums[f"{axis}_hole_source"] += hole
        sums[f"{axis}_combined_source"] += combined
    return sums


def apply_source_scale(sums: dict[str, float], scaled_combined: float | None) -> float:
    unscaled_combined = sums["all_combined_source"]
    if scaled_combined is None or unscaled_combined == 0.0:
        return 1.0
    scale = scaled_combined / unscaled_combined
    for key in list(sums):
        if key.endswith("_source"):
            sums[key] *= scale
    return scale


def replay_row(
    row: dict[str, str],
    items: list[dict[str, Any]],
    consistency_sums: dict[str, float],
    variant: str,
    bias: float | None,
    scaled_combined: float | None,
) -> dict[str, Any]:
    sums = aggregate_sources(items)
    source_scale = apply_source_scale(sums, scaled_combined)
    consistency = empty_consistency_sums()
    consistency.update(consistency_sums)
    for key in CONSISTENCY_SOURCE_KEYS:
        consistency[key] *= source_scale
    all_combined = sums["all_combined_source"]
    electron_required = nonimpact_terms(row, "electron")
    hole_required = nonimpact_terms(row, "hole")
    electron_current = electron_required - all_combined
    hole_current = hole_required - all_combined
    electron_double = electron_required - 2.0 * all_combined
    hole_double = hole_required - 2.0 * all_combined
    electron_plus_e = electron_required - (all_combined + sums["all_electron_source"])
    hole_plus_h = hole_required - (all_combined + sums["all_hole_source"])
    x_combined = sums["x_combined_source"]
    item: dict[str, Any] = {
        "variant": variant,
        "bias_V": bias if bias is not None else row.get("bias_V"),
        "node_id": row.get("node_id"),
        "x": row.get("x", row.get("x_um", "")),
        "y": row.get("y", row.get("y_um", "")),
        "support_class": row.get("support_class", ""),
        "incident_edge_count": len(items),
        "source_scale": source_scale,
        **sums,
        "electron_required_source": electron_required,
        "hole_required_source": hole_required,
        "electron_current_residual": electron_current,
        "hole_current_residual": hole_current,
        "electron_double_all_residual": electron_double,
        "hole_double_all_residual": hole_double,
        "electron_combined_plus_electron_residual": electron_plus_e,
        "hole_combined_plus_hole_residual": hole_plus_h,
        "electron_x_double_residual": electron_required - 2.0 * x_combined,
        "hole_x_double_residual": hole_required - 2.0 * x_combined,
        "electron_required_over_all_combined": safe_ratio(electron_required, all_combined),
        "hole_required_over_all_combined": safe_ratio(hole_required, all_combined),
        "electron_double_all_residual_over_combined": safe_ratio(electron_double, all_combined),
        "hole_combined_plus_hole_residual_over_combined": safe_ratio(hole_plus_h, all_combined),
        **consistency,
    }
    for key in CONSISTENCY_SOURCE_KEYS:
        item[f"electron_required_over_{key}"] = safe_ratio(electron_required, consistency[key])
        item[f"hole_required_over_{key}"] = safe_ratio(hole_required, consistency[key])
    hybrid_electron = electron_required - 2.0 * consistency["current_cxx_endpoint_source"]
    hybrid_hole = hole_required - consistency["sentaurus_edgeavg_current_scaled_source"]
    item.update({
        "hybrid_e2cxx_h_edgeavg_current_electron_residual": hybrid_electron,
        "hybrid_e2cxx_h_edgeavg_current_hole_residual": hybrid_hole,
        "hybrid_e2cxx_h_edgeavg_current_electron_residual_over_combined": safe_ratio(
            hybrid_electron,
            all_combined,
        ),
        "hybrid_e2cxx_h_edgeavg_current_hole_residual_over_combined": safe_ratio(
            hybrid_hole,
            all_combined,
        ),
        "hybrid_e2cxx_h_edgeavg_current_residual_l2": math.sqrt(
            hybrid_electron * hybrid_electron + hybrid_hole * hybrid_hole
        ),
    })
    return item


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    metric_keys = [
        "all_combined_source",
        "x_combined_source",
        "y_combined_source",
        "left_combined_source",
        "right_combined_source",
        "electron_required_over_all_combined",
        "hole_required_over_all_combined",
        "electron_current_residual",
        "hole_current_residual",
        "electron_double_all_residual",
        "hole_combined_plus_hole_residual",
        "electron_double_all_residual_over_combined",
        "hole_combined_plus_hole_residual_over_combined",
    ]
    metric_keys.extend(CONSISTENCY_FIELDS)
    metric_keys.extend(HYBRID_FIELDS)
    item: dict[str, Any] = {"count": len(rows)}
    for key in metric_keys:
        values = finite_values(rows, key)
        item[f"{key}_median"] = statistics.median(values) if values else None
        item[f"{key}_mean"] = statistics.fmean(values) if values else None
        item[f"{key}_min"] = min(values) if values else None
        item[f"{key}_max"] = max(values) if values else None
    hybrid_electron = finite_values(
        rows,
        "hybrid_e2cxx_h_edgeavg_current_electron_residual",
    )
    hybrid_hole = finite_values(
        rows,
        "hybrid_e2cxx_h_edgeavg_current_hole_residual",
    )
    item["hybrid_e2cxx_h_edgeavg_current_electron_l2"] = math.sqrt(
        sum(value * value for value in hybrid_electron)
    )
    item["hybrid_e2cxx_h_edgeavg_current_hole_l2"] = math.sqrt(
        sum(value * value for value in hybrid_hole)
    )
    item["hybrid_e2cxx_h_edgeavg_current_combined_l2"] = math.sqrt(
        item["hybrid_e2cxx_h_edgeavg_current_electron_l2"] ** 2
        + item["hybrid_e2cxx_h_edgeavg_current_hole_l2"] ** 2
    )
    return item


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for support_class in sorted({str(row.get("support_class", "")) for row in rows}):
        result[support_class] = summarize([
            row for row in rows if str(row.get("support_class", "")) == support_class
        ])
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
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    with open_path(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> int:
    args = parse_args()
    contributions = build_edge_contributions(args.sg_edge_csv, args.bias)
    consistency_by_node = build_consistency_source_sums(args.edge_current_csv, args.bias)
    scaled_source_by_node = carrier_source_by_node(args.carrier_term_csv)
    rows = [
        replay_row(
            row,
            contributions.get(int(row["node_id"]), []),
            consistency_by_node.get(int(row["node_id"]), {}),
            args.variant,
            args.bias,
            scaled_source_by_node.get(int(row["node_id"])),
        )
        for row in exact_rows(args.exact_node_csv, args.variant, args.bias)
    ]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "edge_local_source_replay_nodes.csv", rows)
    payload = clean_json({
        "variant": args.variant,
        "bias_V": args.bias,
        "row_count": len(rows),
        "support_classes": class_summary(rows),
    })
    write_text(
        args.out_dir / "edge_local_source_replay_summary.json",
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
