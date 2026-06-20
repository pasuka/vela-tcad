#!/usr/bin/env python3
"""Localize active-edge SG flux ratios and their y trends."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


NODE_FIELDS = [
    "support_node_id",
    "support_class",
    "y_um",
    "edge_count",
    "edge_axes",
    "edge_ids",
    "electron_qf_over_sentaurus_edgeavg_median",
    "hole_qf_over_sentaurus_edgeavg_median",
    "particle_qf_over_sentaurus_edgeavg_median",
    "electron_qf_over_sentaurus_edgeavg_min",
    "hole_qf_over_sentaurus_edgeavg_min",
    "particle_qf_over_sentaurus_edgeavg_min",
    "electron_qf_over_sentaurus_edgeavg_max",
    "hole_qf_over_sentaurus_edgeavg_max",
    "particle_qf_over_sentaurus_edgeavg_max",
]

CORRELATION_FIELDS = [
    "electron_qf_over_sentaurus_edgeavg_median",
    "hole_qf_over_sentaurus_edgeavg_median",
    "particle_qf_over_sentaurus_edgeavg_median",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-state-csv", type=Path, required=True)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--support-class", default="overlap")
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


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


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
    if bias is None or row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def finite_values(values: list[float | None]) -> list[float]:
    return [value for value in values if value is not None and math.isfinite(value)]


def stats(values: list[float | None]) -> tuple[float | None, float | None, float | None]:
    finite = finite_values(values)
    if not finite:
        return None, None, None
    return statistics.median(finite), min(finite), max(finite)


def pearson(xs: list[float], ys: list[float]) -> float | None:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return None
    mean_x = statistics.fmean(x for x, _ in pairs)
    mean_y = statistics.fmean(y for _, y in pairs)
    dx = [x - mean_x for x, _ in pairs]
    dy = [y - mean_y for _, y in pairs]
    denom = math.sqrt(sum(x * x for x in dx) * sum(y * y for y in dy))
    if denom == 0.0:
        return None
    return sum(x * y for x, y in zip(dx, dy)) / denom


def edge_ratios(row: dict[str, str]) -> dict[str, float | None]:
    return {
        "electron": safe_ratio(
            finite_float(row.get("qf_old_slotboom_electron_flux_m2_s")),
            finite_float(row.get("sentaurus_edgeavg_electron_current_flux_m2_s")),
        ),
        "hole": safe_ratio(
            finite_float(row.get("qf_old_slotboom_hole_flux_m2_s")),
            finite_float(row.get("sentaurus_edgeavg_hole_current_flux_m2_s")),
        ),
        "particle": safe_ratio(
            finite_float(row.get("qf_old_slotboom_particle_flux_m2_s")),
            finite_float(row.get("sentaurus_edgeavg_particle_flux_m2_s")),
        ),
    }


def build_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], Counter[str]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    axis_counts: Counter[str] = Counter()
    for row in read_csv(args.edge_state_csv):
        if row.get("variant") != args.variant:
            continue
        if row.get("support_class") != args.support_class:
            continue
        if not matches_bias(row, args.bias):
            continue
        grouped[row["support_node_id"]].append(row)
        axis_counts[row.get("edge_axis", "")] += 1

    rows: list[dict[str, Any]] = []
    for node_id, items in grouped.items():
        axes = sorted({row.get("edge_axis", "") for row in items if row.get("edge_axis", "")})
        edge_ids = sorted(int(row["edge_id"]) for row in items if row.get("edge_id") not in (None, ""))
        ratios = [edge_ratios(row) for row in items]
        e_med, e_min, e_max = stats([row["electron"] for row in ratios])
        h_med, h_min, h_max = stats([row["hole"] for row in ratios])
        p_med, p_min, p_max = stats([row["particle"] for row in ratios])
        rows.append({
            "support_node_id": node_id,
            "support_class": args.support_class,
            "y_um": finite_float(items[0].get("y_um")),
            "edge_count": len(items),
            "edge_axes": ";".join(axes),
            "edge_ids": ";".join(str(edge_id) for edge_id in edge_ids),
            "electron_qf_over_sentaurus_edgeavg_median": e_med,
            "hole_qf_over_sentaurus_edgeavg_median": h_med,
            "particle_qf_over_sentaurus_edgeavg_median": p_med,
            "electron_qf_over_sentaurus_edgeavg_min": e_min,
            "hole_qf_over_sentaurus_edgeavg_min": h_min,
            "particle_qf_over_sentaurus_edgeavg_min": p_min,
            "electron_qf_over_sentaurus_edgeavg_max": e_max,
            "hole_qf_over_sentaurus_edgeavg_max": h_max,
            "particle_qf_over_sentaurus_edgeavg_max": p_max,
        })
    rows.sort(key=lambda row: finite_float(row.get("y_um")))
    return rows, axis_counts


def summarize(rows: list[dict[str, Any]], axis_counts: Counter[str]) -> dict[str, Any]:
    y = [finite_float(row.get("y_um")) for row in rows]
    correlations: dict[str, dict[str, float | None]] = {}
    for field in CORRELATION_FIELDS:
        values = [finite_float(row.get(field), math.nan) for row in rows]
        finite = [value for value in values if math.isfinite(value)]
        correlations[field] = {
            "corr_y": pearson(y, values),
            "median": statistics.median(finite) if finite else None,
            "min": min(finite) if finite else None,
            "max": max(finite) if finite else None,
        }
    return {
        "row_count": len(rows),
        "edge_count": sum(axis_counts.values()),
        "axis_counts": dict(sorted(axis_counts.items())),
        "correlations": correlations,
    }


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> int:
    args = parse_args()
    rows, axis_counts = build_rows(args)
    if not rows:
        raise SystemExit("no edge-state rows matched the requested filters")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "edge_flux_y_trend_nodes.csv", rows)
    with open_path(args.out_dir / "edge_flux_y_trend_summary.json", "w") as handle:
        json.dump(clean_json(summarize(rows, axis_counts)), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
