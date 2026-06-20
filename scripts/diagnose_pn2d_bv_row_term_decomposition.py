#!/usr/bin/env python3
"""Decompose carrier-row non-impact terms and their y trends."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any


NODE_FIELDS = [
    "variant",
    "bias_V",
    "node_id",
    "x",
    "y",
    "support_class",
    "source_policy",
    "source_value",
    "electron_flux",
    "electron_recombination",
    "electron_gauge_boundary",
    "electron_required_source",
    "electron_impact_source",
    "electron_required_over_source",
    "electron_residual_under_policy",
    "electron_residual_under_policy_over_source",
    "hole_flux",
    "hole_recombination",
    "hole_gauge_boundary",
    "hole_required_source",
    "hole_impact_source",
    "hole_required_over_source",
    "hole_residual_under_policy",
    "hole_residual_under_policy_over_source",
]

CORRELATION_FIELDS = [
    "electron_flux",
    "electron_recombination",
    "electron_gauge_boundary",
    "electron_required_source",
    "electron_impact_source",
    "electron_required_over_source",
    "electron_residual_under_policy",
    "electron_residual_under_policy_over_source",
    "hole_flux",
    "hole_recombination",
    "hole_gauge_boundary",
    "hole_required_source",
    "hole_impact_source",
    "hole_required_over_source",
    "hole_residual_under_policy",
    "hole_residual_under_policy_over_source",
]

DOMINANT_FIELDS = {
    "electron": [
        "electron_flux",
        "electron_recombination",
        "electron_gauge_boundary",
        "electron_impact_source",
    ],
    "hole": [
        "hole_flux",
        "hole_recombination",
        "hole_gauge_boundary",
        "hole_impact_source",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-node-csv", type=Path, required=True)
    parser.add_argument("--edge-replay-csv", type=Path, default=None)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    parser.add_argument("--support-class", default="overlap")
    parser.add_argument("--source-policy", default="all_combined_source")
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


def edge_replay_by_node(path: Path | None, variant: str, bias: float | None) -> dict[str, dict[str, str]]:
    if path is None:
        return {}
    result: dict[str, dict[str, str]] = {}
    for row in read_csv(path):
        if row.get("variant") != variant:
            continue
        if not matches_bias(row, bias):
            continue
        node_id = row.get("node_id")
        if node_id not in (None, ""):
            result[node_id] = row
    return result


def carrier_required(row: dict[str, str], carrier: str) -> float:
    return (
        finite_float(row.get(f"{carrier}_flux"))
        + finite_float(row.get(f"{carrier}_recombination"))
        + finite_float(row.get(f"{carrier}_gauge"))
        + finite_float(row.get(f"{carrier}_boundary"))
    )


def gauge_boundary(row: dict[str, str], carrier: str) -> float:
    return finite_float(row.get(f"{carrier}_gauge")) + finite_float(row.get(f"{carrier}_boundary"))


def source_value(
    exact: dict[str, str],
    replay: dict[str, str] | None,
    policy: str,
) -> float:
    if policy == "abs_electron_impact":
        return abs(finite_float(exact.get("electron_impact")))
    if policy == "abs_hole_impact":
        return abs(finite_float(exact.get("hole_impact")))
    if replay is not None and replay.get(policy) not in (None, ""):
        return finite_float(replay.get(policy))
    if exact.get(policy) not in (None, ""):
        return finite_float(exact.get(policy))
    return abs(finite_float(exact.get("electron_impact")))


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    replay_by_node = edge_replay_by_node(args.edge_replay_csv, args.variant, args.bias)
    rows: list[dict[str, Any]] = []
    for exact in read_csv(args.exact_node_csv):
        if exact.get("variant") != args.variant:
            continue
        if exact.get("support_class") != args.support_class:
            continue
        if not matches_bias(exact, args.bias):
            continue
        node_id = exact.get("node_id", "")
        replay = replay_by_node.get(node_id)
        source = source_value(exact, replay, args.source_policy)
        electron_required = carrier_required(exact, "electron")
        hole_required = carrier_required(exact, "hole")
        electron_residual = electron_required - source
        hole_residual = hole_required - source
        item = {
            "variant": args.variant,
            "bias_V": args.bias if args.bias is not None else exact.get("bias_V"),
            "node_id": node_id,
            "x": exact.get("x", exact.get("x_um", "")),
            "y": exact.get("y", exact.get("y_um", "")),
            "support_class": exact.get("support_class", ""),
            "source_policy": args.source_policy,
            "source_value": source,
            "electron_flux": finite_float(exact.get("electron_flux")),
            "electron_recombination": finite_float(exact.get("electron_recombination")),
            "electron_gauge_boundary": gauge_boundary(exact, "electron"),
            "electron_required_source": electron_required,
            "electron_impact_source": abs(finite_float(exact.get("electron_impact"))),
            "electron_required_over_source": safe_ratio(electron_required, source),
            "electron_residual_under_policy": electron_residual,
            "electron_residual_under_policy_over_source": safe_ratio(electron_residual, source),
            "hole_flux": finite_float(exact.get("hole_flux")),
            "hole_recombination": finite_float(exact.get("hole_recombination")),
            "hole_gauge_boundary": gauge_boundary(exact, "hole"),
            "hole_required_source": hole_required,
            "hole_impact_source": abs(finite_float(exact.get("hole_impact"))),
            "hole_required_over_source": safe_ratio(hole_required, source),
            "hole_residual_under_policy": hole_residual,
            "hole_residual_under_policy_over_source": safe_ratio(hole_residual, source),
        }
        rows.append(item)
    rows.sort(key=lambda row: finite_float(row.get("y")))
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


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


def metric_summary(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"median": None, "min": None, "max": None}
    return {
        "median": statistics.median(values),
        "min": min(values),
        "max": max(values),
    }


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    y = finite_values(rows, "y")
    correlations: dict[str, dict[str, Any]] = {}
    for key in CORRELATION_FIELDS:
        values = finite_values(rows, key)
        correlations[key] = {
            "corr_y": pearson(y, values),
            **metric_summary(values),
        }

    dominant: dict[str, str | None] = {}
    for carrier, keys in DOMINANT_FIELDS.items():
        ranked = [
            (abs(corr), key)
            for key in keys
            if (corr := correlations[key]["corr_y"]) is not None
        ]
        dominant[carrier] = max(ranked)[1] if ranked else None
    return {
        "row_count": len(rows),
        "correlations": correlations,
        "dominant_y_correlation": dominant,
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
    rows = build_rows(args)
    if not rows:
        raise SystemExit("no exact rows matched the requested filters")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "row_term_decomposition_nodes.csv", rows)
    with open_path(args.out_dir / "row_term_decomposition_summary.json", "w") as handle:
        json.dump(clean_json(summarize(rows)), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
