#!/usr/bin/env python3
"""Rank carrier-row avalanche source ownership policies from edge replay rows."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from itertools import product
from pathlib import Path
from typing import Any


POLICY_FIELDS = {
    "vela_combined": "all_combined_source",
    "vela_double_combined": "all_combined_source",
    "vela_combined_plus_electron": "all_combined_source",
    "vela_combined_plus_hole": "all_combined_source",
    "vela_electron": "all_electron_source",
    "vela_hole": "all_hole_source",
    "cxx_endpoint": "current_cxx_endpoint_source",
    "sentaurus_edgeavg_current": "sentaurus_edgeavg_current_scaled_source",
    "sentaurus_edgeavg_electron_current": "sentaurus_edgeavg_electron_current_scaled_source",
    "sentaurus_edgeavg_hole_current": "sentaurus_edgeavg_hole_current_scaled_source",
    "sentaurus_edgeavg_abs_electron_current": "sentaurus_edgeavg_abs_electron_current_scaled_source",
    "sentaurus_edgeavg_abs_hole_current": "sentaurus_edgeavg_abs_hole_current_scaled_source",
    "zero": "",
}

POLICIES = list(POLICY_FIELDS)

SUMMARY_FIELDS = [
    "electron_policy",
    "hole_policy",
    "count",
    "electron_l2",
    "hole_l2",
    "combined_carrier_l2",
    "electron_signed_median_over_combined",
    "hole_signed_median_over_combined",
    "electron_abs_median_over_combined",
    "hole_abs_median_over_combined",
    "electron_max_abs_over_combined",
    "hole_max_abs_over_combined",
    "electron_corr_y",
    "hole_corr_y",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edge-replay-csv", type=Path, required=True)
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
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
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


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None or row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def filtered_rows(args: argparse.Namespace) -> list[dict[str, str]]:
    rows = []
    for row in read_csv(args.edge_replay_csv):
        if row.get("variant") != args.variant:
            continue
        if row.get("support_class") != args.support_class:
            continue
        if not matches_bias(row, args.bias):
            continue
        rows.append(row)
    rows.sort(key=lambda row: finite_float(row.get("y")))
    return rows


def policy_source(row: dict[str, str], policy: str) -> float:
    if policy == "zero":
        return 0.0
    if policy == "vela_double_combined":
        return 2.0 * finite_float(row.get("all_combined_source"))
    if policy == "vela_combined_plus_electron":
        return (
            finite_float(row.get("all_combined_source"))
            + finite_float(row.get("all_electron_source"))
        )
    if policy == "vela_combined_plus_hole":
        return (
            finite_float(row.get("all_combined_source"))
            + finite_float(row.get("all_hole_source"))
        )
    return finite_float(row.get(POLICY_FIELDS[policy]))


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if denominator == 0.0:
        return None
    return numerator / denominator


def finite_values(values: list[float | None]) -> list[float]:
    return [value for value in values if value is not None and math.isfinite(value)]


def median(values: list[float | None]) -> float | None:
    finite = finite_values(values)
    return statistics.median(finite) if finite else None


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


def evaluate(rows: list[dict[str, str]], electron_policy: str, hole_policy: str) -> dict[str, Any]:
    electron_residuals: list[float] = []
    hole_residuals: list[float] = []
    electron_normed: list[float | None] = []
    hole_normed: list[float | None] = []
    y_values: list[float] = []
    for row in rows:
        electron_required = finite_float(row.get("electron_required_source"))
        hole_required = finite_float(row.get("hole_required_source"))
        electron_residual = electron_required - policy_source(row, electron_policy)
        hole_residual = hole_required - policy_source(row, hole_policy)
        combined = finite_float(row.get("all_combined_source"))
        electron_residuals.append(electron_residual)
        hole_residuals.append(hole_residual)
        electron_normed.append(safe_ratio(electron_residual, combined))
        hole_normed.append(safe_ratio(hole_residual, combined))
        y_values.append(finite_float(row.get("y")))

    electron_l2 = math.sqrt(sum(value * value for value in electron_residuals))
    hole_l2 = math.sqrt(sum(value * value for value in hole_residuals))
    electron_normed_finite = finite_values(electron_normed)
    hole_normed_finite = finite_values(hole_normed)
    return {
        "electron_policy": electron_policy,
        "hole_policy": hole_policy,
        "count": len(rows),
        "electron_l2": electron_l2,
        "hole_l2": hole_l2,
        "combined_carrier_l2": math.sqrt(electron_l2 * electron_l2 + hole_l2 * hole_l2),
        "electron_signed_median_over_combined": median(electron_normed),
        "hole_signed_median_over_combined": median(hole_normed),
        "electron_abs_median_over_combined": median([abs(value) for value in electron_normed_finite]),
        "hole_abs_median_over_combined": median([abs(value) for value in hole_normed_finite]),
        "electron_max_abs_over_combined": max([abs(value) for value in electron_normed_finite], default=None),
        "hole_max_abs_over_combined": max([abs(value) for value in hole_normed_finite], default=None),
        "electron_corr_y": pearson(y_values, electron_normed_finite),
        "hole_corr_y": pearson(y_values, hole_normed_finite),
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
    rows = filtered_rows(args)
    if not rows:
        raise SystemExit("no edge replay rows matched the requested filters")

    summary_rows = [
        evaluate(rows, electron_policy, hole_policy)
        for electron_policy, hole_policy in product(POLICIES, POLICIES)
    ]
    summary_rows.sort(key=lambda row: (
        finite_float(row.get("combined_carrier_l2"), math.inf),
        POLICIES.index(str(row.get("electron_policy"))),
        POLICIES.index(str(row.get("hole_policy"))),
    ))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "source_ownership_replay_summary.csv", summary_rows)
    summary = {
        "variant": args.variant,
        "bias_V": args.bias,
        "support_class": args.support_class,
        "row_count": len(rows),
        "policies": POLICIES,
        "best_by_combined_l2": summary_rows[0],
        "top5_by_combined_l2": summary_rows[:5],
    }
    with open_path(args.out_dir / "source_ownership_replay_summary.json", "w") as handle:
        json.dump(clean_json(summary), handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
