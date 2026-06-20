#!/usr/bin/env python3
"""Summarize PN2D BV density, mobility, and SG-flux driver ratios."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


FIELDS = [
    "bias_V",
    "driver",
    "carrier",
    "points",
    "log10_p50_vela_over_sentaurus",
    "log10_p95_abs_error",
    "classification",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-edge-csv", type=Path, required=True)
    parser.add_argument("--continuity-feedback-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", required=True)
    parser.add_argument("--mismatch-threshold-decade", type=float, default=0.1)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def value(row: dict[str, str], key: str) -> float | None:
    raw = row.get(key)
    if raw in (None, ""):
        return None
    try:
        result = float(raw)
    except ValueError:
        return None
    if not math.isfinite(result):
        return None
    return result


def log_ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or candidate == 0.0 or reference == 0.0:
        return None
    return math.log10(abs(candidate) / abs(reference))


def percentile(values: list[float], p: float) -> float | None:
    clean = sorted(values)
    if not clean:
        return None
    idx = (len(clean) - 1) * p / 100.0
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return clean[lo]
    return clean[lo] * (hi - idx) + clean[hi] * (idx - lo)


def add_ratio(
    groups: dict[tuple[float, str, str], list[float]],
    bias: float,
    driver: str,
    carrier: str,
    candidate: float | None,
    reference: float | None,
) -> None:
    ratio = log_ratio(candidate, reference)
    if ratio is not None:
        groups[(bias_key(bias), driver, carrier)].append(ratio)


def contact_edge_groups(rows: list[dict[str, str]]) -> dict[tuple[float, str, str], list[float]]:
    groups: dict[tuple[float, str, str], list[float]] = defaultdict(list)
    for row in rows:
        bias = value(row, "bias_V")
        if bias is None:
            continue
        add_ratio(
            groups,
            bias,
            "contact_edge_density",
            "electron",
            value(row, "vela_electron_density_avg_cm3"),
            value(row, "sentaurus_electron_density_cm3"),
        )
        add_ratio(
            groups,
            bias,
            "contact_edge_density",
            "hole",
            value(row, "vela_hole_density_avg_cm3"),
            value(row, "sentaurus_hole_density_cm3"),
        )
        add_ratio(
            groups,
            bias,
            "contact_edge_mobility",
            "electron",
            value(row, "vela_electron_mobility_cm2_V_s"),
            value(row, "sentaurus_electron_mobility_cm2_V_s"),
        )
        add_ratio(
            groups,
            bias,
            "contact_edge_mobility",
            "hole",
            value(row, "vela_hole_mobility_cm2_V_s"),
            value(row, "sentaurus_hole_mobility_cm2_V_s"),
        )
    return groups


def feedback_groups(rows: list[dict[str, str]]) -> dict[tuple[float, str, str], list[float]]:
    groups: dict[tuple[float, str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row.get("relation") not in {"focus", ""}:
            continue
        bias = value(row, "bias_V")
        if bias is None:
            continue
        add_ratio(
            groups,
            bias,
            "interior_sg_flux",
            "electron",
            value(row, "vela_electron_flux_abs_m2_s"),
            value(row, "sentaurus_electron_flux_abs_m2_s"),
        )
        add_ratio(
            groups,
            bias,
            "interior_sg_flux",
            "hole",
            value(row, "vela_hole_flux_abs_m2_s"),
            value(row, "sentaurus_hole_flux_abs_m2_s"),
        )
    return groups


def merge_groups(*items: dict[tuple[float, str, str], list[float]]) -> dict[tuple[float, str, str], list[float]]:
    merged: dict[tuple[float, str, str], list[float]] = defaultdict(list)
    for groups in items:
        for key, values in groups.items():
            merged[key].extend(values)
    return merged


def classify(p50: float | None, threshold: float) -> str:
    if p50 is None:
        return "missing_driver"
    if abs(p50) > threshold:
        return "driver_mismatch"
    return "driver_close"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    biases = [bias_key(bias) for bias in parse_biases(args.biases)]
    groups = merge_groups(
        contact_edge_groups(read_csv_rows(args.contact_edge_csv)),
        feedback_groups(read_csv_rows(args.continuity_feedback_csv)),
    )
    rows: list[dict[str, Any]] = []
    drivers = ["contact_edge_density", "contact_edge_mobility", "interior_sg_flux"]
    carriers = ["electron", "hole"]
    for bias in biases:
        for driver in drivers:
            for carrier in carriers:
                values = groups.get((bias, driver, carrier), [])
                p50 = percentile(values, 50)
                p95 = percentile([abs(item) for item in values], 95)
                rows.append({
                    "bias_V": bias,
                    "driver": driver,
                    "carrier": carrier,
                    "points": len(values),
                    "log10_p50_vela_over_sentaurus": p50,
                    "log10_p95_abs_error": p95,
                    "classification": classify(p50, args.mismatch_threshold_decade),
                })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "transport_driver_summary.csv", rows)
    (args.out_dir / "summary.json").write_text(
        json.dumps({"biases": biases, "rows": rows}, indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
