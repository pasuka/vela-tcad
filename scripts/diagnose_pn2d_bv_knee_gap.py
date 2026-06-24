#!/usr/bin/env python
"""Classify PN2D BV knee gap from voltage-stepped curve evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reference_curves"
    / "pn2d_sentaurus2018_bv_reference.csv"
)


def float_or_none(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def load_curve(path: Path) -> list[tuple[float, float]]:
    rows: dict[float, float] = {}
    with path.resolve().open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            bias = float_or_none(row.get("bias_V"))
            current = float_or_none(row.get("current_total_A_per_um"))
            if current is None:
                current = float_or_none(row.get("current_total"))
            if bias is None or current is None or current == 0.0:
                continue
            previous = rows.get(bias)
            if previous is None or abs(current) > abs(previous):
                rows[bias] = current
    return sorted(rows.items())


def log_abs_at(points: list[tuple[float, float]], bias: float) -> float | None:
    if not points:
        return None
    ordered = sorted(points)
    for b, current in ordered:
        if abs(b - bias) <= 1.0e-9:
            return math.log10(abs(current))
    if bias < ordered[0][0] or bias > ordered[-1][0]:
        return None
    for (b0, i0), (b1, i1) in zip(ordered, ordered[1:]):
        lo = min(b0, b1)
        hi = max(b0, b1)
        if lo <= bias <= hi and b0 != b1:
            t = (bias - b0) / (b1 - b0)
            return (1.0 - t) * math.log10(abs(i0)) + t * math.log10(abs(i1))
    return None


def abs_current_at(points: list[tuple[float, float]], bias: float) -> float | None:
    value = log_abs_at(points, bias)
    return None if value is None else 10.0 ** value


def growth_ratio(points: list[tuple[float, float]], bias: float) -> float | None:
    here = log_abs_at(points, bias)
    prev = log_abs_at(points, bias + 1.0)
    if here is None or prev is None:
        return None
    return 10.0 ** (here - prev)


def local_log_slope(points: list[tuple[float, float]], bias: float) -> float | None:
    here = log_abs_at(points, bias)
    prev = log_abs_at(points, bias + 1.0)
    if here is None or prev is None:
        return None
    return here - prev


def first_growth_bias(points: list[tuple[float, float]], threshold: float,
                      bias_min: float, bias_max: float) -> float | None:
    low = min(bias_min, bias_max)
    high = max(bias_min, bias_max)
    for bias in range(int(math.floor(high)) - 1, int(math.ceil(low)) - 1, -1):
        ratio = growth_ratio(points, float(bias))
        if ratio is not None and ratio > threshold:
            return float(bias)
    return None


def build_per_volt_table(candidate: list[tuple[float, float]],
                         reference: list[tuple[float, float]],
                         bias_min: float,
                         bias_max: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    low = int(math.floor(min(bias_min, bias_max)))
    high = int(math.ceil(max(bias_min, bias_max)))
    for bias in range(high, low - 1, -1):
        b = float(bias)
        cand_abs = abs_current_at(candidate, b)
        ref_abs = abs_current_at(reference, b)
        log_error = None
        current_ratio = None
        if cand_abs is not None and ref_abs is not None and cand_abs > 0.0 and ref_abs > 0.0:
            log_error = math.log10(cand_abs) - math.log10(ref_abs)
            current_ratio = cand_abs / ref_abs
        rows.append({
            "bias_V": b,
            "vela_abs_current": cand_abs,
            "reference_abs_current": ref_abs,
            "vela_over_reference": current_ratio,
            "log10_error_decades": log_error,
            "vela_growth_1V": growth_ratio(candidate, b),
            "reference_growth_1V": growth_ratio(reference, b),
            "vela_dlog10_absI_per_V": local_log_slope(candidate, b),
            "reference_dlog10_absI_per_V": local_log_slope(reference, b),
        })
    return rows


def max_abs_log_error(table: list[dict[str, Any]]) -> float | None:
    values = [abs(row["log10_error_decades"]) for row in table
              if row.get("log10_error_decades") is not None]
    return max(values) if values else None


def summarize_curve(label: str, path: Path, reference: list[tuple[float, float]],
                    bias_min: float, bias_max: float) -> dict[str, Any]:
    points = load_curve(path)
    table = build_per_volt_table(points, reference, bias_min, bias_max)
    converged_rows = 0
    deepest = None
    failure_reason = ""
    if path.exists():
        with path.resolve().open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                bias = float_or_none(row.get("bias_V"))
                if str(row.get("converged", "")).strip() in {"1", "true", "True"}:
                    converged_rows += 1
                    if bias is not None:
                        deepest = bias if deepest is None else min(deepest, bias)
                if not failure_reason:
                    failure_reason = (row.get("failure_reason") or row.get("breakdown_failure_reason") or "").strip()
    vela_growth_1p5 = first_growth_bias(points, 1.5, bias_min, bias_max)
    vela_growth_2p0 = first_growth_bias(points, 2.0, bias_min, bias_max)
    peak_slope = None
    slopes = [row["vela_dlog10_absI_per_V"] for row in table
              if row.get("vela_dlog10_absI_per_V") is not None]
    if slopes:
        peak_slope = max(slopes)
    return {
        "label": label,
        "path": str(path.resolve()),
        "rows_in_window": len(table),
        "converged_rows": converged_rows,
        "deepest_bias_V": deepest,
        "reached_minus20V": deepest is not None and deepest <= -19.99,
        "failure_reason": failure_reason,
        "first_growth_over_1p5_V": vela_growth_1p5,
        "first_growth_over_2p0_V": vela_growth_2p0,
        "peak_dlog10_absI_per_V": peak_slope,
        "max_abs_log10_error_decades": max_abs_log_error(table),
    }


def classify(summaries: list[dict[str, Any]]) -> str:
    # A latent turning point requires an actual steepening signal plus voltage-step failure
    # near that steep region. Reaching -20 V smoothly is treated as physics magnitude/shape gap.
    for row in summaries:
        steep = (row.get("first_growth_over_2p0_V") is not None or
                 (row.get("peak_dlog10_absI_per_V") or 0.0) > math.log10(2.0))
        stopped = not bool(row.get("reached_minus20V")) and bool(row.get("failure_reason"))
        if steep and stopped:
            return "latent_turning_point"
    return "physics_magnitude_gap"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def parse_scan(values: list[str]) -> list[tuple[str, Path]]:
    parsed: list[tuple[str, Path]] = []
    for value in values:
        if "=" not in value:
            raise SystemExit(f"invalid --scan entry, expected LABEL=PATH: {value}")
        label, raw_path = value.split("=", 1)
        parsed.append((label, Path(raw_path).resolve()))
    return parsed


def self_test() -> int:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        ref = root / "ref.csv"
        cand = root / "cand.csv"
        ref.write_text(
            "bias_V,current_total\n-20,-8\n-19,-3\n-18,-1\n-17,-0.5\n-16,-0.25\n-15,-0.125\n-14,-0.0625\n-13,-0.03125\n-12,-0.015625\n-11,-0.0078125\n-10,-0.00390625\n",
            encoding="utf-8",
        )
        cand.write_text(
            "bias_V,current_total,converged\n-20,-2\n-19,-1.4\n-18,-1\n-17,-0.8\n-16,-0.6\n-15,-0.5\n-14,-0.4\n-13,-0.3\n-12,-0.2\n-11,-0.15\n-10,-0.1\n",
            encoding="utf-8",
        )
        reference = load_curve(ref)
        table = build_per_volt_table(load_curve(cand), reference, -20.0, -10.0)
        assert len(table) == 11
        summary = summarize_curve("cand", cand, reference, -20.0, -10.0)
        assert summary["first_growth_over_1p5_V"] is None
        assert classify([summary]) == "physics_magnitude_gap"
    print("self-test passed")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--scan", action="append", default=[], metavar="LABEL=PATH")
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-10.0)
    parser.add_argument("--out-dir", type=Path, required=False)
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if args.baseline is None:
        raise SystemExit("--baseline is required unless --self-test is used")
    if args.out_dir is None:
        raise SystemExit("--out-dir is required unless --self-test is used")

    reference_path = args.reference.resolve()
    baseline_path = args.baseline.resolve()
    out_dir = args.out_dir.resolve()
    reference = load_curve(reference_path)
    baseline = load_curve(baseline_path)
    per_volt = build_per_volt_table(baseline, reference, args.bias_min, args.bias_max)
    write_csv(out_dir / "baseline_knee_table.csv", per_volt)

    summaries = [summarize_curve("baseline", baseline_path, reference, args.bias_min, args.bias_max)]
    for label, path in parse_scan(args.scan):
        summaries.append(summarize_curve(label, path, reference, args.bias_min, args.bias_max))
    write_csv(out_dir / "knee_gap_scan_summary.csv", summaries)
    payload = {
        "classification": classify(summaries),
        "reference": str(reference_path),
        "baseline": str(baseline_path),
        "bias_window_V": {"min": args.bias_min, "max": args.bias_max},
        "summaries": summaries,
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "knee_gap_classification.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())