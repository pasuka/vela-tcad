#!/usr/bin/env python
"""Summarize PN2D BV curve knee shape against the Sentaurus reference."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import NamedTuple


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "reference_curves"
    / "pn2d_sentaurus2018_bv_reference.csv"
)
DEFAULT_CANDIDATE = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "vela"
    / "pn2d_sentaurus2018_bv_minus20_avaljac.csv"
)


class CurveSummary(NamedTuple):
    first_growth_over_1p5: float | None
    first_growth_over_2p0: float | None


def _float_or_none(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def load_curve(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows: dict[float, float] = {}
        for row in reader:
            bias = _float_or_none(row.get("bias_V"))
            current = _float_or_none(row.get("current_total_A_per_um"))
            if current is None:
                current = _float_or_none(row.get("current_total"))
            if bias is None or current is None or current == 0.0:
                continue
            # Multiple Sentaurus rows can exist at 0 V. Keep the largest
            # magnitude current at each bias; the BV knee audit starts at -10 V.
            previous = rows.get(bias)
            if previous is None or abs(current) > abs(previous):
                rows[bias] = current
    return sorted(rows.items())


def _interpolate_log_abs_current(points: list[tuple[float, float]], bias: float) -> float | None:
    if not points:
        return None
    for point_bias, current in points:
        if abs(point_bias - bias) <= 1.0e-9:
            return math.log10(abs(current))

    ordered = sorted(points)
    if bias < ordered[0][0] or bias > ordered[-1][0]:
        return None
    for (b0, i0), (b1, i1) in zip(ordered, ordered[1:]):
        lo = min(b0, b1)
        hi = max(b0, b1)
        if lo <= bias <= hi and b0 != b1:
            t = (bias - b0) / (b1 - b0)
            return (1.0 - t) * math.log10(abs(i0)) + t * math.log10(abs(i1))
    return None


def _current_ratio(points: list[tuple[float, float]], bias: float) -> float | None:
    here = _interpolate_log_abs_current(points, bias)
    previous = _interpolate_log_abs_current(points, bias + 1.0)
    if here is None or previous is None:
        return None
    return 10.0 ** (here - previous)


def first_growth_bias(points: list[tuple[float, float]],
                      threshold: float,
                      bias_min: float,
                      bias_max: float) -> float | None:
    low = min(bias_min, bias_max)
    high = max(bias_min, bias_max)
    start = int(math.floor(high)) - 1
    stop = int(math.ceil(low)) - 1
    for bias in range(start, stop, -1):
        ratio = _current_ratio(points, float(bias))
        if ratio is not None and ratio > threshold:
            return float(bias)
    return None


def summarize_curve(points: list[tuple[float, float]],
                    bias_min: float,
                    bias_max: float) -> CurveSummary:
    return CurveSummary(
        first_growth_over_1p5=first_growth_bias(points, 1.5, bias_min, bias_max),
        first_growth_over_2p0=first_growth_bias(points, 2.0, bias_min, bias_max),
    )


def max_abs_log10_error(candidate: list[tuple[float, float]],
                        reference: list[tuple[float, float]],
                        bias_min: float,
                        bias_max: float) -> float:
    low = min(bias_min, bias_max)
    high = max(bias_min, bias_max)
    max_error = 0.0
    saw_point = False
    for bias, ref_current in reference:
        if bias < low or bias > high or ref_current == 0.0:
            continue
        cand_log = _interpolate_log_abs_current(candidate, bias)
        if cand_log is None:
            continue
        error = abs(cand_log - math.log10(abs(ref_current)))
        max_error = max(max_error, error)
        saw_point = True
    if not saw_point:
        raise ValueError("no overlapping nonzero current points in requested bias window")
    return max_error


def _format_bias(value: float | None) -> str:
    return "none" if value is None else f"{value:.1f} V"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare one-volt PN2D BV current-growth knees."
    )
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-10.0)
    args = parser.parse_args()

    reference = load_curve(args.reference)
    candidate = load_curve(args.candidate)
    reference_summary = summarize_curve(reference, args.bias_min, args.bias_max)
    candidate_summary = summarize_curve(candidate, args.bias_min, args.bias_max)
    max_error = max_abs_log10_error(candidate, reference, args.bias_min, args.bias_max)

    print("PN2D BV knee-shape summary")
    print(f"reference: {args.reference}")
    print(f"candidate: {args.candidate}")
    print(f"bias_window: {args.bias_min:g} V to {args.bias_max:g} V")
    print(
        "Sentaurus first 1V growth ratio > 1.5: "
        f"{_format_bias(reference_summary.first_growth_over_1p5)}"
    )
    print(
        "Sentaurus first 1V growth ratio > 2.0: "
        f"{_format_bias(reference_summary.first_growth_over_2p0)}"
    )
    print(
        "Vela first 1V growth ratio > 1.5: "
        f"{_format_bias(candidate_summary.first_growth_over_1p5)}"
    )
    print(
        "Vela first 1V growth ratio > 2.0: "
        f"{_format_bias(candidate_summary.first_growth_over_2p0)}"
    )
    print(f"max_abs_log10_current_error: {max_error:.6g} decades")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
