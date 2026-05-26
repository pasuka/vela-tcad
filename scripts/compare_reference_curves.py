#!/usr/bin/env python3
"""Compare Vela curves with neutral reference TCAD CSV curves."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


SERIES = {
    "iv": ["current_total", "current_total_A_per_um", "Id", "I"],
    "cv": ["capacitance_F_per_m", "capacitance_F_per_um", "capacitance", "C"],
    "bv": ["max_electric_field_V_per_cm", "max_electric_field_V_per_m", "max_field"],
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def find_column(rows: list[dict[str, str]], candidates: list[str]) -> str | None:
    if not rows:
        return None
    columns = set(rows[0])
    for candidate in candidates:
        if candidate in columns:
            return candidate
    return None


def resolve_column(rows: list[dict[str, str]],
                   candidates: list[str],
                   explicit: str | None) -> str | None:
    if explicit is None:
        return find_column(rows, candidates)
    if not rows:
        return None
    return explicit if explicit in rows[0] else None


def values(rows: list[dict[str, str]], column: str) -> list[float]:
    out = []
    for row in rows:
        text = row.get(column, "")
        if text == "":
            continue
        out.append(float(text))
    return out


def finite_pairs(rows: list[dict[str, str]], value_column: str, scale: float) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for row in rows:
        if "bias_V" not in row:
            continue
        bias_text = row.get("bias_V", "")
        value_text = row.get(value_column, "")
        if bias_text == "" or value_text == "":
            continue
        bias = float(bias_text)
        value = float(value_text) * scale
        if math.isfinite(bias) and math.isfinite(value):
            pairs.append((bias, value))
    return sorted(pairs)


def interpolate_at(pairs: list[tuple[float, float]], bias: float) -> float | None:
    if not pairs:
        return None
    if bias < pairs[0][0] or bias > pairs[-1][0]:
        return None
    for existing_bias, value in pairs:
        if abs(existing_bias - bias) <= max(abs(bias), 1.0) * 1.0e-12:
            return value
    for (b0, v0), (b1, v1) in zip(pairs, pairs[1:]):
        if b0 <= bias <= b1 and b1 != b0:
            t = (bias - b0) / (b1 - b0)
            return v0 + t * (v1 - v0)
    return None


def aligned_values(reference_rows: list[dict[str, str]],
                   candidate_rows: list[dict[str, str]],
                   ref_col: str,
                   cand_col: str,
                   candidate_scale: float,
                   bias_min: float | None,
                   bias_max: float | None) -> tuple[list[float], list[float], list[float]]:
    ref_pairs = finite_pairs(reference_rows, ref_col, 1.0)
    cand_pairs = finite_pairs(candidate_rows, cand_col, candidate_scale)
    if not ref_pairs or not cand_pairs:
        return [], values(reference_rows, ref_col), [
            value * candidate_scale for value in values(candidate_rows, cand_col)
        ]

    biases: list[float] = []
    ref_values: list[float] = []
    cand_values: list[float] = []
    for bias, ref_value in ref_pairs:
        if bias_min is not None and bias < bias_min:
            continue
        if bias_max is not None and bias > bias_max:
            continue
        cand_value = interpolate_at(cand_pairs, bias)
        if cand_value is None:
            continue
        biases.append(bias)
        ref_values.append(ref_value)
        cand_values.append(cand_value)
    return biases, ref_values, cand_values


def trend(values_: list[float]) -> str:
    clean = [value for value in values_ if math.isfinite(value)]
    if len(clean) < 2:
        return "insufficient"
    diffs = [b - a for a, b in zip(clean, clean[1:])]
    eps = max(max(abs(value) for value in clean) * 1.0e-12, 1.0e-300)
    nonnegative = all(diff >= -eps for diff in diffs)
    nonpositive = all(diff <= eps for diff in diffs)
    if nonnegative and not nonpositive:
        return "increasing"
    if nonpositive and not nonnegative:
        return "decreasing"
    if nonnegative and nonpositive:
        return "flat"
    return "mixed"


def orders_of_magnitude(reference: list[float], candidate: list[float]) -> float | None:
    pairs = [
        (abs(ref), abs(cand))
        for ref, cand in zip(reference, candidate)
        if ref != 0.0 and cand != 0.0 and math.isfinite(ref) and math.isfinite(cand)
    ]
    if not pairs:
        return None
    return max(abs(math.log10(cand / ref)) for ref, cand in pairs)


def relative_error(reference: list[float], candidate: list[float]) -> float | None:
    errors = []
    for ref, cand in zip(reference, candidate):
        if not math.isfinite(ref) or not math.isfinite(cand):
            continue
        scale = max(abs(ref), 1.0e-300)
        errors.append(abs(cand - ref) / scale)
    return max(errors) if errors else None


def compare_series(kind: str,
                   reference_rows: list[dict[str, str]],
                   candidate_rows: list[dict[str, str]],
                   candidate_scale: float,
                   bias_min: float | None,
                   bias_max: float | None,
                   reference_column: str | None = None,
                   candidate_column: str | None = None) -> dict[str, Any]:
    ref_col = resolve_column(reference_rows, SERIES[kind], reference_column)
    cand_col = resolve_column(candidate_rows, SERIES[kind], candidate_column)
    if ref_col is None or cand_col is None:
        return {
            "available": False,
            "reference_column": ref_col,
            "candidate_column": cand_col,
            "trend_match": False,
        }

    biases, ref_values, cand_values = aligned_values(
        reference_rows,
        candidate_rows,
        ref_col,
        cand_col,
        candidate_scale,
        bias_min,
        bias_max,
    )
    ref_trend = trend(ref_values)
    cand_trend = trend(cand_values)
    return {
        "available": True,
        "reference_column": ref_col,
        "candidate_column": cand_col,
        "points_compared": min(len(ref_values), len(cand_values)),
        "reference_bias_range": [biases[0], biases[-1]] if biases else None,
        "candidate_scale": candidate_scale,
        "reference_trend": ref_trend,
        "candidate_trend": cand_trend,
        "trend_match": ref_trend == cand_trend,
        "orders_of_magnitude": orders_of_magnitude(ref_values, cand_values),
        "max_relative_error": relative_error(ref_values, cand_values),
    }


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# Reference TCAD Curve Comparison",
        "",
        f"Status: {report.get('status', 'unknown')}",
        "",
        "| Curve | Available | Trend match | Reference trend | Candidate trend | Max relative error | Orders of magnitude |",
        "| --- | --- | --- | --- | --- | ---: | ---: |",
    ]
    for key in report.get("checked_kinds", ("iv", "cv", "bv")):
        item = report[key]
        lines.append(
            "| {key} | {available} | {trend_match} | {ref_trend} | {cand_trend} | {rel} | {oom} |".format(
                key=key.upper(),
                available=item.get("available", False),
                trend_match=item.get("trend_match", False),
                ref_trend=item.get("reference_trend", ""),
                cand_trend=item.get("candidate_trend", ""),
                rel=format_optional(item.get("max_relative_error")),
                oom=format_optional(item.get("orders_of_magnitude")),
            )
        )
    if report.get("failures"):
        lines.append("")
        lines.append("Failures:")
        for failure in report["failures"]:
            lines.append(f"- {failure}")
    lines.append("")
    lines.append("Inputs are explicit CSV/text exports; no proprietary binary formats are parsed.")
    return "\n".join(lines) + "\n"


def format_optional(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    parser.add_argument(
        "--kind",
        choices=["all", "iv", "cv", "bv"],
        default="all",
        help="Curve kind to check. The default preserves the historical all-curve report.",
    )
    parser.add_argument("--require-trend-match", action="store_true")
    parser.add_argument("--min-points", type=int, default=0)
    parser.add_argument("--max-relative-error", type=float)
    parser.add_argument("--max-orders-of-magnitude", type=float)
    parser.add_argument("--candidate-scale", type=float, default=1.0)
    parser.add_argument("--reference-column")
    parser.add_argument("--candidate-column")
    parser.add_argument("--bias-min", type=float)
    parser.add_argument("--bias-max", type=float)
    return parser.parse_args()


def checked_kinds(kind: str) -> list[str]:
    if kind == "all":
        return ["iv", "cv", "bv"]
    return [kind]


def exceeds(value: Any, limit: float | None) -> bool:
    if limit is None or value is None:
        return False
    return float(value) > limit + max(abs(limit), 1.0) * 1.0e-12


def evaluate_failures(report: dict[str, Any], args: argparse.Namespace) -> list[str]:
    failures: list[str] = []
    gate_active = (
        args.kind != "all"
        or args.require_trend_match
        or args.min_points
        or args.max_relative_error is not None
        or args.max_orders_of_magnitude is not None
    )
    for kind in report["checked_kinds"]:
        item = report[kind]
        label = kind.upper()
        if not item.get("available", False):
            if gate_active:
                failures.append(f"{label}: required curve columns are not available")
            continue
        points = int(item.get("points_compared", 0))
        if args.min_points and points < args.min_points:
            failures.append(f"{label}: compared {points} point(s), fewer than required {args.min_points}")
        if args.require_trend_match and not item.get("trend_match", False):
            failures.append(
                f"{label}: trend mismatch "
                f"({item.get('reference_trend')} vs {item.get('candidate_trend')})"
            )
        if exceeds(item.get("max_relative_error"), args.max_relative_error):
            failures.append(
                f"{label}: max relative error {item.get('max_relative_error'):.6g} "
                f"exceeds {args.max_relative_error:.6g}"
            )
        if exceeds(item.get("orders_of_magnitude"), args.max_orders_of_magnitude):
            failures.append(
                f"{label}: orders-of-magnitude delta {item.get('orders_of_magnitude'):.6g} "
                f"exceeds {args.max_orders_of_magnitude:.6g}"
            )
    return failures


def main() -> int:
    args = parse_args()
    reference_rows = read_csv(args.reference)
    candidate_rows = read_csv(args.candidate)
    kinds = checked_kinds(args.kind)
    report = {
        "reference": str(args.reference),
        "candidate": str(args.candidate),
        "checked_kinds": kinds,
        "iv": compare_series(
            "iv", reference_rows, candidate_rows,
            args.candidate_scale, args.bias_min, args.bias_max,
            args.reference_column, args.candidate_column),
        "cv": compare_series(
            "cv", reference_rows, candidate_rows,
            args.candidate_scale, args.bias_min, args.bias_max,
            args.reference_column, args.candidate_column),
        "bv": compare_series(
            "bv", reference_rows, candidate_rows,
            args.candidate_scale, args.bias_min, args.bias_max,
            args.reference_column, args.candidate_column),
    }
    report["failures"] = evaluate_failures(report, args)
    report["status"] = "fail" if report["failures"] else "pass"

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    args.output_md.write_text(markdown_report(report))
    print(json.dumps(report, indent=2))
    return 1 if report["failures"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
