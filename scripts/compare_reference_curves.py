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


def values(rows: list[dict[str, str]], column: str) -> list[float]:
    out = []
    for row in rows:
        text = row.get(column, "")
        if text == "":
            continue
        out.append(float(text))
    return out


def trend(values_: list[float]) -> str:
    clean = [value for value in values_ if math.isfinite(value)]
    if len(clean) < 2:
        return "insufficient"
    diffs = [b - a for a, b in zip(clean, clean[1:])]
    eps = max(max(abs(value) for value in clean), 1.0) * 1.0e-12
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
                   candidate_rows: list[dict[str, str]]) -> dict[str, Any]:
    ref_col = find_column(reference_rows, SERIES[kind])
    cand_col = find_column(candidate_rows, SERIES[kind])
    if ref_col is None or cand_col is None:
        return {
            "available": False,
            "reference_column": ref_col,
            "candidate_column": cand_col,
            "trend_match": False,
        }

    ref_values = values(reference_rows, ref_col)
    cand_values = values(candidate_rows, cand_col)
    ref_trend = trend(ref_values)
    cand_trend = trend(cand_values)
    return {
        "available": True,
        "reference_column": ref_col,
        "candidate_column": cand_col,
        "points_compared": min(len(ref_values), len(cand_values)),
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
        "| Curve | Available | Trend match | Reference trend | Candidate trend | Max relative error | Orders of magnitude |",
        "| --- | --- | --- | --- | --- | ---: | ---: |",
    ]
    for key in ("iv", "cv", "bv"):
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
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    reference_rows = read_csv(args.reference)
    candidate_rows = read_csv(args.candidate)
    report = {
        "reference": str(args.reference),
        "candidate": str(args.candidate),
        "iv": compare_series("iv", reference_rows, candidate_rows),
        "cv": compare_series("cv", reference_rows, candidate_rows),
        "bv": compare_series("bv", reference_rows, candidate_rows),
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n")
    args.output_md.write_text(markdown_report(report))
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
