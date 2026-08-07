#!/usr/bin/env python3
"""Generate machine-readable diffs between Sentaurus VM and Vela exports."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable


def _validate_non_negative_finite(name: str, value: float) -> None:
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"{name} must be finite and >= 0")


def load_series(path: Path, x_col: str, y_col: str) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or x_col not in reader.fieldnames or y_col not in reader.fieldnames:
            raise ValueError(f"CSV {path} missing required columns {x_col}/{y_col}")
        pairs: list[tuple[float, float]] = []
        for row in reader:
            x_value = float(row[x_col])
            y_value = float(row[y_col])
            if not math.isfinite(x_value) or not math.isfinite(y_value):
                raise ValueError(f"CSV {path} contains non-finite values")
            pairs.append((x_value, y_value))
    if not pairs:
        raise ValueError(f"CSV {path} contains no samples")
    if len(pairs) >= 2:
        first_step = pairs[1][0] - pairs[0][0]
        if first_step == 0.0:
            raise ValueError(f"CSV {path} contains duplicate x values")
        direction = 1.0 if first_step > 0.0 else -1.0
        prev_x = pairs[0][0]
        for x_value, _ in pairs[1:]:
            delta = x_value - prev_x
            if delta == 0.0:
                raise ValueError(f"CSV {path} contains duplicate x values")
            if direction > 0.0 and delta <= 0.0:
                raise ValueError(f"CSV {path} x grid changes direction")
            if direction < 0.0 and delta >= 0.0:
                raise ValueError(f"CSV {path} x grid changes direction")
            prev_x = x_value
    return pairs


def max_abs_diff(a: Iterable[tuple[float, float]], b: Iterable[tuple[float, float]]) -> float:
    a_list = list(a)
    b_list = list(b)
    if len(a_list) != len(b_list):
        raise ValueError("series length mismatch")
    if not a_list:
        return 0.0
    diffs = [abs(ay - by) for (_, ay), (_, by) in zip(a_list, b_list)]
    return max(diffs)


def max_abs_x_diff(a: Iterable[tuple[float, float]], b: Iterable[tuple[float, float]]) -> float:
    a_list = list(a)
    b_list = list(b)
    if len(a_list) != len(b_list):
        raise ValueError("series length mismatch")
    if not a_list:
        return 0.0
    diffs = [abs(ax - bx) for (ax, _), (bx, _) in zip(a_list, b_list)]
    return max(diffs)


def compare_series(
    vela_csv: Path,
    oracle_csv: Path,
    x_col: str,
    y_col: str,
    tolerance: float,
    x_tolerance: float = 0.0,
) -> dict:
    _validate_non_negative_finite("tolerance", tolerance)
    _validate_non_negative_finite("x_tolerance", x_tolerance)
    vela = load_series(vela_csv, x_col, y_col)
    oracle = load_series(oracle_csv, x_col, y_col)
    if len(vela) != len(oracle):
        raise ValueError("series length mismatch")
    vela_direction = 0 if len(vela) < 2 else (1 if vela[1][0] > vela[0][0] else -1)
    oracle_direction = 0 if len(oracle) < 2 else (1 if oracle[1][0] > oracle[0][0] else -1)
    if vela_direction != oracle_direction:
        raise ValueError("series direction mismatch")
    x_diff = max_abs_x_diff(vela, oracle)
    diff = max_abs_diff(vela, oracle)
    return {
        "vela_csv": str(vela_csv),
        "oracle_csv": str(oracle_csv),
        "x_col": x_col,
        "y_col": y_col,
        "samples": len(vela),
        "max_abs_x_diff": x_diff,
        "max_abs_diff": diff,
        "tolerance": tolerance,
        "x_tolerance": x_tolerance,
        "x_grid_pass": x_diff <= x_tolerance,
        "y_grid_pass": diff <= tolerance,
        "pass": x_diff <= x_tolerance and diff <= tolerance,
    }


def compare(
    vela_csv: Path,
    oracle_csv: Path,
    x_col: str,
    y_col: str,
    tolerance: float,
    x_tolerance: float = 0.0,
) -> dict:
    return compare_series(vela_csv, oracle_csv, x_col, y_col, tolerance, x_tolerance)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-csv", required=True, type=Path)
    parser.add_argument("--oracle-csv", required=True, type=Path)
    parser.add_argument("--x-col", default="voltage_V")
    parser.add_argument("--y-col", default="current_A_per_um")
    parser.add_argument("--tolerance", required=True, type=float)
    parser.add_argument("--x-tolerance", type=float, default=0.0)
    parser.add_argument("--report-json", type=Path)
    args = parser.parse_args()

    report = compare_series(
        args.vela_csv,
        args.oracle_csv,
        args.x_col,
        args.y_col,
        args.tolerance,
        args.x_tolerance,
    )
    text = json.dumps(report, indent=2)
    if args.report_json:
        args.report_json.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
