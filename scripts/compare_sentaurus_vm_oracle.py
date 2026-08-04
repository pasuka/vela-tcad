#!/usr/bin/env python3
"""Compare Vela outputs against Sentaurus VM oracle CSV curves."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable


def load_series(path: Path, x_col: str, y_col: str) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None or x_col not in reader.fieldnames or y_col not in reader.fieldnames:
            raise ValueError(f"CSV {path} missing required columns {x_col}/{y_col}")
        pairs: list[tuple[float, float]] = []
        for row in reader:
            pairs.append((float(row[x_col]), float(row[y_col])))
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


def compare(vela_csv: Path, oracle_csv: Path, x_col: str, y_col: str, tolerance: float) -> dict:
    vela = load_series(vela_csv, x_col, y_col)
    oracle = load_series(oracle_csv, x_col, y_col)
    diff = max_abs_diff(vela, oracle)
    return {
        "vela_csv": str(vela_csv),
        "oracle_csv": str(oracle_csv),
        "x_col": x_col,
        "y_col": y_col,
        "samples": len(vela),
        "max_abs_diff": diff,
        "tolerance": tolerance,
        "pass": diff <= tolerance,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-csv", required=True, type=Path)
    parser.add_argument("--oracle-csv", required=True, type=Path)
    parser.add_argument("--x-col", default="voltage_V")
    parser.add_argument("--y-col", default="current_A_per_um")
    parser.add_argument("--tolerance", required=True, type=float)
    parser.add_argument("--report-json", type=Path)
    args = parser.parse_args()

    report = compare(args.vela_csv, args.oracle_csv, args.x_col, args.y_col, args.tolerance)
    text = json.dumps(report, indent=2)
    if args.report_json:
        args.report_json.write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

