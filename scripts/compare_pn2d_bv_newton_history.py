#!/usr/bin/env python3
"""Compare PN2D BV Newton-history tails from localization runs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def last_float(rows: list[dict[str, str]], column: str) -> float | None:
    for row in reversed(rows):
        value = row.get(column, "")
        if value:
            return float(value)
    return None


def summarize(path: Path) -> dict[str, object]:
    rows = read_rows(path)
    return {
        "path": str(path),
        "rows": len(rows),
        "last_bias_V": last_float(rows, "bias_V"),
        "last_residual_norm": last_float(rows, "residual_norm"),
        "last_raw_step_norm": last_float(rows, "raw_step_norm"),
        "last_block_psi": last_float(rows, "block_psi"),
        "last_block_phin": last_float(rows, "block_phin"),
        "last_block_phip": last_float(rows, "block_phip"),
        "max_raw_step_norm": max(
            (float(row["raw_step_norm"]) for row in rows if row.get("raw_step_norm")),
            default=None,
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", required=True, type=Path)
    parser.add_argument("--right", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    left = summarize(args.left)
    right = summarize(args.right)
    comparison = {
        "left": left,
        "right": right,
        "raw_step_ratio": (
            right["last_raw_step_norm"] / left["last_raw_step_norm"]
            if left["last_raw_step_norm"] and right["last_raw_step_norm"]
            else None
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(comparison, indent=2), encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
