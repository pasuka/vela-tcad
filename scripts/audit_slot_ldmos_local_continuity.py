#!/usr/bin/env python3
"""Audit every active carrier row from a solved-equation carrier-term probe.

This is a read-only error audit, not a replacement convergence gate. Residuals,
fluxes and sources retain the probe's assembly units; ratios are dimensionless.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def finite(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"nonfinite {key} at node {row['node_id']}")
    return value


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def audit(path: Path, eps: float = 1e-3, source_threshold: float = 1e-13,
          scale_floor: float = 1e-30) -> dict:
    if not all(math.isfinite(value) for value in (eps, source_threshold, scale_floor)) or eps <= 0 or source_threshold < 0 or scale_floor <= 0:
        raise ValueError("invalid audit scales")
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    ids = [int(row["node_id"]) for row in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("duplicate node_id")
    summary = {"input": str(path.resolve()), "node_count": len(rows),
               "eps_row_diagnostic_only": eps, "source_threshold": source_threshold,
               "scale_floor": scale_floor, "carriers": {}}
    for carrier in ("electron", "hole"):
        active = []
        for row in rows:
            flag = row[f"{carrier}_continuity_active"]
            if flag not in ("0", "1"):
                raise ValueError("invalid continuity-active flag")
            values = {key: finite(row, f"{carrier}_{key}") for key in
                      ("flux", "flux_abs_sum", "recombination", "impact", "residual",
                       "gauge", "boundary", "term_sum")}
            if flag == "0":
                continue
            source = max(abs(values["recombination"]), abs(values["impact"]))
            scale = max(abs(values["flux_abs_sum"]), source, scale_floor)
            active.append({"node_id": int(row["node_id"]),
                           "x": finite(row, "x"), "y": finite(row, "y"),
                           **values, "scale": scale, "source_scale": source,
                           "ratio": abs(values["residual"]) / scale})
        if not active:
            raise ValueError(f"no active {carrier} continuity rows")
        residuals = [row["residual"] for row in active]
        source_sum = math.fsum(row["recombination"] + row["impact"] for row in active)
        absolute_source_sum = math.fsum(abs(row["recombination"] + row["impact"])
                                        for row in active)
        absolute_residual_sum = math.fsum(abs(value) for value in residuals)
        ratios = [row["ratio"] for row in active]
        summary["carriers"][carrier] = {
            "active_rows": len(active), "inactive_rows": len(rows) - len(active),
            "source_qualified_rows": sum(row["source_scale"] > 0 and
                                          row["source_scale"] >= source_threshold for row in active),
            "rows_above_eps": sum(value > eps for value in ratios),
            "ratio_p50": percentile(ratios, 0.5), "ratio_p95": percentile(ratios, 0.95),
            "ratio_max": max(ratios), "residual_signed_sum": math.fsum(residuals),
            "residual_absolute_sum": absolute_residual_sum,
            "source_signed_sum": source_sum, "source_absolute_sum": absolute_source_sum,
            "residual_l1_over_source_l1": absolute_residual_sum / max(absolute_source_sum, scale_floor),
            "max_term_residual_disagreement": max(abs(row["term_sum"] - row["residual"])
                                                    for row in active),
            "top_ratio_rows": sorted(active, key=lambda row: row["ratio"], reverse=True)[:20],
            "top_absolute_residual_rows": sorted(active, key=lambda row: abs(row["residual"]),
                                                  reverse=True)[:20],
        }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--eps-row", type=float, default=1e-3)
    parser.add_argument("--source-threshold", type=float, default=1e-13)
    args = parser.parse_args()
    result = audit(args.input, args.eps_row, args.source_threshold)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: {name: value for name, value in carrier.items()
                           if not name.startswith("top_")}
                      for key, carrier in result["carriers"].items()}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
