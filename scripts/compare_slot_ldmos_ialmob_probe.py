#!/usr/bin/env python3
"""Compare production IALMob point evaluations with native TDR mobility exports.

Reports differences, without treating averaged nodal output as an exact model oracle.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def compare(directory: Path) -> dict:
    def indexed(name: str) -> dict:
        rows = json.loads((directory / name).read_text(encoding="utf-8"))
        result = {row["node_id"]: row for row in rows}
        if len(rows) != len(result) or not result:
            raise ValueError(f"Empty or duplicate node support in {name}")
        return result

    points, reference = indexed("points.json"), indexed("reference.json")
    if points.keys() != reference.keys():
        raise ValueError("Reference and point supports differ")
    with (directory / "mobility.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    seen, grouped = set(), {"electron": [], "hole": []}
    for row in rows:
        node, carrier = int(row["node_id"]), row["carrier"]
        if node not in points or carrier not in grouped or (node, carrier) in seen:
            raise ValueError("Unexpected or duplicate mobility support")
        seen.add((node, carrier))
        local = float(row["mobility_cm2_V_s"])
        target = reference[node][carrier + "_mobility_cm2_V_s"]
        if not all(math.isfinite(x) and x > 0 for x in (local, target)):
            raise ValueError("Comparison requires positive finite mobilities")
        grouped[carrier].append({"node_id": node, "relative_error": local / target - 1,
            "local_cm2_V_s": local, "reference_cm2_V_s": target,
            "distance_nm": points[node]["interface_distance_m"] * 1e9})
    if len(seen) != 2 * len(points):
        raise ValueError("Missing mobility rows")
    result = {"schema": "slot_ldmos_ialmob_point_comparison.v1",
        "qualification": "diagnostic only; nodal averaging is not a pointwise formula oracle"}
    for carrier, values in grouped.items():
        result[carrier] = {}
        for name, lower, upper in [("all", -1, math.inf), ("within_30_nm", -1, 30),
                                   ("beyond_100_nm", 100, math.inf)]:
            selected = [v for v in values if lower < v["distance_nm"] <= upper]
            errors = sorted(abs(v["relative_error"]) for v in selected)
            result[carrier][name] = {"count": len(errors)}
            if errors:
                result[carrier][name].update({"median_abs_relative_error": statistics.median(errors),
                    "p95_abs_relative_error": errors[int(.95 * (len(errors) - 1))],
                    "max_abs_relative_error": max(errors),
                    "worst": sorted(selected, key=lambda x: abs(x["relative_error"]), reverse=True)[:10]})
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    report = compare(args.directory)
    output = args.directory / "comparison.json"
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(output)
    for carrier in ("electron", "hole"):
        for name, group in report[carrier].items():
            print(carrier, name, {k: v for k, v in group.items() if k != "worst"})


if __name__ == "__main__":
    main()
