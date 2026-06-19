#!/usr/bin/env python3
"""Compare C++ SG avalanche edge-source dumps with Python reconstruction CSVs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cpp-csv", type=Path, required=True)
    parser.add_argument("--python-csv", type=Path, required=True)
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--node", action="append", type=int, default=[])
    parser.add_argument("--bias-tol", type=float, default=1.0e-9)
    parser.add_argument("--top", type=int, default=20)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in {"", None} else 0.0


def as_int(row: dict[str, str], key: str) -> int:
    return int(float(row[key]))


def contact_class(edge_class: str) -> bool:
    return edge_class in {"contact_edge", "contact_boundary", "contact_adjacent_interior"}


def source_key(row: dict[str, str]) -> float:
    return abs(as_float(row, "edge_source_integral" if "edge_source_integral" in row else "source_integral"))


def log10_ratio(numerator: float, denominator: float) -> float | None:
    if numerator <= 0.0 or denominator <= 0.0:
        return None
    return math.log10(numerator / denominator)


def summarize_cpp(rows: list[dict[str, str]],
                  bias: float,
                  bias_tol: float,
                  nodes: list[int],
                  top: int) -> dict[str, Any]:
    selected = [
        row for row in rows
        if abs(as_float(row, "bias_V") - bias) <= bias_tol
    ]
    selected.sort(key=source_key, reverse=True)
    total = sum(as_float(row, "edge_source_integral") for row in selected)
    contact = sum(
        as_float(row, "edge_source_integral")
        for row in selected
        if contact_class(row.get("edge_class", ""))
    )
    interior_bulk = sum(
        as_float(row, "edge_source_integral")
        for row in selected
        if row.get("edge_class", "") == "interior_bulk"
    )
    node_source = {str(node): 0.0 for node in nodes}
    node_set = set(nodes)
    for row in selected:
        node0 = as_int(row, "node0")
        node1 = as_int(row, "node1")
        if node0 in node_set:
            node_source[str(node0)] += as_float(row, "node0_source_integral")
        if node1 in node_set:
            node_source[str(node1)] += as_float(row, "node1_source_integral")
    return {
        "row_count": len(selected),
        "top_edge_ids": [as_int(row, "edge_id") for row in selected[:top]],
        "total_source_integral": total,
        "contact_edge_source_fraction": contact / total if total > 0.0 else 0.0,
        "interior_bulk_source_fraction": interior_bulk / total if total > 0.0 else 0.0,
        "node_source_integrals": node_source,
    }


def summarize_python(rows: list[dict[str, str]],
                     nodes: list[int],
                     top: int) -> dict[str, Any]:
    selected = list(rows)
    selected.sort(key=source_key, reverse=True)
    total = sum(as_float(row, "source_integral") for row in selected)
    contact = sum(
        as_float(row, "source_integral")
        for row in selected
        if contact_class(row.get("edge_class", ""))
    )
    interior_bulk = sum(
        as_float(row, "source_integral")
        for row in selected
        if row.get("edge_class", "") == "interior_bulk"
    )
    node_source = {str(node): 0.0 for node in nodes}
    node_set = set(nodes)
    for row in selected:
        source = as_float(row, "source_integral")
        half = 0.5 * source
        node0 = as_int(row, "node0")
        node1 = as_int(row, "node1")
        if node0 in node_set:
            node_source[str(node0)] += half
        if node1 in node_set:
            node_source[str(node1)] += half
    return {
        "row_count": len(selected),
        "top_edge_ids": [as_int(row, "edge_id") for row in selected[:top]],
        "total_source_integral": total,
        "contact_edge_source_fraction": contact / total if total > 0.0 else 0.0,
        "interior_bulk_source_fraction": interior_bulk / total if total > 0.0 else 0.0,
        "node_source_integrals": node_source,
    }


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main() -> None:
    args = parse_args()
    cpp = summarize_cpp(read_rows(args.cpp_csv), args.bias, args.bias_tol, args.node, args.top)
    python = summarize_python(read_rows(args.python_csv), args.node, args.top)
    summary = {
        "schema": "pn2d.sg_edge_source_dump_compare.v1",
        "bias_V": args.bias,
        "cpp_csv": str(args.cpp_csv),
        "python_csv": str(args.python_csv),
        "cpp": cpp,
        "python": python,
        "comparison": {
            "total_source_log10_cpp_over_python": log10_ratio(
                cpp["total_source_integral"],
                python["total_source_integral"],
            ),
            "top_edge_order_matches": cpp["top_edge_ids"] == python["top_edge_ids"],
            "contact_edge_fraction_delta": (
                cpp["contact_edge_source_fraction"] - python["contact_edge_source_fraction"]
            ),
            "interior_bulk_fraction_delta": (
                cpp["interior_bulk_source_fraction"] - python["interior_bulk_source_fraction"]
            ),
        },
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(clean_json(summary), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
