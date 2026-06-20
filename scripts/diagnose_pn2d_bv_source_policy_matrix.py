#!/usr/bin/env python3
"""Evaluate carrier-row avalanche source policies from carrier-term probes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import os
from itertools import product
from pathlib import Path
from typing import Any


POLICIES = [
    "combined",
    "double_combined",
    "combined_plus_electron",
    "combined_plus_hole",
    "electron_only",
    "hole_only",
    "zero",
]

SUMMARY_FIELDS = [
    "scope",
    "variant",
    "bias_V",
    "electron_policy",
    "hole_policy",
    "count",
    "electron_l2",
    "hole_l2",
    "combined_carrier_l2",
    "electron_signed_median_over_combined",
    "hole_signed_median_over_combined",
    "electron_abs_median_over_combined",
    "hole_abs_median_over_combined",
    "electron_max_abs_over_combined",
    "hole_max_abs_over_combined",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exact-node-csv", type=Path, required=True)
    parser.add_argument("--carrier-term-csv", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--variant", default="vela_psi_sentaurus_qf")
    parser.add_argument("--bias", type=float, default=None)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def finite_float(raw: Any, default: float = 0.0) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def matches_bias(row: dict[str, str], bias: float | None) -> bool:
    if bias is None:
        return True
    if "bias_V" not in row or row.get("bias_V") in (None, ""):
        return True
    return abs(finite_float(row["bias_V"]) - bias) <= 1.0e-9


def load_exact_rows(path: Path, variant: str, bias: float | None) -> list[dict[str, str]]:
    return [
        row for row in read_csv(path)
        if row.get("variant") == variant and matches_bias(row, bias)
    ]


def merge_source_components(
    exact_rows: list[dict[str, str]],
    carrier_rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    by_node = {row.get("node_id", ""): row for row in carrier_rows}
    component_fields = [
        "impact_electron_source",
        "impact_hole_source",
        "impact_combined_source",
    ]
    merged: list[dict[str, str]] = []
    for row in exact_rows:
        item = dict(row)
        carrier = by_node.get(item.get("node_id", ""))
        if carrier is not None:
            for field in component_fields:
                if item.get(field) in (None, "") and carrier.get(field) not in (None, ""):
                    item[field] = carrier[field]
        merged.append(item)
    return merged


def combined_source(row: dict[str, str]) -> float:
    value = finite_float(row.get("impact_combined_source"))
    if value != 0.0:
        return value
    impact = finite_float(row.get("electron_impact"))
    return abs(impact)


def policy_source(row: dict[str, str], policy: str) -> float:
    electron = finite_float(row.get("impact_electron_source"))
    hole = finite_float(row.get("impact_hole_source"))
    combined = combined_source(row)
    if policy == "combined":
        return combined
    if policy == "double_combined":
        return 2.0 * combined
    if policy == "combined_plus_electron":
        return combined + electron
    if policy == "combined_plus_hole":
        return combined + hole
    if policy == "electron_only":
        return electron
    if policy == "hole_only":
        return hole
    if policy == "zero":
        return 0.0
    raise ValueError(f"unknown source policy: {policy}")


def nonimpact_terms(row: dict[str, str], carrier: str) -> float:
    return (
        finite_float(row.get(f"{carrier}_flux"))
        + finite_float(row.get(f"{carrier}_recombination"))
        + finite_float(row.get(f"{carrier}_gauge"))
        + finite_float(row.get(f"{carrier}_boundary"))
    )


def carrier_residual(row: dict[str, str], carrier: str, policy: str) -> float:
    return nonimpact_terms(row, carrier) - policy_source(row, policy)


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def evaluate_scope(
    *,
    rows: list[dict[str, str]],
    scope: str,
    variant: str,
    bias: float | None,
    electron_policy: str,
    hole_policy: str,
) -> dict[str, Any]:
    electron_residuals: list[float] = []
    hole_residuals: list[float] = []
    electron_normed: list[float] = []
    hole_normed: list[float] = []
    for row in rows:
        electron_residual = carrier_residual(row, "electron", electron_policy)
        hole_residual = carrier_residual(row, "hole", hole_policy)
        electron_residuals.append(electron_residual)
        hole_residuals.append(hole_residual)
        denom = abs(combined_source(row))
        if denom > 0.0:
            electron_normed.append(electron_residual / denom)
            hole_normed.append(hole_residual / denom)

    electron_l2 = math.sqrt(sum(value * value for value in electron_residuals))
    hole_l2 = math.sqrt(sum(value * value for value in hole_residuals))
    return {
        "scope": scope,
        "variant": variant,
        "bias_V": bias,
        "electron_policy": electron_policy,
        "hole_policy": hole_policy,
        "count": len(rows),
        "electron_l2": electron_l2,
        "hole_l2": hole_l2,
        "combined_carrier_l2": math.sqrt(electron_l2 * electron_l2 + hole_l2 * hole_l2),
        "electron_signed_median_over_combined": median(electron_normed),
        "hole_signed_median_over_combined": median(hole_normed),
        "electron_abs_median_over_combined": median([abs(value) for value in electron_normed]),
        "hole_abs_median_over_combined": median([abs(value) for value in hole_normed]),
        "electron_max_abs_over_combined": max([abs(value) for value in electron_normed], default=None),
        "hole_max_abs_over_combined": max([abs(value) for value in hole_normed], default=None),
    }


def policy_matrix_for_scope(
    *,
    rows: list[dict[str, str]],
    scope: str,
    variant: str,
    bias: float | None,
) -> list[dict[str, Any]]:
    return [
        evaluate_scope(
            rows=rows,
            scope=scope,
            variant=variant,
            bias=bias,
            electron_policy=electron_policy,
            hole_policy=hole_policy,
        )
        for electron_policy, hole_policy in product(POLICIES, POLICIES)
    ]


def support_scopes(rows: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        support_class = row.get("support_class", "")
        if not support_class:
            continue
        result.setdefault(f"support:{support_class}", []).append(row)
    return result


def best_row(rows: list[dict[str, Any]]) -> dict[str, Any]:
    if not rows:
        return {}
    return min(rows, key=lambda row: (
        finite_float(row.get("combined_carrier_l2"), math.inf),
        str(row.get("electron_policy")),
        str(row.get("hole_policy")),
    ))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with open_path(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def write_text(path: Path, text: str) -> None:
    with open_path(path, "w", encoding="utf-8") as handle:
        handle.write(text)


def main() -> int:
    args = parse_args()
    exact_rows = load_exact_rows(args.exact_node_csv, args.variant, args.bias)
    if not exact_rows:
        raise SystemExit("no exact rows matched the requested variant/bias")
    carrier_rows = read_csv(args.carrier_term_csv) if args.carrier_term_csv is not None else []
    if carrier_rows:
        exact_rows = merge_source_components(exact_rows, carrier_rows)

    all_summary_rows: list[dict[str, Any]] = []
    scopes = support_scopes(exact_rows)
    if carrier_rows:
        scopes["global"] = carrier_rows

    best_by_scope: dict[str, Any] = {}
    for scope, rows in sorted(scopes.items()):
        matrix = policy_matrix_for_scope(
            rows=rows,
            scope=scope,
            variant=args.variant,
            bias=args.bias,
        )
        all_summary_rows.extend(matrix)
        best_by_scope[scope] = best_row(matrix)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "source_policy_matrix_summary.csv", all_summary_rows)
    payload = clean_json({
        "variant": args.variant,
        "bias_V": args.bias,
        "row_count": len(all_summary_rows),
        "policies": POLICIES,
        "best_by_scope": best_by_scope,
    })
    write_text(
        args.out_dir / "source_policy_matrix_summary.json",
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
