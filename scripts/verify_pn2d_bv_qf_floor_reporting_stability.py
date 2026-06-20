#!/usr/bin/env python3
"""Verify default-vs-QF-floor restart reporting stability for one PN2D BV point."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any


def main() -> int:
    args = parse_args()
    errors: list[str] = []
    restart = read_json(args.restart_summary)

    baseline_current = read_optional_float(restart, "baseline_total_current_A_per_um")
    restart_current = read_optional_float(restart, "restart_drop_total_current_A_per_um")
    reference_current = read_optional_float(restart, "reference_current_A")
    if not math.isfinite(baseline_current):
        errors.append("restart summary is missing finite baseline_total_current_A_per_um")
    if not math.isfinite(restart_current):
        errors.append("restart summary is missing finite restart_drop_total_current_A_per_um")
    if not math.isfinite(reference_current):
        errors.append("restart summary is missing finite reference_current_A")

    default_current = read_iv_current(args.default_iv_csv, args.bias, errors, "default")
    qf_current = read_iv_current(args.qf_floor_iv_csv, args.bias, errors, "qf_floor")
    default_edges = contact_edges(args.default_contact_edge_csv, args.bias, args.contact)
    qf_edges = contact_edges(args.qf_floor_contact_edge_csv, args.bias, args.contact)

    if not default_edges:
        errors.append(f"default contact edge CSV has no {args.contact} rows at {args.bias:g} V")
    if not qf_edges:
        errors.append(f"qf-floor contact edge CSV has no {args.contact} rows at {args.bias:g} V")

    default_override_edges = count_override_edges(default_edges)
    qf_override_edges = count_override_edges(qf_edges)
    if default_override_edges != 0:
        errors.append(f"default probe has {default_override_edges} qf-floor override edges")
    if args.expected_override_edges is not None and qf_override_edges != args.expected_override_edges:
        errors.append(
            f"qf-floor override edge count {qf_override_edges} does not match expected "
            f"{args.expected_override_edges}"
        )
    elif args.expected_override_edges is None and qf_override_edges == 0:
        errors.append("qf-floor probe has no override edges")

    default_edge_sum = sum_edge_current(default_edges)
    qf_edge_sum = sum_edge_current(qf_edges)
    default_edge_delta = abs(default_edge_sum - default_current)
    qf_edge_delta = abs(qf_edge_sum - qf_current)
    if math.isfinite(default_edge_delta) and default_edge_delta > args.current_abs_tol:
        errors.append(
            f"default edge sum delta {default_edge_delta:.6g} A/um exceeds "
            f"{args.current_abs_tol:.6g}"
        )
    if math.isfinite(qf_edge_delta) and qf_edge_delta > args.current_abs_tol:
        errors.append(
            f"qf-floor edge sum delta {qf_edge_delta:.6g} A/um exceeds "
            f"{args.current_abs_tol:.6g}"
        )

    default_delta_baseline = abs(default_current - baseline_current)
    qf_delta_restart = abs(qf_current - restart_current)
    if math.isfinite(default_delta_baseline) and default_delta_baseline > args.current_abs_tol:
        errors.append(
            f"default current delta vs baseline {default_delta_baseline:.6g} A/um "
            f"exceeds {args.current_abs_tol:.6g}"
        )
    if math.isfinite(qf_delta_restart) and qf_delta_restart > args.current_abs_tol:
        errors.append(
            f"qf-floor current delta vs restart {qf_delta_restart:.6g} A/um "
            f"exceeds {args.current_abs_tol:.6g}"
        )

    default_log_error = log10_abs_error(default_current, reference_current)
    qf_log_error = log10_abs_error(qf_current, reference_current)
    if args.require_qf_improves and math.isfinite(default_log_error) and math.isfinite(qf_log_error):
        if qf_log_error >= default_log_error:
            errors.append(
                "qf-floor current does not improve log10 error: "
                f"default={default_log_error:.6g}, qf_floor={qf_log_error:.6g}"
            )

    summary = {
        "passed": not errors,
        "bias_V": args.bias,
        "contact": args.contact,
        "default_current_A_per_um": default_current,
        "qf_floor_current_A_per_um": qf_current,
        "baseline_total_current_A_per_um": baseline_current,
        "restart_drop_total_current_A_per_um": restart_current,
        "reference_current_A": reference_current,
        "default_abs_log10_error": default_log_error,
        "qf_floor_abs_log10_error": qf_log_error,
        "default_delta_vs_baseline_current_A_per_um": default_delta_baseline,
        "qf_floor_delta_vs_restart_current_A_per_um": qf_delta_restart,
        "default_edge_sum_A_per_um": default_edge_sum,
        "qf_floor_edge_sum_A_per_um": qf_edge_sum,
        "default_edge_sum_delta_A_per_um": default_edge_delta,
        "qf_floor_edge_sum_delta_A_per_um": qf_edge_delta,
        "default_override_edge_count": default_override_edges,
        "qf_floor_override_edge_count": qf_override_edges,
        "errors": errors,
    }
    if args.out_json is not None:
        args.out_json.parent.mkdir(parents=True, exist_ok=True)
        args.out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1
    print(
        "qf-floor reporting stability verified: "
        f"bias={args.bias:g} V contact={args.contact} "
        f"default_log10_error={default_log_error:.6g} "
        f"qf_floor_log10_error={qf_log_error:.6g}"
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--default-iv-csv", type=Path, required=True)
    parser.add_argument("--default-contact-edge-csv", type=Path, required=True)
    parser.add_argument("--qf-floor-iv-csv", type=Path, required=True)
    parser.add_argument("--qf-floor-contact-edge-csv", type=Path, required=True)
    parser.add_argument("--restart-summary", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--expected-override-edges", type=int)
    parser.add_argument("--current-abs-tol", type=float, default=1.0e-28)
    parser.add_argument("--require-qf-improves", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--out-json", type=Path)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_optional_float(data: dict[str, Any], key: str) -> float:
    try:
        value = float(data.get(key, math.nan))
    except (TypeError, ValueError):
        return math.nan
    return value if math.isfinite(value) else math.nan


def optional_float(raw: object) -> float | None:
    if raw is None or raw == "":
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def bias_matches(raw: str, bias: float) -> bool:
    try:
        return abs(float(raw) - bias) <= 1.0e-9
    except ValueError:
        return False


def parse_boolish(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def read_iv_current(path: Path, bias: float, errors: list[str], label: str) -> float:
    rows = read_csv(path)
    for row in rows:
        if not bias_matches(row.get("bias_V", ""), bias):
            continue
        converged = row.get("converged", "1")
        if converged != "" and not parse_boolish(converged):
            errors.append(f"{label} IV row at {bias:g} V is not converged")
        value = optional_float(row.get("current_total_A_per_um"))
        if value is None:
            raw_total = optional_float(row.get("current_total"))
            value = None if raw_total is None else raw_total / 1.0e6
        if value is None:
            errors.append(f"{label} IV row has no finite current")
            return math.nan
        return value
    errors.append(f"{label} IV CSV has no row at bias {bias:g} V")
    return math.nan


def contact_edges(path: Path, bias: float, contact: str) -> list[dict[str, str]]:
    return [
        row for row in read_csv(path)
        if row.get("current_contact", row.get("contact", "")) == contact
        and bias_matches(row.get("bias_V", ""), bias)
    ]


def count_override_edges(rows: list[dict[str, str]]) -> int:
    return sum(1 for row in rows if parse_boolish(row.get("hole_qf_drop_override_applied", "0")))


def sum_edge_current(rows: list[dict[str, str]]) -> float:
    total = 0.0
    for row in rows:
        value = optional_float(row.get("current_total_A_per_um"))
        if value is None:
            value_per_m = optional_float(row.get("current_total"))
            value = 0.0 if value_per_m is None else value_per_m / 1.0e6
        total += value
    return total


def log10_abs_error(candidate: float, reference: float) -> float:
    if (
        not math.isfinite(candidate)
        or not math.isfinite(reference)
        or candidate == 0.0
        or reference == 0.0
    ):
        return math.nan
    return abs(math.log10(abs(candidate) / abs(reference)))


if __name__ == "__main__":
    raise SystemExit(main())
