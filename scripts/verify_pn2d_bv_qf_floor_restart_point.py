#!/usr/bin/env python3
"""Verify one PN2D BV qf-floor restart reporting point.

This is intentionally narrow: it checks that a qf-floor-enabled Vela run at one
external Sentaurus-restart bias reports the same terminal current as the
standalone restart-drop diagnostic predicted, and that the expected contact
edges actually used the reporting override.
"""

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

    iv_rows = read_csv(args.iv_csv)
    edge_rows = read_csv(args.contact_edge_csv)
    restart_summary = read_json(args.restart_summary) if args.restart_summary else {}

    iv_row = find_bias_row(iv_rows, args.bias)
    if iv_row is None:
        errors.append(f"missing iv row at bias {args.bias:g} V")
        current = math.nan
    else:
        if not row_is_converged(iv_row):
            errors.append(f"iv row at bias {args.bias:g} V is not converged")
        current = read_current_a_per_um(iv_row, errors)

    reference_current = read_optional_float(
        restart_summary,
        "reference_current_A",
        args.reference_current,
    )
    restart_current = read_optional_float(
        restart_summary,
        "restart_drop_total_current_A_per_um",
        args.restart_current,
    )
    restart_log_error = read_optional_float(
        restart_summary,
        "abs_log10_error_restart_drop",
        math.nan,
    )

    abs_log10_error = (
        log10_abs_error(current, reference_current)
        if math.isfinite(current) and math.isfinite(reference_current)
        else math.nan
    )
    if not math.isfinite(reference_current):
        errors.append("missing reference current; pass --reference-current or --restart-summary")
    elif not math.isfinite(abs_log10_error):
        errors.append("cannot compute finite log10 current error")
    elif abs_log10_error > args.max_log10_error:
        errors.append(
            "log10 current error "
            f"{abs_log10_error:.6g} exceeds limit {args.max_log10_error:.6g}"
        )

    abs_delta_restart = (
        abs(current - restart_current)
        if math.isfinite(current) and math.isfinite(restart_current)
        else math.nan
    )
    if math.isfinite(restart_current):
        if abs_delta_restart > args.current_abs_tol:
            errors.append(
                "current delta vs restart prediction "
                f"{abs_delta_restart:.6g} A/um exceeds limit {args.current_abs_tol:.6g}"
            )
    elif args.restart_summary is not None:
        errors.append("restart summary is missing restart_drop_total_current_A_per_um")

    edge_at_bias = [row for row in edge_rows if bias_matches(row.get("bias_V", ""), args.bias)]
    contact_edges = [
        row for row in edge_at_bias
        if row.get("current_contact", row.get("contact", "")) == args.contact
    ]
    if not contact_edges:
        errors.append(f"missing contact edge rows for {args.contact} at {args.bias:g} V")

    override_column = "hole_qf_drop_override_applied"
    if edge_rows and override_column not in edge_rows[0]:
        errors.append(f"missing {override_column} column in contact edge CSV")

    override_edges = [
        row for row in contact_edges
        if parse_boolish(row.get(override_column, "0"))
    ]
    other_contact_override_edges = [
        row for row in edge_at_bias
        if row.get("current_contact", row.get("contact", "")) != args.contact
        and parse_boolish(row.get(override_column, "0"))
    ]

    if not override_edges:
        errors.append(f"no qf-floor override edges found for {args.contact} at {args.bias:g} V")
    if args.expected_override_edges is not None:
        if len(override_edges) != args.expected_override_edges:
            errors.append(
                f"override edge count {len(override_edges)} does not match expected "
                f"{args.expected_override_edges}"
            )
    if args.forbid_other_contact_overrides and other_contact_override_edges:
        errors.append(
            f"found {len(other_contact_override_edges)} qf-floor override edges on other contacts"
        )

    summary = {
        "passed": not errors,
        "bias_V": args.bias,
        "contact": args.contact,
        "current_total_A_per_um": current,
        "reference_current_A": reference_current,
        "restart_drop_total_current_A_per_um": restart_current,
        "abs_log10_error_vs_reference": abs_log10_error,
        "abs_log10_error_restart_drop": restart_log_error,
        "abs_delta_vs_restart_current_A_per_um": abs_delta_restart,
        "contact_edge_count": len(contact_edges),
        "override_edge_count": len(override_edges),
        "other_contact_override_edge_count": len(other_contact_override_edges),
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
        "qf-floor restart point verified: "
        f"bias={args.bias:g} V contact={args.contact} "
        f"log10_error={abs_log10_error:.6g} "
        f"override_edges={len(override_edges)}"
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iv-csv", type=Path, required=True)
    parser.add_argument("--contact-edge-csv", type=Path, required=True)
    parser.add_argument("--restart-summary", type=Path)
    parser.add_argument("--bias", type=float, default=-12.7)
    parser.add_argument("--contact", default="Anode")
    parser.add_argument("--reference-current", type=float, default=math.nan)
    parser.add_argument("--restart-current", type=float, default=math.nan)
    parser.add_argument("--max-log10-error", type=float, default=0.07)
    parser.add_argument("--current-abs-tol", type=float, default=1.0e-28)
    parser.add_argument("--expected-override-edges", type=int)
    parser.add_argument("--forbid-other-contact-overrides", action="store_true")
    parser.add_argument("--out-json", type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def find_bias_row(rows: list[dict[str, str]], bias: float) -> dict[str, str] | None:
    for row in rows:
        if bias_matches(row.get("bias_V", ""), bias):
            return row
    return None


def bias_matches(value: str, bias: float) -> bool:
    try:
        return abs(float(value) - bias) <= 1.0e-9
    except ValueError:
        return False


def row_is_converged(row: dict[str, str]) -> bool:
    value = row.get("converged")
    if value is None or value == "":
        return True
    return parse_boolish(value)


def parse_boolish(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def read_current_a_per_um(row: dict[str, str], errors: list[str]) -> float:
    if "current_total_A_per_um" in row and row["current_total_A_per_um"] != "":
        return parse_float(row["current_total_A_per_um"], "current_total_A_per_um", errors)
    if "current_total" in row and row["current_total"] != "":
        return parse_float(row["current_total"], "current_total", errors) / 1.0e6
    errors.append("iv row has neither current_total_A_per_um nor current_total")
    return math.nan


def parse_float(value: str, name: str, errors: list[str]) -> float:
    try:
        parsed = float(value)
    except ValueError:
        errors.append(f"cannot parse {name} value {value!r}")
        return math.nan
    if not math.isfinite(parsed):
        errors.append(f"{name} is not finite")
    return parsed


def read_optional_float(
    data: dict[str, Any],
    key: str,
    fallback: float,
) -> float:
    if key not in data:
        return fallback
    try:
        return float(data[key])
    except (TypeError, ValueError):
        return math.nan


def log10_abs_error(candidate: float, reference: float) -> float:
    if candidate == 0.0 or reference == 0.0:
        return math.inf
    return abs(math.log10(abs(candidate)) - math.log10(abs(reference)))


if __name__ == "__main__":
    raise SystemExit(main())
