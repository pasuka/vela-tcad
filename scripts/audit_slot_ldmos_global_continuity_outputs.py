#!/usr/bin/env python3
"""Reclassify existing SLOT-LDMOS IV rows with the two-stage closure gate."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in (None, "") else 0.0


def carrier_gate(contact: float, source: float, source_floor: float,
                 absolute_tolerance: float, relative_tolerance: float) -> dict[str, object]:
    mismatch = contact - source
    physical_scale = max(abs(contact), abs(source))
    qualified = physical_scale >= source_floor
    ratio = abs(mismatch) / max(physical_scale, source_floor, 1.0e-300)
    passed = ratio <= relative_tolerance if qualified else abs(mismatch) <= absolute_tolerance
    return {
        "contact_flux": contact,
        "integrated_source": source,
        "mismatch": mismatch,
        "physical_scale": physical_scale,
        "qualified": qualified,
        "criterion": "relative" if qualified else "absolute",
        "ratio": ratio,
        "passed": passed,
    }


def audit_iv(path: Path, source_floor: float, absolute_tolerance: float,
             relative_tolerance: float) -> dict[str, object] | None:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        return None
    row = rows[-1]
    electron = carrier_gate(
        as_float(row, "global_electron_contact_flux"),
        as_float(row, "global_electron_integrated_source"),
        source_floor, absolute_tolerance, relative_tolerance)
    hole = carrier_gate(
        as_float(row, "global_hole_contact_flux"),
        as_float(row, "global_hole_integrated_source"),
        source_floor, absolute_tolerance, relative_tolerance)
    closure_passed = bool(electron["passed"] and hole["passed"])
    original_converged = row.get("converged", "0") == "1"
    state_exists = (path.parent / "final_state.csv").exists()
    if original_converged and state_exists:
        classification = "existing_checkpoint"
    elif not closure_passed:
        classification = "closure_failure"
    elif state_exists:
        classification = "criterion_candidate_with_state"
    else:
        classification = "criterion_candidate_rerun_required"
    return {
        "case": path.parent.name,
        "iv_csv": str(path.resolve()),
        "bias_V": row.get("bias_V", ""),
        "original_converged": original_converged,
        "state_exists": state_exists,
        "original_failure_reason": row.get("failure_reason", ""),
        "closure_passed": closure_passed,
        "classification": classification,
        "electron": electron,
        "hole": hole,
    }


def write_csv(path: Path, results: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "case", "bias_V", "original_converged", "state_exists",
        "original_failure_reason", "closure_passed", "classification",
        "electron_criterion", "electron_qualified", "electron_mismatch",
        "electron_ratio", "electron_passed", "hole_criterion", "hole_qualified",
        "hole_mismatch", "hole_ratio", "hole_passed", "iv_csv",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for result in results:
            electron = result["electron"]
            hole = result["hole"]
            writer.writerow({
                **{key: result[key] for key in fields if key in result},
                "electron_criterion": electron["criterion"],
                "electron_qualified": electron["qualified"],
                "electron_mismatch": electron["mismatch"],
                "electron_ratio": electron["ratio"],
                "electron_passed": electron["passed"],
                "hole_criterion": hole["criterion"],
                "hole_qualified": hole["qualified"],
                "hole_mismatch": hole["mismatch"],
                "hole_ratio": hole["ratio"],
                "hole_passed": hole["passed"],
            })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--source-floor", type=float, default=1.0e-11)
    parser.add_argument("--absolute-tolerance", type=float, default=1.0e-12)
    parser.add_argument("--relative-tolerance", type=float, default=1.0e-2)
    args = parser.parse_args()
    if args.source_floor < 0 or args.absolute_tolerance < 0 or args.relative_tolerance <= 0:
        parser.error("closure thresholds must be non-negative and relative tolerance positive")

    results = []
    for path in sorted(args.input_root.rglob("iv.csv")):
        result = audit_iv(path, args.source_floor, args.absolute_tolerance,
                          args.relative_tolerance)
        if result is not None:
            results.append(result)
    write_csv(args.output_csv, results)
    summary_path = args.output_csv.with_suffix(".json")
    summary = {
        "schema": "vela.slot_ldmos.global_continuity_audit.v1",
        "source_floor": args.source_floor,
        "absolute_tolerance": args.absolute_tolerance,
        "relative_tolerance": args.relative_tolerance,
        "counts": {
            key: sum(item["classification"] == key for item in results)
            for key in (
                "existing_checkpoint", "criterion_candidate_with_state",
                "criterion_candidate_rerun_required", "closure_failure")
        },
        "results": results,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary["counts"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
