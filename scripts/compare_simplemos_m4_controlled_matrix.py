#!/usr/bin/env python3
"""Compare SimpleMOS M4 curves on the frozen exact bias lattice."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m4_controlled_mobility_contract_v1.json"
)
DEFAULT_REFERENCES = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "controlled_mobility"
)
DEFAULT_VELA = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m4_controlled_mobility" / "vela"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO).as_posix()
    except ValueError:
        return str(resolved)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def exact_curve(path: Path, biases: list[float], bias_column: str,
                current_column: str, tolerance: float,
                require_converged: bool) -> list[tuple[float, float]]:
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"{path}: empty curve")
    required = {bias_column, current_column}
    if require_converged:
        required.add("converged")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    if require_converged and any(row["converged"] != "1" for row in rows):
        raise ValueError(f"{path}: candidate contains a non-converged point")
    observed = [float(row[bias_column]) for row in rows]
    if len(observed) != len(biases):
        raise ValueError(
            f"{path}: exact lattice requires {len(biases)} rows, got {len(observed)}")
    pairs: list[tuple[float, float]] = []
    for target in biases:
        matches = [row for row in rows
                   if abs(float(row[bias_column]) - target) <= tolerance]
        if len(matches) != 1:
            raise ValueError(
                f"{path}: expected exactly one direct row at {target:g} V, "
                f"got {len(matches)}")
        current = float(matches[0][current_column])
        if not math.isfinite(current):
            raise ValueError(f"{path}: non-finite current at {target:g} V")
        pairs.append((target, current))
    return pairs


def trend(values: list[float]) -> str:
    diffs = [right - left for left, right in zip(values, values[1:])]
    tolerance = max(max((abs(value) for value in values), default=0.0) * 1e-12,
                    1e-300)
    if all(value >= -tolerance for value in diffs):
        return "nondecreasing"
    if all(value <= tolerance for value in diffs):
        return "nonincreasing"
    return "mixed"


def compare_case(reference: Path, candidate: Path,
                 contract: dict[str, Any]) -> dict[str, Any]:
    comparison = contract["comparison"]
    biases = [float(value) for value in
              contract["bias_matrix"]["gate_lattice"]["values_V"]]
    tolerance = float(comparison["exact_bias_tolerance_V"])
    reference_points = exact_curve(
        reference, biases, comparison["reference_bias_column"],
        comparison["reference_current_column"], tolerance, False)
    candidate_points = exact_curve(
        candidate, biases, comparison["candidate_bias_column"],
        comparison["candidate_current_column"], tolerance, True)
    floor = float(comparison["current_floor_A_per_um"])
    log_errors: list[float] = []
    relative_errors: list[float] = []
    rows: list[dict[str, float | str]] = []
    for (bias, reference_current), (_, candidate_current) in zip(
            reference_points, candidate_points):
        ref = abs(reference_current)
        cand = abs(candidate_current)
        above_floor = ref >= floor and cand > 0.0
        log_error = abs(math.log10(cand / ref)) if above_floor and ref > 0.0 else math.nan
        relative = abs(cand - ref) / ref if above_floor and ref > 0.0 else math.nan
        if math.isfinite(log_error):
            log_errors.append(log_error)
            relative_errors.append(relative)
        rows.append({
            "gate_voltage_V": bias,
            "sentaurus_current_A_per_um": reference_current,
            "vela_current_A_per_um": candidate_current,
            "above_current_floor": "1" if above_floor else "0",
            "absolute_log10_ratio": log_error,
            "relative_error": relative,
        })
    reference_trend = trend([abs(value) for _, value in reference_points])
    candidate_trend = trend([abs(value) for _, value in candidate_points])
    maximum_log = max(log_errors) if log_errors else None
    median_log = statistics.median(log_errors) if log_errors else None
    maximum_relative = max(relative_errors) if relative_errors else None
    endpoint_ref = abs(reference_points[-1][1])
    endpoint_cand = abs(candidate_points[-1][1])
    endpoint_log = (
        math.log10(endpoint_cand / endpoint_ref)
        if endpoint_ref > 0.0 and endpoint_cand > 0.0 else None)
    numeric = comparison["numeric_parity"]
    failures: list[str] = []
    if maximum_log is None:
        failures.append("no points are above the predeclared current floor")
    elif maximum_log > float(numeric["maximum_absolute_log10_ratio"]):
        failures.append("maximum absolute log10 ratio exceeds contract")
    if maximum_relative is not None and maximum_relative > float(
            numeric["maximum_relative_error"]):
        failures.append("maximum relative error exceeds contract")
    if bool(numeric["trend_match_required"]) and reference_trend != candidate_trend:
        failures.append("trend mismatch")
    return {
        "status": "pass" if not failures else "fail",
        "reference": portable_path(reference),
        "candidate": portable_path(candidate),
        "point_count": len(rows),
        "points_above_floor": len(log_errors),
        "maximum_absolute_log10_ratio_above_floor": maximum_log,
        "median_absolute_log10_ratio_above_floor": median_log,
        "maximum_relative_error_above_floor": maximum_relative,
        "endpoint_log10_ratio": endpoint_log,
        "reference_trend": reference_trend,
        "candidate_trend": candidate_trend,
        "trend_match": reference_trend == candidate_trend,
        "failures": failures,
        "rows": rows,
    }


def write_case_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCES)
    parser.add_argument("--vela-dir", type=Path, default=DEFAULT_VELA)
    parser.add_argument(
        "--variant-vela-dir", action="append", default=[], metavar="ID=PATH",
        help="override one variant directory, for example A3=build/.../vela/a3")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    variant_dirs: dict[str, Path] = {}
    for value in args.variant_vela_dir:
        if "=" not in value:
            parser.error("--variant-vela-dir requires ID=PATH")
        variant_id, path = value.split("=", 1)
        if variant_id not in {"A0", "A1", "A2", "A3"}:
            parser.error(f"unknown M4 variant {variant_id!r}")
        variant_dirs[variant_id] = Path(path).resolve()
    contract = read_json(args.contract.resolve())
    cases: list[dict[str, Any]] = []
    for variant in contract["variants"]:
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            tag = format(float(vd), ".12g").replace(".", "p")
            case_name = f"{str(variant['id']).lower()}_vd_{tag}"
            reference = args.reference_dir.resolve() / f"{case_name}_reference.csv"
            variant_id = str(variant["id"])
            variant_dir = variant_dirs.get(
                variant_id, args.vela_dir.resolve() / variant_id.lower())
            candidate = variant_dir / f"vd_{tag}" / "20_gate_sweep.csv"
            try:
                result = compare_case(reference, candidate, contract)
                comparison_csv = (
                    args.output_dir.resolve() / f"{case_name}_comparison.csv")
                write_case_csv(comparison_csv, result.pop("rows"))
                result["comparison_csv"] = portable_path(comparison_csv)
            except (FileNotFoundError, ValueError) as error:
                result = {
                    "status": "qualification_failed",
                    "reference": portable_path(reference),
                    "candidate": portable_path(candidate),
                    "qualification_error": str(error),
                }
            result.update(variant=variant["id"], drain_voltage_V=float(vd),
                          case=case_name)
            cases.append(result)
    report = {
        "schema": "vela.simplemos.sdevice.m4_comparison.v1",
        "status": "pass" if all(item["status"] == "pass" for item in cases)
        else "fail",
        "contract": portable_path(args.contract),
        "contract_sha256": sha256(args.contract.resolve()),
        "interpolation": "forbidden",
        "cases": cases,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "comparison_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": report["status"], "cases": [
        {"case": item["case"], "status": item["status"]} for item in cases]}))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
