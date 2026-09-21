#!/usr/bin/env python3
"""Score the exact-point Templates/LDMOS G3 Id-Vg qualification curve."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


L2_LIMITS = {
    "median_abs_log_error_dex": 0.10,
    "p95_abs_log_error_dex": 0.20,
    "strong_inversion_endpoint_relative_error": 0.20,
    "vth_absolute_error_V": 0.10,
    "max_gm_relative_error": 0.20,
    "kcl_relative_error": 0.01,
}
FIXED_CURRENT_A_PER_UM = 1.0e-8


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def log_intersection(points: list[tuple[float, float]], target: float) -> float:
    log_target = math.log10(target)
    for (v0, i0), (v1, i1) in zip(points, points[1:]):
        if i0 <= 0.0 or i1 <= 0.0:
            continue
        y0 = math.log10(i0)
        y1 = math.log10(i1)
        if min(y0, y1) <= log_target <= max(y0, y1) and y1 != y0:
            return v0 + (v1 - v0) * (log_target - y0) / (y1 - y0)
    return math.nan


def max_gm(points: list[tuple[float, float]]) -> dict[str, float]:
    candidates = [
        ((v0 + v1) / 2.0, (i1 - i0) / (v1 - v0))
        for (v0, i0), (v1, i1) in zip(points, points[1:])
        if v1 != v0
    ]
    voltage, gm = max(candidates, key=lambda item: item[1])
    return {"value_A_per_um_V": gm, "midpoint_V": voltage}


def markdown(report: dict[str, Any]) -> str:
    metrics = report["metrics"]
    gates = report["gates"]
    lines = [
        "# Templates/LDMOS G3 Id-Vg exact-point qualification",
        "",
        f"Status: {report['status']}",
        "",
        f"Exact shared points: {metrics['exact_shared_points']}; resolved: "
        f"{metrics['resolved_points']}; unresolved biases: "
        f"{metrics['unresolved_biases_V']}.",
        "",
        "| Metric | Value | L2 limit | Pass |",
        "| --- | ---: | ---: | --- |",
    ]
    rows = (
        ("Median abs log error (dex)", "median_abs_log_error_dex"),
        ("P95 abs log error (dex)", "p95_abs_log_error_dex"),
        ("Strong-inversion endpoint relative error", "strong_inversion_endpoint_relative_error"),
        ("Fixed-current Vth absolute error (V)", "vth_absolute_error_V"),
        ("Maximum gm relative error", "max_gm_relative_error"),
        ("Worst resolved KCL relative error", "kcl_relative_error"),
    )
    for label, key in rows:
        lines.append(
            f"| {label} | {metrics[key]:.8g} | {L2_LIMITS[key]:.8g} | "
            f"{gates[key]} |"
        )
    lines.extend([
        "",
        "Curve errors are evaluated only at exact shared CurrentPlot points; no curve interpolation is used.",
        f"Vth uses local log-current interpolation at {FIXED_CURRENT_A_PER_UM:.1e} A/um and remains diagnostic until that current level is formally frozen.",
    ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--balance", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()

    reference_rows = read_csv(args.reference)
    candidate_rows = read_csv(args.candidate)
    balance_rows = read_csv(args.balance)
    if len(reference_rows) != len(candidate_rows) or len(reference_rows) != len(balance_rows):
        raise ValueError("Id-Vg scoring requires equal-length exact bias grids")
    shared: list[float] = []
    reference: dict[float, float] = {}
    candidate: dict[float, float] = {}
    balance: dict[float, dict[str, str]] = {}
    for ref_row, cand_row, balance_row in zip(
        reference_rows, candidate_rows, balance_rows
    ):
        ref_bias = float(ref_row["bias_V"])
        cand_bias = float(cand_row["bias_V"])
        balance_bias = float(balance_row["bias_V"])
        if abs(cand_bias - ref_bias) > 1.0e-12 or abs(balance_bias - ref_bias) > 1.0e-12:
            raise ValueError("Id-Vg scoring requires identical ordered bias grids")
        shared.append(ref_bias)
        reference[ref_bias] = float(ref_row["current_total_A_per_um"])
        candidate[ref_bias] = float(cand_row["current_total_A_per_um"])
        balance[ref_bias] = balance_row

    resolved = [
        bias for bias in shared
        if balance[bias]["numerical_status"] == "resolved"
    ]
    errors = [
        abs(math.log10(abs(candidate[bias]) / abs(reference[bias])))
        for bias in resolved
        if reference[bias] != 0.0 and candidate[bias] != 0.0
    ]
    point_metrics = [
        {
            "bias_V": bias,
            "reference_A_per_um": reference[bias],
            "candidate_A_per_um": candidate[bias],
            "abs_log_error_dex": abs(
                math.log10(abs(candidate[bias]) / abs(reference[bias]))
            ),
            "kcl_relative_error": (
                1.0 / float(balance[bias]["id_to_kcl_residual_ratio"])
            ),
        }
        for bias in resolved
        if reference[bias] != 0.0
        and candidate[bias] != 0.0
        and float(balance[bias]["id_to_kcl_residual_ratio"]) > 0.0
    ]
    reference_points = [(bias, abs(reference[bias])) for bias in shared]
    candidate_points = [(bias, abs(candidate[bias])) for bias in shared]
    reference_vth = log_intersection(reference_points, FIXED_CURRENT_A_PER_UM)
    candidate_vth = log_intersection(candidate_points, FIXED_CURRENT_A_PER_UM)
    reference_gm = max_gm(reference_points)
    candidate_gm = max_gm(candidate_points)
    endpoint = shared[-1]
    kcl_errors = [
        1.0 / float(balance[bias]["id_to_kcl_residual_ratio"])
        for bias in resolved
        if float(balance[bias]["id_to_kcl_residual_ratio"]) > 0.0
    ]

    metrics = {
        "exact_shared_points": len(shared),
        "resolved_points": len(resolved),
        "unresolved_biases_V": [bias for bias in shared if bias not in resolved],
        "median_abs_log_error_dex": statistics.median(errors),
        "p95_abs_log_error_dex": percentile(errors, 0.95),
        "max_abs_log_error_dex": max(errors),
        "strong_inversion_endpoint_relative_error": abs(
            candidate[endpoint] - reference[endpoint]) / abs(reference[endpoint]),
        "reference_endpoint_A_per_um": reference[endpoint],
        "candidate_endpoint_A_per_um": candidate[endpoint],
        "fixed_current_A_per_um": FIXED_CURRENT_A_PER_UM,
        "reference_vth_V": reference_vth,
        "candidate_vth_V": candidate_vth,
        "vth_absolute_error_V": abs(candidate_vth - reference_vth),
        "reference_max_gm": reference_gm,
        "candidate_max_gm": candidate_gm,
        "max_gm_relative_error": abs(
            candidate_gm["value_A_per_um_V"] - reference_gm["value_A_per_um_V"]
        ) / abs(reference_gm["value_A_per_um_V"]),
        "kcl_relative_error": max(kcl_errors),
        "worst_log_error_points": sorted(
            point_metrics, key=lambda row: row["abs_log_error_dex"], reverse=True
        )[:5],
        "worst_kcl_point": max(
            point_metrics, key=lambda row: row["kcl_relative_error"]
        ),
    }
    gates = {
        key: metrics[key] <= limit for key, limit in L2_LIMITS.items()
    }
    report = {
        "schema": "vela.templates_ldmos.g3_idvg_qualification.v1",
        "alignment": "exact shared points; no curve interpolation",
        "reference": str(args.reference.resolve()),
        "candidate": str(args.candidate.resolve()),
        "balance": str(args.balance.resolve()),
        "limits": L2_LIMITS,
        "metrics": metrics,
        "gates": gates,
        "status": "pass" if all(gates.values()) else "fail",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(markdown(report), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
