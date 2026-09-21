#!/usr/bin/env python3
"""Score exact-point Vela D5 Id-Vd curves against the Sentaurus D5 reference."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any, Sequence

from summarize_templates_ldmos_idvd_ablation import align, percentile, read_curve


GATES = ("Vg4", "Vg8")
LIMITS = {
    "engineering": {"median": 15.0, "p95": 25.0, "ron": 20.0,
                    "endpoint": 20.0, "gate_ratio": 15.0, "kcl": 1.0},
    "final": {"median": 5.0, "p95": 12.0, "ron": 10.0,
              "endpoint": 10.0, "gate_ratio": 8.0, "kcl": 0.1},
}


def read_vela_curve(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [(float(row["bias_V"]), float(row["current_total_A_per_um"]))
                for row in csv.DictReader(stream)]


def curve_error(reference: list[tuple[float, float]],
                candidate: list[tuple[float, float]]) -> dict[str, Any]:
    candidate_by_bias = {bias: current for bias, current in candidate}
    if len(candidate_by_bias) != len(candidate):
        raise ValueError("candidate curve contains duplicate bias points")
    missing = [bias for bias, _ in reference if bias not in candidate_by_bias]
    if missing:
        raise ValueError(f"candidate curve misses exact reference points: {missing}")
    aligned_candidate = [(bias, candidate_by_bias[bias]) for bias, _ in reference]
    rows = align(reference, aligned_candidate)
    resolved = [(v, abs(ir), abs(ic)) for v, ir, ic in rows
                if v > 0.0 and abs(ir) >= 1.0e-30 and abs(ic) >= 1.0e-30]
    errors = [abs(ic / ir - 1.0) * 100.0 for _, ir, ic in resolved]
    first, endpoint = resolved[0], resolved[-1]
    return {
        "exact_point_count": len(rows),
        "resolved_point_count": len(resolved),
        "relative_error_percent": {
            "median": statistics.median(errors),
            "p95": percentile(errors, 0.95),
            "max": max(errors),
        },
        "low_vd_differential_resistance_error_percent":
            abs((first[0] / first[2]) / (first[0] / first[1]) - 1.0) * 100.0,
        "endpoint_current_error_percent": abs(endpoint[2] / endpoint[1] - 1.0) * 100.0,
        "endpoint_bias_V": endpoint[0],
        "reference_endpoint_A_per_um": endpoint[1],
        "candidate_endpoint_A_per_um": endpoint[2],
    }


def ratio_error(reference: dict[str, list[tuple[float, float]]],
                candidate: dict[str, list[tuple[float, float]]]) -> dict[str, float]:
    reference_rows = align(reference["Vg4"], reference["Vg8"])
    candidate_maps = {
        gate: {bias: current for bias, current in candidate[gate]} for gate in GATES
    }
    errors: list[float] = []
    endpoint = math.nan
    for vr, ir4, ir8 in reference_rows:
        if vr not in candidate_maps["Vg4"] or vr not in candidate_maps["Vg8"]:
            raise ValueError(f"candidate gate-ratio curves miss exact bias {vr}")
        ic4, ic8 = candidate_maps["Vg4"][vr], candidate_maps["Vg8"][vr]
        if vr > 0.0 and min(abs(ir4), abs(ir8), abs(ic4), abs(ic8)) >= 1.0e-30:
            endpoint = abs((ic8 / ic4) / (ir8 / ir4) - 1.0) * 100.0
            errors.append(endpoint)
    return {"median": statistics.median(errors), "p95": percentile(errors, 0.95),
            "max": max(errors), "endpoint": endpoint}


def kcl_audit(path: Path) -> dict[str, Any]:
    """Keep zero-bias evidence explicit; do not silently exempt it from gates."""
    by_point: dict[int, list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            by_point.setdefault(int(row["point_index"]), []).append(row)
    points = []
    for index, rows in by_point.items():
        biases = {float(row["bias_V"]) for row in rows}
        contacts = {row["contact"] for row in rows}
        if len(biases) != 1 or len(rows) != 4 or contacts != {
                "source", "drain", "gate", "substrate"}:
            raise ValueError(f"incomplete or inconsistent terminal balance at point {index}")
        bias = biases.pop()
        currents = [float(row["current_total_A_per_um"]) for row in rows]
        if not math.isfinite(bias) or not all(math.isfinite(value) for value in currents):
            raise ValueError(f"non-finite terminal balance at point {index}")
        imbalance = abs(sum(currents))
        # Frozen validation plan: terminal KCL / maximum terminal current.
        scale = max(max(abs(value) for value in currents), 1.0e-30)
        points.append({
            "point_index": index, "bias_V": bias,
            "absolute_kcl_A_per_um": imbalance,
            "max_absolute_terminal_current_A_per_um": max(map(abs, currents)),
            "normalized_kcl_percent": imbalance / scale * 100.0,
        })
    if not points:
        raise ValueError("terminal balance contains no points")
    nonzero = [row["normalized_kcl_percent"] for row in points if row["bias_V"] != 0.0]
    return {
        "denominator": "maximum_absolute_terminal_current",
        "max_normalized_kcl_percent": max(row["normalized_kcl_percent"] for row in points),
        "max_nonzero_bias_normalized_kcl_percent": max(nonzero) if nonzero else None,
        "zero_bias_points": [row for row in points if row["bias_V"] == 0.0],
        "zero_bias_policy": "raw ratio retained in gate; no equilibrium exemption or resolution floor inferred",
        "point_count": len(points),
    }


def kcl_error(path: Path) -> float:
    return kcl_audit(path)["max_normalized_kcl_percent"]


def verdict(metrics: dict[str, Any], level: str) -> dict[str, Any]:
    limits = LIMITS[level]
    observed = {
        "median": max(metrics[g]["relative_error_percent"]["median"] for g in GATES),
        "p95": max(metrics[g]["relative_error_percent"]["p95"] for g in GATES),
        "ron": max(metrics[g]["low_vd_differential_resistance_error_percent"] for g in GATES),
        "endpoint": max(metrics[g]["endpoint_current_error_percent"] for g in GATES),
        "gate_ratio": metrics["gate_ratio_error_percent"]["endpoint"],
        "kcl": metrics["max_normalized_kcl_percent"],
    }
    gates = {name: {"observed": observed[name], "limit": limit,
                    "pass": observed[name] <= limit}
             for name, limit in limits.items()}
    return {"status": "pass" if all(item["pass"] for item in gates.values()) else "fail",
            "gates": gates}


def analyze(reference_paths: dict[str, Path], candidate_paths: dict[str, Path],
            balance_paths: dict[str, Path], output: Path, *,
            physics_profile: str = "D5") -> dict[str, Any]:
    if physics_profile not in ("D5", "D4"):
        raise ValueError("Unknown physics profile")
    reference = {gate: read_curve(path) for gate, path in reference_paths.items()}
    candidate = {gate: read_vela_curve(path) for gate, path in candidate_paths.items()}
    metrics: dict[str, Any] = {
        gate: curve_error(reference[gate], candidate[gate]) for gate in GATES
    }
    metrics["gate_ratio_error_percent"] = ratio_error(reference, candidate)
    metrics["kcl_audit"] = {gate: kcl_audit(path) for gate, path in balance_paths.items()}
    metrics["max_normalized_kcl_percent"] = max(
        audit["max_normalized_kcl_percent"] for audit in metrics["kcl_audit"].values())
    result = {
        "schema": f"vela.templates_ldmos.stage4_{physics_profile.lower()}_summary.v1",
        "alignment_policy": "31 exact shared CurrentPlot points selected from the solver path; no curve-score interpolation",
        "metrics": metrics,
        "engineering": verdict(metrics, "engineering"),
        "final": verdict(metrics, "final"),
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    lines = [f"# Templates/LDMOS Stage-4 {physics_profile} Id-Vd qualification", "",
             f"Engineering gate: **{result['engineering']['status']}**; final gate: **{result['final']['status']}**.", "",
             "| Gate | median / P95 / max relative error | Ron error | Vd=40 error |",
             "| --- | ---: | ---: | ---: |"]
    for gate in GATES:
        item = metrics[gate]
        rel = item["relative_error_percent"]
        lines.append(f"| {gate} | {rel['median']:.5g}% / {rel['p95']:.5g}% / {rel['max']:.5g}% | "
                     f"{item['low_vd_differential_resistance_error_percent']:.5g}% | "
                     f"{item['endpoint_current_error_percent']:.5g}% |")
    lines.extend(["", f"Endpoint two-gate current-ratio error: {metrics['gate_ratio_error_percent']['endpoint']:.5g}%.",
                  f"Maximum normalized KCL imbalance: {metrics['max_normalized_kcl_percent']:.5g}%."])
    (output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--physics-profile", choices=("D5", "D4"), default="D5")
    for gate in (4, 8):
        parser.add_argument(f"--reference-vg{gate}", type=Path, required=True)
        parser.add_argument(f"--candidate-vg{gate}", type=Path, required=True)
        parser.add_argument(f"--balance-vg{gate}", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = analyze(
        {"Vg4": args.reference_vg4, "Vg8": args.reference_vg8},
        {"Vg4": args.candidate_vg4, "Vg8": args.candidate_vg8},
        {"Vg4": args.balance_vg4, "Vg8": args.balance_vg8},
        args.output_dir.resolve(),
        physics_profile=args.physics_profile,
    )
    print(json.dumps(result, indent=2))
    return 0 if result["engineering"]["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
