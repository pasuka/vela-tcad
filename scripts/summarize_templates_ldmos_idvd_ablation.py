#!/usr/bin/env python3
"""Normalize and score the Templates/LDMOS Sentaurus Id-Vd D0--D5 chain."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any, Sequence

from run_templates_ldmos_sentaurus_vm import normalize_plt_files


CHAIN = (
    ("D0-original", None, "self-heating increment"),
    ("D1-isothermal", "D0-original", "self-heating increment"),
    ("D2-no-hRecVelocity", "D1-isothermal", "hRecVelocity increment"),
    ("D3-no-hQP", "D2-no-hRecVelocity", "hQP increment"),
    ("D4-classical", "D3-no-hQP", "eQP increment"),
    ("D5-no-IALMob", "D4-classical", "IALMob increment"),
)
GATES = ("Vg4", "Vg8")


def read_curve(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        raw = [(float(row["bias_V"]), float(row["current_total_A_per_um"]))
               for row in csv.DictReader(stream)]
    curve: list[tuple[float, float]] = []
    for point in raw:
        if curve and point[0] == curve[-1][0]:
            curve[-1] = point
        else:
            curve.append(point)
    if len(raw) - len(curve) not in {0, 1}:
        raise ValueError("unexpected repeated Id-Vd bias points")
    return curve


def percentile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * q
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] + (position - low) * (ordered[high] - ordered[low])


def distribution(values: list[float]) -> dict[str, float | int]:
    return {
        "count": len(values),
        "median": statistics.median(values) if values else math.nan,
        "p95": percentile(values, 0.95),
        "max": max(values) if values else math.nan,
    }


def align(left: list[tuple[float, float]], right: list[tuple[float, float]]) -> list[tuple[float, float, float]]:
    if len(left) != len(right):
        raise ValueError(f"curve point count differs: {len(left)} versus {len(right)}")
    rows = []
    for (vl, il), (vr, ir) in zip(left, right):
        if abs(vl - vr) > 1.0e-12:
            raise ValueError(f"exact Id-Vd grids differ: {vl} versus {vr}")
        rows.append((vl, il, ir))
    return rows


def compare_curve(parent: list[tuple[float, float]], child: list[tuple[float, float]]) -> dict[str, Any]:
    rows = align(parent, child)
    resolved = [(v, abs(ip), abs(ic)) for v, ip, ic in rows
                if v > 0.0 and abs(ip) >= 1.0e-30 and abs(ic) >= 1.0e-30]
    relative = [abs(ic / ip - 1.0) * 100.0 for _, ip, ic in resolved]
    first = resolved[0]
    endpoint = resolved[-1]
    return {
        "exact_point_count": len(rows),
        "nonzero_relative_change_percent": distribution(relative),
        "low_vd_differential_resistance_change_percent":
            abs((first[0] / first[2]) / (first[0] / first[1]) - 1.0) * 100.0,
        "endpoint_current_change_percent": abs(endpoint[2] / endpoint[1] - 1.0) * 100.0,
        "endpoint_bias_V": endpoint[0],
    }


def gate_ratio_change(parent: dict[str, list[tuple[float, float]]],
                      child: dict[str, list[tuple[float, float]]]) -> dict[str, float | int]:
    p4p8 = zip(parent["Vg4"], parent["Vg8"])
    c4c8 = zip(child["Vg4"], child["Vg8"])
    values: list[float] = []
    for ((vp4, ip4), (vp8, ip8)), ((vc4, ic4), (vc8, ic8)) in zip(p4p8, c4c8):
        if max(abs(vp4 - vp8), abs(vp4 - vc4), abs(vp4 - vc8)) > 1.0e-12:
            raise ValueError("gate-ratio curves do not share exact drain points")
        if vp4 > 0.0 and min(abs(ip4), abs(ip8), abs(ic4), abs(ic8)) >= 1.0e-30:
            values.append(abs((ic8 / ic4) / (ip8 / ip4) - 1.0) * 100.0)
    return distribution(values)


def hrec_decision(metrics: dict[str, Any],
                  contact_flux: dict[str, Any] | None = None) -> dict[str, Any]:
    limits = {
        "median_nonzero_current_change_percent": 1.0,
        "p95_nonzero_current_change_percent": 2.4,
        "low_vd_resistance_change_percent": 2.0,
        "endpoint_current_change_percent": 2.0,
        "gate_ratio_p95_change_percent": 1.6,
    }
    observed = {
        "median_nonzero_current_change_percent": max(
            metrics[g]["nonzero_relative_change_percent"]["median"] for g in GATES),
        "p95_nonzero_current_change_percent": max(
            metrics[g]["nonzero_relative_change_percent"]["p95"] for g in GATES),
        "low_vd_resistance_change_percent": max(
            metrics[g]["low_vd_differential_resistance_change_percent"] for g in GATES),
        "endpoint_current_change_percent": max(
            metrics[g]["endpoint_current_change_percent"] for g in GATES),
        "gate_ratio_p95_change_percent": metrics["gate_ratio_change_percent"]["p95"],
    }
    gates = {name: {"observed": observed[name], "limit": limit,
                    "pass": observed[name] <= limit}
             for name, limit in limits.items()}
    flux_pass = contact_flux is not None and (
        contact_flux["resolved_sign_change_count"] == 0
        and contact_flux["max_component_delta_over_terminal_current"] <= 1.0e-6
    )
    curve_pass = all(item["pass"] for item in gates.values())
    return {
        "rule": "implement if any primary metric exceeds 20% of its final Stage-4 tolerance",
        "gates": gates,
        "curve_gate_pass": curve_pass,
        "implementation_required_by_curve": not curve_pass,
        "contact_flux_or_kcl_topology_gate": {
            "status": "pass" if flux_pass else ("pending" if contact_flux is None else "fail"),
            "metrics": contact_flux,
        },
        "implementation_required": not (curve_pass and flux_pass),
    }


def contact_flux_topology(parent: Path, child: Path) -> dict[str, Any]:
    with parent.open(newline="", encoding="utf-8") as stream:
        left = list(csv.DictReader(stream))
    with child.open(newline="", encoding="utf-8") as stream:
        right = list(csv.DictReader(stream))
    if len(left) != len(right):
        raise ValueError("hRec contact-component tables differ in row count")
    columns = [f"{contact} {carrier}Current" for contact in ("source", "drain")
               for carrier in ("e", "h")]
    max_normalized = 0.0
    sign_changes = 0
    max_absolute = 0.0
    for lrow, rrow in zip(left, right):
        bias = float(lrow["drain InnerVoltage"])
        if abs(bias - float(rrow["drain InnerVoltage"])) > 1.0e-12 or bias <= 0.0:
            continue
        terminal_scale = max(abs(float(lrow["drain TotalCurrent"])),
                             abs(float(rrow["drain TotalCurrent"])), 1.0e-30)
        resolved_floor = 1.0e-12 * terminal_scale
        for column in columns:
            lv, rv = float(lrow[column]), float(rrow[column])
            delta = abs(rv - lv)
            max_absolute = max(max_absolute, delta)
            max_normalized = max(max_normalized, delta / terminal_scale)
            if lv * rv < 0.0 and max(abs(lv), abs(rv)) > resolved_floor:
                sign_changes += 1
    return {
        "max_component_abs_delta_A_per_um": max_absolute,
        "max_component_delta_over_terminal_current": max_normalized,
        "resolved_sign_change_count": sign_changes,
        "resolved_component_floor": "1e-12 times same-point drain total current",
    }


def normalize_variant(source: Path, destination: Path) -> dict[str, list[tuple[float, float]]]:
    normalize_plt_files(source, destination)
    return {
        "Vg4": read_curve(destination / "IdVd_Vg1_n4_des_drain_curve.csv"),
        "Vg8": read_curve(destination / "IdVd_Vg2_n4_des_drain_curve.csv"),
    }


def summarize(official_raw: Path, run_root: Path, output: Path) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False)
    normalized = output / "normalized"
    curves: dict[str, dict[str, list[tuple[float, float]]]] = {
        "D0-original": normalize_variant(official_raw, normalized / "D0-original")
    }
    for name, _, _ in CHAIN[1:]:
        curves[name] = normalize_variant(run_root / name, normalized / name)
    comparisons: list[dict[str, Any]] = []
    for name, parent, factor in CHAIN[1:]:
        by_gate = {gate: compare_curve(curves[parent][gate], curves[name][gate])
                   for gate in GATES}
        by_gate["gate_ratio_change_percent"] = gate_ratio_change(curves[parent], curves[name])
        comparisons.append({"id": f"{parent}_minus_{name}", "parent": parent,
                            "child": name, "factor": factor, "metrics": by_gate})
    hrec_metrics = next(item["metrics"] for item in comparisons
                        if item["factor"] == "hRecVelocity increment")
    contact_metrics = {
        gate: contact_flux_topology(
            normalized / "D1-isothermal" / f"IdVd_{'Vg1' if gate == 'Vg4' else 'Vg2'}_n4_des.csv",
            normalized / "D2-no-hRecVelocity" / f"IdVd_{'Vg1' if gate == 'Vg4' else 'Vg2'}_n4_des.csv",
        ) for gate in GATES
    }
    contact_metrics = {
        "max_component_abs_delta_A_per_um": max(
            item["max_component_abs_delta_A_per_um"] for item in contact_metrics.values()),
        "max_component_delta_over_terminal_current": max(
            item["max_component_delta_over_terminal_current"] for item in contact_metrics.values()),
        "resolved_sign_change_count": sum(
            item["resolved_sign_change_count"] for item in contact_metrics.values()),
        "by_gate": contact_metrics,
    }
    result = {
        "schema": "vela.templates_ldmos.sentaurus_idvd_ablation_summary.v1",
        "status": "pass",
        "alignment_policy": "exact shared CurrentPlot points; no curve-score interpolation",
        "comparisons": comparisons,
        "hrec_velocity_decision": hrec_decision(hrec_metrics, contact_metrics),
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    lines = ["# Templates/LDMOS Stage-4 Sentaurus Id-Vd ablation", "",
             "| Factor | Gate | Median / P95 / max current delta | Ron delta | Vd=40 delta |",
             "| --- | --- | ---: | ---: | ---: |"]
    for item in comparisons:
        for gate in GATES:
            metric = item["metrics"][gate]
            dist = metric["nonzero_relative_change_percent"]
            lines.append(f"| {item['factor']} | {gate} | {dist['median']:.5g}% / "
                         f"{dist['p95']:.5g}% / {dist['max']:.5g}% | "
                         f"{metric['low_vd_differential_resistance_change_percent']:.5g}% | "
                         f"{metric['endpoint_current_change_percent']:.5g}% |")
    decision = result["hrec_velocity_decision"]
    lines.extend(["", f"hRecVelocity decision: **{'IMPLEMENT' if decision['implementation_required'] else 'DO NOT IMPLEMENT'}**.",
                  f"Contact-flux topology gate: **{decision['contact_flux_or_kcl_topology_gate']['status']}**."])
    (output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official-raw", type=Path, required=True)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(json.dumps(summarize(args.official_raw.resolve(), args.run_root.resolve(),
                               args.output_dir.resolve()), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
