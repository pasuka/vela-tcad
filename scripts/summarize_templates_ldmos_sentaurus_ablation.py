#!/usr/bin/env python3
"""Normalize and compare exact-point Templates/LDMOS Sentaurus ablations."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any, Sequence

from run_templates_ldmos_sentaurus_vm import normalize_plt_files


CHAIN = (
    ("G0-original", None, "official full-physics baseline"),
    ("G1-no-hQP", "G0-original", "hQP increment"),
    ("G2-no-eQP", "G1-no-hQP", "eQP increment"),
    ("G3-no-IALMob", "G2-no-eQP", "IALMob increment"),
    ("G4-no-highfield", "G3-no-IALMob", "high-field increment"),
    ("G-contact-poly-barrier-control", "G0-original", "PolySi/barrier control"),
)


def read_curve(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [
            (float(row["bias_V"]), float(row["current_total_A_per_um"]))
            for row in csv.DictReader(stream)
        ]


def distribution(values: list[float]) -> dict[str, float | int]:
    if not values:
        return {"count": 0, "median": math.nan, "max": math.nan}
    return {"count": len(values), "median": statistics.median(values), "max": max(values)}


def crossing(curve: list[tuple[float, float]], target: float) -> float | None:
    """Interpolate Vg locally in log-current; never interpolate curve scores."""
    log_target = math.log10(target)
    for (v0, i0), (v1, i1) in zip(curve, curve[1:]):
        a, b = abs(i0), abs(i1)
        if a <= 0.0 or b <= 0.0 or a == b:
            continue
        lo, hi = sorted((a, b))
        if lo <= target <= hi:
            fraction = (log_target - math.log10(a)) / (math.log10(b) - math.log10(a))
            return v0 + fraction * (v1 - v0)
    return None


def compare(parent: list[tuple[float, float]], child: list[tuple[float, float]]) -> dict[str, Any]:
    if len(parent) != len(child):
        raise ValueError("ablation curves have different exact-point counts")
    for (vp, _), (vc, _) in zip(parent, child):
        if abs(vp - vc) > 1.0e-12:
            raise ValueError(f"ablation bias grids differ: {vp} versus {vc}")
    resolved = [
        (vp, abs(ip), abs(ic))
        for (vp, ip), (_, ic) in zip(parent, child)
        if abs(ip) >= 1.0e-30 and abs(ic) >= 1.0e-30
    ]
    strong = [(ip, ic) for bias, ip, ic in resolved if bias >= 4.0]
    vth: dict[str, Any] = {}
    for target in (1.0e-12, 1.0e-10, 1.0e-8):
        left, right = crossing(parent, target), crossing(child, target)
        vth[f"{target:.0e}_A_per_um"] = {
            "parent_V": left,
            "child_V": right,
            "absolute_shift_mV": (
                abs(left - right) * 1.0e3 if left is not None and right is not None else None),
            "classification": "diagnostic_until_fixed_current_level_is_frozen",
        }
    return {
        "exact_bias_point_count": len(parent),
        "resolved_log_error_dex": distribution([
            abs(math.log10(ip) - math.log10(ic)) for _, ip, ic in resolved]),
        "strong_inversion_relative_change_percent": distribution([
            abs(ic / ip - 1.0) * 100.0 for ip, ic in strong if ip > 0.0]),
        "fixed_current_vth_diagnostics": vth,
    }


def wall_seconds(path: Path) -> float | None:
    text = path.read_text(encoding="utf-8", errors="ignore")
    match = re.search(
        r"Elapsed \(wall clock\) time \([^\n]+\):\s*([0-9:.]+)\s*$",
        text,
        re.MULTILINE,
    )
    if not match:
        return None
    parts = [float(item) for item in match.group(1).split(":")]
    if len(parts) == 2:
        return parts[0] * 60.0 + parts[1]
    if len(parts) == 3:
        return parts[0] * 3600.0 + parts[1] * 60.0 + parts[2]
    return None


def summarize(run_root: Path, output: Path) -> dict[str, Any]:
    curves: dict[str, list[tuple[float, float]]] = {}
    runs: list[dict[str, Any]] = []
    normalized_root = output / "normalized"
    for index, (name, _, _) in enumerate(CHAIN):
        source = run_root / name
        normalized = normalized_root / f"v{index}"
        files = normalize_plt_files(source, normalized)
        curve_file = normalized / "IdVg_n2_des_drain_curve.csv"
        curves[name] = read_curve(curve_file)
        exit_code = int((source / "exitcode").read_text(encoding="utf-8").strip())
        runs.append({
            "id": name,
            "exit_code": exit_code,
            "status": "pass" if exit_code == 0 else "fail",
            "wall_clock_seconds": wall_seconds(source / "timing.txt"),
            "normalization": files,
        })
    comparisons = []
    for name, parent, factor in CHAIN:
        if parent is None:
            continue
        comparisons.append({
            "id": f"{parent}_minus_{name}",
            "parent": parent,
            "child": name,
            "factor": factor,
            "metrics": compare(curves[parent], curves[name]),
        })
    summary = {
        "schema": "vela.templates_ldmos.sentaurus_ablation_summary.v1",
        "benchmark": "sentaurus_t2022_03_sp2_templates_ldmos",
        "status": "pass" if all(item["status"] == "pass" for item in runs) else "fail",
        "alignment_policy": "exact shared CurrentPlot points; no interpolation in curve scores",
        "runs": runs,
        "comparisons": comparisons,
        "decision_limits": {
            "hqp_strong_inversion_current_delta_percent": 2.0,
            "hqp_fixed_current_vth_shift_mV": 10.0,
            "hqp_vth_status": "not_decidable_until_fixed_current_level_is_frozen",
        },
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "sentaurus_ablation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Templates/LDMOS Sentaurus Id-Vg ablation summary", "",
        f"Status: **{summary['status']}**.", "",
        "| Factor | Pair | Median log delta (dex) | Strong-inversion median/max |",
        "| --- | --- | ---: | ---: |",
    ]
    for item in comparisons:
        log = item["metrics"]["resolved_log_error_dex"]
        strong = item["metrics"]["strong_inversion_relative_change_percent"]
        lines.append(
            f"| {item['factor']} | `{item['parent']}` / `{item['child']}` | "
            f"{log['median']:.6g} | {strong['median']:.6g}% / {strong['max']:.6g}% |")
    lines.extend(["", "Fixed-current Vth shifts are diagnostic only because the current level is not yet frozen."])
    (output / "sentaurus_ablation_summary.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = summarize(args.run_root.resolve(), args.output_dir.resolve())
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
