#!/usr/bin/env python3
"""Run a dense PN2D release BV current-curve comparison."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from PIL import Image, ImageStat


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def load_reference_points(path: Path) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                rows.append((float(row["bias_V"]), float(row["current_total"])))
            except (KeyError, TypeError, ValueError):
                continue
    return rows


def unique_biases(points: list[tuple[float, float]], min_bias: float, max_bias: float) -> list[float]:
    result: list[float] = []
    seen: set[float] = set()
    for bias, _ in points:
        if bias < min_bias - 1.0e-12 or bias > max_bias + 1.0e-12:
            continue
        key = round(bias, 12)
        if key in seen:
            continue
        seen.add(key)
        result.append(bias)
    result.sort(reverse=True)
    return result


def load_vela_points(path: Path) -> list[tuple[float, float]]:
    if not path.exists():
        return []
    points: list[tuple[float, float]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        current_col = "current_total_A_per_um" if "current_total_A_per_um" in fields else "current_total"
        for row in reader:
            if str(row.get("converged", "1")).strip() not in {"1", "true", "True"}:
                continue
            try:
                points.append((float(row["bias_V"]), float(row[current_col])))
            except (KeyError, TypeError, ValueError):
                continue
    return sorted(points, reverse=True)


def interpolate(points: list[tuple[float, float]], bias: float) -> float | None:
    if not points:
        return None
    ordered = sorted(points)
    for point_bias, value in ordered:
        if abs(point_bias - bias) <= 1.0e-10:
            return value
    for (b0, y0), (b1, y1) in zip(ordered, ordered[1:]):
        if b0 <= bias <= b1:
            weight = (bias - b0) / (b1 - b0)
            return y0 + weight * (y1 - y0)
    return None


def build_config(base_config: Path, out_dir: Path, biases: list[float]) -> Path:
    with base_config.open(encoding="utf-8") as handle:
        config = json.load(handle)
    base_dir = base_config.parent
    for key in ("mesh_file", "node_doping_file"):
        config[key] = str((base_dir / config[key]).resolve())
    run_dir = out_dir / "vela_run"
    run_dir.mkdir(parents=True, exist_ok=True)
    config["output_csv"] = "pn2d_bv_dense_curve.csv"
    sweep = config.setdefault("sweep", {})
    sweep["start"] = 0.0
    sweep["stop"] = min(biases)
    sweep["bias_points"] = biases
    sweep["write_vtk"] = False
    sweep["write_state_file"] = "pn2d_bv_dense_curve_last_state.csv"
    config_path = run_dir / "simulation_bv_dense_curve.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path


def write_curve_outputs(
    reference: list[tuple[float, float]],
    vela: list[tuple[float, float]],
    out_dir: Path,
) -> dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_window = [(b, i) for b, i in reference if -20.0000001 <= b <= 0.0000001 and i != 0.0]
    vela_window = [(b, i) for b, i in vela if -20.0000001 <= b <= 0.0000001 and i != 0.0]
    ratio_rows: list[dict[str, Any]] = []
    for bias, vela_current in vela_window:
        sentaurus_current = interpolate(reference_window, bias)
        if sentaurus_current is None or sentaurus_current == 0.0:
            continue
        ratio = abs(vela_current) / abs(sentaurus_current)
        ratio_rows.append({
            "bias_V": bias,
            "reverse_bias_V": -bias,
            "sentaurus_current_A": sentaurus_current,
            "vela_current_A_per_um": vela_current,
            "vela_over_sentaurus_percent": 100.0 * ratio,
            "log10_abs_ratio": math.log10(ratio) if ratio > 0.0 else "",
        })

    points_csv = out_dir / "pn2d_bv_dense_curve_points.csv"
    with points_csv.open("w", newline="") as handle:
        fieldnames = [
            "bias_V",
            "reverse_bias_V",
            "sentaurus_current_A",
            "vela_current_A_per_um",
            "vela_over_sentaurus_percent",
            "log10_abs_ratio",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(ratio_rows)

    fig, (ax_curve, ax_ratio) = plt.subplots(
        2,
        1,
        figsize=(9.8, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )
    ax_curve.semilogy(
        [-b for b, _ in reference_window],
        [abs(i) for _, i in reference_window],
        color="#2A6FBB",
        linewidth=2.0,
        label=f"Sentaurus reference, {len(reference_window)} rows |I| (A)",
    )
    ax_curve.semilogy(
        [-b for b, _ in vela_window],
        [abs(i) for _, i in vela_window],
        color="#D1495B",
        marker=".",
        markersize=4.0,
        linewidth=1.6,
        label=f"Vela release, {len(vela_window)} converged points |I| (A/um)",
    )
    ax_curve.set_ylabel("Absolute current (log scale)")
    ax_curve.set_title("PN2D BV current comparison, Sentaurus-matched bias points")
    ax_curve.grid(True, which="both", linewidth=0.45, alpha=0.34)
    ax_curve.legend(loc="best")

    ax_ratio.plot(
        [row["reverse_bias_V"] for row in ratio_rows],
        [row["vela_over_sentaurus_percent"] for row in ratio_rows],
        color="#4B7F52",
        marker=".",
        markersize=4.0,
        linewidth=1.2,
    )
    ax_ratio.axhline(100.0, color="#444444", linewidth=0.8, linestyle="--", alpha=0.65)
    ax_ratio.set_xlim(0.0, 20.0)
    ax_ratio.set_xlabel("Reverse bias -V_A (V)")
    ax_ratio.set_ylabel("Vela/Sentaurus (%)")
    ax_ratio.grid(True, which="major", linewidth=0.45, alpha=0.34)

    plot_path = out_dir / "pn2d_bv_dense_curve_comparison.png"
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)
    with Image.open(plot_path) as image:
        stat = ImageStat.Stat(image.convert("RGB"))
        image_check = {
            "width_px": image.size[0],
            "height_px": image.size[1],
            "stddev_rgb": stat.stddev,
        }
    return {
        "plot": str(plot_path),
        "points_csv": str(points_csv),
        "reference_rows_in_plot": len(reference_window),
        "vela_points_in_plot": len(vela_window),
        "matched_ratio_rows": len(ratio_rows),
        "image_check": image_check,
    }


def parse_args() -> argparse.Namespace:
    root = repo_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=root / "build-release/reference_tcad/pn2d_sentaurus2018/vela/simulation_bv_minus20_sg_edge_current_vtk_probe.json",
    )
    parser.add_argument("--runner", type=Path, default=root / "build-release/vela_example_runner.exe")
    parser.add_argument(
        "--reference-csv",
        type=Path,
        default=root / "build-release/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_release_dense_curve",
    )
    parser.add_argument("--min-bias", type=float, default=-20.0)
    parser.add_argument("--max-bias", type=float, default=0.0)
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    reference = load_reference_points(args.reference_csv)
    biases = unique_biases(reference, args.min_bias, args.max_bias)
    if not biases or biases[0] != 0.0:
        biases.insert(0, 0.0)
    config_path = build_config(args.base_config, args.out_dir, biases)
    run_dir = args.out_dir / "vela_run"
    vela_csv = run_dir / "pn2d_bv_dense_curve.csv"

    run_return_code = 0
    if not args.skip_run:
        completed = subprocess.run(
            [str(args.runner.resolve()), "--config", str(config_path.resolve())],
            cwd=run_dir,
            check=False,
        )
        run_return_code = completed.returncode

    vela = load_vela_points(vela_csv)
    curve_outputs = write_curve_outputs(reference, vela, args.out_dir)
    manifest = {
        "schema": "pn2d.bv_dense_curve_report.v1",
        "config": str(config_path),
        "reference_csv": str(args.reference_csv),
        "requested_bias_points": len(biases),
        "requested_min_bias_V": min(biases),
        "requested_max_bias_V": max(biases),
        "vela_csv": str(vela_csv),
        "vela_converged_points": len(vela),
        "run_return_code": run_return_code,
        **curve_outputs,
    }
    manifest_path = args.out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return run_return_code


if __name__ == "__main__":
    raise SystemExit(main())
