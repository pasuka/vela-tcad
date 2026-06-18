#!/usr/bin/env python3
"""Run PN2D release BV visual comparison report."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def relative_to(path: Path, base: Path) -> str:
    return str(path.resolve().relative_to(base.resolve()))


def load_curve_points(path: Path, current_col_preference: list[str]) -> list[tuple[float, float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = reader.fieldnames or []
    if not rows:
        return []
    bias_col = next((name for name in fields if name in {"bias_V", "voltage", "V"}), fields[0])
    current_col = next((name for name in current_col_preference if name in fields), fields[-1])
    points: list[tuple[float, float]] = []
    for row in rows:
        if "converged" in row and str(row.get("converged", "")).strip() not in {"1", "true", "True"}:
            continue
        try:
            points.append((float(row[bias_col]), float(row[current_col])))
        except (TypeError, ValueError):
            continue
    return sorted(points, key=lambda item: item[0])


def interpolate(points: list[tuple[float, float]], bias: float) -> float | None:
    if not points:
        return None
    for point_bias, value in points:
        if abs(point_bias - bias) <= 1.0e-10:
            return value
    for (b0, y0), (b1, y1) in zip(points, points[1:]):
        if b0 <= bias <= b1:
            weight = (bias - b0) / (b1 - b0)
            return y0 + weight * (y1 - y0)
    return None


def finite_abs(values: list[float]) -> np.ndarray:
    arr = np.asarray([abs(value) for value in values], dtype=float)
    return arr[np.isfinite(arr) & (arr > 0.0)]


def write_curve_plot(
    reference_csv: Path,
    candidate_csv: Path,
    biases: list[float],
    out_dir: Path,
) -> dict[str, Any]:
    ref = load_curve_points(reference_csv, ["current_total", "current"])
    cand = load_curve_points(candidate_csv, ["current_total_A_per_um", "current_total", "current"])
    min_bias = min(biases)
    ref_window = [(bias, cur) for bias, cur in ref if min_bias <= bias <= 0.0]
    cand_window = [(bias, cur) for bias, cur in cand if min_bias <= bias <= 0.0]
    if not ref_window:
        raise RuntimeError(f"no Sentaurus BV reference points in [{min_bias:g}, 0] V")
    if not cand_window:
        raise RuntimeError(f"no Vela BV candidate points in [{min_bias:g}, 0] V")

    fig, ax = plt.subplots(figsize=(8.4, 5.0), constrained_layout=True)
    ax.semilogy(
        [-bias for bias, _ in ref_window],
        finite_abs([cur for _, cur in ref_window]),
        color="#2A6FBB",
        linewidth=1.9,
        label="Sentaurus reference |I| (A)",
    )
    ax.semilogy(
        [-bias for bias, _ in cand_window],
        finite_abs([cur for _, cur in cand_window]),
        color="#D1495B",
        marker="o",
        linewidth=1.8,
        markersize=4.5,
        label="Vela release |I| (A/um)",
    )
    ax.set_xlabel("Reverse bias -V_A (V)")
    ax.set_ylabel("Absolute current (log scale)")
    ax.set_title("PN2D BV curve comparison, logarithmic current axis")
    ax.grid(True, which="both", linewidth=0.45, alpha=0.34)
    ax.legend(loc="best")
    png_path = out_dir / "pn2d_bv_curve_comparison.png"
    fig.savefig(png_path, dpi=220)
    plt.close(fig)

    point_csv = out_dir / "pn2d_bv_curve_points.csv"
    rows: list[dict[str, Any]] = []
    with point_csv.open("w", newline="") as handle:
        fieldnames = [
            "bias_V",
            "reverse_bias_V",
            "sentaurus_current_A",
            "vela_current_A_per_um",
            "log10_abs_ratio",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for bias in sorted(set([0.0] + biases)):
            s_val = interpolate(ref, bias)
            v_val = interpolate(cand, bias)
            log_ratio = ""
            if (
                s_val is not None
                and v_val is not None
                and abs(s_val) > 0.0
                and abs(v_val) > 0.0
            ):
                log_ratio = math.log10(abs(v_val) / abs(s_val))
            row = {
                "bias_V": bias,
                "reverse_bias_V": -bias,
                "sentaurus_current_A": s_val if s_val is not None else "",
                "vela_current_A_per_um": v_val if v_val is not None else "",
                "log10_abs_ratio": log_ratio,
            }
            rows.append(row)
            writer.writerow(row)

    with Image.open(png_path) as image:
        width, height = image.size
    return {
        "path": str(png_path),
        "points_csv": str(point_csv),
        "width_px": width,
        "height_px": height,
        "rows": rows,
    }


def run_command(cmd: list[str], cwd: Path) -> None:
    print(":: running", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def build_run_config(base_config: Path, run_dir: Path, biases: list[float]) -> Path:
    with base_config.open(encoding="utf-8") as handle:
        config = json.load(handle)
    base_dir = base_config.parent
    for key in ("mesh_file", "node_doping_file"):
        config[key] = str((base_dir / config[key]).resolve())
    config["output_csv"] = "pn2d_bv_release_visual.csv"
    sweep = config.setdefault("sweep", {})
    sweep["bias_points"] = [0.0] + biases
    sweep["start"] = 0.0
    sweep["stop"] = min(biases)
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = "vtk/dc_sweep"
    sweep["write_state_file"] = "pn2d_bv_release_visual_last_state.csv"
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "vtk").mkdir(parents=True, exist_ok=True)
    config_path = run_dir / "simulation_bv_release_visual.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path


def parse_args() -> argparse.Namespace:
    root = repo_root()
    default_base = root / "build-release/reference_tcad/pn2d_sentaurus2018/vela/simulation_bv_minus20_sg_edge_current_vtk_probe.json"
    default_runner = root / "build-release/vela_example_runner.exe"
    default_sentaurus = root / "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias"
    default_curve = root / "build-release/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv"
    default_out = root / "build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_release_curve_heatmaps"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, default=default_base)
    parser.add_argument("--runner", type=Path, default=default_runner)
    parser.add_argument("--sentaurus-root", type=Path, default=default_sentaurus)
    parser.add_argument("--curve-reference", type=Path, default=default_curve)
    parser.add_argument("--out-dir", type=Path, default=default_out)
    parser.add_argument("--biases", default="-1,-2,-4,-8")
    parser.add_argument("--quantities", default="potential,electron_density,hole_density")
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = parse_biases(args.biases)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    run_dir = args.out_dir / "vela_run"
    config_path = build_run_config(args.base_config, run_dir, biases)
    csv_path = run_dir / "pn2d_bv_release_visual.csv"
    vtk_dir = run_dir / "vtk"

    if not args.skip_run:
        run_command([str(args.runner.resolve()), "--config", str(config_path.resolve())], run_dir)

    compare_dir = args.out_dir / "field_compare"
    run_command([
        sys.executable,
        str((repo_root() / "scripts/compare_pn2d_bv_multibias_fields.py").resolve()),
        "--sentaurus-root",
        str(args.sentaurus_root.resolve()),
        "--vela-vtk-root",
        str(vtk_dir.resolve()),
        "--curve-reference",
        str(args.curve_reference.resolve()),
        "--curve-candidate",
        str(csv_path.resolve()),
        "--out-dir",
        str(compare_dir.resolve()),
        "--biases",
        ",".join(f"{bias:g}" for bias in biases),
        "--quantities",
        args.quantities,
    ], repo_root())

    heatmap_dir = args.out_dir / "heatmaps"
    run_command([
        sys.executable,
        str((repo_root() / "scripts/plot_pn2d_multibias_heatmaps.py").resolve()),
        "--sentaurus-root",
        str(args.sentaurus_root.resolve()),
        "--vela-vtk-root",
        str(vtk_dir.resolve()),
        "--out-dir",
        str(heatmap_dir.resolve()),
        "--biases",
        ",".join(f"{bias:g}" for bias in biases),
        "--quantities",
        args.quantities,
    ], repo_root())

    curve = write_curve_plot(args.curve_reference, csv_path, biases, args.out_dir)
    manifest = {
        "schema": "pn2d.release_bv_visual_report.v1",
        "config": str(config_path),
        "vela_csv": str(csv_path),
        "vela_vtk_dir": str(vtk_dir),
        "curve": curve,
        "field_compare_dir": str(compare_dir),
        "heatmap_manifest": str(heatmap_dir / "pn2d_heatmap_manifest.json"),
        "biases_V": biases,
        "quantities": args.quantities.split(","),
    }
    manifest_path = args.out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
