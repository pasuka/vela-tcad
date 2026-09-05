#!/usr/bin/env python3
"""Render the three Genius NPN BJT current-semantics figures used by the daily report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import LogNorm

from diagnose_genius_bjt_sdevice_current_semantics import reconstruct, triangle_areas


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
PROBE_ROOT = BUILD_ROOT / "sdevice_current_semantics_probe"

INK = "#27313A"
MUTED = "#66717C"
GRID = "#D9DEE3"
BLUE = "#2D6CDF"
ORANGE = "#D97706"
OLIVE = "#6B7D2A"
PINK = "#C24E70"
OPEN_BLUE = "#DCE7FA"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Microsoft YaHei", "DejaVu Sans"],
            "font.size": 9.5,
            "axes.titlesize": 11,
            "axes.labelsize": 9.5,
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "axes.titlecolor": INK,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.unicode_minus": False,
        }
    )


def load_mesh(export: Path) -> tuple[np.ndarray, np.ndarray]:
    nodes = read_rows(export / "nodes.csv")
    elements = read_rows(export / "elements.csv")
    coordinates = np.asarray(
        [[float(row["x_um"]), float(row["y_um"])] for row in nodes], dtype=float
    )
    triangles = np.asarray(
        [[int(row["node0"]), int(row["node1"]), int(row["node2"])] for row in elements],
        dtype=int,
    )
    return coordinates, triangles


def load_vector(path: Path, key: str) -> np.ndarray:
    rows = read_rows(path)
    ids = [int(row[key]) for row in rows]
    if ids != list(range(len(rows))):
        raise ValueError(f"unordered vector field: {path}")
    values = np.asarray(
        [[float(row["component0"]), float(row["component1"])] for row in rows],
        dtype=float,
    )
    if not np.all(np.isfinite(values)):
        raise ValueError(f"non-finite vector field: {path}")
    return values


def magnitude(values: np.ndarray) -> np.ndarray:
    return np.linalg.norm(values, axis=1)


def abs_log_error(candidate: np.ndarray, reference: np.ndarray) -> np.ndarray:
    floor = np.finfo(float).tiny
    return np.abs(
        np.log10(np.maximum(magnitude(candidate), floor) / np.maximum(magnitude(reference), floor))
    )


def spatial_format(axis: plt.Axes, zoom_top: bool = False) -> None:
    axis.set_aspect("equal")
    axis.set_xlim(-0.04, 6.04)
    axis.set_ylim(0.88 if zoom_top else 2.04, -0.04)
    axis.set_xlabel("横向位置 x（µm）")
    axis.set_ylabel("深度 y（µm）")


def add_panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        0.015,
        0.985,
        label,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=INK,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 2.0},
        zorder=10,
    )


def plot_figure1(
    triangulation: mtri.Triangulation,
    native: np.ndarray,
    reconstructed: np.ndarray,
    tail_nodes: np.ndarray,
    output: Path,
) -> dict[str, float]:
    native_log = np.log10(np.maximum(magnitude(native), 1.0e-30))
    reconstructed_log = np.log10(np.maximum(magnitude(reconstructed), 1.0e-30))
    error = abs_log_error(reconstructed, native)
    color_min, color_max = -5.0, 2.0

    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.55), constrained_layout=True)
    first = axes[0].tripcolor(
        triangulation,
        native_log,
        shading="gouraud",
        cmap="cividis",
        vmin=color_min,
        vmax=color_max,
        rasterized=True,
    )
    axes[1].tripcolor(
        triangulation,
        reconstructed_log,
        shading="gouraud",
        cmap="cividis",
        vmin=color_min,
        vmax=color_max,
        rasterized=True,
    )
    third = axes[2].tripcolor(
        triangulation,
        np.minimum(error, 0.15),
        shading="gouraud",
        cmap="YlOrBr",
        vmin=0.0,
        vmax=0.15,
        rasterized=True,
    )
    axes[2].scatter(
        triangulation.x[tail_nodes],
        triangulation.y[tail_nodes],
        s=9,
        facecolors="none",
        edgecolors=PINK,
        linewidths=0.55,
        label="原180个弱电流尾部节点",
        zorder=6,
    )
    titles = (
        "SDevice 原始节点电流",
        "单元电流面积投影结果",
        "重构幅值误差（上限显示0.15 decade）",
    )
    for index, (axis, title) in enumerate(zip(axes, titles, strict=True)):
        spatial_format(axis)
        axis.set_title(title)
        add_panel_label(axis, f"({chr(97 + index)})")
    axes[2].legend(loc="lower right", fontsize=8)
    current_bar = figure.colorbar(first, ax=axes[:2], shrink=0.88, pad=0.02)
    current_bar.set_label("log₁₀|空穴电流密度|（A/cm²）")
    error_bar = figure.colorbar(third, ax=axes[2], shrink=0.88, pad=0.02)
    error_bar.set_label("|log₁₀(J重构/J原始)|（decade）")
    figure.suptitle(
        "图1  SDevice单元电流向节点电流的重构验证",
        fontsize=14,
        color=INK,
        fontweight="bold",
    )
    figure.savefig(output, dpi=260, bbox_inches="tight", metadata={"Title": "SDevice element-to-node current reconstruction"})
    plt.close(figure)
    return {
        "all_node_p95_abs_log_error_decade": float(np.quantile(error, 0.95)),
        "tail_p95_abs_log_error_decade": float(np.quantile(error[tail_nodes], 0.95)),
        "tail_max_abs_log_error_decade": float(np.max(error[tail_nodes])),
    }


def group_p95(
    diagnostic_rows: list[dict[str, str]],
    column: str,
    order: list[str],
) -> list[tuple[str, float, int, int]]:
    values: dict[str, list[float]] = defaultdict(list)
    tails: dict[str, int] = defaultdict(int)
    for row in diagnostic_rows:
        if row["selected_by_original_gate_mask"].lower() != "true":
            continue
        key = row[column]
        values[key].append(float(row["absolute_log10_current_error_decade"]))
        tails[key] += int(row["exceeds_0p5_decade"].lower() == "true")
    result = []
    for key in order:
        if key in values:
            result.append((key, float(np.quantile(values[key], 0.95)), tails[key], len(values[key])))
    return result


def plot_group_bars(
    axis: plt.Axes,
    groups: list[tuple[str, float, int, int]],
    labels: dict[str, str],
    title: str,
) -> None:
    y = np.arange(len(groups))
    p95 = np.asarray([row[1] for row in groups])
    colors = [ORANGE if value > 0.5 else BLUE for value in p95]
    axis.barh(y, p95, color=colors, alpha=0.9, height=0.62)
    axis.set_yticks(y, [labels.get(row[0], row[0]) for row in groups])
    axis.invert_yaxis()
    axis.axvline(0.5, color=INK, linestyle="--", linewidth=1.0, label="0.5 decade诊断线")
    axis.set_xlabel("P95幅值误差（decade）")
    axis.set_title(title)
    axis.grid(axis="x", alpha=0.75)
    axis.set_axisbelow(True)
    maximum = max(max(p95, default=0.0), 0.5)
    # Leave enough room for the count annotations used by the report figure.
    axis.set_xlim(0.0, maximum * 1.45 + 0.08)
    for index, (_, value, tail, total) in enumerate(groups):
        axis.text(
            value + maximum * 0.025,
            index,
            f"{tail}/{total}节点超限",
            va="center",
            ha="left",
            fontsize=8,
            color=MUTED,
        )
    axis.legend(loc="upper right", fontsize=8)


def plot_figure2(
    triangulation: mtri.Triangulation,
    sentaurus: np.ndarray,
    vela: np.ndarray,
    tail_nodes: np.ndarray,
    diagnostic_rows: list[dict[str, str]],
    output: Path,
) -> dict[str, float]:
    sentaurus_log = np.log10(np.maximum(magnitude(sentaurus), 1.0e-30))
    vela_log = np.log10(np.maximum(magnitude(vela), 1.0e-30))
    error = abs_log_error(vela, sentaurus)
    figure = plt.figure(figsize=(14.2, 8.8), constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=[1.12, 1.0])
    top = [figure.add_subplot(grid[0, index]) for index in range(3)]
    bottom = [figure.add_subplot(grid[1, index]) for index in range(3)]

    first = top[0].tripcolor(
        triangulation, sentaurus_log, shading="gouraud", cmap="cividis", vmin=-5.0, vmax=2.0, rasterized=True
    )
    top[1].tripcolor(
        triangulation, vela_log, shading="gouraud", cmap="cividis", vmin=-5.0, vmax=2.0, rasterized=True
    )
    third = top[2].tripcolor(
        triangulation,
        np.minimum(error, 2.0),
        shading="gouraud",
        cmap="YlOrBr",
        vmin=0.0,
        vmax=2.0,
        rasterized=True,
    )
    top[2].scatter(
        triangulation.x[tail_nodes],
        triangulation.y[tail_nodes],
        s=9,
        facecolors="none",
        edgecolors=PINK,
        linewidths=0.55,
        zorder=6,
    )
    for index, (axis, title) in enumerate(
        zip(
            top,
            ("SDevice 空穴电流密度", "Vela 当前节点电流密度", "幅值误差与弱电流尾部位置"),
            strict=True,
        )
    ):
        spatial_format(axis)
        axis.set_title(title)
        add_panel_label(axis, f"({chr(97 + index)})")
    current_bar = figure.colorbar(first, ax=top[:2], shrink=0.86, pad=0.02)
    current_bar.set_label("log₁₀|空穴电流密度|（A/cm²）")
    error_bar = figure.colorbar(third, ax=top[2], shrink=0.86, pad=0.02)
    error_bar.set_label("|log₁₀(JVela/JSDevice)|（decade）")

    electrical = group_p95(
        diagnostic_rows,
        "electrical_region",
        ["n_emitter_side", "p_base"],
    )
    depth = group_p95(
        diagnostic_rows,
        "depth_band_um",
        ["[0,0.25)", "[0.25,0.50)", "[0.50,0.75)"],
    )
    current_bin = group_p95(
        diagnostic_rows,
        "reference_magnitude_bin",
        ["[1e-2,1]", "[1e-3,1e-2)", "[1e-4,1e-3)", "[1e-5,1e-4)", "[1e-6,1e-5)"],
    )
    plot_group_bars(
        bottom[0],
        electrical,
        {"n_emitter_side": "发射区侧", "p_base": "P型基区"},
        "按电学区域统计",
    )
    plot_group_bars(
        bottom[1],
        depth,
        {
            "[0,0.25)": "0–0.25 µm",
            "[0.25,0.50)": "0.25–0.50 µm",
            "[0.50,0.75)": "0.50–0.75 µm",
        },
        "按器件深度统计",
    )
    plot_group_bars(
        bottom[2],
        current_bin,
        {
            "[1e-2,1]": "峰值的10⁻²–1",
            "[1e-3,1e-2)": "峰值的10⁻³–10⁻²",
            "[1e-4,1e-3)": "峰值的10⁻⁴–10⁻³",
            "[1e-5,1e-4)": "峰值的10⁻⁵–10⁻⁴",
            "[1e-6,1e-5)": "峰值的10⁻⁶–10⁻⁵",
        },
        "按参考电流强度统计",
    )
    for index, axis in enumerate(bottom, start=3):
        add_panel_label(axis, f"({chr(97 + index)})")
    figure.suptitle(
        "图2  Vela与SDevice空穴电流密度的空间及分区域对比",
        fontsize=14,
        color=INK,
        fontweight="bold",
    )
    figure.savefig(output, dpi=250, bbox_inches="tight", metadata={"Title": "Vela and SDevice hole-current comparison by region"})
    plt.close(figure)
    return {
        "selected_node_count": int(sum(row["selected_by_original_gate_mask"].lower() == "true" for row in diagnostic_rows)),
        "tail_node_count": int(len(tail_nodes)),
        "selected_p95_abs_log_error_decade": float(
            np.quantile(
                [
                    float(row["absolute_log10_current_error_decade"])
                    for row in diagnostic_rows
                    if row["selected_by_original_gate_mask"].lower() == "true"
                ],
                0.95,
            )
        ),
    }


def plot_figure3(
    triangulation: mtri.Triangulation,
    baseline: np.ndarray,
    strict: np.ndarray,
    tail_nodes: np.ndarray,
    output: Path,
) -> dict[str, float]:
    baseline_mag = magnitude(baseline)
    strict_mag = magnitude(strict)
    baseline_log = np.log10(np.maximum(baseline_mag, 1.0e-30))
    strict_log = np.log10(np.maximum(strict_mag, 1.0e-30))
    error = abs_log_error(strict, baseline)
    safe_error = np.maximum(error, 1.0e-16)
    tail_error = np.sort(np.maximum(error[tail_nodes], 1.0e-16))
    cumulative = np.arange(1, len(tail_error) + 1) / len(tail_error)
    tail_p95 = float(np.quantile(error[tail_nodes], 0.95))
    tail_max = float(np.max(error[tail_nodes]))

    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.65), constrained_layout=True)
    axes[0].scatter(
        baseline_log,
        strict_log,
        s=7,
        color=BLUE,
        alpha=0.38,
        linewidths=0,
        rasterized=True,
    )
    lower = min(float(np.min(baseline_log)), float(np.min(strict_log)))
    upper = max(float(np.max(baseline_log)), float(np.max(strict_log)))
    axes[0].plot([lower, upper], [lower, upper], color=INK, linewidth=1.1, linestyle="--", label="完全一致线")
    axes[0].set_xlim(lower, upper)
    axes[0].set_ylim(lower, upper)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel("普通精度 log₁₀|Jh|（A/cm²）")
    axes[0].set_ylabel("Digits=8 + 128位 log₁₀|Jh|（A/cm²）")
    axes[0].set_title("全部5611个节点的一致性")
    axes[0].grid(alpha=0.65)
    axes[0].legend(loc="lower right", fontsize=8)

    spatial = axes[1].tripcolor(
        triangulation,
        safe_error,
        shading="gouraud",
        cmap="YlOrBr",
        norm=LogNorm(vmin=1.0e-14, vmax=1.0e-7),
        rasterized=True,
    )
    spatial_format(axes[1])
    axes[1].set_title("高精度相对普通精度的空间误差")
    colorbar = figure.colorbar(spatial, ax=axes[1], shrink=0.86, pad=0.02)
    colorbar.set_label("|log₁₀(J高精度/J普通)|（decade）")

    axes[2].semilogx(tail_error, cumulative * 100.0, color=BLUE, linewidth=2.0)
    axes[2].axvline(tail_p95, color=ORANGE, linestyle="--", linewidth=1.2, label=f"P95={tail_p95:.2e}")
    axes[2].axvline(tail_max, color=PINK, linestyle=":", linewidth=1.5, label=f"最大值={tail_max:.2e}")
    axes[2].set_xlabel("180个尾部节点的幅值误差（decade，对数轴）")
    axes[2].set_ylabel("累计节点比例（%）")
    axes[2].set_ylim(0.0, 101.0)
    axes[2].set_title("弱电流尾部误差累计分布")
    axes[2].grid(alpha=0.65, which="both")
    axes[2].legend(loc="lower right", fontsize=8)

    for index, axis in enumerate(axes):
        add_panel_label(axis, f"({chr(97 + index)})")
    figure.suptitle(
        "图3  普通精度与高精度SDevice计算结果对比",
        fontsize=14,
        color=INK,
        fontweight="bold",
    )
    figure.savefig(output, dpi=260, bbox_inches="tight", metadata={"Title": "Normal and high-precision SDevice current comparison"})
    plt.close(figure)
    return {
        "tail_p95_abs_log_error_decade": tail_p95,
        "tail_max_abs_log_error_decade": tail_max,
        "all_node_max_abs_log_error_decade": float(np.max(error)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe-root", type=Path, default=PROBE_ROOT)
    parser.add_argument(
        "--tail-diagnostics",
        type=Path,
        default=BUILD_ROOT / "m1_hole_current_tail_diagnosis" / "node_diagnostics.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=FIXTURE / "figures" / "current_semantics")
    args = parser.parse_args()

    configure_style()
    default_export = args.probe_root / "default_export"
    strict_export = args.probe_root / "digits8_ep128_export"
    coordinates, triangles = load_mesh(default_export)
    triangulation = mtri.Triangulation(coordinates[:, 0], coordinates[:, 1], triangles)
    fields = default_export / "fields"
    native = load_vector(fields / "hCurrentDensity_region0.csv", "node_id")
    cell_current = load_vector(fields / "hCurrentDensity_region0_cells.csv", "cell_id")
    strict = load_vector(
        strict_export / "fields" / "hCurrentDensity_region0.csv", "node_id"
    )
    diagnostic_rows = read_rows(args.tail_diagnostics)
    if [int(row["node_id"]) for row in diagnostic_rows] != list(range(len(coordinates))):
        raise ValueError("tail diagnostic rows do not match the mesh node order")
    vela = np.asarray(
        [[float(row["vela_x_A_per_cm2"]), float(row["vela_y_A_per_cm2"])] for row in diagnostic_rows],
        dtype=float,
    )
    tail_nodes = np.asarray(
        [int(row["node_id"]) for row in diagnostic_rows if row["exceeds_0p5_decade"].lower() == "true"],
        dtype=int,
    )
    triangle_ids, areas = triangle_areas(default_export)
    reconstructed = np.asarray(
        reconstruct(
            len(coordinates),
            triangle_ids,
            areas,
            [tuple(value) for value in cell_current],
            "area",
        ),
        dtype=float,
    ) * (1.0e-6 / math.sqrt(2.0))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "figure1": args.output_dir / "figure1_sdevice_element_node_reconstruction.png",
        "figure2": args.output_dir / "figure2_vela_sdevice_hole_current_regions.png",
        "figure3": args.output_dir / "figure3_sdevice_precision_comparison.png",
    }
    metrics = {
        "figure1": plot_figure1(triangulation, native, reconstructed, tail_nodes, outputs["figure1"]),
        "figure2": plot_figure2(triangulation, native, vela, tail_nodes, diagnostic_rows, outputs["figure2"]),
        "figure3": plot_figure3(triangulation, native, strict, tail_nodes, outputs["figure3"]),
    }
    source_paths = {
        "nodes": default_export / "nodes.csv",
        "elements": default_export / "elements.csv",
        "sdevice_node_current": fields / "hCurrentDensity_region0.csv",
        "sdevice_cell_current": fields / "hCurrentDensity_region0_cells.csv",
        "vela_tail_diagnostics": args.tail_diagnostics,
        "sdevice_digits8_ep128_current": strict_export / "fields" / "hCurrentDensity_region0.csv",
    }
    manifest = {
        "schema": "vela.genius_bjt_current_semantics_figures.v1",
        "bias": {"VBE_V": 0.70, "VCE_V": 3.00},
        "mesh": {"nodes": len(coordinates), "triangles": len(triangles)},
        "metrics": metrics,
        "sources_sha256": {name: sha256(path) for name, path in source_paths.items()},
        "outputs_sha256": {name: sha256(path) for name, path in outputs.items()},
    }
    manifest_path = args.output_dir / "current_semantics_figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "outputs": {name: str(path) for name, path in outputs.items()},
                "manifest": str(manifest_path),
                **manifest,
            },
            indent=2,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
