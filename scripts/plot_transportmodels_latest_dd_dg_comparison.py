#!/usr/bin/env python3
"""Regenerate the latest TransportModels DD/DG I-V and spatial comparisons."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm


REPO_ROOT = Path(__file__).resolve().parents[1]
CASE_ROOT = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022"
)
BASELINE_ROOT = CASE_ROOT / "vela_baseline"
LATEST_BASELINE = (
    BASELINE_ROOT / "dd_dg_continuous_contact_basin_kcl_v5_2026-08-24"
)
REFERENCE_CURVES = BASELINE_ROOT / "generated/reference_curves"
REFERENCE_FIELDS = BASELINE_ROOT / "generated/sim_fields"
MESH_PATH = BASELINE_ROOT / "generated/vela/mesh.json"
DEFAULT_OUTPUT = CASE_ROOT / "reports/latest_dd_dg_figures_20260825"

SENTAURUS_COLOR = "#343A40"
VELA_COLOR = "#2878B5"
ERROR_COLOR = "#E07A1F"


@dataclass(frozen=True)
class FieldSpec:
    key: str
    title: str
    sentaurus_file: str
    vela_column: str
    display_label: str
    logarithmic: bool = False
    vela_scale: float = 1.0
    dg_only: bool = False


FIELD_SPECS = (
    FieldSpec(
        "electrostatic_potential",
        "静电势 ψ",
        "ElectrostaticPotential_region3.csv",
        "psi",
        "ψ (V)",
    ),
    FieldSpec(
        "electron_quasi_fermi_potential",
        "电子准费米势 φn",
        "eQuasiFermiPotential_region3.csv",
        "phin",
        "φn (V)",
    ),
    FieldSpec(
        "hole_quasi_fermi_potential",
        "空穴准费米势 φp",
        "hQuasiFermiPotential_region3.csv",
        "phip",
        "φp (V)",
    ),
    FieldSpec(
        "electron_density",
        "电子浓度 n",
        "eDensity_region3.csv",
        "electrons_m3",
        "log10(n / cm⁻³)",
        logarithmic=True,
        vela_scale=1.0e-6,
    ),
    FieldSpec(
        "hole_density",
        "空穴浓度 p",
        "hDensity_region3.csv",
        "holes_m3",
        "log10(p / cm⁻³)",
        logarithmic=True,
        vela_scale=1.0e-6,
    ),
    FieldSpec(
        "electron_quantum_potential",
        "电子量子势 Qn",
        "eQuantumPotential_region3.csv",
        "electron_quantum_potential_V",
        "Qn (V)",
        dg_only=True,
    ),
)


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Microsoft YaHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#56606A",
            "axes.labelcolor": "#252A31",
            "axes.titlecolor": "#20252B",
            "xtick.color": "#4E5661",
            "ytick.color": "#4E5661",
            "grid.color": "#D9DEE5",
            "grid.linewidth": 0.7,
            "grid.alpha": 0.75,
            "legend.frameon": False,
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_curve(path: Path, current_column: str) -> dict[float, float]:
    result: dict[float, float] = {}
    for row in read_rows(path):
        result[float(row["bias_V"])] = abs(float(row[current_column]))
    return result


def aligned_curve(mode: str, sweep: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference = read_curve(
        REFERENCE_CURVES
        / f"transportmodels_sentaurus2022_{mode}_{sweep}_reference.csv",
        "current_total",
    )
    candidate = read_curve(
        LATEST_BASELINE / f"{mode}_{sweep}_curve_comparison_candidate.csv",
        "current_total_A_per_um",
    )
    biases = np.asarray(sorted(reference), dtype=float)
    if len(biases) != 21 or len(candidate) != 21:
        raise RuntimeError(
            f"Expected 21 points for {mode} {sweep}, got "
            f"reference={len(reference)} candidate={len(candidate)}"
        )
    candidate_values = []
    for bias in biases:
        key = min(candidate, key=lambda value: abs(value - bias))
        if abs(key - bias) > 1.0e-10:
            raise RuntimeError(f"Missing candidate bias {bias:g} V")
        candidate_values.append(candidate[key])
    return (
        biases,
        np.asarray([reference[bias] for bias in biases]),
        np.asarray(candidate_values),
    )


def write_curve_csv(
    output_dir: Path,
    mode: str,
    sweep: str,
    bias: np.ndarray,
    reference: np.ndarray,
    candidate: np.ndarray,
    error: np.ndarray,
) -> None:
    path = output_dir / f"transportmodels_{mode}_{sweep}_comparison.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "bias_V",
                "sentaurus_current_A_per_um",
                "vela_current_A_per_um",
                "signed_relative_error_percent",
            ]
        )
        for index in range(len(bias)):
            writer.writerow(
                [
                    f"{bias[index]:.17g}",
                    f"{reference[index]:.17g}",
                    f"{candidate[index]:.17g}",
                    "" if not np.isfinite(error[index]) else f"{error[index]:.17g}",
                ]
            )


def plot_iv_mode(output_dir: Path, mode: str) -> dict[str, object]:
    fig = plt.figure(figsize=(14.4, 8.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[3.2, 1.2])
    result: dict[str, object] = {}
    for column, sweep in enumerate(("idvg", "idvd")):
        bias, reference, candidate = aligned_curve(mode, sweep)
        error = 100.0 * (candidate - reference) / np.abs(reference)
        if sweep == "idvd":
            error[np.isclose(bias, 0.0, atol=1.0e-14)] = np.nan

        curve_ax = fig.add_subplot(grid[0, column])
        error_ax = fig.add_subplot(grid[1, column], sharex=curve_ax)
        plot_method = curve_ax.semilogy if sweep == "idvg" else curve_ax.plot
        plot_method(
            bias,
            reference,
            color=SENTAURUS_COLOR,
            linewidth=2.4,
            marker="o",
            markersize=5.2,
            markerfacecolor="white",
            markeredgewidth=1.3,
            label="Sentaurus 2022",
        )
        plot_method(
            bias,
            candidate,
            color=VELA_COLOR,
            linewidth=2.2,
            linestyle="--",
            marker="s",
            markersize=4.8,
            label=f"Vela {mode.upper()}",
        )
        fixed_bias = "Vd = 1.1 V" if sweep == "idvg" else "Vg = 1.0 V"
        curve_ax.set_title(f"{sweep.replace('idv', 'Id–V').replace('g', 'g').replace('d', 'd')}（{fixed_bias}）")
        curve_ax.set_ylabel("|Id| (A/μm)")
        curve_ax.tick_params(axis="x", labelbottom=False)
        curve_ax.grid(True, which="both")
        curve_ax.legend(loc="best")

        finite = np.isfinite(error)
        error_ax.axhline(0.0, color="#6B747D", linewidth=1.0)
        error_ax.plot(
            bias[finite],
            error[finite],
            color=ERROR_COLOR,
            linewidth=1.9,
            marker="D",
            markersize=4.0,
        )
        error_ax.fill_between(
            bias[finite], error[finite], 0.0, color=ERROR_COLOR, alpha=0.16
        )
        error_ax.set_xlabel("栅极电压 Vg (V)" if sweep == "idvg" else "漏极电压 Vd (V)")
        error_ax.set_ylabel("误差 (%)")
        error_ax.grid(True)
        if sweep == "idvd":
            error_ax.text(
                0.98,
                0.08,
                "Vd=0 时参考电流接近零，百分比误差不定义",
                transform=error_ax.transAxes,
                ha="right",
                fontsize=8.7,
                color="#616A73",
            )

        max_index = int(np.nanargmax(np.abs(error)))
        result[sweep] = {
            "points": int(len(bias)),
            "max_absolute_error_percent": float(abs(error[max_index])),
            "max_error_bias_V": float(bias[max_index]),
            "endpoint_error_percent": float(error[-1]),
        }
        write_curve_csv(output_dir, mode, sweep, bias, reference, candidate, error)

    fig.suptitle(
        f"TransportModels {mode.upper()}：Vela 与 Sentaurus I–V 连续基线对比",
        fontsize=16,
    )
    fig.text(
        0.01,
        -0.015,
        "误差定义：100×(Vela−Sentaurus)/|Sentaurus|；每条曲线均为同一配置、仅由前一点延续的 21 点扫描。",
        fontsize=9.2,
        color="#626A74",
    )
    stem = f"transportmodels_{mode}_iv_comparison"
    for suffix in ("png", "svg"):
        fig.savefig(output_dir / f"{stem}.{suffix}", dpi=240 if suffix == "png" else None)
    plt.close(fig)
    return result


def load_geometry() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mesh = json.loads(MESH_PATH.read_text(encoding="utf-8"))
    nodes = sorted(mesh["nodes"], key=lambda node: int(node["id"]))
    node_ids = [int(node["id"]) for node in nodes]
    if node_ids != list(range(len(nodes))):
        raise RuntimeError("Mesh node IDs must be contiguous")
    lateral = np.asarray([float(node["y"]) for node in nodes])
    depth = np.asarray([float(node["x"]) for node in nodes])
    triangles = np.asarray(
        [
            [int(node_id) for node_id in triangle["node_ids"]]
            for triangle in mesh["triangles"]
            if int(triangle["region_id"]) == 3
        ],
        dtype=int,
    )
    substrate_nodes = np.asarray(sorted(set(triangles.ravel())), dtype=int)
    return lateral, depth, triangles, substrate_nodes


def read_node_field(path: Path, column: str, size: int) -> np.ndarray:
    result = np.full(size, np.nan, dtype=float)
    for row in read_rows(path):
        result[int(row["node_id"])] = float(row[column])
    return result


def bounded_symmetric_percent(candidate: np.ndarray, reference: np.ndarray) -> np.ndarray:
    denominator = np.abs(candidate) + np.abs(reference)
    result = np.zeros_like(denominator)
    np.divide(
        200.0 * (candidate - reference),
        denominator,
        out=result,
        where=denominator > 0.0,
    )
    return result


def write_field_csv(
    output_dir: Path,
    mode: str,
    spec: FieldSpec,
    nodes: np.ndarray,
    lateral: np.ndarray,
    depth: np.ndarray,
    reference: np.ndarray,
    candidate: np.ndarray,
    error: np.ndarray,
) -> None:
    path = output_dir / f"transportmodels_{mode}_{spec.key}_comparison.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "node_id",
                "lateral_um",
                "depth_um",
                "sentaurus_value",
                "vela_value",
                "bounded_symmetric_error_percent",
            ]
        )
        for node in nodes:
            writer.writerow(
                [
                    int(node),
                    f"{lateral[node]:.17g}",
                    f"{depth[node]:.17g}",
                    f"{reference[node]:.17g}",
                    f"{candidate[node]:.17g}",
                    f"{error[node]:.17g}",
                ]
            )


def plot_spatial_field(
    output_dir: Path,
    mode: str,
    spec: FieldSpec,
    triangulation: mtri.Triangulation,
    lateral: np.ndarray,
    depth: np.ndarray,
    substrate_nodes: np.ndarray,
) -> dict[str, object]:
    state_path = LATEST_BASELINE / f"{mode}_idvd_curve_state_bias_2p000000.csv"
    reference_path = REFERENCE_FIELDS / f"{mode}_idvd/fields" / spec.sentaurus_file
    reference = read_node_field(reference_path, "component0", len(lateral))
    candidate = spec.vela_scale * read_node_field(
        state_path, spec.vela_column, len(lateral)
    )
    visible = substrate_nodes[depth[substrate_nodes] <= 0.22]
    if np.any(~np.isfinite(reference[visible])) or np.any(~np.isfinite(candidate[visible])):
        raise RuntimeError(f"Missing {mode} {spec.key} values in the visible silicon region")
    error = bounded_symmetric_percent(candidate, reference)

    if spec.logarithmic:
        reference_display = np.log10(np.maximum(reference, np.finfo(float).tiny))
        candidate_display = np.log10(np.maximum(candidate, np.finfo(float).tiny))
    else:
        reference_display = reference
        candidate_display = candidate
    shared = np.concatenate((reference_display[visible], candidate_display[visible]))
    absolute_norm = Normalize(vmin=float(np.min(shared)), vmax=float(np.max(shared)))
    p99 = float(np.percentile(np.abs(error[visible]), 99.0))
    if p99 <= 0.0:
        p99 = 1.0
    error_norm = TwoSlopeNorm(vmin=-p99, vcenter=0.0, vmax=p99)

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 5.25), constrained_layout=True)
    absolute_map = None
    error_map = None
    for index, (title, values, cmap, norm) in enumerate(
        (
            ("Sentaurus 2022", reference_display, "viridis", absolute_norm),
            (f"Vela {mode.upper()}", candidate_display, "viridis", absolute_norm),
            ("有界对称误差", error, "RdBu_r", error_norm),
        )
    ):
        axis = axes[index]
        mappable = axis.tripcolor(
            triangulation,
            values,
            shading="gouraud",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        axis.triplot(triangulation, color="#5D6670", linewidth=0.1, alpha=0.22)
        axis.axhline(0.0, color="#343A40", linewidth=0.8)
        axis.set_title(title)
        axis.set_xlabel("横向位置 (μm)")
        if index == 0:
            axis.set_ylabel("表面法向深度 (μm)")
        axis.set_xlim(-0.2, 0.2)
        axis.set_ylim(0.22, -0.005)
        axis.set_aspect("equal", adjustable="box")
        if index < 2:
            absolute_map = mappable
        else:
            error_map = mappable
    assert absolute_map is not None and error_map is not None
    cbar = fig.colorbar(absolute_map, ax=axes[:2], location="bottom", shrink=0.82, pad=0.08)
    cbar.set_label(spec.display_label)
    ebar = fig.colorbar(error_map, ax=axes[2], location="bottom", shrink=0.82, pad=0.08)
    ebar.set_label("200×(Vela−Sentaurus)/(|Vela|+|Sentaurus|) (%)")
    fig.suptitle(
        f"TransportModels {mode.upper()} 关键物理量：{spec.title}", fontsize=15.5
    )
    fig.text(
        0.01,
        -0.015,
        f"工作点 Vg=1.0 V、Vd=2.0 V；硅区表面以下 0–0.22 μm；误差色标按 P99={p99:.3g}% 截断。",
        fontsize=9.0,
        color="#626A74",
    )
    stem = f"transportmodels_{mode}_{spec.key}_spatial_comparison"
    for suffix in ("png", "svg"):
        fig.savefig(
            output_dir / f"{stem}.{suffix}",
            dpi=240 if suffix == "png" else None,
        )
    plt.close(fig)
    write_field_csv(
        output_dir,
        mode,
        spec,
        visible,
        lateral,
        depth,
        reference,
        candidate,
        error,
    )
    return {
        "nodes_compared": int(len(visible)),
        "median_absolute_bounded_error_percent": float(np.median(np.abs(error[visible]))),
        "p99_absolute_bounded_error_percent": p99,
        "max_absolute_bounded_error_percent": float(np.max(np.abs(error[visible]))),
        "png": f"{stem}.png",
        "svg": f"{stem}.svg",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()

    summary: dict[str, object] = {
        "baseline": str(LATEST_BASELINE.resolve()),
        "output_dir": str(output_dir),
        "iv_error_definition": "100*(Vela-Sentaurus)/abs(Sentaurus)",
        "spatial_error_definition": "200*(Vela-Sentaurus)/(abs(Vela)+abs(Sentaurus))",
        "iv": {},
        "spatial": {},
    }
    for mode in ("dd", "dg"):
        summary["iv"][mode] = plot_iv_mode(output_dir, mode)  # type: ignore[index]

    lateral, depth, triangles, substrate_nodes = load_geometry()
    triangulation = mtri.Triangulation(lateral, depth, triangles)
    for mode in ("dd", "dg"):
        mode_summary: dict[str, object] = {}
        for spec in FIELD_SPECS:
            if spec.dg_only and mode != "dg":
                continue
            mode_summary[spec.key] = plot_spatial_field(
                output_dir,
                mode,
                spec,
                triangulation,
                lateral,
                depth,
                substrate_nodes,
            )
        summary["spatial"][mode] = mode_summary  # type: ignore[index]

    summary_path = output_dir / "transportmodels_latest_dd_dg_figure_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
