#!/usr/bin/env python3
"""Plot the coarse-grid Genius BJT Sentaurus/Vela comparison suite."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import TwoSlopeNorm


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

from compare_genius_bjt_transport_fields import (  # noqa: E402
    read_vtk_point_data,
    sentaurus_scalar,
    sentaurus_vector,
    vector_metrics,
)
from plot_genius_bjt_sentaurus_vela import (  # noqa: E402
    GRID,
    INK,
    MUTED,
    ORANGE,
    PINK,
    SENTAURUS,
    VELA,
    configure_style,
    load_mesh,
)


FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_VTK = BUILD / "cell_first_recovery_ab" / "coarse" / "vce_030.vtk"
DEFAULT_SENT = BUILD / "m1_p0_multibias" / "sentaurus" / "vce_030"
DEFAULT_OUTPUT = FIXTURE / "figures" / "coarse_m1_comparison"
EPS = 1.0e-300


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtk", type=Path, default=DEFAULT_VTK)
    parser.add_argument("--sentaurus", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def percentile(values: np.ndarray, fraction: float) -> float:
    return float(np.quantile(np.asarray(values, dtype=float), fraction))


def scalar_metrics(
    difference: np.ndarray, mask: np.ndarray | None = None
) -> dict[str, float | int]:
    selected = np.asarray(difference if mask is None else difference[mask], dtype=float)
    selected = selected[np.isfinite(selected)]
    absolute = np.abs(selected)
    return {
        "selected_nodes": int(len(selected)),
        "rmse": float(math.sqrt(np.mean(selected * selected))),
        "p95_absolute_error": percentile(absolute, 0.95),
        "maximum_absolute_error": float(np.max(absolute)),
    }


def spatial_format(axis: plt.Axes) -> None:
    axis.set_aspect("equal")
    axis.set_xlim(-0.03, 6.03)
    axis.set_ylim(2.03, -0.03)
    axis.set_xlabel("x [µm]")
    axis.set_ylabel("y [µm] (depth)")


def robust_limit(values: np.ndarray, mask: np.ndarray | None) -> float:
    selected = values if mask is None else values[mask]
    selected = np.abs(selected[np.isfinite(selected)])
    if len(selected) == 0:
        return 1.0
    return max(float(np.quantile(selected, 0.995)), np.finfo(float).eps)


def plot_scalar_group(
    triangulation: mtri.Triangulation,
    rows: list[dict[str, object]],
    title: str,
    output: Path,
) -> dict[str, dict[str, float | int]]:
    figure, axes = plt.subplots(
        len(rows), 3, figsize=(15.8, 3.55 * len(rows)), constrained_layout=True
    )
    axes = np.atleast_2d(axes)
    metrics: dict[str, dict[str, float | int]] = {}
    for row_index, spec in enumerate(rows):
        name = str(spec["name"])
        sent = np.asarray(spec["sentaurus"], dtype=float)
        vela = np.asarray(spec["vela"], dtype=float)
        mask = spec.get("mask")
        mask_array = None if mask is None else np.asarray(mask, dtype=bool)
        difference = (vela - sent) * float(spec.get("difference_scale", 1.0))
        metrics[name] = scalar_metrics(difference, mask_array)

        display_sent = np.ma.masked_where(~mask_array, sent) if mask_array is not None else sent
        display_vela = np.ma.masked_where(~mask_array, vela) if mask_array is not None else vela
        display_diff = (
            np.ma.masked_where(~mask_array, difference)
            if mask_array is not None
            else difference
        )
        combined = np.concatenate(
            [np.asarray(display_sent.compressed() if np.ma.isMaskedArray(display_sent) else display_sent),
             np.asarray(display_vela.compressed() if np.ma.isMaskedArray(display_vela) else display_vela)]
        )
        vmin = float(np.min(combined))
        vmax = float(np.max(combined))
        limit = robust_limit(difference, mask_array)
        main0 = axes[row_index, 0].tripcolor(
            triangulation,
            display_sent,
            shading="gouraud",
            cmap=str(spec["cmap"]),
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        axes[row_index, 1].tripcolor(
            triangulation,
            display_vela,
            shading="gouraud",
            cmap=str(spec["cmap"]),
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        diff_image = axes[row_index, 2].tripcolor(
            triangulation,
            display_diff,
            shading="gouraud",
            cmap="coolwarm",
            norm=TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit),
            rasterized=True,
        )
        for axis in axes[row_index, :]:
            spatial_format(axis)
        axes[row_index, 0].set_title(f"SDevice · {name}")
        axes[row_index, 1].set_title(f"Vela · {name}")
        axes[row_index, 2].set_title(
            f"Vela − SDevice · P95={metrics[name]['p95_absolute_error']:.3g} "
            f"{spec['difference_unit']}"
        )
        main_bar = figure.colorbar(
            main0, ax=axes[row_index, :2], shrink=0.83, pad=0.018
        )
        main_bar.set_label(str(spec["unit"]))
        diff_bar = figure.colorbar(
            diff_image, ax=axes[row_index, 2], shrink=0.83, pad=0.018
        )
        diff_bar.set_label(
            f"{spec['difference_unit']} (color clipped at 99.5th percentile)"
        )
    figure.suptitle(title, fontsize=14, color=INK)
    figure.savefig(output, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    return metrics


def interpolated_profile(
    triangulation: mtri.Triangulation, values: np.ndarray, x_um: float, y_um: np.ndarray
) -> np.ndarray:
    interpolator = mtri.LinearTriInterpolator(triangulation, values)
    result = interpolator(np.full_like(y_um, x_um), y_um)
    return np.asarray(np.ma.filled(result, np.nan), dtype=float)


def plot_centerline_profiles(
    triangulation: mtri.Triangulation,
    fields: dict[str, tuple[np.ndarray, np.ndarray, str]],
    output: Path,
) -> None:
    y = np.linspace(0.0, 2.0, 801)
    figure, axes = plt.subplots(3, 2, figsize=(12.8, 11.2), constrained_layout=True)
    order = [
        ("potential", "Electrostatic potential", "V"),
        ("electron_density", "Electron density", "log₁₀(cm⁻³)"),
        ("hole_density", "Hole density", "log₁₀(cm⁻³)"),
        ("electron_qf", "Electron quasi-Fermi potential", "V"),
        ("hole_qf", "Hole quasi-Fermi potential", "V"),
    ]
    for axis, (key, title, unit) in zip(axes.flat[:5], order, strict=True):
        sent, vela, _ = fields[key]
        axis.plot(
            interpolated_profile(triangulation, sent, 3.0, y),
            y,
            color=SENTAURUS,
            linewidth=1.8,
            label="SDevice",
        )
        axis.plot(
            interpolated_profile(triangulation, vela, 3.0, y),
            y,
            color=VELA,
            linewidth=1.8,
            linestyle="--",
            label="Vela",
        )
        axis.set_title(title)
        axis.set_xlabel(unit)
        axis.set_ylabel("y [µm] (depth)")
        axis.set_ylim(2.0, 0.0)
        axis.grid(True, alpha=0.75)

    axis = axes.flat[5]
    for key, label, color in (
        ("electron_current", "electron", VELA),
        ("hole_current", "hole", ORANGE),
    ):
        sent, vela, _ = fields[key]
        axis.plot(
            interpolated_profile(triangulation, sent, 3.0, y),
            y,
            color=color,
            linewidth=1.8,
            label=f"SDevice {label}",
        )
        axis.plot(
            interpolated_profile(triangulation, vela, 3.0, y),
            y,
            color=color,
            linewidth=1.8,
            linestyle="--",
            label=f"Vela {label}",
        )
    axis.set_title("Current-density magnitude")
    axis.set_xlabel("log₁₀(A/cm²)")
    axis.set_ylabel("y [µm] (depth)")
    axis.set_ylim(2.0, 0.0)
    axis.grid(True, alpha=0.75)
    axis.legend(fontsize=8)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2)
    figure.suptitle(
        "Coarse-grid centerline profiles · x=3.00 µm · VBE=0.70 V, VCE=3.00 V",
        fontsize=14,
        color=INK,
    )
    figure.savefig(output, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def plot_terminal_curves(path: Path, output: Path) -> None:
    rows = read_rows(path)
    voltage = np.asarray([float(row["VCE_V"]) for row in rows])
    figure, axes = plt.subplots(2, 2, figsize=(12.8, 8.6), constrained_layout=True)
    observables = [
        ("Collector current |Ic|", "Ic", "A/µm", True),
        ("Base current |Ib|", "Ib", "A/µm", True),
        ("Emitter current |Ie|", "Ie", "A/µm", True),
        ("Current gain β=|Ic/Ib|", "beta", "", False),
    ]
    for axis, (title, key, unit, logarithmic) in zip(
        axes.flat, observables, strict=True
    ):
        sent = np.abs(np.asarray([float(row[f"sentaurus_{key}_{'abs' if key == 'beta' else 'A_per_um'}"]) for row in rows]))
        vela = np.abs(np.asarray([float(row[f"vela_{key}_{'abs' if key == 'beta' else 'A_per_um'}"]) for row in rows]))
        axis.plot(voltage, sent, color=SENTAURUS, linewidth=1.8, marker="o", markersize=3.6, markerfacecolor="white", label="SDevice")
        axis.plot(voltage, vela, color=VELA, linewidth=1.8, linestyle="--", marker="s", markersize=3.3, markerfacecolor="white", label="Vela")
        if logarithmic:
            axis.set_yscale("log")
        axis.set_xlim(0.0, 3.0)
        axis.set_xlabel("VCE [V]")
        axis.set_ylabel(unit)
        axis.set_title(title)
        axis.grid(True, which="major", alpha=0.8)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2)
    figure.suptitle(
        "Genius NPN BJT M1 terminal curves · coarse mesh · VBE=0.70 V",
        fontsize=14,
        color=INK,
    )
    figure.savefig(output, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def plot_current_p95(path: Path, output: Path) -> None:
    report = json.loads(path.read_text(encoding="utf-8"))
    rows = [row for row in report["rows"] if row["grid"] == "coarse"]
    figure, axes = plt.subplots(1, 2, figsize=(12.8, 4.7), constrained_layout=True)
    for axis, carrier, gate in zip(axes, ("electron", "hole"), (0.3, 0.5), strict=True):
        for recovery, label, color, marker in (
            ("direct", "Direct nodal fit", SENTAURUS, "o"),
            ("cell_first", "Cell first → node", VELA, "s"),
        ):
            selected = sorted(
                [row for row in rows if row["carrier"] == carrier and row["recovery"] == recovery],
                key=lambda row: row["VCE_V"],
            )
            axis.plot(
                [row["VCE_V"] for row in selected],
                [row["p95_error_decade"] for row in selected],
                color=color,
                linewidth=1.9,
                marker=marker,
                markersize=5,
                markerfacecolor="white",
                label=label,
            )
        axis.axhline(gate, color=PINK, linestyle=":", linewidth=1.6, label=f"Gate {gate:g} decade")
        axis.set_xlim(0.0, 3.0)
        axis.set_xticks([0, 1, 2, 3])
        axis.set_ylim(bottom=0.0)
        axis.set_xlabel("VCE [V]")
        axis.set_ylabel("P95 |Δlog₁₀| [decade]")
        axis.set_title(f"{carrier.capitalize()} current density")
        axis.grid(True, alpha=0.8)
        axis.legend(fontsize=8)
    figure.suptitle(
        "Coarse-grid current-vector recovery A/B · SDevice reference mask",
        fontsize=14,
        color=INK,
    )
    figure.savefig(output, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()

    mesh_path = FIXTURE / "vela" / "input" / "mesh.json"
    coordinates, triangles, _ = load_mesh(mesh_path)
    triangulation = mtri.Triangulation(
        coordinates[:, 0], coordinates[:, 1], triangles
    )
    count, scalars, vectors = read_vtk_point_data(args.vtk.resolve())
    if count != len(coordinates):
        raise ValueError("VTK and coarse mesh node counts differ")
    sent_fields = args.sentaurus.resolve() / "fields"

    sent_potential = np.asarray(sentaurus_scalar(sent_fields / "ElectrostaticPotential_region0.csv"))
    sent_e = np.asarray(sentaurus_scalar(sent_fields / "eDensity_region0.csv"))
    sent_h = np.asarray(sentaurus_scalar(sent_fields / "hDensity_region0.csv"))
    sent_eqf = np.asarray(sentaurus_scalar(sent_fields / "eQuasiFermiPotential_region0.csv"))
    sent_hqf = np.asarray(sentaurus_scalar(sent_fields / "hQuasiFermiPotential_region0.csv"))
    sent_je_vector = np.asarray(sentaurus_vector(sent_fields / "eCurrentDensity_region0.csv"))
    sent_jh_vector = np.asarray(sentaurus_vector(sent_fields / "hCurrentDensity_region0.csv"))

    vela_potential = np.asarray(scalars["Potential"])
    vela_e = np.asarray(scalars["Electrons"])
    vela_h = np.asarray(scalars["Holes"])
    vela_eqf = np.asarray(scalars["ElectronQuasiFermi"])
    vela_hqf = np.asarray(scalars["HoleQuasiFermi"])
    vela_je_vector = np.asarray(vectors["CellFirstSgElectronCurrentDensityVector"])
    vela_jh_vector = np.asarray(vectors["CellFirstSgHoleCurrentDensityVector"])
    sent_je = np.linalg.norm(sent_je_vector, axis=1)
    sent_jh = np.linalg.norm(sent_jh_vector, axis=1)
    vela_je = np.linalg.norm(vela_je_vector[:, :2], axis=1)
    vela_jh = np.linalg.norm(vela_jh_vector[:, :2], axis=1)
    mask_je = sent_je >= np.max(sent_je) * 1.0e-6
    mask_jh = sent_jh >= np.max(sent_jh) * 1.0e-6

    fields = {
        "potential": (sent_potential, vela_potential, "V"),
        "electron_density": (np.log10(np.maximum(sent_e, EPS)), np.log10(np.maximum(vela_e, EPS)), "log10_cm3"),
        "hole_density": (np.log10(np.maximum(sent_h, EPS)), np.log10(np.maximum(vela_h, EPS)), "log10_cm3"),
        "electron_qf": (sent_eqf, vela_eqf, "V"),
        "hole_qf": (sent_hqf, vela_hqf, "V"),
        "electron_current": (np.log10(np.maximum(sent_je, EPS)), np.log10(np.maximum(vela_je, EPS)), "log10_A_cm2"),
        "hole_current": (np.log10(np.maximum(sent_jh, EPS)), np.log10(np.maximum(vela_jh, EPS)), "log10_A_cm2"),
    }
    outputs = {
        "potential": output_dir / "coarse_vce3_potential_comparison.png",
        "density": output_dir / "coarse_vce3_carrier_density_comparison.png",
        "quasi_fermi": output_dir / "coarse_vce3_quasi_fermi_comparison.png",
        "current_density": output_dir / "coarse_vce3_current_density_comparison.png",
        "centerline": output_dir / "coarse_vce3_centerline_profiles.png",
        "terminal_curves": output_dir / "coarse_terminal_curve_comparison.png",
        "current_p95": output_dir / "coarse_current_recovery_p95_curve.png",
    }
    metrics: dict[str, object] = {}
    metrics["potential"] = plot_scalar_group(
        triangulation,
        [{"name": "Electrostatic potential", "sentaurus": sent_potential, "vela": vela_potential, "unit": "V", "difference_scale": 1.0e3, "difference_unit": "mV", "cmap": "viridis"}],
        "Coarse-grid potential comparison · VBE=0.70 V, VCE=3.00 V",
        outputs["potential"],
    )
    metrics["density"] = plot_scalar_group(
        triangulation,
        [
            {"name": "log₁₀ electron density", "sentaurus": fields["electron_density"][0], "vela": fields["electron_density"][1], "unit": "log₁₀(cm⁻³)", "difference_unit": "decade", "cmap": "cividis"},
            {"name": "log₁₀ hole density", "sentaurus": fields["hole_density"][0], "vela": fields["hole_density"][1], "unit": "log₁₀(cm⁻³)", "difference_unit": "decade", "cmap": "magma"},
        ],
        "Coarse-grid carrier-density comparison · VBE=0.70 V, VCE=3.00 V",
        outputs["density"],
    )
    metrics["quasi_fermi"] = plot_scalar_group(
        triangulation,
        [
            {"name": "Electron quasi-Fermi potential", "sentaurus": sent_eqf, "vela": vela_eqf, "unit": "V", "difference_scale": 1.0e3, "difference_unit": "mV", "cmap": "viridis"},
            {"name": "Hole quasi-Fermi potential", "sentaurus": sent_hqf, "vela": vela_hqf, "unit": "V", "difference_scale": 1.0e3, "difference_unit": "mV", "cmap": "viridis"},
        ],
        "Coarse-grid quasi-Fermi comparison · VBE=0.70 V, VCE=3.00 V",
        outputs["quasi_fermi"],
    )
    metrics["current_density"] = plot_scalar_group(
        triangulation,
        [
            {"name": "log₁₀ electron current magnitude", "sentaurus": fields["electron_current"][0], "vela": fields["electron_current"][1], "mask": mask_je, "unit": "log₁₀(A/cm²)", "difference_unit": "decade", "cmap": "plasma"},
            {"name": "log₁₀ hole current magnitude", "sentaurus": fields["hole_current"][0], "vela": fields["hole_current"][1], "mask": mask_jh, "unit": "log₁₀(A/cm²)", "difference_unit": "decade", "cmap": "plasma"},
        ],
        "Coarse-grid current-density comparison · cell-first Vela recovery · VCE=3.00 V",
        outputs["current_density"],
    )
    plot_centerline_profiles(triangulation, fields, outputs["centerline"])
    terminal_path = FIXTURE / "comparison" / "m1_sentaurus_vela.csv"
    ab_path = FIXTURE / "reports" / "cell_first_recovery_ab.json"
    plot_terminal_curves(terminal_path, outputs["terminal_curves"])
    plot_current_p95(ab_path, outputs["current_p95"])

    metrics["current_vector_gate"] = {
        "electron": vector_metrics(sent_je_vector.tolist(), vela_je_vector.tolist(), 1.0e-6),
        "hole": vector_metrics(sent_jh_vector.tolist(), vela_jh_vector.tolist(), 1.0e-6),
    }
    source_paths = {
        "mesh": mesh_path,
        "vela_vtk": args.vtk.resolve(),
        "sentaurus_field_manifest": args.sentaurus.resolve() / "field_manifest.json",
        "terminal_comparison": terminal_path,
        "current_recovery_ab": ab_path,
    }
    manifest = {
        "schema_version": 1,
        "scope": {"mesh": "coarse", "nodes": count, "triangles": len(triangles), "VBE_V": 0.7, "spatial_VCE_V": 3.0},
        "current_field": "CellFirstSg*CurrentDensityVector (opt-in diagnostic)",
        "current_mask": "SDevice magnitude >= 1e-6 of per-carrier peak",
        "metrics": metrics,
        "sources_sha256": {name: sha256(path) for name, path in source_paths.items()},
        "outputs_sha256": {name: sha256(path) for name, path in outputs.items()},
    }
    manifest_path = output_dir / "figure_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({"outputs": {name: str(path) for name, path in outputs.items()}, "manifest": str(manifest_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
