#!/usr/bin/env python3
"""Render reproducible Genius NPN BJT mesh, field, and curve comparisons."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import SymLogNorm, TwoSlopeNorm


REPO = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_RUN = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "genius_bjt_sentaurus2022"
    / "vela_wp3_wp5"
)
DEFAULT_SENT_FIELDS = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "genius_bjt_sentaurus2022"
    / "m1_current_diagnosis"
    / "sentaurus_vce3"
)

INK = "#27313A"
MUTED = "#66717C"
GRID = "#D9DEE3"
SENTAURUS = "#27313A"
VELA = "#2D6CDF"
ORANGE = "#D97706"
OLIVE = "#6B7D2A"
PINK = "#C24E70"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture-root", type=Path, default=DEFAULT_FIXTURE)
    parser.add_argument("--vela-run-root", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--sentaurus-fields-root", type=Path, default=DEFAULT_SENT_FIELDS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIXTURE / "figures")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def load_mesh(path: Path) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    document = json.loads(path.read_text(encoding="utf-8"))
    nodes = sorted(document["nodes"], key=lambda item: int(item["id"]))
    if [int(item["id"]) for item in nodes] != list(range(len(nodes))):
        raise ValueError("mesh node ids must be contiguous and zero-based")
    coordinates = np.asarray([[float(node["x"]), float(node["y"])] for node in nodes])
    triangles = np.asarray(
        [triangle["node_ids"] for triangle in document["triangles"]], dtype=int
    )
    contacts = {
        contact["name"].lower(): np.asarray(contact["node_ids"], dtype=int)
        for contact in document["contacts"]
    }
    return coordinates, triangles, contacts


def load_node_values(path: Path, count: int, column: str | None = None) -> np.ndarray:
    rows = read_csv(path)
    if column is None:
        columns = [name for name in rows[0] if name != "node_id"]
        if len(columns) != 1:
            raise ValueError(f"cannot infer value column in {path}")
        column = columns[0]
    values = np.full(count, np.nan)
    for row in rows:
        values[int(row["node_id"])] = float(row[column])
    if not np.all(np.isfinite(values)):
        raise ValueError(f"missing or non-finite node values in {path}")
    return values


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
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
        }
    )


def contact_xy(
    coordinates: np.ndarray, contacts: dict[str, np.ndarray], name: str
) -> tuple[np.ndarray, np.ndarray]:
    points = coordinates[contacts[name]]
    order = np.argsort(points[:, 0])
    return points[order, 0], points[order, 1]


def plot_mesh(
    coordinates: np.ndarray,
    triangles: np.ndarray,
    contacts: dict[str, np.ndarray],
    net_doping_cm3: np.ndarray,
    output: Path,
) -> None:
    triangulation = mtri.Triangulation(coordinates[:, 0], coordinates[:, 1], triangles)
    figure, axes = plt.subplots(1, 2, figsize=(14.2, 5.1), constrained_layout=True)
    norm = SymLogNorm(linthresh=1.0e15, linscale=0.7, vmin=-1.0e20, vmax=1.0e20)
    contact_colors = {"base": ORANGE, "emitter": PINK, "collector": OLIVE}

    for axis, zoom in zip(axes, (False, True), strict=True):
        field = axis.tripcolor(
            triangulation,
            net_doping_cm3,
            shading="gouraud",
            cmap="RdBu_r",
            norm=norm,
            rasterized=True,
        )
        axis.triplot(triangulation, color=INK, linewidth=0.22, alpha=0.26)
        for name in ("base", "emitter", "collector"):
            x, y = contact_xy(coordinates, contacts, name)
            axis.plot(
                x,
                y,
                color=contact_colors[name],
                linewidth=3.0,
                solid_capstyle="round",
                label=name.capitalize(),
                zorder=5,
            )
        axis.set_aspect("equal")
        axis.set_xlabel("x [µm]")
        axis.set_ylabel("y [µm] (depth)")
        axis.invert_yaxis()
        if zoom:
            axis.set_xlim(1.0, 5.0)
            axis.set_ylim(0.62, -0.03)
            axis.set_title("Top junction and contact mesh detail")
        else:
            axis.set_xlim(-0.03, 6.03)
            axis.set_ylim(2.03, -0.03)
            axis.set_title("Full device mesh and net doping")
        axis.legend(loc="lower right", ncol=3, fontsize=8)

    colorbar = figure.colorbar(field, ax=axes, location="bottom", shrink=0.72, pad=0.08)
    colorbar.set_ticks([-1.0e19, -1.0e17, -1.0e15, 0.0, 1.0e15, 1.0e17, 1.0e19])
    colorbar.set_ticklabels(["−10¹⁹", "−10¹⁷", "−10¹⁵", "0", "10¹⁵", "10¹⁷", "10¹⁹"])
    colorbar.minorticks_off()
    colorbar.set_label("Net doping Nd − Na [cm⁻³], symmetric-log scale")
    figure.suptitle(
        f"Genius NPN BJT mesh · {len(coordinates)} nodes · {len(triangles)} triangles",
        fontsize=14,
        color=INK,
    )
    figure.savefig(output, dpi=240, bbox_inches="tight", metadata={"Title": "Genius NPN BJT mesh"})
    plt.close(figure)


def add_spatial_axis_format(axis: plt.Axes) -> None:
    axis.set_aspect("equal")
    axis.set_xlim(-0.03, 6.03)
    axis.set_ylim(2.03, -0.03)
    axis.set_xlabel("x [µm]")
    axis.set_ylabel("y [µm]")


def plot_field_comparison(
    coordinates: np.ndarray,
    triangles: np.ndarray,
    sentaurus: dict[str, np.ndarray],
    vela: dict[str, np.ndarray],
    output: Path,
) -> dict[str, dict[str, float]]:
    triangulation = mtri.Triangulation(coordinates[:, 0], coordinates[:, 1], triangles)
    fields = [
        ("Potential", sentaurus["potential_V"], vela["potential_V"], "V", "viridis", 1.0e3, "mV"),
        ("log₁₀ electron density", sentaurus["log_e_cm3"], vela["log_e_cm3"], "log₁₀(cm⁻³)", "cividis", 1.0, "dex"),
        ("log₁₀ hole density", sentaurus["log_h_cm3"], vela["log_h_cm3"], "log₁₀(cm⁻³)", "magma", 1.0, "dex"),
    ]
    figure, axes = plt.subplots(3, 3, figsize=(15.5, 10.0), constrained_layout=True)
    metrics: dict[str, dict[str, float]] = {}

    for row, (name, sent, vel, unit, cmap, diff_scale, diff_unit) in enumerate(fields):
        combined = np.concatenate([sent, vel])
        vmin = float(np.min(combined))
        vmax = float(np.max(combined))
        difference = (vel - sent) * diff_scale
        abs_limit = max(float(np.max(np.abs(difference))), np.finfo(float).eps)
        metrics[name] = {
            "rmse": float(math.sqrt(np.mean(difference * difference))),
            "maximum_absolute_error": float(np.max(np.abs(difference))),
            "difference_unit_scale": diff_scale,
        }

        images = [
            axes[row, 0].tripcolor(
                triangulation, sent, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True
            ),
            axes[row, 1].tripcolor(
                triangulation, vel, shading="gouraud", cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True
            ),
            axes[row, 2].tripcolor(
                triangulation,
                difference,
                shading="gouraud",
                cmap="coolwarm",
                norm=TwoSlopeNorm(vcenter=0.0, vmin=-abs_limit, vmax=abs_limit),
                rasterized=True,
            ),
        ]
        for axis in axes[row, :]:
            add_spatial_axis_format(axis)
        axes[row, 0].set_title(f"Sentaurus · {name}")
        axes[row, 1].set_title(f"Vela · {name}")
        axes[row, 2].set_title(
            f"Vela − Sentaurus · max |Δ|={np.max(np.abs(difference)):.3g} {diff_unit}"
        )
        for column in (0, 1):
            colorbar = figure.colorbar(images[column], ax=axes[row, column], shrink=0.82, pad=0.02)
            colorbar.set_label(unit)
        difference_bar = figure.colorbar(images[2], ax=axes[row, 2], shrink=0.82, pad=0.02)
        difference_bar.set_label(diff_unit)

    figure.suptitle(
        "M1 spatial-field comparison at VBE=0.70 V, VCE=3.00 V",
        fontsize=14,
        color=INK,
    )
    figure.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
        metadata={"Title": "Genius BJT M1 spatial field comparison"},
    )
    plt.close(figure)
    return metrics


def comparison_arrays(path: Path) -> dict[str, np.ndarray]:
    rows = read_csv(path)
    columns = (
        "VCE_V",
        "sentaurus_Ic_A_per_um",
        "vela_Ic_A_per_um",
        "sentaurus_Ib_A_per_um",
        "vela_Ib_A_per_um",
        "sentaurus_beta_abs",
        "vela_beta_abs",
        "Ic_absolute_log10_error",
        "Ib_absolute_log10_error",
        "beta_absolute_log10_error",
    )
    return {name: np.asarray([float(row[name]) for row in rows]) for name in columns}


def plot_curves(m0: dict[str, np.ndarray], m1: dict[str, np.ndarray], output: Path) -> None:
    figure, axes = plt.subplots(2, 3, figsize=(14.8, 8.0), sharex=True, constrained_layout=True)
    observables = [
        ("Ic", "Collector current magnitude [A/µm]", "sentaurus_Ic_A_per_um", "vela_Ic_A_per_um", True),
        ("Ib", "Base current magnitude [A/µm]", "sentaurus_Ib_A_per_um", "vela_Ib_A_per_um", True),
        ("β", "Common-emitter gain |Ic/Ib|", "sentaurus_beta_abs", "vela_beta_abs", False),
    ]

    for row, (model, data) in enumerate((("M0", m0), ("M1", m1))):
        for column, (symbol, ylabel, sent_key, vela_key, logarithmic) in enumerate(observables):
            axis = axes[row, column]
            x = data["VCE_V"]
            axis.axvspan(0.5, 3.0, color=VELA, alpha=0.045, linewidth=0)
            axis.plot(
                x,
                np.abs(data[sent_key]),
                color=SENTAURUS,
                linewidth=1.7,
                marker="o",
                markersize=3.5,
                markerfacecolor="white",
                markeredgewidth=0.9,
                label="Sentaurus",
            )
            axis.plot(
                x,
                np.abs(data[vela_key]),
                color=VELA,
                linewidth=1.7,
                linestyle="--",
                marker="s",
                markersize=3.2,
                markerfacecolor="white",
                markeredgewidth=0.9,
                label="Vela",
            )
            if logarithmic:
                axis.set_yscale("log")
            axis.grid(True, which="major", alpha=0.8)
            axis.set_xlim(0.0, 3.0)
            axis.set_xlabel("VCE [V]")
            axis.set_ylabel(ylabel)
            axis.set_title(f"{model} · {symbol}–VCE")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.035))
    figure.suptitle(
        "Genius NPN BJT terminal-curve comparison · VBE=0.70 V",
        fontsize=14,
        color=INK,
        y=1.035,
    )
    figure.savefig(
        output,
        dpi=240,
        bbox_inches="tight",
        metadata={"Title": "Genius BJT current and gain comparison"},
    )
    plt.close(figure)


def plot_m1_parity_error(m1: dict[str, np.ndarray], output: Path) -> None:
    mask = (m1["VCE_V"] >= 0.5) & (m1["VCE_V"] <= 3.0)
    figure, axis = plt.subplots(figsize=(10.5, 5.2), constrained_layout=True)
    series = [
        ("Ic", m1["Ic_absolute_log10_error"], VELA, "o"),
        ("Ib", m1["Ib_absolute_log10_error"], ORANGE, "s"),
        ("β", m1["beta_absolute_log10_error"], OLIVE, "^"),
    ]
    for name, values, color, marker in series:
        axis.plot(
            m1["VCE_V"][mask],
            values[mask],
            color=color,
            linewidth=1.8,
            marker=marker,
            markersize=4.2,
            markerfacecolor="white",
            label=name,
        )
    axis.axhline(0.05, color=INK, linestyle="--", linewidth=1.4, label="Gate: 0.05 decade")
    axis.set_xlim(0.5, 3.0)
    axis.set_ylim(0.0, 0.055)
    axis.set_xlabel("VCE [V]")
    axis.set_ylabel("Absolute log₁₀ magnitude error [decade]")
    axis.set_title("M1 numerical-parity errors in the asserted active region")
    axis.grid(True, alpha=0.85)
    axis.legend(ncol=4, loc="upper center")
    figure.savefig(
        output,
        dpi=240,
        bbox_inches="tight",
        metadata={"Title": "Genius BJT M1 numerical parity error"},
    )
    plt.close(figure)


def main() -> int:
    args = parse_args()
    fixture = args.fixture_root.resolve()
    run_root = args.vela_run_root.resolve()
    sent_root = args.sentaurus_fields_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh_path = fixture / "vela" / "input" / "mesh.json"
    doping_path = fixture / "vela" / "input" / "doping.csv"
    state_path = run_root / "m1_vce300_state.csv"
    m0_comparison_path = fixture / "comparison" / "m0_sentaurus_vela.csv"
    m1_comparison_path = fixture / "comparison" / "m1_sentaurus_vela.csv"
    sent_nodes_path = sent_root / "nodes.csv"

    coordinates, triangles, contacts = load_mesh(mesh_path)
    count = len(coordinates)
    donor = load_node_values(doping_path, count, "donors_cm3")
    acceptor = load_node_values(doping_path, count, "acceptors_cm3")

    state_rows = read_csv(state_path)
    state = {
        name: load_node_values(state_path, count, name)
        for name in ("psi", "electrons_m3", "holes_m3")
    }
    if len(state_rows) != count:
        raise ValueError("Vela state row count does not match the common mesh")

    sent_nodes = read_csv(sent_nodes_path)
    sent_coordinates = np.full_like(coordinates, np.nan)
    for row in sent_nodes:
        sent_coordinates[int(row["id"])] = [float(row["x_um"]), float(row["y_um"])]
    coordinate_error = float(np.max(np.abs(sent_coordinates - coordinates)))
    if coordinate_error > 1.0e-12:
        raise ValueError(f"Sentaurus/Vela node-coordinate mismatch: {coordinate_error:g} µm")

    sentaurus = {
        "potential_V": load_node_values(
            sent_root / "fields" / "ElectrostaticPotential_region0.csv", count
        ),
        "log_e_cm3": np.log10(
            np.maximum(
                load_node_values(sent_root / "fields" / "eDensity_region0.csv", count),
                1.0e-300,
            )
        ),
        "log_h_cm3": np.log10(
            np.maximum(
                load_node_values(sent_root / "fields" / "hDensity_region0.csv", count),
                1.0e-300,
            )
        ),
    }
    vela = {
        "potential_V": state["psi"],
        "log_e_cm3": np.log10(np.maximum(state["electrons_m3"] / 1.0e6, 1.0e-300)),
        "log_h_cm3": np.log10(np.maximum(state["holes_m3"] / 1.0e6, 1.0e-300)),
    }

    outputs = {
        "mesh": output_dir / "genius_bjt_device_mesh.png",
        "fields": output_dir / "genius_bjt_m1_vce3_field_comparison.png",
        "curves": output_dir / "genius_bjt_terminal_curve_comparison.png",
        "parity": output_dir / "genius_bjt_m1_parity_error.png",
    }
    plot_mesh(coordinates, triangles, contacts, donor - acceptor, outputs["mesh"])
    spatial_metrics = plot_field_comparison(
        coordinates, triangles, sentaurus, vela, outputs["fields"]
    )
    m0 = comparison_arrays(m0_comparison_path)
    m1 = comparison_arrays(m1_comparison_path)
    plot_curves(m0, m1, outputs["curves"])
    plot_m1_parity_error(m1, outputs["parity"])

    source_paths = {
        "mesh": mesh_path,
        "doping": doping_path,
        "vela_m1_vce3_state": state_path,
        "sentaurus_m1_vce3_nodes": sent_nodes_path,
        "sentaurus_potential": sent_root / "fields" / "ElectrostaticPotential_region0.csv",
        "sentaurus_electron_density": sent_root / "fields" / "eDensity_region0.csv",
        "sentaurus_hole_density": sent_root / "fields" / "hDensity_region0.csv",
        "m0_comparison": m0_comparison_path,
        "m1_comparison": m1_comparison_path,
    }
    manifest = {
        "schema_version": 1,
        "device": "Genius NPN BJT",
        "bias": {"VBE_V": 0.7, "VCE_V": 3.0},
        "mesh": {"nodes": count, "triangles": len(triangles)},
        "maximum_sentaurus_vela_coordinate_error_um": coordinate_error,
        "spatial_error_metrics": spatial_metrics,
        "sources_sha256": {name: sha256(path) for name, path in source_paths.items()},
        "outputs_sha256": {name: sha256(path) for name, path in outputs.items()},
    }
    manifest_path = output_dir / "figure_manifest.json"
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(manifest, indent=2) + "\n")

    print(json.dumps({"outputs": {name: str(path) for name, path in outputs.items()}, **manifest}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
