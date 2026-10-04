#!/usr/bin/env python3
"""Render SimpleMOS M4 mesh, curve, and nodal-field comparison figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_NEUTRAL = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m2_gate"
)
DEFAULT_COMPARISONS = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "controlled_mobility" / "comparisons"
)
DEFAULT_SENTAURUS_EXPORT = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m4_visuals" / "a0_export"
)
DEFAULT_VELA_STATE = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m4_controlled_mobility_final" / "vela" / "a0" / "vd_0p05"
    / "20_gate_sweep_accepted_state.csv"
)
DEFAULT_OUTPUT = (
    REPO / "docs" / "validation" / "figures" / "simplemos_m4"
)

VARIANT_TITLES = {
    "a0": "A0  Constant mobility",
    "a1": "A1  + Doping dependence",
    "a2": "A2  + High-field saturation",
    "a3": "A3  + Enormal degradation",
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def load_geometry(neutral: Path) -> tuple[
        dict[int, tuple[float, float]], list[dict[str, Any]],
        list[dict[str, Any]], dict[int, float]]:
    nodes = {
        int(row["id"]): (float(row["y_um"]), float(row["x_um"]))
        for row in rows(neutral / "nodes.csv")
    }
    elements = [{
        "id": int(row["id"]),
        "node_ids": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
        "region": row["region"],
        "material": row["material"],
    } for row in rows(neutral / "elements.csv")]
    contacts = [{
        "name": row["name"],
        "node_ids": [int(value) for value in row["node_ids"].split(";") if value],
    } for row in rows(neutral / "contacts.csv")]
    doping = {
        int(row["node_id"]): float(row["donors_cm3"])
        - float(row["acceptors_cm3"])
        for row in rows(neutral / "doping.csv")
    }
    return nodes, elements, contacts, doping


def polygons(nodes: dict[int, tuple[float, float]],
             elements: list[dict[str, Any]], *, silicon_only: bool = False
             ) -> tuple[list[list[tuple[float, float]]], list[dict[str, Any]]]:
    selected = [item for item in elements
                if not silicon_only or item["material"] == "Si"]
    return ([[nodes[node] for node in item["node_ids"]] for item in selected],
            selected)


def style_axis(ax: Any, title: str) -> None:
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_xlabel("Lateral coordinate (um)")
    ax.set_ylabel("Depth coordinate (um)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)


def plot_mesh(neutral: Path, output: Path) -> Path:
    nodes, elements, contacts, _ = load_geometry(neutral)
    material_colors = {"Si": "#6baed6", "SiO2": "#f6c85f", "Nitride": "#9b8ac4"}
    grouped: dict[str, list[list[tuple[float, float]]]] = {}
    segments: list[list[tuple[float, float]]] = []
    for item in elements:
        poly = [nodes[node] for node in item["node_ids"]]
        grouped.setdefault(item["material"], []).append(poly)
        segments.extend(([poly[0], poly[1]], [poly[1], poly[2]], [poly[2], poly[0]]))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
    views = [((-0.52, 0.52), (1.02, -0.11), "Full device mesh"),
             ((-0.28, 0.28), (0.24, -0.11), "Gate and channel detail")]
    contact_colors = {"source": "#2166ac", "drain": "#b2182b",
                      "gate": "#542788", "substrate": "#333333"}
    for ax, (xlim, ylim, title) in zip(axes, views):
        for material, polys in grouped.items():
            ax.add_collection(PolyCollection(
                polys, facecolor=material_colors.get(material, "#cccccc"),
                edgecolor="none", alpha=0.43))
        ax.add_collection(LineCollection(
            segments, colors="#263746", linewidths=0.22, alpha=0.48))
        for contact in contacts:
            points = np.asarray([nodes[node] for node in contact["node_ids"]])
            ax.scatter(points[:, 0], points[:, 1], s=10,
                       color=contact_colors[contact["name"]], zorder=5)
            center = points.mean(axis=0)
            ax.annotate(contact["name"].capitalize(), center,
                        xytext=(0, -11 if contact["name"] == "substrate" else 8),
                        textcoords="offset points", ha="center", fontsize=8,
                        weight="bold", color=contact_colors[contact["name"]])
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        style_axis(ax, title)
    handles = [Patch(facecolor=color, alpha=.55, label=label)
               for label, color in (("Silicon", material_colors["Si"]),
                                    ("Oxide", material_colors["SiO2"]),
                                    ("Nitride spacer", material_colors["Nitride"]))]
    handles.extend(Line2D([], [], marker="o", linestyle="none", markersize=6,
                          color=color, label=name.capitalize())
                   for name, color in contact_colors.items())
    fig.legend(handles=handles, loc="upper center", ncol=7, frameon=False,
               bbox_to_anchor=(0.5, .925), fontsize=9)
    fig.suptitle("SimpleMOS nominal 2-D triangular mesh", fontsize=15,
                 weight="bold", y=.985)
    fig.subplots_adjust(left=.065, right=.985, bottom=.09, top=.82, wspace=.2)
    path = output / "simplemos-m4-device-mesh.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def plot_doping(neutral: Path, output: Path) -> Path:
    nodes, elements, contacts, doping = load_geometry(neutral)
    polys, silicon = polygons(nodes, elements, silicon_only=True)
    values = np.asarray([
        np.mean([doping[node] for node in item["node_ids"]]) for item in silicon])
    norm = SymLogNorm(linthresh=1e14, linscale=.8, vmin=-1.5e20, vmax=1.5e20)
    fig, ax = plt.subplots(figsize=(12.5, 6.6))
    collection = PolyCollection(polys, array=values, cmap="RdBu_r", norm=norm,
                                edgecolor="#2c3e50", linewidth=.16)
    ax.add_collection(collection)
    for contact in contacts:
        points = np.asarray([nodes[node] for node in contact["node_ids"]])
        ax.scatter(points[:, 0], points[:, 1], s=8, color="#202020", zorder=4)
    ax.set_xlim(-.52, .52)
    ax.set_ylim(1.02, -.11)
    style_axis(ax, "Signed net doping on the imported Silicon mesh")
    colorbar = fig.colorbar(collection, ax=ax, pad=.02, shrink=.9)
    colorbar.set_label("Net doping Nd - Na (cm$^{-3}$)")
    fig.suptitle("SimpleMOS nominal doping distribution", fontsize=15,
                 weight="bold", y=.98)
    fig.subplots_adjust(left=.08, right=.9, bottom=.1, top=.88)
    path = output / "simplemos-m4-net-doping.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def load_comparison(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = rows(path)
    return (
        np.asarray([float(row["gate_voltage_V"]) for row in data]),
        np.abs(np.asarray([float(row["sentaurus_current_A_per_um"]) for row in data])),
        np.abs(np.asarray([float(row["vela_current_A_per_um"]) for row in data])),
    )


def plot_curves(comparisons: Path, output: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True)
    colors = {"0p05": "#2166ac", "1": "#d6604d"}
    for ax, variant in zip(axes.flat, VARIANT_TITLES):
        for drain_tag, label in (("0p05", "Vd = 0.05 V"), ("1", "Vd = 1.0 V")):
            gate, sentaurus, vela = load_comparison(
                comparisons / f"{variant}_vd_{drain_tag}_comparison.csv")
            color = colors[drain_tag]
            ax.semilogy(gate, sentaurus, color=color, lw=2.0,
                        label=f"{label}, Sentaurus")
            ax.semilogy(gate, vela, color=color, lw=1.7, ls="--",
                        label=f"{label}, Vela")
        ax.set_title(VARIANT_TITLES[variant], fontsize=11, weight="bold")
        ax.grid(True, which="both", alpha=.22)
        ax.set_xlim(0, 2.5)
        ax.set_ylim(1e-9, 4e3)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage Vg (V)")
    for ax in axes[:, 0]:
        ax.set_ylabel("|Drain current| (A/um)")
    handles = [
        Line2D([], [], color=colors["0p05"], lw=2, label="Vd=0.05 V, Sentaurus"),
        Line2D([], [], color=colors["0p05"], lw=2, ls="--", label="Vd=0.05 V, Vela"),
        Line2D([], [], color=colors["1"], lw=2, label="Vd=1.0 V, Sentaurus"),
        Line2D([], [], color=colors["1"], lw=2, ls="--", label="Vd=1.0 V, Vela"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(.5, .945))
    fig.suptitle("SimpleMOS controlled-mobility Id-Vg comparison",
                 fontsize=15, weight="bold", y=.99)
    fig.subplots_adjust(left=.08, right=.985, bottom=.075, top=.89,
                        hspace=.2, wspace=.12)
    path = output / "simplemos-m4-idvg-comparison.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def nodal_values(path: Path, column: str) -> dict[int, float]:
    return {int(row["node_id"]): float(row[column]) for row in rows(path)}


def cell_average(items: list[dict[str, Any]], values: dict[int, float]) -> np.ndarray:
    return np.asarray([
        np.mean([values[node] for node in item["node_ids"]]) for item in items])


def add_field(ax: Any, polys: list[list[tuple[float, float]]], values: np.ndarray,
              title: str, cmap: str, norm: Any) -> PolyCollection:
    collection = PolyCollection(polys, array=values, cmap=cmap, norm=norm,
                                edgecolor="none")
    ax.add_collection(collection)
    ax.set_xlim(-.52, .52)
    ax.set_ylim(1.02, -.03)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title, fontsize=10, weight="bold")
    ax.set_xlabel("Lateral coordinate (um)")
    ax.set_ylabel("Depth coordinate (um)")
    return collection


def plot_fields(neutral: Path, sentaurus_export: Path, vela_state: Path,
                output: Path) -> Path:
    nodes, elements, _, _ = load_geometry(neutral)
    polys, silicon = polygons(nodes, elements, silicon_only=True)
    sentaurus_psi = nodal_values(
        sentaurus_export / "fields" / "ElectrostaticPotential_region0.csv",
        "component0")
    sentaurus_n = nodal_values(
        sentaurus_export / "fields" / "eDensity_region0.csv", "component0")
    vela_psi = nodal_values(vela_state, "psi")
    vela_n_m3 = nodal_values(vela_state, "electrons_m3")
    common_nodes = set(sentaurus_psi) & set(vela_psi)
    if len(common_nodes) != len(sentaurus_psi):
        raise ValueError("Sentaurus and Vela Silicon node sets do not match")
    vela_n = {node: value / 1e6 for node, value in vela_n_m3.items()}

    s_psi = cell_average(silicon, sentaurus_psi)
    v_psi = cell_average(silicon, vela_psi)
    d_psi = np.abs(v_psi - s_psi)
    s_n = np.maximum(cell_average(silicon, sentaurus_n), 1.0)
    v_n = np.maximum(cell_average(silicon, vela_n), 1.0)
    d_log_n = np.abs(np.log10(v_n / s_n))

    psi_norm = Normalize(vmin=min(s_psi.min(), v_psi.min()),
                         vmax=max(s_psi.max(), v_psi.max()))
    dpsi_norm = Normalize(vmin=0.0, vmax=max(float(d_psi.max()), 1e-12))
    density_norm = LogNorm(vmin=min(s_n.min(), v_n.min()),
                           vmax=max(s_n.max(), v_n.max()))
    density_error_norm = Normalize(vmin=0.0,
                                   vmax=max(float(d_log_n.max()), 1e-12))
    fig, axes = plt.subplots(
        2, 3, figsize=(15.5, 10.2), layout="constrained",
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.08]})
    top0 = add_field(axes[0, 0], polys, s_psi, "Sentaurus potential", "viridis", psi_norm)
    add_field(axes[0, 1], polys, v_psi, "Vela potential", "viridis", psi_norm)
    top2 = add_field(axes[0, 2], polys, d_psi,
                     f"|Difference|, max {d_psi.max():.2e} V", "magma", dpsi_norm)
    bottom0 = add_field(axes[1, 0], polys, s_n, "Sentaurus electron density",
                        "plasma", density_norm)
    add_field(axes[1, 1], polys, v_n, "Vela electron density", "plasma", density_norm)
    bottom2 = add_field(
        axes[1, 2], polys, d_log_n,
        f"|log10 density ratio|, max {d_log_n.max():.3f} dex",
        "magma", density_error_norm)
    fig.colorbar(top0, ax=axes[0, :2], shrink=.82, pad=.015,
                 label="Electrostatic potential (V)")
    fig.colorbar(top2, ax=axes[0, 2], shrink=.82, pad=.015,
                 label="Absolute difference (V)")
    fig.colorbar(bottom0, ax=axes[1, :2], shrink=.82, pad=.015,
                 label="Electron density (cm$^{-3}$)")
    fig.colorbar(bottom2, ax=axes[1, 2], shrink=.82, pad=.015,
                 label="Absolute log10 ratio (dex)")
    fig.suptitle("SimpleMOS A0 nodal fields at Vd=0.05 V, Vg=2.5 V",
                 fontsize=15, weight="bold", y=.985)
    path = output / "simplemos-m4-physical-fields-comparison.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--neutral-dir", type=Path, default=DEFAULT_NEUTRAL)
    parser.add_argument("--comparison-dir", type=Path, default=DEFAULT_COMPARISONS)
    parser.add_argument("--sentaurus-export-dir", type=Path,
                        default=DEFAULT_SENTAURUS_EXPORT)
    parser.add_argument("--vela-state", type=Path, default=DEFAULT_VELA_STATE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    generated = [
        plot_mesh(args.neutral_dir.resolve(), output),
        plot_doping(args.neutral_dir.resolve(), output),
        plot_curves(args.comparison_dir.resolve(), output),
        plot_fields(args.neutral_dir.resolve(), args.sentaurus_export_dir.resolve(),
                    args.vela_state.resolve(), output),
    ]
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
