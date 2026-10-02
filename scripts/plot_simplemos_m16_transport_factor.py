#!/usr/bin/env python3
"""Render static M16 transport-factor evidence figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/transport_factor"
FIGURES = REPO / "docs/validation/figures/simplemos_m16"
BLUE = "#3568A8"
GOLD = "#B08D18"
ORANGE = "#D9822B"
INK = "#252A31"
GREY = "#AAB2BE"
GRID = "#D8DDE5"
KEY_STATE = "n21_vd_0p05_vg_0p8"


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.edgecolor": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "text.color": INK,
        "axes.titleweight": "bold", "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def finish(fig: plt.Figure, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / name, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def grouped_factor_response() -> None:
    data = {row["factor"]: row for row in rows(
        "m16_grouped_factor_contributions.csv") if row["state"] == KEY_STATE}
    names = ["mobility_state", "mobility_drive", "sg_state_kernel"]
    labels = ["Mobility state", "Mobility drive", "SG state kernel"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in names]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    y = np.arange(len(names))
    colors = [GOLD if value < 0 else BLUE for value in values]
    bars = ax.barh(y, values, color=colors, edgecolor=INK, linewidth=0.5)
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_xlim(min(values) - 1.5, max(values) + 2.0)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("Grouped electron-transport factors at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015, "Three-factor Shapley allocation; signed response",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, values):
        label_x = value + 0.18 if value >= 0 else -0.12
        ax.text(label_x,
                bar.get_y() + bar.get_height() / 2,
                f"{value:+.3f}%", va="center",
                ha="left" if value >= 0 else "right")
    finish(fig, "m16_grouped_factor_response_n21_vg0p8.png")


def four_factor_cancellation() -> None:
    data = {row["factor"]: row for row in rows("m16_factor_contributions.csv")
            if row["state"] == KEY_STATE}
    names = ["mobility_state", "mobility_drive", "bernoulli_weights",
             "carrier_population"]
    labels = ["Mobility state", "Mobility drive", "Bernoulli weights",
              "Carrier population"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in names]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    y = np.arange(len(names))
    colors = [GOLD if value < 0 else BLUE for value in values]
    bars = ax.barh(y, values, color=colors, edgecolor=INK, linewidth=0.5)
    ax.set_xscale("symlog", linthresh=0.01, linscale=1.0)
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("adjoint-weighted response / baseline Id (%) — symmetric log scale")
    ax.set_title("Four-factor SG response at n21, Vd=0.05 V, Vg=0.8 V", pad=30)
    ax.text(0.0, 1.015,
            "Independent Bernoulli and population interventions expose their large cancellation",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, values):
        ax.text(0.97 if value >= 0 else 0.03,
                bar.get_y() + bar.get_height() / 2,
                f"{value:+.3g}%", transform=ax.get_yaxis_transform(),
                va="center", ha="right" if value >= 0 else "left")
    finish(fig, "m16_four_factor_cancellation_n21_vg0p8.png")


def grouped_fraction_matrix() -> None:
    data = rows("m16_grouped_factor_contributions.csv")
    states = list(dict.fromkeys(row["state"] for row in data))
    indexed = {(row["state"], row["factor"]): row for row in data}
    names = ["sg_state_kernel", "mobility_state", "mobility_drive"]
    labels = ["SG state kernel", "Mobility state", "Mobility drive"]
    colors = [BLUE, GOLD, ORANGE]
    values = {name: [float(indexed[(state, name)]["absolute_contribution_fraction"])
                     for state in states] for name in names}
    x = np.arange(len(states))
    fig, ax = plt.subplots(figsize=(10.2, 4.8))
    bottom = np.zeros(len(states))
    for name, label, color in zip(names, labels, colors):
        ax.bar(x, values[name], bottom=bottom, color=color, edgecolor=INK,
               linewidth=0.3, label=label)
        bottom += np.array(values[name])
    state_labels = [state.replace("_vd_", "\nVd=").replace("_vg_", ", Vg=")
                    for state in states]
    ax.set_xticks(x, state_labels, rotation=65, ha="right", fontsize=7)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("fraction of absolute grouped contribution")
    ax.set_title("Grouped transport-factor shares across 16 states", pad=62)
    ax.text(0.0, 1.08, "Absolute Shapley magnitudes; signs remain in detailed CSV",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.legend(frameon=False, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, 1.015))
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    finish(fig, "m16_grouped_factor_fraction_16_states.png")


def region_response() -> None:
    data = {row["region"]: row for row in rows(
        "m16_grouped_region_contributions.csv") if row["state"] == KEY_STATE}
    names = ["source", "channel", "drain", "body"]
    labels = ["Source", "Channel", "Drain", "Body"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in names]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    x = np.arange(len(names))
    bars = ax.bar(x, values, color=[BLUE if value >= 0 else GOLD for value in values],
                  edgecolor=INK, linewidth=0.5)
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_ylim(min(values) - 0.8, max(values) + 1.6)
    ax.set_xticks(x, labels)
    ax.set_ylabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("Transport-factor response by region at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015, "Grouped-factor sum; signed M14 node partition",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, values):
        if abs(value) >= 0.001:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    value + (0.22 if value >= 0 else -0.22),
                    f"{value:+.3f}%", ha="center",
                    va="bottom" if value >= 0 else "top")
    ax.text(0.98, 0.92, f"Source = {values[0]:.2e}%\nBody = {values[3]:.2e}%",
            transform=ax.transAxes, ha="right", va="top", color="#5B6573")
    finish(fig, "m16_region_response_n21_vg0p8.png")


def nodal_kernel_response() -> None:
    data = rows("m16_key_state_node_contributions.csv")
    x = np.array([float(row["x_um"]) for row in data])
    y = np.array([float(row["y_um"]) for row in data])
    z = np.array([float(row["grouped_sg_state_kernel_A_per_um"])
                  for row in data])
    nonzero = np.abs(z[np.nonzero(z)])
    linthresh = max(float(np.percentile(nonzero, 20)) if nonzero.size else 1.0e-20,
                    1.0e-20)
    vmax = max(float(np.max(np.abs(z))), linthresh * 10.0)
    cmap = LinearSegmentedColormap.from_list(
        "m16_signed", [GOLD, "#F4F1E4", "#FFFFFF", "#E8EFF8", BLUE])
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    scatter = ax.scatter(y, x, c=z, s=11, cmap=cmap,
                         norm=SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax),
                         linewidths=0.0)
    ax.set_xlabel("lateral coordinate y (um)")
    ax.set_ylabel("depth coordinate x (um)")
    ax.invert_yaxis()
    ax.set_title("Nodal SG state-kernel response at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015, "Adjoint-weighted grouped Shapley value; symmetric-log color",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("response (A/um)")
    finish(fig, "m16_nodal_sg_kernel_response_n21_vg0p8.png")


def main() -> None:
    style()
    grouped_factor_response()
    four_factor_cancellation()
    grouped_fraction_matrix()
    region_response()
    nodal_kernel_response()


if __name__ == "__main__":
    main()
