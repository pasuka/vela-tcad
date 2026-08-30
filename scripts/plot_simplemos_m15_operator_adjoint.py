#!/usr/bin/env python3
"""Render static M15 operator-adjoint evidence figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/operator_adjoint"
FIGURES = REPO / "docs/validation/figures/simplemos_m15"
BLUE = "#3568A8"
ORANGE = "#D9822B"
GOLD = "#B08D18"
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


def component_response() -> None:
    data = {row["component"]: row for row in rows("m15_component_contributions.csv")
            if row["state"] == KEY_STATE}
    components = ["electron_transport", "poisson", "electron_srh",
                  "hole_transport", "hole_srh"]
    labels = ["Electron transport", "Poisson", "Electron SRH",
              "Hole transport", "Hole SRH"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in components]
    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    y = np.arange(len(labels))
    bars = ax.barh(y, values, color=[BLUE, GOLD, ORANGE, GREY, GREY],
                   edgecolor=INK, linewidth=0.5)
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_xlim(min(0.0, min(values) * 1.12), max(values) * 1.14)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("Equation-term response at n21, Vd=0.05 V, Vg=0.8 V", pad=28)
    ax.text(0.0, 1.01, "Vela residual on mapped Sentaurus state; signed contribution",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, values):
        if abs(value) >= 0.001:
            ax.text(value + 0.15, bar.get_y() + bar.get_height() / 2,
                    f"{value:+.3f}%", va="center")
    ax.text(0.99, 0.04,
            f"Electron SRH = {values[2]:.2e}%\nHole transport = {values[3]:.2e}%",
            transform=ax.transAxes, ha="right", va="bottom", color="#5B6573")
    finish(fig, "m15_component_response_n21_vg0p8.png")


def region_response() -> None:
    data = {row["region"]: row for row in rows("m15_region_contributions.csv")
            if row["state"] == KEY_STATE}
    names = ["source", "channel", "drain", "body"]
    labels = ["Source", "Channel", "Drain", "Body"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in names]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    x = np.arange(len(names))
    bars = ax.bar(x, values, color=BLUE, edgecolor=INK, linewidth=0.5)
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_ylim(min(values) - 0.9, max(values) + 1.6)
    ax.set_xticks(x, labels)
    ax.set_ylabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("Spatial operator response at n21, Vd=0.05 V, Vg=0.8 V", pad=28)
    ax.text(0.0, 1.01, "Mutually exclusive M14 partition; signed equation-term sum",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, values):
        if abs(value) >= 0.001:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    value + (0.25 if value >= 0 else -0.25),
                    f"{value:+.3f}%", ha="center",
                    va="bottom" if value >= 0 else "top")
    ax.text(0.98, 0.92, f"Source = {values[0]:.2e}%\nBody = {values[3]:.2e}%",
            transform=ax.transAxes, ha="right", va="top", color="#5B6573")
    finish(fig, "m15_region_response_n21_vg0p8.png")


def component_fraction_matrix() -> None:
    data = rows("m15_component_contributions.csv")
    states = list(dict.fromkeys(row["state"] for row in data))
    indexed = {(row["state"], row["component"]): row for row in data}
    electron = [float(indexed[(state, "electron_transport")]
                      ["absolute_contribution_fraction"]) for state in states]
    poisson = [float(indexed[(state, "poisson")]
                     ["absolute_contribution_fraction"]) for state in states]
    other = [max(0.0, 1.0 - e - p) for e, p in zip(electron, poisson)]
    x = np.arange(len(states))
    fig, ax = plt.subplots(figsize=(10.0, 4.6))
    ax.bar(x, electron, color=BLUE, edgecolor=INK, linewidth=0.3,
           label="Electron transport")
    ax.bar(x, poisson, bottom=electron, color=GOLD, edgecolor=INK,
           linewidth=0.3, label="Poisson")
    ax.bar(x, other, bottom=np.array(electron) + np.array(poisson),
           color=GREY, edgecolor=INK, linewidth=0.3, label="Other")
    labels = [state.replace("_vd_", "\nVd=").replace("_vg_", ", Vg=")
              for state in states]
    ax.set_xticks(x, labels, rotation=65, ha="right", fontsize=7)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("fraction of absolute adjoint contribution")
    ax.set_title("Equation-family contribution fractions across 16 states", pad=72)
    ax.text(0.0, 1.095, "Fractions use absolute component magnitudes; signs are shown in detailed CSV",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(frameon=False, ncol=3, loc="lower center",
              bbox_to_anchor=(0.5, 1.02))
    finish(fig, "m15_component_fraction_16_states.png")


def m11_crosscheck() -> None:
    data = {row["term"]: row for row in rows("m15_m11_mobility_crosscheck.csv")
            if row["state"] == KEY_STATE}
    names = ["phumob", "enormal", "hfs"]
    labels = ["PhuMob", "Enormal", "HFS"]
    values = [abs(float(data[name]["effect"])) for name in names]
    fig, ax = plt.subplots(figsize=(6.8, 4.3))
    bars = ax.bar(labels, values, color=[BLUE, GOLD, ORANGE],
                  edgecolor=INK, linewidth=0.5)
    ax.set_yscale("log")
    ax.set_ylabel("absolute factorial main effect on log10(Id) (dex)")
    ax.set_title("M11 frozen mobility-factor cross-check", pad=28)
    ax.text(0.0, 1.01, "n21, Vd=0.05 V, Vg=0.8 V; same mapped Sentaurus state",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7, which="both")
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value * 1.35,
                f"{value:.3g}", ha="center", va="bottom")
    finish(fig, "m15_m11_mobility_crosscheck_n21_vg0p8.png")


def spatial_map() -> None:
    data = [row for row in rows("m15_key_state_node_contributions.csv")
            if row["state"] == KEY_STATE]
    x = np.array([float(row["y_um"]) for row in data])
    y = np.array([float(row["x_um"]) for row in data])
    value = np.array([float(row["total_A_per_um"]) for row in data])
    maximum = max(abs(value.min()), abs(value.max()))
    cmap = LinearSegmentedColormap.from_list(
        "blue_white_orange", [BLUE, "#F7F8FA", ORANGE])
    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    scatter = ax.scatter(x, y, c=value, s=9, cmap=cmap,
                         norm=SymLogNorm(linthresh=max(maximum * 1.0e-5, 1.0e-30),
                                        vmin=-maximum, vmax=maximum),
                         linewidths=0.0)
    ax.invert_yaxis()
    ax.set_xlabel("lateral coordinate y (um)")
    ax.set_ylabel("depth coordinate x (um)")
    ax.set_title("Nodal adjoint-weighted operator response", pad=28)
    ax.text(0.0, 1.01, "n21, Vd=0.05 V, Vg=0.8 V; symmetric-log color scale",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("nodal contribution (A/um)")
    finish(fig, "m15_nodal_response_n21_vg0p8.png")


def main() -> None:
    style()
    component_response()
    region_response()
    component_fraction_matrix()
    m11_crosscheck()
    spatial_map()


if __name__ == "__main__":
    main()
