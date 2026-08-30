#!/usr/bin/env python3
"""Render static M17 SG secant-factor evidence figures."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/sg_secant"
FIGURES = REPO / "docs/validation/figures/simplemos_m17"
BLUE = "#3568A8"
GOLD = "#B08D18"
ORANGE = "#D9822B"
INK = "#252A31"
GRID = "#D8DDE5"
KEY_STATE = "n21_vd_0p05_vg_0p8"
FACTORS = ("mobility", "sg_secant_conductance", "qf_log_imbalance")
LABELS = ("Mobility", "SG secant conductance", "QF log imbalance")
COLORS = (ORANGE, GOLD, BLUE)


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


def key_factor_response() -> None:
    data = {row["factor"]: row for row in rows("m17_factor_contributions.csv")
            if row["state"] == KEY_STATE}
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in FACTORS]
    fig, ax = plt.subplots(figsize=(7.8, 4.5))
    y = np.arange(len(FACTORS))
    bars = ax.barh(y, values,
                   color=[BLUE if value >= 0 else GOLD for value in values],
                   edgecolor=INK, linewidth=0.5)
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_yticks(y, LABELS)
    ax.invert_yaxis()
    ax.set_xlabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("Stable SG secant factors at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015, "Three-factor Shapley allocation; signed response",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_xlim(min(values) - 2.0, max(values) + 2.5)
    for bar, value in zip(bars, values):
        ax.text(value + (0.25 if value >= 0 else -0.2),
                bar.get_y() + bar.get_height() / 2,
                f"{value:+.3f}%", va="center",
                ha="left" if value >= 0 else "right")
    finish(fig, "m17_factor_response_n21_vg0p8.png")


def factor_fraction_panels() -> None:
    data = rows("m17_factor_contributions.csv")
    indexed = {(row["device"], float(row["drain_voltage_V"]),
                float(row["gate_voltage_V"]), row["factor"]): row
               for row in data}
    gates = [0.0, 0.05, 0.8, 2.5]
    panels = [("n17", 0.05), ("n17", 1.0),
              ("n21", 0.05), ("n21", 1.0)]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.7), sharey=True)
    for ax, (device, drain) in zip(axes.flat, panels):
        bottom = np.zeros(len(gates))
        for factor, label, color in zip(FACTORS, LABELS, COLORS):
            values = [float(indexed[(device, drain, gate, factor)]
                            ["absolute_contribution_fraction"])
                      for gate in gates]
            ax.bar(np.arange(len(gates)), values, bottom=bottom, color=color,
                   edgecolor=INK, linewidth=0.35, label=label)
            bottom += np.array(values)
        ax.set_title(f"{device}, Vd={drain:g} V")
        ax.set_xticks(np.arange(len(gates)), [f"{gate:g}" for gate in gates])
        ax.set_xlabel("Vg (V)")
        ax.set_ylim(0.0, 1.0)
        ax.grid(axis="y", color=GRID, linewidth=0.7)
    axes[0, 0].set_ylabel("fraction of absolute contribution")
    axes[1, 0].set_ylabel("fraction of absolute contribution")
    fig.suptitle("SG secant-factor shares across 16 paired states",
                 y=0.985, fontweight="bold")
    fig.text(0.5, 0.935,
             "Discrete bias comparison; absolute Shapley magnitudes sum to one",
             ha="center", color="#5B6573")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.91))
    fig.subplots_adjust(top=0.82, hspace=0.38, wspace=0.16)
    finish(fig, "m17_factor_fraction_16_states.png")


def cancellation_reduction() -> None:
    report = json.loads((DATA / "m17_sg_secant_report.json")
                        .read_text(encoding="utf-8"))
    values = report["key_state"]["cancellation_amplification"]
    labels = ["M16 independent\nBernoulli/population",
              "M17 stable\nsecant factorization"]
    plotted = [values["m16_ungrouped_four_factor"],
               values["m17_stable_three_factor"]]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    bars = ax.barh(np.arange(2), plotted, color=[GOLD, BLUE],
                   edgecolor=INK, linewidth=0.5)
    ax.set_xscale("log")
    ax.set_yticks(np.arange(2), labels)
    ax.invert_yaxis()
    ax.set_xlabel("sum of absolute factor responses / absolute net response (x)")
    ax.set_title("Factor-cancellation amplification at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015,
            f"Stable reparameterization reduces amplification by {values['reduction_factor']:.1f}x",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="x", color=GRID, linewidth=0.7, which="both")
    for bar, value in zip(bars, plotted):
        ax.text(value * 1.08, bar.get_y() + bar.get_height() / 2,
                f"{value:.3g}x", va="center")
    finish(fig, "m17_cancellation_reduction_n21_vg0p8.png")


def region_response() -> None:
    data = {row["region"]: row for row in rows(
        "m17_region_factor_contributions.csv") if row["state"] == KEY_STATE}
    regions = ["source", "channel", "drain", "body"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in regions]
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    x = np.arange(len(regions))
    bars = ax.bar(x, values, color=[BLUE if value >= 0 else GOLD for value in values],
                  edgecolor=INK, linewidth=0.5)
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_xticks(x, [name.title() for name in regions])
    ax.set_ylabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("M17 factor response by region at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015, "All three stable factors; signed M14 node partition",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.set_ylim(min(values) - 0.8, max(values) + 1.6)
    for bar, value in zip(bars, values):
        if abs(value) >= 0.001:
            ax.text(bar.get_x() + bar.get_width() / 2,
                    value + (0.22 if value >= 0 else -0.22),
                    f"{value:+.3f}%", ha="center",
                    va="bottom" if value >= 0 else "top")
    ax.text(0.98, 0.92, f"Source = {values[0]:.2e}%\nBody = {values[3]:.2e}%",
            transform=ax.transAxes, ha="right", va="top", color="#5B6573")
    finish(fig, "m17_region_response_n21_vg0p8.png")


def nodal_qf_response() -> None:
    data = rows("m17_key_state_node_contributions.csv")
    x = np.array([float(row["x_um"]) for row in data])
    y = np.array([float(row["y_um"]) for row in data])
    z = np.array([float(row["qf_log_imbalance_A_per_um"]) for row in data])
    nonzero = np.abs(z[np.nonzero(z)])
    linthresh = max(float(np.percentile(nonzero, 20)) if nonzero.size else 1.0e-20,
                    1.0e-20)
    vmax = max(float(np.max(np.abs(z))), linthresh * 10.0)
    cmap = LinearSegmentedColormap.from_list(
        "m17_signed", [GOLD, "#F4F1E4", "#FFFFFF", "#E8EFF8", BLUE])
    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    scatter = ax.scatter(y, x, c=z, s=11, cmap=cmap,
                         norm=SymLogNorm(linthresh=linthresh,
                                         vmin=-vmax, vmax=vmax),
                         linewidths=0.0)
    ax.set_xlabel("lateral coordinate y (um)")
    ax.set_ylabel("depth coordinate x (um)")
    ax.invert_yaxis()
    ax.set_title("Nodal QF log-imbalance response at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.015, "Adjoint-weighted Shapley value; symmetric-log color",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("response (A/um)")
    finish(fig, "m17_nodal_qf_imbalance_n21_vg0p8.png")


def main() -> None:
    style()
    key_factor_response()
    factor_fraction_panels()
    cancellation_reduction()
    region_response()
    nodal_qf_response()


if __name__ == "__main__":
    main()
