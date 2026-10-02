#!/usr/bin/env python3
"""Plot the SimpleMOS M30 double-off causal-closure audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/double_off_causal_closure"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m30"
          / "simplemos_m30_double_off_causal_closure.png")


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    decomposition = rows("m30_native_gap_decomposition.csv")
    interactions = rows("m30_nodal_interaction.csv")
    controls = rows("m30_solver_control_summary.csv")
    by_cell = {row["cell"]: row for row in decomposition}
    double_off = by_cell["bgn_off_srh_off"]

    ink = "#243247"
    blue = "#3978A8"
    blue_light = "#DCE8F2"
    orange = "#D88932"
    orange_light = "#F5E3CE"
    grey = "#7C8798"
    grid = "#E5E9EF"
    fig, axes = plt.subplots(2, 2, figsize=(13.4, 9.4))
    fig.patch.set_facecolor("white")
    for ax in axes.flat:
        ax.set_facecolor("white")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(colors=ink, labelsize=9)

    # Exact additive decomposition of the native current gap.
    state = float(double_off["common_operator_state_contribution_dex"])
    operator = float(double_off["sentaurus_state_vela_operator_contribution_dex"])
    total = float(double_off["native_gap_dex"])
    x = np.arange(3)
    bottoms = [0.0, state, 0.0]
    heights = [state, operator, total]
    colors = [blue, orange, blue_light]
    axes[0, 0].bar(x, heights, bottom=bottoms, width=0.58, color=colors,
                   edgecolor=ink, linewidth=0.7)
    axes[0, 0].plot([0, 1], [state, state], color=grey, linewidth=0.9,
                    linestyle="--")
    for i, (bottom, height) in enumerate(zip(bottoms, heights)):
        axes[0, 0].text(i, bottom + height + 0.035, f"{height:.3f}",
                        ha="center", va="bottom", color=ink, fontsize=9)
    axes[0, 0].set_xticks(x, ["State branch\n(common operator)",
                              "Vela operator /\ncontact extraction", "Native gap"])
    axes[0, 0].set_ylabel("log10 current ratio (dex)", color=ink)
    axes[0, 0].set_ylim(0, 1.48)
    axes[0, 0].set_title("Double-off current-gap decomposition", loc="left", color=ink)
    axes[0, 0].grid(axis="y", color=grid, linewidth=0.8)

    # The same decomposition across all four factorial cells.
    order = ["bgn_on_srh_on", "bgn_on_srh_off",
             "bgn_off_srh_on", "bgn_off_srh_off"]
    labels = ["BGN on / SRH on", "BGN on / SRH off",
              "BGN off / SRH on", "BGN off / SRH off"]
    state_values = np.array([
        float(by_cell[cell]["common_operator_state_contribution_dex"])
        for cell in order])
    operator_values = np.array([
        float(by_cell[cell]["sentaurus_state_vela_operator_contribution_dex"])
        for cell in order])
    y = np.arange(len(order))
    height = 0.34
    axes[0, 1].barh(y - height / 2, state_values, height, color=blue,
                    edgecolor=ink, linewidth=0.55, label="State contribution")
    axes[0, 1].barh(y + height / 2, operator_values, height, color=orange_light,
                    edgecolor=orange, linewidth=1.0, label="Operator/extraction")
    axes[0, 1].axvline(0.0, color=ink, linewidth=0.9)
    axes[0, 1].set_yticks(y, labels)
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlabel("Signed contribution (dex)", color=ink)
    axes[0, 1].set_title("Four-cell frozen-state decomposition", loc="left", color=ink)
    axes[0, 1].grid(axis="x", color=grid, linewidth=0.8)
    axes[0, 1].legend(frameon=False, loc="lower right")

    # BGN x SRH interaction in carrier quasi-Fermi potentials.
    for field, color, face, marker, label in (
            ("phin", blue, blue_light, "o", "Electron QF"),
            ("phip", orange, orange_light, "^", "Hole QF")):
        selected = [row for row in interactions if row["field"] == field]
        sx = np.array([float(row["sentaurus_interaction"]) for row in selected])
        sy = np.array([float(row["vela_interaction"]) for row in selected])
        axes[1, 0].scatter(sx, sy, s=14, facecolor=face, edgecolor=color,
                           linewidth=0.45, alpha=0.62, marker=marker, label=label)
    all_values = [float(row[key]) for row in interactions
                  if row["field"] in ("phin", "phip")
                  for key in ("sentaurus_interaction", "vela_interaction")]
    low, high = min(all_values), max(all_values)
    pad = max((high - low) * 0.06, 1e-12)
    axes[1, 0].plot([low - pad, high + pad], [low - pad, high + pad],
                    color=ink, linestyle="--", linewidth=1.0, label="1:1")
    axes[1, 0].set_xlim(low - pad, high + pad)
    axes[1, 0].set_ylim(low - pad, high + pad)
    axes[1, 0].set_xlabel("Sentaurus BGN x SRH interaction (V)", color=ink)
    axes[1, 0].set_ylabel("Vela BGN x SRH interaction (V)", color=ink)
    axes[1, 0].set_title("Nodal quasi-Fermi interaction", loc="left", color=ink)
    axes[1, 0].grid(True, color=grid, linewidth=0.8)
    axes[1, 0].legend(frameon=False, loc="upper left")

    # Failed reclosure and strict path controls expose the numerical floor.
    control_labels = [
        "Baseline\nreclosure", "Tight\nreclosure", "Strict\nreclosure",
        "Row-scaled\nreclosure", "Tight path\nrepeat", "Strict path\nrepeat"]
    residuals = np.array([float(row["failure_residual_norm"]) for row in controls])
    cx = np.arange(len(controls))
    axes[1, 1].bar(cx, residuals, width=0.66,
                   color=[blue_light] * 4 + [orange_light] * 2,
                   edgecolor=[blue] * 4 + [orange] * 2, linewidth=1.0)
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xticks(cx, control_labels)
    axes[1, 1].set_ylabel("Failure residual norm", color=ink)
    axes[1, 1].set_title("Solver-control failure residuals", loc="left", color=ink)
    axes[1, 1].grid(axis="y", which="both", color=grid, linewidth=0.8)
    axes[1, 1].text(0.98, 0.95, "6 / 6 controls: non-converged",
                    transform=axes[1, 1].transAxes, ha="right", va="top",
                    color=grey, fontsize=9)

    fig.suptitle("SimpleMOS M30 — double-off causal closure audit",
                 x=0.065, y=0.985, ha="left", color=ink,
                 fontsize=15, fontweight="bold")
    fig.text(0.065, 0.95,
             "n23, Vd = 0.05 V, Vg = 0.05 V; HFS off; 32 state/operator pairs, 96 probes, 942 silicon nodes",
             ha="left", color=grey, fontsize=10)
    fig.tight_layout(rect=(0.035, 0.035, 0.99, 0.925), h_pad=2.4, w_pad=2.2)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=190, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
