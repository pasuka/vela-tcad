#!/usr/bin/env python3
"""Plot the SimpleMOS M32 scaled electron quasi-Fermi perturbation audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/electron_qf_scaled_perturbation"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m32"
          / "simplemos_m32_electron_qf_scaled_perturbation.png")


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    response = rows("m32_fraction_response.csv")
    cross = rows("m32_cross_block_scaling.csv")
    cut = [row for row in rows("m32_drain_cut_edge_response.csv")
           if float(row["fraction"]) == 1.0]

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

    nonzero = [row for row in response if float(row["fraction"]) > 0.0]
    labels = [f"{100 * float(row['fraction']):g}%" for row in nonzero]
    exact = np.array([float(row["exact_delta_current_A_per_um"])
                      for row in nonzero]) * 1.0e15
    tangent = np.array([float(row["tangent_predicted_delta_current_A_per_um"])
                        for row in nonzero]) * 1.0e15
    x = np.arange(len(nonzero))
    width = 0.36
    axes[0, 0].bar(x - width / 2, exact, width, color=blue,
                   edgecolor=ink, linewidth=0.6, label="Exact frozen SG")
    axes[0, 0].bar(x + width / 2, tangent, width, color=orange_light,
                   edgecolor=orange, linewidth=0.9, label="Baseline tangent")
    axes[0, 0].axhline(0.0, color=ink, linewidth=0.9)
    axes[0, 0].set_xticks(x, labels)
    axes[0, 0].set_ylabel("Drain-current change (fA/um)", color=ink)
    axes[0, 0].set_title("Frozen current: exact response vs tangent",
                         loc="left", color=ink)
    axes[0, 0].grid(axis="y", color=grid, linewidth=0.8)
    axes[0, 0].legend(frameon=False, fontsize=8, loc="lower left")

    errors = np.array([float(row["tangent_relative_error"])
                       for row in nonzero])
    conditions = np.array([float(row["drain_cut_condition"])
                            for row in nonzero])
    axes[0, 1].bar(x, errors, width=0.62, color=blue_light,
                   edgecolor=blue, linewidth=0.9)
    axes[0, 1].axhline(0.1, color=orange, linestyle="--", linewidth=1.1,
                       label="10% linearity criterion")
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_ylabel("Tangent relative error", color=ink)
    axes[0, 1].set_title("Linearization error and drain-cut conditioning",
                         loc="left", color=ink)
    axes[0, 1].grid(axis="y", color=grid, linewidth=0.8)
    axes[0, 1].legend(frameon=False, fontsize=8, loc="upper left")
    for index, (error, condition) in enumerate(zip(errors, conditions)):
        axes[0, 1].text(index, error + 0.035, f"k={condition:.2f}",
                        ha="center", va="bottom", fontsize=8, color=grey)
    axes[0, 1].set_ylim(0.0, max(errors) * 1.2)

    nonzero_cross = [row for row in cross if float(row["fraction"]) > 0.0]
    full_qfp = float(nonzero_cross[-1]["induced_psi_qfp_product_l2"])
    full_psi = float(nonzero_cross[-1]["induced_full_raw_delta_psi_l2_V"])
    expected = np.array([float(row["fraction"]) for row in nonzero_cross])
    qfp_ratio = np.array([
        float(row["induced_psi_qfp_product_l2"]) / full_qfp
        for row in nonzero_cross])
    psi_ratio = np.array([
        float(row["induced_full_raw_delta_psi_l2_V"]) / full_psi
        for row in nonzero_cross])
    axes[1, 0].bar(x - width, expected, width, color="white",
                   edgecolor=ink, linewidth=1.0, hatch="//", label="Ideal fraction")
    axes[1, 0].bar(x, qfp_ratio, width, color=blue,
                   edgecolor=ink, linewidth=0.6, label="Induced Jpsi,qfp response")
    axes[1, 0].bar(x + width, psi_ratio, width, color=orange_light,
                   edgecolor=orange, linewidth=0.9, label="Induced raw dpsi")
    axes[1, 0].set_xticks(x, labels)
    axes[1, 0].set_ylabel("Response / full-fraction response", color=ink)
    axes[1, 0].set_title("Baseline-subtracted Poisson-QFP scaling",
                         loc="left", color=ink)
    axes[1, 0].grid(axis="y", color=grid, linewidth=0.8)
    axes[1, 0].legend(frameon=False, fontsize=8, loc="upper left")

    cut_sorted = sorted(
        cut,
        key=lambda row: abs(float(row["delta_electron_current_A_per_um"])),
        reverse=True)
    edge_labels = [f"edge {row['edge_id']}\n{row['node0']}-{row['node1']}"
                   for row in cut_sorted]
    edge_delta = np.array([
        float(row["delta_electron_current_A_per_um"]) * 1.0e15
        for row in cut_sorted])
    y = np.arange(len(cut_sorted))
    colors = [blue if value >= 0.0 else orange_light for value in edge_delta]
    outlines = [ink if value >= 0.0 else orange for value in edge_delta]
    axes[1, 1].barh(y, edge_delta, color=colors,
                    edgecolor=outlines, linewidth=0.9)
    axes[1, 1].axvline(0.0, color=ink, linewidth=0.9)
    axes[1, 1].set_yticks(y, edge_labels)
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_xlabel("Electron-current change at 100% (fA/um)", color=ink)
    axes[1, 1].set_title("Drain-contact cut edge contributions",
                         loc="left", color=ink)
    axes[1, 1].grid(axis="x", color=grid, linewidth=0.8)

    fig.suptitle("SimpleMOS M32 — scaled electron quasi-Fermi perturbation",
                 x=0.065, y=0.985, ha="left", color=ink,
                 fontsize=15, fontweight="bold")
    fig.text(0.065, 0.95,
             "n23, Vd = 0.05 V, Vg = 0.05 V; BGN/SRH/HFS off; 6 fractions, 18 probes, 7 drain-cut edges",
             ha="left", color=grey, fontsize=10)
    fig.tight_layout(rect=(0.035, 0.035, 0.99, 0.925), h_pad=2.4, w_pad=2.2)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=190, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
