#!/usr/bin/env python3
"""Plot the SimpleMOS M31 minority-hole and Poisson perturbation audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/minority_poisson_perturbation"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m31"
          / "simplemos_m31_minority_poisson_perturbation.png")


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    responses = rows("m31_state_variant_response.csv")
    nodes = rows("m31_poisson_node_audit.csv")
    cross = rows("m31_cross_block_summary.csv")
    by_variant = {row["variant"]: row for row in responses}

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

    # Frozen direct terminal-current response to state substitutions.
    names = ["hole_top20_replace", "hole_full_replace",
             "electron_full_replace", "qf_full_replace",
             "sentaurus_full_replace"]
    labels = ["Hole\ntop 20", "Hole\nfull", "Electron\nfull",
              "Both QF\nfull", "All state\nfull"]
    shifts = np.array([
        float(by_variant[name]["absolute_log10_current_shift_dex"])
        for name in names])
    x = np.arange(len(names))
    axes[0, 0].bar(x, shifts, width=0.66,
                   color=[orange_light, orange, blue, blue, blue_light],
                   edgecolor=[orange, orange, ink, ink, blue], linewidth=0.9)
    axes[0, 0].set_xticks(x, labels)
    axes[0, 0].set_ylabel("Absolute frozen current shift (dex)", color=ink)
    axes[0, 0].set_title("Frozen-state terminal-current response",
                         loc="left", color=ink)
    axes[0, 0].grid(axis="y", color=grid, linewidth=0.8)
    for index, value in enumerate(shifts):
        label = f"{value:.4f}" if value < 0.01 else f"{value:.3f}"
        axes[0, 0].text(index, value + 0.08, label, ha="center",
                        va="bottom", fontsize=8, color=ink)
    axes[0, 0].set_ylim(0.0, max(shifts) * 1.14)

    # Poisson residual concentration at the ten dominant nodes.
    residuals = np.array([float(row["abs_psi_residual"]) for row in nodes])
    node_labels = [f"{row['node_id']}" for row in nodes]
    y = np.arange(len(nodes))
    axes[0, 1].barh(y, residuals, color=blue_light,
                    edgecolor=blue, linewidth=0.9)
    axes[0, 1].set_yticks(y, node_labels)
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_xlabel("Absolute Poisson residual", color=ink)
    axes[0, 1].set_ylabel("Node ID", color=ink)
    axes[0, 1].set_title("Dominant Poisson-residual nodes",
                         loc="left", color=ink)
    axes[0, 1].grid(axis="x", color=grid, linewidth=0.8)
    axes[0, 1].text(0.98, 0.03, "Top 10 = 99.63% of free-Si L2 norm",
                    transform=axes[0, 1].transAxes, ha="right", va="bottom",
                    color=grey, fontsize=9)

    # Residual response to the Poisson-only Newton step.
    step_names = ["vela_baseline", "poisson_top10_step", "poisson_full_step"]
    step_labels = ["Baseline", "Top-10\nPoisson step", "Full\nPoisson step"]
    psi_norms = np.array([
        float(by_variant[name]["psi_residual_norm"]) for name in step_names])
    sx = np.arange(len(step_names))
    axes[1, 0].bar(sx, psi_norms, width=0.62,
                   color=[blue_light, blue, orange_light],
                   edgecolor=[blue, ink, orange], linewidth=0.9)
    axes[1, 0].set_xticks(sx, step_labels)
    axes[1, 0].set_ylabel("Poisson residual norm", color=ink)
    axes[1, 0].set_title("Frozen Poisson-step residual response",
                         loc="left", color=ink)
    axes[1, 0].grid(axis="y", color=grid, linewidth=0.8)
    for index, value in enumerate(psi_norms):
        axes[1, 0].text(index, value + max(psi_norms) * 0.025,
                        f"{value:.3e}", ha="center", va="bottom",
                        fontsize=8, color=ink)
    axes[1, 0].set_ylim(0.0, max(psi_norms) * 1.17)

    # Cross-block indirect coupling response for hole/electron substitutions.
    cross_labels = ["Hole full", "Electron full", "Both QF full"]
    qfp_product = np.array([
        float(row["psi_qfp_product_l2"]) for row in cross])
    raw_psi = np.array([
        float(row["full_raw_delta_psi_l2_V"]) for row in cross])
    cx = np.arange(len(cross))
    width = 0.34
    axes[1, 1].bar(cx - width / 2, qfp_product, width,
                   color=blue, edgecolor=ink, linewidth=0.6,
                   label="||Jpsi,qfp dqfp||2")
    axes[1, 1].bar(cx + width / 2, raw_psi, width,
                   color=orange_light, edgecolor=orange, linewidth=0.9,
                   label="||full raw dpsi||2")
    axes[1, 1].set_xticks(cx, cross_labels)
    axes[1, 1].set_ylabel("L2 response", color=ink)
    axes[1, 1].set_title("Indirect Poisson-QFP coupling",
                         loc="left", color=ink)
    axes[1, 1].grid(axis="y", color=grid, linewidth=0.8)
    axes[1, 1].legend(frameon=False, fontsize=8, loc="upper left")

    fig.suptitle("SimpleMOS M31 — minority-hole and Poisson perturbation audit",
                 x=0.065, y=0.985, ha="left", color=ink,
                 fontsize=15, fontweight="bold")
    fig.text(0.065, 0.95,
             "n23, Vd = 0.05 V, Vg = 0.05 V; BGN/SRH/HFS off; frozen-state and Jacobian diagnostics",
             ha="left", color=grey, fontsize=10)
    fig.tight_layout(rect=(0.035, 0.035, 0.99, 0.925), h_pad=2.4, w_pad=2.2)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=190, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
