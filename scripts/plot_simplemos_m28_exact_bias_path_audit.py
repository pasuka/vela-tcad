#!/usr/bin/env python3
"""Plot the SimpleMOS M28 exact-bias continuation-path audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/exact_bias_path_audit"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m28"
          / "simplemos_m28_exact_bias_path_audit.png")


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    terminal = rows("m28_terminal_route_ledger.csv")
    state = rows("m28_state_route_ledger.csv")
    labels = [
        "M8\nCurrentPlot", "M27 tiny\nendpoint", "Tiny +\nDoZero",
        "Original step\n+ DoZero", "Direct\n+ DoZero",
    ]
    x = np.arange(len(terminal))
    full = np.array([float(row["full_current_A_per_um"]) for row in terminal])
    no_hfs = np.array([float(row["no_hfs_current_A_per_um"]) for row in terminal])
    response = np.array([abs(float(row["hfs_response_A_per_um"])) for row in terminal])
    common = np.array([float(row["common_mode_offset_A_per_um"]) for row in terminal])
    differential = np.array([float(row["differential_offset_A_per_um"]) for row in terminal])

    ink = "#243247"
    blue = "#3978A8"
    orange = "#D88932"
    grey = "#7C8798"
    light = "#DCE8F2"
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.1))
    fig.patch.set_facecolor("white")
    for ax in axes.flat:
        ax.set_facecolor("white")
        ax.grid(axis="y", color="#E5E9EF", linewidth=0.8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(colors=ink, labelsize=9)

    width = 0.36
    axes[0, 0].bar(x - width / 2, full / 1e-16, width, color=blue,
                   edgecolor=ink, linewidth=0.6, label="HFS on")
    axes[0, 0].bar(x + width / 2, no_hfs / 1e-16, width, color=light,
                   edgecolor=blue, linewidth=1.1, label="HFS off")
    axes[0, 0].set_xticks(x, labels)
    axes[0, 0].set_ylabel("Drain current (10⁻¹⁶ A/µm)", color=ink)
    axes[0, 0].set_title("Absolute deep-off current by continuation route", loc="left", color=ink)
    axes[0, 0].legend(frameon=False, ncol=2, loc="upper left")

    axes[0, 1].bar(x, response / 1e-17, color=orange, edgecolor=ink, linewidth=0.6)
    axes[0, 1].axhline(response[0] / 1e-17, color=ink, linestyle="--", linewidth=1.1,
                       label="M8 response")
    axes[0, 1].set_xticks(x, labels)
    axes[0, 1].set_ylabel("|HFS on − off| (10⁻¹⁷ A/µm)", color=ink)
    axes[0, 1].set_title("HFS difference remains stable", loc="left", color=ink)
    axes[0, 1].legend(frameon=False, loc="lower right")

    route_x = x[1:]
    axes[1, 0].bar(route_x, common[1:] / 1e-16, width * 1.35,
                   color=blue, edgecolor=ink, linewidth=0.6, label="Common mode")
    axes[1, 0].axhline(0.0, color=ink, linewidth=0.8)
    axes[1, 0].set_xticks(route_x, labels[1:])
    axes[1, 0].set_ylabel("Common-mode offset (10⁻¹⁶ A/µm)", color=blue)
    axes[1, 0].set_title("Route error is overwhelmingly common-mode", loc="left", color=ink)
    delta_axis = axes[1, 0].twinx()
    delta_axis.plot(route_x, differential[1:] / 1e-19, color=orange,
                    marker="D", markerfacecolor="white", markeredgecolor=orange,
                    linewidth=1.4, label="Differential")
    delta_axis.set_ylabel("Differential offset (10⁻¹⁹ A/µm)", color=orange)
    delta_axis.tick_params(axis="y", colors=orange, labelsize=9)
    delta_axis.spines["top"].set_visible(False)
    handles1, labels1 = axes[1, 0].get_legend_handles_labels()
    handles2, labels2 = delta_axis.get_legend_handles_labels()
    axes[1, 0].legend(handles1 + handles2, labels1 + labels2,
                      frameon=False, loc="upper right")

    chosen = [row for row in state if row["route"] == "original_initial_dozero"]
    sx = np.array([float(row["m27_hfs_delta_phin_V"]) * 1e3 for row in chosen])
    sy = np.array([float(row["hfs_delta_phin_V"]) * 1e3 for row in chosen])
    axes[1, 1].scatter(sx, sy, s=13, facecolor=light, edgecolor=blue,
                       linewidth=0.45, alpha=0.8)
    limits = [min(sx.min(), sy.min()), max(sx.max(), sy.max())]
    axes[1, 1].plot(limits, limits, color=ink, linestyle="--", linewidth=1.0)
    axes[1, 1].set_xlabel("M27 tiny-route Δφn (mV)", color=ink)
    axes[1, 1].set_ylabel("Original-step Δφn (mV)", color=ink)
    axes[1, 1].set_title("HFS state response is route invariant (942 nodes)", loc="left", color=ink)
    axes[1, 1].grid(True, color="#E5E9EF", linewidth=0.8)

    fig.suptitle("SimpleMOS M28 — exact-bias continuation-path audit",
                 x=0.07, y=0.985, ha="left", color=ink, fontsize=15, fontweight="bold")
    fig.text(0.07, 0.95,
             "n23, Vd = 0.05 V, Vg = 0.05 V; same mesh and physics pair, only gate continuation path changes",
             ha="left", color=grey, fontsize=10)
    fig.tight_layout(rect=(0.04, 0.035, 0.99, 0.925), h_pad=2.3, w_pad=2.0)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=190, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
