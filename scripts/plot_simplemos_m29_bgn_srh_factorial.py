#!/usr/bin/env python3
"""Plot the SimpleMOS M29 BGN x SRH causal factorial audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/bgn_srh_factorial"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m29"
          / "simplemos_m29_bgn_srh_factorial.png")


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    matrix = rows("m29_factorial_matrix.csv")
    effects = rows("m29_factor_effects.csv")
    nodes = rows("m29_node_factor_effects.csv")
    labels = ["BGN on\nSRH on", "BGN on\nSRH off",
              "BGN off\nSRH on", "BGN off\nSRH off"]
    x = np.arange(len(matrix))
    sent = np.array([float(row["sentaurus_current_A_per_um"]) for row in matrix])
    vela = np.array([float(row["vela_current_A_per_um"]) for row in matrix])

    ink = "#243247"
    blue = "#3978A8"
    orange = "#D88932"
    grey = "#7C8798"
    light = "#DCE8F2"
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.2))
    fig.patch.set_facecolor("white")
    for ax in axes.flat:
        ax.set_facecolor("white")
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(colors=ink, labelsize=9)

    width = 0.36
    axes[0, 0].bar(x - width / 2, sent, width, color=blue,
                   edgecolor=ink, linewidth=0.6, label="Sentaurus")
    axes[0, 0].bar(x + width / 2, vela, width, color=light,
                   edgecolor=blue, linewidth=1.1, label="Vela")
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_ylim(1e-17, 7e-16)
    axes[0, 0].set_xticks(x, labels)
    axes[0, 0].set_ylabel("|Drain current| (A/um)", color=ink)
    axes[0, 0].set_title("Four-cell terminal-current matrix", loc="left", color=ink)
    axes[0, 0].grid(axis="y", which="both", color="#E5E9EF", linewidth=0.8)
    axes[0, 0].legend(frameon=False, ncol=2, loc="upper left")

    effect_names = ["bgn_main_effect_dex", "srh_main_effect_dex",
                    "bgn_srh_interaction_dex"]
    effect_labels = ["BGN main", "SRH main", "BGN x SRH\ninteraction"]
    chosen = {row["effect"]: row for row in effects}
    ex = np.arange(len(effect_names))
    sent_effect = np.array([float(chosen[name]["sentaurus_dex"])
                            for name in effect_names])
    vela_effect = np.array([float(chosen[name]["vela_dex"])
                            for name in effect_names])
    axes[0, 1].bar(ex - width / 2, sent_effect, width, color=blue,
                   edgecolor=ink, linewidth=0.6, label="Sentaurus")
    axes[0, 1].bar(ex + width / 2, vela_effect, width, color=light,
                   edgecolor=blue, linewidth=1.1, label="Vela")
    axes[0, 1].axhline(0.0, color=ink, linewidth=0.9)
    axes[0, 1].set_xticks(ex, effect_labels)
    axes[0, 1].set_ylabel("Effect on log10(|Id|) (dex)", color=ink)
    axes[0, 1].set_title("Factorial main effects and interaction", loc="left", color=ink)
    axes[0, 1].grid(axis="y", color="#E5E9EF", linewidth=0.8)
    axes[0, 1].legend(frameon=False, ncol=2, loc="upper left")

    for ax, factor, title in (
            (axes[1, 0], "bgn", "BGN nodal electron-density response"),
            (axes[1, 1], "srh", "SRH nodal electron-density response")):
        selected = [row for row in nodes
                    if row["factor"] == factor and row["field"] == "logn"]
        sx = np.array([float(row["sentaurus_effect"]) for row in selected])
        sy = np.array([float(row["vela_effect"]) for row in selected])
        ax.scatter(sx, sy, s=13, facecolor=light, edgecolor=blue,
                   linewidth=0.45, alpha=0.75)
        low = min(float(sx.min()), float(sy.min()))
        high = max(float(sx.max()), float(sy.max()))
        pad = max((high - low) * 0.05, 1e-12)
        limits = [low - pad, high + pad]
        ax.plot(limits, limits, color=ink, linestyle="--", linewidth=1.0,
                label="1:1")
        ax.set_xlim(limits)
        ax.set_ylim(limits)
        ax.set_xlabel("Sentaurus factor response in log10(n)", color=ink)
        ax.set_ylabel("Vela factor response in log10(n)", color=ink)
        ax.set_title(title, loc="left", color=ink)
        ax.grid(True, color="#E5E9EF", linewidth=0.8)
        ax.legend(frameon=False, loc="upper left")

    fig.suptitle("SimpleMOS M29 — BGN x SRH causal factorial audit",
                 x=0.07, y=0.985, ha="left", color=ink,
                 fontsize=15, fontweight="bold")
    fig.text(0.07, 0.95,
             "n23, Vd = 0.05 V, Vg = 0.05 V; HFS off; fixed M28 continuation route; 942 common silicon nodes",
             ha="left", color=grey, fontsize=10)
    fig.tight_layout(rect=(0.04, 0.035, 0.99, 0.925), h_pad=2.3, w_pad=2.0)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=190, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
