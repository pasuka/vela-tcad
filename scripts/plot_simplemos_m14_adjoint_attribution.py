#!/usr/bin/env python3
"""Render static M14 evidence figures from the frozen portable CSV files."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/adjoint_attribution"
FIGURES = REPO / "docs/validation/figures/simplemos_m14"
BLUE = "#3568A8"
ORANGE = "#D9822B"
GOLD = "#B08D18"
PINK = "#B85C82"
INK = "#252A31"
GRID = "#D8DDE5"


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "axes.titleweight": "bold",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def finish(fig: plt.Figure, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / name, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def weak_region() -> None:
    data = rows("m14_weak_region_80_points.csv")
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    palette = {"n21": BLUE, "n22": ORANGE, "n23": GOLD, "n24": PINK}
    markers = {"n21": "o", "n22": "s", "n23": "^", "n24": "D"}
    for device in ("n21", "n22", "n23", "n24"):
        selected = sorted((row for row in data if row["device"] == device),
                          key=lambda row: float(row["gate_voltage_V"]))
        ax.plot([float(row["gate_voltage_V"]) for row in selected],
                [float(row["signed_log10_ratio_dex"]) for row in selected],
                color=palette[device], marker=markers[device], markersize=3.3,
                linewidth=1.5, label=device)
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_title("High-NWell, low-drain terminal-current residual", pad=30)
    ax.text(0.0, 1.01, "Vd = 0.05 V; 80 points; positive means Vela > Sentaurus",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.set_xlabel("Gate voltage (V)")
    ax.set_ylabel("signed log10(Id,Vela / Id,Sentaurus) (dex)")
    ax.set_xlim(0.0, 0.95)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(ncol=4, frameon=False, loc="upper left")
    finish(fig, "m14_weak_region_80_points.png")


def field_cancellation() -> None:
    data = rows("m14_substitution_summary.csv")
    states = ["n21_vd_0p05_vg_0", "n21_vd_0p05_vg_0p8"]
    labels = ["Vg=0 V", "Vg=0.8 V"]
    variants = ["field_psi", "field_phin", "field_phip"]
    variant_labels = ["psi", "phin", "phip"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
    for ax, state, label in zip(axes, states, labels):
        selected = {row["variant"]: row for row in data if row["state"] == state}
        direct = [abs(float(selected[name]["direct_exact_relative_to_current"]))
                  for name in variants]
        relaxation = [abs(float(selected[name]["adjoint_relaxation_relative_to_current"]))
                      for name in variants]
        x = np.arange(len(variants))
        width = 0.36
        ax.bar(x - width / 2, direct, width, color=BLUE, edgecolor=INK,
               linewidth=0.5, label="Frozen direct")
        ax.bar(x + width / 2, relaxation, width, color="white", edgecolor=ORANGE,
               linewidth=1.3, hatch="//", label="Adjoint relaxation")
        ax.set_yscale("log")
        ax.set_xticks(x, variant_labels)
        ax.set_title(label)
        ax.grid(axis="y", color=GRID, linewidth=0.7, which="both")
    axes[0].set_ylabel("absolute current change / |baseline Id|")
    axes[0].legend(frameon=False, loc="upper left")
    fig.suptitle("Field substitution: direct response and self-consistent cancellation",
                 fontweight="bold")
    fig.text(0.5, 0.92, "n21, Vd=0.05 V; deep-off Vg=0.05 V omitted because sub-fA conditioning dominates",
             ha="center", color="#5B6573")
    fig.tight_layout(rect=(0, 0, 1, 0.89))
    finish(fig, "m14_field_cancellation.png")


def region_response() -> None:
    data = rows("m14_substitution_summary.csv")
    selected = {row["variant"]: row for row in data
                if row["state"] == "n21_vd_0p05_vg_0p8"}
    variants = ["region_source", "region_channel", "region_drain", "region_body"]
    labels = ["Source", "Channel", "Drain", "Body"]
    direct = [100.0 * float(selected[name]["direct_exact_relative_to_current"])
              for name in variants]
    relaxation = [100.0 * float(selected[name]["adjoint_relaxation_relative_to_current"])
                  for name in variants]
    x = np.arange(len(variants))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    ax.bar(x - width / 2, direct, width, color=BLUE, edgecolor=INK,
           linewidth=0.5, label="Frozen direct")
    ax.bar(x + width / 2, relaxation, width, color="white", edgecolor=ORANGE,
           linewidth=1.3, hatch="//", label="Adjoint relaxation")
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("signed current change / baseline Id (%)")
    ax.set_title("Spatial substitution response at n21, Vd=0.05 V, Vg=0.8 V",
                 pad=30)
    ax.text(0.0, 1.01, "Mutually exclusive partition; Sentaurus fields substituted into the Vela state",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(frameon=False, loc="lower left")
    finish(fig, "m14_region_response_n21_vg0p8.png")


def adjoint_quality() -> None:
    data = rows("m14_state_summary.csv")
    values = [float(row["adjoint_relative_residual"]) for row in data]
    labels = [row["state"].replace("_vd_", "\nVd=").replace("_vg_", ", Vg=")
              for row in data]
    fig, ax = plt.subplots(figsize=(9.3, 4.2))
    ax.bar(np.arange(len(values)), values, color=BLUE, edgecolor=INK, linewidth=0.4)
    ax.axhline(1.0e-12, color=ORANGE, linewidth=1.2, linestyle="--",
               label="acceptance ceiling 1e-12")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(values)), labels, rotation=65, ha="right", fontsize=7)
    ax.set_ylabel("||J^T lambda - dI/dx||2 / ||dI/dx||2")
    ax.set_title("Terminal-current adjoint linear-solve closure", pad=30)
    ax.text(0.0, 1.01, "16 paired Vela/Sentaurus representative states",
            transform=ax.transAxes, va="bottom", color="#5B6573")
    ax.grid(axis="y", color=GRID, linewidth=0.7, which="both")
    ax.legend(frameon=False, loc="upper left")
    finish(fig, "m14_adjoint_quality.png")


def main() -> None:
    style()
    weak_region()
    field_cancellation()
    region_response()
    adjoint_quality()


if __name__ == "__main__":
    main()
