#!/usr/bin/env python3
"""Render static SimpleMOS M20 drain-adjacent continuity-row figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/drain_adjacent_row_audit"
FIGURES = REPO / "docs/validation/figures/simplemos_m20"
BLUE = "#3568A8"
GOLD = "#B08D18"
ORANGE = "#D9822B"
INK = "#252A31"
MUTED = "#5B6573"
GRID = "#D8DDE5"


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


def key_flux_balance() -> None:
    data = rows("m20_key_state_node_ledger.csv")
    labels = [row["node_id"] for row in data]
    contact = np.array([float(row["delta_contact_flux"]) for row in data])
    internal = np.array([float(row["delta_internal_flux"]) for row in data])
    residual = np.array([float(row["delta_raw_residual"]) for row in data])
    x = np.arange(len(data))
    width = 0.25
    fig, ax = plt.subplots(figsize=(8.4, 4.9))
    ax.bar(x - width, contact, width, color=GOLD, edgecolor=INK,
           linewidth=0.45, label="drain-contact cut flux")
    ax.bar(x, internal, width, color=BLUE, edgecolor=INK,
           linewidth=0.45, label="internal-edge flux")
    ax.bar(x + width, residual, width, color=ORANGE, edgecolor=INK,
           linewidth=0.45, label="resulting row residual")
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_yscale("symlog", linthresh=1.0e-8)
    ax.set_xticks(x, labels)
    ax.set_xlabel("drain-adjacent free node")
    ax.set_ylabel("Sentaurus-state minus Vela-state raw term")
    ax.set_title("Contact and internal SG changes nearly cancel", pad=30)
    ax.text(0.0, 1.015, "n21, Vd=0.05 V, Vg=0.8 V; unchanged Vela operator",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    finish(fig, "m20_key_state_flux_balance.png")


def key_cancellation() -> None:
    data = rows("m20_key_state_node_ledger.csv")
    labels = [row["node_id"] for row in data]
    delta_condition = [
        (abs(float(row["delta_contact_flux"])) + abs(float(row["delta_internal_flux"])))
        / max(abs(float(row["delta_raw_residual"])), 1.0e-300)
        for row in data
    ]
    state_condition = [float(row["sentaurus_flux_condition_number"]) for row in data]
    x = np.arange(len(data))
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.0, 4.7))
    ax.bar(x - width / 2, delta_condition, width, color=ORANGE,
           edgecolor=INK, linewidth=0.45, label="state-difference cancellation")
    ax.bar(x + width / 2, state_condition, width, color=BLUE,
           edgecolor=INK, linewidth=0.45, label="imported-state flux condition")
    ax.set_yscale("log")
    ax.set_xticks(x, labels)
    ax.set_xlabel("drain-adjacent free node")
    ax.set_ylabel("absolute-flux sum / net value")
    ax.set_title("The first drain layer is cancellation-conditioned", pad=30)
    ax.text(0.0, 1.015, "Node 991: 1275x delta cancellation and 7884x state-flux condition",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7, which="both")
    ax.legend(frameon=False)
    finish(fig, "m20_key_state_flux_conditioning.png")


def key_jacobian_partition() -> None:
    data = rows("m20_key_state_node_ledger.csv")
    labels = [row["node_id"] for row in data]
    contact = np.array([float(row["sentaurus_phin_contact_column_share"])
                        for row in data]) * 100.0
    free = np.array([float(row["sentaurus_phin_free_column_share"])
                     for row in data]) * 100.0
    other = 100.0 - contact - free
    x = np.arange(len(data))
    fig, ax = plt.subplots(figsize=(8.0, 4.7))
    ax.bar(x, contact, color=GOLD, edgecolor=INK, linewidth=0.45,
           label="contact-column phin")
    ax.bar(x, free, bottom=contact, color=BLUE, edgecolor=INK, linewidth=0.45,
           label="free-column phin")
    ax.bar(x, other, bottom=contact + free, color=ORANGE, edgecolor=INK,
           linewidth=0.45, label="psi + phip")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.0, 100.0)
    ax.set_xlabel("drain-adjacent free node")
    ax.set_ylabel("absolute production-Jacobian row share (%)")
    ax.set_title("Electron rows are almost purely phin coupled", pad=30)
    ax.text(0.0, 1.015, "Exact production matrix partition; contact/free refers to column constraint",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    finish(fig, "m20_key_state_jacobian_partition.png")


def main() -> None:
    style()
    key_flux_balance()
    key_cancellation()
    key_jacobian_partition()
    print(FIGURES)


if __name__ == "__main__":
    main()
