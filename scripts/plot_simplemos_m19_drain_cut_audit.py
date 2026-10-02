#!/usr/bin/env python3
"""Render static M19 drain-cut audit figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/drain_cut_audit"
FIGURES = REPO / "docs/validation/figures/simplemos_m19"
KEY_STATE = "n21_vd_0p05_vg_0p8"
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


def vg08_gap_closure() -> None:
    data = [row for row in rows("m19_state_summary.csv")
            if abs(float(row["gate_voltage_V"]) - 0.8) < 1.0e-12]
    labels = [f"{row['device']}\nVd={float(row['drain_voltage_V']):g} V"
              for row in data]
    baseline = [100.0 * (10.0 ** float(row["terminal_signed_log10_error_dex"]) - 1.0)
                for row in data]
    replay = [100.0 * (10.0 ** abs(float(row["sentaurus_state_vela_sg_error_dex"])) - 1.0)
              for row in data]
    x = np.arange(len(data))
    width = 0.35
    fig, ax = plt.subplots(figsize=(8.1, 4.7))
    ax.bar(x - width / 2, baseline, width, color=GOLD, edgecolor=INK,
           linewidth=0.45, label="self-consistent Vela")
    ax.bar(x + width / 2, replay, width, color=BLUE, edgecolor=INK,
           linewidth=0.45, label="Sentaurus state through Vela SG")
    ax.set_xticks(x, labels)
    ax.set_ylabel("absolute terminal-current error (%)")
    ax.set_title("Drain-cut state replay closes the Vg=0.8 V current gap", pad=30)
    ax.text(0.0, 1.015,
            "Unchanged Vela SG operator; imported Sentaurus fields use exact TDR node IDs",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(frameon=False)
    finish(fig, "m19_vg0p8_gap_closure.png")


def key_factor_comparison() -> None:
    state = next(row for row in rows("m19_state_summary.csv")
                 if row["state"] == KEY_STATE)
    labels = ["Contact endpoint\nstate", "Adjacent interior\nQF", "Adjacent interior\nother"]
    frozen = np.array([
        float(state["frozen_drain_contact_endpoint_state_relative_to_Id"]),
        float(state["frozen_drain_adjacent_interior_phin_relative_to_Id"]),
        float(state["frozen_drain_adjacent_interior_other_state_relative_to_Id"]),
    ]) * 100.0
    x = np.arange(3)
    fig, ax = plt.subplots(figsize=(8.0, 4.7))
    colors = [BLUE if value >= 0.0 else GOLD for value in frozen]
    bars = ax.bar(x, frozen, color=colors, edgecolor=INK, linewidth=0.5)
    ax.axhline(0.0, color=INK, linewidth=0.8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("exact frozen current change / Vela Id (%)")
    ax.set_title("Frozen drain-cut endpoint-state attribution", pad=30)
    ax.text(0.0, 1.015,
            "n21, Vd=0.05 V, Vg=0.8 V; three-factor Shapley allocation",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, frozen):
        ax.text(bar.get_x() + bar.get_width() / 2,
                value + (0.35 if value >= 0 else -0.35),
                f"{value:+.4f}%", ha="center",
                va="bottom" if value >= 0 else "top")
    ax.text(1.98, min(frozen) * 0.60,
            "M18 self-consistent QF response\n"
            f"= {100.0 * float(state['m18_drain_cut_qf_response_relative_to_Id']):+.4f}%",
            ha="right", va="center", color=ORANGE)
    finish(fig, "m19_key_state_frozen_factors.png")


def key_edge_concentration() -> None:
    data = rows("m19_key_state_drain_cut_edges.csv")
    data.sort(key=lambda row: abs(float(
        row["drain_adjacent_interior_phin_direct_A_per_um"])), reverse=True)
    values = np.array([float(row["drain_adjacent_interior_phin_direct_A_per_um"])
                       for row in data])
    total = sum(abs(value) for value in values)
    support = 100.0 * np.abs(values) / max(total, 1.0e-300)
    labels = [f"edge {row['edge_id']}" for row in data]
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    bars = ax.bar(np.arange(len(data)), support, color=BLUE,
                  edgecolor=INK, linewidth=0.45)
    ax.set_xticks(np.arange(len(data)), labels, rotation=35, ha="right")
    ax.set_ylabel("absolute adjacent-QF contribution share (%)")
    ax.set_title("Drain-cut adjacent-QF response is concentrated on two edges", pad=30)
    ax.text(0.0, 1.015,
            "n21, Vd=0.05 V, Vg=0.8 V; exact frozen cut-current attribution",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, support):
        if value >= 1.0:
            ax.text(bar.get_x() + bar.get_width() / 2, value + 1.0,
                    f"{value:.1f}%", ha="center", va="bottom")
    finish(fig, "m19_key_state_edge_concentration.png")


def contact_boundary_error() -> None:
    data = rows("m19_state_summary.csv")
    values = [max(float(row["max_vela_contact_phin_bias_error_V"]),
                  float(row["max_sentaurus_contact_phin_bias_error_V"]))
              for row in data]
    labels = [row["state"].replace("_vd_", "\nVd=").replace("_vg_", " Vg=")
              for row in data]
    fig, ax = plt.subplots(figsize=(10.2, 4.5))
    ax.semilogy(np.arange(len(data)), np.maximum(values, 1.0e-18), "o-",
                color=BLUE, linewidth=1.1, markersize=4)
    ax.axhline(1.0e-12, color=ORANGE, linestyle="--", linewidth=1.0,
               label="M19 acceptance")
    ax.set_xticks(np.arange(len(data)), labels, rotation=55, ha="right")
    ax.set_ylabel("max |phin(contact)-Vd| (V)")
    ax.set_title("Drain electron quasi-Fermi boundary condition agrees in all states", pad=30)
    ax.text(0.0, 1.015, "Maximum observed error is at floating-point roundoff",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    ax.legend(frameon=False)
    finish(fig, "m19_contact_boundary_error.png")


def main() -> None:
    style()
    vg08_gap_closure()
    key_factor_comparison()
    key_edge_concentration()
    contact_boundary_error()


if __name__ == "__main__":
    main()
