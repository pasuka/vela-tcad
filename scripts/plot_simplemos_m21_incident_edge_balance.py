#!/usr/bin/env python3
"""Render static SimpleMOS M21 incident-edge audit figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/incident_edge_balance"
FIGURES = REPO / "docs/validation/figures/simplemos_m21"
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


def edge_flux_balance() -> None:
    data = rows("m21_key_state_incident_edges.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), sharey=True)
    for ax, node in zip(axes, (991, 992)):
        incident = [row for row in data
                    if node in (int(row["node0"]), int(row["node1"]))]
        incident.sort(key=lambda row: abs(float(row["electron_flux_difference"])),
                      reverse=True)
        values = []
        labels = []
        colors = []
        for row in incident:
            value = float(row["electron_flux_difference"])
            if int(row["node1"]) == node:
                value = -value
            values.append(value)
            labels.append(row["edge_id"])
            colors.append(GOLD if row["edge_class"] == "drain_contact_cut" else BLUE)
        x = np.arange(len(incident))
        ax.bar(x, values, color=colors, edgecolor=INK, linewidth=0.45)
        ax.axhline(0.0, color=INK, linewidth=0.8)
        ax.set_xticks(x, labels, rotation=45, ha="right")
        ax.set_title(f"node {node}")
        ax.grid(axis="y", color=GRID, linewidth=0.7)
    axes[0].set_ylabel("node-oriented imported minus Vela edge flux")
    fig.suptitle("A contact-edge change is nearly cancelled by an internal edge", y=1.05,
                 fontsize=13, fontweight="bold")
    fig.text(0.5, 0.98, "n21, Vd=0.05 V, Vg=0.8 V; gold = drain cut, blue = internal",
             ha="center", color=MUTED)
    finish(fig, "m21_key_edge_flux_balance.png")


def input_changes() -> None:
    selected_ids = ("2107", "2094", "2088", "2114")
    index = {row["edge_id"]: row for row in rows("m21_key_state_incident_edges.csv")}
    data = [index[item] for item in selected_ids]
    labels = [f"{row['edge_id']}\n{row['edge_class'].replace('_', ' ')}" for row in data]
    def percent_from_dex(row: dict[str, str], column: str) -> float:
        return 100.0 * abs(10.0 ** float(row[column]) - 1.0)
    values = {
        "phin drop": [100.0 * abs(float(row["sentaurus_phin_drop_V"])
                                    / float(row["baseline_phin_drop_V"]) - 1.0)
                      for row in data],
        "TCV drive": [percent_from_dex(row, "transport_cell_vector_drive_change_dex")
                      for row in data],
        "mobility": [percent_from_dex(row, "mobility_change_dex") for row in data],
        "density": [percent_from_dex(row, "maximum_endpoint_density_change_dex")
                    for row in data],
        "Bernoulli": [percent_from_dex(row, "maximum_bernoulli_weight_change_dex")
                      for row in data],
    }
    x = np.arange(len(data))
    width = 0.15
    colors = (ORANGE, GOLD, BLUE, "#6C8E5B", "#8B6FA8")
    fig, ax = plt.subplots(figsize=(9.0, 5.0))
    for index_value, ((name, series), color) in enumerate(zip(values.items(), colors)):
        ax.bar(x + (index_value - 2) * width, series, width, color=color,
               edgecolor=INK, linewidth=0.4, label=name)
    ax.set_yscale("log")
    ax.set_xticks(x, labels)
    ax.set_ylabel("absolute imported-state change (%)")
    ax.set_title("QF drop and TCV drive move by ~13%; mobility barely moves", pad=30)
    ax.text(0.0, 1.015, "Production transport_cell_vector drive; four leading incident edges",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7, which="both")
    ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.19))
    finish(fig, "m21_key_edge_input_changes.png")


def factor_shares() -> None:
    data = rows("m21_key_state_incident_edges.csv")
    factors = ("mobility", "sg_secant_conductance", "qf_log_imbalance")
    values = [sum(abs(float(row[f"{factor}_flux_contribution"])) for row in data)
              for factor in factors]
    total = sum(values)
    shares = np.array(values) / total
    labels = ["Mobility", "SG secant\nconductance", "QF log\nimbalance"]
    fig, ax = plt.subplots(figsize=(7.6, 4.7))
    bars = ax.bar(labels, shares * 100.0, color=[BLUE, GOLD, ORANGE],
                  edgecolor=INK, linewidth=0.5)
    ax.set_yscale("log")
    ax.set_ylabel("absolute edge-flux contribution share (%)")
    ax.set_title("Incident-edge response is the QF-imbalance factor", pad=30)
    ax.text(0.0, 1.015, "Eleven edges incident on nodes 991/992; exact three-factor Shapley game",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="y", color=GRID, linewidth=0.7, which="both")
    for bar, share in zip(bars, shares):
        ax.text(bar.get_x() + bar.get_width() / 2, share * 100.0 * 1.25,
                f"{share * 100.0:.5g}%", ha="center", va="bottom")
    finish(fig, "m21_key_factor_shares.png")


def main() -> None:
    style()
    edge_flux_balance()
    input_changes()
    factor_shares()
    print(FIGURES)


if __name__ == "__main__":
    main()
