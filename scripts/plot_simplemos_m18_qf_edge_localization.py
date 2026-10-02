#!/usr/bin/env python3
"""Render static M18 QF edge-localization evidence figures."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, SymLogNorm
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/qf_edge_localization"
FIGURES = REPO / "docs/validation/figures/simplemos_m18"
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


def bucket_support() -> None:
    data = [row for row in rows("m18_bucket_factor_contributions.csv")
            if row["state"] == KEY_STATE
            and row["factor"] == "qf_log_imbalance"]
    order = ["drain_contact_cut", "internal_channel", "internal_drain",
             "internal_source", "source_contact_cut", "substrate_contact_cut",
             "internal_body", "gate_contact_cut"]
    indexed = {row["bucket"]: row for row in data}
    values = [100.0 * float(indexed[name]["absolute_edge_support_fraction"])
              for name in order]
    labels = [name.replace("_", " ").title() for name in order]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    y = np.arange(len(order))
    bars = ax.barh(y, values, color=[BLUE if "contact" in name else GOLD
                                     for name in order],
                   edgecolor=INK, linewidth=0.5)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(0.0, 88.0)
    ax.set_xlabel("fraction of absolute edge QF response (%)")
    ax.set_title("QF edge-response support by finite-volume partition", pad=30)
    ax.text(0.0, 1.015,
            "n21, Vd=0.05 V, Vg=0.8 V; absolute adjoint-weighted Shapley support",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, values):
        label = f"{value:.2f}%" if value >= 0.01 else f"{value:.1e}%"
        ax.text(value + 0.8, bar.get_y() + bar.get_height() / 2,
                label, va="center", ha="left")
    finish(fig, "m18_qf_bucket_support_n21_vg0p8.png")


def edge_map() -> None:
    data = rows("m18_key_state_edge_contributions.csv")
    lateral = np.array([float(row["midpoint_y_um"]) for row in data])
    depth = np.array([float(row["midpoint_x_um"]) for row in data])
    response = np.array([float(row["qf_log_imbalance_A_per_um"])
                         for row in data])
    contact = np.array(["contact_cut" in row["bucket"] for row in data])
    nonzero = np.abs(response[np.nonzero(response)])
    linthresh = max(float(np.percentile(nonzero, 25)) if nonzero.size else 1e-20,
                    1e-20)
    vmax = max(float(np.max(np.abs(response))), 10.0 * linthresh)
    cmap = LinearSegmentedColormap.from_list(
        "m18_signed", [GOLD, "#F4F1E4", "#FFFFFF", "#E8EFF8", BLUE])
    norm = SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax)
    sizes = 7.0 + 45.0 * np.sqrt(np.abs(response) / vmax)
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    internal = ~contact
    scatter = ax.scatter(lateral[internal], depth[internal],
                         c=response[internal], s=sizes[internal], cmap=cmap,
                         norm=norm, marker="o", linewidths=0.0, alpha=0.8)
    ax.scatter(lateral[contact], depth[contact], c=response[contact],
               s=sizes[contact] + 12.0, cmap=cmap, norm=norm, marker="D",
               edgecolors=INK, linewidths=0.35)
    top = sorted(data, key=lambda row: int(row["qf_absolute_rank"]))[:2]
    annotation_offsets = {"1": (-92, 14), "2": (12, 22)}
    for row in top:
        rank = row["qf_absolute_rank"]
        ax.annotate(f"#{row['qf_absolute_rank']} edge {row['edge_id']}",
                    (float(row["midpoint_y_um"]), float(row["midpoint_x_um"])),
                    xytext=annotation_offsets[rank], textcoords="offset points",
                    fontsize=8,
                    arrowprops={"arrowstyle": "-", "color": MUTED,
                                "linewidth": 0.6})
    ax.set_xlabel("lateral coordinate y (um)")
    ax.set_ylabel("depth coordinate x (um)")
    ax.invert_yaxis()
    ax.set_title("Per-edge QF log-imbalance response", pad=30)
    ax.text(0.0, 1.015,
            "n21, Vd=0.05 V, Vg=0.8 V; diamonds are contact-cut edges; symmetric-log color",
            transform=ax.transAxes, va="bottom", color=MUTED)
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("adjoint-weighted response (A/um)")
    finish(fig, "m18_qf_edge_map_n21_vg0p8.png")


def pareto_support() -> None:
    data = rows("m18_key_state_edge_contributions.csv")
    top = data[:30]
    ranks = np.array([int(row["qf_absolute_rank"]) for row in top])
    values = 100.0 * np.array([
        abs(float(row["qf_log_imbalance_A_per_um"])) for row in top])
    total = sum(abs(float(row["qf_log_imbalance_A_per_um"])) for row in data)
    values /= total
    cumulative = 100.0 * np.array([
        float(row["qf_cumulative_absolute_fraction"]) for row in top])
    fig, ax = plt.subplots(figsize=(8.2, 4.7))
    ax.bar(ranks, values, color=BLUE, edgecolor=INK, linewidth=0.35,
           label="individual edge")
    ax.set_xlabel("edge rank by absolute QF response")
    ax.set_ylabel("individual support (%)")
    ax.set_xlim(0.25, 30.75)
    ax.grid(axis="y", color=GRID, linewidth=0.7)
    twin = ax.twinx()
    twin.plot(ranks, cumulative, color=ORANGE, linewidth=1.8, marker="o",
              markersize=3.0, label="cumulative")
    twin.axhline(50.0, color=MUTED, linestyle="--", linewidth=0.8)
    twin.set_ylim(0.0, 105.0)
    twin.set_ylabel("cumulative absolute support (%)")
    ax.set_title("Concentration of QF response on the leading edges", pad=30)
    ax.text(0.0, 1.015,
            "n21, Vd=0.05 V, Vg=0.8 V; edge #1 supplies 66.59%, top two supply 79.41%",
            transform=ax.transAxes, va="bottom", color=MUTED)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = twin.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, frameon=False,
              loc="lower right")
    finish(fig, "m18_qf_edge_pareto_n21_vg0p8.png")


def feedback_ledger() -> None:
    data = {row["component"]: row for row in rows("m18_feedback_ledger.csv")
            if row["state"] == KEY_STATE}
    components = ["mobility", "sg_secant_conductance",
                  "qf_log_imbalance", "poisson"]
    labels = ["Mobility", "SG secant conductance", "QF log imbalance",
              "Poisson feedback"]
    values = [100.0 * float(data[name]["relative_to_baseline_current"])
              for name in components]
    fig, ax = plt.subplots(figsize=(8.0, 4.5))
    y = np.arange(len(components))
    bars = ax.barh(y, values,
                   color=[GOLD if value < 0 else BLUE for value in values],
                   edgecolor=INK, linewidth=0.5)
    ax.axvline(0.0, color=INK, linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("adjoint-weighted response / baseline Id (%)")
    ax.set_title("Self-consistent feedback ledger", pad=30)
    srh_percent = 100.0 * float(data["srh"]["relative_to_baseline_current"])
    ax.text(0.0, 1.015,
            f"n21, Vd=0.05 V, Vg=0.8 V; direct SRH response = {srh_percent:.2e}%",
            transform=ax.transAxes, va="bottom", color=MUTED)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.set_xlim(min(values) - 2.2, max(values) + 2.2)
    for bar, value in zip(bars, values):
        ax.text(value + (0.25 if value >= 0 else -0.2),
                bar.get_y() + bar.get_height() / 2,
                f"{value:+.3f}%", va="center",
                ha="left" if value >= 0 else "right")
    finish(fig, "m18_feedback_ledger_n21_vg0p8.png")


def contact_fraction_matrix() -> None:
    data = rows("m18_bucket_factor_contributions.csv")
    contact = {"source_contact_cut", "drain_contact_cut",
               "substrate_contact_cut", "gate_contact_cut"}
    totals: dict[tuple[str, float, float], float] = {}
    for row in data:
        if row["factor"] != "qf_log_imbalance" or row["bucket"] not in contact:
            continue
        key = (row["device"], float(row["drain_voltage_V"]),
               float(row["gate_voltage_V"]))
        totals[key] = totals.get(key, 0.0) + float(
            row["absolute_edge_support_fraction"])
    row_keys = [("n17", 0.05), ("n17", 1.0),
                ("n21", 0.05), ("n21", 1.0)]
    gates = [0.0, 0.05, 0.8, 2.5]
    matrix = 100.0 * np.array([[totals[(device, drain, gate)]
                                for gate in gates]
                               for device, drain in row_keys])
    fig, ax = plt.subplots(figsize=(7.6, 4.5))
    image = ax.imshow(matrix, cmap="Blues", vmin=0.0, vmax=100.0,
                      aspect="auto")
    ax.set_xticks(np.arange(len(gates)), [f"{gate:g}" for gate in gates])
    ax.set_yticks(np.arange(len(row_keys)),
                  [f"{device}, Vd={drain:g} V" for device, drain in row_keys])
    ax.set_xlabel("Vg (V)")
    ax.set_title("Contact-cut share of absolute QF edge support", pad=30)
    ax.text(0.0, 1.015,
            "16 paired states; percentage of absolute adjoint-weighted QF response",
            transform=ax.transAxes, va="bottom", color=MUTED)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            value = matrix[row, column]
            ax.text(column, row, f"{value:.1f}%", ha="center", va="center",
                    color="white" if value > 55.0 else INK)
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("contact-cut support (%)")
    finish(fig, "m18_qf_contact_fraction_16_states.png")


def main() -> None:
    style()
    bucket_support()
    edge_map()
    pareto_support()
    feedback_ledger()
    contact_fraction_matrix()


if __name__ == "__main__":
    main()
