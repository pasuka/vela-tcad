#!/usr/bin/env python3
"""Plot the SimpleMOS M27 cross-solver HFS response audit."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/cross_solver_response_audit"
FIGURES = REPO / "docs/validation/figures/simplemos_m27"
BLUE = "#2458a6"
ORANGE = "#d97721"
GRAY = "#9aa3af"
DARK = "#343a40"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    report = json.loads(
        (DATA / "m27_cross_solver_response_audit_report.json").read_text(
            encoding="utf-8-sig"))
    node_rows = rows(DATA / "m27_node_response_ledger.csv")
    edge_rows = rows(DATA / "m27_edge_response_ledger.csv")
    terminal = report["terminal_response"]

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.0),
                             constrained_layout=True)

    ax = axes[0, 0]
    solvers = ("vela", "sentaurus")
    x = np.arange(2)
    width = 0.34
    full = [abs(float(terminal[name]["full_current_A_per_um"]))
            for name in solvers]
    off = [abs(float(terminal[name]["no_hfs_current_A_per_um"]))
           for name in solvers]
    ax.bar(x - width / 2, full, width, color=BLUE, label="HFS on")
    ax.bar(x + width / 2, off, width, facecolor="white", edgecolor=ORANGE,
           hatch="//", linewidth=1.4, label="HFS off")
    ax.set_xticks(x, ("Vela", "Sentaurus"))
    ax.set_ylabel("|drain current| (A/um)")
    ax.set_title("(a) Absolute deep-off drain current")
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    response_values = [
        abs(float(terminal["vela"]["hfs_response_A_per_um"])),
        abs(float(terminal["sentaurus"]["hfs_response_A_per_um"])),
        abs(float(terminal["hfs_induced_cross_solver_gap_change_A_per_um"])),
    ]
    labels = ("Vela\nHFS response", "Sentaurus\nHFS response",
              "HFS-induced\ngap change")
    bars = ax.bar(np.arange(3), response_values,
                  color=(BLUE, ORANGE, GRAY))
    bars[1].set_hatch("//")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(3), labels)
    ax.set_ylabel("|current change| (A/um)")
    ax.set_title("(b) HFS response agreement and residual gap")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    for index, value in enumerate(response_values):
        ax.text(index, value * 1.12, f"{value:.3e}", ha="center",
                va="bottom", fontsize=8)

    ax = axes[1, 0]
    max_response = max(
        max(abs(float(row["vela_delta_phin_V"])),
            abs(float(row["sentaurus_delta_phin_V"]))) for row in node_rows)
    selected = [
        row for row in node_rows
        if max(abs(float(row["vela_delta_phin_V"])),
               abs(float(row["sentaurus_delta_phin_V"])))
        > max_response * 1.0e-6
        or row["zone"] == "drain_first_layer"
    ]
    for zone, color, marker, label in (
            ("other", GRAY, "o", "other responsive nodes"),
            ("upstream_operator_source", ORANGE, "s", "M26 upstream source"),
            ("drain_first_layer", BLUE, "*", "drain first layer")):
        subset = [row for row in selected if row["zone"] == zone]
        ax.scatter(
            [1e3 * float(row["vela_delta_phin_V"]) for row in subset],
            [1e3 * float(row["sentaurus_delta_phin_V"]) for row in subset],
            s=55 if zone != "other" else 18, color=color, marker=marker,
            alpha=0.8 if zone != "other" else 0.45,
            edgecolors="white", linewidths=0.35, label=label)
    limit = 1.08 * 1e3 * max_response
    ax.plot([-limit, limit], [-limit, limit], "--", color=DARK,
            linewidth=1, label="1:1")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Vela HFS-on minus HFS-off electron QF (mV)")
    ax.set_ylabel("Sentaurus HFS-on minus HFS-off electron QF (mV)")
    stats = report["cross_solver_response"]["responsive_node_phin"]
    ax.set_title(f"(c) Node QF response (r={stats['pearson']:.5f})")
    ax.grid(True, alpha=0.2)
    ax.legend(fontsize=7, loc="lower right")

    ax = axes[1, 1]
    active = [row for row in edge_rows if row["active_current_edge"] == "True"]
    finite = [
        row for row in active
        if math.isfinite(float(row["vela_mobility_response_dex"]))
        and math.isfinite(float(row["sentaurus_mobility_response_dex"]))
    ]
    maximum = max(
        max(abs(float(row["vela_mobility_response_dex"])),
            abs(float(row["sentaurus_mobility_response_dex"])))
        for row in finite)
    selected_edges = [
        row for row in finite
        if max(abs(float(row["vela_mobility_response_dex"])),
               abs(float(row["sentaurus_mobility_response_dex"])))
        > maximum * 1.0e-6
    ]
    for zone, color, marker, label in (
            ("other", GRAY, "o", "other active edges"),
            ("upstream_operator_source", ORANGE, "s", "M26 upstream-source edges"),
            ("drain_first_layer", BLUE, "*", "drain-first-layer edges")):
        subset = [row for row in selected_edges if row["zone"] == zone]
        ax.scatter(
            [float(row["vela_mobility_response_dex"]) for row in subset],
            [float(row["sentaurus_mobility_response_dex"]) for row in subset],
            s=55 if zone != "other" else 20, color=color, marker=marker,
            alpha=0.8 if zone != "other" else 0.5,
            edgecolors="white", linewidths=0.35, label=label)
    limit = 1.08 * maximum
    ax.plot([-limit, limit], [-limit, limit], "--", color=DARK,
            linewidth=1, label="1:1")
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Vela mobility response (dex)")
    ax.set_ylabel("Sentaurus mobility response (dex)")
    stats = report["cross_solver_response"]["active_edge_mobility"]
    ax.set_title(f"(d) Active-edge mobility response (r={stats['pearson']:.4f})")
    ax.grid(True, alpha=0.2)
    ax.legend(fontsize=7, loc="lower right")

    fig.suptitle(
        "SimpleMOS M27: n23 cross-solver HFS response audit\n"
        "Vd=0.05 V, Vg=0.05 V; shared mesh; HFS is the only switched model",
        fontsize=13.5)
    output = FIGURES / "simplemos_m27_cross_solver_response_audit.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
