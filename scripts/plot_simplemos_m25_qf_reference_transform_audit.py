#!/usr/bin/env python3
"""Plot the SimpleMOS M25 QF-reference transform audit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "qf_reference_transform_audit")
FIGURES = REPO / "docs/validation/figures/simplemos_m25"
BLUE = "#2458a6"
ORANGE = "#d97721"
INK = "#30343b"
LIGHT_BLUE = "#9db8df"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    report = json.loads((DATA / "m25_qf_reference_transform_audit_report.json")
                        .read_text(encoding="utf-8"))
    edge_rows = rows(DATA / "m25_edge_reference_transform_ledger.csv")
    pair_rows = rows(DATA / "m25_variant_edge_delta_ledger.csv")
    summaries = rows(DATA / "m25_class_summary.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.9),
                             constrained_layout=True)

    ax = axes[0, 0]
    topology = report["topology"]
    labels = ["Drain cut", "Reference\ntransition", "Overlap"]
    values = [topology["drain_cut_edge_count"],
              topology["reference_transition_edge_count"],
              topology["drain_cut_reference_transition_edge_count"]]
    bars = ax.bar(labels, values, color=[BLUE, ORANGE, "#d9dde3"],
                  edgecolor=INK, linewidth=0.7)
    ax.bar_label(bars, padding=3)
    ax.set_ylabel("Unique edge count")
    ax.set_title("(a) Audited edge populations")
    ax.set_ylim(0, 58)
    ax.grid(True, axis="y", alpha=0.2)

    ax = axes[0, 1]
    x = np.arange(4)
    labels = ["Endpoint\ntransform", "Physical QF\ndrop",
              "SG log\n(state)", "Factorized\nflux"]
    for offset, variant, color, hatch in (
            (-0.18, "full", BLUE, None),
            (0.18, "no_hfs", ORANGE, "//")):
        group = [row for row in summaries
                 if row["variant"] == variant
                 and row["edge_class"] == "reference_transition"][0]
        values = [float(group["maximum_endpoint_transform_error_V"]),
                  float(group["maximum_physical_drop_closure_error_V"]),
                  float(group["maximum_sg_log_state_closure_error"]),
                  float(group["maximum_factorized_flux_relative_error"])]
        ax.bar(x + offset, values, 0.36, color=color, hatch=hatch,
               edgecolor=INK, linewidth=0.6,
               label="HFS on" if variant == "full" else "HFS off")
    ax.set_yscale("log")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Absolute error or relative flux error")
    ax.set_title("(b) Cross-reference numerical closure")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    transition = [row for row in edge_rows
                  if row["variant"] == "full"
                  and row["reference_transition"].lower() == "true"]
    transition.sort(key=lambda row: abs(float(row["physical_qf_drop_V"])))
    rank = np.arange(1, len(transition) + 1)
    ref_jump = [abs(float(row["reference_jump_V"])) for row in transition]
    increment = [max(abs(float(row["naive_increment_only_drop_V"])), 1e-20)
                 for row in transition]
    physical = [max(abs(float(row["physical_qf_drop_V"])), 1e-20)
                for row in transition]
    ax.plot(rank, ref_jump, color=INK, linestyle="--", linewidth=1.4,
            label="|reference jump|")
    ax.plot(rank, increment, color=ORANGE, linewidth=1.4,
            label="|increment-only drop|")
    ax.plot(rank, physical, color=BLUE, marker="o", markersize=2.6,
            linewidth=1.2, label="|physical QF drop|")
    ax.set_yscale("log")
    ax.set_xlabel("Reference-transition edge (sorted by physical drop)")
    ax.set_ylabel("Magnitude (V)")
    ax.set_title("(c) Reference jump cancellation, HFS on")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    classes = ["drain_cut", "reference_transition"]
    labels = ["Drain cut", "Reference transition"]
    maxima = []
    medians = []
    for edge_class in classes:
        values = [abs(float(row["full_minus_no_hfs_qf_drop_V"]))
                  for row in pair_rows if row["edge_class"] == edge_class]
        maxima.append(max(max(values), 1e-22))
        medians.append(max(float(np.median(values)), 1e-22))
    x = np.arange(2)
    ax.bar(x - 0.18, maxima, 0.36, color=BLUE, edgecolor=INK,
           linewidth=0.6, label="maximum")
    ax.bar(x + 0.18, medians, 0.36, color=LIGHT_BLUE, hatch="//",
           edgecolor=INK, linewidth=0.6, label="median")
    ax.set_yscale("log")
    ax.set_xticks(x, labels)
    ax.set_ylabel("|HFS on - HFS off QF drop| (V)")
    ax.set_title("(d) Self-consistent QF response by edge class")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    fig.suptitle("SimpleMOS M25: QF reference transform and SG log audit",
                 fontsize=14)
    output = FIGURES / "simplemos_m25_qf_reference_transform_audit.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
