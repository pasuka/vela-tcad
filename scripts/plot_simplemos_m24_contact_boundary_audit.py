#!/usr/bin/env python3
"""Plot the SimpleMOS M24 drain-contact boundary audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/contact_boundary_audit"
FIGURES = REPO / "docs/validation/figures/simplemos_m24"
BLUE = "#2458a6"
ORANGE = "#d97721"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    contacts = rows(DATA / "m24_contact_node_ledger.csv")
    edges = rows(DATA / "m24_drain_cut_edge_ledger.csv")
    nodes = rows(DATA / "m24_first_layer_continuity_ledger.csv")
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.8),
                             constrained_layout=True)

    ax = axes[0, 0]
    contact_ids = sorted({int(row["node_id"]) for row in contacts})
    x = np.arange(len(contact_ids))
    for offset, variant, color, hatch in (
            (-0.18, "full", BLUE, None),
            (0.18, "no_hfs", ORANGE, "//")):
        by_node = {int(row["node_id"]): row for row in contacts
                   if row["variant"] == variant}
        ax.bar(x + offset, [float(by_node[node]["electron_qf_bias_error_V"])
                            for node in contact_ids], 0.36, color=color,
               hatch=hatch, label="HFS on" if variant == "full" else "HFS off")
    ax.axhline(1.0e-12, color="#444444", linestyle="--", linewidth=1,
               label="acceptance limit")
    ax.set_yscale("log")
    ax.set_xticks(x, contact_ids)
    ax.set_xlabel("Drain contact node")
    ax.set_ylabel("|electron QF - 0.05 V| (V)")
    ax.set_title("(a) Ohmic contact QF target")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    edge_ids = sorted({int(row["edge_id"]) for row in edges})
    x = np.arange(len(edge_ids))
    for offset, variant, color, hatch in (
            (-0.18, "full", BLUE, None),
            (0.18, "no_hfs", ORANGE, "//")):
        by_edge = {int(row["edge_id"]): row for row in edges
                   if row["variant"] == variant}
        values = [max(abs(float(by_edge[edge]["electron_current_A_per_um"])),
                      1.0e-22) for edge in edge_ids]
        ax.bar(x + offset, values, 0.36, color=color, hatch=hatch,
               label="HFS on" if variant == "full" else "HFS off")
    ax.set_yscale("log")
    ax.set_xticks(x, edge_ids, rotation=45)
    ax.set_xlabel("Drain-cut edge")
    ax.set_ylabel("|electron current| (A/um)")
    ax.set_title("(b) Drain-cut current support")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    node_ids = [int(row["node_id"]) for row in nodes]
    qf_delta = [max(abs(float(row["qf_increment_delta_V"])), 1.0e-22)
                for row in nodes]
    ax.bar(np.arange(len(nodes)), qf_delta, color=BLUE,
           label="|QF increment delta|")
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(nodes)), node_ids)
    ax.set_xlabel("First-layer free node")
    ax.set_ylabel("|HFS-on - HFS-off| (V)")
    ax2 = ax.twinx()
    shares = [float(row["full_contact_column_share"]) for row in nodes]
    ax2.plot(np.arange(len(nodes)), shares, "s--", color=ORANGE,
             label="contact Jacobian-column share")
    ax2.set_ylim(0.0, 0.5)
    ax2.set_ylabel("Contact-column share")
    ax.set_title("(c) First-layer QF state and contact coupling")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2, fontsize=8,
              loc="lower right")

    ax = axes[1, 1]
    x = np.arange(len(nodes))
    contact_delta = [float(row["delta_contact_flux"]) for row in nodes]
    internal_delta = [float(row["delta_internal_flux"]) for row in nodes]
    ax.bar(x - 0.18, contact_delta, 0.36, color=BLUE,
           label="contact-flux delta")
    ax.bar(x + 0.18, internal_delta, 0.36, color=ORANGE, hatch="//",
           label="internal-flux delta")
    ax.axhline(0.0, color="#444444", linewidth=1)
    ax.set_xticks(x, node_ids)
    ax.set_xlabel("First-layer free node")
    ax.set_ylabel("HFS-on - HFS-off raw flux")
    ax.set_title("(d) Continuity-row flux balance")
    ax.grid(True, axis="y", alpha=0.2)
    ax.legend(fontsize=8)

    fig.suptitle("SimpleMOS M24: n23 drain-contact boundary audit",
                 fontsize=14)
    output = FIGURES / "simplemos_m24_contact_boundary_audit.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
