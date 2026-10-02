#!/usr/bin/env python3
"""Plot the SimpleMOS M22 n23 deep-off HFS attribution results."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = (REPO / "reference_tcad/simplemos_sentaurus2022"
        / "n23_hfs_deep_off")
FIGURES = REPO / "docs/validation/figures/simplemos_m22"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    summary = rows(DATA / "m22_summary.csv")
    edges = rows(DATA / "m22_key_state_self_consistent_edges.csv")
    gates = np.array([float(row["gate_voltage_V"]) for row in summary])

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.8), constrained_layout=True)
    ax = axes[0, 0]
    for solver, variant, marker, color in (
            ("sentaurus", "full", "o", "#2458a6"),
            ("vela", "full", "s", "#d34b37"),
            ("sentaurus", "no_hfs", "o", "#4f83cc"),
            ("vela", "no_hfs", "s", "#e58b79")):
        values = [float(row[f"{solver}_{variant}_current_A_per_um"])
                  for row in summary]
        label = f"{solver.capitalize()} {'HFS on' if variant == 'full' else 'HFS off'}"
        ax.semilogy(gates, values, marker=marker,
                    linestyle="-" if variant == "full" else "--",
                    color=color, label=label)
    ax.set_title("(a) n23 deep-off Id-Vg")
    ax.set_xlabel("Gate voltage (V)")
    ax.set_ylabel("|Id| (A/um)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    full_error = [float(row["cross_solver_full_error_dex"]) for row in summary]
    off_error = [float(row["cross_solver_no_hfs_error_dex"]) for row in summary]
    ax.plot(gates, full_error, "o-", label="HFS on", color="#2458a6")
    ax.plot(gates, off_error, "s--", label="HFS off", color="#d34b37")
    ax.axvline(0.05, color="0.5", linewidth=1, linestyle=":")
    ax.set_title("(b) Cross-solver log error")
    ax.set_xlabel("Gate voltage (V)")
    ax.set_ylabel("|log10(Vela/Sentaurus)| (dex)")
    ax.grid(True, alpha=0.25)
    ax.legend()

    ax = axes[1, 0]
    sent_effect = [float(row["sentaurus_hfs_self_consistent_effect_dex"])
                   for row in summary]
    vela_effect = [float(row["vela_hfs_self_consistent_effect_dex"])
                   for row in summary]
    frozen = [float(row["full_state_frozen_hfs_effect_dex"])
              for row in summary]
    ax.plot(gates, sent_effect, "o-", label="Sentaurus self-consistent")
    ax.plot(gates, vela_effect, "s-", label="Vela self-consistent")
    ax.plot(gates, frozen, "^-", label="Vela frozen direct")
    ax.axhline(0.0, color="0.5", linewidth=1)
    ax.set_title("(c) HFS effect: state feedback vs direct")
    ax.set_xlabel("Gate voltage (V)")
    ax.set_ylabel("log10(Id HFS-on / HFS-off) (dex)")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    edge_ids = [row["edge_id"] for row in edges]
    condition = [max(float(row["full_sg_cancellation_condition"]),
                     float(row["no_hfs_sg_cancellation_condition"]))
                 for row in edges]
    current_delta = [abs(float(
        row["electron_current_delta_full_minus_no_hfs_A_per_um"]))
                     for row in edges]
    x = np.arange(len(edges))
    ax.bar(x - 0.18, condition, width=0.36, color="#6a51a3",
           label="SG cancellation condition")
    ax.set_yscale("log")
    ax.set_ylabel("SG cancellation condition")
    ax.set_xticks(x, edge_ids, rotation=45)
    ax.set_xlabel("Drain-cut edge id")
    ax2 = ax.twinx()
    ax2.bar(x + 0.18, current_delta, width=0.36, color="#f28e2b",
            label="|self-consistent current delta|")
    ax2.set_yscale("log")
    ax2.set_ylabel("|Delta edge electron current| (A/um)")
    ax.set_title("(d) Vg=0.05 V drain-cut conditioning")
    ax.grid(True, axis="y", alpha=0.2)
    handles1, labels1 = ax.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(handles1 + handles2, labels1 + labels2,
              loc="lower left", fontsize=7)

    fig.suptitle("SimpleMOS M22: n23 HFS deep-off attribution", fontsize=14)
    output = FIGURES / "simplemos_m22_n23_hfs_deep_off.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
