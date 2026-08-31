#!/usr/bin/env python3
"""Plot the SimpleMOS M23 SG numerical-resolvability audit."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/numerical_resolvability"
FIGURES = REPO / "docs/validation/figures/simplemos_m23"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    edges = rows(DATA / "m23_edge_resolvability.csv")
    states = rows(DATA / "m23_state_summary.csv")
    paths = rows(DATA / "m23_solver_path_reproducibility.csv")
    labels = [f"{row['case'].replace('_stable_weak_inversion', '')}\n{row['state_variant']}"
              for row in states]
    colors = {"m21_stable_weak_inversion": "#2458a6",
              "m22_deep_off": "#d97721"}

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.9),
                             constrained_layout=True)

    ax = axes[0, 0]
    for case in colors:
        selected = [row for row in edges if row["case"] == case]
        ax.scatter([max(float(row["sg_cancellation_condition"]), 1.0)
                    for row in selected],
                   [max(float(row["production_vs_decimal_relative_error"]),
                        1.0e-18) for row in selected],
                   s=34, alpha=0.8, color=colors[case],
                   marker="o" if case.startswith("m21") else "s",
                   label=case.replace("_stable_weak_inversion", ""))
    ax.axhline(1.0e-12, color="#444444", linestyle="--", linewidth=1,
               label="M23 resolved threshold")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title("(a) Production SG evaluation accuracy")
    ax.set_xlabel("Raw subtractive condition number")
    ax.set_ylabel("|production / Decimal-100 - 1|")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    x = np.arange(len(states))
    stable = [max(float(row["maximum_production_vs_decimal_relative_error"]),
                  1.0e-18) for row in states]
    subtractive = [max(float(row[
        "maximum_subtractive_vs_decimal_stable_relative_error"]), 1.0e-18)
                   for row in states]
    ax.bar(x - 0.19, stable, 0.38, color="#2458a6",
           label="Production factorized")
    ax.bar(x + 0.19, subtractive, 0.38, color="#d97721", hatch="//",
           label="Raw subtractive replay")
    ax.set_yscale("log")
    ax.set_xticks(x, labels, fontsize=8)
    ax.set_ylabel("Maximum relative error")
    ax.set_title("(b) Stable and subtractive formulations")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    prod = np.array([max(float(row["production_vs_decimal_relative_error"]),
                         1.0e-18) for row in edges])
    ulp = np.array([max(float(row["one_ulp_flux_relative_envelope"]),
                        1.0e-18) for row in edges])
    order = np.argsort(prod)
    ax.plot(np.arange(len(edges)), prod[order], "o-", markersize=3,
            color="#2458a6", label="Production / Decimal-100")
    ax.plot(np.arange(len(edges)), ulp[order], "s--", markersize=3,
            color="#d97721", label="One-ULP endpoint envelope")
    ax.set_yscale("log")
    ax.set_xlabel("Key edge, sorted by production error (36 edges)")
    ax.set_ylabel("Relative flux change")
    ax.set_title("(c) Stored-state and formula resolution")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[1, 1]
    controls = ["baseline_repeat", "tight_tolerance",
                "same_state_reclosure", "reverse_from_vg_0p1"]
    x = np.arange(len(controls))
    for offset, variant, color, hatch in (
            (-0.18, "full", "#2458a6", None),
            (0.18, "no_hfs", "#d97721", "//")):
        by_control = {row["control"]: row for row in paths
                      if row["variant"] == variant}
        values = [max(abs(float(by_control[name][
            "current_log10_change_dex"])), 1.0e-12) for name in controls]
        ax.bar(x + offset, values, 0.36, color=color, hatch=hatch,
               label="HFS on" if variant == "full" else "HFS off")
    ax.axhline(0.09977223082251109, color="#444444", linestyle="--",
               linewidth=1, label="HFS self-consistent effect")
    ax.set_yscale("log")
    ax.set_xticks(x, ["repeat", "tight tol.", "reclosure", "reverse"],
                  fontsize=8)
    ax.set_ylabel("|current change| (dex)")
    ax.set_title("(d) Vg=0.05 V solver-path reproducibility")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    fig.suptitle("SimpleMOS M23: SG numerical resolvability", fontsize=14)
    output = FIGURES / "simplemos_m23_numerical_resolvability.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
