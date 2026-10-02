#!/usr/bin/env python3
"""Plot M35 signed-AverageBox response and its M33 comparison."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/signed_average_box_assembly"
M33 = REPO / "reference_tcad/simplemos_sentaurus2022/region_resolved_interface_assembly"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m35"
          / "simplemos_m35_signed_average_box_assembly.png")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    response = rows(ROOT / "m35_factorial_response.csv")
    effects = rows(ROOT / "m35_factorial_main_effects.csv")
    comparison = rows(ROOT / "m35_m33_barycentric_comparison.csv")
    m33_effects = rows(M33 / "m33_factorial_main_effects.csv")
    labels = [row["variant"].replace("p", "P").replace("_t", " T").replace("_v", " V")
              for row in response]
    x = list(range(len(response)))

    figure, axes = plt.subplots(1, 3, figsize=(14.4, 4.5))
    axes[0].bar([value - 0.18 for value in x],
                [float(row["m33_barycentric_error_dex"]) for row in comparison],
                width=0.36, color="#9ecae9", label="M33 barycentric")
    axes[0].bar([value + 0.18 for value in x],
                [float(row["m35_signed_average_box_error_dex"]) for row in comparison],
                width=0.36, color="#f58518", label="M35 signed AverageBox")
    axes[0].set_ylabel("absolute Id error (dex)")
    axes[0].set_title("Cross-solver error")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(x, [float(row["log_current_shift_vs_baseline_dex"])
                    for row in response], color="#4c78a8")
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_ylabel("log10(Id / baseline) (dex)")
    axes[1].set_title("M35 self-consistent response")

    old = {row["factor"]: float(row["main_effect_log_current_dex"])
           for row in m33_effects}
    new = {row["factor"]: float(row["main_effect_log_current_dex"])
           for row in effects}
    effect_labels = ["Poisson edge", "transport edge", "node measure"]
    old_values = [old["poisson_edge_coupling"], old["transport_edge_coupling"],
                  old["transport_node_volume"]]
    new_values = [new["poisson_edge_coupling"], new["transport_edge_coupling"],
                  new["transport_signed_average_box_node_volume"]]
    ex = list(range(3))
    axes[2].bar([value - 0.18 for value in ex], old_values, width=0.36,
                color="#9ecae9", label="M33 barycentric")
    axes[2].bar([value + 0.18 for value in ex], new_values, width=0.36,
                color="#f58518", label="M35 signed")
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_ylabel("factorial main effect (dex)")
    axes[2].set_title("Geometry-factor effects")
    axes[2].set_xticks(ex, effect_labels, rotation=22, ha="right")
    axes[2].legend(frameon=False, fontsize=8)

    for axis in axes[:2]:
        axis.set_xticks(x, labels, rotation=50, ha="right")
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
    figure.suptitle("SimpleMOS M35: signed AverageBox shared-node assembly")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=180, bbox_inches="tight")
    print(OUTPUT)


if __name__ == "__main__":
    main()
