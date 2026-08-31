#!/usr/bin/env python3
"""Plot the M33 region-resolved interface-assembly factorial response."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022/region_resolved_interface_assembly"
OUTPUT = (REPO / "docs/validation/figures/simplemos_m33"
          / "simplemos_m33_region_resolved_interface_assembly.png")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    response = rows(ROOT / "m33_factorial_response.csv")
    effects = rows(ROOT / "m33_factorial_main_effects.csv")
    labels = [row["variant"].replace("p", "P").replace("_t", " T").replace("_v", " V")
              for row in response]
    errors = [float(row["absolute_error_dex"]) for row in response]
    shifts = [float(row["log_current_shift_vs_baseline_dex"]) for row in response]
    colors = ["#4c78a8" if row["variant"] != "p1_t1_v1" else "#f58518"
              for row in response]

    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.4))
    axes[0].bar(range(len(response)), errors, color=colors)
    axes[0].axhline(errors[0], color="black", linestyle="--", linewidth=1,
                    label="legacy baseline")
    axes[0].set_ylabel("absolute Id error (dex)")
    axes[0].set_title("Cross-solver error")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].bar(range(len(response)), shifts, color=colors)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set_ylabel("log10(Id / baseline) (dex)")
    axes[1].set_title("Self-consistent current response")

    effect_labels = [
        "Poisson edge", "transport edge", "Si node volume"]
    effect_values = [float(row["main_effect_log_current_dex"])
                     for row in effects]
    axes[2].bar(effect_labels, effect_values,
                color=["#72b7b2", "#54a24b", "#e45756"])
    axes[2].axhline(0.0, color="black", linewidth=0.8)
    axes[2].set_ylabel("factorial main effect (dex)")
    axes[2].set_title("Region-local factor effects")
    axes[2].tick_params(axis="x", rotation=25)

    for axis in axes[:2]:
        axis.set_xticks(range(len(response)), labels, rotation=50, ha="right")
        axis.grid(axis="y", alpha=0.2)
    axes[2].grid(axis="y", alpha=0.2)
    figure.suptitle("SimpleMOS M33: shared-node Si/SiO2 region-resolved assembly")
    figure.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(OUTPUT, dpi=180, bbox_inches="tight")
    print(OUTPUT)


if __name__ == "__main__":
    main()
