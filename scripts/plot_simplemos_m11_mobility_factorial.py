#!/usr/bin/env python3
"""Plot the SimpleMOS M11 three-factor mobility analysis."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                / "m11_mobility_factorial")
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m11"
INK, BLUE, ORANGE, GOLD, PINK, GREY, LIGHT = (
    "#25313C", "#356C9C", "#D08237", "#B58B2A", "#A45C7A",
    "#8B98A5", "#E8ECEF")
LABELS = {
    "phumob": "PhuMob", "enormal": "Enormal", "hfs": "HFS",
    "phumob:enormal": "PhuMob×Enormal",
    "phumob:hfs": "PhuMob×HFS", "enormal:hfs": "Enormal×HFS",
    "phumob:enormal:hfs": "PhuMob×Enormal×HFS",
}


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def configure() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelcolor": INK, "axes.edgecolor": GREY,
                         "axes.titlecolor": INK, "xtick.color": INK,
                         "ytick.color": INK, "text.color": INK,
                         "figure.facecolor": "white", "axes.facecolor": "white"})


def aggregate_effects(self_report, frozen_report, output: Path) -> None:
    terms = list(LABELS)
    self_values = [self_report["aggregate_effects"]
                   ["maximum_absolute_log10_ratio_above_floor"][term]
                   for term in terms]
    frozen_values = [frozen_report["aggregate_effects"]
                     ["mobility_active_p95_error_dex"][term]
                     for term in terms]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.8))
    for ax, values, title, xlabel in (
        (axes[0], self_values, "Self-consistent paired Id–Vg error effects",
         "Effect on maximum current error, dex"),
        (axes[1], frozen_values, "Frozen-state Vela mobility-error effects",
         "Effect on active-edge mobility P95 error, dex"),
    ):
        colors = [ORANGE if value > 0 else BLUE for value in values]
        ax.barh(range(len(terms)), values, color=colors)
        ax.set_yticks(range(len(terms)), [LABELS[term] for term in terms])
        ax.axvline(0.0, color=INK, linewidth=0.9)
        ax.set_xlabel(xlabel)
        ax.set_title(title, fontsize=12, fontweight="bold")
        ax.grid(True, axis="x", color=LIGHT)
    fig.text(0.5, 0.01,
             "Negative values improve cross-solver agreement; positive values worsen it.",
             ha="center")
    fig.tight_layout(rect=(0.01, 0.04, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def condition_heatmap(report, output: Path) -> None:
    terms = list(LABELS)
    conditions = report["conditions"]
    values = np.array([[condition["effects"]
                        ["maximum_absolute_log10_ratio_above_floor"][term]
                        for term in terms] for condition in conditions])
    bound = max(abs(float(values.min())), abs(float(values.max())))
    fig, ax = plt.subplots(figsize=(12.8, 5.0))
    image = ax.imshow(values, aspect="auto", cmap="RdBu_r",
                      vmin=-bound, vmax=bound)
    ax.set_xticks(range(len(terms)), [LABELS[term] for term in terms],
                  rotation=25, ha="right")
    ax.set_yticks(range(len(conditions)), [
        f"{c['device']}, Vd={c['drain_voltage_V']:g} V" for c in conditions])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            ax.text(j, i, f"{values[i, j]:+.4f}", ha="center", va="center",
                    fontsize=8, color="white" if abs(values[i, j]) > 0.6 * bound else INK)
    ax.set_title("Condition-resolved effects on maximum Id–Vg mismatch",
                 fontsize=13, fontweight="bold")
    bar = fig.colorbar(image, ax=ax, pad=0.015)
    bar.set_label("Effect, dex")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def variant_errors(report, output: Path) -> None:
    rows = report["variant_aggregate"]
    labels = [row["variant"].replace("_", " ") for row in rows]
    maximum = [row["maximum_absolute_log10_ratio_above_floor"] for row in rows]
    median = [row["median_absolute_log10_ratio_above_floor"] for row in rows]
    x = np.arange(len(rows))
    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    width = 0.38
    ax.bar(x - width / 2, maximum, width, color=ORANGE, label="Maximum")
    ax.bar(x + width / 2, median, width, color=BLUE, label="Median")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylabel("Absolute log-current error, dex")
    ax.set_title("Aggregate paired error by factorial configuration",
                 fontsize=13, fontweight="bold")
    ax.grid(True, axis="y", color=LIGHT)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def residual_curves(root: Path, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.8, 6.2))
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, 8))
    for color, variant in zip(colors, (
            "p0_e0_h0", "p0_e0_h1", "p0_e1_h0", "p0_e1_h1",
            "p1_e0_h0", "p1_e0_h1", "p1_e1_h0", "p1_e1_h1")):
        rows = read_csv(root / "comparisons"
                        / f"n21_{variant}_vd_0p05_comparison.csv")
        x = [float(row["gate_voltage_V"]) for row in rows]
        y = [math.log10(abs(float(row["vela_current_A_per_um"]))
                        / abs(float(row["sentaurus_current_A_per_um"])))
             for row in rows]
        ax.plot(x, y, color=color, linewidth=1.5,
                label=variant.replace("_", " "))
    ax.axhline(0.0, color=INK, linewidth=0.9)
    ax.set_xlabel("Gate voltage, V")
    ax.set_ylabel("log10(|Id,Vela| / |Id,Sentaurus|), dex")
    ax.set_title("n21, Vd=0.05 V: signed residual across all eight configurations",
                 fontsize=13, fontweight="bold")
    ax.grid(True, color=LIGHT)
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    self_report = read_json(root / "self_consistent_factorial_report.json")
    frozen_report = read_json(root / "frozen_vela/frozen_factorial_report.json")
    configure()
    aggregate_effects(self_report, frozen_report,
                      output / "simplemos_m11_aggregate_effects.png")
    condition_heatmap(self_report,
                      output / "simplemos_m11_condition_effect_heatmap.png")
    variant_errors(self_report,
                   output / "simplemos_m11_variant_errors.png")
    residual_curves(root, output / "simplemos_m11_n21_vd0p05_residuals.png")
    print(json.dumps({"status": "generated", "figures": 4}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
