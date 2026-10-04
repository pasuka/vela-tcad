#!/usr/bin/env python3
"""Plot SimpleMOS M10 fixed-state replay diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                  / "m10_fixed_state_replay/fixed_state_replay_report.json")
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m10"
INK, BLUE, ORANGE, GOLD, GREY, LIGHT = (
    "#25313C", "#356C9C", "#D08237", "#B58B2A", "#8B98A5", "#E8ECEF")


def configure() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelcolor": INK, "axes.edgecolor": GREY,
                         "axes.titlecolor": INK, "xtick.color": INK,
                         "ytick.color": INK, "text.color": INK,
                         "figure.facecolor": "white", "axes.facecolor": "white"})


def label(case: dict) -> str:
    return (f"{case['device']}\nVd={case['drain_voltage_V']:g}\n"
            f"Vg={case['gate_voltage_V']:g}")


def mobility_p95(cases: list[dict], output: Path) -> None:
    x = np.arange(len(cases))
    before = [c["electron_mobility"]["vela_drive"]
              ["active_edges_abs_error_dex"]["p95"] for c in cases]
    after = [c["electron_mobility"]["sentaurus_drive"]
             ["active_edges_abs_error_dex"]["p95"] for c in cases]
    fig, ax = plt.subplots(figsize=(14.2, 6.2))
    ax.plot(x, before, "o-", color=ORANGE, linewidth=1.7,
            label="Vela drive → Vela HFS")
    ax.plot(x, after, "o-", color=BLUE, linewidth=1.7,
            label="Sentaurus drive → Vela HFS")
    ax.set_xticks(x, [label(c) for c in cases], fontsize=7)
    ax.set_ylabel("Active-edge electron mobility P95 error, dex")
    ax.set_title("Frozen-state replay isolates the GradQuasiFermi drive contribution",
                 fontsize=13, fontweight="bold")
    ax.grid(True, axis="y", color=LIGHT)
    ax.legend(frameon=False)
    fig.text(0.5, 0.01, "All 16 states improve at P95; no bias interpolation.",
             ha="center")
    fig.tight_layout(rect=(0.01, 0.04, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def terminal_heatmap(cases: list[dict], output: Path) -> None:
    stages = ["vela_drive_vela_hfs", "sentaurus_drive_vela_hfs",
              "sentaurus_final_mobility"]
    names = ["Vela drive", "Sentaurus drive", "Sentaurus final mobility"]
    values = np.array([[c["terminal_current_replay"][stage]
                        ["absolute_log10_error_dex"] for c in cases]
                       for stage in stages])
    fig, ax = plt.subplots(figsize=(14.2, 4.4))
    image = ax.imshow(values, aspect="auto", cmap="magma", vmin=0,
                      vmax=min(3.0, float(values.max())))
    ax.set_xticks(range(len(cases)), [label(c) for c in cases], fontsize=7)
    ax.set_yticks(range(len(stages)), names)
    ax.set_title("Terminal-current replay: deep-off states are ill-conditioned",
                 fontsize=13, fontweight="bold")
    bar = fig.colorbar(image, ax=ax, pad=0.015)
    bar.set_label("Absolute log-current error, dex")
    fig.text(0.5, 0.01,
             "Vg ≥ 0.8 V is near exact; sub-fA off-current amplifies state-mapping and cancellation errors.",
             ha="center")
    fig.tight_layout(rect=(0.01, 0.05, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def mapping_sensitivity(cases: list[dict], output: Path) -> None:
    primary = [c["electron_mobility"]["sentaurus_drive"]
               ["active_edges_abs_error_dex"]["p95"] for c in cases]
    sensitivity = [c["electron_mobility"]["sentaurus_drive_sensitivity"]
                   ["active_edges_abs_error_dex"]["p95"] for c in cases]
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    ax.scatter(primary, sensitivity, c=[c["gate_voltage_V"] for c in cases],
               cmap="viridis", s=48, edgecolor="white", linewidth=0.5)
    lo, hi = min(primary + sensitivity), max(primary + sensitivity)
    ax.plot([lo, hi], [lo, hi], color=GREY, linestyle="--")
    ax.set_xlabel("Mean endpoint magnitude mapping P95, dex")
    ax.set_ylabel("Magnitude of mean vector mapping P95, dex")
    ax.set_title("Node-to-edge mapping sensitivity", fontsize=13,
                 fontweight="bold")
    ax.grid(True, color=LIGHT)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = json.loads(args.report.resolve().read_text(encoding="utf-8"))
    cases = report["cases"]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    mobility_p95(cases, output / "simplemos_m10_mobility_p95_replay.png")
    terminal_heatmap(cases, output / "simplemos_m10_terminal_error_heatmap.png")
    mapping_sensitivity(cases, output / "simplemos_m10_mapping_sensitivity.png")
    print(json.dumps({"status": "generated", "figures": 3, "states": len(cases)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
