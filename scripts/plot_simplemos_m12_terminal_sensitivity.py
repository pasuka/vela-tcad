#!/usr/bin/env python3
"""Plot SimpleMOS M12 terminal-sensitivity evidence."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
                  / "terminal_sensitivity/m12_terminal_sensitivity_report.json")
DEFAULT_POINTS = (REPO / "reference_tcad/simplemos_sentaurus2022"
                  / "terminal_sensitivity/m12_error_spectrum_points.csv")
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m12"
INK, BLUE, ORANGE, GREEN, RED, GREY, LIGHT = (
    "#25313C", "#356C9C", "#D08237", "#4E8A67", "#B44B4B",
    "#8B98A5", "#E8ECEF")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def configure() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelcolor": INK, "axes.edgecolor": GREY,
                         "axes.titlecolor": INK, "xtick.color": INK,
                         "ytick.color": INK, "text.color": INK,
                         "figure.facecolor": "white", "axes.facecolor": "white"})


def state_label(case: dict) -> str:
    return (f"{case['device']}\nVd={case['drain_voltage_V']:g}\n"
            f"Vg={case['gate_voltage_V']:g}")


def residual_spectrum(points: list[dict[str, str]], output: Path) -> None:
    groups: dict[tuple[str, float], dict[float, list[float]]] = defaultdict(
        lambda: defaultdict(list))
    for row in points:
        key = (row["NWell_cm_minus3"], float(row["drain_voltage_V"]))
        groups[key][float(row["gate_voltage_V"])].append(
            float(row["signed_log10_ratio_dex"]))
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    styles = {
        ("1e17", 0.05): (BLUE, "-", r"$N_{Well}=10^{17}$, $V_D=0.05$ V"),
        ("1e17", 1.0): (BLUE, "--", r"$N_{Well}=10^{17}$, $V_D=1$ V"),
        ("2e17", 0.05): (ORANGE, "-", r"$N_{Well}=2\times10^{17}$, $V_D=0.05$ V"),
        ("2e17", 1.0): (ORANGE, "--", r"$N_{Well}=2\times10^{17}$, $V_D=1$ V"),
    }
    for key, by_gate in sorted(groups.items()):
        gates = sorted(by_gate)
        values = np.array([by_gate[gate] for gate in gates])
        mean = values.mean(axis=1)
        low, high = values.min(axis=1), values.max(axis=1)
        color, line, label = styles[key]
        ax.plot(gates, mean, linestyle=line, color=color, linewidth=1.8,
                label=label)
        ax.fill_between(gates, low, high, color=color, alpha=0.10)
    ax.axhline(0.0, color=GREY, linewidth=0.8)
    ax.axvspan(0.15, 0.60, color=GREEN, alpha=0.05)
    ax.axvspan(0.65, 0.95, color=ORANGE, alpha=0.04)
    ax.set_xlabel("Gate voltage, V")
    ax.set_ylabel(r"Signed $\log_{10}(|I_{Vela}|/|I_{SDevice}|)$, dex")
    ax.set_title("M12 residual spectrum: systematic positive weak-current bias",
                 fontsize=13, fontweight="bold")
    ax.grid(True, color=LIGHT)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def h1_ratio(cases: list[dict], output: Path) -> None:
    x = np.arange(len(cases))
    values = [case["h1"]["contact_to_active_reduction_ratio"] for case in cases]
    threshold = cases[0]["h1"]["threshold"]
    fig, ax = plt.subplots(figsize=(14.2, 5.5))
    ax.bar(x, values, color=[ORANGE if case["drain_voltage_V"] == 0.05
                             else BLUE for case in cases])
    ax.axhline(threshold, color=RED, linestyle="--", linewidth=1.5,
               label=f"H1 threshold = {threshold:g}")
    ax.set_yscale("log")
    ax.set_xticks(x, [state_label(case) for case in cases], fontsize=7)
    ax.set_ylabel("Max(source, drain) weighted reduction / active-edge P95 reduction")
    ax.set_title("H1 passes in all 16 states: mobility improvement is away from contact cuts",
                 fontsize=13, fontweight="bold")
    ax.grid(True, axis="y", color=LIGHT, which="both")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def terminal_shift(cases: list[dict], output: Path) -> None:
    x = np.arange(len(cases))
    values = [case["terminal_current_drive_substitution_signed_shift_dex"]
              for case in cases]
    fig, ax = plt.subplots(figsize=(14.2, 5.5))
    colors = [GREEN if value >= 0.0 else RED for value in values]
    ax.bar(x, values, color=colors)
    ax.axhline(0.0, color=GREY, linewidth=0.8)
    ax.set_xticks(x, [state_label(case) for case in cases], fontsize=7)
    ax.set_ylabel(r"$\log_{10}(|I_{alt}|/|I_{base}|)$, dex")
    ax.set_title("Frozen terminal current is insensitive to Sentaurus drive substitution",
                 fontsize=13, fontweight="bold")
    ax.grid(True, axis="y", color=LIGHT)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def condition_numbers(cases: list[dict], output: Path) -> None:
    contacts = ["source", "drain", "substrate"]
    carriers = ["kappa_electron", "kappa_hole", "kappa_total"]
    values = np.array([[case["contacts"][contact][carrier] or np.nan
                        for case in cases]
                       for contact in contacts for carrier in carriers])
    fig, ax = plt.subplots(figsize=(14.2, 5.5))
    image = ax.imshow(values, aspect="auto", cmap="viridis",
                      norm=LogNorm(vmin=1.0,
                                   vmax=max(3.0, float(np.nanmax(values)))))
    ax.set_xticks(range(len(cases)), [state_label(case) for case in cases],
                  fontsize=7)
    ax.set_yticks(range(len(contacts) * len(carriers)),
                  [f"{contact} {carrier.replace('kappa_', '')}"
                   for contact in contacts for carrier in carriers])
    ax.set_title("Carrier-specific cancellation is localized; total cut current remains conditioned",
                 fontsize=13, fontweight="bold")
    bar = fig.colorbar(image, ax=ax, pad=0.015)
    bar.set_label(r"$\kappa_I$")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = json.loads(args.report.resolve().read_text(encoding="utf-8"))
    points = read_csv(args.points.resolve())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    residual_spectrum(points, output / "simplemos_m12_residual_spectrum.png")
    h1_ratio(report["cases"], output / "simplemos_m12_h1_contact_ratio.png")
    terminal_shift(report["cases"], output / "simplemos_m12_terminal_shift.png")
    condition_numbers(report["cases"], output / "simplemos_m12_contact_conditioning.png")
    print(json.dumps({"status": "generated", "figures": 4,
                      "states": len(report["cases"])}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
