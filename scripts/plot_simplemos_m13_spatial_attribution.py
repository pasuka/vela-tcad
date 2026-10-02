#!/usr/bin/env python3
"""Plot SimpleMOS M13 spatial-attribution evidence."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
                  / "spatial_attribution/m13_spatial_attribution_report.json")
DEFAULT_SUMMARY = (REPO / "reference_tcad/simplemos_sentaurus2022"
                   / "spatial_attribution/m13_state_summary.csv")
DEFAULT_PROFILES = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "spatial_attribution/m13_surface_profiles.csv")
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m13"
INK, BLUE, ORANGE, GREEN, RED, PURPLE, GREY, LIGHT = (
    "#25313C", "#356C9C", "#D08237", "#4E8A67", "#B44B4B",
    "#725A9B", "#8B98A5", "#E8ECEF")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def configure() -> None:
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelcolor": INK, "axes.edgecolor": GREY,
                         "axes.titlecolor": INK, "xtick.color": INK,
                         "ytick.color": INK, "text.color": INK,
                         "figure.facecolor": "white", "axes.facecolor": "white"})


def state_label(row: dict) -> str:
    return (f"{row['device']}\nVd={float(row['drain_voltage_V']):g}\n"
            f"Vg={float(row['gate_voltage_V']):g}")


def drive_discretization(rows: list[dict[str, str]], output: Path) -> None:
    x = np.arange(len(rows))
    width = 0.25
    series = [
        ("transport_active_p95_mobility_error_dex", "transport cell-vector", BLUE),
        ("edge_projection_active_p95_mobility_error_dex", "edge projection", ORANGE),
        ("sentaurus_drive_active_p95_mobility_error_dex", "Sentaurus drive replay", GREEN),
    ]
    fig, ax = plt.subplots(figsize=(14.5, 5.8))
    for offset, (column, label, color) in zip((-1, 0, 1), series):
        ax.bar(x + offset * width, [float(row[column]) for row in rows],
               width=width, label=label, color=color)
    ax.set_xticks(x, [state_label(row) for row in rows], fontsize=7)
    ax.set_ylabel("Active-edge mobility P95 error, dex")
    ax.set_title("M13 drive discretization is regime-dependent",
                 fontsize=13, fontweight="bold")
    ax.grid(True, axis="y", color=LIGHT)
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def barrier_attribution(rows: list[dict[str, str]], output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    for row in rows:
        measured = float(row["terminal_signed_log10_error_dex"])
        predicted = float(row["barrier_predicted_log10_shift_dex"])
        weak = (row["device"] == "n21"
                and float(row["drain_voltage_V"]) == 0.05
                and float(row["gate_voltage_V"]) <= 0.8)
        ax.scatter(predicted, measured, s=70 if weak else 36,
                   color=RED if weak else BLUE, alpha=0.9,
                   edgecolor="white", linewidth=0.6)
        if weak:
            ax.annotate(f"Vg={float(row['gate_voltage_V']):g} V",
                        (predicted, measured), xytext=(5, 5),
                        textcoords="offset points", fontsize=8)
    limits = [-0.18, 0.18]
    ax.plot(limits, limits, color=GREY, linestyle="--", linewidth=1,
            label="barrier proxy = terminal error")
    ax.axhline(0, color=LIGHT, linewidth=1)
    ax.axvline(0, color=LIGHT, linewidth=1)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_xlabel(r"Barrier-predicted $\Delta\log_{10}|I_D|$, dex")
    ax.set_ylabel(r"Measured $\log_{10}(|I_{Vela}|/|I_{SDevice}|)$, dex")
    ax.set_title("Source-barrier shift has the right sign but is incomplete",
                 fontsize=12, fontweight="bold")
    ax.grid(True, color=LIGHT)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def surface_profile(rows: list[dict[str, str]], output: Path) -> None:
    selected = [row for row in rows
                if row["device"] == "n21"
                and float(row["drain_voltage_V"]) == 0.05
                and float(row["gate_voltage_V"]) in (0.0, 0.05, 0.8)]
    lateral_min = min(float(row["lateral_um"]) for row in selected)
    lateral_max = max(float(row["lateral_um"]) for row in selected)
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 8.2), sharex=True)
    colors = {0.0: BLUE, 0.05: ORANGE, 0.8: GREEN}
    for gate in colors:
        group = sorted((row for row in selected
                        if float(row["gate_voltage_V"]) == gate),
                       key=lambda row: float(row["lateral_um"]))
        y = [float(row["lateral_um"]) for row in group]
        axes[0].plot(y, [float(row["psi_difference_mV"]) for row in group],
                     color=colors[gate], label=f"Vg={gate:g} V")
        axes[1].plot(y, [float(row["phin_difference_mV"]) for row in group],
                     color=colors[gate])
        axes[2].plot(y, [float(row["electron_density_log10_ratio_dex"])
                         for row in group], color=colors[gate])
    axes[0].set_ylabel(r"$\Delta\psi$, mV")
    axes[1].set_ylabel(r"$\Delta\phi_n$, mV")
    axes[2].set_ylabel(r"$\log_{10}(n_{Vela}/n_{SDevice})$, dex")
    axes[2].set_xlabel("Si/Ox interface lateral coordinate, µm")
    axes[0].set_title("High-NWell, low-drain surface-state mismatch",
                      fontsize=13, fontweight="bold")
    axes[0].legend(frameon=False, ncol=3)
    for ax in axes:
        ax.axhline(0.0, color=GREY, linewidth=0.8)
        ax.axvspan(lateral_min, 0.0, color=ORANGE, alpha=0.04)
        ax.set_xlim(lateral_min, lateral_max)
        ax.grid(True, color=LIGHT)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def bgn_replay(rows: list[dict[str, str]], output: Path) -> None:
    selected = sorted((row for row in rows
                       if row["state"] == "n21_vd_0p05_vg_0p05"),
                      key=lambda row: float(row["lateral_um"]))
    x = [float(row["lateral_um"]) for row in selected]
    sent = [1e3 * float(row["sentaurus_bgn_eV"]) for row in selected]
    vela = [1e3 * float(row["vela_old_slotboom_total_impurity_eV"])
            for row in selected]
    residual = [float(row["old_slotboom_difference_meV"]) for row in selected]
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.8, 6.2), sharex=True,
                                  gridspec_kw={"height_ratios": [3, 1]})
    ax0.plot(x, sent, color=INK, linewidth=2.4, label="Sentaurus export")
    ax0.plot(x, vela, color=GREEN, linestyle="--", linewidth=1.7,
             label="Vela OldSlotboom replay")
    ax0.set_ylabel("Bandgap narrowing, meV")
    ax0.set_title("OldSlotboom replay matches to floating-point precision",
                  fontsize=13, fontweight="bold")
    ax0.legend(frameon=False)
    ax1.plot(x, residual, color=PURPLE, linewidth=1.5)
    ax1.axhline(0.0, color=GREY, linewidth=0.8)
    ax1.set_ylabel("Residual, meV")
    ax1.set_xlabel("Si/Ox interface lateral coordinate, µm")
    for ax in (ax0, ax1):
        ax.grid(True, color=LIGHT)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def mobility_stage(report: dict, output: Path) -> None:
    cases = report["cases"]
    stages = ["phumob", "phumob_enormal", "phumob_enormal_hfs",
              "sentaurus_drive_hfs", "sentaurus_final"]
    labels = ["PhuMob", "+ Enormal", "+ HFS (Vela drive)",
              "+ HFS (Sentaurus drive)", "Sentaurus final"]
    values = np.array([[case["edge"]["mobility_stage_active_median_cm2_V_s"][stage]
                        for case in cases] for stage in stages])
    fig, ax = plt.subplots(figsize=(14.5, 5.2))
    image = ax.imshow(values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(cases)), [state_label(case) for case in cases],
                  fontsize=7)
    ax.set_yticks(range(len(stages)), labels)
    ax.set_title("Active-edge median mobility stage chain",
                 fontsize=13, fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label(r"Electron mobility, cm$^2$/(V s)")
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--profiles", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    configure()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = json.loads(args.report.read_text(encoding="utf-8"))
    summary = read_csv(args.summary)
    profiles = read_csv(args.profiles)
    drive_discretization(summary, args.output_dir / "m13_drive_discretization_p95.png")
    barrier_attribution(summary, args.output_dir / "m13_barrier_current_attribution.png")
    surface_profile(profiles, args.output_dir / "m13_surface_profile_n21_vd0p05.png")
    bgn_replay(profiles, args.output_dir / "m13_old_slotboom_replay.png")
    mobility_stage(report, args.output_dir / "m13_mobility_stage_chain.png")
    print(json.dumps({"status": "complete", "figure_count": 5,
                      "output_dir": str(args.output_dir)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
