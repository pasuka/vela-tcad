#!/usr/bin/env python3
"""Plot the SimpleMOS M8-A cross-device HFS confirmation evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_confirmation/comparisons/comparison_report.json"
)
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m8a"
COLORS = {"full": "#356C9C", "no_hfs": "#D08237", "phumob_only": "#C49A36"}
LABELS = {"full": "Full", "no_hfs": "No HFS", "phumob_only": "PhuMob only"}
INK = "#25313C"
GREY = "#AAB3BB"
LIGHT_GREY = "#E8ECEF"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [{
            "vg": float(row["gate_voltage_V"]),
            "sent": abs(float(row["sentaurus_current_A_per_um"])),
            "vela": abs(float(row["vela_current_A_per_um"])),
        } for row in csv.DictReader(handle)]


def case_rows(case: dict[str, Any]) -> list[dict[str, float]]:
    path = Path(case["comparison_csv"])
    if not path.is_absolute():
        path = REPO / path
    return rows(path)


def configure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelcolor": INK, "axes.edgecolor": GREY,
        "axes.titlecolor": INK, "xtick.color": INK,
        "ytick.color": INK, "text.color": INK,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def conditions(cases: list[dict[str, Any]]) -> list[tuple[str, float]]:
    return [("n17", 0.05), ("n17", 1.0), ("n21", 0.05), ("n21", 1.0)]


def find_case(cases: list[dict[str, Any]], device: str, vd: float,
              variant: str) -> dict[str, Any]:
    return next(case for case in cases if case["device"] == device
                and math.isclose(float(case["drain_voltage_V"]), vd)
                and case["variant"] == variant)


def full_idvg(cases: list[dict[str, Any]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2), sharex=True, sharey=True)
    for ax, (device, vd) in zip(axes.flat, conditions(cases)):
        case = find_case(cases, device, vd, "full")
        data = case_rows(case)
        x = [row["vg"] for row in data]
        ax.semilogy(x, [row["sent"] for row in data], color="#D08237",
                    linewidth=2.0, label="Sentaurus")
        ax.semilogy(x, [row["vela"] for row in data], color="#356C9C",
                    linewidth=1.7, linestyle="--", label="Vela")
        ax.set_title(f"{device}, Vd={vd:g} V")
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(1e-18, 5e-4)
        ax.grid(True, which="major", color=LIGHT_GREY, linewidth=0.7)
        ax.text(0.03, 0.06,
                f"max error {case['maximum_absolute_log10_ratio_above_floor']:.3f} dex",
                transform=ax.transAxes, fontsize=8)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("|Drain current|, A/um")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("SimpleMOS M8-A full-physics confirmation Id-Vg",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.008,
             "n17/n21, Vd=0.05/1 V, 51 exact gate-bias points; logarithmic current scale",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.93))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def residual(cases: list[dict[str, Any]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.0), sharex=True, sharey=True)
    for ax, (device, vd) in zip(axes.flat, conditions(cases)):
        for variant in ("full", "no_hfs", "phumob_only"):
            case = find_case(cases, device, vd, variant)
            data = case_rows(case)
            x = [row["vg"] for row in data]
            signed = [math.log10(row["vela"] / row["sent"]) for row in data]
            ax.plot(x, signed, color=COLORS[variant], linewidth=1.8,
                    linestyle="--" if variant == "phumob_only" else "-",
                    label=LABELS[variant])
        ax.axhline(0.0, color=INK, linewidth=0.9)
        ax.axvspan(0.65, 0.95, color="#C49A36", alpha=0.1)
        ax.set_title(f"{device}, Vd={vd:g} V")
        ax.set_xlim(0.0, 2.5)
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("log10(|Id,Vela| / |Id,Sentaurus|), dex")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("SimpleMOS M8-A HFS confirmation residuals",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.008,
             "Gold band marks weak inversion (Vg=0.65-0.95 V); zero denotes exact paired agreement",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.93))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def improvement(cases: list[dict[str, Any]], output: Path) -> None:
    condition_list = conditions(cases)
    labels = [f"{device}, Vd={vd:g} V" for device, vd in condition_list]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10.8, 5.2))
    width = 0.34
    for offset, variant in ((-0.5, "no_hfs"), (0.5, "phumob_only")):
        values = [find_case(cases, device, vd, variant)["effect_vs_full"][
            "maximum_log10_ratio_change"] for device, vd in condition_list]
        ax.barh(y + offset * width, values, height=width,
                color=COLORS[variant], edgecolor=INK, linewidth=0.45,
                label=LABELS[variant])
        for yy, value in zip(y + offset * width, values):
            ax.text(value + (0.001 if value >= 0 else -0.001), yy,
                    f"{value:+.4f}", va="center",
                    ha="left" if value >= 0 else "right", fontsize=8)
    ax.axvline(0.0, color=INK, linewidth=1.0)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Change in maximum absolute log-ratio vs full, dex")
    ax.set_title("HFS ablation confirmation across device and drain bias")
    ax.grid(True, axis="x", color=LIGHT_GREY, linewidth=0.7)
    ax.legend(frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, 1.12))
    fig.text(0.5, 0.005,
             "Negative values improve paired agreement; a stable source should improve consistently across controls.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.01, 0.04, 1.0, 0.94))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_summary(cases: list[dict[str, Any]], output: Path) -> None:
    fields = ["device", "drain_voltage_V", "variant", "status",
              "maximum_absolute_log10_ratio_above_floor",
              "maximum_relative_error_above_floor",
              "maximum_log10_ratio_change_vs_full",
              "deep_off_vg_0p05_log10_ratio_change_vs_full",
              "weak_inversion_max_log10_ratio_change_vs_full"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for case in cases:
            effect = case["effect_vs_full"]
            writer.writerow({
                "device": case["device"],
                "drain_voltage_V": case["drain_voltage_V"],
                "variant": case["variant"], "status": case["status"],
                "maximum_absolute_log10_ratio_above_floor": case[
                    "maximum_absolute_log10_ratio_above_floor"],
                "maximum_relative_error_above_floor": case[
                    "maximum_relative_error_above_floor"],
                "maximum_log10_ratio_change_vs_full": effect[
                    "maximum_log10_ratio_change"],
                "deep_off_vg_0p05_log10_ratio_change_vs_full": effect[
                    "deep_off_vg_0p05_log10_ratio_change"],
                "weak_inversion_max_log10_ratio_change_vs_full": effect[
                    "weak_inversion_max_log10_ratio_change"],
            })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-report", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = read_json(args.comparison_report.resolve())
    cases = report["cases"]
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    full_idvg(cases, output / "simplemos_m8a_confirmation_full_idvg.png")
    residual(cases, output / "simplemos_m8a_confirmation_residual.png")
    improvement(cases, output / "simplemos_m8a_confirmation_improvement.png")
    write_summary(cases, output / "simplemos_m8a_confirmation_summary.csv")
    print(json.dumps({"status": "generated", "figures": 3, "cases": len(cases)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
