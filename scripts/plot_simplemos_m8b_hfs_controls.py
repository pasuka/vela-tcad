#!/usr/bin/env python3
"""Plot SimpleMOS M8-B Sentaurus HFS-control diagnostics."""

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
    / "m8b_hfs_controls/comparisons/comparison_report.json"
)
DEFAULT_BASELINE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_confirmation"
)
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m8b"
VARIANTS = [
    "explicit_gradqf", "eparallel", "qf_at_contacts",
    "no_parallel_boundary", "refdens_efield_1e8",
]
LABELS = {
    "explicit_gradqf": "Explicit GradQF (null)",
    "eparallel": "Eparallel",
    "qf_at_contacts": "QF at contacts",
    "no_parallel_boundary": "No parallel boundary",
    "refdens_efield_1e8": "RefDens E-field 1e8",
}
COLORS = {
    "explicit_gradqf": "#8B98A5",
    "eparallel": "#356C9C",
    "qf_at_contacts": "#6F4E7C",
    "no_parallel_boundary": "#D08237",
    "refdens_efield_1e8": "#2F7D6D",
}
INK = "#25313C"
GREY = "#AAB3BB"
LIGHT_GREY = "#E8ECEF"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def resolve(path: str) -> Path:
    value = Path(path)
    return value if value.is_absolute() else REPO / value


def conditions() -> list[tuple[str, float]]:
    return [("n17", 0.05), ("n17", 1.0),
            ("n21", 0.05), ("n21", 1.0)]


def find_case(cases: list[dict[str, Any]], device: str, vd: float,
              variant: str) -> dict[str, Any]:
    return next(case for case in cases if case["device"] == device
                and math.isclose(float(case["drain_voltage_V"]), vd)
                and case["variant"] == variant)


def voltage_tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def configure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelcolor": INK, "axes.edgecolor": GREY,
        "axes.titlecolor": INK, "xtick.color": INK,
        "ytick.color": INK, "text.color": INK,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def response_curves(cases: list[dict[str, Any]], output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2), sharex=True)
    for ax, (device, vd) in zip(axes.flat, conditions()):
        for variant in VARIANTS:
            case = find_case(cases, device, vd, variant)
            data = read_rows(resolve(case["response_csv"]))
            ax.plot(
                [float(row["gate_voltage_V"]) for row in data],
                [float(row["sentaurus_log10_response"]) for row in data],
                color=COLORS[variant], linewidth=1.65,
                linestyle="--" if variant == "explicit_gradqf" else "-",
                label=LABELS[variant])
        ax.axhline(0.0, color=INK, linewidth=0.8)
        ax.set_title(f"{device}, Vd={vd:g} V")
        ax.set_xlim(0.0, 2.5)
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("log10(|Id,control| / |Id,default|), dex")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("SimpleMOS M8-B Sentaurus response to HFS controls",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.008,
             "Each curve changes one Sentaurus control; 51 exact Vg points, no interpolation",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.92))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def error_change(cases: list[dict[str, Any]], output: Path) -> None:
    labels = [f"{device}, Vd={vd:g} V" for device, vd in conditions()]
    y = np.arange(len(labels))
    width = 0.15
    fig, ax = plt.subplots(figsize=(11.4, 6.3))
    offsets = np.linspace(-2, 2, len(VARIANTS)) * width
    for offset, variant in zip(offsets, VARIANTS):
        values = [find_case(cases, device, vd, variant)["effect_vs_default"][
            "maximum_log10_error_change"] for device, vd in conditions()]
        ax.barh(y + offset, values, height=width * 0.9,
                color=COLORS[variant], label=LABELS[variant])
    ax.axvline(0.0, color=INK, linewidth=1.0)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("Change in maximum absolute Vela/Sentaurus error, dex")
    ax.set_title("HFS control effect on the fixed-Vela comparison")
    ax.grid(True, axis="x", color=LIGHT_GREY, linewidth=0.7)
    ax.legend(frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, 1.18))
    fig.text(0.5, 0.01,
             "Negative values improve agreement; consistency across all four conditions is the attribution guard.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.01, 0.045, 1.0, 0.9))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def residual_curves(cases: list[dict[str, Any]], baseline: Path,
                    output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2), sharex=True)
    for ax, (device, vd) in zip(axes.flat, conditions()):
        tag = voltage_tag(vd)
        base = read_rows(baseline / "comparisons"
                         / f"{device}_full_vd_{tag}_comparison.csv")
        ax.plot([float(row["gate_voltage_V"]) for row in base],
                [float(row["absolute_log10_ratio"]) for row in base],
                color=INK, linewidth=2.1, label="Default HFS")
        for variant in VARIANTS[1:]:
            case = find_case(cases, device, vd, variant)
            data = read_rows(resolve(case["response_csv"]))
            ax.plot(
                [float(row["gate_voltage_V"]) for row in data],
                [float(row["controlled_absolute_log10_error"]) for row in data],
                color=COLORS[variant], linewidth=1.45,
                label=LABELS[variant])
        ax.set_title(f"{device}, Vd={vd:g} V")
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(bottom=0.0)
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("|log10(|Id,Vela| / |Id,Sentaurus|)|, dex")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle("SimpleMOS M8-B pointwise residual under HFS controls",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.008,
             "Vela is fixed to the M8-A full-physics candidate; only Sentaurus controls vary",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.92))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def idvg_best_control(report: dict[str, Any], baseline: Path,
                      output: Path) -> str:
    cases = report["cases"]
    candidates = [summary for summary in report["variant_summaries"]
                  if summary["variant"] != "explicit_gradqf"]
    best = min(candidates, key=lambda summary: sum(
        summary["maximum_error_changes_dex"]) / summary["condition_count"])
    variant = str(best["variant"])
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2), sharex=True, sharey=True)
    for ax, (device, vd) in zip(axes.flat, conditions()):
        tag = voltage_tag(vd)
        base = read_rows(baseline / "comparisons"
                         / f"{device}_full_vd_{tag}_comparison.csv")
        control_case = find_case(cases, device, vd, variant)
        control = read_rows(resolve(control_case["response_csv"]))
        x = [float(row["gate_voltage_V"]) for row in base]
        ax.semilogy(x, [abs(float(row["sentaurus_current_A_per_um"]))
                        for row in base], color="#8B98A5", linewidth=1.6,
                    label="Sentaurus default")
        ax.semilogy(x, [abs(float(row[
            "controlled_sentaurus_current_A_per_um"])) for row in control],
                    color=COLORS[variant], linewidth=2.0,
                    label=f"Sentaurus {LABELS[variant]}")
        ax.semilogy(x, [abs(float(row["vela_current_A_per_um"]))
                        for row in control], color="#356C9C", linewidth=1.6,
                    linestyle="--", label="Vela fixed")
        ax.set_title(f"{device}, Vd={vd:g} V")
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(1e-18, 5e-4)
        ax.grid(True, which="major", color=LIGHT_GREY, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("|Drain current|, A/um")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
               bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(f"SimpleMOS M8-B Id-Vg: best mean control ({LABELS[variant]})",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.008,
             "Best is selected by mean change in maximum log error across four conditions",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.92))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return variant


def write_summary(report: dict[str, Any], output: Path) -> None:
    fields = ["device", "drain_voltage_V", "variant", "dimension",
              "maximum_absolute_log10_ratio_above_floor",
              "median_absolute_log10_ratio_above_floor",
              "maximum_error_change_dex", "median_error_change_dex",
              "maximum_sentaurus_response_dex",
              "median_sentaurus_response_dex"]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for case in report["cases"]:
            writer.writerow({
                "device": case["device"],
                "drain_voltage_V": case["drain_voltage_V"],
                "variant": case["variant"], "dimension": case["dimension"],
                "maximum_absolute_log10_ratio_above_floor": case[
                    "maximum_absolute_log10_ratio_above_floor"],
                "median_absolute_log10_ratio_above_floor": case[
                    "median_absolute_log10_ratio_above_floor"],
                "maximum_error_change_dex": case["effect_vs_default"][
                    "maximum_log10_error_change"],
                "median_error_change_dex": case["effect_vs_default"][
                    "median_log10_error_change"],
                "maximum_sentaurus_response_dex": case[
                    "maximum_sentaurus_response_dex"],
                "median_sentaurus_response_dex": case[
                    "median_sentaurus_response_dex"],
            })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-report", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = read_json(args.comparison_report.resolve())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    response_curves(report["cases"], output / "simplemos_m8b_sentaurus_response.png")
    error_change(report["cases"], output / "simplemos_m8b_error_change.png")
    residual_curves(report["cases"], args.baseline_dir.resolve(),
                    output / "simplemos_m8b_residual.png")
    best = idvg_best_control(report, args.baseline_dir.resolve(),
                            output / "simplemos_m8b_best_control_idvg.png")
    write_summary(report, output / "simplemos_m8b_hfs_controls_summary.csv")
    print(json.dumps({"status": "generated", "figures": 4,
                      "cases": len(report["cases"]), "best_control": best}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
