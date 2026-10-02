#!/usr/bin/env python3
"""Plot paired SimpleMOS M8-A model-ablation evidence."""

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
    / "m8a_model_ablation/comparisons/comparison_report.json"
)
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m8a"
VARIANT_LABELS = {
    "full": "Full physics",
    "no_srh": "No SRH",
    "plain_srh": "Plain SRH",
    "no_bgn": "No OldSlotboom",
    "no_enormal": "No Enormal",
    "no_hfs": "No high-field saturation",
    "phumob_only": "PhuMob only",
    "constant_mu": "Constant mobility",
}
BLUE = "#356C9C"
ORANGE = "#D08237"
GOLD = "#C49A36"
INK = "#25313C"
GREY = "#AAB3BB"
LIGHT_GREY = "#E8ECEF"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = []
        for row in csv.DictReader(handle):
            rows.append({
                "gate_voltage_V": float(row["gate_voltage_V"]),
                "sentaurus_current_A_per_um": abs(float(
                    row["sentaurus_current_A_per_um"])),
                "vela_current_A_per_um": abs(float(
                    row["vela_current_A_per_um"])),
            })
    return rows


def resolved(repo: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else repo / path


def load_cases(report_path: Path) -> list[dict[str, Any]]:
    report = read_json(report_path)
    cases = []
    for case in report["cases"]:
        item = dict(case)
        item["rows"] = read_rows(resolved(REPO, case["comparison_csv"]))
        cases.append(item)
    return cases


def configure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelcolor": INK,
        "axes.edgecolor": GREY,
        "axes.titlecolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def paired_idvg(cases: list[dict[str, Any]], output: Path) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.8), sharex=True, sharey=True)
    for ax, case in zip(axes.flat, cases):
        rows = case["rows"]
        x = [row["gate_voltage_V"] for row in rows]
        sent = [row["sentaurus_current_A_per_um"] for row in rows]
        vela = [row["vela_current_A_per_um"] for row in rows]
        ax.semilogy(x, sent, color=ORANGE, linewidth=2.0, label="Sentaurus")
        ax.semilogy(x, vela, color=BLUE, linewidth=1.7, linestyle="--",
                    label="Vela")
        ax.set_title(VARIANT_LABELS[case["variant"]])
        ax.grid(True, which="major", color=LIGHT_GREY, linewidth=0.7)
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(1e-18, 3e-4)
        ax.text(
            0.03, 0.06,
            f"max |log ratio| = {case['maximum_absolute_log10_ratio_above_floor']:.3f} dex",
            transform=ax.transAxes, fontsize=8, color=INK)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("|Drain current|, A/um")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               bbox_to_anchor=(0.5, 0.985))
    fig.suptitle("SimpleMOS M8-A paired Id-Vg model ablations",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.text(0.5, 0.005,
             "n23, Vd=0.05 V, 51 exact Vg points; identical declared physics per solver pair",
             ha="center", fontsize=9, color=INK)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.95))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def signed_residual(cases: list[dict[str, Any]], output: Path) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.6), sharex=True, sharey=True)
    for ax, case in zip(axes.flat, cases):
        rows = case["rows"]
        x = np.array([row["gate_voltage_V"] for row in rows])
        signed = np.array([
            math.log10(row["vela_current_A_per_um"] /
                       row["sentaurus_current_A_per_um"])
            for row in rows
        ])
        ax.axhline(0.0, color=INK, linewidth=1.0)
        ax.axvspan(0.65, 0.95, color=GOLD, alpha=0.12)
        ax.plot(x, signed, color=BLUE, linewidth=1.9)
        ax.scatter(x[[0, 1, -1]], signed[[0, 1, -1]], s=20,
                   color=ORANGE, zorder=3)
        ax.set_title(VARIANT_LABELS[case["variant"]])
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(-0.35, 0.35)
        peak = float(np.max(np.abs(signed)))
        if peak > 0.35:
            ax.text(0.97, 0.93, f"peak {peak:.2f} dex (clipped)",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=8, color=INK)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("log10(|Id,Vela| / |Id,Sentaurus|), dex")
    fig.suptitle("SimpleMOS M8-A signed current residual by model variant",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.text(0.5, 0.005,
             "Gold band: weak inversion window (Vg=0.65-0.95 V); orange markers: Vg=0, 0.05, 2.5 V",
             ha="center", fontsize=9, color=INK)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.96))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def gap_change(cases: list[dict[str, Any]], output: Path) -> None:
    variants = [case for case in cases if case["variant"] != "full"]
    labels = [VARIANT_LABELS[case["variant"]] for case in variants]
    metrics = [
        ("Maximum", "maximum_log10_ratio_change", BLUE),
        ("Weak inversion", "weak_inversion_max_log10_ratio_change", ORANGE),
        ("Deep off at Vg=0.05 V", "deep_off_vg_0p05_log10_ratio_change", GOLD),
    ]
    y = np.arange(len(variants))
    height = 0.22
    fig, (ax, zoom) = plt.subplots(
        1, 2, figsize=(14.5, 6.2), sharey=True,
        gridspec_kw={"width_ratios": [1.15, 1.0]})
    for index, (label, field, color) in enumerate(metrics):
        values = [case["effect_vs_full"][field] for case in variants]
        for target in (ax, zoom):
            target.barh(y + (index - 1) * height, values, height=height,
                        color=color, edgecolor=INK, linewidth=0.45, label=label)
    for target in (ax, zoom):
        target.axvline(0.0, color=INK, linewidth=1.0)
        target.grid(True, axis="x", color=LIGHT_GREY, linewidth=0.7)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlim(-0.04, 1.68)
    zoom.set_xlim(-0.035, 0.055)
    ax.set_xlabel("Full range, dex")
    zoom.set_xlabel("Detail near zero, dex")
    fig.suptitle("Model ablation change in paired current mismatch",
                 fontsize=14, fontweight="bold", y=0.98)
    handles, legend_labels = ax.get_legend_handles_labels()
    fig.legend(handles, legend_labels, frameon=False, ncol=3,
               loc="upper center", bbox_to_anchor=(0.5, 0.94))
    fig.text(0.5, 0.005,
             "Negative values reduce the mismatch; positive values increase it. Diagnostic, not causal, until control cases confirm.",
             ha="center", fontsize=9, color=INK)
    fig.tight_layout(rect=(0.01, 0.04, 1.0, 0.88))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def model_response(cases: list[dict[str, Any]], output: Path) -> None:
    baseline = next(case for case in cases if case["variant"] == "full")
    baseline_rows = baseline["rows"]
    variants = [case for case in cases if case["variant"] != "full"]
    fig, axes = plt.subplots(2, 4, figsize=(15.5, 7.7), sharex=True, sharey=True)
    for ax, case in zip(axes.flat, variants):
        data = case["rows"]
        x = [row["gate_voltage_V"] for row in data]
        sent_response = [math.log10(
            row["sentaurus_current_A_per_um"] /
            base["sentaurus_current_A_per_um"])
            for row, base in zip(data, baseline_rows)]
        vela_response = [math.log10(
            row["vela_current_A_per_um"] / base["vela_current_A_per_um"])
            for row, base in zip(data, baseline_rows)]
        ax.axhline(0.0, color=INK, linewidth=0.9)
        ax.plot(x, sent_response, color=ORANGE, linewidth=1.9,
                label="Sentaurus response")
        ax.plot(x, vela_response, color=BLUE, linewidth=1.7,
                linestyle="--", label="Vela response")
        ax.set_title(VARIANT_LABELS[case["variant"]])
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(-1.05, 1.05)
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
    axes.flat[-1].axis("off")
    for ax in axes[-1, :-1]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("log10(|Id,variant| / |Id,full|), dex")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2, loc="upper center",
               bbox_to_anchor=(0.5, 0.975))
    fig.suptitle("Per-solver current response to each model ablation",
                 fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.006,
             "Matching response curves indicate both solvers react similarly; their separation is the ablation-induced mismatch change.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.94))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_summary(cases: list[dict[str, Any]], output: Path) -> None:
    fields = [
        "variant", "maximum_absolute_log10_ratio_above_floor",
        "maximum_relative_error_above_floor", "endpoint_log10_ratio",
        "deep_off_vg_0_log10_ratio", "deep_off_vg_0p05_log10_ratio",
        "weak_inversion_max_log10_ratio", "weak_inversion_max_vg_V",
        "maximum_log10_ratio_change_vs_full",
        "deep_off_vg_0p05_log10_ratio_change_vs_full",
        "weak_inversion_max_log10_ratio_change_vs_full",
    ]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for case in cases:
            diagnostic = case["diagnostic_metrics"]
            effect = case["effect_vs_full"]
            writer.writerow({
                "variant": case["variant"],
                "maximum_absolute_log10_ratio_above_floor": case[
                    "maximum_absolute_log10_ratio_above_floor"],
                "maximum_relative_error_above_floor": case[
                    "maximum_relative_error_above_floor"],
                "endpoint_log10_ratio": case["endpoint_log10_ratio"],
                "deep_off_vg_0_log10_ratio": diagnostic["deep_off"][0][
                    "absolute_log10_ratio"],
                "deep_off_vg_0p05_log10_ratio": diagnostic["deep_off"][1][
                    "absolute_log10_ratio"],
                "weak_inversion_max_log10_ratio": diagnostic[
                    "weak_inversion_max"]["absolute_log10_ratio"],
                "weak_inversion_max_vg_V": diagnostic[
                    "weak_inversion_max"]["gate_voltage_V"],
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
    report = args.comparison_report.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases = load_cases(report)
    configure()
    paired_idvg(cases, output / "simplemos_m8a_paired_idvg.png")
    signed_residual(cases, output / "simplemos_m8a_signed_residual.png")
    gap_change(cases, output / "simplemos_m8a_gap_change.png")
    model_response(cases, output / "simplemos_m8a_model_response.png")
    write_summary(cases, output / "simplemos_m8a_ablation_summary.csv")
    print(json.dumps({"status": "generated", "figures": 4, "cases": len(cases)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
