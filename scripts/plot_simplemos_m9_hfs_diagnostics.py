#!/usr/bin/env python3
"""Plot SimpleMOS M9 HFS reference-density and Vela edge diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m9_hfs_diagnostics/comparisons/comparison_report.json"
)
DEFAULT_EDGE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m9_hfs_diagnostics/vela_edge_probes/edge_hfs_report.json"
)
DEFAULT_OUTPUT = REPO / "docs/validation/figures/simplemos_m9"
INK = "#25313C"
BLUE = "#356C9C"
GOLD = "#B58B2A"
ORANGE = "#D08237"
PINK = "#A45C7A"
GREY = "#8B98A5"
LIGHT_GREY = "#E8ECEF"
CONDITION_COLORS = {
    ("n17", 0.05): BLUE,
    ("n17", 1.0): GOLD,
    ("n21", 0.05): ORANGE,
    ("n21", 1.0): PINK,
}
CONDITION_MARKERS = {0.05: "o", 1.0: "s"}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def resolve(path: str) -> Path:
    value = Path(path)
    return value if value.is_absolute() else REPO / value


def conditions() -> list[tuple[str, float]]:
    return [("n17", 0.05), ("n17", 1.0),
            ("n21", 0.05), ("n21", 1.0)]


def condition_label(device: str, vd: float) -> str:
    return f"{device}, Vd={vd:g} V"


def configure() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9,
        "axes.labelcolor": INK, "axes.edgecolor": GREY,
        "axes.titlecolor": INK, "xtick.color": INK,
        "ytick.color": INK, "text.color": INK,
        "figure.facecolor": "white", "axes.facecolor": "white",
    })


def positive_scan(report: dict[str, Any]) -> list[dict[str, Any]]:
    return sorted(
        [item for item in report["variant_summaries"]
         if float(item["refdens_cm3"]) > 0.0],
        key=lambda item: float(item["refdens_cm3"]),
    )


def find_case(report: dict[str, Any], device: str, vd: float,
              variant: str) -> dict[str, Any]:
    return next(item for item in report["cases"]
                if item["device"] == device
                and math.isclose(float(item["drain_voltage_V"]), vd)
                and item["variant"] == variant)


def scan_metric(report: dict[str, Any], output: Path, metric: str,
                ylabel: str, title: str, note: str) -> None:
    scan = positive_scan(report)
    fig, ax = plt.subplots(figsize=(10.8, 6.2))
    for device, vd in conditions():
        values = []
        for summary in scan:
            case = find_case(report, device, vd, summary["variant"])
            if metric == "maximum_error_change_dex":
                values.append(float(case["effect_vs_default"][
                    "maximum_log10_error_change"]))
            else:
                values.append(float(case[metric]))
        ax.plot(
            [float(item["refdens_cm3"]) for item in scan], values,
            color=CONDITION_COLORS[(device, vd)],
            marker=CONDITION_MARKERS[vd], linewidth=1.8, markersize=5.5,
            linestyle="-" if device == "n17" else "--",
            label=condition_label(device, vd),
        )
    ax.set_xscale("log")
    ax.axhline(0.0, color=INK, linewidth=0.9)
    ax.set_xlabel("Sentaurus HFS reference density, cm$^{-3}$")
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.grid(True, which="both", color=LIGHT_GREY, linewidth=0.7)
    ax.legend(frameon=False, ncol=2)
    fig.text(0.5, 0.01, note, ha="center", fontsize=9)
    fig.tight_layout(rect=(0.01, 0.045, 1.0, 1.0))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def pointwise_response(report: dict[str, Any], output: Path) -> str:
    positive = [item for item in report["variant_summaries"]
                if float(item["refdens_cm3"]) > 0.0]
    selected = max(positive, key=lambda item: max(
        float(case["maximum_sentaurus_response_dex"])
        for case in report["cases"] if case["variant"] == item["variant"]))
    variant = str(selected["variant"])
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.2), sharex=True)
    for ax, (device, vd) in zip(axes.flat, conditions()):
        case = find_case(report, device, vd, variant)
        rows = read_rows(resolve(case["response_csv"]))
        x = [float(row["gate_voltage_V"]) for row in rows]
        ax.plot(x, [float(row["sentaurus_log10_response"]) for row in rows],
                color=BLUE, linewidth=1.9, label="Sentaurus response")
        ax.plot(x, [float(row["controlled_absolute_log10_error"])
                    for row in rows], color=ORANGE, linewidth=1.6,
                linestyle="--", label="Vela/Sentaurus abs. residual")
        ax.axhline(0.0, color=INK, linewidth=0.8)
        ax.set_title(condition_label(device, vd))
        ax.set_xlim(0.0, 2.5)
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage, V")
    for ax in axes[:, 0]:
        ax.set_ylabel("Log-current difference, dex")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2,
               loc="upper center", bbox_to_anchor=(0.5, 0.965))
    fig.suptitle(
        f"Pointwise HFS RefDens response ({float(selected['refdens_cm3']):.0e} cm$^{{-3}}$)",
        fontsize=14, fontweight="bold", y=1.0)
    fig.text(0.5, 0.008,
             "51 direct Vg points per curve; response is relative to explicit GradQuasiFermi.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.02, 0.035, 1.0, 0.92))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return variant


def edge_limiter(edge_report: dict[str, Any], output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 5.4), sharex=True, sharey=True)
    for ax, carrier in zip(axes, ("electron", "hole")):
        for device, vd in conditions():
            rows = sorted(
                [item for item in edge_report["cases"]
                 if item["device"] == device
                 and math.isclose(float(item["drain_voltage_V"]), vd)],
                key=lambda item: float(item["gate_voltage_V"]),
            )
            ax.plot(
                [float(item["gate_voltage_V"]) for item in rows],
                [float(item[f"minimum_{carrier}_mobility_limiter"])
                 for item in rows],
                color=CONDITION_COLORS[(device, vd)],
                marker=CONDITION_MARKERS[vd], linewidth=1.7,
                linestyle="-" if device == "n17" else "--",
                label=condition_label(device, vd),
            )
        ax.set_title(f"{carrier.capitalize()} effective limiter")
        ax.set_xlabel("Gate voltage, V")
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, color=LIGHT_GREY, linewidth=0.7)
    axes[0].set_ylabel("Minimum edge mobility / low-field mobility")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=2,
               loc="upper center", bbox_to_anchor=(0.5, 0.98))
    fig.suptitle("Vela per-edge HighFieldSaturation limiting strength",
                 fontsize=14, fontweight="bold", y=1.06)
    fig.text(0.5, 0.01,
             "Minimum over transport-capable edges; lower values mean stronger velocity saturation.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0.01, 0.045, 1.0, 0.91))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--edge-report", type=Path, default=DEFAULT_EDGE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = read_json(args.comparison_report.resolve())
    edge_report = read_json(args.edge_report.resolve())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    configure()
    scan_metric(
        report, output / "simplemos_m9_refdens_error_change.png",
        "maximum_error_change_dex",
        "Change in maximum absolute Vela/Sentaurus error, dex",
        "HFS reference-density scan: comparison error change",
        "Negative values improve agreement; all four conditions must move consistently for attribution.",
    )
    scan_metric(
        report, output / "simplemos_m9_refdens_sentaurus_response.png",
        "maximum_sentaurus_response_dex",
        "Maximum |log10(Id,RefDens / Id,default)|, dex",
        "HFS reference-density scan: Sentaurus response",
        "Maximum response over 51 exact Vg points; explicit GradQuasiFermi is the zero reference.",
    )
    selected = pointwise_response(
        report, output / "simplemos_m9_refdens_pointwise_response.png")
    edge_limiter(edge_report, output / "simplemos_m9_vela_edge_limiter.png")
    print(json.dumps({"status": "generated", "figures": 4,
                      "pointwise_variant": selected}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
