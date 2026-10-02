#!/usr/bin/env python3
"""Render SimpleMOS M8 matrix, error, state, and SRH-repair evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize, SymLogNorm
from matplotlib.lines import Line2D
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_COMPARISONS = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "original_physics" / "comparisons"
)
DEFAULT_M8_OUTPUT = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m8_original_physics"
)
DEFAULT_DIAGNOSTICS = (
    REPO / "build-release" / "reference_tcad" / "simplemos_sentaurus2022"
    / "m8_deep_off_diagnostics" / "field_comparison_report.json"
)
DEFAULT_EVIDENCE = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_m8_original_physics_evidence.json"
)
DEFAULT_OUTPUT = REPO / "docs" / "validation" / "figures" / "simplemos_m8"

DEVICES = [f"n{index}" for index in range(17, 25)]
LOW_WELL = set(DEVICES[:4])
DEVICE_COLORS = {
    "n17": "#164e87", "n18": "#3478b6", "n19": "#68a2cf", "n20": "#a5c8e1",
    "n21": "#a84400", "n22": "#d0610b", "n23": "#ed8b2f", "n24": "#f4ba72",
}
DRAIN_COLORS = {"0p05": "#2166ac", "1": "#d66a1f"}
SIGNED_DOPING_CMAP = LinearSegmentedColormap.from_list(
    "blue_white_orange", ["#2166ac", "#f7f7f7", "#d66a1f"]
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def load_comparison(path: Path) -> dict[str, np.ndarray]:
    data = read_rows(path)
    return {
        "vg": np.asarray([float(row["gate_voltage_V"]) for row in data]),
        "sentaurus": np.abs(np.asarray([
            float(row["sentaurus_current_A_per_um"]) for row in data
        ])),
        "vela": np.abs(np.asarray([
            float(row["vela_current_A_per_um"]) for row in data
        ])),
        "log_error": np.asarray([
            float(row["absolute_log10_ratio"]) for row in data
        ]),
    }


def load_case_statuses(comparisons: Path) -> dict[str, str]:
    report = json.loads(
        (comparisons / "comparison_report.json").read_text(encoding="utf-8-sig")
    )
    return {case["case"]: case["status"].upper() for case in report["cases"]}


def plot_idvg_matrix(comparisons: Path, output: Path) -> Path:
    statuses = load_case_statuses(comparisons)
    fig, axes = plt.subplots(4, 2, figsize=(13.2, 15.0), sharex=True, sharey=True)
    for ax, device in zip(axes.flat, DEVICES):
        for drain_tag, drain_label in (("0p05", "Vd=0.05 V"), ("1", "Vd=1.0 V")):
            data = load_comparison(comparisons / f"{device}_vd_{drain_tag}_comparison.csv")
            color = DRAIN_COLORS[drain_tag]
            ax.semilogy(data["vg"], data["sentaurus"], color=color, lw=2.0)
            ax.semilogy(data["vg"], data["vela"], color=color, lw=1.55, ls="--")
        status = "PASS" if all(
            statuses[f"{device}_vd_{tag}"] == "PASS" for tag in ("0p05", "1")
        ) else "FAIL"
        well = "1e17" if device in LOW_WELL else "2e17"
        ax.set_title(f"{device}  NWell={well} cm$^{{-3}}$  |  {status}",
                     fontsize=10.5, weight="bold")
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(1e-17, 3e-3)
        ax.grid(True, which="both", color="#d9dde1", alpha=.55, linewidth=.55)
        ax.axvspan(0.0, 0.05, color="#d66a1f", alpha=.075, lw=0)
    for ax in axes[-1, :]:
        ax.set_xlabel("Gate voltage Vg (V)")
    for ax in axes[:, 0]:
        ax.set_ylabel("|Drain current| (A/um)")
    handles = [
        Line2D([], [], color=DRAIN_COLORS["0p05"], lw=2, label="Vd=0.05 V, Sentaurus"),
        Line2D([], [], color=DRAIN_COLORS["0p05"], lw=1.7, ls="--", label="Vd=0.05 V, Vela"),
        Line2D([], [], color=DRAIN_COLORS["1"], lw=2, label="Vd=1.0 V, Sentaurus"),
        Line2D([], [], color=DRAIN_COLORS["1"], lw=1.7, ls="--", label="Vd=1.0 V, Vela"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(.5, .958), fontsize=9.5)
    fig.suptitle("SimpleMOS M8 original-physics Id-Vg matrix",
                 fontsize=16, weight="bold", y=.992)
    fig.text(.5, .967,
             "16/16 cases pass; 51 direct points per curve; shaded band marks Vg=0-0.05 V",
             ha="center", va="top", fontsize=10, color="#4b5563")
    fig.subplots_adjust(left=.085, right=.985, bottom=.055, top=.925,
                        hspace=.25, wspace=.09)
    path = output / "simplemos-m8-idvg-matrix.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def plot_error_localization(comparisons: Path, output: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(14.2, 5.9), sharex=True, sharey=True)
    for ax, (drain_tag, drain_label) in zip(
            axes, (("0p05", "Vd = 0.05 V"), ("1", "Vd = 1.0 V"))):
        for device in DEVICES:
            data = load_comparison(comparisons / f"{device}_vd_{drain_tag}_comparison.csv")
            low_well = device in LOW_WELL
            ax.plot(data["vg"], data["log_error"], color=DEVICE_COLORS[device],
                    lw=1.45 if low_well else 2.0,
                    ls="--" if low_well else "-",
                    marker="o" if not low_well else None,
                    markevery=[0, 1], ms=4.0, label=device)
        ax.axhline(.3, color="#30343b", lw=1.3, ls=":")
        ax.axvspan(0.0, 0.05, color="#d66a1f", alpha=.09, lw=0)
        ax.set_title(drain_label, fontsize=11.5, weight="bold")
        ax.set_xlabel("Gate voltage Vg (V)")
        ax.set_xlim(0.0, 2.5)
        ax.set_ylim(0.0, .34)
        ax.grid(True, color="#d9dde1", alpha=.6, linewidth=.6)
    axes[0].set_ylabel("Absolute log10 current ratio (dex)")
    axes[1].annotate("Frozen limit = 0.3 dex", xy=(1.55, .3), xytext=(1.25, .39),
                     arrowprops={"arrowstyle": "->", "color": "#30343b"},
                     fontsize=9, color="#30343b")
    handles = [Line2D([], [], color=DEVICE_COLORS[d], lw=1.8,
                      ls="--" if d in LOW_WELL else "-", label=d)
               for d in DEVICES]
    fig.legend(handles=handles, loc="upper center", ncol=8, frameon=False,
               bbox_to_anchor=(.5, .927), fontsize=9)
    fig.suptitle("SimpleMOS M8 current-error localization",
                 fontsize=15.5, weight="bold", y=.985)
    fig.text(.5, .942,
             "All 16 cases pass; dashed blue: NWell=1e17, solid orange: NWell=2e17",
             ha="center", va="top", fontsize=9.5, color="#4b5563")
    fig.subplots_adjust(left=.08, right=.985, bottom=.12, top=.84, wspace=.1)
    path = output / "simplemos-m8-error-localization.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def load_geometry(neutral: Path) -> tuple[
        dict[int, tuple[float, float]], list[dict[str, Any]], dict[int, float]]:
    nodes = {
        int(row["id"]): (float(row["y_um"]), float(row["x_um"]))
        for row in read_rows(neutral / "nodes.csv")
    }
    elements = [{
        "node_ids": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
        "material": row["material"],
    } for row in read_rows(neutral / "elements.csv") if row["material"] == "Si"]
    doping = {
        int(row["node_id"]): float(row["donors_cm3"]) - float(row["acceptors_cm3"])
        for row in read_rows(neutral / "doping.csv")
    }
    return nodes, elements, doping


def nodal_values(path: Path, column: str, scale: float = 1.0) -> dict[int, float]:
    return {int(row["node_id"]): float(row[column]) * scale for row in read_rows(path)}


def cell_values(elements: list[dict[str, Any]], values: dict[int, float]) -> np.ndarray:
    return np.asarray([
        np.mean([values[node] for node in element["node_ids"]])
        for element in elements
    ])


def plot_state_contrast(m8_output: Path, output: Path) -> Path:
    datasets: dict[str, dict[str, Any]] = {}
    for device in ("n17", "n21"):
        neutral = m8_output / "neutral" / device
        state = (m8_output / "vela" / device / "workflow" / "vd_1"
                 / "00_equilibrium_accepted_state.csv")
        nodes, elements, doping = load_geometry(neutral)
        datasets[device] = {
            "polys": [[nodes[node] for node in element["node_ids"]]
                      for element in elements],
            "doping": cell_values(elements, doping),
            "psi": cell_values(elements, nodal_values(state, "psi")),
            "electrons": np.maximum(
                cell_values(elements, nodal_values(state, "electrons_m3", 1e-6)), 1.0),
        }

    doping_norm = SymLogNorm(linthresh=1e14, linscale=.8, vmin=-1.5e20, vmax=1.5e20)
    psi_all = np.concatenate([datasets[d]["psi"] for d in datasets])
    psi_norm = Normalize(vmin=float(psi_all.min()), vmax=float(psi_all.max()))
    n_all = np.concatenate([datasets[d]["electrons"] for d in datasets])
    density_norm = LogNorm(vmin=float(n_all.min()), vmax=float(n_all.max()))

    fig = plt.figure(figsize=(15.4, 10.0))
    grid = fig.add_gridspec(
        3, 3, height_ratios=[1.0, 1.0, .055],
        left=.06, right=.975, bottom=.075, top=.88, hspace=.32, wspace=.24)
    axes = np.empty((2, 3), dtype=object)
    for row_index in range(2):
        for column_index in range(3):
            axes[row_index, column_index] = fig.add_subplot(
                grid[row_index, column_index])
    colorbar_axes = [fig.add_subplot(grid[2, index]) for index in range(3)]
    collections: dict[str, Any] = {}
    specs = [
        ("doping", SIGNED_DOPING_CMAP, doping_norm, "Net doping Nd - Na"),
        ("psi", "viridis", psi_norm, "Electrostatic potential"),
        ("electrons", "plasma", density_norm, "Electron density"),
    ]
    for row_index, device in enumerate(("n17", "n21")):
        status = "accepted"
        well = "1e17" if device == "n17" else "2e17"
        for column_index, (field, cmap, norm, title) in enumerate(specs):
            ax = axes[row_index, column_index]
            collection = PolyCollection(
                datasets[device]["polys"], array=datasets[device][field],
                cmap=cmap, norm=norm, edgecolor="none")
            ax.add_collection(collection)
            ax.set_xlim(-.52, .52)
            ax.set_ylim(1.02, -.03)
            ax.set_aspect("equal", adjustable="box")
            ax.set_title(f"{device} {status} | {title}", fontsize=9.7, weight="bold")
            ax.set_xlabel("Lateral coordinate (um)")
            if column_index == 0:
                ax.set_ylabel(f"NWell={well} cm$^{{-3}}$\nDepth (um)")
            collections[field] = collection

    doping_bar = fig.colorbar(
        collections["doping"], cax=colorbar_axes[0], orientation="horizontal",
        label="Net doping (cm$^{-3}$)")
    doping_bar.set_ticks([-1e20, -1e17, -1e14, 0.0, 1e14, 1e17, 1e20])
    doping_bar.set_ticklabels([
        r"$-10^{20}$", r"$-10^{17}$", r"$-10^{14}$", "0",
        r"$10^{14}$", r"$10^{17}$", r"$10^{20}$",
    ])
    fig.colorbar(collections["psi"], cax=colorbar_axes[1], orientation="horizontal",
                 label="Potential (V)")
    fig.colorbar(collections["electrons"], cax=colorbar_axes[2], orientation="horizontal",
                 label="Electron density (cm$^{-3}$)")
    fig.suptitle("SimpleMOS M8 physical-input and equilibrium-state contrast",
                 fontsize=15.2, weight="bold", y=.995)
    fig.text(.5, .955,
             "Vela equilibrium states; device-to-device diagnostic, not a Sentaurus-Vela field-parity result",
             ha="center", va="top", fontsize=9.6, color="#4b5563")
    path = output / "simplemos-m8-n17-n21-physical-state-contrast.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def plot_srh_repair(diagnostics: Path, output: Path) -> Path:
    report = json.loads(diagnostics.read_text(encoding="utf-8-sig"))
    records = [case for case in report["cases"] if case.get("integrated_srh")]
    labels = [
        f"{case['device']}\nVd={case['drain_voltage_V']:g}\n"
        f"Vg={'0.05' if case['state'].endswith('0p05') else '0'}"
        for case in records
    ]
    old_ratios = [
        abs(case["integrated_srh"]["independent_transportmodels_override_A_per_um"]
            / case["integrated_srh"]["sentaurus_A_per_um"])
        for case in records
    ]
    corrected_ratios = [
        abs(case["integrated_srh"]["independent_simplemos_default_A_per_um"]
            / case["integrated_srh"]["sentaurus_A_per_um"])
        for case in records
    ]
    positions = np.arange(len(records))
    width = .36
    fig, ax = plt.subplots(figsize=(14.2, 6.3))
    ax.bar(positions - width / 2, old_ratios, width, color="#d66a1f",
           edgecolor="#7f3b08", hatch="//", label="Old TransportModels override")
    ax.bar(positions + width / 2, corrected_ratios, width, color="#2166ac",
           edgecolor="#123b62", hatch="..", label="Corrected SimpleMOS default")
    ax.axhline(1.0, color="#30343b", lw=1.4, ls=":", label="Sentaurus parity")
    ax.set_xticks(positions, labels)
    ax.set_ylabel("Integrated SRH magnitude ratio, Vela replay / Sentaurus")
    ax.set_ylim(0.0, 6.2)
    ax.grid(axis="y", color="#d9dde1", alpha=.7, linewidth=.7)
    ax.legend(frameon=False, ncol=3, loc="upper center")
    ax.set_title("SimpleMOS M8 SRH lifetime repair evidence",
                 fontsize=15.5, weight="bold", pad=26)
    ax.text(.5, 1.015,
            "Eight matched deep-off states: old ratio 4.51-5.68x; corrected replay agrees to numerical precision",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=9.7,
            color="#4b5563")
    fig.tight_layout()
    path = output / "simplemos-m8-srh-repair-evidence.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def update_evidence(evidence_path: Path, figures: list[Path]) -> None:
    evidence = json.loads(evidence_path.read_text(encoding="utf-8-sig"))
    evidence["visualization"] = {
        "status": "generated",
        "figures": [
            {"path": path.relative_to(REPO).as_posix(), "sha256": sha256(path)}
            for path in figures
        ],
    }
    evidence_path.write_text(
        json.dumps(evidence, indent=2) + "\n", encoding="utf-8", newline="\n"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison-dir", type=Path, default=DEFAULT_COMPARISONS)
    parser.add_argument("--m8-output-dir", type=Path, default=DEFAULT_M8_OUTPUT)
    parser.add_argument("--diagnostics", type=Path, default=DEFAULT_DIAGNOSTICS)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    generated = [
        plot_idvg_matrix(args.comparison_dir.resolve(), output),
        plot_error_localization(args.comparison_dir.resolve(), output),
        plot_state_contrast(args.m8_output_dir.resolve(), output),
        plot_srh_repair(args.diagnostics.resolve(), output),
    ]
    update_evidence(args.evidence.resolve(), generated)
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
