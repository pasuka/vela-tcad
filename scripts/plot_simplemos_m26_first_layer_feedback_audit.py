#!/usr/bin/env python3
"""Plot the SimpleMOS M26 first-layer self-consistent feedback audit."""

from __future__ import annotations

import csv
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "reference_tcad/simplemos_sentaurus2022/first_layer_feedback_audit"
FIGURES = REPO / "docs/validation/figures/simplemos_m26"
BLUE = "#2458a6"
ORANGE = "#d97721"
GRAY = "#697386"
GREEN = "#27864a"


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    report = json.loads(
        (DATA / "m26_first_layer_feedback_audit_report.json").read_text(
            encoding="utf-8-sig"
        )
    )
    operator = csv_rows(DATA / "m26_operator_adjoint_ledger.csv")

    fig, axes = plt.subplots(
        2, 2, figsize=(12.4, 9.0), constrained_layout=True
    )

    ax = axes[0, 0]
    directions = report["operator_directions"]
    x = np.arange(len(directions))
    quantities = (
        ("actual_self_consistent_current_change_A_per_um", "actual self-consistent", BLUE, None),
        ("direct_frozen_operator_current_change_A_per_um", "direct frozen operator", GRAY, "xx"),
        ("adjoint_relaxation_current_change_A_per_um", "adjoint relaxation", ORANGE, "//"),
        ("first_order_total_current_change_A_per_um", "first-order total", GREEN, ".."),
    )
    offsets = np.linspace(-0.27, 0.27, len(quantities))
    for offset, (field, label, color, hatch) in zip(offsets, quantities):
        values = [abs(float(row[field])) for row in directions]
        ax.bar(x + offset, values, 0.18, color=color, hatch=hatch, label=label)
    ax.set_yscale("log")
    ax.set_xticks(x, ["HFS off -> on", "HFS on -> off"])
    ax.set_ylabel("|drain-current change| (A/um)")
    ax.set_title("(a) Direct contribution vs self-consistent feedback")
    ax.grid(True, axis="y", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    ax = axes[0, 1]
    fractions = report["findings"][
        "first_layer_cut_current_fraction_by_node_no_hfs_linearization"
    ]
    node_ids = [int(node) for node in fractions]
    shares = [100.0 * abs(float(fractions[str(node)])) for node in node_ids]
    bars = ax.bar(np.arange(len(node_ids)), shares, color=BLUE)
    bars[2].set_color(ORANGE)
    bars[2].set_hatch("//")
    ax.set_xticks(np.arange(len(node_ids)), node_ids)
    ax.set_xlabel("Drain first-layer node")
    ax.set_ylabel("Terminal cut-current share (%)")
    ax.set_title("(b) First-layer terminal-current conduit")
    ax.grid(True, axis="y", alpha=0.2)
    for index, value in enumerate(shares):
        ax.text(index, value + 1.2, f"{value:.2f}%", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0.0, max(shares) * 1.16)

    ax = axes[1, 0]
    by_node: dict[int, dict[str, float]] = defaultdict(
        lambda: {"x": 0.0, "y": 0.0, "signed": 0.0, "absolute": 0.0}
    )
    first_layer = set(report["execution"]["first_layer_nodes"])
    for row in operator:
        if row["direction"] != "turn_on":
            continue
        node = int(row["node_id"])
        contribution = float(row["relaxation_current_contribution_A_per_um"])
        by_node[node]["x"] = float(row["x_m"]) * 1.0e9
        by_node[node]["y"] = float(row["y_m"]) * 1.0e9
        by_node[node]["signed"] += contribution
        by_node[node]["absolute"] += abs(contribution)
    active = [(node, values) for node, values in by_node.items() if values["absolute"] > 0.0]
    vmax = max(values["absolute"] for _, values in active)
    sizes = [18.0 + 240.0 * np.sqrt(values["absolute"] / vmax) for _, values in active]
    colors = [ORANGE if values["signed"] >= 0.0 else BLUE for _, values in active]
    ax.scatter(
        [values["x"] for _, values in active],
        [values["y"] for _, values in active],
        s=sizes,
        c=colors,
        alpha=0.58,
        edgecolors="white",
        linewidths=0.35,
        label="adjoint-weighted operator source",
    )
    layer_points = [(node, by_node[node]) for node in first_layer]
    ax.scatter(
        [values["x"] for _, values in layer_points],
        [values["y"] for _, values in layer_points],
        s=105,
        marker="*",
        facecolors="none",
        edgecolors="black",
        linewidths=1.2,
        label="drain first layer",
    )
    label_offsets = ((8, 18), (8, 3), (8, -12))
    for item, offset in zip(
        report["findings"]["top_turn_on_operator_relaxation_nodes"][:3],
        label_offsets,
    ):
        ax.annotate(
            str(item["node_id"]),
            (item["x_m"] * 1.0e9, item["y_m"] * 1.0e9),
            xytext=offset,
            textcoords="offset points",
            fontsize=7,
            arrowprops={"arrowstyle": "-", "color": "#555555", "lw": 0.6},
        )
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title("(c) HFS operator origin is upstream, not at drain cut")
    ax.grid(True, alpha=0.18)
    ax.legend(fontsize=8, loc="upper right")

    ax = axes[1, 1]
    closure = report["closure"]
    metrics = [
        ("adjoint residual", closure["maximum_adjoint_relative_residual"], 1.0e-8),
        ("FD / adjoint gradient", closure["maximum_direct_gradient_relative_disagreement"], 1.0e-4),
        ("functional / cut", closure["maximum_functional_cut_gradient_relative_disagreement"], 1.0e-8),
        ("term closure", closure["maximum_fd_term_closure_relative"], 1.0e-8),
        ("first-order model", closure["maximum_first_order_relative_error"], 1.0e-1),
        ("cut prediction", closure["maximum_cut_prediction_relative_error"], 1.0e-3),
        ("current extraction", closure["maximum_reference_current_relative_disagreement"], 1.0e-3),
    ]
    ratio = [max(value / limit, 1.0e-12) for _, value, limit in metrics]
    colors = [GREEN if value <= 1.0 else ORANGE for value in ratio]
    ax.barh(np.arange(len(metrics)), ratio, color=colors)
    ax.axvline(1.0, color="#333333", linestyle="--", linewidth=1, label="acceptance limit")
    ax.set_xscale("log")
    ax.set_yticks(np.arange(len(metrics)), [name for name, _, _ in metrics])
    ax.invert_yaxis()
    ax.set_xlabel("measured error / acceptance limit")
    ax.set_title("(d) Numerical closure and acceptance")
    ax.grid(True, axis="x", which="both", alpha=0.2)
    ax.legend(fontsize=8)

    fig.suptitle(
        "SimpleMOS M26: drain first-layer feedback audit (n23, Vd=0.05 V, deep off)",
        fontsize=13.5,
    )
    output = FIGURES / "simplemos_m26_first_layer_feedback_audit.png"
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
