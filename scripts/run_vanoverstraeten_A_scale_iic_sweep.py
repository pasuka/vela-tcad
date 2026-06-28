#!/usr/bin/env python3
"""Post-process a fixed PN2D solution with VanOverstraeten A-scale IIC sweeps."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
E_CHARGE_C = 1.602176634e-19
DEPTH_CM = 1.0e-4

DEFAULT_SCALES = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 24.0, 32.0]
WINDOWS = {
    "full": None,
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}

# Built-in 300 K VanOverstraeten/de Man defaults, converted from the C++ SI values.
ELECTRON_A_LOW_CM_INV = 7.03e5
ELECTRON_A_HIGH_CM_INV = 7.03e5
ELECTRON_B_LOW_V_CM = 1.231e6
ELECTRON_B_HIGH_V_CM = 1.231e6
HOLE_A_LOW_CM_INV = 1.582e6
HOLE_A_HIGH_CM_INV = 6.71e5
HOLE_B_LOW_V_CM = 2.036e6
HOLE_B_HIGH_V_CM = 1.693e6
SWITCH_FIELD_V_CM = 4.0e5


def clean_row(row: dict[str, Any]) -> dict[str, str]:
    return {
        str(key).strip(): str(value).strip()
        for key, value in row.items()
        if key is not None
    }


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [clean_row(row) for row in csv.DictReader(handle, skipinitialspace=True)]


def write_rows(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite_float(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def fmt(value: Any) -> str:
    number = finite_float(value)
    if number is None:
        return ""
    return f"{number:.12g}"


def fmt_scale(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def parse_scales(text: str) -> list[float]:
    scales = [float(item.strip()) for item in text.split(",") if item.strip()]
    if not scales:
        raise ValueError("--scales must contain at least one value")
    return scales


class NodeRecord:
    def __init__(self, node_id: int, x_um: float, y_um: float) -> None:
        self.node_id = node_id
        self.x_um = x_um
        self.y_um = y_um
        self.values: dict[str, tuple[float | None, float | None]] = {}

    def set_value(self, quantity: str, sentaurus: float | None, vela: float | None) -> None:
        self.values[quantity] = (sentaurus, vela)

    def sentaurus(self, quantity: str) -> float:
        value = self.values.get(quantity, (0.0, 0.0))[0]
        return value if value is not None else 0.0

    def vela(self, quantity: str) -> float:
        value = self.values.get(quantity, (0.0, 0.0))[1]
        return value if value is not None else 0.0


def load_compare(path: Path) -> dict[float, dict[int, NodeRecord]]:
    by_bias: dict[float, dict[int, NodeRecord]] = defaultdict(dict)
    for row in read_rows(path):
        bias = finite_float(row.get("bias_V"))
        node_id_f = finite_float(row.get("node_id"))
        x_um = finite_float(row.get("x_um"))
        y_um = finite_float(row.get("y_um"))
        quantity = row.get("quantity", "")
        if bias is None or node_id_f is None or x_um is None or y_um is None or not quantity:
            continue
        node_id = int(node_id_f)
        nodes = by_bias[bias]
        if node_id not in nodes:
            nodes[node_id] = NodeRecord(node_id, x_um, y_um)
        nodes[node_id].set_value(
            quantity,
            finite_float(row.get("sentaurus_value")),
            finite_float(row.get("vela_value_scaled_to_sentaurus_units")),
        )
    return dict(by_bias)


def load_elements(path: Path) -> list[tuple[int, int, int]]:
    elements: list[tuple[int, int, int]] = []
    for row in read_rows(path):
        try:
            elements.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
        except (KeyError, ValueError):
            continue
    return elements


def neighbor_map(elements: list[tuple[int, int, int]]) -> dict[int, set[int]]:
    neighbors: dict[int, set[int]] = defaultdict(set)
    for a, b, c in elements:
        for i, j in ((a, b), (a, c), (b, a), (b, c), (c, a), (c, b)):
            neighbors[i].add(j)
    return neighbors


def recover_gradients_v_per_cm(
    nodes: dict[int, NodeRecord],
    elements: list[tuple[int, int, int]],
    quantity: str,
) -> dict[int, float]:
    neighbors = neighbor_map(elements)
    out: dict[int, float] = {}
    for node_id, node in nodes.items():
        center = node.vela(quantity)
        ata00 = ata01 = ata11 = atb0 = atb1 = 0.0
        for other_id in neighbors.get(node_id, set()):
            other = nodes.get(other_id)
            if other is None:
                continue
            dx = other.x_um - node.x_um
            dy = other.y_um - node.y_um
            dist = math.hypot(dx, dy)
            if dist <= 0.0:
                continue
            weight = 1.0 / dist
            dv = other.vela(quantity) - center
            ata00 += weight * dx * dx
            ata01 += weight * dx * dy
            ata11 += weight * dy * dy
            atb0 += weight * dx * dv
            atb1 += weight * dy * dv
        det = ata00 * ata11 - ata01 * ata01
        if abs(det) <= 1.0e-30:
            out[node_id] = 0.0
            continue
        grad_x_v_um = (atb0 * ata11 - atb1 * ata01) / det
        grad_y_v_um = (ata00 * atb1 - ata01 * atb0) / det
        out[node_id] = math.hypot(grad_x_v_um, grad_y_v_um) * 1.0e4
    return out


def alpha_vanoverstraeten_cm_inv(field_v_cm: float, carrier: str, scale: float) -> float:
    field = abs(field_v_cm)
    if field <= 0.0:
        return 0.0
    low = field < SWITCH_FIELD_V_CM
    if carrier == "electron":
        a = ELECTRON_A_LOW_CM_INV if low else ELECTRON_A_HIGH_CM_INV
        b = ELECTRON_B_LOW_V_CM if low else ELECTRON_B_HIGH_V_CM
    else:
        a = HOLE_A_LOW_CM_INV if low else HOLE_A_HIGH_CM_INV
        b = HOLE_B_LOW_V_CM if low else HOLE_B_HIGH_V_CM
    exponent = max(min(-b / field, 0.0), -700.0)
    return scale * a * math.exp(exponent)


def triangle_area_um2(nodes: dict[int, NodeRecord], element: tuple[int, int, int]) -> float:
    n0, n1, n2 = (nodes[element[0]], nodes[element[1]], nodes[element[2]])
    return abs(
        (n1.x_um - n0.x_um) * (n2.y_um - n0.y_um)
        - (n2.x_um - n0.x_um) * (n1.y_um - n0.y_um)
    ) * 0.5


def centroid_x_um(nodes: dict[int, NodeRecord], element: tuple[int, int, int]) -> float:
    return sum(nodes[node_id].x_um for node_id in element) / 3.0


def in_window(x_um: float, window: tuple[float, float] | None) -> bool:
    return window is None or window[0] <= x_um <= window[1]


def integrate_by_window(
    nodes: dict[int, NodeRecord],
    elements: list[tuple[int, int, int]],
    values: dict[int, float],
) -> dict[str, float]:
    totals = {name: 0.0 for name in WINDOWS}
    for element in elements:
        if any(node_id not in nodes or node_id not in values for node_id in element):
            continue
        area_cm2 = triangle_area_um2(nodes, element) * 1.0e-8
        volume_cm3 = area_cm2 * DEPTH_CM
        avg_value = sum(values[node_id] for node_id in element) / 3.0
        x_um = centroid_x_um(nodes, element)
        for name, window in WINDOWS.items():
            if in_window(x_um, window):
                totals[name] += avg_value * volume_cm3
    return totals


def ratio_or_blank(value: float, reference: float, floor: float) -> float | None:
    if abs(reference) < floor:
        return None
    return value / reference


def region_error_metric(ratios: list[float | None]) -> float | None:
    scores = [abs(math.log10(abs(ratio))) for ratio in ratios if ratio is not None and ratio > 0.0]
    if not scores:
        return None
    return statistics.median(scores)


def compute_rows(
    data: dict[float, dict[int, NodeRecord]],
    elements: list[tuple[int, int, int]],
    scales: list[float],
    sentaurus_qg_floor: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias in sorted(data):
        nodes = data[bias]
        sentaurus_g = {
            node_id: max(node.sentaurus("avalanche"), 0.0)
            for node_id, node in nodes.items()
        }
        default_self_g = {
            node_id: max(node.vela("avalanche"), 0.0)
            for node_id, node in nodes.items()
        }
        sentaurus_qg = {
            node_id: sentaurus_g[node_id] * E_CHARGE_C
            for node_id in nodes
        }
        default_self_qg = {
            node_id: default_self_g[node_id] * E_CHARGE_C
            for node_id in nodes
        }
        sentaurus_integrals = integrate_by_window(nodes, elements, sentaurus_qg)
        default_self_integrals = integrate_by_window(nodes, elements, default_self_qg)
        max_g_sentaurus = max(sentaurus_g.values(), default=0.0)
        max_g_default_self = max(default_self_g.values(), default=0.0)

        for scale in scales:
            self_integrals = {
                name: scale * value
                for name, value in default_self_integrals.items()
            }
            ratios = {
                name: ratio_or_blank(self_integrals[name], sentaurus_integrals[name], sentaurus_qg_floor)
                for name in WINDOWS
            }
            row = {
                "bias_V": fmt(bias),
                "A_scale": fmt_scale(scale),
                "qG_full_self_A_per_um": fmt(self_integrals["full"]),
                "qG_full_sentaurus_A_per_um": fmt(sentaurus_integrals["full"]),
                "qG_full_ratio": fmt(ratios["full"]),
                "qG_junction_ratio": fmt(ratios["junction"]),
                "qG_center_ratio": fmt(ratios["center"]),
                "qG_left_shoulder_ratio": fmt(ratios["left_shoulder"]),
                "qG_right_shoulder_ratio": fmt(ratios["right_shoulder"]),
                "max_G_self": fmt(scale * max_g_default_self),
                "max_G_sentaurus": fmt(max_g_sentaurus),
                "best_region_error_metric": fmt(region_error_metric(list(ratios.values()))),
                "signal_status": (
                    "low_signal"
                    if abs(sentaurus_integrals["full"]) < sentaurus_qg_floor
                    else "ok"
                ),
            }
            rows.append(row)
    return rows


def best_scale_for_bias(rows: list[dict[str, Any]], bias: float) -> tuple[str, float | None]:
    candidates: list[tuple[float, str]] = []
    bias_text = fmt(bias)
    for row in rows:
        if row.get("bias_V") != bias_text:
            continue
        score = finite_float(row.get("best_region_error_metric"))
        if score is not None:
            candidates.append((score, row["A_scale"]))
    if not candidates:
        return "insufficient_signal", None
    score, scale = min(candidates, key=lambda item: item[0])
    return scale, score


def rank_scales(rows: list[dict[str, Any]], target_biases: list[float]) -> list[tuple[float, str]]:
    scores_by_scale: dict[str, list[float]] = defaultdict(list)
    bias_texts = {fmt(bias) for bias in target_biases}
    for row in rows:
        if row.get("bias_V") not in bias_texts:
            continue
        score = finite_float(row.get("best_region_error_metric"))
        if score is not None:
            scores_by_scale[row["A_scale"]].append(score)
    ranked = [
        (statistics.median(scores), scale)
        for scale, scores in scores_by_scale.items()
        if scores
    ]
    return sorted(ranked, key=lambda item: item[0])


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    target_biases = [-5.0, -10.0, -16.0, -18.0, -20.0]
    best = {bias: best_scale_for_bias(rows, bias) for bias in target_biases}
    ranked = rank_scales(rows, target_biases)
    recommendations = [scale for _, scale in ranked[:3]]
    best_scales = {
        scale for scale, score in best.values()
        if scale != "insufficient_signal" and score is not None
    }
    single_scale = len(best_scales) == 1 and len(best_scales) > 0
    low_signal_biases = sorted({
        row["bias_V"] for row in rows
        if row.get("signal_status") == "low_signal"
    }, key=lambda text: float(text))

    lines = [
        "# VanOverstraeten A-Scale IIC Sweep",
        "",
        "- This is post-processing/IIC only: avalanche generation is not fed back into the continuity equations.",
        "- A_scale=1 integrates the existing Vela default AvalancheGeneration node output from the aligned CSV.",
        "- Other A_scale values linearly rescale that fixed default source, preserving the default B values, default low/high branch selection, default switchField=4.0e5 V/cm, and current source mapping.",
        "- qG columns are q integral AvalancheGeneration dV over a 1 um depth and are reported in A/um.",
        "- Rows marked low-signal have Sentaurus full-domain qG below the configured floor, so ratios are blanked for low-signal regions.",
        "",
        "## Best A_scale By Bias",
        "",
        "| bias | best A_scale | median abs log10 region error |",
        "|---:|---:|---:|",
    ]
    for bias in target_biases:
        scale, score = best[bias]
        lines.append(f"| {bias:g}V | {scale} | {fmt(score)} |")

    lines.extend([
        "",
        "## Answers",
        "",
        f"1. -5V q integral Gava is closest at A_scale={best[-5.0][0]}.",
        f"2. -10V q integral Gava is closest at A_scale={best[-10.0][0]}.",
        (
            "3. High-bias best A_scale values: "
            f"-16V={best[-16.0][0]}, -18V={best[-18.0][0]}, -20V={best[-20.0][0]}."
        ),
        (
            "4. A single A_scale covers -5V to -20V: yes."
            if single_scale
            else "4. A single A_scale does not cleanly cover -5V to -20V."
        ),
    ])
    if not single_scale:
        lines.append(
            "5. Low-bias and high-bias points prefer different A_scale values, so B, cutoff, or low-field smoothing still needs correction."
        )
    else:
        lines.append(
            "5. The same A_scale is preferred across the target biases in this IIC sweep; B/cutoff/low-field smoothing is less strongly indicated by this metric."
        )
    lines.append(
        "6. Recommended full coupled BV A_scale candidates (at most three): "
        + (", ".join(recommendations) if recommendations else "insufficient signal")
        + "."
    )
    lines.append(
        "7. Do not directly use sentaurus_fit_A_B for full coupled BV: prior BV matrix results showed it already makes -5V q integral Gava too large by about 1e4."
    )
    if low_signal_biases:
        lines.append("")
        lines.append("Low-signal biases: " + ", ".join(f"{bias}V" for bias in low_signal_biases) + ".")
    lines.append("")
    lines.append("## Candidate Ranking")
    lines.append("")
    lines.append("| rank | A_scale | median abs log10 region error |")
    lines.append("|---:|---:|---:|")
    for rank, (score, scale) in enumerate(ranked, start=1):
        lines.append(f"| {rank} | {scale} | {fmt(score)} |")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def output_fields() -> list[str]:
    return [
        "bias_V",
        "A_scale",
        "qG_full_self_A_per_um",
        "qG_full_sentaurus_A_per_um",
        "qG_full_ratio",
        "qG_junction_ratio",
        "qG_center_ratio",
        "qG_left_shoulder_ratio",
        "qG_right_shoulder_ratio",
        "max_G_self",
        "max_G_sentaurus",
        "best_region_error_metric",
        "signal_status",
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare/coarse_node_field_compare_aligned.csv",
    )
    parser.add_argument(
        "--elements-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_vm_vector/elements.csv",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=REPO / "build/diagnostics/vanoverstraeten_A_scale_iic_sweep.csv",
    )
    parser.add_argument(
        "--summary-md",
        type=Path,
        default=REPO / "build/diagnostics/vanoverstraeten_A_scale_iic_sweep_summary.md",
    )
    parser.add_argument(
        "--scales",
        default=",".join(fmt_scale(scale) for scale in DEFAULT_SCALES),
    )
    parser.add_argument(
        "--sentaurus-qg-floor",
        type=float,
        default=1.0e-30,
        help="Blank ratios when a Sentaurus qG integral is below this A/um floor.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data = load_compare(args.input_csv)
    elements = load_elements(args.elements_csv)
    if not data:
        raise ValueError(f"no rows loaded from {args.input_csv}")
    if not elements:
        raise ValueError(f"no triangular elements loaded from {args.elements_csv}")
    rows = compute_rows(data, elements, parse_scales(args.scales), args.sentaurus_qg_floor)
    write_rows(args.output_csv, rows, output_fields())
    write_summary(args.summary_md, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
