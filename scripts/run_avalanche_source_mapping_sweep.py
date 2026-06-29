#!/usr/bin/env python3
"""Offline AvalancheGeneration source mapping sweep for the PN2D A2/B1.05 BV state."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
DEPTH_M = 1.0e-6
MODES = [
    "node_F_node_alpha_node_G",
    "edge_F_edge_alpha_edge_G_to_node",
    "cell_F_cell_alpha_cell_G_to_node",
    "cell_F_cell_alpha_cell_G_integral_only",
]


def default_case_dir() -> Path:
    return (
        REPO / "build/diagnostics/vanoverstraeten_A2_B105_spatial_validation_case"
        / "BV-A2-B1p05-spatial"
    )


def parse_args() -> argparse.Namespace:
    case_dir = default_case_dir()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, default=case_dir)
    parser.add_argument("--sweep-csv", type=Path, default=case_dir / "BV-A2-B1p05-spatial.csv")
    parser.add_argument(
        "--sentaurus-multibias-dir",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/sentaurus_multibias",
    )
    parser.add_argument(
        "--reference-current-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/sentaurus_coarse_bv_reference_aligned.csv",
    )
    parser.add_argument("--out-dir", type=Path, default=REPO / "build/diagnostics")
    parser.add_argument("--bias-points", default="0,-5,-10,-16,-18,-20")
    return parser.parse_args()


def ratio(value: Any, reference: Any) -> float | None:
    lhs = bv.finite_float(value)
    rhs = bv.finite_float(reference)
    if lhs is None or rhs is None or rhs == 0.0:
        return None
    return lhs / rhs


def point_field(data: bv.VtkData, names: list[str]) -> list[float]:
    return bv.field_values(data, names)


def reconstructed_node_generation(data: bv.VtkData) -> list[float]:
    e_alpha = point_field(data, ["ElectronAlphaAvalanche"])
    h_alpha = point_field(data, ["HoleAlphaAvalanche"])
    jn = point_field(data, ["ElectronCurrentDensityVector", "J_n_total"])
    jp = point_field(data, ["HoleCurrentDensityVector", "J_p_total"])
    count = min(len(e_alpha), len(h_alpha), len(jn), len(jp), len(data.points))
    if count == 0:
        return point_field(data, ["AvalancheGeneration"])
    return [
        (abs(e_alpha[i]) * abs(jn[i]) + abs(h_alpha[i]) * abs(jp[i])) / bv.E_CHARGE_C
        for i in range(count)
    ]


def edge_mapped_generation(data: bv.VtkData, node_g: list[float]) -> list[float]:
    sums = [0.0] * len(data.points)
    weights = [0.0] * len(data.points)
    seen: set[tuple[int, int]] = set()
    for cell in data.cells:
        for a, b in ((cell[0], cell[1]), (cell[1], cell[2]), (cell[2], cell[0])):
            edge = (min(a, b), max(a, b))
            if edge in seen or max(edge) >= len(node_g):
                continue
            seen.add(edge)
            length = distance_m(data.points[a], data.points[b])
            value = 0.5 * (node_g[a] + node_g[b])
            for node in edge:
                sums[node] += value * length
                weights[node] += length
    return [sums[i] / weights[i] if weights[i] > 0.0 else 0.0 for i in range(len(data.points))]


def cell_mapped_generation(data: bv.VtkData, node_g: list[float]) -> tuple[list[float], list[float]]:
    sums = [0.0] * len(data.points)
    weights = [0.0] * len(data.points)
    cell_values: list[float] = []
    for cell in data.cells:
        if max(cell) >= len(node_g):
            cell_values.append(0.0)
            continue
        area = bv.triangle_area_m2(data.points, cell)
        value = sum(node_g[node] for node in cell) / 3.0
        cell_values.append(value)
        for node in cell:
            sums[node] += value * area / 3.0
            weights[node] += area / 3.0
    return [sums[i] / weights[i] if weights[i] > 0.0 else 0.0 for i in range(len(data.points))], cell_values

def scale_node_values_to_full_integral(
    data: bv.VtkData,
    values: list[float],
    target_full_integral: float | None,
) -> list[float]:
    raw_full = integrate_node_generation(data.points, data.cells, values, None)
    factor = scale_factor(raw_full, target_full_integral)
    return [value * factor for value in values]


def scale_cell_values_to_full_integral(
    data: bv.VtkData,
    values: list[float],
    target_full_integral: float | None,
) -> list[float]:
    raw_full = integrate_cell_generation(data.points, data.cells, values, None)
    factor = scale_factor(raw_full, target_full_integral)
    return [value * factor for value in values]


def scale_factor(raw_full: float | None, target_full: float | None) -> float:
    if raw_full is None or target_full is None or raw_full == 0.0:
        return 1.0
    return target_full / raw_full


def distance_m(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def integrate_node_generation(
    points: list[tuple[float, float]],
    cells: list[tuple[int, int, int]],
    values: list[float],
    window: tuple[float, float] | None,
) -> float:
    return bv.integrate_generation(points, cells, values, window)


def integrate_cell_generation(
    points: list[tuple[float, float]],
    cells: list[tuple[int, int, int]],
    cell_values: list[float],
    window: tuple[float, float] | None,
) -> float:
    total = 0.0
    for idx, cell in enumerate(cells):
        if idx >= len(cell_values):
            continue
        cx = bv.centroid_x_um(points, cell)
        if window is not None and not (window[0] <= cx <= window[1]):
            continue
        total += cell_values[idx] * bv.triangle_area_m2(points, cell) * DEPTH_M
    return total


def node_max(points: list[tuple[float, float]], values: list[float]) -> tuple[float | None, float | None, float | None]:
    if not points or not values:
        return None, None, None
    indexed = [(abs(value), idx) for idx, value in enumerate(values) if idx < len(points)]
    if not indexed:
        return None, None, None
    value, idx = max(indexed, key=lambda item: item[0])
    x_m, y_m = points[idx]
    return value, x_m * 1.0e6, y_m * 1.0e6


def cell_max(
    points: list[tuple[float, float]],
    cells: list[tuple[int, int, int]],
    values: list[float],
) -> tuple[float | None, float | None, float | None]:
    indexed = [(abs(value), idx) for idx, value in enumerate(values) if idx < len(cells)]
    if not indexed:
        return None, None, None
    value, idx = max(indexed, key=lambda item: item[0])
    cell = cells[idx]
    x_um = sum(points[node][0] for node in cell) / 3.0 * 1.0e6
    y_um = sum(points[node][1] for node in cell) / 3.0 * 1.0e6
    return value, x_um, y_um


def load_sentaurus_reference(root: Path) -> dict[float, dict[str, Any]]:
    result: dict[float, dict[str, Any]] = {}
    if not root.exists():
        return result
    for case_dir in root.iterdir():
        if not case_dir.is_dir() or not case_dir.name.startswith("sentaurus_"):
            continue
        token = case_dir.name.removeprefix("sentaurus_").removesuffix("v").replace("p", ".")
        try:
            bias = float(token)
        except ValueError:
            continue
        nodes = bv.load_sentaurus_nodes(case_dir / "nodes.csv")
        cells = bv.load_sentaurus_elements(case_dir / "elements.csv")
        values_cm3_s = bv.load_sentaurus_scalar(case_dir / "fields" / "ImpactIonization_region0.csv")
        values_m3_s = [value * 1.0e6 for value in values_cm3_s]
        item: dict[str, Any] = {}
        for name, window in bv.WINDOWS.items():
            item[f"sentaurus_integral_{name}"] = bv.integrate_generation(nodes, cells, values_m3_s, window)
        max_value, x_um, y_um = node_max(nodes, values_cm3_s)
        item["max_Gava_sentaurus"] = max_value
        item["max_Gava_sentaurus_x_um"] = x_um
        item["max_Gava_sentaurus_y_um"] = y_um
        result[bv.bias_key(bias)] = item
    return result


def discover_vtks(case_dir: Path) -> dict[float, Path]:
    return bv.discover_vtks(case_dir)


def load_curve_rows(path: Path) -> dict[float, dict[str, str]]:
    return {
        bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0): row
        for row in bv.read_rows(path)
    }


def make_mapping_rows(
    bias: float,
    data: bv.VtkData,
    sentaurus: dict[str, Any],
    curve_row: dict[str, str],
    terminal_ref: float | None,
) -> list[dict[str, Any]]:
    node_g = reconstructed_node_generation(data)
    existing_g = point_field(data, ["AvalancheGeneration"])
    target_full_integral = (
        integrate_node_generation(data.points, data.cells, existing_g, None)
        if existing_g else None
    )
    node_g = scale_node_values_to_full_integral(data, node_g, target_full_integral)
    edge_g = edge_mapped_generation(data, node_g)
    cell_node_g, cell_values = cell_mapped_generation(data, node_g)
    edge_g = scale_node_values_to_full_integral(data, edge_g, target_full_integral)
    cell_node_g = scale_node_values_to_full_integral(data, cell_node_g, target_full_integral)
    cell_values = scale_cell_values_to_full_integral(data, cell_values, target_full_integral)
    mappings = {
        "node_F_node_alpha_node_G": ("node", node_g),
        "edge_F_edge_alpha_edge_G_to_node": ("node", edge_g),
        "cell_F_cell_alpha_cell_G_to_node": ("node", cell_node_g),
        "cell_F_cell_alpha_cell_G_integral_only": ("cell", cell_values),
    }
    rows = []
    terminal_current = bv.finite_float(
        curve_row.get("current_total_A_per_um") or curve_row.get("current_total"))
    for mode, (kind, values) in mappings.items():
        row: dict[str, Any] = {
            "bias_V": bias,
            "mapping_mode": mode,
            "terminal_current_ratio": ratio(terminal_current, terminal_ref),
            "convergence_status": "converged" if curve_row.get("converged") == "1" else "not_converged",
        }
        if kind == "cell":
            max_value, x_um, y_um = cell_max(data.points, data.cells, values)
        else:
            max_value, x_um, y_um = node_max(data.points, values)
        row["max_Gava_self"] = cm3_per_s(max_value)
        row["max_Gava_self_x_um"] = x_um
        row["max_Gava_self_y_um"] = y_um
        row["max_Gava_sentaurus"] = sentaurus.get("max_Gava_sentaurus")
        row["max_Gava_sentaurus_x_um"] = sentaurus.get("max_Gava_sentaurus_x_um")
        row["max_Gava_sentaurus_y_um"] = sentaurus.get("max_Gava_sentaurus_y_um")
        row["max_Gava_distance_um"] = distance_um(
            x_um, y_um, sentaurus.get("max_Gava_sentaurus_x_um"), sentaurus.get("max_Gava_sentaurus_y_um"))
        for name, window in bv.WINDOWS.items():
            if kind == "cell":
                integral = integrate_cell_generation(data.points, data.cells, values, window)
            else:
                integral = integrate_node_generation(data.points, data.cells, values, window)
            sent_integral = sentaurus.get(f"sentaurus_integral_{name}")
            row[f"{short_window(name)}_qG_ratio"] = ratio(integral, sent_integral)
        rows.append(row)
    return rows


def short_window(name: str) -> str:
    return {
        "full_semiconductor": "full",
        "junction_window": "junction",
        "center_window": "center",
        "left_shoulder": "left_shoulder",
        "right_shoulder": "right_shoulder",
    }[name]


def cm3_per_s(value: Any) -> float | None:
    parsed = bv.finite_float(value)
    return None if parsed is None else parsed * 1.0e-6


def distance_um(x0: Any, y0: Any, x1: Any, y1: Any) -> float | None:
    values = [bv.finite_float(value) for value in (x0, y0, x1, y1)]
    if any(value is None for value in values):
        return None
    assert values[0] is not None and values[1] is not None
    assert values[2] is not None and values[3] is not None
    return math.hypot(values[0] - values[2], values[1] - values[3])


def output_fields() -> list[str]:
    return [
        "bias_V", "mapping_mode",
        "full_qG_ratio", "junction_qG_ratio", "center_qG_ratio",
        "left_shoulder_qG_ratio", "right_shoulder_qG_ratio",
        "max_Gava_self", "max_Gava_sentaurus",
        "max_Gava_self_x_um", "max_Gava_self_y_um",
        "max_Gava_sentaurus_x_um", "max_Gava_sentaurus_y_um",
        "max_Gava_distance_um",
        "terminal_current_ratio", "convergence_status",
    ]


def score_mode(rows: list[dict[str, Any]], mode: str, columns: list[str]) -> float | None:
    errors: list[float] = []
    for row in rows:
        if row.get("mapping_mode") != mode:
            continue
        bias = bv.finite_float(row.get("bias_V"))
        if bias == 0.0:
            continue
        for column in columns:
            value = bv.finite_float(row.get(column))
            if value is not None and value != 0.0:
                errors.append(abs(math.log10(abs(value))))
    return median(errors)


def median(values: list[float]) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def best_mode(rows: list[dict[str, Any]], columns: list[str]) -> str:
    best = ""
    best_score = math.inf
    for mode in MODES:
        score = score_mode(rows, mode, columns)
        if score is not None and score < best_score:
            best = mode
            best_score = score
    return best or "insufficient_data"


def best_position_mode(rows: list[dict[str, Any]]) -> str:
    scores = {}
    for mode in MODES:
        values = [
            value for value in (bv.finite_float(row.get("max_Gava_distance_um")) for row in rows
                                if row.get("mapping_mode") == mode and bv.finite_float(row.get("bias_V")) != 0.0)
            if value is not None
        ]
        if values:
            scores[mode] = median(values)
    if not scores:
        return "insufficient_data"
    return min(scores, key=lambda mode: scores[mode] if scores[mode] is not None else math.inf)


def mode_pair_difference(rows: list[dict[str, Any]], mode_a: str, mode_b: str) -> float | None:
    by_key = {(row["mapping_mode"], bv.bias_key(float(row["bias_V"]))): row for row in rows}
    diffs: list[float] = []
    for bias in {key[1] for key in by_key if key[1] != 0.0}:
        a = by_key.get((mode_a, bias), {})
        b = by_key.get((mode_b, bias), {})
        for column in ("full_qG_ratio", "junction_qG_ratio", "center_qG_ratio",
                       "left_shoulder_qG_ratio", "right_shoulder_qG_ratio"):
            va = bv.finite_float(a.get(column))
            vb = bv.finite_float(b.get(column))
            if va is not None and vb is not None:
                diffs.append(abs(math.log10(abs(va) + 1.0e-300) - math.log10(abs(vb) + 1.0e-300)))
    return median(diffs)


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    best_pos = best_position_mode(rows)
    best_shoulder = best_mode(rows, ["left_shoulder_qG_ratio", "right_shoulder_qG_ratio"])
    best_total = best_mode(rows, ["full_qG_ratio", "junction_qG_ratio"])
    cell_vs_node = mode_pair_difference(
        rows, "cell_F_cell_alpha_cell_G_integral_only", "cell_F_cell_alpha_cell_G_to_node")
    edge_or_cell_improves = best_shoulder in {
        "edge_F_edge_alpha_edge_G_to_node",
        "cell_F_cell_alpha_cell_G_to_node",
        "cell_F_cell_alpha_cell_G_integral_only",
    }

    lines = [
        "# Avalanche Source Mapping Sweep",
        "",
        "- Diagnostic option only; default coupled solve path is not changed.",
        "- Fixed physics context: VanOverstraeten/de Man, A_scale=2.0, B_scale=1.05, parameter_set=default.",
        "- Driving force context: GradQuasiFermi from the converged BV state.",
        "- Each bias is full-domain normalized to the current A2_B105 Vela AvalancheGeneration integral, so this isolates spatial redistribution rather than global source magnitude.",
        "- box_face_F_alpha_G_to_node was not emitted because this offline export has no box/control-volume face table.",
        "",
        "## Mode Scores",
        "",
        "| mapping mode | position median distance um | shoulder score | total score |",
        "|---|---:|---:|---:|",
    ]
    for mode in MODES:
        distances = [
            value for value in (bv.finite_float(row.get("max_Gava_distance_um")) for row in rows
                                if row.get("mapping_mode") == mode and bv.finite_float(row.get("bias_V")) != 0.0)
            if value is not None
        ]
        lines.append(
            f"| {mode} | {bv.fmt(median(distances))} | "
            f"{bv.fmt(score_mode(rows, mode, ['left_shoulder_qG_ratio', 'right_shoulder_qG_ratio']))} | "
            f"{bv.fmt(score_mode(rows, mode, ['full_qG_ratio', 'junction_qG_ratio']))} |"
        )
    lines.extend([
        "",
        "## Answers",
        "",
        f"1. The mapping mode with max Gava position closest to Sentaurus is {best_pos}.",
        f"2. The mapping mode with left/right shoulder qG closest to Sentaurus is {best_shoulder}.",
        f"3. cell integral only vs node mapped qG difference is {bv.fmt(cell_vs_node)} median log10-ratio units; {'large' if cell_vs_node is not None and cell_vs_node > 0.3 else 'not large'} by this diagnostic threshold.",
    ])
    if edge_or_cell_improves:
        lines.append(
            "4. An edge/cell mapping is best for shoulders, so the primary suspect is source mapping rather than global A/B parameters."
        )
    else:
        lines.append(
            "4. Edge/cell mapping does not clearly improve shoulders over node mapping in this sweep; source mapping remains the diagnostic target."
        )
    lines.extend([
        f"5. If all mapping modes remain poor after detailed review, run B_low/B_high split plus cutoff/RefDens checks; current best total mode is {best_total}.",
        "6. Keep this as a diagnostic option and do not alter the default solve path.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'avalanche_source_mapping_sweep.csv'}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    biases = bv.parse_biases(args.bias_points)
    vtks = discover_vtks(args.case_dir)
    curve_rows = load_curve_rows(args.sweep_csv)
    sentaurus = load_sentaurus_reference(args.sentaurus_multibias_dir)
    terminal_refs = bv.load_reference_current(args.reference_current_csv)

    rows: list[dict[str, Any]] = []
    for bias in biases:
        vtk = vtks.get(bv.bias_key(bias))
        if not vtk:
            continue
        data = bv.parse_vtk(vtk)
        rows.extend(make_mapping_rows(
            bias,
            data,
            sentaurus.get(bv.bias_key(bias), {}),
            curve_rows.get(bv.bias_key(bias), {}),
            terminal_refs.get(bv.bias_key(bias)),
        ))

    csv_path = args.out_dir / "avalanche_source_mapping_sweep.csv"
    summary_path = args.out_dir / "avalanche_source_mapping_sweep_summary.md"
    bv.write_rows(csv_path, rows, output_fields())
    write_summary(summary_path, rows)
    print(csv_path)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
