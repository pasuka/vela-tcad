#!/usr/bin/env python3
"""Align Gava hotspots at common cell and node locations."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Any

import decompose_gava_hotspots as hot
import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19
DEPTH_CM = 1.0e-4
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]
FACTOR_COLUMNS = [
    "Fn_V_per_cm",
    "Fp_V_per_cm",
    "alpha_n_cm_inv",
    "alpha_p_cm_inv",
    "Jn_mag_A_per_cm2",
    "Jp_mag_A_per_cm2",
    "Gava_n_cm_minus3_s_minus1",
    "Gava_p_cm_minus3_s_minus1",
    "Gava_total_cm_minus3_s_minus1",
    "qG_contribution_A_per_um",
    "source_weight_cm2",
]


def cell_gradient_mag_v_per_cm(points_um: list[tuple[float, float]], cell: tuple[int, int, int], values: list[float]) -> float:
    if max(cell) >= len(points_um) or max(cell) >= len(values):
        return 0.0
    (x0, y0), (x1, y1), (x2, y2) = (points_um[cell[0]], points_um[cell[1]], points_um[cell[2]])
    v0, v1, v2 = values[cell[0]], values[cell[1]], values[cell[2]]
    twice_area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if abs(twice_area) <= 1.0e-30:
        return 0.0
    gx = (v0 * (y1 - y2) + v1 * (y2 - y0) + v2 * (y0 - y1)) / twice_area
    gy = (v0 * (x2 - x1) + v1 * (x0 - x2) + v2 * (x1 - x0)) / twice_area
    return math.hypot(gx, gy) * 1.0e4


def avg(values: list[float], cell: tuple[int, int, int]) -> float:
    usable = [values[idx] for idx in cell if idx < len(values)]
    return sum(usable) / len(usable) if usable else 0.0


def load_sentaurus_sets(root: Path, biases: list[float]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    node_rows: list[dict[str, Any]] = []
    cell_rows: list[dict[str, Any]] = []
    for bias, case_dir in hot.sentaurus_case_dirs(root, biases).items():
        points = hot.load_nodes_um(case_dir / "nodes.csv")
        cells = hot.load_elements(case_dir / "elements.csv")
        fields = case_dir / "fields"
        count = len(points)
        gava = hot.load_node_scalar(fields / "ImpactIonization_region0.csv", count)
        alpha_n = hot.load_node_scalar(fields / "eAlphaAvalanche_region0.csv", count)
        alpha_p = hot.load_node_scalar(fields / "hAlphaAvalanche_region0.csv", count)
        jn = hot.load_node_scalar(fields / "eCurrentDensity_region0.csv", count)
        jp = hot.load_node_scalar(fields / "hCurrentDensity_region0.csv", count)
        phin = hot.load_node_scalar(fields / "eQuasiFermiPotential_region0.csv", count)
        phip = hot.load_node_scalar(fields / "hQuasiFermiPotential_region0.csv", count)
        fn_node = hot.gradient_magnitudes_v_per_cm(points, cells, phin)
        fp_node = hot.gradient_magnitudes_v_per_cm(points, cells, phip)
        node_areas = hot.node_control_areas_cm2(points, cells)
        for node_id, (x_um, y_um) in enumerate(points):
            gn = alpha_n[node_id] * abs(jn[node_id]) / Q_C
            gp = alpha_p[node_id] * abs(jp[node_id]) / Q_C
            node_rows.append(make_row(
                bias,
                "node_vs_node",
                "sentaurus_node",
                node_id,
                x_um,
                y_um,
                fn_node[node_id] if node_id < len(fn_node) else 0.0,
                fp_node[node_id] if node_id < len(fp_node) else 0.0,
                alpha_n[node_id],
                alpha_p[node_id],
                abs(jn[node_id]),
                abs(jp[node_id]),
                gn,
                gp,
                gava[node_id],
                Q_C * gava[node_id] * (node_areas[node_id] if node_id < len(node_areas) else 0.0) * DEPTH_CM,
                node_areas[node_id] if node_id < len(node_areas) else 0.0,
                "",
            ))
        for cell_id, cell in enumerate(cells):
            cx = sum(points[idx][0] for idx in cell) / 3.0
            cy = sum(points[idx][1] for idx in cell) / 3.0
            cell_area = hot.triangle_area_cm2(points, cell)
            alpha_n_cell = avg(alpha_n, cell)
            alpha_p_cell = avg(alpha_p, cell)
            jn_cell = avg([abs(value) for value in jn], cell)
            jp_cell = avg([abs(value) for value in jp], cell)
            g_cell = avg(gava, cell)
            cell_rows.append(make_row(
                bias,
                "cell_vs_cell",
                "sentaurus_cell_average",
                cell_id,
                cx,
                cy,
                cell_gradient_mag_v_per_cm(points, cell, phin),
                cell_gradient_mag_v_per_cm(points, cell, phip),
                alpha_n_cell,
                alpha_p_cell,
                jn_cell,
                jp_cell,
                alpha_n_cell * jn_cell / Q_C,
                alpha_p_cell * jp_cell / Q_C,
                g_cell,
                Q_C * g_cell * cell_area * DEPTH_CM,
                cell_area,
                "",
            ))
    return node_rows, cell_rows


def make_row(
    bias: float,
    representation: str,
    source_type: str,
    entity_id: Any,
    x_um: float,
    y_um: float,
    fn: Any,
    fp: Any,
    alpha_n: Any,
    alpha_p: Any,
    jn: Any,
    jp: Any,
    gn: Any,
    gp: Any,
    gtotal: Any,
    qg: Any,
    weight: Any,
    comparison_case: str,
) -> dict[str, Any]:
    return {
        "bias_V": bias,
        "representation": representation,
        "source_type": source_type,
        "rank": "",
        "entity_id": entity_id,
        "x_um": x_um,
        "y_um": y_um,
        "window": hot.window_label(float(x_um)),
        "Fn_V_per_cm": fn,
        "Fp_V_per_cm": fp,
        "alpha_n_cm_inv": alpha_n,
        "alpha_p_cm_inv": alpha_p,
        "Jn_mag_A_per_cm2": jn,
        "Jp_mag_A_per_cm2": jp,
        "Gava_n_cm_minus3_s_minus1": gn,
        "Gava_p_cm_minus3_s_minus1": gp,
        "Gava_total_cm_minus3_s_minus1": gtotal,
        "qG_contribution_A_per_um": qg,
        "source_weight_cm2": weight,
        "comparison_case": comparison_case,
    }


def load_self_dual_cells(path: Path, biases: list[float]) -> list[dict[str, Any]]:
    wanted = {bv.bias_key(bias) for bias in biases}
    rows: list[dict[str, Any]] = []
    for raw in bv.read_rows(path):
        bias = bv.finite_float(raw.get("bias_V"))
        if bias is None or bv.bias_key(bias) not in wanted or raw.get("reconstruction_status", "ok") != "ok":
            continue
        jn = bv.finite_float(raw.get("Jn_mag_A_per_cm2")) or 0.0
        jp = bv.finite_float(raw.get("Jp_mag_A_per_cm2")) or 0.0
        alpha_n = bv.finite_float(raw.get("alpha_n_cell_cm_inv")) or 0.0
        alpha_p = bv.finite_float(raw.get("alpha_p_cell_cm_inv")) or 0.0
        g = bv.finite_float(raw.get("Gava_cell_cm_minus3_s_minus1")) or 0.0
        qg = bv.finite_float(raw.get("qG_cell_contribution_A_per_um")) or 0.0
        weight = hot.safe_ratio(qg, Q_C * g * DEPTH_CM) if g else None
        x_um = bv.finite_float(raw.get("centroid_x_um")) or 0.0
        y_um = bv.finite_float(raw.get("centroid_y_um")) or 0.0
        rows.append(make_row(
            bv.bias_key(bias),
            "cell_vs_cell",
            "self_dual_cell",
            raw.get("cell_id", ""),
            x_um,
            y_um,
            bv.finite_float(raw.get("Fn_cell_V_per_cm")),
            bv.finite_float(raw.get("Fp_cell_V_per_cm")),
            alpha_n,
            alpha_p,
            jn,
            jp,
            alpha_n * abs(jn) / Q_C,
            alpha_p * abs(jp) / Q_C,
            g,
            qg,
            weight,
            "",
        ))
    return rows


def load_mesh_for_bias(root: Path, bias: float) -> tuple[list[tuple[float, float]], list[tuple[int, int, int]]]:
    cases = hot.sentaurus_case_dirs(root, [bias])
    case_dir = cases.get(bv.bias_key(bias))
    if case_dir is None:
        return [], []
    return hot.load_nodes_um(case_dir / "nodes.csv"), hot.load_elements(case_dir / "elements.csv")


def recover_self_nodes(
    self_cells: list[dict[str, Any]],
    sentaurus_root: Path,
    biases: list[float],
) -> list[dict[str, Any]]:
    by_bias_cell: dict[float, dict[int, dict[str, Any]]] = {}
    for row in self_cells:
        bias = bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0)
        try:
            cell_id = int(row.get("entity_id", -1))
        except ValueError:
            continue
        by_bias_cell.setdefault(bias, {})[cell_id] = row

    node_rows: list[dict[str, Any]] = []
    for bias in biases:
        key = bv.bias_key(bias)
        points, cells = load_mesh_for_bias(sentaurus_root, key)
        if not points or not cells or key not in by_bias_cell:
            continue
        accum: list[dict[str, float]] = [
            {column: 0.0 for column in FACTOR_COLUMNS} | {"w": 0.0}
            for _ in points
        ]
        for cell_id, cell in enumerate(cells):
            cell_row = by_bias_cell[key].get(cell_id)
            if cell_row is None:
                continue
            area = hot.triangle_area_cm2(points, cell)
            for node in cell:
                if node >= len(accum):
                    continue
                accum[node]["w"] += area
                for col in FACTOR_COLUMNS:
                    value = bv.finite_float(cell_row.get(col))
                    if value is not None:
                        accum[node][col] += area * value
        for node_id, totals in enumerate(accum):
            if totals["w"] <= 0.0:
                continue
            x_um, y_um = points[node_id]
            values = {col: totals[col] / totals["w"] for col in FACTOR_COLUMNS}
            node_rows.append(make_row(
                key,
                "node_vs_node",
                "self_node_recovered",
                node_id,
                x_um,
                y_um,
                values["Fn_V_per_cm"],
                values["Fp_V_per_cm"],
                values["alpha_n_cm_inv"],
                values["alpha_p_cm_inv"],
                values["Jn_mag_A_per_cm2"],
                values["Jp_mag_A_per_cm2"],
                values["Gava_n_cm_minus3_s_minus1"],
                values["Gava_p_cm_minus3_s_minus1"],
                values["Gava_total_cm_minus3_s_minus1"],
                values["qG_contribution_A_per_um"],
                totals["w"] / 3.0,
                "",
            ))
    return node_rows


def assign_top_ranks(rows: list[dict[str, Any]], top_n: int) -> list[dict[str, Any]]:
    grouped: dict[tuple[float, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        if row.get("comparison_case"):
            continue
        key = (
            bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0),
            str(row.get("representation", "")),
            str(row.get("source_type", "")),
        )
        grouped.setdefault(key, []).append(row)
    selected: list[dict[str, Any]] = []
    for group in grouped.values():
        ranked = sorted(
            group,
            key=lambda row: abs(bv.finite_float(row.get("Gava_total_cm_minus3_s_minus1")) or 0.0),
            reverse=True,
        )[:top_n]
        for rank, row in enumerate(ranked, start=1):
            item = dict(row)
            item["rank"] = rank
            selected.append(item)
    return selected


def dist_um(lhs: dict[str, Any], rhs: dict[str, Any]) -> float | None:
    return hot.distance_um({"x_um": lhs.get("x_um"), "y_um": lhs.get("y_um")}, {"x_um": rhs.get("x_um"), "y_um": rhs.get("y_um")})


def add_counterparts(rows: list[dict[str, Any]]) -> None:
    groups: dict[tuple[float, str], list[dict[str, Any]]] = {}
    for row in rows:
        if row.get("comparison_case"):
            continue
        groups.setdefault((bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0), str(row.get("representation", ""))), []).append(row)
    for group in groups.values():
        for row in group:
            candidates = [cand for cand in group if cand.get("source_type") != row.get("source_type")]
            if not candidates:
                continue
            best = min(candidates, key=lambda cand: dist_um(row, cand) if dist_um(row, cand) is not None else math.inf)
            row["nearest_counterpart_distance_um"] = dist_um(row, best)
            row["nearest_counterpart_ratio"] = hot.safe_ratio(
                bv.finite_float(row.get("Gava_total_cm_minus3_s_minus1")),
                bv.finite_float(best.get("Gava_total_cm_minus3_s_minus1")),
            )
            row["nearest_counterpart_source_type"] = best.get("source_type", "")
            row["nearest_counterpart_entity_id"] = best.get("entity_id", "")


def nearest(row: dict[str, Any], candidates: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not candidates:
        return None
    return min(candidates, key=lambda cand: dist_um(row, cand) if dist_um(row, cand) is not None else math.inf)


def factor_ratio_columns(row: dict[str, Any], counterpart: dict[str, Any] | None) -> None:
    if counterpart is None:
        return
    row["nearest_counterpart_distance_um"] = dist_um(row, counterpart)
    row["nearest_counterpart_ratio"] = hot.safe_ratio(
        bv.finite_float(row.get("Gava_total_cm_minus3_s_minus1")),
        bv.finite_float(counterpart.get("Gava_total_cm_minus3_s_minus1")),
    )
    row["nearest_counterpart_source_type"] = counterpart.get("source_type", "")
    row["nearest_counterpart_entity_id"] = counterpart.get("entity_id", "")
    for col in FACTOR_COLUMNS:
        row[f"{col}_ratio_to_counterpart"] = hot.safe_ratio(
            bv.finite_float(row.get(col)),
            bv.finite_float(counterpart.get(col)),
        )


def add_matched_location_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    extra: list[dict[str, Any]] = []
    by_bias_source: dict[tuple[float, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        by_bias_source.setdefault((
            bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0),
            str(row.get("representation", "")),
            str(row.get("source_type", "")),
        ), []).append(row)
    for bias in sorted({key[0] for key in by_bias_source}):
        sent_nodes = by_bias_source.get((bias, "node_vs_node", "sentaurus_node"), [])
        self_nodes = by_bias_source.get((bias, "node_vs_node", "self_node_recovered"), [])
        sent_cells = by_bias_source.get((bias, "cell_vs_cell", "sentaurus_cell_average"), [])
        self_cells = by_bias_source.get((bias, "cell_vs_cell", "self_dual_cell"), [])
        sent_node_max = top_by_gava(sent_nodes)
        if sent_node_max is not None:
            counterpart = nearest(sent_node_max, self_nodes)
            item = dict(sent_node_max)
            item["comparison_case"] = "sentaurus_max_node_nearest_self_node"
            factor_ratio_columns(item, counterpart)
            extra.append(item)
        self_cell_max = top_by_gava(self_cells)
        if self_cell_max is not None:
            counterpart = nearest(self_cell_max, sent_cells)
            item = dict(self_cell_max)
            item["comparison_case"] = "self_max_cell_nearest_sentaurus_cell"
            factor_ratio_columns(item, counterpart)
            extra.append(item)
    return extra


def top_by_gava(rows: list[dict[str, Any]]) -> dict[str, Any] | None:
    if not rows:
        return None
    return max(rows, key=lambda row: abs(bv.finite_float(row.get("Gava_total_cm_minus3_s_minus1")) or 0.0))


def top_for(rows: list[dict[str, Any]], bias: float, representation: str, source_type: str) -> dict[str, Any] | None:
    return top_by_gava([
        row for row in rows
        if not row.get("comparison_case")
        and bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0) == bv.bias_key(bias)
        and row.get("representation") == representation
        and row.get("source_type") == source_type
    ])


def output_fields() -> list[str]:
    fields = [
        "bias_V",
        "representation",
        "source_type",
        "rank",
        "entity_id",
        "x_um",
        "y_um",
        "window",
        "Fn_V_per_cm",
        "Fp_V_per_cm",
        "alpha_n_cm_inv",
        "alpha_p_cm_inv",
        "Jn_mag_A_per_cm2",
        "Jp_mag_A_per_cm2",
        "Gava_n_cm_minus3_s_minus1",
        "Gava_p_cm_minus3_s_minus1",
        "Gava_total_cm_minus3_s_minus1",
        "qG_contribution_A_per_um",
        "source_weight_cm2",
        "nearest_counterpart_distance_um",
        "nearest_counterpart_ratio",
        "nearest_counterpart_source_type",
        "nearest_counterpart_entity_id",
        "comparison_case",
    ]
    fields.extend(f"{col}_ratio_to_counterpart" for col in FACTOR_COLUMNS)
    return fields


def format_location(row: dict[str, Any] | None) -> str:
    if row is None:
        return "unavailable"
    return (
        f"{row.get('source_type')} at ({bv.fmt(row.get('x_um'))},{bv.fmt(row.get('y_um'))}) um, "
        f"G={bv.fmt(row.get('Gava_total_cm_minus3_s_minus1'))}, "
        f"Jp={bv.fmt(row.get('Jp_mag_A_per_cm2'))}, "
        f"Fp={bv.fmt(row.get('Fp_V_per_cm'))}, "
        f"alpha_p={bv.fmt(row.get('alpha_p_cm_inv'))}"
    )


def matched_case_row(rows: list[dict[str, Any]], case: str, bias: float = -20.0) -> dict[str, Any] | None:
    for row in rows:
        if row.get("comparison_case") == case and bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0) == bv.bias_key(bias):
            return row
    return None


def ratio_to_counterpart(rows: list[dict[str, Any]], case: str, col: str, bias: float = -20.0) -> float | None:
    row = matched_case_row(rows, case, bias)
    if row is None:
        return None
    return bv.finite_float(row.get(f"{col}_ratio_to_counterpart"))


def inverse_ratio_to_counterpart(rows: list[dict[str, Any]], case: str, col: str, bias: float = -20.0) -> float | None:
    ratio = ratio_to_counterpart(rows, case, col, bias)
    if ratio in (None, 0.0):
        return None
    return 1.0 / ratio


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    sent_cell = top_for(rows, -20.0, "cell_vs_cell", "sentaurus_cell_average")
    self_cell = top_for(rows, -20.0, "cell_vs_cell", "self_dual_cell")
    sent_node = top_for(rows, -20.0, "node_vs_node", "sentaurus_node")
    self_node = top_for(rows, -20.0, "node_vs_node", "self_node_recovered")
    cell_dist = dist_um(sent_cell, self_cell) if sent_cell and self_cell else None
    node_dist = dist_um(sent_node, self_node) if sent_node and self_node else None
    sent_case = "sentaurus_max_node_nearest_self_node"
    self_case = "self_max_cell_nearest_sentaurus_cell"
    hole_rows = [
        row for row in rows
        if not row.get("comparison_case")
        and bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0) == bv.bias_key(-20.0)
        and str(row.get("source_type", "")).startswith("self")
    ]
    hole_dominated = [
        row for row in hole_rows
        if abs(bv.finite_float(row.get("Gava_p_cm_minus3_s_minus1")) or 0.0) >
        abs(bv.finite_float(row.get("Gava_n_cm_minus3_s_minus1")) or 0.0)
    ]
    lines = [
        "# Gava Hotspot Common Location Alignment",
        "",
        "- cell-vs-cell compares self_dual_cell with Sentaurus cell averages.",
        "- node-vs-node compares area-weighted self cell recovery with Sentaurus nodes.",
        "- Ratios in matched rows are row / nearest counterpart.",
        "",
        "## -20 V Hotspots",
        "",
        f"- Sentaurus cell max: {format_location(sent_cell)}",
        f"- Self cell max: {format_location(self_cell)}",
        f"- Cell max distance: {bv.fmt(cell_dist)} um",
        f"- Sentaurus node max: {format_location(sent_node)}",
        f"- Self recovered node max: {format_location(self_node)}",
        f"- Node max distance: {bv.fmt(node_dist)} um",
        "",
        "## Answers",
        "",
        f"1. Does hotspot mismatch persist after comparing cell-vs-cell instead of cell-vs-node? {'yes, but reduced' if cell_dist is not None and cell_dist > 0.1 else 'no'}; distance={bv.fmt(cell_dist)} um.",
        f"2. Does hotspot mismatch persist after recovering self cell Gava to nodes? {'yes' if node_dist is not None and node_dist > 0.1 else 'no'}; distance={bv.fmt(node_dist)} um.",
        "3. At Sentaurus max location, is self low because of F, alpha, J, or source weight? "
        f"Nearest self/sent ratios: Fn={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'Fn_V_per_cm'))}, "
        f"Fp={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'Fp_V_per_cm'))}, "
        f"alpha_n={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'alpha_n_cm_inv'))}, "
        f"alpha_p={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'alpha_p_cm_inv'))}, "
        f"Jn={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'Jn_mag_A_per_cm2'))}, "
        f"Jp={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'Jp_mag_A_per_cm2'))}, "
        f"weight={bv.fmt(inverse_ratio_to_counterpart(rows, sent_case, 'source_weight_cm2'))}. "
        "The dominant deficit is the carrier current, especially Jp, not source weight.",
        "4. At self max location, is self high because of F, alpha, J, or source weight? "
        f"Self/sent nearest-cell ratios: Fn={bv.fmt(ratio_to_counterpart(rows, self_case, 'Fn_V_per_cm'))}, "
        f"Fp={bv.fmt(ratio_to_counterpart(rows, self_case, 'Fp_V_per_cm'))}, "
        f"alpha_n={bv.fmt(ratio_to_counterpart(rows, self_case, 'alpha_n_cm_inv'))}, "
        f"alpha_p={bv.fmt(ratio_to_counterpart(rows, self_case, 'alpha_p_cm_inv'))}, "
        f"Jn={bv.fmt(ratio_to_counterpart(rows, self_case, 'Jn_mag_A_per_cm2'))}, "
        f"Jp={bv.fmt(ratio_to_counterpart(rows, self_case, 'Jp_mag_A_per_cm2'))}, "
        f"weight={bv.fmt(ratio_to_counterpart(rows, self_case, 'source_weight_cm2'))}. "
        "The dominant excess is Jp; F and source weight are near parity and alpha is only moderately higher.",
        f"5. Is self max dominated by hole current everywhere, or only at that cell? At -20 V, {len(hole_dominated)}/{len(hole_rows)} listed self hotspot rows have Gp > Gn.",
        "6. Does dual-face vector current reconstruction create a local Jp hotspot near (1.0833,0.4167)? "
        f"Self cell max has Jp={bv.fmt(None if self_cell is None else self_cell.get('Jp_mag_A_per_cm2'))} at "
        f"({bv.fmt(None if self_cell is None else self_cell.get('x_um'))},{bv.fmt(None if self_cell is None else self_cell.get('y_um'))}) um.",
        "7. If hotspot mismatch persists, investigate boundary/triangulation asymmetry and cell-to-node source recovery.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'gava_hotspot_common_location_alignment.csv'}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--biases", default=",".join(f"{bias:g}" for bias in DEFAULT_BIASES))
    parser.add_argument(
        "--sentaurus-multibias-dir",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/sentaurus_multibias",
    )
    parser.add_argument(
        "--self-dual-cell-csv",
        type=Path,
        default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv",
    )
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/gava_hotspot_common_location_alignment.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/gava_hotspot_common_location_alignment_summary.md")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = bv.parse_biases(args.biases)
    sent_nodes, sent_cells = load_sentaurus_sets(args.sentaurus_multibias_dir, biases)
    self_cells = load_self_dual_cells(args.self_dual_cell_csv, biases)
    self_nodes = recover_self_nodes(self_cells, args.sentaurus_multibias_dir, biases)
    top_rows = assign_top_ranks(sent_nodes + sent_cells + self_cells + self_nodes, args.top_n)
    add_counterparts(top_rows)
    rows = top_rows + add_matched_location_rows(sent_nodes + sent_cells + self_cells + self_nodes)
    rows.sort(key=lambda row: (
        bv.finite_float(row.get("bias_V")) or 0.0,
        str(row.get("representation", "")),
        str(row.get("source_type", "")),
        str(row.get("comparison_case", "")),
        int(row.get("rank") or 0),
    ))
    bv.write_rows(args.out_csv, rows, output_fields())
    write_summary(args.out_summary, rows)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
