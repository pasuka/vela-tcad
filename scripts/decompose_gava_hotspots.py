#!/usr/bin/env python3
"""Decompose top AvalancheGeneration hotspots for Sentaurus and Vela self data."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path
from typing import Any

import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19
DEPTH_CM = 1.0e-4
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]


def window_label(x_um: float) -> str:
    if 0.9 <= x_um <= 1.1:
        return "center"
    if 0.7 <= x_um <= 0.85:
        return "left_shoulder"
    if 1.15 <= x_um <= 1.3:
        return "right_shoulder"
    if 0.7 <= x_um <= 1.3:
        return "junction"
    return "full_semiconductor"


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator in (None, 0.0):
        return None
    return numerator / denominator


def fmt(value: Any) -> str:
    return bv.fmt(value)


def distance_um(a: dict[str, Any], b: dict[str, Any]) -> float | None:
    ax = bv.finite_float(a.get("x_um"))
    ay = bv.finite_float(a.get("y_um"))
    bx = bv.finite_float(b.get("x_um"))
    by = bv.finite_float(b.get("y_um"))
    if ax is None or ay is None or bx is None or by is None:
        return None
    return math.hypot(ax - bx, ay - by)


def triangle_area_cm2(points_um: list[tuple[float, float]], cell: tuple[int, int, int]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = (points_um[cell[0]], points_um[cell[1]], points_um[cell[2]])
    area_um2 = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) * 0.5
    return area_um2 * 1.0e-8


def node_control_areas_cm2(points_um: list[tuple[float, float]], cells: list[tuple[int, int, int]]) -> list[float]:
    areas = [0.0] * len(points_um)
    for cell in cells:
        area = triangle_area_cm2(points_um, cell)
        for node in cell:
            if node < len(areas):
                areas[node] += area / 3.0
    return areas


def load_nodes_um(path: Path) -> list[tuple[float, float]]:
    rows = bv.read_rows(path)
    points: list[tuple[float, float]] = []
    for row in rows:
        points.append((bv.finite_float(row.get("x_um")) or 0.0, bv.finite_float(row.get("y_um")) or 0.0))
    return points


def load_elements(path: Path) -> list[tuple[int, int, int]]:
    cells: list[tuple[int, int, int]] = []
    for row in bv.read_rows(path):
        try:
            cells.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
        except (KeyError, ValueError):
            continue
    return cells


def load_node_scalar(path: Path, count: int) -> list[float]:
    values = [0.0] * count
    for row in bv.read_rows(path):
        try:
            idx = int(row.get("node_id", row.get("id", "0")))
        except ValueError:
            continue
        if 0 <= idx < count:
            values[idx] = bv.finite_float(row.get("component0") or row.get("value")) or 0.0
    return values


def gradient_magnitudes_v_per_cm(
    points_um: list[tuple[float, float]],
    cells: list[tuple[int, int, int]],
    values_v: list[float],
) -> list[float]:
    neighbors = [set() for _ in points_um]
    for a, b, c in cells:
        for u, v in ((a, b), (b, c), (c, a)):
            if u < len(neighbors) and v < len(neighbors):
                neighbors[u].add(v)
                neighbors[v].add(u)
    out = [0.0] * len(points_um)
    for i, nbrs in enumerate(neighbors):
        if not nbrs:
            continue
        x0, y0 = points_um[i]
        v0 = values_v[i] if i < len(values_v) else 0.0
        a00 = a01 = a11 = b0 = b1 = 0.0
        for j in nbrs:
            if j >= len(points_um):
                continue
            dx = points_um[j][0] - x0
            dy = points_um[j][1] - y0
            dist = math.hypot(dx, dy)
            if dist <= 1.0e-30:
                continue
            weight = 1.0 / dist
            dv = (values_v[j] if j < len(values_v) else 0.0) - v0
            a00 += weight * dx * dx
            a01 += weight * dx * dy
            a11 += weight * dy * dy
            b0 += weight * dx * dv
            b1 += weight * dy * dv
        det = a00 * a11 - a01 * a01
        if abs(det) <= 1.0e-30:
            continue
        gx_v_per_um = (b0 * a11 - b1 * a01) / det
        gy_v_per_um = (a00 * b1 - a01 * b0) / det
        out[i] = math.hypot(gx_v_per_um, gy_v_per_um) * 1.0e4
    return out


def sentaurus_case_dirs(root: Path, biases: list[float]) -> dict[float, Path]:
    wanted = {bv.bias_key(bias): bias for bias in biases}
    result: dict[float, Path] = {}
    pattern = re.compile(r"sentaurus_(?P<bias>[-+0-9p.]+)v$")
    for path in root.iterdir() if root.exists() else []:
        if not path.is_dir():
            continue
        match = pattern.match(path.name)
        if not match:
            continue
        try:
            bias = float(match.group("bias").replace("p", "."))
        except ValueError:
            continue
        key = bv.bias_key(bias)
        if key in wanted:
            result[key] = path
    return result


def load_sentaurus_hotspots(root: Path, biases: list[float], top_n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias_key, case_dir in sentaurus_case_dirs(root, biases).items():
        points = load_nodes_um(case_dir / "nodes.csv")
        cells = load_elements(case_dir / "elements.csv")
        count = len(points)
        fields = case_dir / "fields"
        gava = load_node_scalar(fields / "ImpactIonization_region0.csv", count)
        alpha_n = load_node_scalar(fields / "eAlphaAvalanche_region0.csv", count)
        alpha_p = load_node_scalar(fields / "hAlphaAvalanche_region0.csv", count)
        jn = load_node_scalar(fields / "eCurrentDensity_region0.csv", count)
        jp = load_node_scalar(fields / "hCurrentDensity_region0.csv", count)
        phin = load_node_scalar(fields / "eQuasiFermiPotential_region0.csv", count)
        phip = load_node_scalar(fields / "hQuasiFermiPotential_region0.csv", count)
        fn = gradient_magnitudes_v_per_cm(points, cells, phin)
        fp = gradient_magnitudes_v_per_cm(points, cells, phip)
        areas = node_control_areas_cm2(points, cells)
        ranked = sorted(range(count), key=lambda idx: abs(gava[idx]), reverse=True)[:top_n]
        for rank, idx in enumerate(ranked, start=1):
            x_um, y_um = points[idx]
            g_n = alpha_n[idx] * abs(jn[idx]) / Q_C
            g_p = alpha_p[idx] * abs(jp[idx]) / Q_C
            rows.append({
                "bias_V": bias_key,
                "source": "sentaurus_node",
                "rank": rank,
                "entity_id": idx,
                "x_um": x_um,
                "y_um": y_um,
                "window_label": window_label(x_um),
                "Fn_V_per_cm": fn[idx],
                "Fp_V_per_cm": fp[idx],
                "alpha_n_cm_inv": alpha_n[idx],
                "alpha_p_cm_inv": alpha_p[idx],
                "Jn_mag_A_per_cm2": abs(jn[idx]),
                "Jp_mag_A_per_cm2": abs(jp[idx]),
                "Gava_n_cm_minus3_s_minus1": g_n,
                "Gava_p_cm_minus3_s_minus1": g_p,
                "Gava_total_cm_minus3_s_minus1": gava[idx],
                "qG_contribution_A_per_um": Q_C * gava[idx] * (areas[idx] if idx < len(areas) else 0.0) * DEPTH_CM,
                "support_weight_cm2": areas[idx] if idx < len(areas) else 0.0,
            })
    return rows


def discover_vtk_for_bias(case_dir: Path) -> dict[float, Path]:
    return bv.discover_vtks(case_dir)


def vector_field(data: bv.VtkData, name: str) -> list[float]:
    return data.fields.get(name, [])


def load_self_vtk_hotspots(case_dir: Path, top_n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not case_dir.exists():
        return rows
    for bias, vtk_path in discover_vtk_for_bias(case_dir).items():
        data = bv.parse_vtk(vtk_path)
        points_um = [(x * 1.0e6, y * 1.0e6) for x, y in data.points]
        areas = node_control_areas_cm2(points_um, data.cells)
        gava_m3 = data.fields.get("AvalancheGeneration", [])
        alpha_n_m = data.fields.get("ElectronAlphaAvalanche", [])
        alpha_p_m = data.fields.get("HoleAlphaAvalanche", [])
        fn_m = data.fields.get("ElectronImpactIonizationDrive") or data.fields.get("ElectronHighFieldDrive", [])
        fp_m = data.fields.get("HoleImpactIonizationDrive") or data.fields.get("HoleHighFieldDrive", [])
        jn_m2 = vector_field(data, "ElectronCurrentDensityVector")
        jp_m2 = vector_field(data, "HoleCurrentDensityVector")
        ranked = sorted(range(min(len(points_um), len(gava_m3))), key=lambda idx: abs(gava_m3[idx]), reverse=True)[:top_n]
        for rank, idx in enumerate(ranked, start=1):
            x_um, y_um = points_um[idx]
            alpha_n = (alpha_n_m[idx] if idx < len(alpha_n_m) else 0.0) * 0.01
            alpha_p = (alpha_p_m[idx] if idx < len(alpha_p_m) else 0.0) * 0.01
            jn = (abs(jn_m2[idx]) if idx < len(jn_m2) else 0.0) * 1.0e-4
            jp = (abs(jp_m2[idx]) if idx < len(jp_m2) else 0.0) * 1.0e-4
            g_total = gava_m3[idx] * 1.0e-6
            rows.append({
                "bias_V": bias,
                "source": "self_vtk_node",
                "rank": rank,
                "entity_id": idx,
                "x_um": x_um,
                "y_um": y_um,
                "window_label": window_label(x_um),
                "Fn_V_per_cm": (fn_m[idx] if idx < len(fn_m) else 0.0) * 0.01,
                "Fp_V_per_cm": (fp_m[idx] if idx < len(fp_m) else 0.0) * 0.01,
                "alpha_n_cm_inv": alpha_n,
                "alpha_p_cm_inv": alpha_p,
                "Jn_mag_A_per_cm2": jn,
                "Jp_mag_A_per_cm2": jp,
                "Gava_n_cm_minus3_s_minus1": alpha_n * jn / Q_C,
                "Gava_p_cm_minus3_s_minus1": alpha_p * jp / Q_C,
                "Gava_total_cm_minus3_s_minus1": g_total,
                "qG_contribution_A_per_um": Q_C * g_total * (areas[idx] if idx < len(areas) else 0.0) * DEPTH_CM,
                "support_weight_cm2": areas[idx] if idx < len(areas) else 0.0,
            })
    return rows


def load_self_edge_hotspots(path: Path, top_n: int, biases: list[float]) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    wanted = {bv.bias_key(bias) for bias in biases}
    by_bias: dict[float, list[dict[str, str]]] = {}
    for row in bv.read_rows(path):
        bias = bv.finite_float(row.get("bias_V"))
        if bias is None:
            continue
        key = bv.bias_key(bias)
        if key in wanted:
            by_bias.setdefault(key, []).append(row)
    for bias, items in by_bias.items():
        ranked = sorted(
            items,
            key=lambda row: abs(bv.finite_float(row.get("edge_source_integral")) or 0.0) /
            max(abs(bv.finite_float(row.get("edge_area_proxy_m2")) or 0.0), 1.0e-300),
            reverse=True,
        )[:top_n]
        for rank, row in enumerate(ranked, start=1):
            area_m2 = bv.finite_float(row.get("edge_area_proxy_m2")) or 0.0
            x_um = 0.5 * ((bv.finite_float(row.get("x0_um")) or 0.0) + (bv.finite_float(row.get("x1_um")) or 0.0))
            y_um = 0.5 * ((bv.finite_float(row.get("y0_um")) or 0.0) + (bv.finite_float(row.get("y1_um")) or 0.0))
            en_source = bv.finite_float(row.get("electron_source_integral")) or 0.0
            hp_source = bv.finite_float(row.get("hole_source_integral")) or 0.0
            total_source = bv.finite_float(row.get("edge_source_integral")) or 0.0
            rows.append({
                "bias_V": bias,
                "source": "self_internal_edge",
                "rank": rank,
                "entity_id": row.get("edge_id", ""),
                "x_um": x_um,
                "y_um": y_um,
                "window_label": window_label(x_um),
                "Fn_V_per_cm": (bv.finite_float(row.get("electron_impact_field_V_per_m")) or 0.0) * 0.01,
                "Fp_V_per_cm": (bv.finite_float(row.get("hole_impact_field_V_per_m")) or 0.0) * 0.01,
                "alpha_n_cm_inv": (bv.finite_float(row.get("electron_alpha_m_inv")) or 0.0) * 0.01,
                "alpha_p_cm_inv": (bv.finite_float(row.get("hole_alpha_m_inv")) or 0.0) * 0.01,
                "Jn_mag_A_per_cm2": (bv.finite_float(row.get("electron_flux_proxy")) or 0.0) * 1.0e-4,
                "Jp_mag_A_per_cm2": (bv.finite_float(row.get("hole_flux_proxy")) or 0.0) * 1.0e-4,
                "Gava_n_cm_minus3_s_minus1": en_source / max(area_m2, 1.0e-300) * 1.0e-6,
                "Gava_p_cm_minus3_s_minus1": hp_source / max(area_m2, 1.0e-300) * 1.0e-6,
                "Gava_total_cm_minus3_s_minus1": total_source / max(area_m2, 1.0e-300) * 1.0e-6,
                "qG_contribution_A_per_um": Q_C * total_source * 1.0e-6,
                "support_weight_cm2": area_m2 * 1.0e4,
            })
    return rows


def load_self_dual_cell_hotspots(path: Path, top_n: int, biases: list[float]) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    wanted = {bv.bias_key(bias) for bias in biases}
    by_bias: dict[float, list[dict[str, str]]] = {}
    for row in bv.read_rows(path):
        bias = bv.finite_float(row.get("bias_V"))
        if bias is None:
            continue
        key = bv.bias_key(bias)
        if key in wanted and row.get("reconstruction_status", "ok") == "ok":
            by_bias.setdefault(key, []).append(row)
    rows: list[dict[str, Any]] = []
    for bias, items in by_bias.items():
        ranked = sorted(
            items,
            key=lambda row: abs(bv.finite_float(row.get("Gava_cell_cm_minus3_s_minus1")) or 0.0),
            reverse=True,
        )[:top_n]
        for rank, row in enumerate(ranked, start=1):
            jn = bv.finite_float(row.get("Jn_mag_A_per_cm2")) or 0.0
            jp = bv.finite_float(row.get("Jp_mag_A_per_cm2")) or 0.0
            alpha_n = bv.finite_float(row.get("alpha_n_cell_cm_inv")) or 0.0
            alpha_p = bv.finite_float(row.get("alpha_p_cell_cm_inv")) or 0.0
            x_um = bv.finite_float(row.get("centroid_x_um")) or 0.0
            y_um = bv.finite_float(row.get("centroid_y_um")) or 0.0
            rows.append({
                "bias_V": bias,
                "source": "self_dual_cell",
                "rank": rank,
                "entity_id": row.get("cell_id", ""),
                "x_um": x_um,
                "y_um": y_um,
                "window_label": row.get("window_label") or window_label(x_um),
                "Fn_V_per_cm": bv.finite_float(row.get("Fn_cell_V_per_cm")),
                "Fp_V_per_cm": bv.finite_float(row.get("Fp_cell_V_per_cm")),
                "alpha_n_cm_inv": alpha_n,
                "alpha_p_cm_inv": alpha_p,
                "Jn_mag_A_per_cm2": jn,
                "Jp_mag_A_per_cm2": jp,
                "Gava_n_cm_minus3_s_minus1": alpha_n * abs(jn) / Q_C,
                "Gava_p_cm_minus3_s_minus1": alpha_p * abs(jp) / Q_C,
                "Gava_total_cm_minus3_s_minus1": bv.finite_float(row.get("Gava_cell_cm_minus3_s_minus1")),
                "qG_contribution_A_per_um": bv.finite_float(row.get("qG_cell_contribution_A_per_um")),
                "support_weight_cm2": safe_ratio(
                    bv.finite_float(row.get("qG_cell_contribution_A_per_um")),
                    Q_C * (bv.finite_float(row.get("Gava_cell_cm_minus3_s_minus1")) or 0.0) * DEPTH_CM,
                ),
            })
    return rows


def add_nearest_counterparts(rows: list[dict[str, Any]]) -> None:
    by_bias: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        bias = bv.finite_float(row.get("bias_V"))
        if bias is not None:
            by_bias.setdefault(bv.bias_key(bias), []).append(row)
    for subset in by_bias.values():
        sentaurus = [row for row in subset if row.get("source") == "sentaurus_node"]
        self_rows = [row for row in subset if str(row.get("source", "")).startswith("self_")]
        for row in subset:
            candidates = self_rows if row.get("source") == "sentaurus_node" else sentaurus
            best = min(
                ((distance_um(row, cand), cand) for cand in candidates),
                key=lambda item: math.inf if item[0] is None else item[0],
                default=(None, None),
            )
            row["nearest_counterpart_source"] = "" if best[1] is None else best[1].get("source", "")
            row["nearest_counterpart_rank"] = "" if best[1] is None else best[1].get("rank", "")
            row["nearest_counterpart_distance_um"] = best[0]


def output_fields() -> list[str]:
    return [
        "bias_V",
        "source",
        "rank",
        "entity_id",
        "x_um",
        "y_um",
        "window_label",
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
        "support_weight_cm2",
        "nearest_counterpart_source",
        "nearest_counterpart_rank",
        "nearest_counterpart_distance_um",
    ]


def top_row(rows: list[dict[str, Any]], bias: float, source: str, metric: str = "Gava_total_cm_minus3_s_minus1") -> dict[str, Any] | None:
    candidates = [
        row for row in rows
        if bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0) == bv.bias_key(bias) and row.get("source") == source
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda row: abs(bv.finite_float(row.get(metric)) or 0.0))


def max_metric_row(rows: list[dict[str, Any]], bias: float, source: str, metrics: list[str]) -> dict[str, Any] | None:
    candidates = [
        row for row in rows
        if bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0) == bv.bias_key(bias) and row.get("source") == source
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda row: max(abs(bv.finite_float(row.get(metric)) or 0.0) for metric in metrics))


def describe_row(row: dict[str, Any] | None) -> str:
    if row is None:
        return "unavailable"
    return (
        f"{row.get('source')} rank {row.get('rank')} at "
        f"({fmt(row.get('x_um'))},{fmt(row.get('y_um'))}) um, "
        f"window={row.get('window_label')}, G={fmt(row.get('Gava_total_cm_minus3_s_minus1'))}"
    )


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    sent20 = top_row(rows, -20.0, "sentaurus_node")
    self20 = top_row(rows, -20.0, "self_dual_cell") or top_row(rows, -20.0, "self_internal_edge") or top_row(rows, -20.0, "self_vtk_node")
    self_source = "" if self20 is None else str(self20.get("source"))
    sent_f_alpha = max_metric_row(rows, -20.0, "sentaurus_node", ["Fn_V_per_cm", "Fp_V_per_cm", "alpha_n_cm_inv", "alpha_p_cm_inv"])
    self_f_alpha = max_metric_row(rows, -20.0, self_source, ["Fn_V_per_cm", "Fp_V_per_cm", "alpha_n_cm_inv", "alpha_p_cm_inv"]) if self_source else None
    sent_j = max_metric_row(rows, -20.0, "sentaurus_node", ["Jn_mag_A_per_cm2", "Jp_mag_A_per_cm2"])
    self_j = max_metric_row(rows, -20.0, self_source, ["Jn_mag_A_per_cm2", "Jp_mag_A_per_cm2"]) if self_source else None
    self_qg = top_row(rows, -20.0, self_source, "qG_contribution_A_per_um") if self_source else None
    self_to_sent = distance_um(self20, sent20) if self20 and sent20 else None
    self_to_f = distance_um(self20, self_f_alpha) if self20 and self_f_alpha else None
    self_to_j = distance_um(self20, self_j) if self20 and self_j else None
    self_to_qg = distance_um(self20, self_qg) if self20 and self_qg else None
    lines = [
        "# Gava Hotspot Decomposition",
        "",
        "- Rows list top Gava hotspots per bias for Sentaurus and available self sources.",
        "- Units: F in V/cm, alpha in cm^-1, J in A/cm^2, Gava in cm^-3 s^-1, qG in A/um.",
        "",
        "## -20 V Focus",
        "",
        f"- Sentaurus max Gava: {describe_row(sent20)}",
        f"- Self max Gava: {describe_row(self20)}",
        f"- Self-to-Sentaurus max distance: {fmt(self_to_sent)} um",
        "",
        "## Answers",
        "",
        "1. Is self hotspot shifted because F/alpha peak shifts? Within the listed top-Gava hotspots, "
        f"Self F/alpha peak is {describe_row(self_f_alpha)}; distance from self Gava max is {fmt(self_to_f)} um. "
        f"Sentaurus F/alpha peak is {describe_row(sent_f_alpha)}.",
        "2. Is it shifted because J peak shifts? Within the listed top-Gava hotspots, "
        f"Self J peak is {describe_row(self_j)}; distance from self Gava max is {fmt(self_to_j)} um. "
        f"Sentaurus J peak is {describe_row(sent_j)}.",
        "3. Is it shifted because source weight/mapping differs? Within the listed top-Gava hotspots, "
        f"Self qG-contribution peak is {describe_row(self_qg)}; distance from self Gava max is {fmt(self_to_qg)} um. "
        "If this distance is small while Sentaurus max is elsewhere, the mapping/weight is not moving the self maximum by itself; the local J/F/alpha product is.",
        "4. At -20V, why does max remain at (1.08333,0.416667) um? "
        f"The selected self source is {self_source}; its max row is {describe_row(self20)} with "
        f"Fn={fmt(None if self20 is None else self20.get('Fn_V_per_cm'))}, "
        f"Fp={fmt(None if self20 is None else self20.get('Fp_V_per_cm'))}, "
        f"alpha_n={fmt(None if self20 is None else self20.get('alpha_n_cm_inv'))}, "
        f"alpha_p={fmt(None if self20 is None else self20.get('alpha_p_cm_inv'))}, "
        f"Jn={fmt(None if self20 is None else self20.get('Jn_mag_A_per_cm2'))}, "
        f"Jp={fmt(None if self20 is None else self20.get('Jp_mag_A_per_cm2'))}.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'gava_hotspot_decomposition.csv'}",
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
        "--self-case-dir",
        type=Path,
        default=REPO / "build/diagnostics/avalanche_current_magnitude_mode_coupled_matrix_cases/ava-current-edge-scalar",
    )
    parser.add_argument(
        "--self-edge-csv",
        type=Path,
        default=REPO / "build/diagnostics/avalanche_current_magnitude_mode_coupled_matrix_cases/ava-current-edge-scalar/ava-current-edge-scalar_sg_avalanche_edges.csv",
    )
    parser.add_argument(
        "--self-dual-cell-csv",
        type=Path,
        default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv",
    )
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/gava_hotspot_decomposition.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/gava_hotspot_decomposition_summary.md")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = bv.parse_biases(args.biases)
    rows: list[dict[str, Any]] = []
    rows.extend(load_sentaurus_hotspots(args.sentaurus_multibias_dir, biases, args.top_n))
    rows.extend(load_self_vtk_hotspots(args.self_case_dir, args.top_n))
    rows.extend(load_self_edge_hotspots(args.self_edge_csv, args.top_n, biases))
    rows.extend(load_self_dual_cell_hotspots(args.self_dual_cell_csv, args.top_n, biases))
    rows.sort(key=lambda row: (
        bv.finite_float(row.get("bias_V")) or 0.0,
        str(row.get("source", "")),
        int(row.get("rank", 0) or 0),
    ))
    add_nearest_counterparts(rows)
    bv.write_rows(args.out_csv, rows, output_fields())
    write_summary(args.out_summary, rows)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
