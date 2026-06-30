#!/usr/bin/env python3
"""Reconstruct cell vector current from internal edge scalar current diagnostics."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19
WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}
QF_DENSITY_MEANS = ("arithmetic", "geometric", "logarithmic", "upwind")
BASE_MODES = (
    "current_edge_scalar_existing",
    "current_cell_vector_ls_from_edge_tangent",
    "current_cell_vector_ls_from_box_normal",
)


@dataclass(frozen=True)
class Edge:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    length_cm: float
    raw_e_sign: float = 1.0
    raw_h_sign: float = 1.0

    @property
    def key(self) -> tuple[int, int]:
        return tuple(sorted((self.node0, self.node1)))

    @property
    def tangent(self) -> tuple[float, float]:
        dx = self.x1_um - self.x0_um
        dy = self.y1_um - self.y0_um
        length = math.hypot(dx, dy)
        if length <= 0.0:
            return (0.0, 0.0)
        return (dx / length, dy / length)

    @property
    def normal(self) -> tuple[float, float]:
        tx, ty = self.tangent
        return (-ty, tx)


@dataclass(frozen=True)
class Element:
    cell_id: int
    nodes: tuple[int, int, int]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/avalanche_internal_source_current_audit_default.csv")
    parser.add_argument(
        "--self-edge-topology",
        type=Path,
        default=(
            REPO
            / "build/diagnostics/avalanche_internal_source_current_audit_case"
            / "BV-A1-B1p00-internal-source-current-audit"
            / "BV-A1-B1p00-internal-source-current-audit_sg_avalanche_edges.csv"
        ),
    )
    parser.add_argument(
        "--elements-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_fields_20260629_155228/elements.csv",
    )
    parser.add_argument(
        "--node-compare",
        type=Path,
        default=(
            REPO
            / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
            / "coarse_previous_full20_fields_20260629_155228/coarse_node_field_compare_aligned.csv"
        ),
    )
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/edge_to_cell_vector_current_reconstruction.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/edge_to_cell_vector_current_reconstruction_summary.md")
    return parser.parse_args()


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items() if key is not None}


def f(raw: Any, default: float = 0.0) -> float:
    try:
        value = float(str(raw).strip())
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def safe_ratio(value: float | None, reference: float | None) -> float | None:
    if value is None or reference is None or abs(reference) <= 1.0e-300:
        return None
    ratio = value / reference
    return ratio if math.isfinite(ratio) else None


def parse_biases(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def nearest_bias(value: float, targets: list[float], tol: float) -> float | None:
    best = min(targets, key=lambda item: abs(item - value))
    return best if abs(best - value) <= tol else None


def sign_from(value: float) -> float:
    if value > 0.0:
        return 1.0
    if value < 0.0:
        return -1.0
    return 1.0


def load_edges(path: Path) -> dict[int, Edge]:
    edges: dict[int, Edge] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            edge_id = int(f(row["edge_id"]))
            if edge_id in edges:
                continue
            edges[edge_id] = Edge(
                edge_id=edge_id,
                node0=int(f(row["node0"])),
                node1=int(f(row["node1"])),
                x0_um=f(row["x0_um"]),
                y0_um=f(row["y0_um"]),
                x1_um=f(row["x1_um"]),
                y1_um=f(row["y1_um"]),
                length_cm=f(row["edge_length_m"]) * 100.0,
                raw_e_sign=sign_from(f(row.get("electron_raw_signed_flux_proxy", "1"))),
                raw_h_sign=sign_from(f(row.get("hole_raw_signed_flux_proxy", "1"))),
            )
    return edges


def load_elements(path: Path) -> list[Element]:
    elements: list[Element] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            elements.append(Element(
                cell_id=int(f(row.get("id", len(elements)))),
                nodes=(int(f(row["node0"])), int(f(row["node1"])), int(f(row["node2"]))),
            ))
    return elements


def load_node_compare(path: Path) -> dict[tuple[float, int], dict[str, dict[str, float | None]]]:
    nodes: dict[tuple[float, int], dict[str, dict[str, float | None]]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            bias = f(row.get("bias_V"))
            node_id = int(f(row.get("node_id")))
            quantity = row.get("quantity", "")
            entry = nodes.setdefault((bias, node_id), {})
            entry[quantity] = {
                "sent": f(row.get("sentaurus_value"), math.nan),
                "self": f(row.get("vela_value_scaled_to_sentaurus_units"), math.nan),
                "x_um": f(row.get("x_um"), math.nan),
                "y_um": f(row.get("y_um"), math.nan),
            }
    return nodes


def node_value(nodes: dict[tuple[float, int], dict[str, dict[str, float | None]]], bias: float, node_id: int, quantity: str, side: str) -> float:
    value = nodes.get((bias, node_id), {}).get(quantity, {}).get(side)
    return 0.0 if value is None or not math.isfinite(value) else float(value)


def node_xy(nodes: dict[tuple[float, int], dict[str, dict[str, float | None]]], bias: float, node_id: int) -> tuple[float, float]:
    entry = nodes.get((bias, node_id), {})
    for quantity in entry.values():
        x = quantity.get("x_um")
        y = quantity.get("y_um")
        if x is not None and y is not None and math.isfinite(x) and math.isfinite(y):
            return float(x), float(y)
    return 0.0, 0.0


def triangle_area_um2(coords: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = coords
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))


def triangle_grad_v_per_cm(coords: list[tuple[float, float]], values: list[float]) -> tuple[float, float]:
    (x1, y1), (x2, y2), (x3, y3) = coords
    v1, v2, v3 = values
    double_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if abs(double_area) <= 1.0e-300:
        return 0.0, 0.0
    grad_x_v_per_um = (v1 * (y2 - y3) + v2 * (y3 - y1) + v3 * (y1 - y2)) / double_area
    grad_y_v_per_um = (v1 * (x3 - x2) + v2 * (x1 - x3) + v3 * (x2 - x1)) / double_area
    return grad_x_v_per_um * 1.0e4, grad_y_v_per_um * 1.0e4


def ls_vector(rows: list[tuple[tuple[float, float], float]]) -> tuple[float, float]:
    a11 = a12 = a22 = b1 = b2 = 0.0
    for (ux, uy), value in rows:
        a11 += ux * ux
        a12 += ux * uy
        a22 += uy * uy
        b1 += ux * value
        b2 += uy * value
    det = a11 * a22 - a12 * a12
    if abs(det) <= 1.0e-300:
        return 0.0, 0.0
    return (b1 * a22 - b2 * a12) / det, (a11 * b2 - a12 * b1) / det


def geometric(values: list[float]) -> float:
    clean = [max(value, 0.0) for value in values]
    product = 1.0
    for value in clean:
        product *= value
    return product ** (1.0 / len(clean))


def logarithmic_pair(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.5 * (a + b)
    if abs(a - b) <= 1.0e-14 * max(abs(a), abs(b), 1.0):
        return 0.5 * (a + b)
    return (b - a) / (math.log(b) - math.log(a))


def density_mean(values: list[float], mode: str, grad: tuple[float, float], coords: list[tuple[float, float]]) -> float:
    if mode == "arithmetic":
        return sum(values) / len(values)
    if mode == "geometric":
        return geometric(values)
    if mode == "logarithmic":
        return sum(logarithmic_pair(values[i], values[(i + 1) % len(values)]) for i in range(len(values))) / len(values)
    gx, gy = grad
    if mode == "upwind":
        cx = sum(x for x, _ in coords) / 3.0
        cy = sum(y for _, y in coords) / 3.0
        scores = [((x - cx) * gx + (y - cy) * gy, idx) for idx, (x, y) in enumerate(coords)]
        return values[min(scores, key=lambda item: item[0])[1]]
    return sum(values) / len(values)


def load_internal_rows(path: Path, requested_biases: list[float], tol: float) -> dict[tuple[float, int], dict[str, float]]:
    rows: dict[tuple[float, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            if row.get("source_location_type", "edge") != "edge":
                continue
            bias = nearest_bias(f(row.get("bias_V")), requested_biases, tol)
            if bias is None:
                continue
            edge_id = int(f(row.get("source_entity_id")))
            rows[(bias, edge_id)] = {
                "Fn": f(row.get("Fn_used_V_per_cm")),
                "Fp": f(row.get("Fp_used_V_per_cm")),
                "alpha_n": f(row.get("alpha_n_used_cm_inv")),
                "alpha_p": f(row.get("alpha_p_used_cm_inv")),
                "Jn": f(row.get("Jn_mag_used_A_per_cm2")),
                "Jp": f(row.get("Jp_mag_used_A_per_cm2")),
                "Gava": f(row.get("Gava_total_used_cm_minus3_s_minus1")),
                "area_cm2": f(row.get("source_weight_or_volume_cm2_for_2D") or row.get("contribution_volume_cm3_or_area_cm2_for_2D")),
                "qG": f(row.get("qG_contribution_A_per_um")),
            }
    return rows


def build_edge_cell_counts(elements: list[Element]) -> dict[tuple[int, int], int]:
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for elem in elements:
        a, b, c = elem.nodes
        for pair in ((a, b), (b, c), (c, a)):
            counts[tuple(sorted(pair))] += 1
    return counts


def element_edges(elem: Element, edge_by_key: dict[tuple[int, int], Edge]) -> list[Edge]:
    a, b, c = elem.nodes
    return [
        edge_by_key[key]
        for key in (tuple(sorted(pair)) for pair in ((a, b), (b, c), (c, a)))
        if key in edge_by_key
    ]


def cell_sentaurus_values(nodes: dict[tuple[float, int], dict[str, dict[str, float | None]]], bias: float, elem: Element) -> dict[str, float]:
    nids = elem.nodes
    def avg(quantity: str, side: str = "sent") -> float:
        return sum(node_value(nodes, bias, node, quantity, side) for node in nids) / 3.0
    jnx = avg("electron_current_density_x")
    jny = avg("electron_current_density_y")
    jpx = avg("hole_current_density_x")
    jpy = avg("hole_current_density_y")
    return {
        "Jn_sent_cell_vector_mag_A_per_cm2": math.hypot(jnx, jny),
        "Jp_sent_cell_vector_mag_A_per_cm2": math.hypot(jpx, jpy),
        "alpha_n_sent_cell_cm_inv": avg("electron_alpha_avalanche"),
        "alpha_p_sent_cell_cm_inv": avg("hole_alpha_avalanche"),
        "Gava_sent_cell_cm_minus3_s_minus1": avg("avalanche"),
    }


def sentaurus_cell_edge_projection_mean(
    nodes: dict[tuple[float, int], dict[str, dict[str, float | None]]],
    bias: float,
    cell_edges: list[Edge],
    prefix: str,
) -> float:
    values: list[float] = []
    for edge in cell_edges:
        jx0 = node_value(nodes, bias, edge.node0, f"{prefix}_x", "sent")
        jy0 = node_value(nodes, bias, edge.node0, f"{prefix}_y", "sent")
        jx1 = node_value(nodes, bias, edge.node1, f"{prefix}_x", "sent")
        jy1 = node_value(nodes, bias, edge.node1, f"{prefix}_y", "sent")
        tx, ty = edge.tangent
        values.append(abs((0.5 * (jx0 + jx1)) * tx + (0.5 * (jy0 + jy1)) * ty))
    return sum(values) / len(values) if values else 0.0
def make_output_row(
    *,
    bias: float,
    elem: Element,
    coords: list[tuple[float, float]],
    mode: str,
    jn_vec: tuple[float, float],
    jp_vec: tuple[float, float],
    alpha_n: float,
    alpha_p: float,
    fn: float,
    fp: float,
    area_cm2: float,
    sent: dict[str, float],
    baseline_qg: float | None = None,
) -> dict[str, Any]:
    jn_mag = math.hypot(*jn_vec)
    jp_mag = math.hypot(*jp_vec)
    gava = (alpha_n * jn_mag + alpha_p * jp_mag) / Q_C
    qg = baseline_qg if baseline_qg is not None else Q_C * gava * area_cm2 / 1.0e4
    sent_qg = Q_C * sent["Gava_sent_cell_cm_minus3_s_minus1"] * area_cm2 / 1.0e4
    cx = sum(x for x, _ in coords) / 3.0
    cy = sum(y for _, y in coords) / 3.0
    return {
        "bias_V": bias,
        "mode": mode,
        "cell_id": elem.cell_id,
        "centroid_x_um": cx,
        "centroid_y_um": cy,
        "Jn_x_A_per_cm2": jn_vec[0],
        "Jn_y_A_per_cm2": jn_vec[1],
        "Jn_mag_A_per_cm2": jn_mag,
        "Jp_x_A_per_cm2": jp_vec[0],
        "Jp_y_A_per_cm2": jp_vec[1],
        "Jp_mag_A_per_cm2": jp_mag,
        "alpha_n_cell_cm_inv": alpha_n,
        "alpha_p_cell_cm_inv": alpha_p,
        "Fn_cell_V_per_cm": fn,
        "Fp_cell_V_per_cm": fp,
        "Gava_cell_cm_minus3_s_minus1": gava,
        "qG_cell_contribution_A_per_um": qg,
        **sent,
        "qG_sent_cell_contribution_A_per_um": sent_qg,
        "qG_ratio_to_sentaurus": safe_ratio(qg, sent_qg),
        "window_label": window_label(cx),
    }


def window_label(x_um: float) -> str:
    for name in ("left_shoulder", "center", "right_shoulder"):
        lo, hi = WINDOWS[name]
        if lo <= x_um <= hi:
            return name
    lo, hi = WINDOWS["junction"]
    if lo <= x_um <= hi:
        return "junction"
    return "outside_junction"


def window_members(x_um: float) -> list[str]:
    return [name for name, (lo, hi) in WINDOWS.items() if lo <= x_um <= hi]


def run(args: argparse.Namespace) -> None:
    requested_biases = parse_biases(args.biases)
    edges = load_edges(args.self_edge_topology)
    elements = load_elements(args.elements_csv)
    nodes = load_node_compare(args.node_compare)
    internal = load_internal_rows(args.internal_audit, requested_biases, args.bias_tol)
    edge_by_key = {edge.key: edge for edge in edges.values()}
    edge_counts = build_edge_cell_counts(elements)
    rows: list[dict[str, Any]] = []

    for bias in requested_biases:
        for elem in elements:
            coords = [node_xy(nodes, bias, node) for node in elem.nodes]
            area_cm2 = triangle_area_um2(coords) * 1.0e-8
            if area_cm2 <= 0.0:
                continue
            cell_edges = element_edges(elem, edge_by_key)
            sent = cell_sentaurus_values(nodes, bias, elem)

            baseline_qg = 0.0
            weighted_jn = weighted_jp = weighted_fn = weighted_fp = weighted_an = weighted_ap = weight_sum = 0.0
            tangent_e: list[tuple[tuple[float, float], float]] = []
            tangent_h: list[tuple[tuple[float, float], float]] = []
            normal_e: list[tuple[tuple[float, float], float]] = []
            normal_h: list[tuple[tuple[float, float], float]] = []
            for edge in cell_edges:
                edge_row = internal.get((bias, edge.edge_id))
                if edge_row is None:
                    continue
                share = max(edge_counts.get(edge.key, 1), 1)
                baseline_qg += edge_row["qG"] / share
                weight = max(edge.length_cm, 1.0e-300)
                weighted_jn += edge_row["Jn"] * weight
                weighted_jp += edge_row["Jp"] * weight
                weighted_fn += edge_row["Fn"] * weight
                weighted_fp += edge_row["Fp"] * weight
                weighted_an += edge_row["alpha_n"] * weight
                weighted_ap += edge_row["alpha_p"] * weight
                weight_sum += weight
                tx, ty = edge.tangent
                nx, ny = edge.normal
                tangent_e.append(((tx, ty), edge.raw_e_sign * edge_row["Jn"]))
                tangent_h.append(((tx, ty), edge.raw_h_sign * edge_row["Jp"]))
                normal_e.append(((nx, ny), edge.raw_e_sign * edge_row["Jn"]))
                normal_h.append(((nx, ny), edge.raw_h_sign * edge_row["Jp"]))

            sent["Jn_sent_cell_edge_projection_mean_A_per_cm2"] = sentaurus_cell_edge_projection_mean(nodes, bias, cell_edges, "electron_current_density")
            sent["Jp_sent_cell_edge_projection_mean_A_per_cm2"] = sentaurus_cell_edge_projection_mean(nodes, bias, cell_edges, "hole_current_density")
            sent["Jn_sent_vector_mag_over_edge_projection"] = safe_ratio(sent["Jn_sent_cell_vector_mag_A_per_cm2"], sent["Jn_sent_cell_edge_projection_mean_A_per_cm2"])
            sent["Jp_sent_vector_mag_over_edge_projection"] = safe_ratio(sent["Jp_sent_cell_vector_mag_A_per_cm2"], sent["Jp_sent_cell_edge_projection_mean_A_per_cm2"])

            if weight_sum <= 0.0:
                continue
            avg_jn = weighted_jn / weight_sum
            avg_jp = weighted_jp / weight_sum
            alpha_n = weighted_an / weight_sum
            alpha_p = weighted_ap / weight_sum
            fn = weighted_fn / weight_sum
            fp = weighted_fp / weight_sum
            rows.append(make_output_row(
                bias=bias, elem=elem, coords=coords, mode="current_edge_scalar_existing",
                jn_vec=(avg_jn, 0.0), jp_vec=(avg_jp, 0.0),
                alpha_n=alpha_n, alpha_p=alpha_p, fn=fn, fp=fp,
                area_cm2=area_cm2, sent=sent, baseline_qg=baseline_qg,
            ))
            rows.append(make_output_row(
                bias=bias, elem=elem, coords=coords, mode="current_cell_vector_ls_from_edge_tangent",
                jn_vec=ls_vector(tangent_e), jp_vec=ls_vector(tangent_h),
                alpha_n=alpha_n, alpha_p=alpha_p, fn=fn, fp=fp,
                area_cm2=area_cm2, sent=sent,
            ))
            rows.append(make_output_row(
                bias=bias, elem=elem, coords=coords, mode="current_cell_vector_ls_from_box_normal",
                jn_vec=ls_vector(normal_e), jp_vec=ls_vector(normal_h),
                alpha_n=alpha_n, alpha_p=alpha_p, fn=fn, fp=fp,
                area_cm2=area_cm2, sent=sent,
            ))

            phin = [node_value(nodes, bias, node, "electron_qf", "self") for node in elem.nodes]
            phip = [node_value(nodes, bias, node, "hole_qf", "self") for node in elem.nodes]
            grad_n = triangle_grad_v_per_cm(coords, phin)
            grad_p = triangle_grad_v_per_cm(coords, phip)
            fn_cell = math.hypot(*grad_n)
            fp_cell = math.hypot(*grad_p)
            n_vals = [node_value(nodes, bias, node, "electron_density", "self") for node in elem.nodes]
            p_vals = [node_value(nodes, bias, node, "hole_density", "self") for node in elem.nodes]
            mu_n_vals = [node_value(nodes, bias, node, "electron_mobility", "self") for node in elem.nodes]
            mu_p_vals = [node_value(nodes, bias, node, "hole_mobility", "self") for node in elem.nodes]
            mu_n = sum(mu_n_vals) / 3.0
            mu_p = sum(mu_p_vals) / 3.0
            for mean in QF_DENSITY_MEANS:
                n_eff = density_mean(n_vals, mean, grad_n, coords)
                p_eff = density_mean(p_vals, mean, grad_p, coords)
                jn_vec = (-Q_C * mu_n * n_eff * grad_n[0], -Q_C * mu_n * n_eff * grad_n[1])
                jp_vec = (-Q_C * mu_p * p_eff * grad_p[0], -Q_C * mu_p * p_eff * grad_p[1])
                rows.append(make_output_row(
                    bias=bias, elem=elem, coords=coords,
                    mode=f"current_cell_vector_from_qf_gradient_{mean}",
                    jn_vec=jn_vec, jp_vec=jp_vec,
                    alpha_n=alpha_n, alpha_p=alpha_p, fn=fn_cell, fp=fp_cell,
                    area_cm2=area_cm2, sent=sent,
                ))

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "bias_V", "mode", "cell_id", "centroid_x_um", "centroid_y_um",
        "Jn_x_A_per_cm2", "Jn_y_A_per_cm2", "Jn_mag_A_per_cm2",
        "Jp_x_A_per_cm2", "Jp_y_A_per_cm2", "Jp_mag_A_per_cm2",
        "alpha_n_cell_cm_inv", "alpha_p_cell_cm_inv", "Fn_cell_V_per_cm", "Fp_cell_V_per_cm",
        "Gava_cell_cm_minus3_s_minus1", "qG_cell_contribution_A_per_um",
        "Jn_sent_cell_vector_mag_A_per_cm2", "Jp_sent_cell_vector_mag_A_per_cm2",
        "Jn_sent_cell_edge_projection_mean_A_per_cm2", "Jp_sent_cell_edge_projection_mean_A_per_cm2",
        "Jn_sent_vector_mag_over_edge_projection", "Jp_sent_vector_mag_over_edge_projection",
        "alpha_n_sent_cell_cm_inv", "alpha_p_sent_cell_cm_inv",
        "Gava_sent_cell_cm_minus3_s_minus1", "qG_sent_cell_contribution_A_per_um",
        "qG_ratio_to_sentaurus", "window_label",
    ]
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_value(row.get(key)) for key in header})
    write_summary(args.out_summary, rows, internal)


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.17g}"
    return str(value)


def write_summary(path: Path, rows: list[dict[str, Any]], internal: dict[tuple[float, int], dict[str, float]]) -> None:
    by_mode_window: dict[tuple[float, str, str], dict[str, float]] = defaultdict(lambda: {"self": 0.0, "sent": 0.0})
    max_by_mode: dict[tuple[float, str], tuple[float, float, float]] = {}
    for row in rows:
        bias = float(row["bias_V"])
        mode = str(row["mode"])
        x = float(row["centroid_x_um"])
        for window in window_members(x):
            bucket = by_mode_window[(bias, mode, window)]
            bucket["self"] += float(row["qG_cell_contribution_A_per_um"])
            bucket["sent"] += float(row["qG_sent_cell_contribution_A_per_um"])
        key = (bias, mode)
        g = abs(float(row["Gava_cell_cm_minus3_s_minus1"]))
        if key not in max_by_mode or g > max_by_mode[key][0]:
            max_by_mode[key] = (g, float(row["centroid_x_um"]), float(row["centroid_y_um"]))

    internal_qg = sum(value["qG"] for value in internal.values())
    baseline_qg = sum(float(row["qG_cell_contribution_A_per_um"]) for row in rows if row["mode"] == "current_edge_scalar_existing")
    baseline_ok = abs(baseline_qg - internal_qg) / max(abs(internal_qg), 1.0e-300) <= 1.0e-10

    def ratio_for(bias: float, mode: str, window: str) -> float | None:
        bucket = by_mode_window.get((bias, mode, window))
        if not bucket:
            return None
        return safe_ratio(bucket["self"], bucket["sent"])

    def score(mode: str, window: str) -> float | None:
        values = []
        for bias in (-20.0,):
            ratio = ratio_for(bias, mode, window)
            if ratio is not None and ratio > 0.0:
                values.append(abs(math.log10(abs(ratio))))
        return median(values) if values else None

    modes = sorted({str(row["mode"]) for row in rows})
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Edge-to-Cell Vector Current Reconstruction\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver default path changed: no\n")
        out.write("- physics: VanOverstraeten/de Man, GradQuasiFermi, parameter_set=default\n")
        out.write("- A/B source: Sentaurus par/default; no A/B tuning in this diagnostic\n")
        out.write("- Sentaurus current reference includes both cell vector magnitude and cell edge-projection mean\n")
        out.write("- current-density unit: A/cm^2\n")
        out.write("- alpha unit: cm^-1\n")
        out.write("- Gava unit: cm^-3 s^-1\n")
        out.write("- qG unit: A/um\n")
        out.write(f"- edge scalar baseline reproduces current internal qG: {'yes' if baseline_ok else 'no'}\n\n")
        out.write("## -20V Window Ratios\n\n")
        out.write("| mode | full | junction | center | left_shoulder | right_shoulder |\n")
        out.write("| --- | --- | --- | --- | --- | --- |\n")
        for mode in modes:
            values = [ratio_for(-20.0, mode, window) for window in ("full", "junction", "center", "left_shoulder", "right_shoulder")]
            out.write("| " + mode + " | " + " | ".join("" if value is None else f"{value:.6g}" for value in values) + " |\n")
        out.write("\n## Required Answers\n\n")
        for window in ("full", "junction", "center", "left_shoulder", "right_shoulder"):
            ranked = [(mode, score(mode, window)) for mode in modes]
            ranked = [(mode, value) for mode, value in ranked if value is not None]
            best = min(ranked, key=lambda item: item[1])[0] if ranked else "insufficient_data"
            out.write(f"1.{window} best qG mode: {best}.\n")
        base = ratio_for(-20.0, "current_edge_scalar_existing", "full")
        best_full = min(
            [(mode, score(mode, "full")) for mode in modes if score(mode, "full") is not None],
            key=lambda item: item[1],
            default=("insufficient_data", None),
        )[0]
        best_ratio = ratio_for(-20.0, best_full, "full") if best_full != "insufficient_data" else None
        out.write(f"2. Cell vector |J| changes -20V full qG ratio from {format_ratio(base)} to {format_ratio(best_ratio)} for best mode {best_full}.\n")
        base_max = max_by_mode.get((-20.0, "current_edge_scalar_existing"))
        best_max = max_by_mode.get((-20.0, best_full))
        out.write(f"3. max Gava position baseline={format_pos(base_max)}, best={format_pos(best_max)}.\n")
        out.write(f"4. edge scalar baseline vs best full-window qG ratio delta: {format_delta(base, best_ratio)}.\n")
        t = score("current_cell_vector_ls_from_edge_tangent", "full")
        n = score("current_cell_vector_ls_from_box_normal", "full")
        if t is not None and n is not None:
            better = "tangent" if t < n else "box-normal"
        else:
            better = "insufficient_data"
        out.write(f"5. Tangent projection vs box-normal projection: {better} is closer by -20V full-window qG score.\n")
        if best_full not in {"insufficient_data", "current_edge_scalar_existing"}:
            out.write("6. Cell/vector reconstruction improves or changes the qG comparison; consider a guarded optional source mode using reconstructed vector magnitude.\n")
        else:
            out.write("6. Cell/vector reconstruction does not beat the edge scalar baseline; do not change source current magnitude from this evidence alone.\n")
        out.write("7. If no vector mode is satisfactory, continue SG formula/density averaging/sign convention investigation.\n")
        out.write("8. Sentaurus comparison distinguishes cell vector magnitude from edge projection via Jn/Jp_sent_vector_mag_over_edge_projection columns.\n")


def format_ratio(value: float | None) -> str:
    return "insufficient_data" if value is None else f"{value:.6g}"


def format_delta(a: float | None, b: float | None) -> str:
    if a is None or b is None:
        return "insufficient_data"
    return f"{(b - a):.6g}"


def format_pos(value: tuple[float, float, float] | None) -> str:
    if value is None:
        return "insufficient_data"
    return f"({value[1]:.6g}, {value[2]:.6g}) um"


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
