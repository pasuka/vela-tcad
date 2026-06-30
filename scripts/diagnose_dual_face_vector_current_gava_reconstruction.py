#!/usr/bin/env python3
"""Reconstruct cell current vectors from true dual-face normal fluxes."""

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
PREVIOUS_FULL_QG_RATIO = 0.246
PREVIOUS_RIGHT_SHOULDER_QG_RATIO = 0.0667
PREVIOUS_MAX_GAVA_X_UM = 1.08333
PREVIOUS_MAX_GAVA_Y_UM = 0.416667
WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}


@dataclass(frozen=True)
class Element:
    cell_id: int
    nodes: tuple[int, int, int]


@dataclass(frozen=True)
class DualFace:
    dual_face_id: int
    cell_id: int
    edge_id: int
    mid_x_um: float
    mid_y_um: float
    normal_x: float | None
    normal_y: float | None
    length_cm: float | None
    source_area_cm2: float | None
    orientation_sign: float
    material_region: str
    boundary_type: str
    dual_type: str

    @property
    def has_true_normal(self) -> bool:
        return (
            self.dual_type != "edge_coupled_only"
            and self.normal_x is not None
            and self.normal_y is not None
            and math.hypot(self.normal_x, self.normal_y) > 0.0
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dual-face-csv", type=Path, default=REPO / "build/diagnostics/dual_face_geometry_audit.csv")
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/a1_b1_internal_source_current_audit/avalanche_internal_source_current_audit.csv")
    parser.add_argument("--elements-csv", type=Path, default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_fields_20260629_155228/elements.csv")
    parser.add_argument("--node-compare", type=Path, default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20_fields_20260629_155228/coarse_node_field_compare_aligned.csv")
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--weight-mode", choices=("dual_face_length_cm", "source_area_cm2"), default="dual_face_length_cm")
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction_summary.md")
    return parser.parse_args()


def clean(row: dict[str, str]) -> dict[str, str]:
    return {str(k).strip(): str(v).strip() for k, v in row.items() if k is not None}


def f(value: Any, default: float = 0.0) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def optional_f(value: Any) -> float | None:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def i(value: Any, default: int = -1) -> int:
    try:
        return int(float(str(value).strip()))
    except (TypeError, ValueError):
        return default


def parse_biases(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def nearest_bias(value: float, targets: list[float], tol: float) -> float | None:
    best = min(targets, key=lambda item: abs(item - value))
    return best if abs(best - value) <= tol else None


def safe_ratio(value: float | None, reference: float | None) -> float | None:
    if value is None or reference is None or abs(reference) <= 1.0e-300:
        return None
    ratio = value / reference
    return ratio if math.isfinite(ratio) else None


def sign(value: float) -> float:
    if value < 0.0:
        return -1.0
    return 1.0


def load_dual_faces(path: Path) -> list[DualFace]:
    faces: list[DualFace] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean(raw)
            owner_area = optional_f(row.get("owner_subcv_area_cm2"))
            neighbor_area = optional_f(row.get("neighbor_subcv_area_cm2"))
            source_area = None
            if owner_area is not None and neighbor_area is not None:
                source_area = 0.5 * (owner_area + neighbor_area)
            elif owner_area is not None:
                source_area = owner_area
            faces.append(DualFace(
                dual_face_id=i(row.get("dual_face_id")),
                cell_id=i(row.get("cell_id")),
                edge_id=i(row.get("associated_primal_edge_id")),
                mid_x_um=f(row.get("dual_face_mid_x_um")),
                mid_y_um=f(row.get("dual_face_mid_y_um")),
                normal_x=optional_f(row.get("dual_face_normal_x")),
                normal_y=optional_f(row.get("dual_face_normal_y")),
                length_cm=optional_f(row.get("dual_face_length_cm")),
                source_area_cm2=source_area,
                orientation_sign=f(row.get("orientation_sign"), 1.0),
                material_region=row.get("material_region", ""),
                boundary_type=row.get("boundary_type", ""),
                dual_type=row.get("dual_type", "edge_coupled_only"),
            ))
    return faces


def load_elements(path: Path) -> dict[int, Element]:
    elements: dict[int, Element] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean(raw)
            cell_id = i(row.get("id", row.get("cell_id")))
            elements[cell_id] = Element(cell_id, (i(row.get("node0")), i(row.get("node1")), i(row.get("node2"))))
    return elements


def load_internal(path: Path, biases: list[float], tol: float) -> dict[tuple[float, int], dict[str, float]]:
    rows: dict[tuple[float, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean(raw)
            if row.get("source_location_type", "edge") != "edge":
                continue
            bias = nearest_bias(f(row.get("bias_V")), biases, tol)
            if bias is None:
                continue
            edge_id = i(row.get("source_entity_id"))
            rows[(bias, edge_id)] = {
                "Fn": f(row.get("Fn_used_V_per_cm")),
                "Fp": f(row.get("Fp_used_V_per_cm")),
                "alpha_n": f(row.get("alpha_n_used_cm_inv")),
                "alpha_p": f(row.get("alpha_p_used_cm_inv")),
                "Jn": f(row.get("Jn_mag_used_A_per_cm2")),
                "Jp": f(row.get("Jp_mag_used_A_per_cm2")),
                "area": f(row.get("source_weight_or_volume_cm2_for_2D") or row.get("contribution_volume_cm3_or_area_cm2_for_2D")),
                "e_sign": sign(f(row.get("electron_raw_signed_flux_proxy", "1"))),
                "h_sign": sign(f(row.get("hole_raw_signed_flux_proxy", "1"))),
            }
    return rows


def load_node_compare(path: Path) -> dict[tuple[float, int], dict[str, dict[str, float]]]:
    nodes: dict[tuple[float, int], dict[str, dict[str, float]]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean(raw)
            bias = f(row.get("bias_V"))
            node_id = i(row.get("node_id"))
            quantity = row.get("quantity", "")
            nodes.setdefault((bias, node_id), {})[quantity] = {
                "sent": f(row.get("sentaurus_value"), math.nan),
                "self": f(row.get("vela_value_scaled_to_sentaurus_units"), math.nan),
                "x_um": f(row.get("x_um"), math.nan),
                "y_um": f(row.get("y_um"), math.nan),
            }
    return nodes


def node_value(nodes: dict[tuple[float, int], dict[str, dict[str, float]]], bias: float, node_id: int, quantity: str, side: str = "sent") -> float:
    value = nodes.get((bias, node_id), {}).get(quantity, {}).get(side, math.nan)
    return value if math.isfinite(value) else 0.0


def node_xy(nodes: dict[tuple[float, int], dict[str, dict[str, float]]], bias: float, node_id: int) -> tuple[float, float]:
    entry = nodes.get((bias, node_id), {})
    for values in entry.values():
        x = values.get("x_um", math.nan)
        y = values.get("y_um", math.nan)
        if math.isfinite(x) and math.isfinite(y):
            return x, y
    return 0.0, 0.0


def triangle_area_cm2(coords: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = coords
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) * 1.0e-8


def ls_weighted(constraints: list[tuple[float, float, float, float]]) -> tuple[float, float] | None:
    a11 = a12 = a22 = b1 = b2 = 0.0
    for nx, ny, value, weight in constraints:
        norm = math.hypot(nx, ny)
        if norm <= 0.0:
            continue
        ux = nx / norm
        uy = ny / norm
        w = max(weight, 1.0e-300)
        a11 += w * ux * ux
        a12 += w * ux * uy
        a22 += w * uy * uy
        b1 += w * ux * value
        b2 += w * uy * value
    det = a11 * a22 - a12 * a12
    scale = max(abs(a11 * a22), abs(a12 * a12), 1.0)
    if abs(det) <= 1.0e-24 * scale:
        return None
    return (b1 * a22 - b2 * a12) / det, (a11 * b2 - a12 * b1) / det


def cell_sentaurus(nodes: dict[tuple[float, int], dict[str, dict[str, float]]], bias: float, elem: Element) -> dict[str, float]:
    def avg(quantity: str) -> float:
        return sum(node_value(nodes, bias, node, quantity, "sent") for node in elem.nodes) / 3.0
    return {
        "alpha_n_sent_cell_cm_inv": avg("electron_alpha_avalanche"),
        "alpha_p_sent_cell_cm_inv": avg("hole_alpha_avalanche"),
        "Gava_sent_cell_cm_minus3_s_minus1": avg("avalanche"),
    }


def window_members(x_um: float) -> list[str]:
    return [name for name, (lo, hi) in WINDOWS.items() if lo <= x_um <= hi]


def make_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    biases = parse_biases(args.biases)
    faces = load_dual_faces(args.dual_face_csv)
    elements = load_elements(args.elements_csv)
    internal = load_internal(args.internal_audit, biases, args.bias_tol)
    nodes = load_node_compare(args.node_compare)
    by_cell: dict[int, list[DualFace]] = defaultdict(list)
    for face in faces:
        by_cell[face.cell_id].append(face)
    true_face_count = sum(1 for face in faces if face.has_true_normal)
    rows: list[dict[str, Any]] = []

    for bias in biases:
        for cell_id, elem in elements.items():
            cell_faces = [face for face in by_cell.get(cell_id, []) if face.has_true_normal]
            coords = [node_xy(nodes, bias, node) for node in elem.nodes]
            area_cm2 = triangle_area_cm2(coords)
            sent = cell_sentaurus(nodes, bias, elem)
            cx = sum(x for x, _ in coords) / 3.0
            cy = sum(y for _, y in coords) / 3.0
            if not cell_faces:
                rows.append(base_row(bias, cell_id, cx, cy, "no_true_dual_faces", area_cm2, sent))
                continue
            constraints_n: list[tuple[float, float, float, float]] = []
            constraints_p: list[tuple[float, float, float, float]] = []
            alpha_n_sum = alpha_p_sum = fn_sum = fp_sum = weight_sum = 0.0
            for face in cell_faces:
                edge = internal.get((bias, face.edge_id))
                if edge is None:
                    continue
                weight = face.length_cm if args.weight_mode == "dual_face_length_cm" else face.source_area_cm2
                if weight is None or weight <= 0.0:
                    weight = 1.0
                orientation = face.orientation_sign
                nx = face.normal_x or 0.0
                ny = face.normal_y or 0.0
                constraints_n.append((nx, ny, orientation * edge["e_sign"] * edge["Jn"], weight))
                constraints_p.append((nx, ny, orientation * edge["h_sign"] * edge["Jp"], weight))
                alpha_n_sum += edge["alpha_n"] * weight
                alpha_p_sum += edge["alpha_p"] * weight
                fn_sum += edge["Fn"] * weight
                fp_sum += edge["Fp"] * weight
                weight_sum += weight
            jn_vec = ls_weighted(constraints_n)
            jp_vec = ls_weighted(constraints_p)
            if jn_vec is None or jp_vec is None or weight_sum <= 0.0 or area_cm2 <= 0.0:
                rows.append(base_row(bias, cell_id, cx, cy, "insufficient_dual_constraints", area_cm2, sent))
                continue
            alpha_n = alpha_n_sum / weight_sum
            alpha_p = alpha_p_sum / weight_sum
            jn_mag = math.hypot(*jn_vec)
            jp_mag = math.hypot(*jp_vec)
            gava = (alpha_n * jn_mag + alpha_p * jp_mag) / Q_C
            qg = Q_C * gava * area_cm2 / 1.0e4
            sent_qg = Q_C * sent["Gava_sent_cell_cm_minus3_s_minus1"] * area_cm2 / 1.0e4
            rows.append({
                "bias_V": bias,
                "cell_id": cell_id,
                "centroid_x_um": cx,
                "centroid_y_um": cy,
                "reconstruction_status": "ok",
                "mode": "current_cell_vector_ls_from_true_dual_face_normals",
                "dual_faces_used": len(cell_faces),
                "weight_mode": args.weight_mode,
                "Jn_x_A_per_cm2": jn_vec[0],
                "Jn_y_A_per_cm2": jn_vec[1],
                "Jn_mag_A_per_cm2": jn_mag,
                "Jp_x_A_per_cm2": jp_vec[0],
                "Jp_y_A_per_cm2": jp_vec[1],
                "Jp_mag_A_per_cm2": jp_mag,
                "alpha_n_cell_cm_inv": alpha_n,
                "alpha_p_cell_cm_inv": alpha_p,
                "Fn_cell_V_per_cm": fn_sum / weight_sum,
                "Fp_cell_V_per_cm": fp_sum / weight_sum,
                "Gava_cell_cm_minus3_s_minus1": gava,
                "qG_cell_contribution_A_per_um": qg,
                **sent,
                "qG_sent_cell_contribution_A_per_um": sent_qg,
                "qG_ratio_to_sentaurus": safe_ratio(qg, sent_qg),
                "window_label": window_label(cx),
            })
    return rows, {"true_face_count": true_face_count, "total_face_count": len(faces)}


def base_row(bias: float, cell_id: int, cx: float, cy: float, status: str, area_cm2: float, sent: dict[str, float]) -> dict[str, Any]:
    sent_qg = Q_C * sent["Gava_sent_cell_cm_minus3_s_minus1"] * area_cm2 / 1.0e4
    return {
        "bias_V": bias,
        "cell_id": cell_id,
        "centroid_x_um": cx,
        "centroid_y_um": cy,
        "reconstruction_status": status,
        "mode": "current_cell_vector_ls_from_true_dual_face_normals",
        "dual_faces_used": 0,
        "weight_mode": "",
        "Jn_x_A_per_cm2": None,
        "Jn_y_A_per_cm2": None,
        "Jn_mag_A_per_cm2": None,
        "Jp_x_A_per_cm2": None,
        "Jp_y_A_per_cm2": None,
        "Jp_mag_A_per_cm2": None,
        "alpha_n_cell_cm_inv": None,
        "alpha_p_cell_cm_inv": None,
        "Fn_cell_V_per_cm": None,
        "Fp_cell_V_per_cm": None,
        "Gava_cell_cm_minus3_s_minus1": None,
        "qG_cell_contribution_A_per_um": None,
        **sent,
        "qG_sent_cell_contribution_A_per_um": sent_qg,
        "qG_ratio_to_sentaurus": None,
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


HEADER = [
    "bias_V", "cell_id", "centroid_x_um", "centroid_y_um", "reconstruction_status", "mode",
    "dual_faces_used", "weight_mode",
    "Jn_x_A_per_cm2", "Jn_y_A_per_cm2", "Jn_mag_A_per_cm2",
    "Jp_x_A_per_cm2", "Jp_y_A_per_cm2", "Jp_mag_A_per_cm2",
    "alpha_n_cell_cm_inv", "alpha_p_cell_cm_inv", "Fn_cell_V_per_cm", "Fp_cell_V_per_cm",
    "Gava_cell_cm_minus3_s_minus1", "qG_cell_contribution_A_per_um",
    "alpha_n_sent_cell_cm_inv", "alpha_p_sent_cell_cm_inv",
    "Gava_sent_cell_cm_minus3_s_minus1", "qG_sent_cell_contribution_A_per_um",
    "qG_ratio_to_sentaurus", "window_label",
]


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.17g}" if math.isfinite(value) else ""
    return str(value)


def write_csv_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in HEADER})


def aggregate_ratios(rows: list[dict[str, Any]]) -> dict[tuple[float, str], float | None]:
    buckets: dict[tuple[float, str], dict[str, float]] = defaultdict(lambda: {"self": 0.0, "sent": 0.0})
    for row in rows:
        if row.get("reconstruction_status") != "ok":
            continue
        bias = float(row["bias_V"])
        x = float(row["centroid_x_um"])
        for window in window_members(x):
            bucket = buckets[(bias, window)]
            bucket["self"] += float(row["qG_cell_contribution_A_per_um"])
            bucket["sent"] += float(row["qG_sent_cell_contribution_A_per_um"])
    return {key: safe_ratio(value["self"], value["sent"]) for key, value in buckets.items()}


def max_gava(rows: list[dict[str, Any]], bias: float) -> tuple[float, float] | None:
    candidates = [
        row for row in rows
        if row.get("reconstruction_status") == "ok" and abs(float(row["bias_V"]) - bias) <= 1.0e-9
    ]
    if not candidates:
        return None
    best = max(candidates, key=lambda row: abs(float(row["Gava_cell_cm_minus3_s_minus1"])))
    return float(best["centroid_x_um"]), float(best["centroid_y_um"])


def write_summary(path: Path, rows: list[dict[str, Any]], meta: dict[str, Any]) -> None:
    ratios = aggregate_ratios(rows)
    full = ratios.get((-20.0, "full"))
    right = ratios.get((-20.0, "right_shoulder"))
    ok_count = sum(1 for row in rows if row.get("reconstruction_status") == "ok")
    no_true = meta["true_face_count"] == 0
    full_improved = full is not None and full > PREVIOUS_FULL_QG_RATIO
    right_improved = right is not None and right > PREVIOUS_RIGHT_SHOULDER_QG_RATIO
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Dual Face Vector Current Gava Reconstruction\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver path changed: no\n")
        out.write("- current-density unit: A/cm^2\n")
        out.write("- alpha unit: cm^-1\n")
        out.write("- Gava unit: cm^-3 s^-1\n")
        out.write("- qG unit: A/um\n")
        out.write(f"- total dual-face rows: {meta['total_face_count']}\n")
        out.write(f"- true dual-face rows: {meta['true_face_count']}\n")
        out.write(f"- reconstructed cell rows: {ok_count}\n\n")
        out.write("## -20V Window Ratios\n\n")
        out.write("| window | qG ratio to Sentaurus |\n")
        out.write("| --- | --- |\n")
        for window in ("full", "junction", "center", "left_shoulder", "right_shoulder"):
            out.write(f"| {window} | {format_ratio(ratios.get((-20.0, window)))} |\n")
        out.write("\n## Required Answers\n\n")
        if no_true:
            out.write("1. Does true dual-face vector reconstruction improve -20V full qG beyond previous 0.246? no; no true dual-face normals are available.\n")
            out.write("2. Does it improve right_shoulder beyond previous 0.0667? no; no true dual-face normals are available.\n")
            out.write("3. Does max Gava location improve? no; reconstruction was not possible.\n")
            out.write("4. recommend guarded AvalancheCurrentMagnitudeMode=dual_face_vector_mag: no.\n")
            out.write("5. Return to SG formula / current magnitude model, or first extend mesh geometry to export true dual-face normals.\n")
        else:
            out.write(f"1. Does true dual-face vector reconstruction improve -20V full qG beyond previous 0.246? {'yes' if full_improved else 'no'}; ratio={format_ratio(full)}.\n")
            out.write(f"2. Does it improve right_shoulder beyond previous 0.0667? {'yes' if right_improved else 'no'}; ratio={format_ratio(right)}.\n")
            max_pos = max_gava(rows, -20.0)
            max_changed = max_pos is not None and math.hypot(max_pos[0] - PREVIOUS_MAX_GAVA_X_UM, max_pos[1] - PREVIOUS_MAX_GAVA_Y_UM) > 1.0e-4
            out.write(f"3. Does max Gava location improve? {'yes' if max_changed else 'no'}; reconstructed max at {format_pos(max_pos)}.\n")
            recommend = full_improved and right_improved
            out.write(f"4. recommend guarded AvalancheCurrentMagnitudeMode=dual_face_vector_mag: {'yes' if recommend else 'no'} for qG magnitude; hotspot location still needs separate validation.\n")
            if recommend:
                out.write("5. Keep this as a guarded optional avalanche current magnitude mode before coupled solve experiments, but continue source mapping/hotspot-position checks.\n")
            else:
                out.write("5. Return to SG formula / current magnitude model.\n")


def format_ratio(value: float | None) -> str:
    return "insufficient_data" if value is None else f"{value:.6g}"


def format_pos(value: tuple[float, float] | None) -> str:
    return "insufficient_data" if value is None else f"({value[0]:.6g}, {value[1]:.6g}) um"


def run(args: argparse.Namespace) -> None:
    rows, meta = make_rows(args)
    write_csv_rows(args.out_csv, rows)
    write_summary(args.out_summary, rows, meta)


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


