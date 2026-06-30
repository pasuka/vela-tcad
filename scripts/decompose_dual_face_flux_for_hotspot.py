#!/usr/bin/env python3
"""Decompose true-dual face fluxes for a local self Gava hotspot cell."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any

import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19


def clean(row: dict[str, str]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items() if key is not None}


def f(value: Any, default: float = 0.0) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def opt_f(value: Any) -> float | None:
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


def sign(value: float) -> float:
    return -1.0 if value < 0.0 else 1.0


def safe_ratio(value: float | None, ref: float | None) -> float | None:
    if value is None or ref is None or abs(ref) <= 1.0e-300:
        return None
    return value / ref


def read_rows(path: Path) -> list[dict[str, str]]:
    return [clean(row) for row in bv.read_rows(path)]


def nearest_bias(value: float, target: float, tol: float) -> bool:
    return abs(value - target) <= tol


def load_dual_faces(path: Path, cell_id: int) -> list[dict[str, Any]]:
    faces: list[dict[str, Any]] = []
    for row in read_rows(path):
        if i(row.get("cell_id")) != cell_id:
            continue
        faces.append({
            "face_id": i(row.get("dual_face_id")),
            "cell_id": cell_id,
            "owner_node_id": i(row.get("owner_node_id")),
            "neighbor_node_id": i(row.get("neighbor_node_id")),
            "edge_id": i(row.get("associated_primal_edge_id")),
            "face_mid_x_um": f(row.get("dual_face_mid_x_um")),
            "face_mid_y_um": f(row.get("dual_face_mid_y_um")),
            "normal_x": f(row.get("dual_face_normal_x")),
            "normal_y": f(row.get("dual_face_normal_y")),
            "face_length_cm": f(row.get("dual_face_length_cm")),
            "orientation_sign": f(row.get("orientation_sign"), 1.0),
            "boundary_type": row.get("boundary_type", ""),
            "dual_type": row.get("dual_type", ""),
        })
    return sorted(faces, key=lambda row: row["face_id"])


def load_internal(path: Path, bias: float, tol: float) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    for row in read_rows(path):
        if row.get("source_location_type", "edge") != "edge":
            continue
        row_bias = f(row.get("bias_V"))
        if not nearest_bias(row_bias, bias, tol):
            continue
        edge_id = i(row.get("source_entity_id"))
        rows[edge_id] = {
            "Fn": f(row.get("Fn_used_V_per_cm")),
            "Fp": f(row.get("Fp_used_V_per_cm")),
            "alpha_n": f(row.get("alpha_n_used_cm_inv")),
            "alpha_p": f(row.get("alpha_p_used_cm_inv")),
            "Jn": f(row.get("Jn_mag_used_A_per_cm2")),
            "Jp": f(row.get("Jp_mag_used_A_per_cm2")),
            "e_sign": sign(f(row.get("electron_raw_signed_flux_proxy", "1"))),
            "h_sign": sign(f(row.get("hole_raw_signed_flux_proxy", "1"))),
        }
    return rows


def load_node_compare(path: Path, bias: float, tol: float) -> dict[int, dict[str, float]]:
    nodes: dict[int, dict[str, float]] = {}
    for row in read_rows(path):
        if not nearest_bias(f(row.get("bias_V")), bias, tol):
            continue
        node_id = i(row.get("node_id"))
        quantity = row.get("quantity", "")
        value = f(row.get("vela_value_scaled_to_sentaurus_units"))
        item = nodes.setdefault(node_id, {})
        item[quantity] = value
        item["x_um"] = f(row.get("x_um"), item.get("x_um", 0.0))
        item["y_um"] = f(row.get("y_um"), item.get("y_um", 0.0))
    return nodes


def choose_cell(path: Path, bias: float, target_x: float, target_y: float, cell_id: int | None) -> dict[str, Any]:
    best: dict[str, Any] | None = None
    best_dist = math.inf
    for row in read_rows(path):
        if not nearest_bias(f(row.get("bias_V")), bias, 1.0e-6):
            continue
        if row.get("reconstruction_status", "ok") != "ok":
            continue
        current_id = i(row.get("cell_id"))
        x = f(row.get("centroid_x_um"))
        y = f(row.get("centroid_y_um"))
        dist = math.hypot(x - target_x, y - target_y)
        if cell_id is not None and current_id != cell_id:
            continue
        if dist < best_dist:
            best = row
            best_dist = dist
    if best is None:
        raise SystemExit("target hotspot cell not found")
    return best


def invert_2x2(a00: float, a01: float, a11: float) -> tuple[float, float, float] | None:
    det = a00 * a11 - a01 * a01
    scale = max(abs(a00 * a11), abs(a01 * a01), 1.0)
    if abs(det) <= 1.0e-24 * scale:
        return None
    return a11 / det, -a01 / det, a00 / det


def ls_matrix(constraints: list[dict[str, Any]]) -> tuple[tuple[float, float, float], tuple[float, float], tuple[float, float] | None]:
    a00 = a01 = a11 = b0 = b1 = 0.0
    for row in constraints:
        w = row["weight"]
        nx = row["normal_x"]
        ny = row["normal_y"]
        scalar = row["scalar"]
        a00 += w * nx * nx
        a01 += w * nx * ny
        a11 += w * ny * ny
        b0 += w * nx * scalar
        b1 += w * ny * scalar
    inv = invert_2x2(a00, a01, a11)
    if inv is None:
        return (a00, a01, a11), (b0, b1), None
    i00, i01, i11 = inv
    return (a00, a01, a11), (b0, b1), (i00 * b0 + i01 * b1, i01 * b0 + i11 * b1)


def contribution(inv: tuple[float, float, float] | None, normal_x: float, normal_y: float, scalar: float, weight: float) -> tuple[float | None, float | None]:
    if inv is None:
        return None, None
    i00, i01, i11 = inv
    bx = weight * normal_x * scalar
    by = weight * normal_y * scalar
    return i00 * bx + i01 * by, i01 * bx + i11 * by


def node_value(nodes: dict[int, dict[str, float]], node_id: int, quantity: str) -> float:
    return nodes.get(node_id, {}).get(quantity, 0.0)


def face_node_quantities(face: dict[str, Any], nodes: dict[int, dict[str, float]]) -> dict[str, Any]:
    owner = face["owner_node_id"]
    neighbor = face["neighbor_node_id"]
    phip0 = node_value(nodes, owner, "hole_qf")
    phip1 = node_value(nodes, neighbor, "hole_qf")
    p0 = node_value(nodes, owner, "hole_density")
    p1 = node_value(nodes, neighbor, "hole_density")
    mu0 = node_value(nodes, owner, "hole_mobility")
    mu1 = node_value(nodes, neighbor, "hole_mobility")
    x0 = nodes.get(owner, {}).get("x_um", 0.0)
    y0 = nodes.get(owner, {}).get("y_um", 0.0)
    x1 = nodes.get(neighbor, {}).get("x_um", 0.0)
    y1 = nodes.get(neighbor, {}).get("y_um", 0.0)
    dist_cm = math.hypot(x1 - x0, y1 - y0) * 1.0e-4
    return {
        "associated_nodes": f"{owner};{neighbor}",
        "phip_node0_V": phip0,
        "phip_node1_V": phip1,
        "Fp_across_face_V_per_cm": abs(phip1 - phip0) / dist_cm if dist_cm > 0.0 else 0.0,
        "p_density_edge_cm_minus3": 0.5 * (p0 + p1),
        "mu_p_edge_cm2_per_Vs": 0.5 * (mu0 + mu1),
    }


def build_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    target = choose_cell(
        args.self_dual_cell_csv,
        args.bias,
        args.target_x_um,
        args.target_y_um,
        args.cell_id,
    )
    cell_id = i(target.get("cell_id"))
    faces = load_dual_faces(args.dual_face_csv, cell_id)
    internal = load_internal(args.internal_audit, args.bias, args.bias_tol)
    nodes = load_node_compare(args.node_compare, args.bias, args.bias_tol)

    enriched: list[dict[str, Any]] = []
    constraints_n: list[dict[str, Any]] = []
    constraints_p: list[dict[str, Any]] = []
    for face in faces:
        edge = internal.get(face["edge_id"])
        status = "ok" if edge is not None else "no_associated_primal_edge"
        weight = face["face_length_cm"] if face["face_length_cm"] > 0.0 else 1.0
        jn_scalar = face["orientation_sign"] * (edge["e_sign"] if edge else 1.0) * (edge["Jn"] if edge else 0.0)
        jp_scalar = face["orientation_sign"] * (edge["h_sign"] if edge else 1.0) * (edge["Jp"] if edge else 0.0)
        item = dict(face)
        item.update(face_node_quantities(face, nodes))
        item.update({
            "edge_lookup_status": status,
            "weight": weight,
            "Jn_face_scalar_A_per_cm2": jn_scalar,
            "Jp_face_scalar_A_per_cm2": jp_scalar,
            "Fn_face_V_per_cm": edge["Fn"] if edge else None,
            "Fp_face_V_per_cm": edge["Fp"] if edge else None,
            "alpha_p_cm_inv": edge["alpha_p"] if edge else None,
            "Gava_p_face_cm_minus3_s_minus1": (edge["alpha_p"] * abs(jp_scalar) / Q_C) if edge else None,
        })
        enriched.append(item)
        if edge is not None:
            constraints_n.append({"normal_x": face["normal_x"], "normal_y": face["normal_y"], "scalar": jn_scalar, "weight": weight, "face_id": face["face_id"]})
            constraints_p.append({"normal_x": face["normal_x"], "normal_y": face["normal_y"], "scalar": jp_scalar, "weight": weight, "face_id": face["face_id"]})

    a_n, _, jn_vec = ls_matrix(constraints_n)
    a_p, _, jp_vec = ls_matrix(constraints_p)
    inv_n = invert_2x2(*a_n)
    inv_p = invert_2x2(*a_p)

    rows: list[dict[str, Any]] = []
    for item in enriched:
        cnx, cny = contribution(inv_n, item["normal_x"], item["normal_y"], item["Jn_face_scalar_A_per_cm2"], item["weight"])
        cpx, cpy = contribution(inv_p, item["normal_x"], item["normal_y"], item["Jp_face_scalar_A_per_cm2"], item["weight"])
        residual_n = None if jn_vec is None else item["normal_x"] * jn_vec[0] + item["normal_y"] * jn_vec[1] - item["Jn_face_scalar_A_per_cm2"]
        residual_p = None if jp_vec is None else item["normal_x"] * jp_vec[0] + item["normal_y"] * jp_vec[1] - item["Jp_face_scalar_A_per_cm2"]
        rows.append({
            "bias_V": args.bias,
            "cell_id": cell_id,
            "cell_centroid_x_um": f(target.get("centroid_x_um")),
            "cell_centroid_y_um": f(target.get("centroid_y_um")),
            "face_id": item["face_id"],
            "associated_edge_id": item["edge_id"] if item["edge_id"] >= 0 else "",
            "edge_lookup_status": item["edge_lookup_status"],
            "face_mid_x_um": item["face_mid_x_um"],
            "face_mid_y_um": item["face_mid_y_um"],
            "normal_x": item["normal_x"],
            "normal_y": item["normal_y"],
            "face_length_cm": item["face_length_cm"],
            "boundary_type": item["boundary_type"],
            "dual_type": item["dual_type"],
            "Jp_face_scalar_A_per_cm2": item["Jp_face_scalar_A_per_cm2"],
            "Jn_face_scalar_A_per_cm2": item["Jn_face_scalar_A_per_cm2"],
            "Jp_vector_reconstruction_contribution_x_A_per_cm2": cpx,
            "Jp_vector_reconstruction_contribution_y_A_per_cm2": cpy,
            "Jn_vector_reconstruction_contribution_x_A_per_cm2": cnx,
            "Jn_vector_reconstruction_contribution_y_A_per_cm2": cny,
            "weight": item["weight"],
            "residual_p_A_per_cm2": residual_p,
            "residual_n_A_per_cm2": residual_n,
            "associated_nodes": item["associated_nodes"],
            "phip_node0_V": item["phip_node0_V"],
            "phip_node1_V": item["phip_node1_V"],
            "Fp_across_face_V_per_cm": item["Fp_across_face_V_per_cm"],
            "p_density_edge_cm_minus3": item["p_density_edge_cm_minus3"],
            "mu_p_edge_cm2_per_Vs": item["mu_p_edge_cm2_per_Vs"],
            "alpha_p_cm_inv": item["alpha_p_cm_inv"],
            "Gava_p_face_cm_minus3_s_minus1": item["Gava_p_face_cm_minus3_s_minus1"],
            "mesh_diagonal_alignment_abs": max(abs(item["normal_x"] + item["normal_y"]), abs(item["normal_x"] - item["normal_y"])),
        })
    meta = {
        "cell_id": cell_id,
        "centroid_x": f(target.get("centroid_x_um")),
        "centroid_y": f(target.get("centroid_y_um")),
        "Jp_x": f(target.get("Jp_x_A_per_cm2")),
        "Jp_y": f(target.get("Jp_y_A_per_cm2")),
        "Jp_mag": f(target.get("Jp_mag_A_per_cm2")),
        "Gava": f(target.get("Gava_cell_cm_minus3_s_minus1")),
    }
    return rows, meta


FIELDS = [
    "bias_V", "cell_id", "cell_centroid_x_um", "cell_centroid_y_um",
    "face_id", "associated_edge_id", "edge_lookup_status",
    "face_mid_x_um", "face_mid_y_um", "normal_x", "normal_y", "face_length_cm",
    "boundary_type", "dual_type",
    "Jp_face_scalar_A_per_cm2", "Jn_face_scalar_A_per_cm2",
    "Jp_vector_reconstruction_contribution_x_A_per_cm2",
    "Jp_vector_reconstruction_contribution_y_A_per_cm2",
    "Jn_vector_reconstruction_contribution_x_A_per_cm2",
    "Jn_vector_reconstruction_contribution_y_A_per_cm2",
    "weight", "residual_p_A_per_cm2", "residual_n_A_per_cm2",
    "associated_nodes", "phip_node0_V", "phip_node1_V", "Fp_across_face_V_per_cm",
    "p_density_edge_cm_minus3", "mu_p_edge_cm2_per_Vs", "alpha_p_cm_inv",
    "Gava_p_face_cm_minus3_s_minus1", "mesh_diagonal_alignment_abs",
]


def dominant(rows: list[dict[str, Any]], key: str) -> dict[str, Any] | None:
    usable = [row for row in rows if row.get("edge_lookup_status") == "ok"]
    if not usable:
        return None
    return max(usable, key=lambda row: abs(f(row.get(key))))


def tied_dominant_faces(rows: list[dict[str, Any]], key: str, rel_tol: float = 1.0e-9) -> list[dict[str, Any]]:
    usable = [row for row in rows if row.get("edge_lookup_status") == "ok"]
    if not usable:
        return []
    top = max(abs(f(row.get(key))) for row in usable)
    if top <= 0.0:
        return []
    return [row for row in usable if abs(abs(f(row.get(key))) - top) <= rel_tol * top]


def fmt(value: Any) -> str:
    return bv.fmt(value)


def write_summary(path: Path, rows: list[dict[str, Any]], meta: dict[str, Any]) -> None:
    dom_jp = dominant(rows, "Jp_face_scalar_A_per_cm2")
    tied_jp = tied_dominant_faces(rows, "Jp_face_scalar_A_per_cm2")
    tied_face_ids = ";".join(str(row.get("face_id")) for row in tied_jp)
    tied_edge_ids = ";".join(sorted({str(row.get("associated_edge_id")) for row in tied_jp}))
    boundary = dom_jp is not None and "boundary" in str(dom_jp.get("boundary_type", "")).lower()
    top_edge = dom_jp is not None and f(dom_jp.get("face_mid_y_um")) >= f(meta.get("centroid_y"))
    diagonal = dom_jp is not None and f(dom_jp.get("mesh_diagonal_alignment_abs")) > 1.3
    missing = sum(1 for row in rows if row.get("edge_lookup_status") != "ok")
    lines = [
        "# Dual-Face Flux Decomposition For Self Max Gava Cell",
        "",
        f"- bias_V: {fmt(rows[0]['bias_V'] if rows else '')}",
        f"- cell_id: {meta.get('cell_id')}",
        f"- centroid_um: ({fmt(meta.get('centroid_x'))}, {fmt(meta.get('centroid_y'))})",
        f"- reconstructed Jp vector: ({fmt(meta.get('Jp_x'))}, {fmt(meta.get('Jp_y'))}) A/cm^2, |Jp|={fmt(meta.get('Jp_mag'))}",
        f"- Gava_cell: {fmt(meta.get('Gava'))} cm^-3 s^-1",
        f"- dual_faces: {len(rows)}, missing_flux_constraints: {missing}",
        "",
        "## Answers",
        "",
        "1. Which face dominates reconstructed Jp magnitude? "
        f"face_id={'' if dom_jp is None else dom_jp.get('face_id')} (tied_faces={tied_face_ids}), "
        f"edge={'' if dom_jp is None else dom_jp.get('associated_edge_id')} (tied_edges={tied_edge_ids}), "
        f"Jp_face_scalar={fmt(None if dom_jp is None else dom_jp.get('Jp_face_scalar_A_per_cm2'))}, "
        f"contribution=({fmt(None if dom_jp is None else dom_jp.get('Jp_vector_reconstruction_contribution_x_A_per_cm2'))}, "
        f"{fmt(None if dom_jp is None else dom_jp.get('Jp_vector_reconstruction_contribution_y_A_per_cm2'))}) A/cm^2.",
        f"2. Is Jp dominated by a boundary/top-edge face? boundary={boundary}, top_or_above_centroid={top_edge}.",
        f"3. Does the face align with mesh diagonal orientation? {'yes' if diagonal else 'no'}; normal=({fmt(None if dom_jp is None else dom_jp.get('normal_x'))}, {fmt(None if dom_jp is None else dom_jp.get('normal_y'))}).",
        "4. Would changing source mapping from cell centroid to node average reduce the hotspot? "
        "Likely yes for visualization/location if the dominant cell source is spread to its nodes; it will not remove the local Jp flux excess itself.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'dual_face_flux_decomposition_self_max_cell.csv'}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bias", type=float, default=-20.0)
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--target-x-um", type=float, default=1.0833333333333333)
    parser.add_argument("--target-y-um", type=float, default=0.4166666666666667)
    parser.add_argument("--cell-id", type=int, default=None)
    parser.add_argument("--dual-face-csv", type=Path, default=REPO / "build/diagnostics/dual_face_geometry_audit.csv")
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/a1_b1_internal_source_current_audit/avalanche_internal_source_current_audit.csv")
    parser.add_argument("--node-compare", type=Path, default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/coarse_node_field_compare_aligned.csv")
    parser.add_argument("--self-dual-cell-csv", type=Path, default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv")
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/dual_face_flux_decomposition_self_max_cell.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/dual_face_flux_decomposition_self_max_cell_summary.md")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows, meta = build_rows(args)
    bv.write_rows(args.out_csv, rows, FIELDS)
    write_summary(args.out_summary, rows, meta)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
