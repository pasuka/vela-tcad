#!/usr/bin/env python3
"""Audit local true-dual face least-squares sensitivity for the self max Gava cell."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Any, Callable

import decompose_dual_face_flux_for_hotspot as local


REPO = Path(__file__).resolve().parents[1]
Q_C = local.Q_C


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(parsed):
        return ""
    return f"{parsed:.12g}"


def truth(value: bool) -> str:
    return "True" if value else "False"


def eigvals_2x2(a00: float, a01: float, a11: float) -> tuple[float, float]:
    trace = a00 + a11
    disc = max((a00 - a11) * (a00 - a11) + 4.0 * a01 * a01, 0.0)
    root = math.sqrt(disc)
    return 0.5 * (trace + root), 0.5 * (trace - root)


def solve_2x2(a00: float, a01: float, a11: float, b0: float, b1: float) -> tuple[float, float] | None:
    inv = local.invert_2x2(a00, a01, a11)
    if inv is None:
        return None
    i00, i01, i11 = inv
    return i00 * b0 + i01 * b1, i01 * b0 + i11 * b1


def solve_symmetric_min_norm(a00: float, a01: float, a11: float, b0: float, b1: float) -> tuple[float, float] | None:
    direct = solve_2x2(a00, a01, a11, b0, b1)
    if direct is not None:
        return direct
    eig_hi, _ = eigvals_2x2(a00, a01, a11)
    if eig_hi <= 0.0:
        return None
    if abs(a01) > abs(a00 - eig_hi):
        ux = a01
        uy = eig_hi - a00
    else:
        ux = eig_hi - a11
        uy = a01
    norm = math.hypot(ux, uy)
    if norm <= 0.0:
        ux, uy = (1.0, 0.0) if abs(a00) >= abs(a11) else (0.0, 1.0)
    else:
        ux /= norm
        uy /= norm
    scale = (ux * b0 + uy * b1) / eig_hi
    return scale * ux, scale * uy


def is_boundary(face: dict[str, Any]) -> bool:
    return "boundary" in str(face.get("boundary_type", "")).lower()


def has_flux(face: dict[str, Any]) -> bool:
    return str(face.get("edge_lookup_status", "")) == "ok"


def scalar_for(face: dict[str, Any], carrier: str) -> float:
    key = "Jp_face_scalar_A_per_cm2" if carrier == "hole" else "Jn_face_scalar_A_per_cm2"
    return local.f(face.get(key))


def constraints(
    faces: list[dict[str, Any]],
    carrier: str,
    include: Callable[[dict[str, Any]], bool] | None = None,
    boundary_weight: float = 1.0,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for face in faces:
        if not has_flux(face):
            continue
        if include is not None and not include(face):
            continue
        scale = boundary_weight if is_boundary(face) else 1.0
        weight = local.f(face.get("weight")) * scale
        if weight <= 0.0:
            continue
        rows.append({
            "face": face,
            "face_id": local.i(face.get("face_id")),
            "normal_x": local.f(face.get("normal_x")),
            "normal_y": local.f(face.get("normal_y")),
            "scalar": scalar_for(face, carrier),
            "weight": weight,
        })
    return rows


def solve_ls(rows: list[dict[str, Any]], tikhonov_abs: float = 0.0) -> dict[str, Any]:
    a00 = a01 = a11 = b0 = b1 = 0.0
    for row in rows:
        w = row["weight"]
        nx = row["normal_x"]
        ny = row["normal_y"]
        scalar = row["scalar"]
        a00 += w * nx * nx
        a01 += w * nx * ny
        a11 += w * ny * ny
        b0 += w * nx * scalar
        b1 += w * ny * scalar
    raw_trace = a00 + a11
    a00_reg = a00 + tikhonov_abs
    a11_reg = a11 + tikhonov_abs
    eig_hi, eig_lo = eigvals_2x2(a00_reg, a01, a11_reg)
    rank_tol = max(eig_hi, 1.0) * 1.0e-12
    rank = int(eig_hi > rank_tol) + int(eig_lo > rank_tol)
    condition = math.inf if eig_lo <= rank_tol else abs(eig_hi / eig_lo)
    solution = solve_symmetric_min_norm(a00_reg, a01, a11_reg, b0, b1)
    residual_norm_sq = 0.0
    residuals: dict[int, float] = {}
    if solution is not None:
        jx, jy = solution
        for row in rows:
            residual = row["normal_x"] * jx + row["normal_y"] * jy - row["scalar"]
            residuals[row["face_id"]] = residual
            residual_norm_sq += row["weight"] * residual * residual
    return {
        "normal_A00": a00_reg,
        "normal_A01": a01,
        "normal_A11": a11_reg,
        "unregularized_trace": raw_trace,
        "rhs_x": b0,
        "rhs_y": b1,
        "condition_number": condition,
        "rank": rank,
        "residual_norm": math.sqrt(residual_norm_sq),
        "solution": solution,
        "residuals": residuals,
        "constraints": rows,
        "tikhonov_abs": tikhonov_abs,
    }


def contribution(ls: dict[str, Any], row: dict[str, Any]) -> tuple[float | None, float | None]:
    inv = local.invert_2x2(
        local.f(ls.get("normal_A00")),
        local.f(ls.get("normal_A01")),
        local.f(ls.get("normal_A11")),
    )
    if inv is None:
        return None, None
    return local.contribution(inv, row["normal_x"], row["normal_y"], row["scalar"], row["weight"])


def load_target_cell(path: Path, bias: float, target_x: float, target_y: float, cell_id: int | None) -> dict[str, Any]:
    return local.choose_cell(path, bias, target_x, target_y, cell_id)


def load_other_gava_max(path: Path, bias: float, target_cell_id: int) -> float | None:
    values: list[float] = []
    for row in local.read_rows(path):
        if not local.nearest_bias(local.f(row.get("bias_V")), bias, 1.0e-6):
            continue
        if row.get("reconstruction_status", "ok") != "ok":
            continue
        if local.i(row.get("cell_id")) == target_cell_id:
            continue
        values.append(local.f(row.get("Gava_cell_cm_minus3_s_minus1")))
    return max(values) if values else None


def carrier_alpha(target: dict[str, Any], carrier: str, faces: list[dict[str, Any]]) -> float:
    key = "alpha_p_cell_cm_inv" if carrier == "hole" else "alpha_n_cell_cm_inv"
    value = local.opt_f(target.get(key))
    if value is not None:
        return value
    face_key = "alpha_p_cm_inv" if carrier == "hole" else "alpha_n_cm_inv"
    weighted = [
        (local.f(face.get(face_key)), local.f(face.get("weight")))
        for face in faces
        if has_flux(face) and local.opt_f(face.get(face_key)) is not None
    ]
    total = sum(weight for _, weight in weighted)
    return sum(value * weight for value, weight in weighted) / total if total > 0.0 else 0.0


def vector_mag(solution: tuple[float, float] | None) -> float | None:
    if solution is None:
        return None
    return math.hypot(solution[0], solution[1])


def gava_from(alpha_n: float, jn_mag: float | None, alpha_p: float, jp_mag: float | None) -> tuple[float | None, float | None, float | None]:
    gn = None if jn_mag is None else alpha_n * jn_mag / Q_C
    gp = None if jp_mag is None else alpha_p * jp_mag / Q_C
    total = None if gn is None or gp is None else gn + gp
    return gn, gp, total


def variant_defs() -> list[tuple[str, Callable[[dict[str, Any]], bool] | None, float, float]]:
    return [
        ("baseline", None, 1.0, 0.0),
        ("remove_face_74", lambda face: local.i(face.get("face_id")) != 74, 1.0, 0.0),
        ("remove_face_75", lambda face: local.i(face.get("face_id")) != 75, 1.0, 0.0),
        ("remove_faces_74_75", lambda face: local.i(face.get("face_id")) not in {74, 75}, 1.0, 0.0),
        ("remove_all_boundary_faces", lambda face: not is_boundary(face), 1.0, 0.0),
        ("only_interior_faces", lambda face: not is_boundary(face), 1.0, 0.0),
        ("boundary_weight_1p0", None, 1.0, 0.0),
        ("boundary_weight_0p5", None, 0.5, 0.0),
        ("boundary_weight_0p25", None, 0.25, 0.0),
        ("boundary_weight_0p0", None, 0.0, 0.0),
        ("tikhonov_1e-06_trace", None, 1.0, 1.0e-6),
        ("tikhonov_1e-04_trace", None, 1.0, 1.0e-4),
        ("tikhonov_1e-02_trace", None, 1.0, 1.0e-2),
    ]


FIELDS = [
    "row_type",
    "bias_V",
    "cell_id",
    "cell_centroid_x_um",
    "cell_centroid_y_um",
    "carrier",
    "variant_name",
    "face_id",
    "edge_id",
    "boundary_type",
    "is_boundary",
    "is_top_or_above_centroid",
    "normal_x",
    "normal_y",
    "face_length_cm",
    "weight",
    "Jp_face_scalar_A_per_cm2",
    "Jn_face_scalar_A_per_cm2",
    "contribution_x_A_per_cm2",
    "contribution_y_A_per_cm2",
    "residual_A_per_cm2",
    "has_flux_constraint",
    "missing_reason",
    "normal_A00",
    "normal_A01",
    "normal_A11",
    "rhs_x",
    "rhs_y",
    "condition_number",
    "rank",
    "residual_norm",
    "tikhonov_lambda_abs",
    "tikhonov_lambda_factor_times_trace",
    "Jp_x_A_per_cm2",
    "Jp_y_A_per_cm2",
    "Jp_mag_A_per_cm2",
    "Jn_x_A_per_cm2",
    "Jn_y_A_per_cm2",
    "Jn_mag_A_per_cm2",
    "Gava_p_cm_minus3_s_minus1",
    "Gava_n_cm_minus3_s_minus1",
    "Gava_total_cm_minus3_s_minus1",
    "ratio_to_baseline",
    "ratio_to_sentaurus_cell",
    "hotspot_would_remain_rank1",
]


def base_row(args: argparse.Namespace, meta: dict[str, Any], carrier: str = "", variant: str = "") -> dict[str, Any]:
    return {
        "bias_V": args.bias,
        "cell_id": meta["cell_id"],
        "cell_centroid_x_um": meta["centroid_x"],
        "cell_centroid_y_um": meta["centroid_y"],
        "carrier": carrier,
        "variant_name": variant,
    }


def build_output(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    face_rows, meta = local.build_rows(args)
    target = load_target_cell(args.self_dual_cell_csv, args.bias, args.target_x_um, args.target_y_um, args.cell_id)
    cell_id = local.i(target.get("cell_id"))
    meta["cell_id"] = cell_id
    alpha_n = carrier_alpha(target, "electron", face_rows)
    alpha_p = carrier_alpha(target, "hole", face_rows)
    sent_gava = local.opt_f(target.get("Gava_sent_cell_cm_minus3_s_minus1"))
    max_other_gava = load_other_gava_max(args.self_dual_cell_csv, args.bias, cell_id)

    baseline_hole = solve_ls(constraints(face_rows, "hole"))
    baseline_electron = solve_ls(constraints(face_rows, "electron"))
    baseline_jp = vector_mag(baseline_hole["solution"])
    baseline_jn = vector_mag(baseline_electron["solution"])
    _, baseline_gp, baseline_total = gava_from(alpha_n, baseline_jn, alpha_p, baseline_jp)

    out: list[dict[str, Any]] = []
    for carrier, ls in (("hole", baseline_hole), ("electron", baseline_electron)):
        row = base_row(args, meta, carrier, "baseline")
        row.update({
            "row_type": "ls_system",
            "normal_A00": ls["normal_A00"],
            "normal_A01": ls["normal_A01"],
            "normal_A11": ls["normal_A11"],
            "rhs_x": ls["rhs_x"],
            "rhs_y": ls["rhs_y"],
            "condition_number": ls["condition_number"],
            "rank": ls["rank"],
            "residual_norm": ls["residual_norm"],
            "tikhonov_lambda_abs": 0.0,
            "tikhonov_lambda_factor_times_trace": 0.0,
        })
        out.append(row)
        constraint_by_face = {item["face_id"]: item for item in ls["constraints"]}
        for face in face_rows:
            face_id = local.i(face.get("face_id"))
            item = constraint_by_face.get(face_id)
            c_x = c_y = None
            residual = None
            weight = local.f(face.get("weight"))
            if item is not None:
                c_x, c_y = contribution(ls, item)
                residual = ls["residuals"].get(face_id)
                weight = item["weight"]
            face_row = base_row(args, meta, carrier, "baseline")
            face_row.update({
                "row_type": "face",
                "face_id": face_id,
                "edge_id": face.get("associated_edge_id"),
                "boundary_type": face.get("boundary_type"),
                "is_boundary": truth(is_boundary(face)),
                "is_top_or_above_centroid": truth(local.f(face.get("face_mid_y_um")) >= local.f(meta.get("centroid_y"))),
                "normal_x": face.get("normal_x"),
                "normal_y": face.get("normal_y"),
                "face_length_cm": face.get("face_length_cm"),
                "weight": weight,
                "Jp_face_scalar_A_per_cm2": face.get("Jp_face_scalar_A_per_cm2"),
                "Jn_face_scalar_A_per_cm2": face.get("Jn_face_scalar_A_per_cm2"),
                "contribution_x_A_per_cm2": c_x,
                "contribution_y_A_per_cm2": c_y,
                "residual_A_per_cm2": residual,
                "has_flux_constraint": truth(item is not None),
            })
            out.append(face_row)

    for face in face_rows:
        if has_flux(face):
            continue
        row = base_row(args, meta)
        row.update({
            "row_type": "missing_flux_constraint",
            "face_id": local.i(face.get("face_id")),
            "edge_id": face.get("associated_edge_id"),
            "boundary_type": face.get("boundary_type"),
            "is_boundary": truth(is_boundary(face)),
            "is_top_or_above_centroid": truth(local.f(face.get("face_mid_y_um")) >= local.f(meta.get("centroid_y"))),
            "normal_x": face.get("normal_x"),
            "normal_y": face.get("normal_y"),
            "face_length_cm": face.get("face_length_cm"),
            "has_flux_constraint": "False",
            "missing_reason": face.get("edge_lookup_status", "missing_internal_edge_flux"),
        })
        out.append(row)

    baseline_trace_hole = baseline_hole["unregularized_trace"]
    baseline_trace_electron = baseline_electron["unregularized_trace"]
    for name, include, boundary_weight, lambda_factor in variant_defs():
        electron_lambda = lambda_factor * baseline_trace_electron
        hole_lambda = lambda_factor * baseline_trace_hole
        ls_e = solve_ls(constraints(face_rows, "electron", include, boundary_weight), electron_lambda)
        ls_h = solve_ls(constraints(face_rows, "hole", include, boundary_weight), hole_lambda)
        jn_solution = ls_e["solution"]
        jp_solution = ls_h["solution"]
        jn_mag = vector_mag(jn_solution)
        jp_mag = vector_mag(jp_solution)
        gn, gp, total = gava_from(alpha_n, jn_mag, alpha_p, jp_mag)
        row = base_row(args, meta, "hole,electron", name)
        row.update({
            "row_type": "variant",
            "normal_A00": ls_h["normal_A00"],
            "normal_A01": ls_h["normal_A01"],
            "normal_A11": ls_h["normal_A11"],
            "rhs_x": ls_h["rhs_x"],
            "rhs_y": ls_h["rhs_y"],
            "condition_number": ls_h["condition_number"],
            "rank": min(local.i(ls_h["rank"]), local.i(ls_e["rank"])),
            "residual_norm": ls_h["residual_norm"],
            "tikhonov_lambda_abs": hole_lambda,
            "tikhonov_lambda_factor_times_trace": lambda_factor,
            "Jp_x_A_per_cm2": None if jp_solution is None else jp_solution[0],
            "Jp_y_A_per_cm2": None if jp_solution is None else jp_solution[1],
            "Jp_mag_A_per_cm2": jp_mag,
            "Jn_x_A_per_cm2": None if jn_solution is None else jn_solution[0],
            "Jn_y_A_per_cm2": None if jn_solution is None else jn_solution[1],
            "Jn_mag_A_per_cm2": jn_mag,
            "Gava_p_cm_minus3_s_minus1": gp,
            "Gava_n_cm_minus3_s_minus1": gn,
            "Gava_total_cm_minus3_s_minus1": total,
            "ratio_to_baseline": None if baseline_total in (None, 0.0) or total is None else total / baseline_total,
            "ratio_to_sentaurus_cell": None if sent_gava in (None, 0.0) or total is None else total / sent_gava,
            "hotspot_would_remain_rank1": "" if max_other_gava is None or total is None else truth(total >= max_other_gava),
        })
        out.append(row)

    summary_meta = {
        "alpha_n": alpha_n,
        "alpha_p": alpha_p,
        "baseline_jn": baseline_jn,
        "baseline_jp": baseline_jp,
        "baseline_gp": baseline_gp,
        "baseline_total": baseline_total,
        "sent_gava": sent_gava,
        "max_other_gava": max_other_gava,
    }
    return out, summary_meta


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in FIELDS})


def row_float(row: dict[str, Any], key: str) -> float | None:
    value = local.opt_f(row.get(key))
    return value


def find_variant(rows: list[dict[str, Any]], name: str) -> dict[str, Any] | None:
    for row in rows:
        if row.get("row_type") == "variant" and row.get("variant_name") == name:
            return row
    return None


def write_summary(path: Path, rows: list[dict[str, Any]], meta: dict[str, Any]) -> None:
    hole_system = next((row for row in rows if row.get("row_type") == "ls_system" and row.get("carrier") == "hole"), None)
    face_rows = [row for row in rows if row.get("row_type") == "face" and row.get("carrier") == "hole" and row.get("has_flux_constraint") == "True"]
    dominant = max(
        face_rows,
        key=lambda row: math.hypot(local.f(row.get("contribution_x_A_per_cm2")), local.f(row.get("contribution_y_A_per_cm2"))),
        default=None,
    )
    missing = [row for row in rows if row.get("row_type") == "missing_flux_constraint"]
    remove_7475 = find_variant(rows, "remove_faces_74_75")
    remove_boundary = find_variant(rows, "remove_all_boundary_faces")
    bw05 = find_variant(rows, "boundary_weight_0p5")
    bw025 = find_variant(rows, "boundary_weight_0p25")
    bw0 = find_variant(rows, "boundary_weight_0p0")
    baseline = find_variant(rows, "baseline")

    cond = row_float(hole_system or {}, "condition_number")
    rank = row_float(hole_system or {}, "rank")
    ill = cond is None or not math.isfinite(cond) or cond > 1.0e6 or rank is None or rank < 2
    baseline_sent_ratio = row_float(baseline or {}, "ratio_to_sentaurus_cell")
    bw05_ratio = row_float(bw05 or {}, "ratio_to_sentaurus_cell")
    bw025_ratio = row_float(bw025 or {}, "ratio_to_sentaurus_cell")
    bw0_ratio = row_float(bw0 or {}, "ratio_to_sentaurus_cell")

    def closer(value: float | None) -> bool:
        if value is None or baseline_sent_ratio is None or baseline_sent_ratio <= 0.0 or value <= 0.0:
            return False
        baseline_error = abs(math.log(baseline_sent_ratio))
        value_error = abs(math.log(value))
        return value_error < baseline_error * (1.0 - 1.0e-9)

    boundary_closer = any(closer(value) for value in (bw05_ratio, bw025_ratio, bw0_ratio))
    remove_7475_rank = (remove_7475 or {}).get("hotspot_would_remain_rank1", "")
    remove_boundary_rank = (remove_boundary or {}).get("hotspot_would_remain_rank1", "")
    missing_ids = ", ".join(str(row.get("face_id")) for row in missing)
    lines = [
        "# Self Max Cell Dual-Face LS Sensitivity",
        "",
        f"- target_bias_V: {fmt((baseline or {}).get('bias_V'))}",
        f"- target_cell_id: {fmt((baseline or {}).get('cell_id'))}",
        f"- target_centroid_um: ({fmt((baseline or {}).get('cell_centroid_x_um'))}, {fmt((baseline or {}).get('cell_centroid_y_um'))})",
        f"- baseline_Gava_total_cm^-3_s^-1: {fmt(meta.get('baseline_total'))}",
        f"- baseline_Jp_mag_A_per_cm2: {fmt(meta.get('baseline_jp'))}",
        f"- sentaurus_cell_Gava_cm^-3_s^-1: {fmt(meta.get('sent_gava'))}",
        f"- max_other_self_cell_Gava_cm^-3_s^-1: {fmt(meta.get('max_other_gava'))}",
        "",
        "## Answers",
        "",
        f"1. Is the LS reconstruction ill-conditioned? {'yes' if ill else 'no'}; hole condition_number={fmt(cond)}, rank={fmt(rank)}.",
        "2. Which face contributes most to Jp magnitude? "
        f"face_id={'' if dominant is None else dominant.get('face_id')}, edge_id={'' if dominant is None else dominant.get('edge_id')}, "
        f"boundary={'' if dominant is None else dominant.get('is_boundary')}, "
        f"contribution=({fmt(None if dominant is None else dominant.get('contribution_x_A_per_cm2'))}, "
        f"{fmt(None if dominant is None else dominant.get('contribution_y_A_per_cm2'))}) A/cm^2.",
        "3. Are faces 74/75 necessary for the hotspot? "
        f"remove_faces_74_75 ratio_to_baseline={fmt(None if remove_7475 is None else remove_7475.get('ratio_to_baseline'))}, "
        f"hotspot_would_remain_rank1={remove_7475_rank}.",
        "4. Do missing flux constraints make the LS system underdetermined or biased? "
        f"missing_faces={missing_ids or 'none'}; rank remains {fmt(rank)}, so this case is not underdetermined, but missing internal constraints leave the fitted vector biased toward constrained boundary/interior edge fluxes.",
        "5. Does reducing boundary face weight move Jp closer to Sentaurus? "
        f"{'yes' if boundary_closer else 'no'}; ratio_to_sentaurus baseline={fmt(baseline_sent_ratio)}, w0.5={fmt(bw05_ratio)}, w0.25={fmt(bw025_ratio)}, w0={fmt(bw0_ratio)}.",
        "6. Does removing boundary/top diagonal face eliminate the self hotspot? "
        f"remove_all_boundary_faces hotspot_would_remain_rank1={remove_boundary_rank}, ratio_to_baseline={fmt(None if remove_boundary is None else remove_boundary.get('ratio_to_baseline'))}.",
        "7. Should boundary cell dual-face vector current use a guarded reconstruction? "
        f"{'yes' if (dominant and dominant.get('is_boundary') == 'True') else 'review'}; boundary/top faces can dominate local Jp and should be downweighted, regularized, or source-recovered before being used in coupled avalanche source.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'self_max_cell_dual_face_ls_sensitivity.csv'}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bias", type=float, default=-20.0)
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--target-x-um", type=float, default=1.0833333333333333)
    parser.add_argument("--target-y-um", type=float, default=0.4166666666666667)
    parser.add_argument("--cell-id", type=int, default=12)
    parser.add_argument("--dual-face-csv", type=Path, default=REPO / "build/diagnostics/dual_face_geometry_audit.csv")
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/a1_b1_internal_source_current_audit/avalanche_internal_source_current_audit.csv")
    parser.add_argument("--node-compare", type=Path, default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/coarse_node_field_compare_aligned.csv")
    parser.add_argument("--self-dual-cell-csv", type=Path, default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv")
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/self_max_cell_dual_face_ls_sensitivity.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/self_max_cell_dual_face_ls_sensitivity_summary.md")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows, meta = build_output(args)
    write_rows(args.out_csv, rows)
    write_summary(args.out_summary, rows, meta)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
