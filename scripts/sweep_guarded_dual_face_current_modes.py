#!/usr/bin/env python3
"""Sweep guarded current-magnitude modes for true-dual avalanche post-processing."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
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
MODES = (
    "edge_scalar_abs",
    "dual_face_vector_mag",
    "guarded_fallback_edge_scalar",
    "guarded_blend_beta_0.25",
    "guarded_blend_beta_0.5",
    "guarded_ratio_limiter_2",
    "guarded_ratio_limiter_3",
    "guarded_ratio_limiter_5",
)


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


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{parsed:.12g}" if math.isfinite(parsed) else ""


def safe_ratio(value: float | None, ref: float | None) -> float | None:
    if value is None or ref is None or abs(ref) <= 1.0e-300:
        return None
    ratio = value / ref
    return ratio if math.isfinite(ratio) else None


def sign(value: float) -> float:
    return -1.0 if value < 0.0 else 1.0


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [clean(row) for row in csv.DictReader(handle, skipinitialspace=True)]


def is_boundary(boundary_type: str) -> bool:
    return "boundary" in boundary_type.lower()


def windows_for_x(x_um: float) -> list[str]:
    return [name for name, (lo, hi) in WINDOWS.items() if lo <= x_um <= hi]


def solve_min_norm(constraints: list[tuple[float, float, float, float]]) -> tuple[float, float] | None:
    a00 = a01 = a11 = b0 = b1 = 0.0
    for nx, ny, scalar, weight in constraints:
        a00 += weight * nx * nx
        a01 += weight * nx * ny
        a11 += weight * ny * ny
        b0 += weight * nx * scalar
        b1 += weight * ny * scalar
    det = a00 * a11 - a01 * a01
    scale = max(abs(a00 * a11), abs(a01 * a01), 1.0)
    if abs(det) > 1.0e-24 * scale:
        return ((a11 * b0 - a01 * b1) / det, (-a01 * b0 + a00 * b1) / det)
    trace = a00 + a11
    disc = max((a00 - a11) * (a00 - a11) + 4.0 * a01 * a01, 0.0)
    eig = 0.5 * (trace + math.sqrt(disc))
    if eig <= 0.0:
        return None
    if abs(a01) > abs(a00 - eig):
        ux, uy = a01, eig - a00
    else:
        ux, uy = eig - a11, a01
    norm = math.hypot(ux, uy)
    if norm <= 0.0:
        ux, uy = (1.0, 0.0) if abs(a00) >= abs(a11) else (0.0, 1.0)
    else:
        ux, uy = ux / norm, uy / norm
    projection = (ux * b0 + uy * b1) / eig
    return projection * ux, projection * uy


def vector_contributions(constraints: list[dict[str, Any]]) -> dict[int, tuple[float, float]]:
    tuples = [(row["nx"], row["ny"], row["scalar"], row["weight"]) for row in constraints]
    solution = solve_min_norm(tuples)
    if solution is None:
        return {}
    a00 = a01 = a11 = 0.0
    for nx, ny, _, weight in tuples:
        a00 += weight * nx * nx
        a01 += weight * nx * ny
        a11 += weight * ny * ny
    det = a00 * a11 - a01 * a01
    if abs(det) <= 1.0e-24 * max(abs(a00 * a11), abs(a01 * a01), 1.0):
        return {}
    i00, i01, i11 = a11 / det, -a01 / det, a00 / det
    out: dict[int, tuple[float, float]] = {}
    for row in constraints:
        bx = row["weight"] * row["nx"] * row["scalar"]
        by = row["weight"] * row["ny"] * row["scalar"]
        out[row["face_id"]] = (i00 * bx + i01 * by, i01 * bx + i11 * by)
    return out


def load_cell_rows(path: Path, mode_filter: str | None = None) -> dict[tuple[float, int], dict[str, Any]]:
    cells: dict[tuple[float, int], dict[str, Any]] = {}
    for row in read_rows(path):
        if row.get("reconstruction_status", "ok") != "ok":
            continue
        if mode_filter is not None and row.get("mode") != mode_filter:
            continue
        bias = f(row.get("bias_V"))
        cell_id = i(row.get("cell_id"))
        gava = f(row.get("Gava_cell_cm_minus3_s_minus1"))
        qg = f(row.get("qG_cell_contribution_A_per_um"))
        sent_gava = f(row.get("Gava_sent_cell_cm_minus3_s_minus1"))
        sent_qg = f(row.get("qG_sent_cell_contribution_A_per_um"))
        area_factor = safe_ratio(qg, gava) or safe_ratio(sent_qg, sent_gava) or 0.0
        cells[(bias, cell_id)] = {
            "bias": bias,
            "cell_id": cell_id,
            "x": f(row.get("centroid_x_um")),
            "y": f(row.get("centroid_y_um")),
            "jn": f(row.get("Jn_mag_A_per_cm2")),
            "jp": f(row.get("Jp_mag_A_per_cm2")),
            "alpha_n": f(row.get("alpha_n_cell_cm_inv")),
            "alpha_p": f(row.get("alpha_p_cell_cm_inv")),
            "gava": gava,
            "qg": qg,
            "sent_gava": sent_gava,
            "sent_qg": sent_qg,
            "area_factor": area_factor,
            "window_label": row.get("window_label", ""),
        }
    return cells


def load_faces(path: Path) -> dict[int, list[dict[str, Any]]]:
    by_cell: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in read_rows(path):
        face = {
            "face_id": i(row.get("dual_face_id")),
            "cell_id": i(row.get("cell_id")),
            "edge_id": i(row.get("associated_primal_edge_id")),
            "nx": f(row.get("dual_face_normal_x")),
            "ny": f(row.get("dual_face_normal_y")),
            "length": f(row.get("dual_face_length_cm"), 1.0),
            "orientation": f(row.get("orientation_sign"), 1.0),
            "boundary_type": row.get("boundary_type", ""),
        }
        by_cell[face["cell_id"]].append(face)
    return by_cell


def load_internal(path: Path) -> dict[tuple[float, int], dict[str, float]]:
    edges: dict[tuple[float, int], dict[str, float]] = {}
    for row in read_rows(path):
        if row.get("source_location_type", "edge") != "edge":
            continue
        bias = f(row.get("bias_V"))
        edge_id = i(row.get("source_entity_id"))
        edges[(bias, edge_id)] = {
            "jn": f(row.get("Jn_mag_used_A_per_cm2")),
            "jp": f(row.get("Jp_mag_used_A_per_cm2")),
            "e_sign": sign(f(row.get("electron_raw_signed_flux_proxy", "1"))),
            "h_sign": sign(f(row.get("hole_raw_signed_flux_proxy", "1"))),
        }
    return edges


def load_local_ls_flags(path: Path) -> dict[tuple[float, int], bool]:
    flags: dict[tuple[float, int], bool] = {}
    for row in read_rows(path):
        if row.get("row_type") != "variant":
            continue
        if row.get("variant_name") != "remove_faces_74_75":
            continue
        ratio = opt_f(row.get("ratio_to_baseline"))
        rank1 = row.get("hotspot_would_remain_rank1", "").lower() == "true"
        if ratio is not None and ratio < 0.5 and not rank1:
            flags[(f(row.get("bias_V")), i(row.get("cell_id")))] = True
    return flags


def cell_flag_metrics(
    bias: float,
    cell_id: int,
    faces_by_cell: dict[int, list[dict[str, Any]]],
    internal: dict[tuple[float, int], dict[str, float]],
) -> dict[str, Any]:
    faces = faces_by_cell.get(cell_id, [])
    constraints: list[dict[str, Any]] = []
    missing = 0
    boundary_count = 0
    for face in faces:
        edge = internal.get((bias, face["edge_id"]))
        if edge is None:
            missing += 1
            continue
        boundary = is_boundary(face["boundary_type"])
        if boundary:
            boundary_count += 1
        constraints.append({
            "face_id": face["face_id"],
            "nx": face["nx"],
            "ny": face["ny"],
            "scalar": face["orientation"] * edge["h_sign"] * edge["jp"],
            "weight": face["length"] if face["length"] > 0.0 else 1.0,
            "boundary": boundary,
        })
    solution = solve_min_norm([(row["nx"], row["ny"], row["scalar"], row["weight"]) for row in constraints])
    base_mag = math.hypot(*solution) if solution is not None else 0.0
    contributions = vector_contributions(constraints)
    contribution_mags = {face_id: math.hypot(*vec) for face_id, vec in contributions.items()}
    total_contrib = sum(contribution_mags.values())
    dominant_fraction = max(contribution_mags.values(), default=0.0) / total_contrib if total_contrib > 0.0 else 0.0
    boundary_fraction = (
        sum(contribution_mags.get(row["face_id"], 0.0) for row in constraints if row["boundary"]) / total_contrib
        if total_contrib > 0.0 else 0.0
    )
    loo_ratio = 1.0
    for row in constraints:
        trial = [item for item in constraints if item["face_id"] != row["face_id"]]
        trial_solution = solve_min_norm([(item["nx"], item["ny"], item["scalar"], item["weight"]) for item in trial])
        trial_mag = math.hypot(*trial_solution) if trial_solution is not None else 0.0
        if base_mag > 0.0 and trial_mag > 0.0:
            loo_ratio = max(loo_ratio, base_mag / trial_mag, trial_mag / base_mag)
        elif base_mag > 0.0:
            loo_ratio = math.inf
    no_boundary = [row for row in constraints if not row["boundary"]]
    no_boundary_solution = solve_min_norm([(row["nx"], row["ny"], row["scalar"], row["weight"]) for row in no_boundary])
    no_boundary_mag = math.hypot(*no_boundary_solution) if no_boundary_solution is not None else 0.0
    remove_boundary_change = abs(base_mag - no_boundary_mag) / base_mag if base_mag > 0.0 else 0.0
    flagged = (
        missing > 0
        or (boundary_count > 0 and dominant_fraction > 0.4)
        or boundary_fraction > 0.5
        or loo_ratio > 2.0
        or remove_boundary_change > 0.5
    )
    return {
        "missing_flux_constraints": missing,
        "boundary_face_count": boundary_count,
        "dominant_face_contribution_fraction": dominant_fraction,
        "boundary_contribution_fraction": boundary_fraction,
        "leave_one_out_max_change_ratio": loo_ratio,
        "remove_boundary_change_fraction": remove_boundary_change,
        "flagged": flagged,
    }


def choose_j(edge_value: float, dual_value: float, flagged: bool, mode: str) -> float:
    if mode == "edge_scalar_abs":
        return edge_value
    if mode == "dual_face_vector_mag" or not flagged:
        return dual_value
    if mode == "guarded_fallback_edge_scalar":
        return edge_value
    if mode == "guarded_blend_beta_0.25":
        return 0.75 * edge_value + 0.25 * dual_value
    if mode == "guarded_blend_beta_0.5":
        return 0.5 * edge_value + 0.5 * dual_value
    if mode == "guarded_ratio_limiter_2":
        return min(dual_value, 2.0 * edge_value)
    if mode == "guarded_ratio_limiter_3":
        return min(dual_value, 3.0 * edge_value)
    if mode == "guarded_ratio_limiter_5":
        return min(dual_value, 5.0 * edge_value)
    raise ValueError(mode)


def build_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    dual_cells = load_cell_rows(args.dual_cell_csv)
    edge_cells = load_cell_rows(args.edge_cell_csv, "current_edge_scalar_existing")
    faces_by_cell = load_faces(args.dual_face_csv)
    internal = load_internal(args.internal_audit)
    local_flags = load_local_ls_flags(args.ls_sensitivity_csv)
    biases = sorted({bias for bias, _ in dual_cells})
    output: list[dict[str, Any]] = []

    flag_metrics = {
        key: cell_flag_metrics(key[0], key[1], faces_by_cell, internal)
        for key in dual_cells
    }
    for key, flagged in local_flags.items():
        flag_metrics.setdefault(key, {})["flagged"] = flag_metrics.get(key, {}).get("flagged", False) or flagged

    for bias in biases:
        sent_candidates = [cell for (row_bias, _), cell in dual_cells.items() if abs(row_bias - bias) <= 1.0e-9]
        sent_max = max(sent_candidates, key=lambda cell: cell["sent_gava"], default=None)
        for mode in MODES:
            mode_cells: list[dict[str, Any]] = []
            for key, dual in dual_cells.items():
                if abs(key[0] - bias) > 1.0e-9:
                    continue
                edge = edge_cells.get(key, dual)
                metrics = flag_metrics.get(key, {})
                flagged = bool(metrics.get("flagged", False))
                jn = choose_j(edge["jn"], dual["jn"], flagged, mode)
                jp = choose_j(edge["jp"], dual["jp"], flagged, mode)
                if mode == "edge_scalar_abs":
                    gava = edge["gava"]
                    qg = edge["qg"]
                elif mode == "dual_face_vector_mag":
                    gava = dual["gava"]
                    qg = dual["qg"]
                elif flagged and mode == "guarded_fallback_edge_scalar":
                    gava = edge["gava"]
                    qg = edge["qg"]
                else:
                    gava = (dual["alpha_n"] * jn + dual["alpha_p"] * jp) / Q_C
                    qg = gava * dual["area_factor"]
                mode_cell = {
                    **dual,
                    "mode": mode,
                    "jn_mode": jn,
                    "jp_mode": jp,
                    "gava_mode": gava,
                    "qg_mode": qg,
                    "flagged": flagged,
                    **metrics,
                }
                mode_cells.append(mode_cell)
            buckets = {name: {"self": 0.0, "sent": 0.0} for name in WINDOWS}
            for cell in mode_cells:
                for window in windows_for_x(cell["x"]):
                    buckets[window]["self"] += cell["qg_mode"]
                    buckets[window]["sent"] += cell["sent_qg"]
            max_cell = max(mode_cells, key=lambda cell: cell["gava_mode"], default=None)
            cell12 = next((cell for cell in mode_cells if cell["cell_id"] == args.target_cell_id), None)
            distance = None
            if max_cell is not None and sent_max is not None:
                distance = math.hypot(max_cell["x"] - sent_max["x"], max_cell["y"] - sent_max["y"])
            row = {
                "bias_V": bias,
                "current_magnitude_mode": mode,
                "qG_full_ratio": safe_ratio(buckets["full"]["self"], buckets["full"]["sent"]),
                "qG_junction_ratio": safe_ratio(buckets["junction"]["self"], buckets["junction"]["sent"]),
                "qG_center_ratio": safe_ratio(buckets["center"]["self"], buckets["center"]["sent"]),
                "qG_left_shoulder_ratio": safe_ratio(buckets["left_shoulder"]["self"], buckets["left_shoulder"]["sent"]),
                "qG_right_shoulder_ratio": safe_ratio(buckets["right_shoulder"]["self"], buckets["right_shoulder"]["sent"]),
                "max_Gava_cm_minus3_s_minus1": None if max_cell is None else max_cell["gava_mode"],
                "max_Gava_x_um": None if max_cell is None else max_cell["x"],
                "max_Gava_y_um": None if max_cell is None else max_cell["y"],
                "max_Gava_distance_to_sentaurus_um": distance,
                "self_max_cell_id": None if max_cell is None else max_cell["cell_id"],
                "cell12_remains_rank1": "" if cell12 is None or max_cell is None else str(max_cell["cell_id"] == args.target_cell_id),
                "cell12_Jp_A_per_cm2": None if cell12 is None else cell12["jp_mode"],
                "cell12_Gava_cm_minus3_s_minus1": None if cell12 is None else cell12["gava_mode"],
                "flagged_cell_count": sum(1 for cell in mode_cells if cell["flagged"]),
                "cell12_flagged": "" if cell12 is None else str(cell12["flagged"]),
                "cell12_missing_flux_constraints": None if cell12 is None else cell12.get("missing_flux_constraints"),
                "cell12_boundary_face_count": None if cell12 is None else cell12.get("boundary_face_count"),
                "cell12_boundary_contribution_fraction": None if cell12 is None else cell12.get("boundary_contribution_fraction"),
                "cell12_remove_boundary_change_fraction": None if cell12 is None else cell12.get("remove_boundary_change_fraction"),
            }
            output.append(row)
    return output


FIELDS = [
    "bias_V",
    "current_magnitude_mode",
    "qG_full_ratio",
    "qG_junction_ratio",
    "qG_center_ratio",
    "qG_left_shoulder_ratio",
    "qG_right_shoulder_ratio",
    "max_Gava_cm_minus3_s_minus1",
    "max_Gava_x_um",
    "max_Gava_y_um",
    "max_Gava_distance_to_sentaurus_um",
    "self_max_cell_id",
    "cell12_remains_rank1",
    "cell12_Jp_A_per_cm2",
    "cell12_Gava_cm_minus3_s_minus1",
    "flagged_cell_count",
    "cell12_flagged",
    "cell12_missing_flux_constraints",
    "cell12_boundary_face_count",
    "cell12_boundary_contribution_fraction",
    "cell12_remove_boundary_change_fraction",
]


def write_csv_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in FIELDS})


def score_ratio(value: float | None) -> float:
    if value is None or value <= 0.0:
        return math.inf
    return abs(math.log(value))


def row_value(row: dict[str, Any], key: str) -> float | None:
    return opt_f(row.get(key))


def rows_at(rows: list[dict[str, Any]], bias: float = -20.0) -> list[dict[str, Any]]:
    return [row for row in rows if abs(f(row.get("bias_V")) - bias) <= 1.0e-9]


def best(rows: list[dict[str, Any]], key_fn: Any) -> dict[str, Any] | None:
    return min(rows, key=key_fn, default=None)


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    r20 = rows_at(rows)
    guarded = [row for row in r20 if str(row.get("current_magnitude_mode", "")).startswith("guarded_")]
    edge = next((row for row in r20 if row.get("current_magnitude_mode") == "edge_scalar_abs"), None)
    dual = next((row for row in r20 if row.get("current_magnitude_mode") == "dual_face_vector_mag"), None)
    best_hotspot = best(guarded, lambda row: row_value(row, "cell12_Gava_cm_minus3_s_minus1") or math.inf)
    best_full_junction = best(
        r20,
        lambda row: score_ratio(row_value(row, "qG_full_ratio")) + score_ratio(row_value(row, "qG_junction_ratio")),
    )
    no_hotspot_guarded = [row for row in guarded if row.get("cell12_remains_rank1") == "False"]
    best_right = best(no_hotspot_guarded, lambda row: score_ratio(row_value(row, "qG_right_shoulder_ratio")))
    best_right_any_guarded = best(guarded, lambda row: score_ratio(row_value(row, "qG_right_shoulder_ratio")))
    edge_full = row_value(edge or {}, "qG_full_ratio")
    improving_guarded = [
        row for row in no_hotspot_guarded
        if edge_full is not None and (row_value(row, "qG_full_ratio") or -math.inf) > edge_full
    ]
    recommended = best(
        improving_guarded,
        lambda row: score_ratio(row_value(row, "qG_full_ratio")) + score_ratio(row_value(row, "qG_junction_ratio")),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Guarded Dual-Face Current Mode Sweep",
        "",
        "- diagnostic: post-processing only",
        "- default coupled solver changed: no",
        "- current unit: A/cm^2",
        "- Gava unit: cm^-3 s^-1",
        "",
        "## -20V Mode Table",
        "",
        "| mode | full | junction | center | right_shoulder | max_cell | cell12_rank1 | cell12_Gava |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in r20:
        lines.append(
            f"| {row.get('current_magnitude_mode')} | {fmt(row.get('qG_full_ratio'))} | "
            f"{fmt(row.get('qG_junction_ratio'))} | {fmt(row.get('qG_center_ratio'))} | "
            f"{fmt(row.get('qG_right_shoulder_ratio'))} | {fmt(row.get('self_max_cell_id'))} | "
            f"{row.get('cell12_remains_rank1')} | {fmt(row.get('cell12_Gava_cm_minus3_s_minus1'))} |"
        )
    lines.extend([
        "",
        "## Answers",
        "",
        "1. Which guarded mode best reduces the artificial hotspot at cell 12? "
        f"{'' if best_hotspot is None else best_hotspot.get('current_magnitude_mode')}; "
        f"cell12_Gava={fmt(None if best_hotspot is None else best_hotspot.get('cell12_Gava_cm_minus3_s_minus1'))}.",
        "2. Which mode best matches -20V full/junction qG? "
        f"{'' if best_full_junction is None else best_full_junction.get('current_magnitude_mode')}; "
        f"full={fmt(None if best_full_junction is None else best_full_junction.get('qG_full_ratio'))}, "
        f"junction={fmt(None if best_full_junction is None else best_full_junction.get('qG_junction_ratio'))}.",
        "3. Which mode improves right_shoulder without creating a new hotspot? "
        f"{'none' if best_right is None else best_right.get('current_magnitude_mode')}; "
        f"best_guarded_right_shoulder={'' if best_right_any_guarded is None else best_right_any_guarded.get('current_magnitude_mode')} "
        f"({fmt(None if best_right_any_guarded is None else best_right_any_guarded.get('qG_right_shoulder_ratio'))}) "
        f"still has cell12_rank1={'' if best_right_any_guarded is None else best_right_any_guarded.get('cell12_remains_rank1')}.",
        "4. Does any guarded mode keep qG improvement over edge_scalar_abs while avoiding dual_face_vector_mag instability? "
        f"{'yes' if recommended is not None else 'no'}"
        f"{'' if recommended is None else '; candidate=' + str(recommended.get('current_magnitude_mode'))}.",
        "5. Should a guarded mode be tested in coupled BV? "
        f"{'yes, as diagnostic candidate ' + str(recommended.get('current_magnitude_mode')) if recommended is not None else 'not yet; current guarded modes trade off qG recovery too strongly or leave instability.'}",
        "6. Do not recommend unguarded dual_face_vector_mag as default. "
        f"Confirmed; unguarded cell12_rank1={'' if dual is None else dual.get('cell12_remains_rank1')}, "
        f"edge_scalar_full={fmt(None if edge is None else edge.get('qG_full_ratio'))}.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'guarded_dual_face_current_mode_sweep.csv'}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dual-cell-csv", type=Path, default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv")
    parser.add_argument("--edge-cell-csv", type=Path, default=REPO / "build/diagnostics/edge_to_cell_vector_current_reconstruction.csv")
    parser.add_argument("--dual-face-csv", type=Path, default=REPO / "build/diagnostics/dual_face_geometry_audit.csv")
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/a1_b1_internal_source_current_audit/avalanche_internal_source_current_audit.csv")
    parser.add_argument("--ls-sensitivity-csv", type=Path, default=REPO / "build/diagnostics/self_max_cell_dual_face_ls_sensitivity.csv")
    parser.add_argument("--target-cell-id", type=int, default=12)
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/guarded_dual_face_current_mode_sweep.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/guarded_dual_face_current_mode_sweep_summary.md")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = build_rows(args)
    write_csv_rows(args.out_csv, rows)
    write_summary(args.out_summary, rows)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
