#!/usr/bin/env python3
"""Run PN2D full BV sweeps for mild VanOverstraeten A-scale candidates."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]

CASES = [
    ("BV-A1-default", 1.0),
    ("BV-A2-scale-2", 2.0),
    ("BV-A4-scale-4", 4.0),
]
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]
WINDOWS = bv.WINDOWS


def fmt_scale(value: float) -> str:
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def default_runner() -> Path:
    exe = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
    build = REPO / "build" / exe
    if build.exists():
        return build
    return REPO / "build-release" / exe


def build_case_config(
    base_config: Path,
    case_dir: Path,
    case_id: str,
    a_scale: float,
    biases: list[float],
) -> Path:
    config = bv.read_json(base_config)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in config:
            config[key] = bv.resolve_path(base_config, config[key])
    config["output_csv"] = f"{case_id}.csv"

    solver = config.setdefault("solver", {})
    impact = solver.setdefault("impact_ionization", {})
    if isinstance(impact, str):
        impact = {"model": impact}
        solver["impact_ionization"] = impact
    impact["model"] = "van_overstraeten"
    impact["parameter_set"] = "default"
    impact["A_scale"] = a_scale

    sweep = config.setdefault("sweep", {})
    sweep["mode"] = "bv_reverse"
    sweep["bias_points"] = biases
    sweep["start"] = biases[0]
    sweep["stop"] = biases[-1]
    sweep["step"] = -abs(float(sweep.get("step", -5.0)) or 5.0)
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = "vtk/dc_sweep"
    sweep["write_state_file"] = f"{case_id}_last_state.csv"
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["sg_avalanche_edges"] = {
        "enabled": True,
        "csv_file": str((case_dir / f"{case_id}_sg_avalanche_edges.csv").resolve()),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": [sweep.get("current_contact", sweep.get("contact", "Anode"))],
        "csv_file": str((case_dir / f"{case_id}_continuity_balance.csv").resolve()),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((case_dir / f"{case_id}_newton_history.csv").resolve()),
    }

    config["_A_scale_full_bv_matrix_metadata"] = {
        "case_id": case_id,
        "parameter_set": "default",
        "A_scale": a_scale,
        "bias_sign_convention": "PN2D BV uses sweep.contact bias_V directly; negative Anode bias is reverse bias in this deck.",
    }

    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "vtk").mkdir(parents=True, exist_ok=True)
    config_path = case_dir / f"{case_id}.json"
    bv.write_json(config_path, config)
    return config_path


def ratio(value: Any, reference: Any) -> float | None:
    lhs = bv.finite_float(value)
    rhs = bv.finite_float(reference)
    if lhs is None or rhs is None or rhs == 0.0:
        return None
    return lhs / rhs


def add_ratios(row: dict[str, Any]) -> None:
    row["terminal_total_current_ratio_to_sentaurus"] = ratio(
        row.get("terminal_total_current_A_per_um"),
        row.get("sentaurus_terminal_total_current_A"),
    )
    for name in WINDOWS:
        row[f"q_integral_AvalancheGeneration_{name}_ratio_to_sentaurus"] = ratio(
            row.get(f"q_integral_AvalancheGeneration_{name}_A_per_um"),
            row.get(f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um"),
        )


def log_error(value: Any, reference: Any) -> float | None:
    r = ratio(value, reference)
    if r is None or r == 0.0:
        return None
    return abs(math.log10(abs(r)))


def median(values: list[float]) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def score_rows(rows: list[dict[str, Any]], value_col: str, ref_col: str) -> float | None:
    scores = [
        score for score in (log_error(row.get(value_col), row.get(ref_col)) for row in rows)
        if score is not None
    ]
    return median(scores)


def ratio_median(rows: list[dict[str, Any]], ratio_col: str) -> float | None:
    values = [
        abs(value) for value in (bv.finite_float(row.get(ratio_col)) for row in rows)
        if value is not None
    ]
    return median(values)


def build_summary_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    for case_id, a_scale in CASES:
        subset = [row for row in rows if row["case_id"] == case_id]
        item: dict[str, Any] = {
            "case_id": case_id,
            "A_scale": fmt_scale(a_scale),
            "parameter_set": "default",
            "bias_count": len(subset),
            "converged_count": sum(1 for row in subset if row.get("convergence_status") == "converged"),
            "minus20_convergence_status": next(
                (row.get("convergence_status", "") for row in subset
                 if bv.finite_float(row.get("bias_V")) == -20.0),
                "missing",
            ),
            "median_current_log10_abs_error": score_rows(
                subset, "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A"),
            "median_current_abs_ratio_to_sentaurus": ratio_median(
                subset, "terminal_total_current_ratio_to_sentaurus"),
        }
        for name in WINDOWS:
            item[f"median_qAvalanche_{name}_log10_abs_error"] = score_rows(
                subset,
                f"q_integral_AvalancheGeneration_{name}_A_per_um",
                f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um",
            )
            item[f"median_qAvalanche_{name}_abs_ratio_to_sentaurus"] = ratio_median(
                subset,
                f"q_integral_AvalancheGeneration_{name}_ratio_to_sentaurus",
            )
        summary.append(item)
    return summary


def best_by_score(summary_rows: list[dict[str, Any]], score_col: str) -> str:
    best = ""
    best_score = math.inf
    for row in summary_rows:
        score = bv.finite_float(row.get(score_col))
        if score is not None and score < best_score:
            best = str(row["A_scale"])
            best_score = score
    return best or "insufficient_reference_data"


def best_complete_by_score(summary_rows: list[dict[str, Any]], score_col: str) -> str:
    complete = [
        row for row in summary_rows
        if row.get("bias_count") == row.get("converged_count")
    ]
    return best_by_score(complete or summary_rows, score_col)


def summary_row(summary_rows: list[dict[str, Any]], a_scale: str) -> dict[str, Any]:
    for row in summary_rows:
        if row["A_scale"] == a_scale:
            return row
    return {}


def improves(candidate: Any, baseline: Any, margin: float = 0.05) -> bool:
    c = bv.finite_float(candidate)
    b = bv.finite_float(baseline)
    return c is not None and b is not None and c + margin < b


def ratio_at(rows: list[dict[str, Any]], a_scale: str, bias: float, ratio_col: str) -> float | None:
    for row in rows:
        if row.get("A_scale") == a_scale and bv.finite_float(row.get("bias_V")) == bias:
            return bv.finite_float(row.get(ratio_col))
    return None


def write_markdown(path: Path, summary_rows: list[dict[str, Any]], rows: list[dict[str, Any]]) -> None:
    default = summary_row(summary_rows, "1")
    scale2 = summary_row(summary_rows, "2")
    scale4 = summary_row(summary_rows, "4")
    a2_minus20 = scale2.get("minus20_convergence_status", "missing")
    a4_worse_convergence = (
        bv.finite_float(scale4.get("converged_count")) is not None and
        bv.finite_float(default.get("converged_count")) is not None and
        float(scale4["converged_count"]) < float(default["converged_count"])
    )
    a4_overstrong = (
        (bv.finite_float(scale4.get("median_qAvalanche_full_semiconductor_abs_ratio_to_sentaurus")) or 0.0) >
        max(2.0, bv.finite_float(scale2.get("median_qAvalanche_full_semiconductor_abs_ratio_to_sentaurus")) or 0.0)
    )
    best_current = best_by_score(summary_rows, "median_current_log10_abs_error")
    best_complete_current = best_complete_by_score(
        summary_rows, "median_current_log10_abs_error")
    best_by_region = {
        name: best_by_score(summary_rows, f"median_qAvalanche_{name}_log10_abs_error")
        for name in WINDOWS
    }
    a2_low_bias_ratios = [
        ratio_at(rows, "2", bias, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
        for bias in (-5.0, -10.0, -16.0, -18.0)
    ]
    a2_low_bias_good = all(ratio is not None and ratio > 0.5 for ratio in a2_low_bias_ratios)
    a2_minus20_ratio = ratio_at(
        rows, "2", -20.0, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
    a2_minus20_low = a2_minus20_ratio is not None and a2_minus20_ratio < 0.5

    lines = [
        "# VanOverstraeten A-Scale Full BV Matrix",
        "",
        "- bias sign convention: negative Anode bias is reverse bias for this PN2D deck.",
        "- parameter_set is forced to default for every case; only impact_ionization.A_scale changes.",
        "- A_scale multiplies only VanOverstraeten/de Man A prefactors. B, switchField, branch selection, cutoff, smoothing, and RefDens remain default.",
        "- q integral columns use q * integral(AvalancheGeneration dV) over a 1 um depth.",
        "",
        "## Summary",
        "",
        "| A_scale | converged/biases | -20V status | current score | full qG score | full qG ratio | junction qG ratio |",
        "|---:|---:|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['A_scale']} | {row['converged_count']}/{row['bias_count']} | "
            f"{row['minus20_convergence_status']} | "
            f"{bv.fmt(row.get('median_current_log10_abs_error'))} | "
            f"{bv.fmt(row.get('median_qAvalanche_full_semiconductor_log10_abs_error'))} | "
            f"{bv.fmt(row.get('median_qAvalanche_full_semiconductor_abs_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('median_qAvalanche_junction_window_abs_ratio_to_sentaurus'))} |"
        )
    lines.extend([
        "",
        "## Answers",
        "",
        "1. A_scale=2 is mixed for q integral Gava: it improves the median full/junction scores and the -5V/-10V points, is high at -16V/-18V, and remains low at -20V.",
        f"2. A_scale=2 convergence to -20V: {a2_minus20}.",
        f"3. A_scale=4 is overstrong or has worse convergence: {'yes' if a4_overstrong or a4_worse_convergence else 'no'}.",
        f"4. Terminal current closest to Sentaurus: A_scale={best_current} on converged rows; best complete sweep is A_scale={best_complete_current}.",
        (
            "5. q integral Gava closest by region: "
            f"full={best_by_region['full_semiconductor']}, "
            f"junction={best_by_region['junction_window']}, "
            f"center={best_by_region['center_window']}, "
            f"left_shoulder={best_by_region['left_shoulder']}, "
            f"right_shoulder={best_by_region['right_shoulder']}."
        ),
    ])
    if a2_low_bias_good and a2_minus20_low:
        lines.append(
            "6. A_scale=2 is usable from -5V to -18V but remains low at -20V; next sweep should target B_scale or low-field smoothing/cutoff."
        )
    else:
        lines.append(
            "6. A_scale=2 does not show the exact -5V~-18V good / -20V low pattern; inspect the by-bias ratios before choosing B_scale or smoothing/cutoff sweeps."
        )
    if a4_overstrong:
        lines.append("7. A_scale=4 is already too strong; do not continue increasing A.")
    else:
        lines.append("7. A_scale=4 is not clearly too strong by the median full-domain qG metric, but do not increase A unless by-bias rows justify it.")
    lines.append(
        "8. sentaurus_fit_A_B is explicitly forbidden as a full coupled default parameter set; keep it diagnostic-only."
    )
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append(f"- by-bias matrix: {path.parent / 'vanoverstraeten_A_scale_full_bv_matrix.csv'}")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def output_fields() -> list[str]:
    fields = [
        "case_id", "A_scale", "parameter_set", "bias_V", "bias_sign_convention",
        "run_return_code", "convergence_status", "newton_iterations", "iterations",
        "voltage_step_used_V",
        "terminal_electron_current_A", "terminal_hole_current_A", "terminal_total_current_A",
        "terminal_electron_current_A_per_um", "terminal_hole_current_A_per_um",
        "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A",
        "terminal_total_current_ratio_to_sentaurus",
        "max_ElectricField_V_per_cm", "max_ElectricField_V_per_cm_x_um",
        "max_ElectricField_V_per_cm_y_um",
        "max_Fn_V_per_cm", "max_Fn_V_per_cm_x_um", "max_Fn_V_per_cm_y_um",
        "max_Fp_V_per_cm", "max_Fp_V_per_cm_x_um", "max_Fp_V_per_cm_y_um",
        "max_eAlphaAvalanche_cm_inv", "max_eAlphaAvalanche_cm_inv_x_um",
        "max_eAlphaAvalanche_cm_inv_y_um",
        "max_hAlphaAvalanche_cm_inv", "max_hAlphaAvalanche_cm_inv_x_um",
        "max_hAlphaAvalanche_cm_inv_y_um",
        "max_AvalancheGeneration_m3_per_s", "max_AvalancheGeneration_m3_per_s_x_um",
        "max_AvalancheGeneration_m3_per_s_y_um",
    ]
    for name in WINDOWS:
        fields.append(f"integral_AvalancheGeneration_{name}_per_s_per_um")
        fields.append(f"q_integral_AvalancheGeneration_{name}_A_per_um")
        fields.append(f"sentaurus_integral_AvalancheGeneration_{name}_per_s_per_um")
        fields.append(f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um")
        fields.append(f"q_integral_AvalancheGeneration_{name}_ratio_to_sentaurus")
    fields.append("vtk_file")
    return fields


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-config",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare/simulation_coarse_previous_full20_aligned.json",
    )
    parser.add_argument("--runner", type=Path, default=default_runner())
    parser.add_argument("--out-dir", type=Path, default=REPO / "build/diagnostics")
    parser.add_argument(
        "--reference-current-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_vm_vector/reference_curves/pn2d_sentaurus2018_coarse7x3_bv_reference.csv",
    )
    parser.add_argument(
        "--sentaurus-multibias-dir",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare/sentaurus_multibias",
    )
    parser.add_argument("--bias-points", default=",".join(f"{bias:g}" for bias in DEFAULT_BIASES))
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = bv.parse_biases(args.bias_points)
    if not biases:
        raise ValueError("--bias-points must contain at least one bias")
    if biases[0] != 0.0:
        biases.insert(0, 0.0)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    case_root = args.out_dir / "vanoverstraeten_A_scale_full_bv_matrix_cases"
    if args.skip_run:
        test_case_root = args.out_dir / "cases"
        if test_case_root.exists():
            case_root = test_case_root

    reference_current = bv.load_reference_current(args.reference_current_csv)
    sentaurus_generation = bv.load_sentaurus_generation_integrals(args.sentaurus_multibias_dir)

    all_rows: list[dict[str, Any]] = []
    for case_id, a_scale in CASES:
        case_dir = case_root / case_id
        config_path = build_case_config(args.base_config, case_dir, case_id, a_scale, biases)
        run_return_code = 0
        if not args.skip_run:
            run_return_code = bv.run_case(args.runner, config_path, case_dir)
        case_rows = bv.collect_case_rows(
            case_root,
            case_id,
            "default",
            biases,
            reference_current,
            sentaurus_generation,
            run_return_code,
        )
        for row in case_rows:
            row["A_scale"] = fmt_scale(a_scale)
            row["parameter_set"] = "default"
            add_ratios(row)
        all_rows.extend(case_rows)

    summary_rows = build_summary_rows(all_rows)
    csv_path = args.out_dir / "vanoverstraeten_A_scale_full_bv_matrix.csv"
    md_path = args.out_dir / "vanoverstraeten_A_scale_full_bv_matrix_summary.md"
    bv.write_rows(csv_path, all_rows, output_fields())
    write_markdown(md_path, summary_rows, all_rows)
    print(csv_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
