#!/usr/bin/env python3
"""Run PN2D full BV sweeps for A_scale=2 global VanOverstraeten B-scale candidates."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
A_SCALE = 2.0
CASES = [
    ("BV-A2-B0p85", 0.85),
    ("BV-A2-B0p90", 0.90),
    ("BV-A2-B0p95", 0.95),
    ("BV-A2-B1p00", 1.00),
    ("BV-A2-B1p05", 1.05),
    ("BV-A2-B1p10", 1.10),
]
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]
WINDOWS = bv.WINDOWS


def fmt_scale(value: float) -> str:
    return f"{value:.2f}"


def fmt_a_scale(value: float) -> str:
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
    b_scale: float,
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
    impact["A_scale"] = A_SCALE
    impact["B_scale"] = b_scale

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

    config["_A2_B_scale_full_bv_matrix_metadata"] = {
        "case_id": case_id,
        "parameter_set": "default",
        "A_scale": A_SCALE,
        "B_scale": b_scale,
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
    row["actual_voltage_step_V"] = row.get("voltage_step_used_V")
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
    for case_id, b_scale in CASES:
        subset = [row for row in rows if row["case_id"] == case_id]
        item: dict[str, Any] = {
            "case_id": case_id,
            "A_scale": fmt_a_scale(A_SCALE),
            "B_scale": fmt_scale(b_scale),
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
            best = str(row["B_scale"])
            best_score = score
    return best or "insufficient_reference_data"


def row_for(summary_rows: list[dict[str, Any]], b_scale: str) -> dict[str, Any]:
    for row in summary_rows:
        if row.get("B_scale") == b_scale:
            return row
    return {}


def ratio_at(rows: list[dict[str, Any]], b_scale: str, bias: float, ratio_col: str) -> float | None:
    for row in rows:
        if row.get("B_scale") == b_scale and bv.finite_float(row.get("bias_V")) == bias:
            return bv.finite_float(row.get(ratio_col))
    return None


def current_ratio_at(rows: list[dict[str, Any]], b_scale: str, bias: float) -> float | None:
    return ratio_at(rows, b_scale, bias, "terminal_total_current_ratio_to_sentaurus")


def improves_minus20_without_hurting_mid(rows: list[dict[str, Any]], b_scale: str) -> bool:
    base20 = ratio_at(rows, "1.00", -20.0, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
    cand20 = ratio_at(rows, b_scale, -20.0, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
    if base20 is None or cand20 is None or abs(cand20) <= abs(base20):
        return False
    for bias in (-16.0, -18.0):
        base = ratio_at(rows, "1.00", bias, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
        cand = ratio_at(rows, b_scale, bias, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
        if base is not None and cand is not None and abs(cand) > 1.5 * max(abs(base), 1.0):
            return False
    return True


def all_target_ratios_reasonable(rows: list[dict[str, Any]], b_scale: str) -> bool:
    for bias in (-5.0, -10.0, -16.0, -18.0, -20.0):
        q_ratio = ratio_at(rows, b_scale, bias, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
        i_ratio = current_ratio_at(rows, b_scale, bias)
        if q_ratio is None or i_ratio is None:
            return False
        if not (0.5 <= abs(q_ratio) <= 2.0 and 0.5 <= abs(i_ratio) <= 2.0):
            return False
    return True


def write_markdown(path: Path, summary_rows: list[dict[str, Any]], rows: list[dict[str, Any]]) -> None:
    best_current = best_by_score(summary_rows, "median_current_log10_abs_error")
    best_full = best_by_score(summary_rows, "median_qAvalanche_full_semiconductor_log10_abs_error")
    best_junction = best_by_score(summary_rows, "median_qAvalanche_junction_window_log10_abs_error")
    low_scales = ["0.85", "0.90", "0.95"]
    low_improves_20 = [scale for scale in low_scales if improves_minus20_without_hurting_mid(rows, scale)]
    low_makes_mid_stronger = []
    for scale in low_scales:
        for bias in (-16.0, -18.0):
            base = ratio_at(rows, "1.00", bias, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
            cand = ratio_at(rows, scale, bias, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
            if base is not None and cand is not None and abs(cand) > abs(base):
                low_makes_mid_stronger.append(scale)
                break
    covering = [row["B_scale"] for row in summary_rows if all_target_ratios_reasonable(rows, row["B_scale"])]
    all_minus20 = [
        ratio_at(rows, row["B_scale"], -20.0, "q_integral_AvalancheGeneration_full_semiconductor_ratio_to_sentaurus")
        for row in summary_rows
    ]
    solved_minus20 = any(value is not None and 0.5 <= abs(value) <= 2.0 for value in all_minus20)
    recommended = [scale for scale in ("0.95", "0.90") if improves_minus20_without_hurting_mid(rows, scale)]

    lines = [
        "# VanOverstraeten A2 B-Scale Full BV Matrix",
        "",
        "- bias sign convention: negative Anode bias is reverse bias for this PN2D deck.",
        "- parameter_set is forced to default for every case; sentaurus_fit_A_B is not used.",
        "- A_scale is fixed at 2. B_scale multiplies only VanOverstraeten/de Man electron/hole low/high B critical fields.",
        "- A, switchField, branch selection, cutoff, smoothing, and RefDens remain default.",
        "- q integral columns use q * integral(AvalancheGeneration dV) over a 1 um depth.",
        "",
        "## Summary",
        "",
        "| B_scale | converged/biases | -20V status | current score | full qG score | full qG ratio | junction qG ratio |",
        "|---:|---:|---|---:|---:|---:|---:|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['B_scale']} | {row['converged_count']}/{row['bias_count']} | "
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
        f"1. Terminal current closest to Sentaurus: B_scale={best_current}.",
        f"2. Full qG score best: B_scale={best_full}.",
        f"3. Junction qG score best: B_scale={best_junction}.",
        f"4. B_scale<1 improves -20V low qG without significant -16V/-18V worsening: {'yes, ' + ', '.join(low_improves_20) if low_improves_20 else 'no'}.",
        f"5. B_scale<1 makes -16V/-18V stronger: {'yes, ' + ', '.join(sorted(set(low_makes_mid_stronger))) if low_makes_mid_stronger else 'no'}.",
        f"6. Single B_scale covering -5V to -20V by current and full qG 0.5..2 ratio gate: {'yes, ' + ', '.join(covering) if covering else 'no'}.",
    ])
    if not covering:
        lines.append("7. No single global B_scale cleanly covers -5V to -20V; next sweep should split B_low_scale and B_high_scale.")
    else:
        lines.append("7. A single global B_scale passes the coarse 0.5..2 gate, but low/high branch split remains the next diagnostic if shoulder errors persist.")
    if not solved_minus20:
        lines.append("8. All B_scale cases leave -20V outside the target qG range; inspect cutoff/smoothing/RefDens or source mapping.")
    else:
        lines.append("8. At least one B_scale brings -20V qG into the coarse target range; still check cutoff/smoothing/RefDens or source mapping if spatial windows disagree.")
    if recommended:
        lines.append(f"9. Recommended next full BV candidate from this sweep: B_scale={recommended[0]}.")
    else:
        lines.append("9. No B_scale=0.95 or 0.90 candidate simultaneously improves -20V without materially worsening -16V/-18V by this gate.")
    lines.extend([
        "10. sentaurus_fit_A_B remains explicitly forbidden as a full coupled default parameter set.",
        "",
        "## Files",
        "",
        f"- by-bias matrix: {path.parent / 'vanoverstraeten_A2_B_scale_full_bv_matrix.csv'}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def output_fields() -> list[str]:
    fields = [
        "case_id", "A_scale", "B_scale", "parameter_set", "bias_V", "bias_sign_convention",
        "run_return_code", "convergence_status", "newton_iterations", "iterations",
        "voltage_step_used_V", "actual_voltage_step_V",
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
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/simulation_coarse_previous_full20_aligned.json",
    )
    parser.add_argument("--runner", type=Path, default=default_runner())
    parser.add_argument("--out-dir", type=Path, default=REPO / "build/diagnostics")
    parser.add_argument(
        "--reference-current-csv",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/sentaurus_coarse_bv_reference_aligned.csv",
    )
    parser.add_argument(
        "--sentaurus-multibias-dir",
        type=Path,
        default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_previous_full20/sentaurus_multibias",
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
    case_root = args.out_dir / "vanoverstraeten_A2_B_scale_full_bv_matrix_cases"
    if args.skip_run:
        test_case_root = args.out_dir / "cases"
        if test_case_root.exists():
            case_root = test_case_root

    reference_current = bv.load_reference_current(args.reference_current_csv)
    sentaurus_generation = bv.load_sentaurus_generation_integrals(args.sentaurus_multibias_dir)

    all_rows: list[dict[str, Any]] = []
    for case_id, b_scale in CASES:
        case_dir = case_root / case_id
        config_path = build_case_config(args.base_config, case_dir, case_id, b_scale, biases)
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
            row["A_scale"] = fmt_a_scale(A_SCALE)
            row["B_scale"] = fmt_scale(b_scale)
            row["parameter_set"] = "default"
            add_ratios(row)
        all_rows.extend(case_rows)

    summary_rows = build_summary_rows(all_rows)
    csv_path = args.out_dir / "vanoverstraeten_A2_B_scale_full_bv_matrix.csv"
    md_path = args.out_dir / "vanoverstraeten_A2_B_scale_full_bv_matrix_summary.md"
    bv.write_rows(csv_path, all_rows, output_fields())
    write_markdown(md_path, summary_rows, all_rows)
    print(csv_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())