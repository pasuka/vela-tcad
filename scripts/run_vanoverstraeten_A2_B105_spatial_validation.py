#!/usr/bin/env python3
"""Run and summarize the PN2D A=2, B=1.05 spatial BV validation case."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

import run_vanoverstraeten_A2_B_scale_full_bv_matrix as bmatrix
import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
CASE_ID = "BV-A2-B1p05-spatial"
A_SCALE = 2.0
B_SCALE = 1.05
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]
WINDOWS = bv.WINDOWS


def default_runner() -> Path:
    exe = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
    build = REPO / "build" / exe
    if build.exists():
        return build
    return REPO / "build-release" / exe


def ratio(value: Any, reference: Any) -> float | None:
    lhs = bv.finite_float(value)
    rhs = bv.finite_float(reference)
    if lhs is None or rhs is None or rhs == 0.0:
        return None
    return lhs / rhs


def add_ratios(row: dict[str, Any]) -> None:
    row["contact_electron_current_A"] = row.get("terminal_electron_current_A")
    row["contact_hole_current_A"] = row.get("terminal_hole_current_A")
    row["contact_total_current_A"] = row.get("terminal_total_current_A")
    row["contact_electron_current_A_per_um"] = row.get("terminal_electron_current_A_per_um")
    row["contact_hole_current_A_per_um"] = row.get("terminal_hole_current_A_per_um")
    row["contact_total_current_A_per_um"] = row.get("terminal_total_current_A_per_um")
    row["terminal_total_current_ratio_to_sentaurus"] = ratio(
        row.get("terminal_total_current_A_per_um"),
        row.get("sentaurus_terminal_total_current_A"),
    )
    row["actual_voltage_step_V"] = row.get("voltage_step_used_V")
    for name in WINDOWS:
        row[f"qG_{name}_ratio_to_sentaurus"] = ratio(
            row.get(f"q_integral_AvalancheGeneration_{name}_A_per_um"),
            row.get(f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um"),
        )


def load_sentaurus_field_maxima(root: Path) -> dict[float, dict[str, Any]]:
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
        fields = case_dir / "fields"
        item: dict[str, Any] = {}
        add_sentaurus_max(item, nodes, fields / "ImpactIonization_region0.csv",
                          "sentaurus_max_Gava_cm3_per_s", 1.0)
        add_sentaurus_max(item, nodes, fields / "ElectricField_region0.csv",
                          "sentaurus_max_ElectricField_V_per_cm", 1.0)
        add_sentaurus_max(item, nodes, fields / "eAlphaAvalanche_region0.csv",
                          "sentaurus_max_eAlphaAvalanche_cm_inv", 1.0)
        add_sentaurus_max(item, nodes, fields / "hAlphaAvalanche_region0.csv",
                          "sentaurus_max_hAlphaAvalanche_cm_inv", 1.0)
        result[bv.bias_key(bias)] = item
    return result


def add_sentaurus_max(
    out: dict[str, Any],
    nodes: list[tuple[float, float]],
    path: Path,
    label: str,
    scale: float,
) -> None:
    values = bv.load_sentaurus_scalar(path)
    value, x_um, y_um = max_sentaurus_with_coord(nodes, values, scale)
    out[label] = value
    out[f"{label}_x_um"] = x_um
    out[f"{label}_y_um"] = y_um


def max_sentaurus_with_coord(
    nodes: list[tuple[float, float]],
    values: list[float],
    scale: float,
) -> tuple[float | None, float | None, float | None]:
    indexed = [(abs(value * scale), idx) for idx, value in enumerate(values) if idx < len(nodes)]
    if not indexed:
        return None, None, None
    value, idx = max(indexed, key=lambda item: item[0])
    x_m, y_m = nodes[idx]
    return value, x_m * 1.0e6, y_m * 1.0e6


def add_sentaurus_spatial_rows(rows: list[dict[str, Any]], sentaurus_maxima: dict[float, dict[str, Any]]) -> None:
    for row in rows:
        row.update(sentaurus_maxima.get(bv.bias_key(float(row["bias_V"])), {}))
        row["max_Gava_ratio_to_sentaurus"] = ratio(
            cm3_per_s(row.get("max_AvalancheGeneration_m3_per_s")),
            row.get("sentaurus_max_Gava_cm3_per_s"),
        )
        row["max_Gava_distance_to_sentaurus_um"] = distance_um(
            row.get("max_AvalancheGeneration_m3_per_s_x_um"),
            row.get("max_AvalancheGeneration_m3_per_s_y_um"),
            row.get("sentaurus_max_Gava_cm3_per_s_x_um"),
            row.get("sentaurus_max_Gava_cm3_per_s_y_um"),
        )


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


def load_b100_rows(path: Path) -> dict[float, dict[str, str]]:
    if not path.exists():
        return {}
    return {
        bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0): row
        for row in bv.read_rows(path)
        if row.get("B_scale") == "1.00"
    }


def add_b100_comparison(rows: list[dict[str, Any]], b100_rows: dict[float, dict[str, str]]) -> None:
    for row in rows:
        base = b100_rows.get(bv.bias_key(float(row["bias_V"])), {})
        for name in WINDOWS:
            row[f"A2_B100_qG_{name}_ratio_to_sentaurus"] = base.get(
                f"q_integral_AvalancheGeneration_{name}_ratio_to_sentaurus", "")
        row["A2_B100_terminal_total_current_ratio_to_sentaurus"] = base.get(
            "terminal_total_current_ratio_to_sentaurus", "")


def output_fields() -> list[str]:
    fields = [
        "case_id", "A_scale", "B_scale", "parameter_set", "bias_V", "bias_sign_convention",
        "run_return_code", "convergence_status", "newton_iterations", "iterations",
        "voltage_step_used_V", "actual_voltage_step_V",
        "terminal_electron_current_A", "terminal_hole_current_A", "terminal_total_current_A",
        "terminal_electron_current_A_per_um", "terminal_hole_current_A_per_um",
        "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A",
        "contact_electron_current_A", "contact_hole_current_A", "contact_total_current_A",
        "contact_electron_current_A_per_um", "contact_hole_current_A_per_um",
        "contact_total_current_A_per_um",
        "terminal_total_current_ratio_to_sentaurus",
        "A2_B100_terminal_total_current_ratio_to_sentaurus",
        "max_ElectricField_V_per_cm", "max_ElectricField_V_per_cm_x_um",
        "max_ElectricField_V_per_cm_y_um", "sentaurus_max_ElectricField_V_per_cm",
        "sentaurus_max_ElectricField_V_per_cm_x_um", "sentaurus_max_ElectricField_V_per_cm_y_um",
        "max_Fn_V_per_cm", "max_Fn_V_per_cm_x_um", "max_Fn_V_per_cm_y_um",
        "max_Fp_V_per_cm", "max_Fp_V_per_cm_x_um", "max_Fp_V_per_cm_y_um",
        "max_eAlphaAvalanche_cm_inv", "max_eAlphaAvalanche_cm_inv_x_um",
        "max_eAlphaAvalanche_cm_inv_y_um", "sentaurus_max_eAlphaAvalanche_cm_inv",
        "sentaurus_max_eAlphaAvalanche_cm_inv_x_um", "sentaurus_max_eAlphaAvalanche_cm_inv_y_um",
        "max_hAlphaAvalanche_cm_inv", "max_hAlphaAvalanche_cm_inv_x_um",
        "max_hAlphaAvalanche_cm_inv_y_um", "sentaurus_max_hAlphaAvalanche_cm_inv",
        "sentaurus_max_hAlphaAvalanche_cm_inv_x_um", "sentaurus_max_hAlphaAvalanche_cm_inv_y_um",
        "max_AvalancheGeneration_m3_per_s", "max_AvalancheGeneration_m3_per_s_x_um",
        "max_AvalancheGeneration_m3_per_s_y_um", "sentaurus_max_Gava_cm3_per_s",
        "sentaurus_max_Gava_cm3_per_s_x_um", "sentaurus_max_Gava_cm3_per_s_y_um",
        "max_Gava_ratio_to_sentaurus", "max_Gava_distance_to_sentaurus_um",
    ]
    for name in WINDOWS:
        fields.extend([
            f"q_integral_AvalancheGeneration_{name}_A_per_um",
            f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um",
            f"qG_{name}_ratio_to_sentaurus",
            f"A2_B100_qG_{name}_ratio_to_sentaurus",
        ])
    fields.append("vtk_file")
    return fields


def median_abs_ratio(rows: list[dict[str, Any]], column: str) -> float | None:
    values = [abs(value) for value in (bv.finite_float(row.get(column)) for row in rows)
              if value is not None and math.isfinite(value)]
    if not values:
        return None
    return bmatrix.median(values)


def ratio_gate(rows: list[dict[str, Any]], columns: list[str]) -> bool:
    for row in rows:
        bias = bv.finite_float(row.get("bias_V"))
        if bias == 0.0:
            continue
        for column in columns:
            value = bv.finite_float(row.get(column))
            if value is None or not (0.5 <= abs(value) <= 2.0):
                return False
    return True


def b105_better_than_b100(rows: list[dict[str, Any]], window: str) -> bool | None:
    better = False
    for row in rows:
        bias = bv.finite_float(row.get("bias_V"))
        if bias == 0.0:
            continue
        b105 = bv.finite_float(row.get(f"qG_{window}_ratio_to_sentaurus"))
        b100 = bv.finite_float(row.get(f"A2_B100_qG_{window}_ratio_to_sentaurus"))
        if b105 is None or b100 is None or b100 == 0.0:
            continue
        if abs(math.log10(abs(b105))) > abs(math.log10(abs(b100))):
            return False
        if abs(math.log10(abs(b105))) < abs(math.log10(abs(b100))):
            better = True
    return better if better else None


def max_gava_position_consistent(rows: list[dict[str, Any]], tolerance_um: float) -> bool | None:
    distances = [
        value for value in (bv.finite_float(row.get("max_Gava_distance_to_sentaurus_um")) for row in rows)
        if value is not None
    ]
    if not distances:
        return None
    return max(distances) <= tolerance_um


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    full_better = b105_better_than_b100(rows, "full_semiconductor")
    junction_better = b105_better_than_b100(rows, "junction_window")
    windows_ok = ratio_gate(rows, [
        "qG_center_window_ratio_to_sentaurus",
        "qG_left_shoulder_ratio_to_sentaurus",
        "qG_right_shoulder_ratio_to_sentaurus",
    ])
    gava_position_ok = max_gava_position_consistent(rows, tolerance_um=0.25)
    qg_close = ratio_gate(rows, [
        "qG_full_semiconductor_ratio_to_sentaurus",
        "qG_junction_window_ratio_to_sentaurus",
    ])
    current_close = ratio_gate(rows, ["terminal_total_current_ratio_to_sentaurus"])
    shoulder_ok = ratio_gate(rows, [
        "qG_left_shoulder_ratio_to_sentaurus",
        "qG_right_shoulder_ratio_to_sentaurus",
    ])

    lines = [
        "# VanOverstraeten A2 B1.05 Spatial Validation",
        "",
        "- bias sign convention: negative Anode bias is reverse bias for this PN2D deck.",
        "- impact_ionization.A_scale = 2.0",
        "- impact_ionization.B_scale = 1.05",
        "- parameter_set = default; sentaurus_fit_A_B is not used.",
        "",
        "## By-Bias Summary",
        "",
        "| bias V | status | current ratio | full qG | junction qG | center qG | left qG | right qG | max Gava distance um |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {bv.fmt(row.get('bias_V'))} | {row.get('convergence_status', '')} | "
            f"{bv.fmt(row.get('terminal_total_current_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('qG_full_semiconductor_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('qG_junction_window_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('qG_center_window_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('qG_left_shoulder_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('qG_right_shoulder_ratio_to_sentaurus'))} | "
            f"{bv.fmt(row.get('max_Gava_distance_to_sentaurus_um'))} |"
        )
    lines.extend([
        "",
        "## Aggregate Ratios",
        "",
        f"- median terminal current abs ratio: {bv.fmt(median_abs_ratio(rows, 'terminal_total_current_ratio_to_sentaurus'))}",
        f"- median full qG abs ratio: {bv.fmt(median_abs_ratio(rows, 'qG_full_semiconductor_ratio_to_sentaurus'))}",
        f"- median junction qG abs ratio: {bv.fmt(median_abs_ratio(rows, 'qG_junction_window_ratio_to_sentaurus'))}",
        f"- median center qG abs ratio: {bv.fmt(median_abs_ratio(rows, 'qG_center_window_ratio_to_sentaurus'))}",
        f"- median left shoulder qG abs ratio: {bv.fmt(median_abs_ratio(rows, 'qG_left_shoulder_ratio_to_sentaurus'))}",
        f"- median right shoulder qG abs ratio: {bv.fmt(median_abs_ratio(rows, 'qG_right_shoulder_ratio_to_sentaurus'))}",
        "",
        "## Answers",
        "",
        f"1. A2_B105 is stable better than A2_B100 on full qG: {yes_no_unknown(full_better)}; on junction qG: {yes_no_unknown(junction_better)}.",
        f"2. A2_B105 center/left/right windows all stay in 0.5..2 for nonzero biases: {'yes' if windows_ok else 'no'}.",
        f"3. Max Gava position matches Sentaurus within 0.25 um: {yes_no_unknown(gava_position_ok)}.",
    ])
    if qg_close and not windows_ok:
        lines.append("4. Total qG is close but spatial windows disagree; next step should inspect source mapping.")
    else:
        lines.append("4. Total qG/spatial-window agreement does not reduce to a pure total-vs-window conflict by this gate.")
    if not shoulder_ok:
        lines.append("5. Shoulder windows still disagree; next step should split B_low_scale and B_high_scale.")
    else:
        lines.append("5. Shoulder windows pass the coarse 0.5..2 gate; low/high branch split remains optional if finer spatial errors matter.")
    if qg_close and not current_close:
        lines.append("6. qG is close but terminal current is not; next step should inspect current density and continuity source feedback.")
    else:
        lines.append("6. Current and qG gates do not indicate a qG-close/current-far split; if that split appears in a later run, inspect current density and continuity source feedback.")
    lines.extend([
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'vanoverstraeten_A2_B105_spatial_validation.csv'}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def yes_no_unknown(value: bool | None) -> str:
    if value is None:
        return "insufficient data"
    return "yes" if value else "no"


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
    parser.add_argument(
        "--b-scale-matrix-csv",
        type=Path,
        default=REPO / "build/diagnostics/vanoverstraeten_A2_B_scale_full_bv_matrix.csv",
    )
    parser.add_argument("--bias-points", default=",".join(f"{bias:g}" for bias in DEFAULT_BIASES))
    parser.add_argument("--skip-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    biases = bv.parse_biases(args.bias_points)
    if 0.0 not in [bv.bias_key(bias) for bias in biases]:
        biases.insert(0, 0.0)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    case_root = args.out_dir / "vanoverstraeten_A2_B105_spatial_validation_case"
    case_dir = case_root / CASE_ID

    reference_current = bv.load_reference_current(args.reference_current_csv)
    sentaurus_generation = bv.load_sentaurus_generation_integrals(args.sentaurus_multibias_dir)
    sentaurus_maxima = load_sentaurus_field_maxima(args.sentaurus_multibias_dir)

    config_path = bmatrix.build_case_config(args.base_config, case_dir, CASE_ID, B_SCALE, biases)
    run_return_code = 0
    if not args.skip_run:
        run_return_code = bv.run_case(args.runner, config_path, case_dir)

    rows = bv.collect_case_rows(
        case_root,
        CASE_ID,
        "default",
        biases,
        reference_current,
        sentaurus_generation,
        run_return_code,
    )
    for row in rows:
        row["A_scale"] = "2.0"
        row["B_scale"] = "1.05"
        row["parameter_set"] = "default"
        add_ratios(row)
    add_sentaurus_spatial_rows(rows, sentaurus_maxima)
    add_b100_comparison(rows, load_b100_rows(args.b_scale_matrix_csv))

    csv_path = args.out_dir / "vanoverstraeten_A2_B105_spatial_validation.csv"
    md_path = args.out_dir / "vanoverstraeten_A2_B105_spatial_validation_summary.md"
    bv.write_rows(csv_path, rows, output_fields())
    write_summary(md_path, rows)
    print(csv_path)
    print(md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
