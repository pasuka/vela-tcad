#!/usr/bin/env python3
"""Run coupled BV sweeps for A2/B1.05 avalanche source mapping modes."""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from typing import Any

import run_vanoverstraeten_A2_B105_spatial_validation as spatial
import run_vanoverstraeten_A2_B_scale_full_bv_matrix as bmatrix
import run_vanoverstraeten_bv_parameter_matrix as bv


REPO = Path(__file__).resolve().parents[1]
A_SCALE = 2.0
B_SCALE = 1.05
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]
CASES = [
    ("A2B105-map-node", "node_F_node_alpha_node_G"),
    ("A2B105-map-edge", "edge_F_edge_alpha_edge_G_to_node"),
    ("A2B105-map-cell", "cell_F_cell_alpha_cell_G_to_node"),
]


def default_runner() -> Path:
    exe = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
    build = REPO / "build" / exe
    if build.exists():
        return build
    return REPO / "build-release" / exe


def read_json(path: Path) -> dict[str, Any]:
    return bv.read_json(path)


def write_json(path: Path, data: dict[str, Any]) -> None:
    bv.write_json(path, data)


def build_case_config(
    base_config: Path,
    case_dir: Path,
    case_id: str,
    source_mapping_mode: str,
    biases: list[float],
) -> Path:
    config_path = bmatrix.build_case_config(base_config, case_dir, case_id, B_SCALE, biases)
    config = read_json(config_path)
    impact = config.setdefault("solver", {}).setdefault("impact_ionization", {})
    impact["model"] = "van_overstraeten"
    impact["parameter_set"] = "default"
    impact["driving_force"] = "quasi_fermi_gradient"
    impact["generation"] = "current_density"
    impact["current_approximation"] = "grad_qf"
    impact["A_scale"] = A_SCALE
    impact["B_scale"] = B_SCALE
    impact["source_mapping_mode"] = source_mapping_mode
    config["_A2_B105_source_mapping_coupled_matrix_metadata"] = {
        "case_id": case_id,
        "source_mapping_mode": source_mapping_mode,
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "A_scale": A_SCALE,
        "B_scale": B_SCALE,
        "parameter_set": "default",
        "sentaurus_fit_A_B": "forbidden_for_full_coupled_default",
    }
    write_json(config_path, config)
    return config_path


def add_common_ratios(row: dict[str, Any]) -> None:
    bmatrix.add_ratios(row)
    for name in bv.WINDOWS:
        row[f"qG_{name}_ratio_to_sentaurus"] = row.get(
            f"q_integral_AvalancheGeneration_{name}_ratio_to_sentaurus")


def add_reference_reproduction(rows: list[dict[str, Any]], reference_csv: Path) -> None:
    refs = {
        bv.bias_key(bv.finite_float(row.get("bias_V")) or 0.0): row
        for row in bv.read_rows(reference_csv)
    }
    for row in rows:
        if row.get("source_mapping_mode") != "node_F_node_alpha_node_G":
            continue
        ref = refs.get(bv.bias_key(float(row["bias_V"])), {})
        row["node_reproduce_A2_B105_current_ratio_delta"] = delta(
            row.get("terminal_total_current_ratio_to_sentaurus"),
            ref.get("terminal_total_current_ratio_to_sentaurus"),
        )
        for name in bv.WINDOWS:
            row[f"node_reproduce_A2_B105_{name}_qG_ratio_delta"] = delta(
                row.get(f"qG_{name}_ratio_to_sentaurus"),
                ref.get(f"qG_{name}_ratio_to_sentaurus"),
            )


def delta(value: Any, reference: Any) -> float | None:
    lhs = bv.finite_float(value)
    rhs = bv.finite_float(reference)
    if lhs is None or rhs is None:
        return None
    return lhs - rhs


def output_fields() -> list[str]:
    fields = [
        "case_id", "source_mapping_mode", "A_scale", "B_scale", "parameter_set",
        "bias_V", "bias_sign_convention", "run_return_code",
        "convergence_status", "newton_iterations", "iterations",
        "voltage_step_used_V", "actual_voltage_step_V",
        "terminal_electron_current_A", "terminal_hole_current_A", "terminal_total_current_A",
        "terminal_electron_current_A_per_um", "terminal_hole_current_A_per_um",
        "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A",
        "terminal_total_current_ratio_to_sentaurus",
    ]
    for name in bv.WINDOWS:
        fields.extend([
            f"q_integral_AvalancheGeneration_{name}_A_per_um",
            f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um",
            f"qG_{name}_ratio_to_sentaurus",
            f"node_reproduce_A2_B105_{name}_qG_ratio_delta",
        ])
    fields.extend([
        "node_reproduce_A2_B105_current_ratio_delta",
        "max_AvalancheGeneration_m3_per_s", "max_AvalancheGeneration_m3_per_s_x_um",
        "max_AvalancheGeneration_m3_per_s_y_um", "sentaurus_max_Gava_cm3_per_s",
        "sentaurus_max_Gava_cm3_per_s_x_um", "sentaurus_max_Gava_cm3_per_s_y_um",
        "max_Gava_distance_to_sentaurus_um",
        "max_eAlphaAvalanche_cm_inv", "max_eAlphaAvalanche_cm_inv_x_um",
        "max_eAlphaAvalanche_cm_inv_y_um",
        "max_hAlphaAvalanche_cm_inv", "max_hAlphaAvalanche_cm_inv_x_um",
        "max_hAlphaAvalanche_cm_inv_y_um",
        "max_Fn_V_per_cm", "max_Fn_V_per_cm_x_um", "max_Fn_V_per_cm_y_um",
        "max_Fp_V_per_cm", "max_Fp_V_per_cm_x_um", "max_Fp_V_per_cm_y_um",
        "max_ElectricField_V_per_cm", "max_ElectricField_V_per_cm_x_um",
        "max_ElectricField_V_per_cm_y_um",
        "vtk_file",
    ])
    return fields


def median(values: list[float]) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def score_mode(rows: list[dict[str, Any]], mode: str, columns: list[str]) -> float | None:
    errors: list[float] = []
    for row in rows:
        if row.get("source_mapping_mode") != mode or bv.finite_float(row.get("bias_V")) == 0.0:
            continue
        for col in columns:
            value = bv.finite_float(row.get(col))
            if value is not None and value != 0.0:
                errors.append(abs(math.log10(abs(value))))
    return median(errors)


def best_mode(rows: list[dict[str, Any]], columns: list[str]) -> str:
    best = "insufficient_data"
    best_score = math.inf
    for _, mode in CASES:
        score = score_mode(rows, mode, columns)
        if score is not None and score < best_score:
            best = mode
            best_score = score
    return best


def best_position_mode(rows: list[dict[str, Any]]) -> str:
    best = "insufficient_data"
    best_score = math.inf
    for _, mode in CASES:
        values = [
            value for value in (
                bv.finite_float(row.get("max_Gava_distance_to_sentaurus_um"))
                for row in rows
                if row.get("source_mapping_mode") == mode and bv.finite_float(row.get("bias_V")) != 0.0
            )
            if value is not None
        ]
        score = median(values)
        if score is not None and score < best_score:
            best = mode
            best_score = score
    return best


def mode_score_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = [
        "| source_mapping_mode | converged/biases | position score | shoulder score | full/junction/center score |",
        "|---|---:|---:|---:|---:|",
    ]
    for _, mode in CASES:
        subset = [row for row in rows if row.get("source_mapping_mode") == mode]
        nonzero = [row for row in subset if bv.finite_float(row.get("bias_V")) != 0.0]
        distances = [
            value for value in (bv.finite_float(row.get("max_Gava_distance_to_sentaurus_um")) for row in nonzero)
            if value is not None
        ]
        lines.append(
            f"| {mode} | "
            f"{sum(1 for row in subset if row.get('convergence_status') == 'converged')}/{len(subset)} | "
            f"{bv.fmt(median(distances))} | "
            f"{bv.fmt(score_mode(rows, mode, ['qG_left_shoulder_ratio_to_sentaurus', 'qG_right_shoulder_ratio_to_sentaurus']))} | "
            f"{bv.fmt(score_mode(rows, mode, ['qG_full_semiconductor_ratio_to_sentaurus', 'qG_junction_window_ratio_to_sentaurus', 'qG_center_window_ratio_to_sentaurus']))} |"
        )
    return lines


def node_reproduction_max_abs_delta(rows: list[dict[str, Any]]) -> float | None:
    values: list[float] = []
    for row in rows:
        if row.get("source_mapping_mode") != "node_F_node_alpha_node_G":
            continue
        for key, value in row.items():
            if key.startswith("node_reproduce_A2_B105_") and key.endswith("_delta"):
                parsed = bv.finite_float(value)
                if parsed is not None:
                    values.append(abs(parsed))
    return max(values) if values else None


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    best_pos = best_position_mode(rows)
    best_shoulder = best_mode(rows, [
        "qG_left_shoulder_ratio_to_sentaurus",
        "qG_right_shoulder_ratio_to_sentaurus",
    ])
    best_core = best_mode(rows, [
        "qG_full_semiconductor_ratio_to_sentaurus",
        "qG_junction_window_ratio_to_sentaurus",
        "qG_center_window_ratio_to_sentaurus",
    ])
    node_delta = node_reproduction_max_abs_delta(rows)
    edge_or_cell_improves = best_shoulder != "node_F_node_alpha_node_G"
    lines = [
        "# VanOverstraeten A2 B1.05 Coupled Source Mapping Matrix",
        "",
        "- model = van_overstraeten",
        "- driving_force = quasi_fermi_gradient",
        "- A_scale = 2.0",
        "- B_scale = 1.05",
        "- parameter_set = default",
        "- sentaurus_fit_A_B is forbidden as a full coupled default.",
        "- Each row records source_mapping_mode metadata; default parameter values are not changed.",
        "",
        "## Mode Scores",
        "",
        *mode_score_table(rows),
        "",
        "## Answers",
        "",
        f"1. Max Gava position closest to Sentaurus: {best_pos}.",
        f"2. Left/right shoulder qG closest to Sentaurus: {best_shoulder}.",
        f"3. Full/junction/center qG closest to Sentaurus: {best_core}.",
        f"4. Edge/cell mapping improves spatial windows over node mapping: {'yes' if edge_or_cell_improves else 'no'}.",
    ]
    if best_shoulder == "cell_F_cell_alpha_cell_G_to_node" or best_core == "cell_F_cell_alpha_cell_G_to_node":
        lines.append("5. cell_F_cell_alpha_cell_G_to_node is a strong optional coupled source mapping candidate, but should not be promoted to global default yet.")
    else:
        lines.append("5. cell_F_cell_alpha_cell_G_to_node is not clearly dominant in this coupled matrix.")
    if not edge_or_cell_improves:
        lines.append("6. Coupled mapping did not improve; next run should do B_low/B_high split plus cutoff/RefDens checks.")
    else:
        lines.append("6. Since coupled mapping changes the spatial windows, defer B_low/B_high and cutoff/RefDens until this source mapping branch is reviewed.")
    lines.extend([
        f"7. node_F_node_alpha_node_G reproduction max abs ratio delta vs current A2_B105: {bv.fmt(node_delta)}.",
        "8. No default parameter set was changed; all mapping modes are diagnostic metadata in the generated configs.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'vanoverstraeten_A2_B105_source_mapping_coupled_matrix.csv'}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


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
        "--a2-b105-reference-csv",
        type=Path,
        default=REPO / "build/diagnostics/vanoverstraeten_A2_B105_spatial_validation.csv",
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
    case_root = args.out_dir / "vanoverstraeten_A2_B105_source_mapping_coupled_matrix_cases"
    reference_current = bv.load_reference_current(args.reference_current_csv)
    sentaurus_generation = bv.load_sentaurus_generation_integrals(args.sentaurus_multibias_dir)
    sentaurus_maxima = spatial.load_sentaurus_field_maxima(args.sentaurus_multibias_dir)

    rows: list[dict[str, Any]] = []
    for case_id, mode in CASES:
        case_dir = case_root / case_id
        config_path = build_case_config(args.base_config, case_dir, case_id, mode, biases)
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
            row["A_scale"] = "2.0"
            row["B_scale"] = "1.05"
            row["source_mapping_mode"] = mode
            row["parameter_set"] = "default"
            add_common_ratios(row)
        spatial.add_sentaurus_spatial_rows(case_rows, sentaurus_maxima)
        rows.extend(case_rows)
    add_reference_reproduction(rows, args.a2_b105_reference_csv)

    csv_path = args.out_dir / "vanoverstraeten_A2_B105_source_mapping_coupled_matrix.csv"
    summary_path = args.out_dir / "vanoverstraeten_A2_B105_source_mapping_coupled_matrix_summary.md"
    bv.write_rows(csv_path, rows, output_fields())
    write_summary(summary_path, rows)
    print(csv_path)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
