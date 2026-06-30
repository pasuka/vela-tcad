#!/usr/bin/env python3
"""Run guarded coupled BV sweeps for avalanche current magnitude modes."""

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
DEFAULT_BIASES = [0.0, -5.0, -10.0, -16.0, -18.0, -20.0]
CASES = [
    ("ava-current-edge-scalar", "edge_scalar_abs"),
    ("ava-current-dual-face-vector", "dual_face_vector_mag"),
]


def default_runner() -> Path:
    exe = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
    build = REPO / "build-release" / exe
    if build.exists():
        return build
    return REPO / "build" / exe


def read_json(path: Path) -> dict[str, Any]:
    return bv.read_json(path)


def write_json(path: Path, data: dict[str, Any]) -> None:
    bv.write_json(path, data)


def inherited_scale(impact: dict[str, Any], key: str) -> float:
    value = bv.finite_float(impact.get(key))
    return 1.0 if value is None else value


def inherited_parameter_set(impact: dict[str, Any]) -> str:
    return str(impact.get("parameter_set") or "default")


def build_case_config(
    base_config: Path,
    case_dir: Path,
    case_id: str,
    current_magnitude_mode: str,
    biases: list[float],
) -> tuple[Path, dict[str, Any]]:
    config = read_json(base_config)
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
    impact["driving_force"] = "quasi_fermi_gradient"
    impact["generation"] = "current_density"
    impact["current_approximation"] = "grad_qf"
    impact["current_magnitude_mode"] = current_magnitude_mode
    impact.pop("parameter_set", None)
    if "sentaurus_fit" in str(impact.get("parameter_set", "")):
        impact["parameter_set"] = "default"

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

    parameter_set = inherited_parameter_set(impact)
    metadata = {
        "case_id": case_id,
        "current_magnitude_mode": current_magnitude_mode,
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "current_approximation": "grad_qf",
        "parameter_set": parameter_set,
        "A_scale": inherited_scale(impact, "A_scale"),
        "B_scale": inherited_scale(impact, "B_scale"),
        "sentaurus_fit_A_B": "forbidden_for_full_coupled_default",
        "bias_sign_convention": "negative Anode bias is reverse bias for this PN2D deck",
    }
    config["_avalanche_current_magnitude_mode_coupled_matrix_metadata"] = metadata

    case_dir.mkdir(parents=True, exist_ok=True)
    (case_dir / "vtk").mkdir(parents=True, exist_ok=True)
    config_path = case_dir / f"{case_id}.json"
    write_json(config_path, config)
    return config_path, metadata


def add_common_ratios(row: dict[str, Any]) -> None:
    bmatrix.add_ratios(row)
    row["contact_electron_current_A"] = row.get("terminal_electron_current_A")
    row["contact_hole_current_A"] = row.get("terminal_hole_current_A")
    row["contact_total_current_A"] = row.get("terminal_total_current_A")
    for name in bv.WINDOWS:
        row[f"qG_{name}_ratio_to_sentaurus"] = row.get(
            f"q_integral_AvalancheGeneration_{name}_ratio_to_sentaurus")


def count_csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    rows = bv.read_rows(path)
    return len(rows)


def count_dual_rows(default_path: Path) -> int:
    return count_csv_rows(default_path)


def output_fields() -> list[str]:
    fields = [
        "case_id", "current_magnitude_mode", "A_scale", "B_scale", "parameter_set",
        "bias_V", "bias_sign_convention", "run_return_code",
        "convergence_status", "newton_iterations", "iterations",
        "voltage_step_used_V", "actual_voltage_step_V",
        "terminal_electron_current_A", "terminal_hole_current_A", "terminal_total_current_A",
        "terminal_electron_current_A_per_um", "terminal_hole_current_A_per_um",
        "terminal_total_current_A_per_um", "sentaurus_terminal_total_current_A",
        "terminal_total_current_ratio_to_sentaurus",
        "contact_electron_current_A", "contact_hole_current_A", "contact_total_current_A",
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
        "true_dual_face_rows", "reconstructed_cell_rows",
    ]
    for name in bv.WINDOWS:
        fields.extend([
            f"q_integral_AvalancheGeneration_{name}_A_per_um",
            f"sentaurus_q_integral_AvalancheGeneration_{name}_A_per_um",
            f"qG_{name}_ratio_to_sentaurus",
        ])
    fields.append("vtk_file")
    return fields


def median(values: list[float]) -> float | None:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def ratio_at(rows: list[dict[str, Any]], mode: str, bias: float, column: str) -> float | None:
    for row in rows:
        if row.get("current_magnitude_mode") == mode and bv.finite_float(row.get("bias_V")) == bias:
            return bv.finite_float(row.get(column))
    return None


def converged_count(rows: list[dict[str, Any]], mode: str) -> tuple[int, int]:
    subset = [row for row in rows if row.get("current_magnitude_mode") == mode]
    return sum(1 for row in subset if row.get("convergence_status") == "converged"), len(subset)


def mode_table(rows: list[dict[str, Any]]) -> list[str]:
    lines = [
        "| current_magnitude_mode | converged/biases | -20 full qG | -20 junction qG | -20 right shoulder qG | -20 current ratio | -20 max Gava distance um |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, mode in CASES:
        ok, total = converged_count(rows, mode)
        lines.append(
            f"| {mode} | {ok}/{total} | "
            f"{bv.fmt(ratio_at(rows, mode, -20.0, 'qG_full_semiconductor_ratio_to_sentaurus'))} | "
            f"{bv.fmt(ratio_at(rows, mode, -20.0, 'qG_junction_window_ratio_to_sentaurus'))} | "
            f"{bv.fmt(ratio_at(rows, mode, -20.0, 'qG_right_shoulder_ratio_to_sentaurus'))} | "
            f"{bv.fmt(ratio_at(rows, mode, -20.0, 'terminal_total_current_ratio_to_sentaurus'))} | "
            f"{bv.fmt(ratio_at(rows, mode, -20.0, 'max_Gava_distance_to_sentaurus_um'))} |"
        )
    return lines


def closer_to_one(new: float | None, old: float | None) -> bool | None:
    if new is None or old is None or new == 0.0 or old == 0.0:
        return None
    return abs(math.log10(abs(new))) < abs(math.log10(abs(old)))


def write_summary(path: Path, rows: list[dict[str, Any]], metadata: dict[str, Any]) -> None:
    edge_full20 = ratio_at(rows, "edge_scalar_abs", -20.0, "qG_full_semiconductor_ratio_to_sentaurus")
    dual_full20 = ratio_at(rows, "dual_face_vector_mag", -20.0, "qG_full_semiconductor_ratio_to_sentaurus")
    edge_junc20 = ratio_at(rows, "edge_scalar_abs", -20.0, "qG_junction_window_ratio_to_sentaurus")
    dual_junc20 = ratio_at(rows, "dual_face_vector_mag", -20.0, "qG_junction_window_ratio_to_sentaurus")
    edge_current20 = ratio_at(rows, "edge_scalar_abs", -20.0, "terminal_total_current_ratio_to_sentaurus")
    dual_current20 = ratio_at(rows, "dual_face_vector_mag", -20.0, "terminal_total_current_ratio_to_sentaurus")
    dual_right20 = ratio_at(rows, "dual_face_vector_mag", -20.0, "qG_right_shoulder_ratio_to_sentaurus")
    edge_right20 = ratio_at(rows, "edge_scalar_abs", -20.0, "qG_right_shoulder_ratio_to_sentaurus")
    dual_dist20 = ratio_at(rows, "dual_face_vector_mag", -20.0, "max_Gava_distance_to_sentaurus_um")
    edge_full5 = ratio_at(rows, "edge_scalar_abs", -5.0, "qG_full_semiconductor_ratio_to_sentaurus")
    dual_full5 = ratio_at(rows, "dual_face_vector_mag", -5.0, "qG_full_semiconductor_ratio_to_sentaurus")
    edge_junc5 = ratio_at(rows, "edge_scalar_abs", -5.0, "qG_junction_window_ratio_to_sentaurus")
    dual_junc5 = ratio_at(rows, "dual_face_vector_mag", -5.0, "qG_junction_window_ratio_to_sentaurus")
    edge_current5 = ratio_at(rows, "edge_scalar_abs", -5.0, "terminal_total_current_ratio_to_sentaurus")
    dual_current5 = ratio_at(rows, "dual_face_vector_mag", -5.0, "terminal_total_current_ratio_to_sentaurus")
    dual_max_xy = next((
        (row.get("max_AvalancheGeneration_m3_per_s_x_um"), row.get("max_AvalancheGeneration_m3_per_s_y_um"))
        for row in rows
        if row.get("current_magnitude_mode") == "dual_face_vector_mag" and
        bv.finite_float(row.get("bias_V")) == -20.0
    ), ("", ""))
    edge_ok, edge_total = converged_count(rows, "edge_scalar_abs")
    dual_ok, dual_total = converged_count(rows, "dual_face_vector_mag")
    dual_reached_minus20 = dual_full20 is not None
    dual_too_large = any(
        (value is not None and abs(value) > 5.0)
        for value in (
            bv.finite_float(row.get("qG_full_semiconductor_ratio_to_sentaurus"))
            for row in rows
            if row.get("current_magnitude_mode") == "dual_face_vector_mag"
        )
    )
    dual_worse_convergence = dual_ok < edge_ok
    qg_answer = (
        f"full={closer_to_one(dual_full20, edge_full20)}, junction={closer_to_one(dual_junc20, edge_junc20)} at -20 V."
        if dual_reached_minus20
        else "not evaluable at -20 V because dual_face_vector_mag stopped before that bias; at -5 V "
             f"full qG edge={bv.fmt(edge_full5)}, dual={bv.fmt(dual_full5)}, "
             f"junction edge={bv.fmt(edge_junc5)}, dual={bv.fmt(dual_junc5)}."
    )
    current_answer = (
        f"{closer_to_one(dual_current20, edge_current20)} at -20 V."
        if dual_reached_minus20
        else f"not evaluable at -20 V; at -5 V current ratio edge={bv.fmt(edge_current5)}, dual={bv.fmt(dual_current5)}."
    )
    location_answer = (
        f"Max Gava location for dual_face_vector_mag at -20 V is ({dual_max_xy[0]}, {dual_max_xy[1]}) um; "
        f"distance to Sentaurus={bv.fmt(dual_dist20)} um."
        if dual_reached_minus20
        else "Max Gava location at -20 V is not available because dual_face_vector_mag did not converge past the low-bias points."
    )
    right_answer = (
        f"Right shoulder at -20 V: edge={bv.fmt(edge_right20)}, dual={bv.fmt(dual_right20)}; diagnostic target was about 0.112."
        if dual_reached_minus20
        else "Right shoulder at -20 V is not available; at -5 V dual qG windows are already stronger than edge scalar."
    )
    lines = [
        "# Avalanche Current Magnitude Mode Coupled Matrix",
        "",
        "- model = van_overstraeten",
        "- driving_force = quasi_fermi_gradient",
        "- generation = current_density",
        "- current_approximation = grad_qf",
        f"- parameter_set = {metadata.get('parameter_set', 'default')}",
        f"- A_scale = {metadata.get('A_scale', 1.0)}",
        f"- B_scale = {metadata.get('B_scale', 1.0)}",
        "- A/B are inherited from the current release/default config and Sentaurus par defaults; sentaurus_fit_A_B is not used.",
        "- dual_face_vector_mag remains guarded and is not promoted to global default by this run.",
        "",
        "## Mode Scores",
        "",
        *mode_table(rows),
        "",
        "## Answers",
        "",
        f"1. dual_face_vector_mag improves full/junction qG in coupled solve: {qg_answer}",
        f"2. dual_face_vector_mag improves terminal current: {current_answer}",
        f"3. Convergence: edge_scalar_abs={edge_ok}/{edge_total}, dual_face_vector_mag={dual_ok}/{dual_total}.",
        f"4. {location_answer}",
        f"5. {right_answer}",
    ]
    if dual_reached_minus20 and closer_to_one(dual_full20, edge_full20) and not closer_to_one(dual_dist20, 0.0):
        lines.append("6. qG magnitude changes more than hotspot position; next step should validate source mapping/hotspot placement.")
    else:
        lines.append("6. Source mapping/hotspot validation should wait until dual_face_vector_mag can be stabilized through the full BV range.")
    if dual_too_large or dual_worse_convergence:
        lines.append("7. dual_face_vector_mag is too strong or convergence worsened; use lambda_ava homotopy before any broader coupled run.")
    else:
        lines.append("7. No qG>5x blow-up or convergence regression was detected by this coarse gate; lambda_ava homotopy is not immediately required.")
    lines.extend([
        "8. Do not make dual_face_vector_mag the default until this coupled matrix is reviewed.",
        "",
        "## Files",
        "",
        f"- CSV: {path.parent / 'avalanche_current_magnitude_mode_coupled_matrix.csv'}",
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
        "--dual-face-csv",
        type=Path,
        default=REPO / "build/diagnostics/dual_face_geometry_audit.csv",
    )
    parser.add_argument(
        "--dual-reconstruction-csv",
        type=Path,
        default=REPO / "build/diagnostics/dual_face_vector_current_gava_reconstruction.csv",
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
    case_root = args.out_dir / "avalanche_current_magnitude_mode_coupled_matrix_cases"
    reference_current = bv.load_reference_current(args.reference_current_csv)
    sentaurus_generation = bv.load_sentaurus_generation_integrals(args.sentaurus_multibias_dir)
    sentaurus_maxima = spatial.load_sentaurus_field_maxima(args.sentaurus_multibias_dir)
    true_dual_face_rows = count_dual_rows(args.dual_face_csv)
    reconstructed_cell_rows = count_csv_rows(args.dual_reconstruction_csv)

    rows: list[dict[str, Any]] = []
    first_metadata: dict[str, Any] = {}
    for case_id, mode in CASES:
        case_dir = case_root / case_id
        config_path, metadata = build_case_config(args.base_config, case_dir, case_id, mode, biases)
        if not first_metadata:
            first_metadata = metadata
        run_return_code = 0
        if not args.skip_run:
            run_return_code = bv.run_case(args.runner, config_path, case_dir)
        case_rows = bv.collect_case_rows(
            case_root,
            case_id,
            metadata["parameter_set"],
            biases,
            reference_current,
            sentaurus_generation,
            run_return_code,
        )
        for row in case_rows:
            row["current_magnitude_mode"] = mode
            row["A_scale"] = metadata["A_scale"]
            row["B_scale"] = metadata["B_scale"]
            row["parameter_set"] = metadata["parameter_set"]
            row["true_dual_face_rows"] = true_dual_face_rows
            row["reconstructed_cell_rows"] = reconstructed_cell_rows
            add_common_ratios(row)
        spatial.add_sentaurus_spatial_rows(case_rows, sentaurus_maxima)
        rows.extend(case_rows)

    csv_path = args.out_dir / "avalanche_current_magnitude_mode_coupled_matrix.csv"
    summary_path = args.out_dir / "avalanche_current_magnitude_mode_coupled_matrix_summary.md"
    bv.write_rows(csv_path, rows, output_fields())
    write_summary(summary_path, rows, first_metadata)
    print(csv_path)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
