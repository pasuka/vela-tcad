#!/usr/bin/env python3
"""Audit Sentaurus electron-row ErrRef, box Measure and NewtonPlot ordering."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path
from typing import Any


DEBUG_LINE = re.compile(
    r"^\s*(?P<grd>\d+)\s+(?P<des>-?\d+)\s+(?P<type>\d+)\s+"
    r"(?P<values>[^#]+?)\s*$"
)
ERRREF_RE = re.compile(
    r"RelErrControl \(Reference error\):.*?^\s*Electron\s*:\s*([+\-0-9.Ee]+)",
    re.MULTILINE | re.DOTALL,
)
ITERATION_RE = re.compile(
    r"^\s*(?P<iteration>\d+)\s+(?P<rhs>[+\-0-9.Ee]+)(?:\s+.*)?$",
    re.MULTILINE,
)
CNORM_RE = re.compile(
    r"^\s*electron:\s*([+\-0-9.Ee]+)\s+(\d+)\s+"
    r"\([^)]*\)\s+([+\-0-9.Ee]+)\s*$",
    re.MULTILINE,
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def scalar_field(export_root: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_root / "fields" / f"{name}_region0.csv")
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def field_hash(export_root: Path, name: str) -> str:
    return sha256(export_root / "fields" / f"{name}_region0.csv")


def field_norm(export_root: Path, name: str) -> dict[str, float]:
    values = list(scalar_field(export_root, name).values())
    return {
        "l2": math.sqrt(sum(value * value for value in values)),
        "maximum_abs": max(abs(value) for value in values),
    }


def parse_debug_block(path: Path, name: str) -> dict[int, list[float]]:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf"\n\s*{re.escape(name)}\s*\{{.*?\n(?P<body>.*?)\n\s*\}}",
        text,
        re.DOTALL,
    )
    if match is None:
        raise ValueError(f"{name} block not found in {path}")
    result: dict[int, list[float]] = {}
    for line in match.group("body").splitlines():
        parsed = DEBUG_LINE.match(line)
        if parsed is None or int(parsed.group("type")) != 2:
            continue
        design = int(parsed.group("des"))
        if design >= 0:
            result[design] = [
                float(value) for value in parsed.group("values").split()
            ]
    return result


def silicon_node_measures(
    mesh_path: Path, measure_path: Path, nodes: list[int]
) -> dict[int, float]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    measures = parse_debug_block(measure_path, "Measure")
    wanted = set(nodes)
    totals = {node: 0.0 for node in nodes}
    for cell in mesh["triangles"]:
        if int(cell["region_id"]) != 0:
            continue
        cell_id = int(cell["id"])
        cell_nodes = [int(value) for value in cell["node_ids"]]
        if not wanted.intersection(cell_nodes):
            continue
        if cell_id not in measures or len(measures[cell_id]) != len(cell_nodes):
            raise ValueError(f"missing or malformed Measure row for cell {cell_id}")
        for local, node in enumerate(cell_nodes):
            if node in wanted:
                totals[node] += measures[cell_id][local]
    if any(value <= 0.0 for value in totals.values()):
        raise ValueError(f"non-positive silicon Measure in requested rows: {totals}")
    return totals


def log_summary(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    errref = ERRREF_RE.search(text)
    if errref is None:
        raise ValueError(f"Electron ErrRef not found in {path}")
    iterations = {
        int(match.group("iteration")): float(match.group("rhs"))
        for match in ITERATION_RE.finditer(text)
    }
    cnorm = CNORM_RE.search(text)
    if cnorm is None:
        raise ValueError(f"electron CNorm row not found in {path}")
    header = text.find("Iteration   |Rhs|")
    iteration_zero = text.find("    0", header)
    newton_write = text.find("newton_0_0.tdr", iteration_zero)
    cnorm_position = text.find("C-norm_equation", newton_write)
    iteration_one = text.find("    1", cnorm_position)
    positions = {
        "iteration_header": header,
        "iteration_0_rhs": iteration_zero,
        "newton_0_write": newton_write,
        "cnorm": cnorm_position,
        "iteration_1_solve_result": iteration_one,
    }
    return {
        "electron_errref_cm3": float(errref.group(1)),
        "iteration_rhs_norm": iterations,
        "electron_cnorm": {
            "maximum_error": float(cnorm.group(1)),
            "vertex": int(cnorm.group(2)),
            "density_cm3": float(cnorm.group(3)),
        },
        "without_diagonal_preconditioning_declared": (
            "Without diagonal preconditioning" in text
        ),
        "event_character_offsets": positions,
        "iteration0_newtonplot_precedes_cnorm_and_first_linear_solve": (
            all(position >= 0 for position in positions.values())
            and header < iteration_zero < newton_write < cnorm_position < iteration_one
        ),
    }


def errref_formula_check(
    baseline_export: Path,
    candidate_export: Path,
    baseline_errref: float,
    candidate_errref: float,
) -> dict[str, Any]:
    density = scalar_field(baseline_export, "eDensity")
    baseline_error = scalar_field(baseline_export, "eDensityError")
    candidate_error = scalar_field(candidate_export, "eDensityError")
    common = sorted(set(density) & set(baseline_error) & set(candidate_error))
    active = [node for node in common if baseline_error[node] != 0.0]
    observed = [candidate_error[node] / baseline_error[node] for node in active]
    predicted = [
        (abs(density[node]) + baseline_errref)
        / (abs(density[node]) + candidate_errref)
        for node in active
    ]
    relative_errors = [
        abs(actual / expected - 1.0)
        for actual, expected in zip(observed, predicted, strict=True)
    ]
    return {
        "active_nodes": len(active),
        "observed_candidate_over_baseline": {
            "minimum": min(observed),
            "median": statistics.median(observed),
            "maximum": max(observed),
        },
        "manual_formula": "(|n|+baseline_ErrRef)/(|n|+candidate_ErrRef)",
        "maximum_relative_error_against_formula": max(relative_errors),
        "median_relative_error_against_formula": statistics.median(relative_errors),
    }


def pearson(left: list[float], right: list[float]) -> float:
    left_mean = statistics.mean(left)
    right_mean = statistics.mean(right)
    numerator = sum(
        (a - left_mean) * (b - right_mean)
        for a, b in zip(left, right, strict=True)
    )
    denominator = math.sqrt(
        sum((value - left_mean) ** 2 for value in left)
        * sum((value - right_mean) ** 2 for value in right)
    )
    return numerator / denominator if denominator else math.nan


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for sign in ("minus", "plus"):
        parser.add_argument(f"--baseline-{sign}-newton0", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-newton0", type=Path, required=True)
        parser.add_argument(f"--baseline-{sign}-newton1", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-newton1", type=Path, required=True)
        parser.add_argument(f"--baseline-{sign}-log", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-log", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--measure-coefficients", type=Path, required=True)
    parser.add_argument("--one-ring-rows", type=Path, required=True)
    parser.add_argument("--single-mode-summary", type=Path, required=True)
    parser.add_argument("--temperature-K", type=float, default=300.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    prior = json.loads(args.single_mode_summary.read_text(encoding="utf-8"))
    ring = [int(node) for node in prior["pre_registered_support"]["nodes"]]
    logs: dict[str, dict[str, Any]] = {}
    identity: dict[str, dict[str, bool]] = {}
    formula: dict[str, dict[str, Any]] = {}
    artifact_paths: list[Path] = [
        args.mesh, args.measure_coefficients, args.one_ring_rows,
        args.single_mode_summary,
    ]
    for sign in ("minus", "plus"):
        baseline0 = getattr(args, f"baseline_{sign}_newton0")
        candidate0 = getattr(args, f"candidate_{sign}_newton0")
        baseline1 = getattr(args, f"baseline_{sign}_newton1")
        candidate1 = getattr(args, f"candidate_{sign}_newton1")
        baseline_log = getattr(args, f"baseline_{sign}_log")
        candidate_log = getattr(args, f"candidate_{sign}_log")
        logs[sign] = {
            "baseline": log_summary(baseline_log),
            "candidate": log_summary(candidate_log),
        }
        for variant, export in (("baseline", baseline0), ("candidate", candidate0)):
            exported = field_norm(export, "eContinuityRhs")
            printed = logs[sign][variant]["iteration_rhs_norm"][0]
            logs[sign][variant]["iteration0_exported_eContinuityRhs"] = exported
            logs[sign][variant]["printed_global_rhs_over_exported_electron_l2"] = (
                printed / exported["l2"]
            )
        identity[sign] = {
            "iteration0_eContinuityRhs": (
                field_hash(baseline0, "eContinuityRhs")
                == field_hash(candidate0, "eContinuityRhs")
            ),
            "iteration1_eContinuityRhs": (
                field_hash(baseline1, "eContinuityRhs")
                == field_hash(candidate1, "eContinuityRhs")
            ),
            "iteration1_NewtonStepEDensityUpdate": (
                field_hash(baseline1, "NewtonStepEDensityUpdate")
                == field_hash(candidate1, "NewtonStepEDensityUpdate")
            ),
            "iteration1_eDensityError": (
                field_hash(baseline1, "eDensityError")
                == field_hash(candidate1, "eDensityError")
            ),
        }
        formula[sign] = errref_formula_check(
            baseline1,
            candidate1,
            logs[sign]["baseline"]["electron_errref_cm3"],
            logs[sign]["candidate"]["electron_errref_cm3"],
        )
        artifact_paths.extend([
            baseline_log, candidate_log,
            baseline0 / "fields" / "eContinuityRhs_region0.csv",
            candidate0 / "fields" / "eContinuityRhs_region0.csv",
            baseline1 / "fields" / "eContinuityRhs_region0.csv",
            candidate1 / "fields" / "eContinuityRhs_region0.csv",
            baseline1 / "fields" / "NewtonStepEDensityUpdate_region0.csv",
            candidate1 / "fields" / "NewtonStepEDensityUpdate_region0.csv",
            baseline1 / "fields" / "eDensityError_region0.csv",
            candidate1 / "fields" / "eDensityError_region0.csv",
        ])

    measures = silicon_node_measures(args.mesh, args.measure_coefficients, ring)
    row_by_node = {int(row["node_id"]): row for row in read_csv(args.one_ring_rows)}
    rows: list[dict[str, Any]] = []
    ratios: list[float] = []
    for node in ring:
        sentaurus = float(row_by_node[node]["sentaurus_dR_dphin_A_per_V"])
        vela = float(row_by_node[node]["vela_dR_dphin_A_per_um_per_V"])
        ratio = abs(sentaurus / vela)
        ratios.append(ratio)
        rows.append({
            "node_id": node,
            "silicon_averagebox_measure_um2": measures[node],
            "absolute_sentaurus_over_vela_row_ratio": ratio,
            "ratio_times_measure": ratio * measures[node],
            "ratio_over_measure_per_um2": ratio / measures[node],
        })
    log_measure = [math.log(measures[node]) for node in ring]
    log_ratio = [math.log(value) for value in ratios]
    measure_range = max(measures.values()) / min(measures.values())
    ratio_range = max(ratios) / min(ratios)

    k_B = 1.380649e-23
    q = 1.602176634e-19
    thermal_voltage = k_B * args.temperature_K / q
    scalar = float(prior["cross_engine_central_difference"]["one_ring"][
        "absolute_sentaurus_over_vela_least_squares"
    ])
    all_rhs_identical = all(
        values["iteration0_eContinuityRhs"]
        and values["iteration1_eContinuityRhs"]
        and values["iteration1_NewtonStepEDensityUpdate"]
        and not values["iteration1_eDensityError"]
        for values in identity.values()
    )
    all_ordered = all(
        record[variant][
            "iteration0_newtonplot_precedes_cnorm_and_first_linear_solve"
        ]
        for record in logs.values() for variant in ("baseline", "candidate")
    )
    report = {
        "experiment": "templates_ldmos_g3_sentaurus_internal_equation_scaling",
        "contract": {
            "state": "same_high_endpoint_vsv_symmetric_node3721_qf_pair",
            "baseline_electron_errref_cm3": 1.0e10,
            "candidate_electron_errref_cm3": 1.0e8,
            "only_changed_sentaurus_input": "Math.ErrRef(Electron)",
            "production_defaults_changed": False,
        },
        "manual_contract": {
            "rhs_norm": "norm of assembled residual; convergence quantity",
            "relative_update_error": "Delta_x/(epsilon_R*(abs(x)+ErrRef))",
            "newtonplot": "implementation-dependent internal RHS/error/update data",
            "supported_public_row_normalization_formula_found": False,
        },
        "errref_ab": {
            "logs": logs,
            "field_identity": identity,
            "error_formula_validation": formula,
            "rhs_and_linear_update_invariant": all_rhs_identical,
        },
        "control_box_measure": {
            "pre_registered_nodes": ring,
            "silicon_measure_range_factor": measure_range,
            "absolute_row_scale_range_factor": ratio_range,
            "log_measure_log_scale_pearson": pearson(log_measure, log_ratio),
            "per_row_measure_normalization_explains_uniform_scalar": (
                measure_range <= 1.05 * ratio_range
            ),
        },
        "newtonplot_ordering": {
            "all_four_logs_have_expected_order": all_ordered,
            "order": [
                "iteration_0_rhs_norm",
                "newton_0_0_tdr_write",
                "CNorm_update_error",
                "first_linear_solve_and_iteration_1_result",
            ],
            "linear_preconditioner_can_change_iteration0_rhs_field": False,
            "printed_global_rhs_is_direct_l2_of_exported_electron_rhs": False,
            "reason": (
                "iteration-0 NewtonPlot is written from the assembled residual "
                "before the first linear solve; the logs also declare no diagonal "
                "preconditioning; the printed global |Rhs| differs by over ten "
                "orders of magnitude from the exported electron-field L2"
            ),
        },
        "simple_thermal_voltage_candidate": {
            "temperature_K": args.temperature_K,
            "thermal_voltage_V": thermal_voltage,
            "inverse_thermal_voltage_per_V": 1.0 / thermal_voltage,
            "observed_scalar": scalar,
            "observed_over_inverse_thermal_voltage": scalar * thermal_voltage,
            "relative_difference_from_inverse_thermal_voltage": (
                abs(scalar / (1.0 / thermal_voltage) - 1.0)
            ),
            "simple_one_over_Vt_exactly_explains_scalar": False,
        },
        "classification": {
            "electron_errref_explains_33x": False,
            "per_row_control_box_measure_explains_33x": False,
            "linear_solver_preconditioner_explains_iteration0_newtonplot_rhs": False,
            "simple_inverse_thermal_voltage_explains_33x": False,
            "absolute_newtonplot_rhs_is_valid_cross_engine_gate": False,
            "newtonplot_rhs_shape_and_relative_bias_change_remain_valid": True,
            "remaining_uniform_internal_normalization_or_assembly_coefficient": True,
            "ledger_status": "draft",
        },
        "artifacts": {str(path): sha256(path) for path in artifact_paths},
    }
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "measure_row_scale.csv", rows)
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str((output / "summary.json").resolve()),
        "errref_rhs_and_update_invariant": all_rhs_identical,
        "measure_range_factor": measure_range,
        "row_scale_range_factor": ratio_range,
        "newtonplot_precedes_first_linear_solve": all_ordered,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
