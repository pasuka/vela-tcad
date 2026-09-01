#!/usr/bin/env python3
"""Align native Sentaurus electron equation balance with Vela frozen rows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


DEFAULT_NODES = (3721, 3974, 3973, 4091, 3727, 4021, 3771)
ELEMENTARY_CHARGE_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty alignment table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def cosine(left: list[float], right: list[float]) -> float | None:
    denominator = l2(left) * l2(right)
    if denominator == 0.0:
        return None
    return sum(a * b for a, b in zip(left, right, strict=True)) / denominator


def field_by_node(export_root: Path, field: str) -> dict[int, float]:
    rows = read_csv(export_root / "fields" / f"{field}_region0.csv")
    return {int(row["node_id"]): float(row["component0"]) for row in rows}


def coordinates(export_root: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(export_root / "nodes.csv")
    }


def carrier_rows(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            "x_um": float(row["x"]),
            "y_um": float(row["y"]),
            "electron_residual": float(row["electron_residual"]),
            "electron_flux_abs_sum": float(row["electron_flux_abs_sum"]),
            "electron_recombination": float(row["electron_recombination"]),
        }
        for row in read_csv(path)
    }


def diagnostic_state(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            "psi": float(row["psi"]),
            "phin": float(row["phin"]),
            "phip": float(row["phip"]),
        }
        for row in read_csv(path)
    }


def fitted_scale(target: list[float], source: list[float]) -> float:
    denominator = sum(value * value for value in source)
    if denominator == 0.0:
        raise ValueError("cannot fit a scale to a zero source vector")
    return sum(a * b for a, b in zip(target, source, strict=True)) / denominator


def current_scale_from_sg(path: Path) -> dict[str, float]:
    ratios: list[float] = []
    for row in read_csv(path):
        scaled = float(row["electron_flux"])
        physical = float(row["electron_particle_line_flux_per_m_s"])
        if abs(scaled) <= 1.0e-250 or abs(physical) <= 1.0e-250:
            continue
        ratio = physical / scaled
        if math.isfinite(ratio) and ratio > 0.0:
            ratios.append(ratio)
    if not ratios:
        raise ValueError(f"no finite electron flux scale in {path}")
    continuity_scale = statistics.median(ratios)
    relative_spread = max(
        abs(value / continuity_scale - 1.0) for value in ratios
    )
    return {
        "continuity_particle_scale_per_m_s": continuity_scale,
        "current_scale_A_per_um_per_scaled_residual": (
            continuity_scale * ELEMENTARY_CHARGE_C * 1.0e-6
        ),
        "sample_count": len(ratios),
        "maximum_relative_spread": relative_spread,
    }


def point_from_summary(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    points = payload.get("points", [])
    if len(points) != 1:
        raise ValueError(f"expected exactly one endpoint in {path}")
    return points[0]


def analyze_endpoint(
    *,
    name: str,
    sentaurus_export: Path,
    vela_carrier_csv: Path,
    vela_sg_csv: Path,
    state_feedback_summary: Path,
    selected_nodes: tuple[int, ...],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    sent_coordinates = coordinates(sentaurus_export)
    sent_rhs = field_by_node(sentaurus_export, "eContinuityRhs")
    vela = carrier_rows(vela_carrier_csv)
    point = point_from_summary(state_feedback_summary)
    scale = current_scale_from_sg(vela_sg_csv)
    current_scale = scale["current_scale_A_per_um_per_scaled_residual"]
    terminal = float(point["sentaurus_terminal_A_per_um"])

    missing = [node for node in selected_nodes if node not in sent_rhs or node not in vela]
    if missing:
        raise ValueError(f"{name}: missing selected nodes {missing}")

    rows: list[dict[str, Any]] = []
    sent_selected: list[float] = []
    vela_selected: list[float] = []
    maximum_coordinate_error = 0.0
    for node in selected_nodes:
        sx, sy = sent_coordinates[node]
        coordinate_error = math.hypot(sx - vela[node]["x_um"], sy - vela[node]["y_um"])
        maximum_coordinate_error = max(maximum_coordinate_error, coordinate_error)
        sent_value = sent_rhs[node]
        vela_scaled = vela[node]["electron_residual"]
        vela_physical = vela_scaled * current_scale
        sent_selected.append(sent_value)
        vela_selected.append(vela_scaled)
        rows.append({
            "bias": name,
            "node_id": node,
            "x_um": sx,
            "y_um": sy,
            "coordinate_error_um": coordinate_error,
            "sentaurus_e_continuity_rhs_A": sent_value,
            "vela_e_continuity_residual_scaled": vela_scaled,
            "vela_e_continuity_residual_A_per_um": vela_physical,
            "vela_flux_abs_sum_scaled": vela[node]["electron_flux_abs_sum"],
            "vela_row_cancellation_fraction": abs(vela_scaled) / max(
                vela[node]["electron_flux_abs_sum"], 1.0e-300
            ),
            "vela_recombination_scaled": vela[node]["electron_recombination"],
        })

    sent_all = list(sent_rhs.values())
    vela_all = [row["electron_residual"] for row in vela.values()]
    sent_l2 = l2(sent_selected)
    vela_l2_scaled = l2(vela_selected)
    vela_l2_physical = vela_l2_scaled * current_scale
    vela_global_l2_scaled = l2(vela_all)
    summary = {
        "bias": name,
        "selected_nodes": list(selected_nodes),
        "maximum_coordinate_error_um": maximum_coordinate_error,
        "sentaurus_native": {
            "seven_node_l2_A": sent_l2,
            "seven_node_max_abs_A": max(abs(value) for value in sent_selected),
            "silicon_l2_A": l2(sent_all),
            "silicon_max_abs_A": max(abs(value) for value in sent_all),
            "seven_node_l2_over_terminal": sent_l2 / abs(terminal),
        },
        "vela_same_sentaurus_state": {
            "seven_node_l2_scaled": vela_l2_scaled,
            "seven_node_l2_A_per_um": vela_l2_physical,
            "seven_node_max_abs_A_per_um": max(
                abs(value * current_scale) for value in vela_selected
            ),
            "silicon_l2_scaled": vela_global_l2_scaled,
            "silicon_l2_A_per_um": vela_global_l2_scaled * current_scale,
            "seven_node_energy_fraction_of_silicon": (
                vela_l2_scaled / max(vela_global_l2_scaled, 1.0e-300)
            ) ** 2,
            "seven_node_l2_over_terminal": vela_l2_physical / abs(terminal),
            "terminal_current_A_per_um": float(
                point["variants"]["SSS"]["current_A_per_um"]
            ),
            "terminal_ratio_to_sentaurus": float(
                point["variants"]["SSS"]["ratio_to_sentaurus"]
            ),
        },
        "sentaurus_terminal_A_per_um": terminal,
        "scaling": scale,
        "artifacts": {
            "sentaurus_rhs": str(
                sentaurus_export / "fields" / "eContinuityRhs_region0.csv"
            ),
            "vela_carrier_terms": str(vela_carrier_csv),
            "vela_sg_edges": str(vela_sg_csv),
        },
    }
    return summary, rows


def analyze_vsv_endpoint(
    *,
    name: str,
    sentaurus_newton_export: Path,
    sentaurus_loaded_export: Path,
    target_state_csv: Path,
    vela_carrier_csv: Path,
    vela_sg_csv: Path,
    selected_nodes: tuple[int, ...],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    sent_coordinates = coordinates(sentaurus_newton_export)
    sent_rhs = field_by_node(sentaurus_newton_export, "eContinuityRhs")
    loaded = {
        "psi": field_by_node(sentaurus_loaded_export, "ElectrostaticPotential"),
        "phin": field_by_node(sentaurus_loaded_export, "eQuasiFermiPotential"),
        "phip": field_by_node(sentaurus_loaded_export, "hQuasiFermiPotential"),
    }
    target = diagnostic_state(target_state_csv)
    vela = carrier_rows(vela_carrier_csv)
    scale = current_scale_from_sg(vela_sg_csv)
    current_scale = scale["current_scale_A_per_um_per_scaled_residual"]

    missing = [
        node for node in selected_nodes
        if node not in sent_rhs or node not in vela or node not in target
        or any(node not in values for values in loaded.values())
    ]
    if missing:
        raise ValueError(f"{name}: missing selected VSV nodes {missing}")

    sent_selected = [sent_rhs[node] for node in selected_nodes]
    vela_selected = [
        vela[node]["electron_residual"] * current_scale for node in selected_nodes
    ]
    scale_fit = fitted_scale(sent_selected, vela_selected)
    fit_residual = [
        sent - scale_fit * vela
        for sent, vela in zip(sent_selected, vela_selected, strict=True)
    ]
    sent_l2 = l2(sent_selected)
    vela_l2 = l2(vela_selected)

    rows: list[dict[str, Any]] = []
    maximum_coordinate_error = 0.0
    for node, sent_value, vela_value in zip(
        selected_nodes, sent_selected, vela_selected, strict=True
    ):
        sx, sy = sent_coordinates[node]
        coordinate_error = math.hypot(
            sx - vela[node]["x_um"], sy - vela[node]["y_um"]
        )
        maximum_coordinate_error = max(maximum_coordinate_error, coordinate_error)
        rows.append({
            "bias": name,
            "node_id": node,
            "x_um": sx,
            "y_um": sy,
            "coordinate_error_um": coordinate_error,
            "loaded_psi_error_V": loaded["psi"][node] - target[node]["psi"],
            "loaded_phin_error_V": loaded["phin"][node] - target[node]["phin"],
            "loaded_phip_error_V": loaded["phip"][node] - target[node]["phip"],
            "sentaurus_e_continuity_rhs_A": sent_value,
            "vela_e_continuity_residual_A_per_um": vela_value,
            "signed_sentaurus_over_vela": (
                sent_value / vela_value if vela_value != 0.0 else None
            ),
            "sentaurus_after_scalar_fit_residual_A": (
                sent_value - scale_fit * vela_value
            ),
        })

    sent_all = list(sent_rhs.values())
    vela_all = [
        vela[node]["electron_residual"] * current_scale
        for node in sent_rhs if node in vela
    ]
    sent_common = [sent_rhs[node] for node in sent_rhs if node in vela]
    global_fit = fitted_scale(sent_common, vela_all)
    global_fit_residual = [
        sent - global_fit * value
        for sent, value in zip(sent_common, vela_all, strict=True)
    ]
    state_errors = {
        field: max(
            abs(loaded[field][node] - target[node][field])
            for node in selected_nodes
        )
        for field in ("psi", "phin", "phip")
    }
    summary = {
        "bias": name,
        "selected_nodes": list(selected_nodes),
        "maximum_coordinate_error_um": maximum_coordinate_error,
        "loaded_state_seven_node_max_abs_error_V": state_errors,
        "sentaurus_vsv": {
            "seven_node_l2_A": sent_l2,
            "seven_node_max_abs_A": max(abs(value) for value in sent_selected),
            "silicon_l2_A": l2(sent_all),
        },
        "vela_vsv": {
            "seven_node_l2_A_per_um": vela_l2,
            "seven_node_max_abs_A_per_um": max(abs(value) for value in vela_selected),
            "silicon_l2_A_per_um": l2(vela_all),
        },
        "cross_engine_row_mode": {
            "cosine": cosine(sent_selected, vela_selected),
            "signed_sentaurus_over_vela_least_squares": scale_fit,
            "absolute_sentaurus_over_vela_least_squares": abs(scale_fit),
            "sentaurus_l2_over_vela_l2": sent_l2 / vela_l2,
            "relative_l2_after_scalar_fit": l2(fit_residual) / sent_l2,
            "maximum_abs_after_scalar_fit_over_sentaurus_l2": (
                max(abs(value) for value in fit_residual) / sent_l2
            ),
        },
        "whole_silicon_context": {
            "cosine": cosine(sent_common, vela_all),
            "signed_sentaurus_over_vela_least_squares": global_fit,
            "relative_l2_after_scalar_fit": l2(global_fit_residual) / l2(sent_common),
        },
        "scaling": scale,
        "artifacts": {
            "sentaurus_rhs": str(
                sentaurus_newton_export / "fields" / "eContinuityRhs_region0.csv"
            ),
            "sentaurus_loaded_state": str(sentaurus_loaded_export),
            "target_state": str(target_state_csv),
            "vela_carrier_terms": str(vela_carrier_csv),
            "vela_sg_edges": str(vela_sg_csv),
        },
    }
    return summary, rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for endpoint in ("low", "high"):
        parser.add_argument(f"--{endpoint}-sentaurus-export", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-vela-carrier", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-vela-sg", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-state-feedback-summary", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-vsv-sentaurus-export", type=Path)
        parser.add_argument(f"--{endpoint}-vsv-loaded-export", type=Path)
        parser.add_argument(f"--{endpoint}-vsv-state", type=Path)
        parser.add_argument(f"--{endpoint}-vsv-vela-carrier", type=Path)
        parser.add_argument(f"--{endpoint}-vsv-vela-sg", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nodes", default=",".join(map(str, DEFAULT_NODES)))
    args = parser.parse_args()

    selected = tuple(int(value) for value in args.nodes.split(",") if value.strip())
    endpoints: dict[str, dict[str, Any]] = {}
    rows: list[dict[str, Any]] = []
    input_paths: list[Path] = []
    for name in ("low", "high"):
        sentaurus_export = getattr(args, f"{name}_sentaurus_export").resolve()
        vela_carrier = getattr(args, f"{name}_vela_carrier").resolve()
        vela_sg = getattr(args, f"{name}_vela_sg").resolve()
        feedback = getattr(args, f"{name}_state_feedback_summary").resolve()
        endpoint, endpoint_rows = analyze_endpoint(
            name=name,
            sentaurus_export=sentaurus_export,
            vela_carrier_csv=vela_carrier,
            vela_sg_csv=vela_sg,
            state_feedback_summary=feedback,
            selected_nodes=selected,
        )
        endpoints[name] = endpoint
        rows.extend(endpoint_rows)
        input_paths.extend([
            sentaurus_export / "nodes.csv",
            sentaurus_export / "fields" / "eContinuityRhs_region0.csv",
            vela_carrier,
            vela_sg,
            feedback,
        ])

    low_vela = [
        row["vela_e_continuity_residual_scaled"] for row in rows if row["bias"] == "low"
    ]
    high_vela = [
        row["vela_e_continuity_residual_scaled"] for row in rows if row["bias"] == "high"
    ]
    native_roundoff = all(
        endpoint["sentaurus_native"]["seven_node_l2_over_terminal"] <= 1.0e-8
        for endpoint in endpoints.values()
    )
    vela_mismatch = all(
        endpoint["vela_same_sentaurus_state"]["seven_node_l2_over_terminal"] >= 1.0e-5
        for endpoint in endpoints.values()
    )
    report = {
        "schema": "vela.templates_ldmos.g3_sentaurus_equation_balance.v2",
        "selected_nodes": list(selected),
        "endpoints": endpoints,
        "bias_growth": {
            "vela_seven_node_l2_high_over_low": l2(high_vela) / l2(low_vela),
            "vela_seven_node_residual_cosine": cosine(low_vela, high_vela),
        },
        "classification": {
            "sentaurus_native_balance_is_roundoff": native_roundoff,
            "vela_same_state_local_balance_mismatch": vela_mismatch,
            "same_state_local_operator_mismatch_confirmed": (
                native_roundoff and vela_mismatch
            ),
            "exact_qf_drive_rule_identified": False,
            "reason": (
                "A converged equation balance constrains only each assembled row sum. "
                "Because the native Sentaurus RHS is at roundoff and the TDR current "
                "density is a vertex post-processing field, these data prove a local "
                "same-state operator mismatch but do not uniquely recover element-edge "
                "QF-drive or density averaging."
            ),
            "next_decisive_probe": (
                "Evaluate Sentaurus NewtonPlot on the identical VSV state or obtain a "
                "supported element-edge assembly-current export; do not tune mobility, "
                "threshold, IALMob, predictor, or production defaults."
            ),
            "ledger_status": "draft",
        },
        "provenance": {
            str(path): sha256(path) for path in input_paths
        },
    }

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    vsv_argument_names = [
        f"{endpoint}_vsv_{suffix}"
        for endpoint in ("low", "high")
        for suffix in (
            "sentaurus_export", "loaded_export", "state", "vela_carrier", "vela_sg"
        )
    ]
    vsv_values = [getattr(args, name) for name in vsv_argument_names]
    if any(value is not None for value in vsv_values):
        if not all(value is not None for value in vsv_values):
            missing = [
                name.replace("_", "-") for name, value in zip(
                    vsv_argument_names, vsv_values, strict=True
                ) if value is None
            ]
            raise ValueError(f"incomplete VSV argument group; missing {missing}")
        vsv_endpoints: dict[str, dict[str, Any]] = {}
        vsv_rows: list[dict[str, Any]] = []
        for name in ("low", "high"):
            endpoint, endpoint_rows = analyze_vsv_endpoint(
                name=name,
                sentaurus_newton_export=getattr(
                    args, f"{name}_vsv_sentaurus_export"
                ).resolve(),
                sentaurus_loaded_export=getattr(
                    args, f"{name}_vsv_loaded_export"
                ).resolve(),
                target_state_csv=getattr(args, f"{name}_vsv_state").resolve(),
                vela_carrier_csv=getattr(
                    args, f"{name}_vsv_vela_carrier"
                ).resolve(),
                vela_sg_csv=getattr(args, f"{name}_vsv_vela_sg").resolve(),
                selected_nodes=selected,
            )
            vsv_endpoints[name] = endpoint
            vsv_rows.extend(endpoint_rows)
            input_paths.extend([
                getattr(args, f"{name}_vsv_sentaurus_export").resolve()
                / "fields" / "eContinuityRhs_region0.csv",
                getattr(args, f"{name}_vsv_loaded_export").resolve()
                / "fields" / "ElectrostaticPotential_region0.csv",
                getattr(args, f"{name}_vsv_state").resolve(),
                getattr(args, f"{name}_vsv_vela_carrier").resolve(),
                getattr(args, f"{name}_vsv_vela_sg").resolve(),
            ])
        scales = [
            endpoint["cross_engine_row_mode"][
                "absolute_sentaurus_over_vela_least_squares"
            ]
            for endpoint in vsv_endpoints.values()
        ]
        maximum_state_error = max(
            value
            for endpoint in vsv_endpoints.values()
            for value in endpoint["loaded_state_seven_node_max_abs_error_V"].values()
        )
        scale_spread = max(scales) / min(scales) - 1.0
        row_modes_aligned = all(
            abs(endpoint["cross_engine_row_mode"]["cosine"]) >= 0.999
            and endpoint["cross_engine_row_mode"]["relative_l2_after_scalar_fit"] <= 0.02
            for endpoint in vsv_endpoints.values()
        )
        report["vsv_endpoints"] = vsv_endpoints
        report["vsv_classification"] = {
            "seven_node_state_exact_within_1e_12_V": maximum_state_error <= 1.0e-12,
            "assembled_seven_node_row_mode_aligned_up_to_scalar": row_modes_aligned,
            "absolute_scale_bias_stable_within_1_percent": scale_spread <= 0.01,
            "absolute_scale_high_over_low_minus_one": scale_spread,
            "individual_element_edge_rule_identified": False,
            "interpretation": (
                "On two nonzero identical VSV states, the seven-node Sentaurus and "
                "Vela electron-row vectors are antiparallel by sign convention and "
                "agree in shape after one bias-stable scalar. This strongly excludes "
                "a different local QF-drive spatial mode as the dominant cause, but "
                "does not by itself identify whether the remaining scalar belongs to "
                "NewtonPlot RHS normalization, 2-D width/current convention, or a "
                "uniform assembly coefficient."
            ),
            "next_decisive_probe": (
                "Calibrate the Sentaurus NewtonPlot eContinuityRhs scalar with a "
                "controlled single-mode QF perturbation or a supported assembly-edge "
                "current export before changing SG, mobility, IALMob, predictor, or "
                "production defaults."
            ),
            "ledger_status": "draft",
        }
        report["classification"]["vsv_row_mode_aligned_up_to_scalar"] = (
            row_modes_aligned
        )
        report["classification"]["next_decisive_probe"] = report[
            "vsv_classification"
        ]["next_decisive_probe"]
        write_csv(output / "vsv_seven_node_equation_balance.csv", vsv_rows)

    report["provenance"] = {str(path): sha256(path) for path in input_paths}
    write_csv(output / "seven_node_equation_balance.csv", rows)
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str(output / "summary.json"),
        "classification": report["classification"],
        "vsv_classification": report.get("vsv_classification"),
        "bias_growth": report["bias_growth"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
