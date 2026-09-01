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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for endpoint in ("low", "high"):
        parser.add_argument(f"--{endpoint}-sentaurus-export", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-vela-carrier", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-vela-sg", type=Path, required=True)
        parser.add_argument(f"--{endpoint}-state-feedback-summary", type=Path, required=True)
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
        "schema": "vela.templates_ldmos.g3_sentaurus_equation_balance.v1",
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
    write_csv(output / "seven_node_equation_balance.csv", rows)
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str(output / "summary.json"),
        "classification": report["classification"],
        "bias_growth": report["bias_growth"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
