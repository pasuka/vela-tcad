#!/usr/bin/env python3
"""Decompose the production SG state on selected LDMOS interface rows."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any


DEFAULT_NODES = (3721, 3974, 3973, 4091, 3727, 4021, 3771)
FACTOR_NAMES = (
    "geometry",
    "mobility",
    "einstein",
    "bernoulli",
    "right_density",
    "qf_drive",
)


def runner_environment() -> dict[str, str]:
    environment = dict(os.environ)
    if os.name == "nt":
        environment["PATH"] = os.pathsep.join([
            r"D:\msys64\ucrt64\bin",
            r"D:\msys64\usr\bin",
            environment.get("PATH", ""),
        ])
    environment["VELA_LINEAR_SOLVER"] = "sparselu"
    return environment


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty edge-state table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run_edge_probe(runner: Path, source_config: Path, output: Path) -> Path:
    config = json.loads(source_config.read_text(encoding="utf-8"))
    if config.get("simulation_type") != "sg_edge_flux_probe":
        raise ValueError(f"not an SG edge probe: {source_config}")
    mobility = config.get("solver", {}).get("mobility", {})
    mobility_model = str(mobility.get("model", "")).lower()
    predictor = config.get("sweep", {}).get("predictor", False)
    if "ialmob" in mobility_model or predictor not in (False, None, "off"):
        raise ValueError("interface edge audit requires IALMob and predictor disabled")
    if config.get("mesh_geometry", {}).get("carrier_transport_couple_profile") != (
        "templates_ldmos_external_averagebox"
    ):
        raise ValueError("interface edge audit requires external AverageBox")
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "sg_edges_extended.csv"
    config_path = output / "sg_probe_extended.json"
    config["output_csv"] = str(csv_path.resolve())
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=output,
        env=runner_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    (output / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"SG edge probe failed ({completed.returncode}): {completed.stderr}"
        )
    return csv_path


def edge_factors(row: dict[str, str]) -> dict[str, float]:
    length = float(row["length_m"])
    if length <= 0.0:
        raise ValueError(f"non-positive edge length for edge {row['edge_id']}")
    return {
        "geometry": float(row["couple_m"]) / length,
        "mobility": float(row["electron_mobility_m2_V_s"]),
        "einstein": float(row["electron_generalized_einstein_factor"]),
        "bernoulli": float(row["electron_bernoulli_plus"]),
        "right_density": float(row["electron_density1_m3"]),
        "qf_drive": math.expm1(float(row["electron_quasi_fermi_argument"])),
    }


def factor_product(factors: dict[str, float]) -> float:
    result = 1.0
    for name in FACTOR_NAMES:
        result *= factors[name]
    return result


def signed_at_node(row: dict[str, str], node: int, flux: float) -> float:
    if int(row["node0"]) == node:
        return flux
    if int(row["node1"]) == node:
        return -flux
    raise ValueError(f"edge {row['edge_id']} is not incident to node {node}")


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def analyze_edges(
    low_rows: list[dict[str, str]],
    high_rows: list[dict[str, str]],
    selected_nodes: tuple[int, ...],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    low_by_edge = {int(row["edge_id"]): row for row in low_rows}
    high_by_edge = {int(row["edge_id"]): row for row in high_rows}
    if set(low_by_edge) != set(high_by_edge):
        raise ValueError("low/high SG edge sets differ")

    selected = set(selected_nodes)
    incident_ids = sorted(
        edge_id for edge_id, row in low_by_edge.items()
        if int(row["node0"]) in selected or int(row["node1"]) in selected
    )
    records: list[dict[str, Any]] = []
    node_contributions: dict[int, dict[str, list[float]]] = {
        node: {"low": [], "high": [], **{name: [] for name in FACTOR_NAMES}}
        for node in selected_nodes
    }
    maximum_ratio_error = 0.0
    inactive_zero_flux_edges = 0

    for edge_id in incident_ids:
        low = low_by_edge[edge_id]
        high = high_by_edge[edge_id]
        for invariant in ("node0", "node1", "length_m", "couple_m"):
            if float(low[invariant]) != float(high[invariant]):
                raise ValueError(f"edge {edge_id} changes invariant {invariant}")
        low_flux = float(low["electron_flux"])
        high_flux = float(high["electron_flux"])
        low_factors = edge_factors(low)
        high_factors = edge_factors(high)
        factorization_finite = all(
            math.isfinite(value)
            for value in (*low_factors.values(), *high_factors.values())
        )
        low_kernel = factor_product(low_factors) if factorization_finite else 0.0
        high_kernel = factor_product(high_factors) if factorization_finite else 0.0
        predicted_ratio = None
        actual_ratio = None
        ratio_error = None
        if abs(low_flux) > 1.0e-300 and abs(low_kernel) > 1.0e-300:
            actual_ratio = high_flux / low_flux
            predicted_ratio = high_kernel / low_kernel
            ratio_error = abs(predicted_ratio - actual_ratio) / max(
                abs(actual_ratio), 1.0e-300
            )
            maximum_ratio_error = max(maximum_ratio_error, ratio_error)
        if not factorization_finite and (
            abs(low_flux) > 1.0e-300 or abs(high_flux) > 1.0e-300
        ):
            raise ValueError(f"active edge {edge_id} has non-finite SG factors")
        if not factorization_finite:
            inactive_zero_flux_edges += 1
        if abs(low_kernel) > 1.0e-300:
            scale = low_flux / low_kernel
        elif abs(high_kernel) > 1.0e-300:
            scale = high_flux / high_kernel
        else:
            scale = 0.0
        held_low_flux: dict[str, float] = {}
        for name in FACTOR_NAMES:
            if factorization_finite:
                factors = dict(high_factors)
                factors[name] = low_factors[name]
                held_low_flux[name] = scale * factor_product(factors)
            else:
                held_low_flux[name] = 0.0

        endpoints = [int(low["node0"]), int(low["node1"])]
        for node in endpoints:
            if node not in selected:
                continue
            low_signed = signed_at_node(low, node, low_flux)
            high_signed = signed_at_node(high, node, high_flux)
            node_contributions[node]["low"].append(low_signed)
            node_contributions[node]["high"].append(high_signed)
            for name in FACTOR_NAMES:
                node_contributions[node][name].append(
                    signed_at_node(high, node, held_low_flux[name])
                )
            record = {
                "node_id": node,
                "edge_id": edge_id,
                "node0": endpoints[0],
                "node1": endpoints[1],
                "other_node": endpoints[1] if node == endpoints[0] else endpoints[0],
                "x0_m": float(low["x0"]),
                "y0_m": float(low["y0"]),
                "x1_m": float(low["x1"]),
                "y1_m": float(low["y1"]),
                "length_m": float(low["length_m"]),
                "couple_m": float(low["couple_m"]),
                "low_signed_electron_flux": low_signed,
                "high_signed_electron_flux": high_signed,
                "signed_flux_change": high_signed - low_signed,
                "low_psi_drop_V": float(low["psi1_V"]) - float(low["psi0_V"]),
                "high_psi_drop_V": float(high["psi1_V"]) - float(high["psi0_V"]),
                "low_phin_drop_V": float(low["phin1_V"]) - float(low["phin0_V"]),
                "high_phin_drop_V": float(high["phin1_V"]) - float(high["phin0_V"]),
                "low_mobility_field_V_m": float(low["electron_mobility_field_V_m"]),
                "high_mobility_field_V_m": float(high["electron_mobility_field_V_m"]),
                "low_bernoulli_argument": float(low["electron_bernoulli_argument"]),
                "high_bernoulli_argument": float(high["electron_bernoulli_argument"]),
                "low_qf_argument": float(low["electron_quasi_fermi_argument"]),
                "high_qf_argument": float(high["electron_quasi_fermi_argument"]),
                "actual_flux_ratio": actual_ratio,
                "factor_product_flux_ratio": predicted_ratio,
                "factor_ratio_relative_error": ratio_error,
            }
            for name in FACTOR_NAMES:
                denominator = low_factors[name]
                record[f"low_{name}"] = (
                    low_factors[name] if math.isfinite(low_factors[name]) else None
                )
                record[f"high_{name}"] = (
                    high_factors[name] if math.isfinite(high_factors[name]) else None
                )
                record[f"{name}_ratio"] = (
                    high_factors[name] / denominator
                    if factorization_finite and abs(denominator) > 1.0e-300
                    else None
                )
                record[f"high_flux_with_{name}_held_low"] = signed_at_node(
                    high, node, held_low_flux[name]
                )
            records.append(record)

    node_summaries: list[dict[str, Any]] = []
    low_residuals: list[float] = []
    high_residuals: list[float] = []
    held_low_residuals: dict[str, list[float]] = {name: [] for name in FACTOR_NAMES}
    for node in selected_nodes:
        contribution = node_contributions[node]
        low_residual = sum(contribution["low"])
        high_residual = sum(contribution["high"])
        low_residuals.append(low_residual)
        high_residuals.append(high_residual)
        held = {name: sum(contribution[name]) for name in FACTOR_NAMES}
        for name in FACTOR_NAMES:
            held_low_residuals[name].append(held[name])
        node_records = [row for row in records if row["node_id"] == node]
        node_summaries.append({
            "node_id": node,
            "incident_edges": len(node_records),
            "low_residual": low_residual,
            "high_residual": high_residual,
            "absolute_growth": abs(high_residual) / max(abs(low_residual), 1.0e-300),
            "low_cancellation_condition": sum(abs(value) for value in contribution["low"])
                / max(abs(low_residual), 1.0e-300),
            "high_cancellation_condition": sum(abs(value) for value in contribution["high"])
                / max(abs(high_residual), 1.0e-300),
            "high_residual_with_factor_held_low": held,
            "largest_high_edges": [
                {
                    "edge_id": row["edge_id"],
                    "other_node": row["other_node"],
                    "signed_flux": row["high_signed_electron_flux"],
                    "signed_flux_change": row["signed_flux_change"],
                }
                for row in sorted(
                    node_records,
                    key=lambda item: abs(item["high_signed_electron_flux"]),
                    reverse=True,
                )[:4]
            ],
        })

    low_norm = l2(low_residuals)
    high_norm = l2(high_residuals)
    held_norms = {name: l2(values) for name, values in held_low_residuals.items()}
    summary = {
        "schema": "vela.templates_ldmos.g3_interface_edge_state.v1",
        "selected_nodes": list(selected_nodes),
        "unique_incident_edges": len(incident_ids),
        "row_edge_records": len(records),
        "inactive_zero_flux_incident_edges": inactive_zero_flux_edges,
        "low_residual_l2": low_norm,
        "high_residual_l2": high_norm,
        "residual_l2_growth_high_over_low": high_norm / max(low_norm, 1.0e-300),
        "high_residual_l2_with_factor_held_low": held_norms,
        "held_low_l2_over_high": {
            name: value / max(high_norm, 1.0e-300)
            for name, value in held_norms.items()
        },
        "held_low_l2_over_low": {
            name: value / max(low_norm, 1.0e-300)
            for name, value in held_norms.items()
        },
        "maximum_factor_product_flux_ratio_relative_error": maximum_ratio_error,
        "nodes": node_summaries,
        "interpretation": (
            "Each held-low result is a diagnostic high-bias counterfactual; "
            "it does not change the production solver."
        ),
    }
    return summary, records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--low-config", type=Path, required=True)
    parser.add_argument("--high-config", type=Path, required=True)
    parser.add_argument("--low-target-config", type=Path)
    parser.add_argument("--high-target-config", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nodes", default=",".join(map(str, DEFAULT_NODES)))
    args = parser.parse_args()

    output = args.output_dir.resolve()
    low_csv = run_edge_probe(args.runner, args.low_config, output / "low")
    high_csv = run_edge_probe(args.runner, args.high_config, output / "high")
    selected = tuple(int(value) for value in args.nodes.split(",") if value.strip())
    bias_summary, bias_rows = analyze_edges(
        read_csv(low_csv), read_csv(high_csv), selected
    )
    write_csv(output / "bias_growth_interface_edge_state.csv", bias_rows)
    report: dict[str, Any] = {
        "schema": "vela.templates_ldmos.g3_interface_edge_state.v2",
        "bias_growth_vsv": bias_summary,
    }
    target_artifacts: dict[str, str] = {}
    if (args.low_target_config is None) != (args.high_target_config is None):
        raise ValueError("low/high target configs must be supplied together")
    if args.low_target_config is not None and args.high_target_config is not None:
        low_target_csv = run_edge_probe(
            args.runner, args.low_target_config, output / "low_target"
        )
        high_target_csv = run_edge_probe(
            args.runner, args.high_target_config, output / "high_target"
        )
        low_reclose, low_reclose_rows = analyze_edges(
            read_csv(low_csv), read_csv(low_target_csv), selected
        )
        high_reclose, high_reclose_rows = analyze_edges(
            read_csv(high_csv), read_csv(high_target_csv), selected
        )
        report["same_bias_vsv_to_vvv"] = {
            "vg_1p0": low_reclose,
            "vg_1p166667": high_reclose,
        }
        write_csv(output / "low_reclose_interface_edge_state.csv", low_reclose_rows)
        write_csv(output / "high_reclose_interface_edge_state.csv", high_reclose_rows)
        target_artifacts = {
            "low_target_extended_edges": str(low_target_csv),
            "high_target_extended_edges": str(high_target_csv),
            "low_reclose_row_edge_table": str(
                output / "low_reclose_interface_edge_state.csv"
            ),
            "high_reclose_row_edge_table": str(
                output / "high_reclose_interface_edge_state.csv"
            ),
        }
    report["contract"] = {
        "bias_points_V": [1.0, 1.1666666666666667],
        "bias_growth_state": "VSV (Vela psi, Sentaurus phin, Vela phip)",
        "same_bias_transition": "VSV to VVV when target configs are supplied",
        "ialmob": "disabled",
        "predictor": "disabled",
        "carrier_couples": "external AverageBox",
        "production_parameters_changed": False,
    }
    report["artifacts"] = {
        "low_extended_edges": str(low_csv),
        "high_extended_edges": str(high_csv),
        "bias_growth_row_edge_table": str(
            output / "bias_growth_interface_edge_state.csv"
        ),
        **target_artifacts,
    }
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    console = {
        "summary": str(output / "summary.json"),
        "bias_growth_vsv": {
            key: bias_summary[key]
            for key in (
                "residual_l2_growth_high_over_low",
                "held_low_l2_over_high",
                "maximum_factor_product_flux_ratio_relative_error",
            )
        },
    }
    if "same_bias_vsv_to_vvv" in report:
        console["same_bias_vsv_to_vvv"] = {
            name: {
                "initial_l2": value["low_residual_l2"],
                "target_l2": value["high_residual_l2"],
                "held_initial_l2_over_initial": value["held_low_l2_over_low"],
            }
            for name, value in report["same_bias_vsv_to_vvv"].items()
        }
    print(json.dumps(console, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
