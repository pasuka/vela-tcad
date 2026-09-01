#!/usr/bin/env python3
"""Audit Fermi generalized-Einstein edge averages on frozen LDMOS rows."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

DEFAULT_NODES = (3721, 3974, 3973, 4091, 3727, 4021, 3771)
KB_J_K = 1.380649e-23
Q_C = 1.602176634e-19
BEDNARCZYK_COEFFICIENT = 0.75 * math.sqrt(math.pi)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fermi_half(eta: float) -> float:
    if not math.isfinite(eta):
        return 0.0 if eta < 0.0 else float("inf")
    if eta < -40.0:
        return math.exp(eta)
    exponential = math.exp(max(-eta, -700.0))
    shifted = eta + 1.0
    gaussian = math.exp(-0.17 * shifted * shifted)
    value = eta**4 + 50.0 + 33.6 * eta * (1.0 - 0.68 * gaussian)
    denominator = exponential + BEDNARCZYK_COEFFICIENT * value**-0.375
    return 1.0 / denominator


def fermi_half_derivative(eta: float) -> float:
    if not math.isfinite(eta):
        return 0.0
    if eta < -40.0:
        return math.exp(eta)
    exponential = math.exp(max(-eta, -700.0))
    shifted = eta + 1.0
    gaussian = math.exp(-0.17 * shifted * shifted)
    bracket = 1.0 - 0.68 * gaussian
    bracket_derivative = 0.2312 * shifted * gaussian
    value = eta**4 + 50.0 + 33.6 * eta * bracket
    value_derivative = 4.0 * eta**3 + 33.6 * (
        bracket + eta * bracket_derivative
    )
    power = value**-0.375
    denominator = exponential + BEDNARCZYK_COEFFICIENT * power
    denominator_derivative = -exponential - (
        0.375 * BEDNARCZYK_COEFFICIENT * value**-1.375 * value_derivative
    )
    return max(0.0, -denominator_derivative / denominator**2)


def local_generalized_factor(eta: float) -> float:
    derivative = fermi_half_derivative(eta)
    return max(1.0, fermi_half(eta) / derivative) if derivative > 0.0 else 1.0


def inverse_fermi_half(value: float) -> float:
    if not (value > 0.0 and math.isfinite(value)):
        raise ValueError("inverse Fermi argument must be positive and finite")
    lower = -500.0
    upper = max(2.0, (value / 0.752252778063675) ** (2.0 / 3.0) * 1.5)
    while fermi_half(upper) < value:
        upper *= 2.0
    for _ in range(100):
        midpoint = 0.5 * (lower + upper)
        if fermi_half(midpoint) < value:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def bernoulli(value: float) -> float:
    if abs(value) < 1.0e-4:
        value2 = value * value
        return 1.0 - 0.5 * value + value2 / 12.0 - value2 * value2 / 720.0
    if value > 50.0:
        return value * math.exp(-value)
    if value < -50.0:
        return -value
    return value / math.expm1(value)


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def factor_variants(row: dict[str, str]) -> dict[str, float]:
    eta0 = float(row["electron_eta0"])
    eta1 = float(row["electron_eta1"])
    density0 = float(row["electron_density0_m3"])
    density1 = float(row["electron_density1_m3"])
    production = float(row["electron_generalized_einstein_factor"])
    local0 = local_generalized_factor(eta0)
    local1 = local_generalized_factor(eta1)
    midpoint = local_generalized_factor(0.5 * (eta0 + eta1))
    arithmetic = 0.5 * (local0 + local1)
    geometric = math.sqrt(local0 * local1)
    harmonic = 2.0 / (1.0 / local0 + 1.0 / local1)
    density_weighted = (
        density0 * local0 + density1 * local1
    ) / max(density0 + density1, 1.0e-300)
    log_density_eta = inverse_fermi_half(
        math.sqrt(fermi_half(eta0) * fermi_half(eta1))
    )
    log_density_midpoint = local_generalized_factor(log_density_eta)
    # Simpson integration of the local factor over eta. It is symmetric and
    # provides a smooth integral-average diagnostic distinct from the secant.
    integral_average = (local0 + 4.0 * midpoint + local1) / 6.0
    return {
        "production_secant": production,
        "classical": 1.0,
        "eta_midpoint": midpoint,
        "endpoint_arithmetic": arithmetic,
        "endpoint_geometric": geometric,
        "endpoint_harmonic": harmonic,
        "density_weighted": density_weighted,
        "log_density_midpoint": log_density_midpoint,
        "eta_integral_simpson": integral_average,
        "endpoint_min_bound": min(local0, local1),
        "endpoint_max_bound": max(local0, local1),
    }


def candidate_flux(row: dict[str, str], factor: float, thermal_voltage: float) -> float:
    production_flux = float(row["electron_flux"])
    production_factor = float(row["electron_generalized_einstein_factor"])
    production_argument = float(row["electron_bernoulli_argument"])
    production_qf_argument = float(row["electron_quasi_fermi_argument"])
    production_kernel = (
        production_factor
        * bernoulli(production_argument)
        * math.expm1(production_qf_argument)
    )
    if abs(production_kernel) <= 1.0e-300:
        return 0.0
    drift = float(row["electron_drift_potential_V"])
    qf_drop = float(row["phin1_V"]) - float(row["phin0_V"])
    candidate_kernel = (
        factor
        * bernoulli(drift / (thermal_voltage * factor))
        * math.expm1(qf_drop / (thermal_voltage * factor))
    )
    return production_flux * candidate_kernel / production_kernel


def signed_at_node(row: dict[str, str], node: int, flux: float) -> float:
    if int(row["node0"]) == node:
        return flux
    if int(row["node1"]) == node:
        return -flux
    raise ValueError(f"edge {row['edge_id']} is not incident to node {node}")


def norm(values: list[float]) -> dict[str, float]:
    return {
        "l2": math.sqrt(sum(value * value for value in values)),
        "maximum_abs": max(map(abs, values), default=0.0),
    }


def audit_state(
    rows: list[dict[str, str]], selected_nodes: tuple[int, ...], temperature_K: float,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    selected = set(selected_nodes)
    incident = [
        row for row in rows
        if int(row["node0"]) in selected or int(row["node1"]) in selected
    ]
    thermal_voltage = KB_J_K * temperature_K / Q_C
    variant_names = list(factor_variants(next(
        row for row in incident if abs(float(row["electron_flux"])) > 1.0e-300
    )))
    contributions: dict[str, dict[int, list[float]]] = {
        variant: {node: [] for node in selected_nodes}
        for variant in variant_names
    }
    factor_ratios: dict[str, list[float]] = {name: [] for name in variant_names}
    edge_output: list[dict[str, Any]] = []
    maximum_production_factor_error = 0.0
    maximum_production_flux_error = 0.0

    for row in incident:
        production_flux = float(row["electron_flux"])
        active = abs(production_flux) > 1.0e-300
        variants = factor_variants(row) if active else {
            name: 1.0 for name in variant_names
        }
        if active:
            eta0 = float(row["electron_eta0"])
            eta1 = float(row["electron_eta1"])
            n0 = float(row["electron_density0_m3"])
            n1 = float(row["electron_density1_m3"])
            delta_log_n = math.log(n1 / n0)
            reconstructed = (
                (eta1 - eta0) / delta_log_n
                if abs(delta_log_n) > 1.0e-10
                else local_generalized_factor(0.5 * (eta0 + eta1))
            )
            maximum_production_factor_error = max(
                maximum_production_factor_error,
                abs(reconstructed - variants["production_secant"])
                / max(abs(variants["production_secant"]), 1.0e-300),
            )
        fluxes = {
            name: candidate_flux(row, value, thermal_voltage) if active else 0.0
            for name, value in variants.items()
        }
        if active:
            maximum_production_flux_error = max(
                maximum_production_flux_error,
                abs(fluxes["production_secant"] - production_flux)
                / max(abs(production_flux), 1.0e-300),
            )
        production_factor = variants["production_secant"]
        if active:
            for name, value in variants.items():
                factor_ratios[name].append(value / production_factor)
        for node in (int(row["node0"]), int(row["node1"])):
            if node not in selected:
                continue
            for name, flux in fluxes.items():
                contributions[name][node].append(signed_at_node(row, node, flux))
            edge_output.append({
                "node_id": node,
                "edge_id": int(row["edge_id"]),
                "other_node": int(row["node1"]) if node == int(row["node0"])
                    else int(row["node0"]),
                "production_signed_flux": signed_at_node(row, node, production_flux),
                **{
                    f"{name}_factor": variants[name]
                    for name in variant_names
                },
                **{
                    f"{name}_signed_flux": signed_at_node(row, node, fluxes[name])
                    for name in variant_names
                },
            })

    residuals = {
        name: [sum(contributions[name][node]) for node in selected_nodes]
        for name in variant_names
    }
    norms = {name: norm(values) for name, values in residuals.items()}
    baseline = norms["production_secant"]
    variants_summary: dict[str, Any] = {}
    for name in variant_names:
        ratios = factor_ratios[name]
        node_ratios = [
            abs(candidate) / max(abs(production), 1.0e-300)
            for candidate, production in zip(
                residuals[name], residuals["production_secant"], strict=True
            )
        ]
        variants_summary[name] = {
            "residual": norms[name],
            "residual_over_production": {
                key: norms[name][key] / max(baseline[key], 1.0e-300)
                for key in baseline
            },
            "maximum_node_residual_ratio": max(node_ratios, default=0.0),
            "factor_over_production": {
                "minimum": min(ratios, default=1.0),
                "median": percentile(ratios, 0.5),
                "p95": percentile(ratios, 0.95),
                "maximum": max(ratios, default=1.0),
            },
            "node_residuals": [
                {"node_id": node, "value": value}
                for node, value in zip(selected_nodes, residuals[name], strict=True)
            ],
        }
    return {
        "selected_nodes": list(selected_nodes),
        "unique_incident_edges": len({int(row["edge_id"]) for row in incident}),
        "active_incident_edges": sum(
            abs(float(row["electron_flux"])) > 1.0e-300 for row in incident
        ),
        "maximum_production_factor_reconstruction_relative_error":
            maximum_production_factor_error,
        "maximum_production_flux_reconstruction_relative_error":
            maximum_production_flux_error,
        "variants": variants_summary,
    }, edge_output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, action="append", required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--temperature-K", type=float, default=300.0)
    parser.add_argument("--nodes", default=",".join(map(str, DEFAULT_NODES)))
    args = parser.parse_args()
    if len(args.state) != len(args.label) or len(args.state) != 2:
        raise ValueError("exactly two aligned states and labels are required")
    selected = tuple(int(value) for value in args.nodes.split(",") if value.strip())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    states: dict[str, Any] = {}
    for label, path in zip(args.label, args.state, strict=True):
        summary, edges = audit_state(read_csv(path), selected, args.temperature_K)
        states[label] = summary
        write_csv(output / f"{label}_edge_variants.csv", edges)

    variants = list(next(iter(states.values()))["variants"])
    combined: dict[str, Any] = {}
    for name in variants:
        l2_ratios = [
            states[label]["variants"][name]["residual_over_production"]["l2"]
            for label in args.label
        ]
        max_ratios = [
            states[label]["variants"][name]["residual_over_production"]["maximum_abs"]
            for label in args.label
        ]
        node_ratios = [
            states[label]["variants"][name]["maximum_node_residual_ratio"]
            for label in args.label
        ]
        combined[name] = {
            "l2_ratios": l2_ratios,
            "maximum_abs_ratios": max_ratios,
            "maximum_node_ratios": node_ratios,
            "worst_l2_ratio": max(l2_ratios),
            "worst_maximum_abs_ratio": max(max_ratios),
            "worst_node_ratio": max(node_ratios),
            "passes_fixed_state_gate": (
                name != "production_secant"
                and max(l2_ratios) <= 0.5
                and max(max_ratios) <= 0.5
                and max(node_ratios) <= 0.5
            ),
        }
    ranked = sorted(
        (name for name in variants if name != "production_secant"),
        key=lambda name: (
            combined[name]["worst_l2_ratio"],
            combined[name]["worst_maximum_abs_ratio"],
        ),
    )
    report = {
        "schema": "vela.templates_ldmos.g3_fermi_edge_average_ab.v1",
        "contract": {
            "mode": "read_only_frozen_state",
            "bias_points_V": [1.0, 1.1666666666666667],
            "state": "VSV (Vela psi, Sentaurus phin, Vela phip)",
            "carrier_couple": "qualified external AverageBox",
            "ialmob": "disabled",
            "predictor": "disabled",
            "production_solver_changed": False,
            "gate": "both endpoints L2, maximum and every selected-node residual <=0.5x production",
        },
        "states": states,
        "combined": combined,
        "ranked_nonproduction_variants": ranked,
        "passing_variants": [
            name for name in ranked if combined[name]["passes_fixed_state_gate"]
        ],
    }
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "ranked": [
            {"variant": name, **combined[name]} for name in ranked
        ],
        "passing_variants": report["passing_variants"],
        "summary": str((output / "summary.json").resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
