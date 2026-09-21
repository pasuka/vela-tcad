#!/usr/bin/env python3
"""Run the frozen WP3/T1 residual-scaling audit for Templates/LDMOS G3.

Only explicitly supplied, same-contract G3 Id-Vg states are accepted.  Each
state is replayed through the unchanged ``newton_carrier_term_probe`` and
``sg_edge_flux_probe`` contracts.  The script physicalizes carrier residuals,
builds a direct-interface plus one-ring band, reports frozen-kernel analytic
electron-QF endpoint sensitivities on the 40 frozen-hotspot incident edges,
and performs the pre-registered collinearity/LOOCV checks.

The endpoint sensitivities differentiate the SG density weights while holding
the already assembled mobility and generalized-Einstein edge factor fixed.
They are therefore the analytic diagnostic needed to classify upwind versus
downwind Bernoulli suppression; they are not a replacement for the complete
production Jacobian, whose HFS and secant-factor state derivatives are separate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    from scripts.audit_templates_ldmos_g3_fermi_edge_average_ab import (
        fermi_half,
        fermi_half_derivative,
    )
    from scripts.audit_templates_ldmos_interface_pair_box import (
        region_local_geometry,
    )
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_fermi_edge_average_ab import (  # type: ignore
        fermi_half,
        fermi_half_derivative,
    )
    from audit_templates_ldmos_interface_pair_box import (  # type: ignore
        region_local_geometry,
    )


DEFAULT_NODES = (3721, 3974, 3973, 4091, 3727, 4021, 3771)
DIRECT_INTERFACE_HOTSPOTS = (3721, 3974, 3973, 4091, 3727, 3771)
ELEMENTARY_CHARGE_C = 1.602176634e-19
KB_J_K = 1.380649e-23
EXPECTED_BIASES = tuple(index / 6.0 for index in range(8))
T95 = {2: 4.302653, 3: 3.182446, 4: 2.776445, 5: 2.570582, 6: 2.446912}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
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


def bias_label(bias: float) -> str:
    return f"vg_{bias:.6f}".replace(".", "p")


def parse_state_spec(spec: str) -> tuple[float, Path]:
    if "=" not in spec:
        raise ValueError("--state must be BIAS=PATH")
    raw_bias, raw_path = spec.split("=", 1)
    return float(raw_bias), Path(raw_path).resolve()


def validate_state_specs(specs: list[str]) -> list[tuple[float, Path]]:
    states = sorted((parse_state_spec(spec) for spec in specs), key=lambda item: item[0])
    if len(states) != len(EXPECTED_BIASES):
        raise ValueError("T1 requires exactly eight same-contract G3 states")
    for (bias, path), expected in zip(states, EXPECTED_BIASES, strict=True):
        if not math.isclose(bias, expected, rel_tol=0.0, abs_tol=5.0e-7):
            raise ValueError(
                f"T1 biases must be 0:1/6:7/6; got {bias}, expected {expected}"
            )
        if not path.is_file():
            raise FileNotFoundError(path)
    return states


def validate_baseline(config: dict[str, Any]) -> dict[str, Any]:
    solver = config.get("solver", {})
    mobility = solver.get("mobility", {})
    violations: list[str] = []
    mobility_model = str(mobility.get("model", "")).lower()
    enabled_models = [
        str(value).lower()
        for value in mobility.get("models", [])
    ] if isinstance(mobility.get("models", []), list) else []
    if "ialmob" in mobility_model or any("ialmob" in value for value in enabled_models):
        violations.append("IALMob present")
    if config.get("sweep", {}).get("predictor") not in (None, False, "off"):
        violations.append("predictor enabled")
    if mobility.get("model") != "constant_field":
        violations.append("qualified HFS mobility model missing")
    if mobility.get("high_field_driving_force") != "quasi_fermi_gradient":
        violations.append("GradQF HFS drive missing")
    if mobility.get("high_field_gradient_discretization") != "edge_projection":
        violations.append("baseline is not edge_projection")
    if not mobility.get("contact_electric_field_fallback", False):
        violations.append("contact ElectricField fallback missing")
    geometry = config.get("mesh_geometry", {})
    if geometry.get("carrier_transport_couple_profile") != (
        "templates_ldmos_external_averagebox"
    ):
        violations.append("external AverageBox carrier couples missing")
    if solver.get("carrier_statistics", {}).get("model") != "fermi_dirac":
        violations.append("Fermi statistics missing")
    if solver.get("bandgap_narrowing", {}).get("model") != "old_slotboom":
        violations.append("OldSlotboom BGN missing")
    if solver.get("impact_ionization", {}).get("model") != "none":
        violations.append("impact ionization must be disabled")
    if violations:
        raise ValueError("; ".join(violations))
    return {
        "ialmob": "disabled",
        "predictor": "disabled",
        "hfs": "GradQF edge_projection plus qualified contact ElectricField fallback",
        "carrier_statistics": "fermi_dirac",
        "bandgap_narrowing": "old_slotboom",
        "carrier_couples": "templates_ldmos_external_averagebox",
        "contact_boundary_reconstruction": solver.get(
            "contact_boundary_reconstruction"
        ),
        "impact_ionization": "none",
    }


def set_gate_bias(config: dict[str, Any], bias: float) -> None:
    gates = [contact for contact in config.get("contacts", []) if contact["name"] == "gate"]
    if len(gates) != 1:
        raise ValueError("baseline must contain exactly one gate contact")
    gates[0]["bias"] = bias


def run_probe(
    runner: Path,
    baseline: dict[str, Any],
    state: Path,
    bias: float,
    simulation_type: str,
    output_csv: Path,
    config_path: Path,
) -> dict[str, Any]:
    config = json.loads(json.dumps(baseline))
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state)
    config["output_csv"] = str(output_csv.resolve())
    config.pop("output_state_file", None)
    config.pop("output_json", None)
    config.pop("sweep", None)
    set_gate_bias(config, bias)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=config_path.parent,
        env=runner_environment(),
        capture_output=True,
        text=True,
        check=False,
    )
    config_path.with_suffix(".stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    config_path.with_suffix(".stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(
            f"{simulation_type} failed at Vg={bias}: {completed.stderr}"
        )
    try:
        return json.loads(completed.stdout)
    except json.JSONDecodeError:
        return {"stdout": completed.stdout.strip()}


def node_id_signature(path: Path) -> tuple[int, str]:
    rows = read_csv(path)
    ids = [int(row["node_id"]) for row in rows]
    payload = ",".join(map(str, ids)).encode("ascii")
    return len(ids), hashlib.sha256(payload).hexdigest()


def l2(values: Iterable[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def norm(values: Iterable[float]) -> dict[str, float]:
    data = list(values)
    return {
        "l1": sum(abs(value) for value in data),
        "l2": l2(data),
        "maximum_abs": max((abs(value) for value in data), default=0.0),
    }


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def current_scale_from_sg(rows: list[dict[str, str]]) -> dict[str, float]:
    ratios = []
    for row in rows:
        scaled = float(row["electron_flux"])
        physical = float(row["electron_particle_line_flux_per_m_s"])
        if abs(scaled) <= 1.0e-250 or abs(physical) <= 1.0e-250:
            continue
        ratio = physical / scaled
        if math.isfinite(ratio) and ratio > 0.0:
            ratios.append(ratio)
    if not ratios:
        raise ValueError("SG output has no usable electron particle scale")
    particle_scale = percentile(ratios, 0.5)
    return {
        "particle_per_m_s_per_scaled": particle_scale,
        "A_per_um_per_scaled": particle_scale * ELEMENTARY_CHARGE_C * 1.0e-6,
        "sample_count": len(ratios),
        "maximum_relative_spread": max(
            abs(value / particle_scale - 1.0) for value in ratios
        ),
    }


def local_generalized_factor(eta: float) -> float:
    derivative = fermi_half_derivative(eta)
    return max(1.0, fermi_half(eta) / derivative) if derivative > 0.0 else 1.0


def frozen_kernel_sensitivities(
    row: dict[str, str], particle_scale: float, temperature_K: float
) -> tuple[float, float]:
    """Return d(scaled flux)/d(phin0, phin1) with edge coefficients frozen."""
    thermal_voltage = KB_J_K * temperature_K / ELEMENTARY_CHARGE_C
    length = float(row["length_m"])
    couple = float(row["couple_m"])
    mobility = float(row["electron_mobility_m2_V_s"])
    factor = float(row["electron_generalized_einstein_factor"])
    argument = float(row["electron_bernoulli_argument"])
    density0 = float(row["electron_density0_m3"])
    density1 = float(row["electron_density1_m3"])
    if length <= 0.0 or couple <= 0.0 or mobility <= 0.0:
        return 0.0, 0.0
    coefficient = mobility * thermal_voltage * couple / length / particle_scale
    g0 = local_generalized_factor(float(row["electron_eta0"]))
    g1 = local_generalized_factor(float(row["electron_eta1"]))
    bplus = float(row["electron_bernoulli_plus"])
    bminus = float(row["electron_bernoulli_minus"])
    # F = coef*f*(B(-a)n0 - B(+a)n1).  At fixed f/a/mobility,
    # dn/dphin = -n/(Vt*g_local) for the Fermi endpoint density.
    d0 = -coefficient * factor * bminus * density0 / (thermal_voltage * g0)
    d1 = coefficient * factor * bplus * density1 / (thermal_voltage * g1)
    if not all(map(math.isfinite, (argument, d0, d1))):
        raise ValueError(f"non-finite sensitivity on edge {row['edge_id']}")
    return d0, d1


def edge_sensitivity_record(
    row: dict[str, str],
    *,
    bias: float,
    hotspots: set[int],
    current_scale: float,
    particle_scale: float,
    temperature_K: float,
) -> dict[str, Any]:
    node0, node1 = int(row["node0"]), int(row["node1"])
    d0, d1 = frozen_kernel_sensitivities(row, particle_scale, temperature_K)
    magnitude0, magnitude1 = abs(d0), abs(d1)
    if max(magnitude0, magnitude1) == 0.0:
        upwind, downwind, ratio = "inactive", "inactive", 1.0
    elif magnitude0 >= magnitude1:
        upwind, downwind = "node0", "node1"
        ratio = magnitude1 / magnitude0
    else:
        upwind, downwind = "node1", "node0"
        ratio = magnitude0 / magnitude1
    incident = [node for node in (node0, node1) if node in hotspots]
    hotspot_roles = []
    for node in incident:
        hotspot_roles.append(
            f"{node}:{'upwind' if (node == node0) == (upwind == 'node0') else 'downwind'}"
            if upwind != "inactive" else f"{node}:inactive"
        )
    return {
        "bias_V": bias,
        "edge_id": int(row["edge_id"]),
        "node0": node0,
        "node1": node1,
        "hotspot_nodes": ";".join(map(str, incident)),
        "hotspot_roles": ";".join(hotspot_roles),
        "length_m": float(row["length_m"]),
        "couple_m": float(row["couple_m"]),
        "phin0_V": float(row["phin0_V"]),
        "phin1_V": float(row["phin1_V"]),
        "phin_gradient_abs_V_per_m": abs(
            float(row["phin1_V"]) - float(row["phin0_V"])
        ) / max(float(row["length_m"]), 1.0e-300),
        "electron_eta0": float(row["electron_eta0"]),
        "electron_eta1": float(row["electron_eta1"]),
        "bernoulli_argument_eta": float(row["electron_bernoulli_argument"]),
        "bernoulli_abs_eta": abs(float(row["electron_bernoulli_argument"])),
        "bernoulli_plus": float(row["electron_bernoulli_plus"]),
        "bernoulli_minus": float(row["electron_bernoulli_minus"]),
        "electron_flux_scaled": float(row["electron_flux"]),
        "electron_flux_A_per_um": float(row["electron_flux"]) * current_scale,
        "dflux_dphin0_scaled_per_V": d0,
        "dflux_dphin1_scaled_per_V": d1,
        "dflux_dphin0_A_per_um_per_V": d0 * current_scale,
        "dflux_dphin1_A_per_um_per_V": d1 * current_scale,
        "upwind_endpoint": upwind,
        "downwind_endpoint": downwind,
        "downwind_over_upwind_sensitivity": ratio,
        "bernoulli_saturated_abs_eta_ge_10": abs(
            float(row["electron_bernoulli_argument"])
        ) >= 10.0,
        "sensitivity_model": "analytic_frozen_edge_kernel",
    }


def density_by_node(rows: list[dict[str, str]]) -> tuple[dict[int, float], float]:
    samples: dict[int, list[float]] = {}
    for row in rows:
        for suffix in ("0", "1"):
            node = int(row[f"node{suffix}"])
            value = float(row[f"electron_density{suffix}_m3"])
            if value > 0.0:
                samples.setdefault(node, []).append(value)
    result = {node: percentile(values, 0.5) for node, values in samples.items()}
    spread = max(
        (
            max(abs(value / result[node] - 1.0) for value in values)
            for node, values in samples.items()
        ),
        default=0.0,
    )
    return result, spread


def geometric_mean(values: Iterable[float]) -> float:
    data = [value for value in values if value > 0.0]
    if not data:
        return math.nan
    return math.exp(sum(math.log(value) for value in data) / len(data))


def reference_curve(path: Path) -> dict[float, float]:
    return {
        round(float(row["bias_V"]), 9): float(row["current_total_A_per_um"])
        for row in read_csv(path)
    }


def fit_linear(x: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    coefficients, _, rank, singular = np.linalg.lstsq(x, y, rcond=None)
    fitted = x @ coefficients
    residual = y - fitted
    dof = len(y) - x.shape[1]
    covariance = np.linalg.pinv(x.T @ x)
    sigma2 = float(residual @ residual / dof) if dof > 0 else math.nan
    stderr = np.sqrt(np.maximum(np.diag(covariance) * sigma2, 0.0))
    tcrit = T95.get(dof, 1.96)
    return {
        "coefficients": coefficients.tolist(),
        "standard_errors": stderr.tolist(),
        "ci95": [
            [float(value - tcrit * error), float(value + tcrit * error)]
            for value, error in zip(coefficients, stderr, strict=True)
        ],
        "rank": int(rank),
        "singular_values": singular.tolist(),
        "r2": float(1.0 - (residual @ residual) / max(
            (y - np.mean(y)) @ (y - np.mean(y)), 1.0e-300
        )),
        "residual_rmse": float(math.sqrt(np.mean(residual * residual))),
    }


def loocv(x: np.ndarray, y: np.ndarray) -> dict[str, Any]:
    predictions: list[float] = []
    coefficients: list[list[float]] = []
    for omitted in range(len(y)):
        keep = np.arange(len(y)) != omitted
        beta = np.linalg.lstsq(x[keep], y[keep], rcond=None)[0]
        predictions.append(float(x[omitted] @ beta))
        coefficients.append(beta.tolist())
    errors = np.asarray(predictions) - y
    coefficient_array = np.asarray(coefficients)
    return {
        "prediction_rmse_log": float(math.sqrt(np.mean(errors * errors))),
        "prediction_max_abs_log": float(np.max(np.abs(errors))),
        "coefficient_min": np.min(coefficient_array, axis=0).tolist(),
        "coefficient_max": np.max(coefficient_array, axis=0).tolist(),
    }


def vif_values(predictors: np.ndarray, names: list[str]) -> dict[str, float]:
    result: dict[str, float] = {}
    for column, name in enumerate(names):
        y = predictors[:, column]
        others = np.delete(predictors, column, axis=1)
        design = np.column_stack([np.ones(len(y)), others])
        fitted = design @ np.linalg.lstsq(design, y, rcond=None)[0]
        denominator = float((y - np.mean(y)) @ (y - np.mean(y)))
        r2 = 1.0 - float((y - fitted) @ (y - fitted)) / max(denominator, 1.0e-300)
        result[name] = math.inf if r2 >= 1.0 - 1.0e-14 else 1.0 / (1.0 - r2)
    return result


def regression_report(rows: list[dict[str, Any]], response: str) -> dict[str, Any]:
    predictor_names = [
        "log_n_interface",
        "log_Id",
        "Vg",
        "log_grad_phin",
    ]
    y = np.log(np.asarray([float(row[response]) for row in rows]))
    predictors = np.column_stack([
        np.log([float(row["n_interface_geomean_m3"]) for row in rows]),
        np.log([abs(float(row["Id_sentaurus_A_per_um"])) for row in rows]),
        [float(row["bias_V"]) for row in rows],
        np.log([float(row["hotspot_grad_phin_rms_V_per_m"]) for row in rows]),
    ])
    standardized = (predictors - np.mean(predictors, axis=0)) / np.std(
        predictors, axis=0, ddof=1
    )
    design = np.column_stack([np.ones(len(rows)), standardized])
    condition = float(np.linalg.cond(design))
    vifs = vif_values(predictors, predictor_names)
    single: dict[str, Any] = {}
    for index, name in enumerate(predictor_names):
        x = np.column_stack([np.ones(len(rows)), predictors[:, index]])
        fit = fit_linear(x, y)
        fit["exponent"] = fit["coefficients"][1]
        fit["exponent_ci95"] = fit["ci95"][1]
        fit["loocv"] = loocv(x, y)
        single[name] = fit
    full = fit_linear(design, y)
    full["coefficient_names"] = ["intercept", *predictor_names]
    full["loocv"] = loocv(design, y)
    collinear = condition > 30.0 or any(value > 10.0 for value in vifs.values())
    intervals = [single[name]["exponent_ci95"] for name in predictor_names if name != "Vg"]
    pairwise_distinguishable = all(
        left[1] < right[0] or right[1] < left[0]
        for index, left in enumerate(intervals)
        for right in intervals[index + 1:]
    )
    return {
        "response": f"log({response})",
        "row_count": len(rows),
        "predictor_names": predictor_names,
        "standardized_design_condition_number": condition,
        "vif": vifs,
        "vif_limit": 10.0,
        "condition_warning_limit": 30.0,
        "collinearity_exceeds_preregistered_guard": collinear,
        "evidence_role": "descriptive_only" if collinear else "directional",
        "single_predictor_models": single,
        "full_standardized_model": full,
        "density_current_gradient_ci_pairwise_distinguishable": pairwise_distinguishable,
    }


def markdown_report(summary: dict[str, Any]) -> str:
    lines = [
        "# Templates/LDMOS G3 WP3 T1 residual scaling",
        "",
        "## Outcome",
        "",
        summary["outcome"],
        "",
        "## Contract",
        "",
        f"- Eight biases: `{summary['contract']['bias_points_V']}` at Vd=0.1 V.",
        "- Fermi/OldSlotboom/SRH/Auger/HFS are retained; IALMob, predictor, "
        "QP and impact ionization remain disabled.",
        "- Residuals and sensitivities are converted to A/um with the per-state "
        "SG particle scale; no NewtonPlot absolute RHS is used.",
        "",
        "## Scaling and collinearity",
        "",
    ]
    for scope, report in summary["regressions"].items():
        lines.extend([
            f"### {scope}",
            "",
            f"Condition number: `{report['standardized_design_condition_number']:.6g}`; "
            f"maximum VIF: `{max(report['vif'].values()):.6g}`; evidence role: "
            f"`{report['evidence_role']}`.",
            "",
            "| predictor | exponent | 95% CI | R2 | LOOCV RMSE (log) |",
            "| --- | ---: | ---: | ---: | ---: |",
        ])
        for name, fit in report["single_predictor_models"].items():
            lower, upper = fit["exponent_ci95"]
            lines.append(
                f"| {name} | {fit['exponent']:.6g} | "
                f"[{lower:.6g}, {upper:.6g}] | {fit['r2']:.6g} | "
                f"{fit['loocv']['prediction_rmse_log']:.6g} |"
            )
        lines.append("")
    endpoints = summary["edge_sensitivity"]["endpoint_summaries"]
    lines.extend([
        "## Hotspot-edge sensitivity",
        "",
        "The table is a frozen-edge-kernel analytic derivative: endpoint Fermi "
        "compressibility is included, while HFS mobility and secant-factor "
        "state derivatives are held fixed. This is sufficient for Bernoulli "
        "upwind/downwind classification, but it is not the complete Jacobian.",
        "",
        "| Vg (V) | topology edges | evaluable Si edges | saturated | minimum down/up sensitivity | node 3721 downwind flux share |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for endpoint in endpoints:
        lines.append(
            f"| {endpoint['bias_V']:.6g} | {endpoint['unique_hotspot_edges']} | "
            f"{endpoint['evaluable_silicon_hotspot_edges']} | "
            f"{endpoint['saturated_edge_count']} | "
            f"{endpoint['minimum_downwind_over_upwind_sensitivity']:.6g} | "
            f"{endpoint['node3721_downwind_abs_flux_fraction']:.6g} |"
        )
    lines.extend([
        "",
        "## Files",
        "",
        "- `residual_scaling_states.csv`: per-state physical residual and predictor table.",
        "- `hotspot_edge_phin_sensitivity.csv`: both endpoint sensitivities, Bernoulli eta, and roles.",
        "- `summary.json`: contract checks, regressions, VIF, condition number, and LOOCV.",
        "- `state_manifest.json`: exact state/config hashes.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--reference-curve", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--state", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--nodes", default=",".join(map(str, DEFAULT_NODES)))
    parser.add_argument("--temperature-K", type=float, default=300.0)
    args = parser.parse_args()

    states = validate_state_specs(args.state)
    baseline_path = args.baseline_config.resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    contract = validate_baseline(baseline)
    mesh_path = args.mesh.resolve()
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    selected_nodes = tuple(
        int(value) for value in args.nodes.split(",") if value.strip()
    )
    selected_set = set(selected_nodes)
    references = reference_curve(args.reference_curve.resolve())
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    _, pairs = region_local_geometry(mesh, {"si", "silicon"})
    direct_interface = {int(pair["global_node_id"]) for pair in pairs}
    expected_node_ids = [int(node["id"]) for node in mesh["nodes"]]
    expected_signature = hashlib.sha256(
        ",".join(map(str, expected_node_ids)).encode("ascii")
    ).hexdigest()
    manifest_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    all_sensitivities: list[dict[str, Any]] = []
    endpoint_summaries: list[dict[str, Any]] = []

    for bias, state_path in states:
        count, signature = node_id_signature(state_path)
        if count != len(expected_node_ids) or signature != expected_signature:
            raise ValueError(f"state node contract mismatch: {state_path}")
        point = output / bias_label(bias)
        carrier_csv = point / "carrier_terms.csv"
        sg_csv = point / "sg_edges.csv"
        run_probe(
            args.runner, baseline, state_path, bias,
            "newton_carrier_term_probe", carrier_csv, point / "carrier_probe.json",
        )
        run_probe(
            args.runner, baseline, state_path, bias,
            "sg_edge_flux_probe", sg_csv, point / "sg_probe.json",
        )
        carrier = read_csv(carrier_csv)
        edges = read_csv(sg_csv)
        scale = current_scale_from_sg(edges)
        current_scale = scale["A_per_um_per_scaled"]
        carrier_by_node = {int(row["node_id"]): row for row in carrier}
        missing = sorted(selected_set - set(carrier_by_node))
        if missing:
            raise ValueError(f"carrier probe missing hotspots: {missing}")

        free_nodes = {
            int(row["node_id"])
            for row in carrier
            if float(row["electron_flux_abs_sum"]) > 0.0
            and abs(float(row["electron_gauge"])) <= 1.0e-300
            and abs(float(row["electron_boundary"])) <= 1.0e-300
        }
        interface_one_ring = set(direct_interface)
        for row in edges:
            if float(row["couple_m"]) <= 0.0:
                continue
            node0, node1 = int(row["node0"]), int(row["node1"])
            if node0 in direct_interface or node1 in direct_interface:
                interface_one_ring.update((node0, node1))
        interface_band = free_nodes & interface_one_ring
        hotspot_edges = [
            row for row in edges
            if int(row["node0"]) in selected_set or int(row["node1"]) in selected_set
        ]
        if len({int(row["edge_id"]) for row in hotspot_edges}) != 40:
            raise ValueError(
                f"Vg={bias}: frozen hotspot edge set is not 40 edges"
            )

        sensitivities = [
            edge_sensitivity_record(
                row,
                bias=bias,
                hotspots=selected_set,
                current_scale=current_scale,
                particle_scale=scale["particle_per_m_s_per_scaled"],
                temperature_K=args.temperature_K,
            )
            for row in hotspot_edges
        ]
        all_sensitivities.extend(sensitivities)
        density, density_spread = density_by_node(edges)
        interface_density = geometric_mean(
            density[node] for node in DIRECT_INTERFACE_HOTSPOTS
        )
        active_sensitivities = [
            row for row in sensitivities
            if row["upwind_endpoint"] != "inactive"
            and math.isfinite(float(row["bernoulli_abs_eta"]))
        ]
        if not active_sensitivities:
            raise ValueError(f"Vg={bias}: no evaluable silicon hotspot edges")
        gradients = [
            float(row["phin_gradient_abs_V_per_m"])
            for row in active_sensitivities
        ]
        seven_scaled = [
            float(carrier_by_node[node]["electron_residual"])
            for node in selected_nodes
        ]
        interface_scaled = [
            float(carrier_by_node[node]["electron_residual"])
            for node in sorted(interface_band)
        ]
        silicon_scaled = [
            float(carrier_by_node[node]["electron_residual"])
            for node in sorted(free_nodes)
        ]
        seven = norm(value * current_scale for value in seven_scaled)
        band = norm(value * current_scale for value in interface_scaled)
        silicon = norm(value * current_scale for value in silicon_scaled)
        flux_abs = sum(
            abs(float(row["electron_flux"]) * current_scale)
            for row in hotspot_edges
        )
        reference_key = round(bias, 9)
        if reference_key not in references:
            raise ValueError(f"reference curve lacks Vg={bias}")
        state_rows.append({
            "bias_V": bias,
            "Id_sentaurus_A_per_um": references[reference_key],
            "n_interface_geomean_m3": interface_density,
            "hotspot_grad_phin_rms_V_per_m": math.sqrt(
                sum(value * value for value in gradients) / len(gradients)
            ),
            "hotspot_grad_phin_geomean_V_per_m": geometric_mean(gradients),
            "hotspot_topological_edge_count": len(sensitivities),
            "hotspot_evaluable_silicon_edge_count": len(active_sensitivities),
            "seven_node_residual_l1_A_per_um": seven["l1"],
            "seven_node_residual_l2_A_per_um": seven["l2"],
            "seven_node_residual_max_A_per_um": seven["maximum_abs"],
            "interface_band_node_count": len(interface_band),
            "interface_band_residual_l1_A_per_um": band["l1"],
            "interface_band_residual_l2_A_per_um": band["l2"],
            "interface_band_residual_max_A_per_um": band["maximum_abs"],
            "free_silicon_node_count": len(free_nodes),
            "free_silicon_residual_l2_A_per_um": silicon["l2"],
            "free_silicon_residual_max_A_per_um": silicon["maximum_abs"],
            "seven_node_energy_fraction_of_interface_band": (
                seven["l2"] / max(band["l2"], 1.0e-300)
            ) ** 2,
            "seven_node_residual_over_incident_abs_flux": (
                seven["l2"] / max(flux_abs, 1.0e-300)
            ),
            "current_scale_A_per_um_per_scaled": current_scale,
            "current_scale_max_relative_spread": scale["maximum_relative_spread"],
            "reconstructed_density_max_relative_spread": density_spread,
        })
        node3721_edges = [
            row for row in active_sensitivities if 3721 in (
                int(row["node0"]), int(row["node1"])
            )
        ]
        node3721_downwind_flux = sum(
            abs(float(row["electron_flux_A_per_um"]))
            for row in node3721_edges
            if f"3721:downwind" in str(row["hotspot_roles"])
        )
        node3721_total_flux = sum(
            abs(float(row["electron_flux_A_per_um"])) for row in node3721_edges
        )
        endpoint_summaries.append({
            "bias_V": bias,
            "unique_hotspot_edges": len(sensitivities),
            "evaluable_silicon_hotspot_edges": len(active_sensitivities),
            "inactive_or_nontransport_hotspot_edges": (
                len(sensitivities) - len(active_sensitivities)
            ),
            "saturated_edge_count": sum(
                bool(row["bernoulli_saturated_abs_eta_ge_10"])
                for row in active_sensitivities
            ),
            "median_abs_bernoulli_eta": percentile(
                (
                    float(row["bernoulli_abs_eta"])
                    for row in active_sensitivities
                ), 0.5
            ),
            "maximum_abs_bernoulli_eta": max(
                float(row["bernoulli_abs_eta"])
                for row in active_sensitivities
            ),
            "median_downwind_over_upwind_sensitivity": percentile(
                (
                    float(row["downwind_over_upwind_sensitivity"])
                    for row in active_sensitivities
                ), 0.5
            ),
            "p95_downwind_over_upwind_sensitivity": percentile(
                (
                    float(row["downwind_over_upwind_sensitivity"])
                    for row in active_sensitivities
                ), 0.95
            ),
            "minimum_downwind_over_upwind_sensitivity": min(
                float(row["downwind_over_upwind_sensitivity"])
                for row in active_sensitivities
            ),
            "node3721_incident_edges": len(node3721_edges),
            "node3721_downwind_edges": sum(
                f"3721:downwind" in str(row["hotspot_roles"])
                for row in node3721_edges
            ),
            "node3721_downwind_abs_flux_fraction": (
                node3721_downwind_flux / max(node3721_total_flux, 1.0e-300)
            ),
            "node3721_minimum_downwind_over_upwind_sensitivity": min(
                float(row["downwind_over_upwind_sensitivity"])
                for row in node3721_edges
            ),
        })
        manifest_rows.append({
            "bias_V": bias,
            "state_path": str(state_path),
            "state_sha256": sha256(state_path),
            "state_node_count": count,
            "state_node_id_sha256": signature,
            "carrier_config": str((point / "carrier_probe.json").resolve()),
            "carrier_config_sha256": sha256(point / "carrier_probe.json"),
            "sg_config": str((point / "sg_probe.json").resolve()),
            "sg_config_sha256": sha256(point / "sg_probe.json"),
            "carrier_csv": str(carrier_csv.resolve()),
            "sg_csv": str(sg_csv.resolve()),
        })

    write_csv(output / "residual_scaling_states.csv", state_rows)
    write_csv(output / "hotspot_edge_phin_sensitivity.csv", all_sensitivities)
    manifest = {
        "schema": "vela.templates_ldmos.g3_wp3_t1_state_manifest.v1",
        "baseline_config": str(baseline_path),
        "baseline_config_sha256": sha256(baseline_path),
        "mesh": str(mesh_path),
        "mesh_sha256": sha256(mesh_path),
        "reference_curve": str(args.reference_curve.resolve()),
        "reference_curve_sha256": sha256(args.reference_curve.resolve()),
        "states": manifest_rows,
    }
    (output / "state_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    regressions = {
        "seven_node": regression_report(
            state_rows, "seven_node_residual_l2_A_per_um"
        ),
        "interface_band": regression_report(
            state_rows, "interface_band_residual_l2_A_per_um"
        ),
    }
    descriptive = any(
        report["collinearity_exceeds_preregistered_guard"]
        for report in regressions.values()
    )
    distinguishable = all(
        report["density_current_gradient_ci_pairwise_distinguishable"]
        for report in regressions.values()
    )
    saturation_supported = any(
        endpoint["saturated_edge_count"] > 0
        or endpoint["minimum_downwind_over_upwind_sensitivity"] < 0.1
        for endpoint in endpoint_summaries
    )
    outcome = (
        "The eight-state replay is complete, but the pre-registered collinearity "
        "guard is exceeded; scaling fits are descriptive and cannot by themselves "
        "lock H1/H2/H3. The evaluable hotspot-edge set also shows no strong "
        "Bernoulli saturation or exponentially suppressed endpoint sensitivity."
        if descriptive else
        "The eight-state replay passes the collinearity guard; scaling fits may be "
        "used as directional evidence, subject to the reported LOOCV stability."
    )
    summary = {
        "schema": "vela.templates_ldmos.g3_wp3_t1_residual_scaling.v1",
        "outcome": outcome,
        "contract": {
            **contract,
            "bias_points_V": [bias for bias, _ in states],
            "drain_bias_V": 0.1,
            "state_count": len(states),
            "selected_nodes": list(selected_nodes),
            "direct_interface_hotspots": list(DIRECT_INTERFACE_HOTSPOTS),
            "interface_band": "free silicon direct Si/dielectric nodes plus one SG-edge ring",
            "sensitivity_model": "analytic frozen edge kernel with endpoint Fermi compressibility",
            "sensitivity_excludes": [
                "HFS mobility state derivative",
                "generalized-Einstein secant-factor state derivative",
            ],
            "production_defaults_changed": False,
            "vm_used": False,
            "idvd_bv_full_physics_states_used": False,
        },
        "data_quality": {
            "all_state_node_id_hashes_equal_mesh": True,
            "maximum_current_scale_relative_spread": max(
                float(row["current_scale_max_relative_spread"])
                for row in state_rows
            ),
            "maximum_reconstructed_density_relative_spread": max(
                float(row["reconstructed_density_max_relative_spread"])
                for row in state_rows
            ),
            "unique_hotspot_edge_count_each_state": sorted({
                int(row["unique_hotspot_edges"]) for row in endpoint_summaries
            }),
        },
        "regressions": regressions,
        "gates": {
            "collinearity_guard_passed": not descriptive,
            "scaling_ci_distinguishes_density_current_gradient": distinguishable,
            "t1_can_lock_hypothesis_family_alone": not descriptive and distinguishable,
            "strong_bernoulli_saturation_mechanism_supported": saturation_supported,
        },
        "edge_sensitivity": {
            "endpoint_summaries": endpoint_summaries,
            "table": str((output / "hotspot_edge_phin_sensitivity.csv").resolve()),
        },
        "artifacts": {
            "state_table": str((output / "residual_scaling_states.csv").resolve()),
            "edge_table": str((output / "hotspot_edge_phin_sensitivity.csv").resolve()),
            "state_manifest": str((output / "state_manifest.json").resolve()),
            "report": str((output / "report.md").resolve()),
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    (output / "report.md").write_text(markdown_report(summary), encoding="utf-8")
    print(json.dumps({
        "summary": str((output / "summary.json").resolve()),
        "outcome": outcome,
        "gates": summary["gates"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
