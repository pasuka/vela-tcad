#!/usr/bin/env python3
"""Compare Sentaurus and Vela central-difference electron-QF row responses."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


ELEMENTARY_CHARGE_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def scalar_field(export_root: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_root / "fields" / f"{name}_region0.csv")
    }


def state_rows(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            "psi": float(row["psi"]),
            "phin": float(row["phin"]),
            "phip": float(row["phip"]),
        }
        for row in read_csv(path)
    }


def carrier_residuals(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["electron_residual"])
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
    return {
        "continuity_particle_scale_per_m_s": continuity_scale,
        "current_scale_A_per_um_per_scaled_residual": (
            continuity_scale * ELEMENTARY_CHARGE_C * 1.0e-6
        ),
        "sample_count": len(ratios),
        "maximum_relative_spread": max(
            abs(value / continuity_scale - 1.0) for value in ratios
        ),
    }


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def cosine(left: list[float], right: list[float]) -> float | None:
    denominator = l2(left) * l2(right)
    if denominator == 0.0:
        return None
    return sum(a * b for a, b in zip(left, right, strict=True)) / denominator


def fitted_scale(target: list[float], source: list[float]) -> float:
    denominator = sum(value * value for value in source)
    if denominator == 0.0:
        raise ValueError("cannot fit scale to a zero vector")
    return sum(a * b for a, b in zip(target, source, strict=True)) / denominator


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def one_ring(mesh_path: Path, target: int) -> list[int]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    triangles = [
        triangle for triangle in mesh["triangles"]
        if int(triangle["region_id"]) == 0
        and target in [int(node) for node in triangle["node_ids"]]
    ]
    if not triangles:
        raise ValueError(f"target node {target} has no silicon triangle")
    return sorted({int(node) for triangle in triangles for node in triangle["node_ids"]})


def state_validation(
    loaded_export: Path,
    target_state: dict[int, dict[str, float]],
) -> dict[str, Any]:
    loaded = {
        "psi": scalar_field(loaded_export, "ElectrostaticPotential"),
        "phin": scalar_field(loaded_export, "eQuasiFermiPotential"),
        "phip": scalar_field(loaded_export, "hQuasiFermiPotential"),
    }
    common = sorted(set(target_state).intersection(*[set(field) for field in loaded.values()]))
    return {
        "common_silicon_nodes": len(common),
        "max_abs_error_V": {
            name: max(abs(values[node] - target_state[node][name]) for node in common)
            for name, values in loaded.items()
        },
    }


def loaded_pair_validation(
    minus_export: Path,
    plus_export: Path,
    target: int,
    amplitude: float,
) -> dict[str, Any]:
    minus = {
        "psi": scalar_field(minus_export, "ElectrostaticPotential"),
        "phin": scalar_field(minus_export, "eQuasiFermiPotential"),
        "phip": scalar_field(minus_export, "hQuasiFermiPotential"),
    }
    plus = {
        "psi": scalar_field(plus_export, "ElectrostaticPotential"),
        "phin": scalar_field(plus_export, "eQuasiFermiPotential"),
        "phip": scalar_field(plus_export, "hQuasiFermiPotential"),
    }
    common = sorted(set(minus["psi"]) & set(plus["psi"]))
    return {
        "target_phin_delta_V": plus["phin"][target] - minus["phin"][target],
        "target_phin_delta_error_V": (
            plus["phin"][target] - minus["phin"][target] - 2.0 * amplitude
        ),
        "maximum_abs_unintended_pair_delta_V": {
            "psi": max(abs(plus["psi"][node] - minus["psi"][node]) for node in common),
            "phip": max(abs(plus["phip"][node] - minus["phip"][node]) for node in common),
            "phin_excluding_target": max(
                abs(plus["phin"][node] - minus["phin"][node])
                for node in common if node != target
            ),
        },
    }
def vector_metrics(sent: list[float], vela: list[float]) -> dict[str, Any]:
    scale = fitted_scale(sent, vela)
    error = [a - scale * b for a, b in zip(sent, vela, strict=True)]
    sent_norm = l2(sent)
    return {
        "sentaurus_l2_A_per_V": sent_norm,
        "vela_l2_A_per_um_per_V": l2(vela),
        "cosine": cosine(sent, vela),
        "signed_sentaurus_over_vela_least_squares": scale,
        "absolute_sentaurus_over_vela_least_squares": abs(scale),
        "relative_l2_after_scalar_fit": l2(error) / max(sent_norm, 1.0e-300),
        "maximum_abs_after_scalar_fit_over_sentaurus_l2": (
            max(abs(value) for value in error) / max(sent_norm, 1.0e-300)
        ),
    }


def nonlinearity(
    minus: dict[int, float],
    baseline: dict[int, float],
    plus: dict[int, float],
    nodes: list[int],
) -> dict[str, float]:
    odd = [plus[node] - minus[node] for node in nodes]
    even = [plus[node] + minus[node] - 2.0 * baseline[node] for node in nodes]
    return {
        "odd_difference_l2": l2(odd),
        "even_second_difference_l2": l2(even),
        "even_over_odd_l2": l2(even) / max(l2(odd), 1.0e-300),
    }


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--baseline-sentaurus-export", type=Path, required=True)
    parser.add_argument("--equation-balance-summary", type=Path, required=True)
    parser.add_argument("--minus-loaded-export", type=Path, required=True)
    parser.add_argument("--plus-loaded-export", type=Path, required=True)
    parser.add_argument("--minus-newton-export", type=Path, required=True)
    parser.add_argument("--plus-newton-export", type=Path, required=True)
    parser.add_argument("--baseline-vela-carrier", type=Path, required=True)
    parser.add_argument("--minus-vela-carrier", type=Path, required=True)
    parser.add_argument("--plus-vela-carrier", type=Path, required=True)
    parser.add_argument("--vela-sg", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    target = int(manifest["contract"]["node_id"])
    amplitude = float(manifest["contract"]["amplitude_V"])
    if amplitude <= 0.0:
        raise ValueError("manifest amplitude must be positive")
    ring = one_ring(args.mesh, target)
    minus_state = state_rows(Path(manifest["variants"]["minus"]["state_file"]))
    plus_state = state_rows(Path(manifest["variants"]["plus"]["state_file"]))

    sent_minus = scalar_field(args.minus_newton_export, "eContinuityRhs")
    sent_plus = scalar_field(args.plus_newton_export, "eContinuityRhs")
    sent_base = scalar_field(args.baseline_sentaurus_export, "eContinuityRhs")
    vela_minus = carrier_residuals(args.minus_vela_carrier)
    vela_plus = carrier_residuals(args.plus_vela_carrier)
    vela_base = carrier_residuals(args.baseline_vela_carrier)
    scale = current_scale_from_sg(args.vela_sg)
    vela_physical_scale = scale["current_scale_A_per_um_per_scaled_residual"]

    common = sorted(
        set(sent_minus) & set(sent_plus) & set(sent_base)
        & set(vela_minus) & set(vela_plus) & set(vela_base)
    )
    missing_ring = sorted(set(ring) - set(common))
    if missing_ring:
        raise ValueError(f"one-ring nodes missing from engine outputs: {missing_ring}")

    sent_derivative = {
        node: (sent_plus[node] - sent_minus[node]) / (2.0 * amplitude)
        for node in common
    }
    vela_derivative = {
        node: (vela_plus[node] - vela_minus[node])
        * vela_physical_scale / (2.0 * amplitude)
        for node in common
    }
    maximum = max(max(abs(value) for value in sent_derivative.values()), 1.0e-300)
    active = sorted({
        node for node in common
        if abs(sent_derivative[node]) >= 1.0e-8 * maximum
        or abs(vela_derivative[node]) >= 1.0e-8
        * max(max(abs(value) for value in vela_derivative.values()), 1.0e-300)
    })

    def values(mapping: dict[int, float], nodes: list[int]) -> list[float]:
        return [mapping[node] for node in nodes]

    ring_metrics = vector_metrics(values(sent_derivative, ring), values(vela_derivative, ring))
    active_metrics = vector_metrics(
        values(sent_derivative, active), values(vela_derivative, active)
    )
    equation_balance = json.loads(args.equation_balance_summary.read_text(encoding="utf-8"))
    baseline_scalar = abs(float(
        equation_balance["vsv_endpoints"]["high"]["cross_engine_row_mode"]
        ["signed_sentaurus_over_vela_least_squares"]
    ))
    rows = [{
        "node_id": node,
        "in_pre_registered_one_ring": node in ring,
        "sentaurus_minus_rhs_A": sent_minus[node],
        "sentaurus_plus_rhs_A": sent_plus[node],
        "sentaurus_dR_dphin_A_per_V": sent_derivative[node],
        "vela_minus_residual_scaled": vela_minus[node],
        "vela_plus_residual_scaled": vela_plus[node],
        "vela_dR_dphin_A_per_um_per_V": vela_derivative[node],
        "sentaurus_after_ring_scalar_fit_A_per_V": (
            sent_derivative[node]
            - ring_metrics["signed_sentaurus_over_vela_least_squares"]
            * vela_derivative[node]
        ),
    } for node in ring]

    inputs = [
        args.manifest, args.mesh, args.equation_balance_summary,
        args.baseline_sentaurus_export / "fields" / "eContinuityRhs_region0.csv",
        args.minus_loaded_export / "fields" / "eQuasiFermiPotential_region0.csv",
        args.plus_loaded_export / "fields" / "eQuasiFermiPotential_region0.csv",
        args.minus_newton_export / "fields" / "eContinuityRhs_region0.csv",
        args.plus_newton_export / "fields" / "eContinuityRhs_region0.csv",
        args.baseline_vela_carrier, args.minus_vela_carrier, args.plus_vela_carrier,
        args.vela_sg,
    ]
    report = {
        "experiment": "templates_ldmos_g3_qf_single_mode_central_difference",
        "contract": manifest["contract"],
        "pre_registered_support": {
            "definition": "target_node_plus_all_nodes_sharing_a_silicon_triangle",
            "nodes": ring,
        },
        "loaded_state_validation": {
            "minus": state_validation(args.minus_loaded_export, minus_state),
            "plus": state_validation(args.plus_loaded_export, plus_state),
            "symmetric_pair": loaded_pair_validation(
                args.minus_loaded_export, args.plus_loaded_export, target, amplitude
            ),
        },
        "cross_engine_central_difference": {
            "one_ring": ring_metrics,
            "active_support": {"nodes": active, **active_metrics},
            "calibration_against_unperturbed_row_scale": {
                "unperturbed_high_endpoint_absolute_scale": baseline_scalar,
                "single_mode_absolute_scale": ring_metrics[
                    "absolute_sentaurus_over_vela_least_squares"
                ],
                "relative_drift": abs(
                    ring_metrics["absolute_sentaurus_over_vela_least_squares"]
                    / baseline_scalar - 1.0
                ),
            },
        },
        "linearity": {
            "sentaurus_one_ring": nonlinearity(
                sent_minus, sent_base, sent_plus, ring
            ),
            "vela_one_ring_scaled": nonlinearity(
                vela_minus, vela_base, vela_plus, ring
            ),
        },
        "scaling": scale,
        "classification": {
            "stable_uniform_scalar_supported": (
                ring_metrics["cosine"] is not None
                and abs(ring_metrics["cosine"]) >= 0.999
                and ring_metrics["relative_l2_after_scalar_fit"] <= 0.05
            ),
            "rhs_or_width_normalization_remains_primary_candidate": (
                ring_metrics["cosine"] is not None
                and abs(ring_metrics["cosine"]) >= 0.999
                and 25.0 <= ring_metrics["absolute_sentaurus_over_vela_least_squares"] <= 45.0
            ),
            "qf_drive_spatial_shape_is_primary_cause": (
                ring_metrics["cosine"] is None
                or abs(ring_metrics["cosine"]) < 0.99
            ),
            "ledger_status": "draft",
        },
        "artifacts": {str(path): sha256(path) for path in inputs},
    }
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(output_dir / "one_ring_rows.csv", rows)
    (output_dir / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
