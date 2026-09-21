#!/usr/bin/env python3
"""Audit Sentaurus NewtonPlot continuity-RHS AreaFactor scaling."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any


def read_field(export_root: Path, name: str = "eContinuityRhs") -> dict[int, float]:
    path = export_root / "fields" / f"{name}_region0.csv"
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(stream)
        }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def vector_scale(target: list[float], source: list[float]) -> dict[str, float]:
    denominator = sum(value * value for value in source)
    if denominator == 0.0:
        raise ValueError("cannot fit a zero vector")
    fitted = sum(a * b for a, b in zip(target, source, strict=True)) / denominator
    source_norm = l2(source)
    target_norm = l2(target)
    residual = [a - fitted * b for a, b in zip(target, source, strict=True)]
    return {
        "least_squares_scale": fitted,
        "cosine": (
            sum(a * b for a, b in zip(target, source, strict=True))
            / (target_norm * source_norm)
        ),
        "relative_l2_after_scale": l2(residual) / max(target_norm, 1.0e-300),
    }


def central_difference(
    minus: dict[int, float], plus: dict[int, float], amplitude_V: float
) -> dict[int, float]:
    common = sorted(set(minus) & set(plus))
    return {
        node: (plus[node] - minus[node]) / (2.0 * amplitude_V)
        for node in common
    }


CONTACT_RE = re.compile(
    r"^\s*drain\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+"
    r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s*$",
    re.MULTILINE,
)


def drain_current(log_path: Path) -> dict[str, float]:
    matches = CONTACT_RE.findall(log_path.read_text(encoding="utf-8", errors="replace"))
    if not matches:
        raise ValueError(f"no final drain-current row in {log_path}")
    voltage, electron, hole, conduction = map(float, matches[-1])
    return {
        "voltage_V": voltage,
        "electron_current_A": electron,
        "hole_current_A": hole,
        "conduction_current_A": conduction,
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-minus-export", type=Path, required=True)
    parser.add_argument("--baseline-plus-export", type=Path, required=True)
    parser.add_argument("--candidate-minus-export", type=Path, required=True)
    parser.add_argument("--candidate-plus-export", type=Path, required=True)
    parser.add_argument("--baseline-minus-log", type=Path, required=True)
    parser.add_argument("--baseline-plus-log", type=Path, required=True)
    parser.add_argument("--candidate-minus-log", type=Path, required=True)
    parser.add_argument("--candidate-plus-log", type=Path, required=True)
    parser.add_argument("--single-mode-summary", type=Path, required=True)
    parser.add_argument("--baseline-area-factor", type=float, default=1.0)
    parser.add_argument("--candidate-area-factor", type=float, default=2.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    prior = json.loads(args.single_mode_summary.read_text(encoding="utf-8"))
    amplitude = float(prior["contract"]["amplitude_V"])
    ring = [int(node) for node in prior["pre_registered_support"]["nodes"]]

    baseline_minus = read_field(args.baseline_minus_export)
    baseline_plus = read_field(args.baseline_plus_export)
    candidate_minus = read_field(args.candidate_minus_export)
    candidate_plus = read_field(args.candidate_plus_export)
    common = sorted(
        set(baseline_minus) & set(baseline_plus)
        & set(candidate_minus) & set(candidate_plus)
    )
    if sorted(set(ring) - set(common)):
        raise ValueError("pre-registered one-ring is incomplete")

    baseline_derivative = central_difference(baseline_minus, baseline_plus, amplitude)
    candidate_derivative = central_difference(candidate_minus, candidate_plus, amplitude)
    maximum = max(abs(value) for value in baseline_derivative.values())
    active = [
        node for node in common
        if abs(baseline_derivative[node]) >= 1.0e-8 * maximum
    ]
    scaling = vector_scale(
        [candidate_derivative[node] for node in active],
        [baseline_derivative[node] for node in active],
    )

    baseline_logs = {
        "minus": drain_current(args.baseline_minus_log),
        "plus": drain_current(args.baseline_plus_log),
    }
    candidate_logs = {
        "minus": drain_current(args.candidate_minus_log),
        "plus": drain_current(args.candidate_plus_log),
    }
    terminal_ratios = {
        sign: {
            quantity: candidate_logs[sign][quantity] / baseline_logs[sign][quantity]
            for quantity in (
                "electron_current_A", "hole_current_A", "conduction_current_A"
            )
        }
        for sign in ("minus", "plus")
    }

    field_paths = {
        "baseline_minus": args.baseline_minus_export / "fields" / "eContinuityRhs_region0.csv",
        "baseline_plus": args.baseline_plus_export / "fields" / "eContinuityRhs_region0.csv",
        "candidate_minus": args.candidate_minus_export / "fields" / "eContinuityRhs_region0.csv",
        "candidate_plus": args.candidate_plus_export / "fields" / "eContinuityRhs_region0.csv",
    }
    hashes = {name: sha256(path) for name, path in field_paths.items()}
    exact_pairwise_identity = {
        "minus": hashes["baseline_minus"] == hashes["candidate_minus"],
        "plus": hashes["baseline_plus"] == hashes["candidate_plus"],
    }
    rows = [{
        "node_id": node,
        "baseline_dR_dphin_A_per_V": baseline_derivative[node],
        "candidate_dR_dphin_A_per_V": candidate_derivative[node],
        "candidate_over_baseline": (
            candidate_derivative[node] / baseline_derivative[node]
        ),
    } for node in ring]

    expected_terminal_scale = args.candidate_area_factor / args.baseline_area_factor
    terminal_scales = [
        terminal_ratios[sign][quantity]
        for sign in terminal_ratios
        for quantity in terminal_ratios[sign]
    ]
    terminal_follows_area_factor = all(
        abs(value / expected_terminal_scale - 1.0) <= 5.0e-4
        for value in terminal_scales
    )
    rhs_is_area_factor_invariant = (
        all(exact_pairwise_identity.values())
        and abs(scaling["least_squares_scale"] - 1.0) <= 1.0e-12
        and scaling["relative_l2_after_scale"] <= 1.0e-12
    )
    report = {
        "experiment": "templates_ldmos_g3_rhs_dimension_chain_area_factor_ab",
        "contract": {
            "state": "same_high_endpoint_vsv_symmetric_node3721_qf_pair",
            "amplitude_V": amplitude,
            "baseline_area_factor": args.baseline_area_factor,
            "candidate_area_factor": args.candidate_area_factor,
            "only_changed_sentaurus_input": "global_Physics.AreaFactor",
            "production_defaults_changed": False,
        },
        "sentaurus_newtonplot_rhs": {
            "declared_tdr_unit": "A",
            "common_silicon_nodes": len(common),
            "active_derivative_nodes": len(active),
            "exact_csv_identity": exact_pairwise_identity,
            "candidate_over_baseline_derivative": scaling,
            "pre_registered_one_ring_nodes": ring,
        },
        "sentaurus_terminal_current": {
            "baseline": baseline_logs,
            "candidate": candidate_logs,
            "candidate_over_baseline": terminal_ratios,
            "expected_area_factor_scale": expected_terminal_scale,
        },
        "vela_dimension_chain": {
            "mesh_length_internal_unit": "um",
            "carrier_particle_line_flux_unit": "m^-1 s^-1",
            "terminal_current_unit": "A/um",
            "particle_line_flux_to_A_per_um_factor": 1.602176634e-25,
            "continuity_particle_scale_per_m_s": prior["scaling"][
                "continuity_particle_scale_per_m_s"
            ],
            "scaled_residual_to_A_per_um": prior["scaling"][
                "current_scale_A_per_um_per_scaled_residual"
            ],
        },
        "classification": {
            "terminal_current_follows_area_factor": terminal_follows_area_factor,
            "newtonplot_rhs_is_area_factor_invariant": rhs_is_area_factor_invariant,
            "two_dimensional_width_or_area_factor_explains_33x": False,
            "newtonplot_rhs_must_be_treated_as_implementation_dependent_internal_quantity": True,
            "remaining_uniform_assembly_or_internal_normalization_candidate": True,
            "ledger_status": "draft",
        },
        "artifacts": {
            str(path): sha256(path)
            for path in [*field_paths.values(), args.baseline_minus_log,
                         args.baseline_plus_log, args.candidate_minus_log,
                         args.candidate_plus_log, args.single_mode_summary]
        },
    }
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / "one_ring_rows.csv", rows)
    (output_dir / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
