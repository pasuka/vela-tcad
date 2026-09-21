#!/usr/bin/env python3
"""Run a three-level fixed-state interface audit on selected LDMOS nodes.

Level A is the current diagnostic contract: one shared geometric node, direct
Silicon-only AverageBox transport couples, and Vela's global barycentric source
volume.  Level B keeps the transport rows fixed and changes only the local
continuity source volume to the Silicon-side AverageBox Measure.  Level C adds
an algebraic Silicon-master/oxide-slave potential pair with an exact continuity
constraint.  With the constraint eliminated, its carrier row must equal B;
the region-resolved Sentaurus potentials independently check that assumption.

This is deliberately a read-only frozen-state audit.  It neither adds solver
degrees of freedom nor changes a production default.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.audit_templates_ldmos_averagebox_full_mesh import averagebox_geometry
from scripts.audit_templates_ldmos_averagebox_node4492 import (
    norm,
    parse_debug_block,
    read_csv,
    write_csv,
)
from scripts.audit_templates_ldmos_interface_pair_box import region_local_geometry


DEFAULT_NODES = (3721, 3974, 3973, 4091, 3727, 4021, 3771)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_scalar_field(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            result[int(row["node_id"])] = float(row["component0"])
    return result


def _ratios(candidate: dict[str, float], baseline: dict[str, float]) -> dict[str, float]:
    return {
        key: candidate[key] / max(baseline[key], 1.0e-300)
        for key in baseline
    }


def audit_endpoint(
    *,
    bias_v: float,
    mesh: dict[str, Any],
    carrier_rows: list[dict[str, str]],
    measures: dict[int, list[float]],
    coefficients: dict[int, list[float]],
    transport_materials: set[str],
    selected_nodes: tuple[int, ...],
    silicon_potential: dict[int, float],
    oxide_potential: dict[int, float],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    _, silicon_measure, global_volume, _, geometry_metadata = averagebox_geometry(
        mesh, measures, coefficients, transport_materials
    )
    _, pairs = region_local_geometry(mesh, transport_materials)
    interface_nodes = {int(pair["global_node_id"]) for pair in pairs}
    carriers = {int(row["node_id"]): row for row in carrier_rows}
    missing = [node for node in selected_nodes if node not in carriers]
    if missing:
        raise ValueError(f"carrier probe lacks selected nodes: {missing}")

    rows: list[dict[str, Any]] = []
    level_a_values: list[float] = []
    level_b_values: list[float] = []
    level_c_values: list[float] = []
    constraint_values: list[float] = []
    for node in selected_nodes:
        raw = carriers[node]
        flux = float(raw["electron_flux"])
        recombination = float(raw["electron_recombination"])
        impact = float(raw["electron_impact"])
        gauge = float(raw["electron_gauge"])
        boundary = float(raw["electron_boundary"])
        source = recombination + impact
        level_a = float(raw["electron_residual"])
        volume_scale = silicon_measure.get(node, 0.0) / max(
            global_volume.get(node, 0.0), 1.0e-300
        )
        level_b = flux + volume_scale * source + gauge + boundary

        # Exact psi_slave = psi_master elimination leaves the master carrier row
        # unchanged.  This is the Schur-complement form of the constrained pair.
        level_c = level_b
        has_pair = node in interface_nodes
        si_psi = silicon_potential.get(node)
        ox_psi = oxide_potential.get(node)
        potential_split = (
            si_psi - ox_psi
            if has_pair and si_psi is not None and ox_psi is not None
            else None
        )
        if potential_split is not None:
            constraint_values.append(potential_split)

        row = {
            "bias_V": bias_v,
            "node_id": node,
            "x_um": float(raw["x"]),
            "y_um": float(raw["y"]),
            "is_direct_si_dielectric_interface_node": has_pair,
            "has_sentaurus_silicon_potential": si_psi is not None,
            "has_sentaurus_oxide_potential": ox_psi is not None,
            "sentaurus_silicon_potential_V": si_psi,
            "sentaurus_oxide_potential_V": ox_psi,
            "sentaurus_region_potential_split_V": potential_split,
            "global_barycentric_volume_um2": global_volume.get(node, 0.0),
            "silicon_averagebox_measure_um2": silicon_measure.get(node, 0.0),
            "silicon_measure_over_global_volume": volume_scale,
            "electron_flux": flux,
            "electron_source_global_volume": source,
            "electron_gauge": gauge,
            "electron_boundary": boundary,
            "level_a_shared_node_global_source_residual": level_a,
            "level_b_region_local_source_residual": level_b,
            "level_c_constrained_pair_residual": level_c,
            "level_b_minus_a": level_b - level_a,
            "level_c_minus_b": level_c - level_b,
        }
        rows.append(row)
        level_a_values.append(level_a)
        level_b_values.append(level_b)
        level_c_values.append(level_c)

    level_a_norm = norm(level_a_values)
    level_b_norm = norm(level_b_values)
    level_c_norm = norm(level_c_values)
    constraint_norm = norm(constraint_values) if constraint_values else {
        "l1": 0.0, "l2": 0.0, "maximum_abs": 0.0
    }
    summary = {
        "bias_V": bias_v,
        "selected_nodes": list(selected_nodes),
        "direct_interface_nodes": sum(
            row["is_direct_si_dielectric_interface_node"] for row in rows
        ),
        "interface_one_ring_only_nodes": [
            int(row["node_id"]) for row in rows
            if not row["is_direct_si_dielectric_interface_node"]
        ],
        "level_a_residual": level_a_norm,
        "level_b_residual": level_b_norm,
        "level_b_over_level_a": _ratios(level_b_norm, level_a_norm),
        "level_c_residual": level_c_norm,
        "level_c_over_level_b": _ratios(level_c_norm, level_b_norm),
        "sentaurus_pair_constraint_residual_V": constraint_norm,
        "maximum_absolute_level_c_minus_b": max(
            abs(float(row["level_c_minus_b"])) for row in rows
        ),
        "geometry": geometry_metadata,
    }
    return summary, rows


def audit(
    *,
    mesh: dict[str, Any],
    measures: dict[int, list[float]],
    coefficients: dict[int, list[float]],
    transport_materials: set[str],
    selected_nodes: tuple[int, ...],
    endpoints: list[dict[str, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    endpoint_summaries: list[dict[str, Any]] = []
    all_rows: list[dict[str, Any]] = []
    for endpoint in endpoints:
        summary, rows = audit_endpoint(
            bias_v=float(endpoint["bias_V"]),
            mesh=mesh,
            carrier_rows=endpoint["carrier_rows"],
            measures=measures,
            coefficients=coefficients,
            transport_materials=transport_materials,
            selected_nodes=selected_nodes,
            silicon_potential=endpoint["silicon_potential"],
            oxide_potential=endpoint["oxide_potential"],
        )
        endpoint_summaries.append(summary)
        all_rows.extend(rows)

    growth: dict[str, float] = {}
    if len(endpoint_summaries) == 2:
        low, high = sorted(endpoint_summaries, key=lambda item: item["bias_V"])
        for level in ("level_a", "level_b", "level_c"):
            growth[level] = (
                high[f"{level}_residual"]["l2"]
                / max(low[f"{level}_residual"]["l2"], 1.0e-300)
            )

    max_split = max(
        (
            abs(float(row["sentaurus_region_potential_split_V"]))
            for row in all_rows
            if row["sentaurus_region_potential_split_V"] is not None
        ),
        default=0.0,
    )
    max_b_change = max(
        (
            abs(float(row["level_b_minus_a"]))
            / max(abs(float(row["level_a_shared_node_global_source_residual"])),
                  1.0e-300)
            for row in all_rows
        ),
        default=0.0,
    )
    max_source = max(
        abs(float(row["electron_source_global_volume"])) for row in all_rows
    )
    max_residual = max(
        abs(float(row["level_a_shared_node_global_source_residual"]))
        for row in all_rows
    )
    report = {
        "schema": "vela.templates_ldmos.g3_seven_node_interface_ablation.v1",
        "contract": {
            "mode": "read_only_fixed_state",
            "state_variant": "VSV (Vela psi, Sentaurus phin, Vela phip)",
            "level_a": "shared node + silicon-only external AverageBox couples + global barycentric source volume",
            "level_b": "level A + silicon-side AverageBox source Measure",
            "level_c": "level B + exact psi master-slave continuity constraint, statically eliminated",
            "production_solver_modified": False,
            "full_self_consistent_double_node_solve_performed": False,
        },
        "selected_nodes": list(selected_nodes),
        "endpoints": endpoint_summaries,
        "endpoint_l2_growth_high_over_low": growth,
        "causal_gates": {
            "maximum_relative_row_change_from_region_local_source_volume": max_b_change,
            "maximum_absolute_electron_source_term": max_source,
            "maximum_source_over_maximum_level_a_residual": (
                max_source / max(max_residual, 1.0e-300)
            ),
            "maximum_sentaurus_si_oxide_potential_split_V": max_split,
            "level_b_material_if_any_row_changes_by_at_least_1_percent":
                max_b_change >= 0.01,
            "level_c_has_independent_frozen_carrier_effect": max_split > 1.0e-12,
        },
        "interpretation_limits": [
            "Level C proves only the exact-constraint Schur equivalence at frozen state.",
            "A nonzero double-node effect would require relaxed/discontinuous potential or different side-local Poisson assembly.",
            "A full self-consistent explicit-pair test requires new C++ degrees of freedom and interface equations.",
        ],
    }
    return report, all_rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--measure-coefficients", type=Path, required=True)
    parser.add_argument("--low-carrier-terms", type=Path, required=True)
    parser.add_argument("--low-sentaurus-export", type=Path, required=True)
    parser.add_argument("--low-bias", type=float, default=1.0)
    parser.add_argument("--high-carrier-terms", type=Path, required=True)
    parser.add_argument("--high-sentaurus-export", type=Path, required=True)
    parser.add_argument("--high-bias", type=float, default=1.1666666666666667)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--transport-materials", default="Si,Silicon")
    parser.add_argument(
        "--nodes", default=",".join(str(node) for node in DEFAULT_NODES)
    )
    args = parser.parse_args()

    debug = args.measure_coefficients.resolve()
    endpoints = []
    for bias, carrier, export in (
        (args.low_bias, args.low_carrier_terms, args.low_sentaurus_export),
        (args.high_bias, args.high_carrier_terms, args.high_sentaurus_export),
    ):
        fields = export.resolve() / "fields"
        endpoints.append({
            "bias_V": bias,
            "carrier_rows": read_csv(carrier.resolve()),
            "silicon_potential": read_scalar_field(
                fields / "ElectrostaticPotential_region0.csv"
            ),
            "oxide_potential": read_scalar_field(
                fields / "ElectrostaticPotential_region1.csv"
            ),
        })

    report, rows = audit(
        mesh=json.loads(args.mesh.read_text(encoding="utf-8")),
        measures=parse_debug_block(debug, "Measure"),
        coefficients=parse_debug_block(debug, "Coefficients"),
        transport_materials={
            value.strip().lower()
            for value in args.transport_materials.split(",") if value.strip()
        },
        selected_nodes=tuple(
            int(value) for value in args.nodes.split(",") if value.strip()
        ),
        endpoints=endpoints,
    )
    report["oracle"] = {"path": str(debug), "sha256": sha256(debug)}
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "seven_node_ablation.csv"
    json_path = output / "summary.json"
    write_csv(csv_path, rows)
    report["artifacts"] = {
        "node_table": str(csv_path),
        "summary": str(json_path),
    }
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "summary": str(json_path),
        "growth": report["endpoint_l2_growth_high_over_low"],
        "causal_gates": report["causal_gates"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
