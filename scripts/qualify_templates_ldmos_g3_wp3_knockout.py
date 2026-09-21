#!/usr/bin/env python3
"""Apply the frozen R4 qualification and residual-reduction gates to one T4 endpoint."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Sequence


NP_P95_LIMIT_DEX = 1.0e-3
NP_MAX_LIMIT_DEX = 1.0e-2
CONTACT_QF_LIMIT_V = 1.0e-9
NEUTRALITY_LIMIT_REL = 1.0e-6
PORT_LIMIT_REL = 1.0e-3
DENOMINATOR_PORT_FRACTION = 1.0e-6
FIXED_HOTSPOTS = (3721, 3974, 3973, 4091, 3727, 4021, 3771)
VALID_VARIANTS = {"k1_hfs_off", "k2_boltzmann", "k3_bgn_off"}


def percentile95(values: Sequence[float]) -> float:
    ordered = sorted(abs(float(value)) for value in values)
    if not ordered:
        raise ValueError("cannot compute P95 of an empty sequence")
    rank = 0.95 * (len(ordered) - 1)
    low = math.floor(rank)
    high = math.ceil(rank)
    if low == high:
        return ordered[low]
    return ordered[low] + (rank - low) * (ordered[high] - ordered[low])


def l2(values: Sequence[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in values))


def require_sha256_map(raw: dict[str, Any]) -> dict[str, str]:
    required = {
        "deck", "grid", "parameter", "sentaurus_tdr", "imported_state",
        "vela_replay_config", "transport_edges",
    }
    if set(raw) != required:
        raise ValueError(f"artifact_sha256 keys must be exactly {sorted(required)}")
    for name, digest in raw.items():
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in "0123456789abcdef" for c in digest):
            raise ValueError(f"invalid SHA-256 for {name}")
    return {name: raw[name] for name in sorted(raw)}


def qualify(raw: dict[str, Any]) -> dict[str, Any]:
    variant = str(raw["variant"])
    if variant not in VALID_VARIANTS:
        raise ValueError(f"unregistered variant: {variant}")
    bias = float(raw["bias_V"])
    if min(abs(bias - 1.0), abs(bias - 1.166666666666667)) > 1.0e-12:
        raise ValueError(f"unregistered endpoint bias: {bias}")

    electron_errors = [abs(float(value)) for value in raw["np_log10_errors_dex"]["electron"]]
    hole_errors = [abs(float(value)) for value in raw["np_log10_errors_dex"]["hole"]]
    np_p95 = max(percentile95(electron_errors), percentile95(hole_errors))
    np_max = max(electron_errors + hole_errors)

    contact = raw["contact"]
    qf_max = max(
        max((abs(float(value)) for value in contact["phin_minus_contact_V"]), default=0.0),
        max((abs(float(value)) for value in contact["phip_minus_contact_V"]), default=0.0),
    )
    neutrality_max = max((abs(float(value)) for value in contact["neutrality_relative_residual"]), default=0.0)

    port = raw["port_current_A_per_um"]
    sentaurus_current = float(port["sentaurus"])
    vela_current = float(port["vela"])
    if sentaurus_current == 0.0:
        raise ValueError("Sentaurus port current must be nonzero")
    port_rel = abs(vela_current / sentaurus_current - 1.0)

    mesh = raw["mesh"]
    mesh_equal = bool(
        mesh["vertex_count_equal"]
        and mesh["coordinate_sha256_equal"]
        and mesh["transport_edge_sha256_equal"]
    )

    residual = raw["residual"]
    fixed_rows = [float(value) for value in residual["fixed_seven_rows_A_per_um"]]
    if len(fixed_rows) != 7:
        raise ValueError("fixed_seven_rows_A_per_um must contain exactly seven rows")
    top_rows = [float(value) for value in residual["variant_top7_rows_A_per_um"]]
    top_ids = [int(value) for value in residual["variant_top7_node_ids"]]
    if len(top_rows) != 7 or len(top_ids) != 7 or len(set(top_ids)) != 7:
        raise ValueError("variant top-7 must contain seven unique nodes and seven rows")
    edge_ids = [int(value) for value in residual["incident_silicon_edge_ids"]]
    fluxes = [abs(float(value)) for value in residual["incident_silicon_edge_flux_A_per_um"]]
    if len(edge_ids) != 40 or len(fluxes) != 40 or len(set(edge_ids)) != 40:
        raise ValueError("fixed-seven denominator must contain the 40 unique incident silicon edges")
    r7_l2 = l2(fixed_rows)
    sum_abs_phi = sum(fluxes)
    # R4's I_port^variant is frozen here to the Sentaurus terminal current of
    # the independently converged variant endpoint.  The preceding 1e-3 port
    # gate ensures the matching Vela fixed-state current is equivalent.
    floor = DENOMINATOR_PORT_FRACTION * abs(sentaurus_current)
    floor_hit = sum_abs_phi < floor
    denominator = max(sum_abs_phi, floor)
    normalized = r7_l2 / denominator if denominator > 0.0 else math.nan
    baseline = float(raw["baseline_normalized_residual"])
    if not math.isfinite(baseline) or baseline <= 0.0:
        raise ValueError("baseline_normalized_residual must be finite and positive")
    ratio = normalized / baseline
    reduction = baseline / normalized if normalized > 0.0 else sys.float_info.max

    gates = {
        "np_reconstruction": np_p95 <= NP_P95_LIMIT_DEX and np_max <= NP_MAX_LIMIT_DEX,
        "contact_rows": qf_max <= CONTACT_QF_LIMIT_V and neutrality_max <= NEUTRALITY_LIMIT_REL,
        "fixed_state_port_current": port_rel <= PORT_LIMIT_REL,
        "mesh_and_transport_edges": mesh_equal,
    }
    if not all(gates.values()):
        verdict = "contract_not_equivalent"
    elif floor_hit:
        verdict = "uncertain"
    elif reduction > 10.0:
        verdict = "discretization_carrier"
    elif reduction < 2.0:
        verdict = "family_ruled_out"
    else:
        verdict = "uncertain"

    return {
        "schema": "vela.templates_ldmos.g3_wp3_knockout_summary.v1",
        "benchmark": "sentaurus_t2022_03_sp2_templates_ldmos",
        "variant": variant,
        "bias_V": bias,
        "gate_np_p95_dex": np_p95,
        "gate_np_max_dex": np_max,
        "gate_contact_qf_max_V": qf_max,
        "gate_neutrality_max_rel": neutrality_max,
        "gate_port_rel": port_rel,
        "gate_mesh_equal": mesh_equal,
        "qualification_gates": gates,
        "r7_l2": r7_l2,
        "sum_abs_phi": sum_abs_phi,
        "incident_edge_count": len(edge_ids),
        "denominator": denominator,
        "denominator_floor": floor,
        "denominator_port_current_source": "sentaurus_variant_endpoint",
        "residual_and_flux_unit": "A/um",
        "denominator_floor_hit": floor_hit,
        "normalized_residual": normalized,
        "baseline_normalized_residual": baseline,
        "ratio_vs_baseline": ratio,
        "reduction_factor_vs_baseline": reduction,
        "hotspot_set": {
            "fixed_seven_node_ids": list(FIXED_HOTSPOTS),
            "variant_top7_node_ids": top_ids,
        },
        "residual_sets": {
            "fixed_seven_l2_A_per_um": r7_l2,
            "fixed_seven_max_abs_A_per_um": max(abs(value) for value in fixed_rows),
            "variant_top7_l2_A_per_um": l2(top_rows),
            "variant_top7_max_abs_A_per_um": max(abs(value) for value in top_rows),
        },
        "artifact_sha256": require_sha256_map(raw["artifact_sha256"]),
        "verdict": verdict,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    summary = qualify(json.loads(args.input.read_text(encoding="utf-8")))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["verdict"] != "contract_not_equivalent" else 2


if __name__ == "__main__":
    raise SystemExit(main())
