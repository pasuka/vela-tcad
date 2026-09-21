#!/usr/bin/env python3
"""Replay G4 Sentaurus states through LDMOS SG coefficient controls.

This is a fixed-state, single-factor audit.  It does not enable IALMob or a
predictor, and it deliberately does not inherit the PN2D atomic profile.  The
two edge-coefficient variants are the existing positive barycentric fallback
and truncated Voronoi (zero negative local cotangent).  Mixed-Voronoi node
volume is included only as an orthogonal invariance control.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from copy import deepcopy
from pathlib import Path
from typing import Any

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import (
        BIAS_POINTS_V,
        drain_nodes,
        read_csv,
        run_sg_probe,
        write_sentaurus_state_csv,
    )
except ModuleNotFoundError:  # Direct ``python scripts/<name>.py`` execution.
    from audit_templates_ldmos_g3_idvg_shift_kcl import (
        BIAS_POINTS_V,
        drain_nodes,
        read_csv,
        run_sg_probe,
        write_sentaurus_state_csv,
    )


VARIANTS: dict[str, dict[str, Any]] = {
    "barycentric_fallback": {
        "node_volume_policy": "barycentric",
        "fallback_negative_cotangent": True,
    },
    "truncated_voronoi": {
        "node_volume_policy": "barycentric",
        "fallback_negative_cotangent": False,
    },
    "mixed_volume_control": {
        "node_volume_policy": "mixed_voronoi",
        "fallback_negative_cotangent": True,
    },
    "mixed_truncated_control": {
        "node_volume_policy": "mixed_voronoi",
        "fallback_negative_cotangent": False,
    },
}


def reference_currents(path: Path) -> dict[float, float]:
    return {
        round(float(row["bias_V"]), 12): float(row["current_total_A_per_um"])
        for row in read_csv(path)
    }


def variant_config(baseline: dict[str, Any], name: str) -> dict[str, Any]:
    config = deepcopy(baseline)
    geometry = config.setdefault("mesh_geometry", {})
    geometry.update(VARIANTS[name])
    config["_comment"] = (
        "Templates/LDMOS G4 fixed-state coefficient qualification; "
        f"variant={name}; no IALMob, predictor, or PN2D profile inheritance."
    )
    return config


def edge_record(row: dict[str, str]) -> dict[str, float | int]:
    length = float(row["length_m"])
    return {
        "edge_id": int(row["edge_id"]),
        "node0": int(row["node0"]),
        "node1": int(row["node1"]),
        "length_m": length,
        "couple_m": float(row["couple_m"]),
        "couple_over_length": float(row["couple_m"]) / max(length, 1.0e-300),
        "net_doping_avg_m3": float(row["net_doping_avg_m3"]),
        "electron_density0_m3": float(row["electron_density0_m3"]),
        "electron_density1_m3": float(row["electron_density1_m3"]),
        "electron_mobility_m2_V_s": float(row["electron_mobility_m2_V_s"]),
        "electron_particle_line_flux_per_m_s": float(
            row["electron_particle_line_flux_per_m_s"]
        ),
    }


def edge_delta_summary(
    baseline_rows: list[dict[str, str]], candidate_rows: list[dict[str, str]],
) -> dict[str, float | int]:
    base = {int(row["edge_id"]): edge_record(row) for row in baseline_rows}
    candidate = {
        int(row["edge_id"]): edge_record(row) for row in candidate_rows
    }
    common = sorted(set(base) & set(candidate))
    couple_changed = [
        edge for edge in common
        if not math.isclose(
            float(base[edge]["couple_m"]), float(candidate[edge]["couple_m"]),
            rel_tol=1.0e-14, abs_tol=1.0e-30,
        )
    ]
    flux_relative = [
        abs(
            float(candidate[edge]["electron_particle_line_flux_per_m_s"])
            - float(base[edge]["electron_particle_line_flux_per_m_s"])
        ) / max(
            abs(float(base[edge]["electron_particle_line_flux_per_m_s"])),
            1.0e-300,
        )
        for edge in common
    ]
    return {
        "common_edges": len(common),
        "couple_changed_edges": len(couple_changed),
        "maximum_couple_absolute_change_m": max(
            (
                abs(float(candidate[edge]["couple_m"]) - float(base[edge]["couple_m"]))
                for edge in common
            ),
            default=0.0,
        ),
        "maximum_flux_relative_change": max(flux_relative, default=0.0),
    }


def write_aligned_edges(
    output: Path, bias: float, rows_by_variant: dict[str, list[dict[str, str]]],
    contact_nodes: set[int],
) -> None:
    indexed = {
        name: {int(row["edge_id"]): edge_record(row) for row in rows}
        for name, rows in rows_by_variant.items()
    }
    edges = sorted(set.intersection(*(set(rows) for rows in indexed.values())))
    fields = [
        "bias_V", "edge_id", "node0", "node1", "crosses_drain",
        "variant", "length_m", "couple_m", "couple_over_length",
        "net_doping_avg_m3", "electron_density0_m3",
        "electron_density1_m3", "electron_mobility_m2_V_s",
        "electron_particle_line_flux_per_m_s",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for edge in edges:
            for name in VARIANTS:
                record = indexed[name][edge]
                writer.writerow({
                    "bias_V": format(bias, ".17g"),
                    "edge_id": edge,
                    "node0": record["node0"],
                    "node1": record["node1"],
                    "crosses_drain": int(
                        (int(record["node0"]) in contact_nodes)
                        != (int(record["node1"]) in contact_nodes)
                    ),
                    "variant": name,
                    **{key: record[key] for key in fields[6:]},
                })


def classify_branch(
    g4_ratios: list[float], g3_ratios: list[float],
    g4_curve_ratios: list[float],
) -> dict[str, Any]:
    log_g4 = statistics.median(math.log10(abs(value)) for value in g4_ratios)
    log_g3 = statistics.median(math.log10(abs(value)) for value in g3_ratios)
    log_curve = statistics.median(
        math.log10(abs(value)) for value in g4_curve_ratios
    )
    distance_to_g3 = abs(log_g4 - log_g3)
    distance_to_curve = abs(log_g4 - log_curve)
    return {
        "classification": (
            "plot_time_supported"
            if distance_to_g3 < distance_to_curve
            else "assembly_real_supported"
        ),
        "median_g4_fixed_replay_log10_ratio": log_g4,
        "median_g3_fixed_replay_log10_ratio": log_g3,
        "median_g4_self_consistent_curve_log10_ratio": log_curve,
        "distance_to_g3_fixed_replay_dex": distance_to_g3,
        "distance_to_g4_self_consistent_curve_dex": distance_to_curve,
        "interpretation": (
            "Nearest-anchor classification only: it separates whether removal "
            "of HFS in the reference state collapses the replay mismatch."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--sentaurus-export-root", type=Path, required=True)
    parser.add_argument("--reference-curve", type=Path, required=True)
    parser.add_argument("--g3-summary", type=Path, required=True)
    parser.add_argument("--g4-control-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    reference = reference_currents(args.reference_curve)
    mesh_path = Path(baseline["mesh_file"])
    contacts = drain_nodes(mesh_path)
    summary: dict[str, Any] = {
        "schema": "vela.templates_ldmos_g4_fixed_state_coefficients.v1",
        "scope": {
            "bias_indices": [1, 3, 5],
            "ialmob_enabled": False,
            "predictor_enabled": False,
            "pn2d_atomic_profile_inherited": False,
            "variants": VARIANTS,
        },
        "points": [],
    }
    for index in (1, 3, 5):
        bias = BIAS_POINTS_V[index]
        point_dir = output / f"vg_{bias:.6f}".replace(".", "p")
        state = write_sentaurus_state_csv(
            args.sentaurus_export_root / f"export_{index:04d}", mesh_path,
            point_dir / "sentaurus_state.csv",
        )
        rows_by_variant: dict[str, list[dict[str, str]]] = {}
        point: dict[str, Any] = {"bias_V": bias, "variants": {}}
        for name in VARIANTS:
            case_dir = point_dir / name
            probe = run_sg_probe(
                args.runner, variant_config(baseline, name), state, bias, case_dir,
            )
            rows = read_csv(Path(probe["artifacts"]["sg_edges"]))
            rows_by_variant[name] = rows
            current = float(probe["drain_cut"]["total_A_per_um"])
            ref = reference[round(bias, 12)]
            point["variants"][name] = {
                "sentaurus_terminal_A_per_um": ref,
                "vela_sg_on_sentaurus_state_A_per_um": current,
                "signed_ratio": current / ref,
                "magnitude_error_dex": abs(
                    math.log10(max(abs(current), 1.0e-300))
                    - math.log10(max(abs(ref), 1.0e-300))
                ),
                "drain_cut": probe["drain_cut"],
                "dominant_drain_edges": probe["dominant_drain_edges"],
                "artifacts": probe["artifacts"],
            }
        base_rows = rows_by_variant["barycentric_fallback"]
        point["edge_deltas_from_barycentric_fallback"] = {
            name: edge_delta_summary(base_rows, rows)
            for name, rows in rows_by_variant.items()
            if name != "barycentric_fallback"
        }
        write_aligned_edges(
            point_dir / "all_edge_coefficient_ab.csv", bias,
            rows_by_variant, contacts,
        )
        summary["points"].append(point)

    g3 = json.loads(args.g3_summary.read_text(encoding="utf-8"))
    g4_control = json.loads(args.g4_control_summary.read_text(encoding="utf-8"))
    g4_ratios = [
        point["variants"]["barycentric_fallback"]["signed_ratio"]
        for point in summary["points"]
    ]
    g3_ratios = [
        point["signed_ratio"]
        for point in g3["sentaurus_state_sg_replay"]["points"]
    ]
    selected = {round(BIAS_POINTS_V[index], 12) for index in (1, 3, 5)}
    g4_curve_ratios = [
        point["candidate_A_per_um"] / point["reference_A_per_um"]
        for point in g4_control["high_field_off_control"]["point_metrics"]
        if round(float(point["bias_V"]), 12) in selected
    ]
    summary["branch_decision"] = classify_branch(
        g4_ratios, g3_ratios, g4_curve_ratios,
    )
    summary["acceptance"] = {
        "mixed_volume_fixed_sg_invariant": all(
            point["edge_deltas_from_barycentric_fallback"]
            ["mixed_volume_control"]["maximum_flux_relative_change"] <= 1.0e-14
            for point in summary["points"]
        ),
        "ledger_status": "draft",
        "ledger_approval_allowed": False,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8",
    )
    print(json.dumps({
        "branch_decision": summary["branch_decision"],
        "acceptance": summary["acceptance"],
        "summary": str((output / "summary.json").resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
