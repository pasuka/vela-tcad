#!/usr/bin/env python3
"""Decompose Templates/LDMOS G3 current error by frozen state-family swaps.

The production G3 operator is held fixed with IALMob and predictor disabled.
For every selected exact bias, all eight combinations of Sentaurus/Vela
``psi``, electron QF and hole QF are replayed through the same SG operator.
This separates a read-only operator/coefficient mismatch from the nonlinear
self-consistent state feedback that remains after the contact-HFS repair.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import (
        drain_nodes,
        drain_sg_cut,
        read_csv,
        run_probe,
        run_sg_probe,
        write_sentaurus_state_csv,
    )
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_idvg_shift_kcl import (  # type: ignore[no-redef]
        drain_nodes,
        drain_sg_cut,
        read_csv,
        run_probe,
        run_sg_probe,
        write_sentaurus_state_csv,
    )


VARIANTS = {
    "VVV": ("V", "V", "V"),
    "SSS": ("S", "S", "S"),
    "SVV": ("S", "V", "V"),
    "VSV": ("V", "S", "V"),
    "VVS": ("V", "V", "S"),
    "SSV": ("S", "S", "V"),
    "SVS": ("S", "V", "S"),
    "VSS": ("V", "S", "S"),
}


def state_rows(path: Path) -> dict[int, dict[str, float]]:
    rows: dict[int, dict[str, float]] = {}
    for raw in read_csv(path):
        node = int(raw["node_id"])
        rows[node] = {
            key: float(raw[key])
            for key in (
                "psi", "phin", "phip", "electrons_m3", "holes_m3"
            )
        }
    return rows


def merge_state_rows(
    sentaurus: dict[int, dict[str, float]],
    vela: dict[int, dict[str, float]],
    owners: tuple[str, str, str],
) -> list[dict[str, float | int]]:
    if set(sentaurus) != set(vela):
        raise ValueError("Sentaurus and Vela states do not share the same nodes")
    psi_owner, phin_owner, phip_owner = owners
    sources = {"S": sentaurus, "V": vela}
    rows: list[dict[str, float | int]] = []
    for node in sorted(sentaurus):
        rows.append({
            "node_id": node,
            "psi": sources[psi_owner][node]["psi"],
            "phin": sources[phin_owner][node]["phin"],
            "phip": sources[phip_owner][node]["phip"],
            # The coupled SG diagnostic reconstructs n/p from psi and QFs.
            # Keep finite provenance values for the restart CSV contract.
            "electrons_m3": sources[phin_owner][node]["electrons_m3"],
            "holes_m3": sources[phip_owner][node]["holes_m3"],
        })
    return rows


def write_state(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                key: row[key] if key == "node_id" else format(float(row[key]), ".17g")
                for key in fields
            })


def find_export(roots: list[Path], index: int) -> Path:
    name = f"export_{index:04d}"
    matches = [root / name for root in roots if (root / name).is_dir()]
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected one {name} below export roots, found {len(matches)}"
        )
    return matches[0]


def curve_map(path: Path) -> dict[float, float]:
    result: dict[float, float] = {}
    for row in read_csv(path):
        bias = round(float(row["bias_V"]), 12)
        current_key = (
            "current_total_A_per_um"
            if "current_total_A_per_um" in row
            else "drain_current_A_per_um"
        )
        result[bias] = float(row[current_key])
    return result


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def state_delta(
    sentaurus: dict[int, dict[str, float]],
    vela: dict[int, dict[str, float]],
    field: str,
) -> dict[str, float]:
    values = [vela[node][field] - sentaurus[node][field] for node in sentaurus]
    center = statistics.median(values)
    centered = [value - center for value in values]
    return {
        "median_V": center,
        "centered_p95_abs_V": percentile([abs(value) for value in centered], 0.95),
        "max_abs_V": max(abs(value) for value in values),
    }


def attribution(variant_currents: dict[str, float], reference: float) -> dict[str, Any]:
    def error(name: str) -> float:
        return math.log10(max(abs(variant_currents[name] / reference), 1.0e-300))

    full_feedback = error("VVV") - error("SSS")
    single_sentaurus = {
        "psi": error("VVV") - error("SVV"),
        "phin": error("VVV") - error("VSV"),
        "phip": error("VVV") - error("VVS"),
    }
    qf_pair = error("VVV") - error("VSS")
    return {
        "operator_log10_ratio_dex": error("SSS"),
        "self_consistent_log10_ratio_dex": error("VVV"),
        "feedback_amplification_dex": full_feedback,
        "single_sentaurus_family_error_recovery_dex": single_sentaurus,
        "joint_qf_error_recovery_dex": qf_pair,
        "dominant_single_family": max(single_sentaurus, key=single_sentaurus.get),
    }


def carrier_term_summary(path: Path, limit: int = 10) -> dict[str, Any]:
    rows = read_csv(path)
    ranked = sorted(
        rows, key=lambda row: abs(float(row["electron_residual"])), reverse=True
    )
    residuals = [float(row["electron_residual"]) for row in rows]
    fluxes = [float(row["electron_flux"]) for row in rows]
    recombination = [float(row["electron_recombination"]) for row in rows]
    return {
        "electron_residual_l2": math.sqrt(sum(value * value for value in residuals)),
        "electron_residual_l1": sum(abs(value) for value in residuals),
        "electron_flux_l1": sum(abs(value) for value in fluxes),
        "electron_recombination_l1": sum(abs(value) for value in recombination),
        "top_electron_residual_rows": [
            {
                key: (int(row[key]) if key == "node_id" else float(row[key]))
                for key in (
                    "node_id", "x", "y", "electron_residual", "electron_flux",
                    "electron_flux_abs_sum", "electron_recombination",
                    "electron_boundary", "net_doping_m3",
                )
            }
            for row in ranked[:limit]
        ],
    }


def edge_feedback_summary(
    vela_edges_path: Path, sentaurus_phin_edges_path: Path, limit: int = 12,
) -> dict[str, Any]:
    vela = {int(row["edge_id"]): row for row in read_csv(vela_edges_path)}
    swapped = {
        int(row["edge_id"]): row for row in read_csv(sentaurus_phin_edges_path)
    }
    records = []
    for edge in sorted(set(vela) & set(swapped)):
        vrow, srow = vela[edge], swapped[edge]
        vdrop = float(vrow["phin1_V"]) - float(vrow["phin0_V"])
        sdrop = float(srow["phin1_V"]) - float(srow["phin0_V"])
        vflux = float(vrow["electron_particle_line_flux_per_m_s"])
        sflux = float(srow["electron_particle_line_flux_per_m_s"])
        records.append({
            "edge_id": edge,
            "node0": int(vrow["node0"]),
            "node1": int(vrow["node1"]),
            "x_mid": 0.5 * (float(vrow["x0"]) + float(vrow["x1"])),
            "y_mid": 0.5 * (float(vrow["y0"]) + float(vrow["y1"])),
            "couple_over_length": float(vrow["couple_m"]) / max(
                float(vrow["length_m"]), 1.0e-300
            ),
            "vela_phin_drop_V": vdrop,
            "sentaurus_phin_drop_V": sdrop,
            "phin_drop_difference_V": vdrop - sdrop,
            "vela_electron_line_flux": vflux,
            "sentaurus_phin_electron_line_flux": sflux,
            "electron_line_flux_difference": vflux - sflux,
        })
    drop_differences = [abs(row["phin_drop_difference_V"]) for row in records]
    return {
        "common_edges": len(records),
        "phin_drop_difference_abs_V": {
            "median": statistics.median(drop_differences),
            "p95": percentile(drop_differences, 0.95),
            "maximum": max(drop_differences),
        },
        "top_flux_feedback_edges": sorted(
            records,
            key=lambda row: abs(row["electron_line_flux_difference"]),
            reverse=True,
        )[:limit],
        "top_phin_drop_difference_edges": sorted(
            records,
            key=lambda row: abs(row["phin_drop_difference_V"]),
            reverse=True,
        )[:limit],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument(
        "--sentaurus-export-root", type=Path, action="append", required=True
    )
    parser.add_argument("--vela-state-root", type=Path, required=True)
    parser.add_argument("--reference-curve", type=Path, required=True)
    parser.add_argument("--vela-curve", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--indices", default="1,2,3,4")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    mobility_text = json.dumps(baseline["solver"]["mobility"]).lower()
    if "ialmob" in mobility_text:
        raise ValueError("state feedback audit requires IALMob disabled")
    if not baseline["solver"]["mobility"].get(
        "contact_electric_field_fallback", False
    ):
        raise ValueError("state feedback audit requires qualified contact HFS fallback")
    if "predictor" in json.dumps(baseline.get("sweep", {})).lower():
        raise ValueError("state feedback audit requires predictor disabled")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    contact_nodes = drain_nodes(args.mesh)
    references = curve_map(args.reference_curve)
    vela_curve = curve_map(args.vela_curve)
    indices = [int(item) for item in args.indices.split(",") if item.strip()]
    points: list[dict[str, Any]] = []

    for index in indices:
        bias = index / 6.0
        label = f"vg_{bias:.6f}".replace(".", "p")
        point_dir = output / label
        sentaurus_path = write_sentaurus_state_csv(
            find_export(args.sentaurus_export_root, index),
            args.mesh,
            point_dir / "sentaurus_state.csv",
        )
        vela_path = args.vela_state_root / (
            "g3_contact_hfs_point_bias_" + f"{bias:.6f}".replace(".", "p") + ".csv"
        )
        if not vela_path.is_file():
            raise FileNotFoundError(vela_path)
        sentaurus = state_rows(sentaurus_path)
        vela = state_rows(vela_path)
        currents: dict[str, float] = {}
        variants: dict[str, Any] = {}
        for name, owners in VARIANTS.items():
            state = point_dir / name / "hybrid_state.csv"
            if not state.is_file():
                write_state(state, merge_state_rows(sentaurus, vela, owners))
            variant_dir = point_dir / name
            edge_path = variant_dir / "sg_edges.csv"
            if edge_path.is_file():
                cut = drain_sg_cut(read_csv(edge_path), contact_nodes)
                probe = {
                    "drain_cut": cut,
                    "artifacts": {
                        "config": str((variant_dir / "sg_probe.json").resolve()),
                        "sg_edges": str(edge_path.resolve()),
                    },
                }
            else:
                probe = run_sg_probe(
                    args.runner, baseline, state, bias, variant_dir
                )
            current = float(probe["drain_cut"]["total_A_per_um"])
            currents[name] = current
            variants[name] = {
                "owners": {"psi": owners[0], "phin": owners[1], "phip": owners[2]},
                "current_A_per_um": current,
                "ratio_to_sentaurus": current / references[round(bias, 12)],
                "artifacts": probe["artifacts"],
            }
        carrier_terms: dict[str, Any] = {}
        for name in ("SSS", "VVV", "VSV"):
            variant_dir = point_dir / name
            carrier_path = variant_dir / "carrier_terms.csv"
            if not carrier_path.is_file():
                run_probe(
                    args.runner, baseline,
                    variant_dir / "hybrid_state.csv", bias, variant_dir,
                )
            carrier_terms[name] = {
                **carrier_term_summary(carrier_path),
                "artifact": str(carrier_path.resolve()),
            }
        point = {
            "bias_V": bias,
            "sentaurus_terminal_A_per_um": references[round(bias, 12)],
            "vela_curve_A_per_um": vela_curve[round(bias, 12)],
            "vela_curve_to_probe_relative_difference": (
                currents["VVV"] / vela_curve[round(bias, 12)] - 1.0
            ),
            "state_delta_vela_minus_sentaurus": {
                field: state_delta(sentaurus, vela, field)
                for field in ("psi", "phin", "phip")
            },
            "variants": variants,
            "attribution": attribution(currents, references[round(bias, 12)]),
            "carrier_continuity": carrier_terms,
            "edge_feedback": edge_feedback_summary(
                point_dir / "VVV" / "sg_edges.csv",
                point_dir / "VSV" / "sg_edges.csv",
            ),
        }
        points.append(point)

    feedback = [point["attribution"]["feedback_amplification_dex"] for point in points]
    summary = {
        "schema": "vela.templates_ldmos.g3_state_feedback.v1",
        "contracts": {
            "ialmob": "disabled",
            "predictor": "disabled",
            "operator": "qualified_contact_hfs_g3",
            "state_substitution": "read_only_psi_phin_phip_factorial",
            "pn2d_atomic_profile_inherited": False,
        },
        "points": points,
        "aggregate": {
            "median_feedback_amplification_dex": statistics.median(feedback),
            "median_feedback_amplification_factor": 10.0 ** statistics.median(feedback),
            "dominant_single_family_counts": {
                field: sum(
                    point["attribution"]["dominant_single_family"] == field
                    for point in points
                )
                for field in ("psi", "phin", "phip")
            },
            "ledger_status": "draft",
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), **summary["aggregate"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
