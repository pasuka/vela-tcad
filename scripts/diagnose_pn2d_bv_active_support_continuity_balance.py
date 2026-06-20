#!/usr/bin/env python3
"""Summarize continuity-balance terms on active PN2D BV support nodes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "active_edge_count",
    "active_edge_ids",
    "node_volume_m2",
    "vela_psi_minus_phin_V",
    "sentaurus_psi_minus_phin_V",
    "delta_psi_minus_phin_V",
    "vela_phip_minus_psi_V",
    "sentaurus_phip_minus_psi_V",
    "delta_phip_minus_psi_V",
    "vela_electron_density_cm3",
    "sentaurus_electron_density_cm3",
    "vela_electron_density_over_sentaurus",
    "vela_hole_density_cm3",
    "sentaurus_hole_density_cm3",
    "vela_hole_density_over_sentaurus",
    "vela_electron_transport_integral_s_inv",
    "vela_hole_transport_integral_s_inv",
    "vela_edge_avalanche_node_integral_s_inv",
    "vela_vtk_avalanche_node_integral_s_inv",
    "sentaurus_generation_node_integral_s_inv",
    "vela_edge_avalanche_over_sentaurus_generation",
    "vela_vtk_avalanche_over_sentaurus_generation",
    "vela_srh_node_integral_s_inv",
    "sentaurus_srh_node_integral_s_inv",
    "vela_srh_over_sentaurus_srh",
    "vela_electron_residual_estimate_s_inv",
    "vela_hole_residual_estimate_s_inv",
    "vela_electron_residual_over_edge_avalanche",
    "vela_hole_residual_over_edge_avalanche",
]


SUMMARY_FIELDS = [
    "delta_psi_minus_phin_V",
    "delta_phip_minus_psi_V",
    "vela_electron_density_over_sentaurus",
    "vela_hole_density_over_sentaurus",
    "vela_edge_avalanche_over_sentaurus_generation",
    "vela_vtk_avalanche_over_sentaurus_generation",
    "vela_srh_over_sentaurus_srh",
    "vela_electron_residual_over_edge_avalanche",
    "vela_hole_residual_over_edge_avalanche",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--active-axis", choices=["x", "y"], default="x")
    parser.add_argument("--active-source-relative-threshold", type=float, default=1.0e-12)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize_group(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {"count": len(rows)}
    for key in SUMMARY_FIELDS:
        values = finite_values(rows, key)
        if values:
            summary[f"{key}_median"] = statistics.median(values)
            summary[f"{key}_min"] = min(values)
            summary[f"{key}_max"] = max(values)
    return summary


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def active_edge_ids(active_by_node: dict[int, list[dict[str, Any]]], node_id: int) -> str:
    ids = sorted(int(item["edge_id"]) for item in active_by_node.get(node_id, []))
    return ";".join(str(edge_id) for edge_id in ids)


def support_nodes(path: Path) -> list[dict[str, str]]:
    rows = []
    for row in read_csv_rows(path):
        if row.get("support_class") in (None, "", "inactive"):
            continue
        rows.append(row)
    return rows


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, contacts = cont.sgdiag.read_mesh(args.mesh)
    edges = cont.sgdiag.build_edges(nodes, triangles)
    volumes = cont.sgdiag.node_volumes(nodes, triangles)
    vt = cont.sgdiag.K_B_OVER_Q * args.temperature_k
    vtk = fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    scalars = cont.sgdiag.parse_vtk_scalars(vtk)
    ni = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    source_rows, node_edge_avalanche = cont.source_rows_with_signed_flux(
        nodes,
        edges,
        contacts,
        scalars,
        ni,
        vt,
    )
    electron_transport, hole_transport = cont.compute_signed_terms(source_rows)
    active_by_node = fluxfac.load_active_cxx_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    sent_nodes = cont.load_sentaurus_nodes(args.sentaurus_dir)
    sent_fields = cont.load_sentaurus_fields(args.sentaurus_dir)
    srh = scalars.get("SRHRecombination", [0.0 for _ in nodes])
    avalanche = scalars.get("AvalancheGeneration", [0.0 for _ in nodes])

    rows: list[dict[str, Any]] = []
    for support in support_nodes(args.support_csv):
        node_id = int(support["node_id"])
        node = nodes[node_id]
        sent_id, sent_distance = cont.nearest_node(sent_nodes, node["x_um"], node["y_um"])
        volume = volumes[node_id]
        vela_psi = scalars["Potential"][node_id]
        vela_phin = scalars["ElectronQuasiFermi"][node_id]
        vela_phip = scalars["HoleQuasiFermi"][node_id]
        sent_psi = cont.scalar_at(sent_fields, "ElectrostaticPotential", sent_id)
        sent_phin = cont.scalar_at(sent_fields, "eQuasiFermiPotential", sent_id)
        sent_phip = cont.scalar_at(sent_fields, "hQuasiFermiPotential", sent_id)
        vela_n = scalars["Electrons"][node_id] * 1.0e-6
        vela_p = scalars["Holes"][node_id] * 1.0e-6
        sent_n = cont.scalar_at(sent_fields, "eDensity", sent_id)
        sent_p = cont.scalar_at(sent_fields, "hDensity", sent_id)
        vela_edge_avalanche = node_edge_avalanche[node_id]
        vela_vtk_avalanche = avalanche[node_id] * volume
        sent_generation = cont.scalar_at(sent_fields, "ImpactIonization", sent_id) * 1.0e6 * volume
        vela_srh = srh[node_id] * volume
        sent_srh = cont.scalar_at(sent_fields, "srhRecombination", sent_id) * 1.0e6 * volume
        electron_residual = electron_transport[node_id] + vela_srh - vela_edge_avalanche
        hole_residual = hole_transport[node_id] + vela_srh - vela_edge_avalanche
        rows.append({
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "support_class": support["support_class"],
            "sentaurus_node_id": sent_id,
            "sentaurus_distance_um": sent_distance,
            "active_edge_count": len(active_by_node.get(node_id, [])),
            "active_edge_ids": active_edge_ids(active_by_node, node_id),
            "node_volume_m2": volume,
            "vela_psi_minus_phin_V": vela_psi - vela_phin,
            "sentaurus_psi_minus_phin_V": sent_psi - sent_phin,
            "delta_psi_minus_phin_V": (vela_psi - vela_phin) - (sent_psi - sent_phin),
            "vela_phip_minus_psi_V": vela_phip - vela_psi,
            "sentaurus_phip_minus_psi_V": sent_phip - sent_psi,
            "delta_phip_minus_psi_V": (vela_phip - vela_psi) - (sent_phip - sent_psi),
            "vela_electron_density_cm3": vela_n,
            "sentaurus_electron_density_cm3": sent_n,
            "vela_electron_density_over_sentaurus": ratio(vela_n, sent_n),
            "vela_hole_density_cm3": vela_p,
            "sentaurus_hole_density_cm3": sent_p,
            "vela_hole_density_over_sentaurus": ratio(vela_p, sent_p),
            "vela_electron_transport_integral_s_inv": electron_transport[node_id],
            "vela_hole_transport_integral_s_inv": hole_transport[node_id],
            "vela_edge_avalanche_node_integral_s_inv": vela_edge_avalanche,
            "vela_vtk_avalanche_node_integral_s_inv": vela_vtk_avalanche,
            "sentaurus_generation_node_integral_s_inv": sent_generation,
            "vela_edge_avalanche_over_sentaurus_generation": ratio(vela_edge_avalanche, sent_generation),
            "vela_vtk_avalanche_over_sentaurus_generation": ratio(vela_vtk_avalanche, sent_generation),
            "vela_srh_node_integral_s_inv": vela_srh,
            "sentaurus_srh_node_integral_s_inv": sent_srh,
            "vela_srh_over_sentaurus_srh": ratio(vela_srh, sent_srh),
            "vela_electron_residual_estimate_s_inv": electron_residual,
            "vela_hole_residual_estimate_s_inv": hole_residual,
            "vela_electron_residual_over_edge_avalanche": ratio(electron_residual, vela_edge_avalanche),
            "vela_hole_residual_over_edge_avalanche": ratio(hole_residual, vela_edge_avalanche),
        })
    return rows


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "active_support_continuity_balance_nodes.csv", rows)
    by_class: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_class.setdefault(str(row["support_class"]), []).append(row)
    summary = {
        "bias": args.bias,
        "row_count": len(rows),
        "support_classes": {
            support_class: summarize_group(items)
            for support_class, items in sorted(by_class.items())
        },
        "all": summarize_group(rows),
    }
    (args.out_dir / "active_support_continuity_balance_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
