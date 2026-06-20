#!/usr/bin/env python3
"""Compare qF-shifted active-support continuity balance against exported current."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_mixed_state_replay as mixed
import diagnose_pn2d_bv_active_support_residual_proxy as residual
import diagnose_pn2d_bv_edge_source_current_consistency as consistency
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "bias_V",
    "variant",
    "node_id",
    "x_um",
    "y_um",
    "support_class",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "active_edge_count",
    "active_edge_ids",
    "active_endpoint_area_sum_m2",
    "node_volume_m2",
    "transport_model",
    "sentaurus_support_particle_flux_m2_s",
    "sentaurus_active_current_integral_s_inv",
    "particle_flux_avg_m2_s",
    "particle_flux_over_sentaurus_support_current",
    "electron_transport_s_inv",
    "hole_transport_s_inv",
    "cxx_edge_source_s_inv",
    "replayed_edge_source_s_inv",
    "sentaurus_node_source_s_inv",
    "srh_source_s_inv",
    "cxx_edge_source_over_sentaurus_node_source",
    "replayed_edge_source_over_sentaurus_node_source",
    "cxx_edge_source_over_sentaurus_active_current_integral",
    "replayed_edge_source_over_sentaurus_active_current_integral",
    "electron_transport_over_sentaurus_node_source",
    "hole_transport_over_sentaurus_node_source",
    "electron_transport_over_cxx_edge_source",
    "hole_transport_over_cxx_edge_source",
    "electron_residual_cxx_source_s_inv",
    "hole_residual_cxx_source_s_inv",
    "electron_residual_replayed_source_s_inv",
    "hole_residual_replayed_source_s_inv",
    "electron_residual_sentaurus_source_s_inv",
    "hole_residual_sentaurus_source_s_inv",
    "electron_residual_over_impact",
    "hole_residual_over_impact",
    "electron_residual_replayed_over_replayed_impact",
    "hole_residual_replayed_over_replayed_impact",
    "electron_residual_sentaurus_over_sentaurus_impact",
    "hole_residual_sentaurus_over_sentaurus_impact",
]

SUMMARY_KEYS = [
    "particle_flux_over_sentaurus_support_current",
    "cxx_edge_source_over_sentaurus_node_source",
    "replayed_edge_source_over_sentaurus_node_source",
    "electron_transport_over_sentaurus_node_source",
    "hole_transport_over_sentaurus_node_source",
    "electron_transport_over_cxx_edge_source",
    "hole_transport_over_cxx_edge_source",
    "electron_residual_over_impact",
    "hole_residual_over_impact",
    "electron_residual_replayed_over_replayed_impact",
    "hole_residual_replayed_over_replayed_impact",
    "electron_residual_sentaurus_over_sentaurus_impact",
    "hole_residual_sentaurus_over_sentaurus_impact",
]

VARIANT_MAP = {
    "vela_baseline": "vela_baseline",
    "vela_qf_shift": "vela_qf_shift",
    "vela_psi_sentaurus_qf": "vela_psi_sentaurus_qf",
}


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
    parser.add_argument("--transport-model", choices=sorted(residual.VALID_TRANSPORT_MODES),
                        default="qf_old_slotboom_ni")
    parser.add_argument("--electron-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--hole-qf-shift-v", type=float, default=0.0)
    return parser.parse_args()


def avg2(a: float, b: float) -> float:
    return 0.5 * (a + b)


def active_endpoint_area(active_items: list[dict[str, Any]]) -> float:
    return sum(float(item["endpoint_area_m2"]) for item in active_items)


def active_edge_flux_avg(
    state: dict[str, list[float]],
    active_items: list[dict[str, Any]],
    mesh_edges: dict[int, dict[str, Any]],
    nodes: list[dict[str, Any]],
    contacts: dict[int, str],
    ni_model: list[float],
    vt: float,
    bias: float,
) -> float | None:
    numerator = 0.0
    denominator = 0.0
    for item in active_items:
        edge = mesh_edges.get(int(item["edge_id"]))
        if edge is None:
            continue
        area = float(item["endpoint_area_m2"])
        edge_row = fluxforms.edge_row(bias, "state", edge, nodes, contacts, state, ni_model, vt)
        particle = abs(float(edge_row["electron_flux_qf_model_abs"])) + abs(float(edge_row["hole_flux_qf_model_abs"]))
        numerator += particle * area
        denominator += area
    if denominator == 0.0:
        return None
    return numerator / denominator


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = consistency.optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {"count": len(rows)}
    for key in SUMMARY_KEYS:
        values = finite_values(rows, key)
        summary[f"{key}_median"] = statistics.median(values) if values else None
        summary[f"{key}_mean"] = statistics.fmean(values) if values else None
    return summary


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, contacts = fluxforms.sgdiag.read_mesh(args.mesh)
    edges = fluxforms.sgdiag.build_edges(nodes, triangles)
    mesh_edges = {int(edge["edge_id"]): edge for edge in edges}
    volumes = fluxforms.sgdiag.node_volumes(nodes, triangles)
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    ni_model = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    vela = fluxforms.load_vela_state(
        residual.fluxfac.discover_vtk(args.vela_vtk_root, args.bias),
        len(nodes),
    )
    sentaurus = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
    all_states = mixed.make_variants(
        vela,
        sentaurus,
        vt,
        args.electron_qf_shift_v,
        args.hole_qf_shift_v,
    )
    states = {name: all_states[source_name] for name, source_name in VARIANT_MAP.items()}
    transports = {
        name: residual.transport_terms(state, edges, vt, args.transport_model, ni_model)
        for name, state in states.items()
    }
    active_by_node = residual.load_active_edges_with_source(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    sent_nodes = residual.cont.load_sentaurus_nodes(args.sentaurus_dir)
    sentaurus_impact = fluxforms.read_scalar_csv(args.sentaurus_dir / "fields" / "ImpactIonization_region0.csv")
    sent_e_current = consistency.load_field_magnitude(args.sentaurus_dir, "eCurrentDensity")
    sent_h_current = consistency.load_field_magnitude(args.sentaurus_dir, "hCurrentDensity")
    vela_vtk = residual.fluxfac.discover_vtk(args.vela_vtk_root, args.bias)
    vela_scalars = fluxforms.sgdiag.parse_vtk_scalars(vela_vtk)
    vela_srh = vela_scalars.get("SRHRecombination", [0.0 for _ in nodes])

    rows: list[dict[str, Any]] = []
    for support in residual.support_nodes(args.support_csv):
        node_id = int(support["node_id"])
        node = nodes[node_id]
        sent_id, sent_distance = residual.cont.nearest_node(sent_nodes, node["x_um"], node["y_um"])
        active_items = active_by_node.get(node_id, [])
        endpoint_area = active_endpoint_area(active_items)
        support_e_flux = consistency.particle_flux_from_current(sent_e_current.get(node_id, 0.0))
        support_h_flux = consistency.particle_flux_from_current(sent_h_current.get(node_id, 0.0))
        support_particle_flux = support_e_flux + support_h_flux
        sentaurus_current_integral = support_particle_flux * endpoint_area
        volume = volumes[node_id]
        cxx_source = residual.cxx_edge_source(active_items)
        sentaurus_source = sentaurus_impact[sent_id] * 1.0e6 * volume
        srh_source = vela_srh[node_id] * volume
        for variant, state in states.items():
            electron_transport, hole_transport = transports[variant]
            electron_t = electron_transport[node_id]
            hole_t = hole_transport[node_id]
            replayed_source = residual.replayed_edge_source(
                state,
                active_items,
                mesh_edges,
                vt,
                args.transport_model,
                ni_model,
            )
            particle_flux = active_edge_flux_avg(
                state, active_items, mesh_edges, nodes, contacts, ni_model, vt, args.bias)
            electron_residual_cxx = electron_t + srh_source - cxx_source
            hole_residual_cxx = hole_t + srh_source - cxx_source
            electron_residual_replayed = electron_t + srh_source - replayed_source
            hole_residual_replayed = hole_t + srh_source - replayed_source
            electron_residual_sentaurus = electron_t + srh_source - sentaurus_source
            hole_residual_sentaurus = hole_t + srh_source - sentaurus_source
            rows.append({
                "bias_V": args.bias,
                "variant": variant,
                "node_id": node_id,
                "x_um": node["x_um"],
                "y_um": node["y_um"],
                "support_class": support.get("support_class", ""),
                "sentaurus_node_id": sent_id,
                "sentaurus_distance_um": sent_distance,
                "active_edge_count": len(active_items),
                "active_edge_ids": residual.active_edge_ids(active_by_node, node_id),
                "active_endpoint_area_sum_m2": endpoint_area,
                "node_volume_m2": volume,
                "transport_model": args.transport_model,
                "sentaurus_support_particle_flux_m2_s": support_particle_flux,
                "sentaurus_active_current_integral_s_inv": sentaurus_current_integral,
                "particle_flux_avg_m2_s": particle_flux,
                "particle_flux_over_sentaurus_support_current": consistency.ratio(
                    particle_flux, support_particle_flux),
                "electron_transport_s_inv": electron_t,
                "hole_transport_s_inv": hole_t,
                "cxx_edge_source_s_inv": cxx_source,
                "replayed_edge_source_s_inv": replayed_source,
                "sentaurus_node_source_s_inv": sentaurus_source,
                "srh_source_s_inv": srh_source,
                "cxx_edge_source_over_sentaurus_node_source": consistency.ratio(cxx_source, sentaurus_source),
                "replayed_edge_source_over_sentaurus_node_source": consistency.ratio(
                    replayed_source, sentaurus_source),
                "cxx_edge_source_over_sentaurus_active_current_integral": consistency.ratio(
                    cxx_source, sentaurus_current_integral),
                "replayed_edge_source_over_sentaurus_active_current_integral": consistency.ratio(
                    replayed_source, sentaurus_current_integral),
                "electron_transport_over_sentaurus_node_source": consistency.ratio(electron_t, sentaurus_source),
                "hole_transport_over_sentaurus_node_source": consistency.ratio(hole_t, sentaurus_source),
                "electron_transport_over_cxx_edge_source": consistency.ratio(electron_t, cxx_source),
                "hole_transport_over_cxx_edge_source": consistency.ratio(hole_t, cxx_source),
                "electron_residual_cxx_source_s_inv": electron_residual_cxx,
                "hole_residual_cxx_source_s_inv": hole_residual_cxx,
                "electron_residual_replayed_source_s_inv": electron_residual_replayed,
                "hole_residual_replayed_source_s_inv": hole_residual_replayed,
                "electron_residual_sentaurus_source_s_inv": electron_residual_sentaurus,
                "hole_residual_sentaurus_source_s_inv": hole_residual_sentaurus,
                "electron_residual_over_impact": consistency.ratio(electron_residual_cxx, cxx_source),
                "hole_residual_over_impact": consistency.ratio(hole_residual_cxx, cxx_source),
                "electron_residual_replayed_over_replayed_impact": consistency.ratio(
                    electron_residual_replayed, replayed_source),
                "hole_residual_replayed_over_replayed_impact": consistency.ratio(
                    hole_residual_replayed, replayed_source),
                "electron_residual_sentaurus_over_sentaurus_impact": consistency.ratio(
                    electron_residual_sentaurus, sentaurus_source),
                "hole_residual_sentaurus_over_sentaurus_impact": consistency.ratio(
                    hole_residual_sentaurus, sentaurus_source),
            })
    return rows


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for support_class in sorted({str(row["support_class"]) for row in rows}):
        class_rows = [row for row in rows if str(row["support_class"]) == support_class]
        variants: dict[str, Any] = {}
        for variant in sorted({str(row["variant"]) for row in class_rows}):
            variants[variant] = summarize_rows([
                row for row in class_rows if str(row["variant"]) == variant
            ])
        result[support_class] = {"variants": variants}
    return result


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


def main() -> int:
    args = parse_args()
    rows = make_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "qf_shift_continuity_balance_nodes.csv", rows)
    payload = clean_json({
        "bias_V": args.bias,
        "row_count": len(rows),
        "transport_model": args.transport_model,
        "electron_qf_shift_v": args.electron_qf_shift_v,
        "hole_qf_shift_v": args.hole_qf_shift_v,
        "support_classes": class_summary(rows),
    })
    (args.out_dir / "qf_shift_continuity_balance_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
