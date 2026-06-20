#!/usr/bin/env python3
"""Replay active-edge state substitutions against Sentaurus exported current."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_active_edge_mixed_state_replay as mixed
import diagnose_pn2d_bv_edge_source_current_consistency as consistency
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


CSV_FIELDS = [
    "bias_V",
    "variant",
    "support_node_id",
    "support_class",
    "x_um",
    "y_um",
    "edge_id",
    "edge_node0",
    "edge_node1",
    "edge_axis",
    "sentaurus_support_electron_current_flux_m2_s",
    "sentaurus_support_hole_current_flux_m2_s",
    "sentaurus_support_particle_flux_m2_s",
    "sentaurus_edgeavg_electron_current_flux_m2_s",
    "sentaurus_edgeavg_hole_current_flux_m2_s",
    "sentaurus_edgeavg_particle_flux_m2_s",
    "density_electron_flux_m2_s",
    "density_hole_flux_m2_s",
    "density_particle_flux_m2_s",
    "qf_old_slotboom_electron_flux_m2_s",
    "qf_old_slotboom_hole_flux_m2_s",
    "qf_old_slotboom_particle_flux_m2_s",
    "density_particle_flux_over_sentaurus_support_current",
    "qf_old_slotboom_particle_flux_over_sentaurus_support_current",
    "density_particle_flux_over_sentaurus_edgeavg_current",
    "qf_old_slotboom_particle_flux_over_sentaurus_edgeavg_current",
]

SUMMARY_METRICS = [
    "density_particle_flux_over_sentaurus_support_current",
    "qf_old_slotboom_particle_flux_over_sentaurus_support_current",
    "density_particle_flux_over_sentaurus_edgeavg_current",
    "qf_old_slotboom_particle_flux_over_sentaurus_edgeavg_current",
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
    parser.add_argument("--electron-qf-shift-v", type=float, default=0.0)
    parser.add_argument("--hole-qf-shift-v", type=float, default=0.0)
    return parser.parse_args()


def avg2(a: float, b: float) -> float:
    return 0.5 * (a + b)


def flux_bundle(edge_row: dict[str, Any]) -> dict[str, dict[str, float]]:
    density_e = abs(float(edge_row["electron_flux_density_abs"]))
    density_h = abs(float(edge_row["hole_flux_density_abs"]))
    qf_e = abs(float(edge_row["electron_flux_qf_model_abs"]))
    qf_h = abs(float(edge_row["hole_flux_qf_model_abs"]))
    return {
        "density": {
            "electron": density_e,
            "hole": density_h,
            "particle": density_e + density_h,
        },
        "qf_old_slotboom": {
            "electron": qf_e,
            "hole": qf_h,
            "particle": qf_e + qf_h,
        },
    }


def make_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    nodes, triangles, contacts = fluxforms.sgdiag.read_mesh(args.mesh)
    mesh_edges = {int(edge["edge_id"]): edge for edge in fluxforms.sgdiag.build_edges(nodes, triangles)}
    vt = fluxforms.sgdiag.K_B_OVER_Q * args.temperature_k
    ni_model = fluxforms.effective_ni_from_doping(
        args.doping_csv,
        len(nodes),
        args.material_ni_m3,
        vt,
        args.bandgap_narrowing,
    )
    vela = fluxforms.load_vela_state(fluxfac.discover_vtk(args.vela_vtk_root, args.bias), len(nodes))
    sentaurus = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
    variants = mixed.make_variants(vela, sentaurus, vt, args.electron_qf_shift_v, args.hole_qf_shift_v)
    active_by_node = consistency.load_active_endpoint_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    sent_e_current = consistency.load_field_magnitude(args.sentaurus_dir, "eCurrentDensity")
    sent_h_current = consistency.load_field_magnitude(args.sentaurus_dir, "hCurrentDensity")

    rows: list[dict[str, Any]] = []
    for support in consistency.support_rows(args.support_csv):
        support_node = int(support["node_id"])
        support_e = consistency.particle_flux_from_current(sent_e_current.get(support_node, 0.0))
        support_h = consistency.particle_flux_from_current(sent_h_current.get(support_node, 0.0))
        support_particle = support_e + support_h
        for item in active_by_node.get(support_node, []):
            edge_id = int(item["edge_id"])
            edge = mesh_edges[edge_id]
            i = int(edge["node0"])
            j = int(edge["node1"])
            edgeavg_e = avg2(
                consistency.particle_flux_from_current(sent_e_current.get(i, 0.0)),
                consistency.particle_flux_from_current(sent_e_current.get(j, 0.0)),
            )
            edgeavg_h = avg2(
                consistency.particle_flux_from_current(sent_h_current.get(i, 0.0)),
                consistency.particle_flux_from_current(sent_h_current.get(j, 0.0)),
            )
            edgeavg_particle = edgeavg_e + edgeavg_h
            for name, state in variants.items():
                edge_flux = flux_bundle(
                    fluxforms.edge_row(args.bias, name, edge, nodes, contacts, state, ni_model, vt)
                )
                density = edge_flux["density"]
                qf = edge_flux["qf_old_slotboom"]
                rows.append({
                    "bias_V": args.bias,
                    "variant": name,
                    "support_node_id": support_node,
                    "support_class": support.get("support_class", ""),
                    "x_um": support.get("x_um", ""),
                    "y_um": support.get("y_um", ""),
                    "edge_id": edge_id,
                    "edge_node0": i,
                    "edge_node1": j,
                    "edge_axis": item["axis"],
                    "sentaurus_support_electron_current_flux_m2_s": support_e,
                    "sentaurus_support_hole_current_flux_m2_s": support_h,
                    "sentaurus_support_particle_flux_m2_s": support_particle,
                    "sentaurus_edgeavg_electron_current_flux_m2_s": edgeavg_e,
                    "sentaurus_edgeavg_hole_current_flux_m2_s": edgeavg_h,
                    "sentaurus_edgeavg_particle_flux_m2_s": edgeavg_particle,
                    "density_electron_flux_m2_s": density["electron"],
                    "density_hole_flux_m2_s": density["hole"],
                    "density_particle_flux_m2_s": density["particle"],
                    "qf_old_slotboom_electron_flux_m2_s": qf["electron"],
                    "qf_old_slotboom_hole_flux_m2_s": qf["hole"],
                    "qf_old_slotboom_particle_flux_m2_s": qf["particle"],
                    "density_particle_flux_over_sentaurus_support_current": consistency.ratio(
                        density["particle"], support_particle),
                    "qf_old_slotboom_particle_flux_over_sentaurus_support_current": consistency.ratio(
                        qf["particle"], support_particle),
                    "density_particle_flux_over_sentaurus_edgeavg_current": consistency.ratio(
                        density["particle"], edgeavg_particle),
                    "qf_old_slotboom_particle_flux_over_sentaurus_edgeavg_current": consistency.ratio(
                        qf["particle"], edgeavg_particle),
                })
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = consistency.optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize_variant(rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "count": len(rows),
        "active_edge_ids": sorted({int(row["edge_id"]) for row in rows}),
    }
    for key in SUMMARY_METRICS:
        values = finite_values(rows, key)
        summary[f"{key}_median"] = statistics.median(values) if values else None
        summary[f"{key}_mean"] = statistics.fmean(values) if values else None
    return summary


def class_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for support_class in sorted({str(row["support_class"]) for row in rows}):
        class_rows = [row for row in rows if str(row["support_class"]) == support_class]
        variants: dict[str, Any] = {}
        for variant in sorted({str(row["variant"]) for row in class_rows}):
            variants[variant] = summarize_variant([
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
    write_csv(args.out_dir / "edge_state_substitution_current_edges.csv", rows)
    summary = clean_json({
        "bias_V": args.bias,
        "row_count": len(rows),
        "electron_qf_shift_v": args.electron_qf_shift_v,
        "hole_qf_shift_v": args.hole_qf_shift_v,
        "support_classes": class_summary(rows),
    })
    (args.out_dir / "edge_state_substitution_current_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
