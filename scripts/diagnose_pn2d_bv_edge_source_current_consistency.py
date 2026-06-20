#!/usr/bin/env python3
"""Compare active-edge source proxies with density/qF/Sentaurus current fluxes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


Q_C = 1.602176634e-19

CSV_FIELDS = [
    "bias_V",
    "support_node_id",
    "support_class",
    "x_um",
    "y_um",
    "edge_id",
    "edge_node0",
    "edge_node1",
    "edge_axis",
    "endpoint_area_m2",
    "edge_area_m2",
    "cxx_electron_flux_m2_s",
    "cxx_hole_flux_m2_s",
    "cxx_particle_flux_m2_s",
    "cxx_generation_density_m3_s",
    "cxx_endpoint_source_integral_s_inv",
    "cxx_edge_source_integral_s_inv",
    "sentaurus_support_electron_current_flux_m2_s",
    "sentaurus_support_hole_current_flux_m2_s",
    "sentaurus_support_particle_flux_m2_s",
    "sentaurus_edgeavg_electron_current_flux_m2_s",
    "sentaurus_edgeavg_hole_current_flux_m2_s",
    "sentaurus_edgeavg_particle_flux_m2_s",
    "sentaurus_support_generation_m3_s",
    "sentaurus_edgeavg_generation_m3_s",
    "vela_density_electron_flux_m2_s",
    "vela_density_hole_flux_m2_s",
    "vela_density_particle_flux_m2_s",
    "vela_qf_inferred_electron_flux_m2_s",
    "vela_qf_inferred_hole_flux_m2_s",
    "vela_qf_inferred_particle_flux_m2_s",
    "vela_qf_old_slotboom_electron_flux_m2_s",
    "vela_qf_old_slotboom_hole_flux_m2_s",
    "vela_qf_old_slotboom_particle_flux_m2_s",
    "sentaurus_density_electron_flux_m2_s",
    "sentaurus_density_hole_flux_m2_s",
    "sentaurus_density_particle_flux_m2_s",
    "sentaurus_qf_inferred_electron_flux_m2_s",
    "sentaurus_qf_inferred_hole_flux_m2_s",
    "sentaurus_qf_inferred_particle_flux_m2_s",
    "sentaurus_qf_old_slotboom_electron_flux_m2_s",
    "sentaurus_qf_old_slotboom_hole_flux_m2_s",
    "sentaurus_qf_old_slotboom_particle_flux_m2_s",
    "cxx_particle_flux_over_sentaurus_support_current",
    "cxx_particle_flux_over_sentaurus_edgeavg_current",
    "cxx_generation_over_sentaurus_support_generation",
    "cxx_generation_over_sentaurus_edgeavg_generation",
    "cxx_particle_flux_over_vela_density_flux",
    "cxx_particle_flux_over_vela_qf_old_slotboom_flux",
    "vela_density_flux_over_sentaurus_support_current",
    "vela_qf_old_slotboom_flux_over_sentaurus_support_current",
    "sentaurus_density_flux_over_sentaurus_support_current",
    "sentaurus_qf_old_slotboom_flux_over_sentaurus_support_current",
]

SUMMARY_METRICS = [
    "cxx_particle_flux_over_sentaurus_support_current",
    "cxx_particle_flux_over_sentaurus_edgeavg_current",
    "cxx_generation_over_sentaurus_support_generation",
    "cxx_generation_over_sentaurus_edgeavg_generation",
    "cxx_particle_flux_over_vela_density_flux",
    "cxx_particle_flux_over_vela_qf_old_slotboom_flux",
    "vela_density_flux_over_sentaurus_support_current",
    "vela_qf_old_slotboom_flux_over_sentaurus_support_current",
    "sentaurus_density_flux_over_sentaurus_support_current",
    "sentaurus_qf_old_slotboom_flux_over_sentaurus_support_current",
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


def support_rows(path: Path) -> list[dict[str, str]]:
    return [
        row for row in read_csv_rows(path)
        if row.get("support_class") not in (None, "", "inactive")
    ]


def numeric_components(row: dict[str, str]) -> list[float]:
    values: list[float] = []
    for key, raw in row.items():
        if key == "node_id":
            continue
        value = optional_float(raw)
        if value is not None:
            values.append(value)
    return values


def load_field_magnitude(sentaurus_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted((sentaurus_dir / "fields").glob(f"{name}_region*.csv")):
        for row in read_csv_rows(path):
            values = numeric_components(row)
            result[int(row["node_id"])] = math.sqrt(sum(value * value for value in values)) if values else 0.0
    return result


def load_scalar_field(sentaurus_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted((sentaurus_dir / "fields").glob(f"{name}_region*.csv")):
        for row in read_csv_rows(path):
            values = numeric_components(row)
            result[int(row["node_id"])] = values[0] if values else 0.0
    return result


def particle_flux_from_current(current_a_cm2: float) -> float:
    return abs(current_a_cm2) * 1.0e4 / Q_C


def edge_axis(row: dict[str, str]) -> str:
    dx = abs((optional_float(row.get("x1_um")) or 0.0) - (optional_float(row.get("x0_um")) or 0.0))
    dy = abs((optional_float(row.get("y1_um")) or 0.0) - (optional_float(row.get("y0_um")) or 0.0))
    return "x" if dx >= dy else "y"


def edge_source_integral(row: dict[str, str], alpha_flux: float, edge_area: float) -> float:
    explicit = optional_float(row.get("edge_source_integral"))
    if explicit is not None:
        return explicit
    electron = optional_float(row.get("electron_source_integral")) or 0.0
    hole = optional_float(row.get("hole_source_integral")) or 0.0
    if electron != 0.0 or hole != 0.0:
        return electron + hole
    return alpha_flux * edge_area


def load_active_endpoint_edges(
    edge_csv: Path,
    bias: float,
    active_axis: str,
    threshold_scale: float,
) -> dict[int, list[dict[str, Any]]]:
    by_node: dict[int, list[dict[str, Any]]] = {}
    for row in read_csv_rows(edge_csv):
        if not fluxfac.matches_bias(row, bias):
            continue
        electron_flux = abs(optional_float(row.get("electron_flux_proxy")) or 0.0)
        hole_flux = abs(optional_float(row.get("hole_flux_proxy")) or 0.0)
        electron_alpha = optional_float(row.get("electron_alpha_m_inv")) or 0.0
        hole_alpha = optional_float(row.get("hole_alpha_m_inv")) or 0.0
        alpha_flux = electron_alpha * electron_flux + hole_alpha * hole_flux
        edge_area = optional_float(row.get("edge_area_proxy_m2")) or 0.0
        edge_source = edge_source_integral(row, alpha_flux, edge_area)
        item = {
            "edge_id": int(row["edge_id"]),
            "axis": edge_axis(row),
            "endpoint_area_m2": 0.5 * edge_area,
            "edge_area_m2": edge_area,
            "electron_flux_m2_s": electron_flux,
            "hole_flux_m2_s": hole_flux,
            "generation_density_m3_s": alpha_flux,
            "edge_source_integral_s_inv": edge_source,
            "endpoint_source_integral_s_inv": 0.5 * edge_source,
        }
        for endpoint in ["node0", "node1"]:
            if row.get(endpoint) not in (None, ""):
                by_node.setdefault(int(row[endpoint]), []).append(item)
    active: dict[int, list[dict[str, Any]]] = {}
    for node_id, items in by_node.items():
        max_alpha_flux = max((abs(float(item["generation_density_m3_s"])) for item in items), default=0.0)
        threshold = max_alpha_flux * threshold_scale
        active[node_id] = [
            item for item in items
            if item["axis"] == active_axis and abs(float(item["generation_density_m3_s"])) > threshold
        ]
    return active


def flux_triplet(edge_row: dict[str, Any], prefix: str) -> dict[str, float]:
    electron = abs(float(edge_row[f"{prefix}_electron_flux_m2_s"]))
    hole = abs(float(edge_row[f"{prefix}_hole_flux_m2_s"]))
    return {"electron": electron, "hole": hole, "particle": electron + hole}


def state_fluxes(edge_row: dict[str, Any]) -> dict[str, dict[str, float]]:
    density_e = abs(float(edge_row["electron_flux_density_abs"]))
    density_h = abs(float(edge_row["hole_flux_density_abs"]))
    qf_inferred_e = abs(float(edge_row["electron_flux_qf_inferred_abs"]))
    qf_inferred_h = abs(float(edge_row["hole_flux_qf_inferred_abs"]))
    qf_model_e = abs(float(edge_row["electron_flux_qf_model_abs"]))
    qf_model_h = abs(float(edge_row["hole_flux_qf_model_abs"]))
    return {
        "density": {"electron": density_e, "hole": density_h, "particle": density_e + density_h},
        "qf_inferred": {"electron": qf_inferred_e, "hole": qf_inferred_h, "particle": qf_inferred_e + qf_inferred_h},
        "qf_old_slotboom": {"electron": qf_model_e, "hole": qf_model_h, "particle": qf_model_e + qf_model_h},
    }


def avg2(a: float, b: float) -> float:
    return 0.5 * (a + b)


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
    vela_state = fluxforms.load_vela_state(fluxfac.discover_vtk(args.vela_vtk_root, args.bias), len(nodes))
    sentaurus_state = fluxforms.load_sentaurus_state(args.sentaurus_dir, len(nodes))
    active_by_node = load_active_endpoint_edges(
        args.sg_edge_csv,
        args.bias,
        args.active_axis,
        args.active_source_relative_threshold,
    )
    sent_e_current = load_field_magnitude(args.sentaurus_dir, "eCurrentDensity")
    sent_h_current = load_field_magnitude(args.sentaurus_dir, "hCurrentDensity")
    sent_impact = load_scalar_field(args.sentaurus_dir, "ImpactIonization")

    needed_edges = sorted({
        int(item["edge_id"])
        for support in support_rows(args.support_csv)
        for item in active_by_node.get(int(support["node_id"]), [])
    })
    vela_rows = {
        edge_id: fluxforms.edge_row(args.bias, "vela", mesh_edges[edge_id], nodes, contacts, vela_state, ni_model, vt)
        for edge_id in needed_edges
    }
    sent_rows = {
        edge_id: fluxforms.edge_row(
            args.bias, "sentaurus", mesh_edges[edge_id], nodes, contacts, sentaurus_state, ni_model, vt)
        for edge_id in needed_edges
    }

    rows: list[dict[str, Any]] = []
    for support in support_rows(args.support_csv):
        support_node = int(support["node_id"])
        support_particle = (
            particle_flux_from_current(sent_e_current.get(support_node, 0.0))
            + particle_flux_from_current(sent_h_current.get(support_node, 0.0))
        )
        support_generation = sent_impact.get(support_node, 0.0) * 1.0e6
        for item in active_by_node.get(support_node, []):
            edge_id = int(item["edge_id"])
            edge = mesh_edges[edge_id]
            i = int(edge["node0"])
            j = int(edge["node1"])
            edgeavg_e = avg2(
                particle_flux_from_current(sent_e_current.get(i, 0.0)),
                particle_flux_from_current(sent_e_current.get(j, 0.0)),
            )
            edgeavg_h = avg2(
                particle_flux_from_current(sent_h_current.get(i, 0.0)),
                particle_flux_from_current(sent_h_current.get(j, 0.0)),
            )
            edgeavg_particle = edgeavg_e + edgeavg_h
            edgeavg_generation = avg2(sent_impact.get(i, 0.0), sent_impact.get(j, 0.0)) * 1.0e6
            vela_flux = state_fluxes(vela_rows[edge_id])
            sent_flux = state_fluxes(sent_rows[edge_id])
            cxx_e = float(item["electron_flux_m2_s"])
            cxx_h = float(item["hole_flux_m2_s"])
            cxx_particle = cxx_e + cxx_h
            row = {
                "bias_V": args.bias,
                "support_node_id": support_node,
                "support_class": support.get("support_class", ""),
                "x_um": support.get("x_um", ""),
                "y_um": support.get("y_um", ""),
                "edge_id": edge_id,
                "edge_node0": i,
                "edge_node1": j,
                "edge_axis": item["axis"],
                "endpoint_area_m2": item["endpoint_area_m2"],
                "edge_area_m2": item["edge_area_m2"],
                "cxx_electron_flux_m2_s": cxx_e,
                "cxx_hole_flux_m2_s": cxx_h,
                "cxx_particle_flux_m2_s": cxx_particle,
                "cxx_generation_density_m3_s": item["generation_density_m3_s"],
                "cxx_endpoint_source_integral_s_inv": item["endpoint_source_integral_s_inv"],
                "cxx_edge_source_integral_s_inv": item["edge_source_integral_s_inv"],
                "sentaurus_support_electron_current_flux_m2_s": particle_flux_from_current(
                    sent_e_current.get(support_node, 0.0)),
                "sentaurus_support_hole_current_flux_m2_s": particle_flux_from_current(
                    sent_h_current.get(support_node, 0.0)),
                "sentaurus_support_particle_flux_m2_s": support_particle,
                "sentaurus_edgeavg_electron_current_flux_m2_s": edgeavg_e,
                "sentaurus_edgeavg_hole_current_flux_m2_s": edgeavg_h,
                "sentaurus_edgeavg_particle_flux_m2_s": edgeavg_particle,
                "sentaurus_support_generation_m3_s": support_generation,
                "sentaurus_edgeavg_generation_m3_s": edgeavg_generation,
                "vela_density_electron_flux_m2_s": vela_flux["density"]["electron"],
                "vela_density_hole_flux_m2_s": vela_flux["density"]["hole"],
                "vela_density_particle_flux_m2_s": vela_flux["density"]["particle"],
                "vela_qf_inferred_electron_flux_m2_s": vela_flux["qf_inferred"]["electron"],
                "vela_qf_inferred_hole_flux_m2_s": vela_flux["qf_inferred"]["hole"],
                "vela_qf_inferred_particle_flux_m2_s": vela_flux["qf_inferred"]["particle"],
                "vela_qf_old_slotboom_electron_flux_m2_s": vela_flux["qf_old_slotboom"]["electron"],
                "vela_qf_old_slotboom_hole_flux_m2_s": vela_flux["qf_old_slotboom"]["hole"],
                "vela_qf_old_slotboom_particle_flux_m2_s": vela_flux["qf_old_slotboom"]["particle"],
                "sentaurus_density_electron_flux_m2_s": sent_flux["density"]["electron"],
                "sentaurus_density_hole_flux_m2_s": sent_flux["density"]["hole"],
                "sentaurus_density_particle_flux_m2_s": sent_flux["density"]["particle"],
                "sentaurus_qf_inferred_electron_flux_m2_s": sent_flux["qf_inferred"]["electron"],
                "sentaurus_qf_inferred_hole_flux_m2_s": sent_flux["qf_inferred"]["hole"],
                "sentaurus_qf_inferred_particle_flux_m2_s": sent_flux["qf_inferred"]["particle"],
                "sentaurus_qf_old_slotboom_electron_flux_m2_s": sent_flux["qf_old_slotboom"]["electron"],
                "sentaurus_qf_old_slotboom_hole_flux_m2_s": sent_flux["qf_old_slotboom"]["hole"],
                "sentaurus_qf_old_slotboom_particle_flux_m2_s": sent_flux["qf_old_slotboom"]["particle"],
            }
            row.update({
                "cxx_particle_flux_over_sentaurus_support_current": ratio(cxx_particle, support_particle),
                "cxx_particle_flux_over_sentaurus_edgeavg_current": ratio(cxx_particle, edgeavg_particle),
                "cxx_generation_over_sentaurus_support_generation": ratio(
                    item["generation_density_m3_s"], support_generation),
                "cxx_generation_over_sentaurus_edgeavg_generation": ratio(
                    item["generation_density_m3_s"], edgeavg_generation),
                "cxx_particle_flux_over_vela_density_flux": ratio(cxx_particle, vela_flux["density"]["particle"]),
                "cxx_particle_flux_over_vela_qf_old_slotboom_flux": ratio(
                    cxx_particle, vela_flux["qf_old_slotboom"]["particle"]),
                "vela_density_flux_over_sentaurus_support_current": ratio(
                    vela_flux["density"]["particle"], support_particle),
                "vela_qf_old_slotboom_flux_over_sentaurus_support_current": ratio(
                    vela_flux["qf_old_slotboom"]["particle"], support_particle),
                "sentaurus_density_flux_over_sentaurus_support_current": ratio(
                    sent_flux["density"]["particle"], support_particle),
                "sentaurus_qf_old_slotboom_flux_over_sentaurus_support_current": ratio(
                    sent_flux["qf_old_slotboom"]["particle"], support_particle),
            })
            rows.append(row)
    return rows


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
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
        result[support_class] = summarize_rows([
            row for row in rows if str(row["support_class"]) == support_class
        ])
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
    write_csv(args.out_dir / "edge_source_current_consistency_edges.csv", rows)
    summary = {
        "bias_V": args.bias,
        "row_count": len(rows),
        "support_classes": class_summary(rows),
    }
    payload = clean_json(summary)
    (args.out_dir / "edge_source_current_consistency_summary.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
