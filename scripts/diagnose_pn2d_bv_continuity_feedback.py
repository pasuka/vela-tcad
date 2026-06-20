#!/usr/bin/env python3
"""Localize PN2D BV carrier/continuity feedback near a selected SG edge."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag


Q = 1.602176634e-19


EDGE_FIELDS = [
    "bias_V",
    "focus_edge_id",
    "edge_id",
    "source_rank",
    "relation",
    "node0",
    "node1",
    "sentaurus_node0",
    "sentaurus_node1",
    "edge_mid_x_um",
    "edge_mid_y_um",
    "edge_class",
    "edge_area_m2",
    "vela_dpsi_V",
    "sentaurus_dpsi_V",
    "vela_dphin_V",
    "sentaurus_dphin_V",
    "vela_dphip_V",
    "sentaurus_dphip_V",
    "vela_electron_qf_field_V_m",
    "sentaurus_electron_qf_field_V_m",
    "vela_hole_qf_field_V_m",
    "sentaurus_hole_qf_field_V_m",
    "vela_electron_flux_abs_m2_s",
    "sentaurus_electron_flux_abs_m2_s",
    "vela_hole_flux_abs_m2_s",
    "sentaurus_hole_flux_abs_m2_s",
    "vela_source_density_m3_s",
    "sentaurus_generation_avg_m3_s",
    "log10_vela_over_sentaurus_generation",
    "vela_electron_density_avg_cm3",
    "sentaurus_electron_density_avg_cm3",
    "log10_vela_over_sentaurus_electron_density",
    "vela_hole_density_avg_cm3",
    "sentaurus_hole_density_avg_cm3",
    "log10_vela_over_sentaurus_hole_density",
]


NODE_FIELDS = [
    "bias_V",
    "focus_edge_id",
    "node_id",
    "relation",
    "x_um",
    "y_um",
    "sentaurus_node_id",
    "sentaurus_distance_um",
    "node_volume_m2",
    "vela_psi_V",
    "sentaurus_psi_V",
    "vela_phin_V",
    "sentaurus_phin_V",
    "vela_phip_V",
    "sentaurus_phip_V",
    "vela_psi_minus_phin_V",
    "sentaurus_psi_minus_phin_V",
    "delta_psi_minus_phin_V",
    "vela_phip_minus_psi_V",
    "sentaurus_phip_minus_psi_V",
    "delta_phip_minus_psi_V",
    "vela_electron_density_cm3",
    "sentaurus_electron_density_cm3",
    "log10_vela_over_sentaurus_electron_density",
    "vela_hole_density_cm3",
    "sentaurus_hole_density_cm3",
    "log10_vela_over_sentaurus_hole_density",
    "vela_ni_eff_cm3",
    "sentaurus_ni_eff_from_electron_cm3",
    "sentaurus_ni_eff_from_hole_cm3",
    "vela_electron_density_reconstructed_cm3",
    "vela_hole_density_reconstructed_cm3",
    "vela_electron_density_from_sentaurus_qf_cm3",
    "vela_hole_density_from_sentaurus_qf_cm3",
    "vela_electron_transport_integral_s_inv",
    "vela_hole_transport_integral_s_inv",
    "vela_avalanche_node_integral_s_inv",
    "vela_srh_node_integral_s_inv",
    "vela_electron_residual_estimate_s_inv",
    "vela_hole_residual_estimate_s_inv",
    "sentaurus_generation_node_integral_s_inv",
    "sentaurus_srh_node_integral_s_inv",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-10,-13.2")
    parser.add_argument("--edge-id", type=int, default=2886)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def bias_token(bias: float) -> str:
    return f"{bias:g}".replace("-", "m").replace("+", "").replace(".", "p")


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if match:
            result[bias_key(float(match.group("bias")))] = path
    return result


def sentaurus_dir(root: Path, bias: float) -> Path:
    candidates = [
        root / f"sentaurus_{signed_bias_token(bias)}v",
        root / f"sentaurus_{bias_token(bias)}v",
        root / f"sentaurus_{bias_token(abs(bias))}v",
        root / f"sentaurus_{bias:g}v",
    ]
    for candidate in candidates:
        if (candidate / "nodes.csv").exists():
            return candidate
    raise FileNotFoundError(f"missing Sentaurus export for {bias:g} V under {root}")


def load_sentaurus_nodes(root: Path) -> dict[int, dict[str, float]]:
    with (root / "nodes.csv").open(newline="") as handle:
        return {
            int(row["id"]): {
                "id": int(row["id"]),
                "x_um": float(row["x_um"]),
                "y_um": float(row["y_um"]),
            }
            for row in csv.DictReader(handle)
        }


def load_components(root: Path, name: str) -> dict[int, list[float]]:
    path = root / "fields" / f"{name}_region0.csv"
    if not path.exists():
        raise FileNotFoundError(path)
    values: dict[int, list[float]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        columns = [column for column in reader.fieldnames or [] if column != "node_id"]
        for row in reader:
            comps = [float(row[column]) for column in columns if row.get(column, "") != ""]
            values[int(row["node_id"])] = comps
    return values


def scalar_at(fields: dict[str, dict[int, list[float]]], name: str, node_id: int) -> float:
    return fields[name][node_id][0]


def magnitude_at(fields: dict[str, dict[int, list[float]]], name: str, node_id: int) -> float:
    values = fields[name][node_id]
    return math.sqrt(sum(value * value for value in values))


def nearest_node(nodes: dict[int, dict[str, float]], x_um: float, y_um: float) -> tuple[int, float]:
    best = min(
        nodes.values(),
        key=lambda item: (item["x_um"] - x_um) ** 2 + (item["y_um"] - y_um) ** 2,
    )
    distance_um = math.hypot(best["x_um"] - x_um, best["y_um"] - y_um)
    return int(best["id"]), distance_um


def log10_ratio(candidate: float, reference: float) -> float | None:
    if candidate == 0.0 or reference == 0.0:
        return None
    return math.log10(abs(candidate) / abs(reference))


def limited_exp(value: float) -> float:
    return math.exp(max(-500.0, min(500.0, value)))


def density_reconstruction(ni_m3: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if carrier == "electron":
        return ni_m3 * limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return ni_m3 * limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def inferred_ni_cm3(density_cm3: float, psi: float, qf: float, vt: float, carrier: str) -> float:
    if density_cm3 <= 0.0:
        return 0.0
    if carrier == "electron":
        return density_cm3 / limited_exp((psi - qf) / vt)
    if carrier == "hole":
        return density_cm3 / limited_exp((qf - psi) / vt)
    raise ValueError(carrier)


def selected_edges(edges: list[dict[str, Any]], focus_edge_id: int) -> list[int]:
    focus = next(edge for edge in edges if int(edge["edge_id"]) == focus_edge_id)
    focus_nodes = {int(focus["node0"]), int(focus["node1"])}
    ids = [
        int(edge["edge_id"])
        for edge in edges
        if int(edge["node0"]) in focus_nodes or int(edge["node1"]) in focus_nodes
    ]
    return sorted(set(ids))


def compute_signed_terms(rows: list[dict[str, Any]]) -> tuple[dict[int, float], dict[int, float]]:
    electron: dict[int, float] = defaultdict(float)
    hole: dict[int, float] = defaultdict(float)
    for row in rows:
        i = int(row["node0"])
        j = int(row["node1"])
        couple = float(row["couple_m"])
        electron_flux = float(row["electron_flux_signed_m2_s"])
        hole_flux = float(row["hole_flux_signed_m2_s"])
        electron[i] += electron_flux * couple
        electron[j] -= electron_flux * couple
        hole[i] += hole_flux * couple
        hole[j] -= hole_flux * couple
    return electron, hole


def source_rows_with_signed_flux(
    nodes: dict[int, dict[str, float]],
    edges: list[dict[str, Any]],
    contacts: dict[int, set[str]],
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
) -> tuple[list[dict[str, Any]], list[float]]:
    rows, node_source = sgdiag.source_rows(
        nodes, edges, contacts, scalars, vt, "quasi_fermi_gradient",
        "quasi_fermi_variable_ni", ni)
    by_edge = {int(row["edge_id"]): row for row in rows}
    psi = scalars["Potential"]
    phin = scalars["ElectronQuasiFermi"]
    phip = scalars["HoleQuasiFermi"]
    for edge in edges:
        edge_id = int(edge["edge_id"])
        row = by_edge.get(edge_id)
        if row is None:
            continue
        i = int(edge["node0"])
        j = int(edge["node1"])
        h = float(edge["length_m"])
        if h <= 1.0e-30:
            continue
        mun = float(row["electron_mobility_m2_V_s"])
        mup = float(row["hole_mobility_m2_V_s"])
        row["electron_flux_signed_m2_s"] = sgdiag.sg_electron_flux_qf_variable_ni(
            ni[i], ni[j], psi[i], psi[j], phin[i], phin[j], vt, mun * vt / h)
        row["hole_flux_signed_m2_s"] = sgdiag.sg_hole_flux_qf_variable_ni(
            ni[i], ni[j], psi[i], psi[j], phip[i], phip[j], vt, mup * vt / h)
    return rows, node_source


def load_sentaurus_fields(root: Path) -> dict[str, dict[int, list[float]]]:
    names = [
        "ElectrostaticPotential",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eDensity",
        "hDensity",
        "eCurrentDensity",
        "hCurrentDensity",
        "ImpactIonization",
        "srhRecombination",
    ]
    return {name: load_components(root, name) for name in names}


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def append_bias_rows(
    bias: float,
    focus_edge_id: int,
    mesh: Path,
    doping_csv: Path,
    vtk: Path,
    sentaurus: Path,
    vt: float,
    edge_rows: list[dict[str, Any]],
    node_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    nodes, triangles, contacts = sgdiag.read_mesh(mesh)
    edges = sgdiag.build_edges(nodes, triangles)
    scalars = sgdiag.parse_vtk_scalars(vtk)
    ni = sgdiag.effective_ni_from_doping(doping_csv, len(nodes), 1.0e16, vt, "old_slotboom")
    rows, node_avalanche = source_rows_with_signed_flux(nodes, edges, contacts, scalars, ni, vt)
    source_by_edge = {int(row["edge_id"]): row for row in rows}
    selected = [
        edge_id for edge_id in selected_edges(edges, focus_edge_id)
        if edge_id in source_by_edge
    ]
    selected_nodes = sorted({
        int(edge["node0"])
        for edge in edges
        if int(edge["edge_id"]) in selected
    } | {
        int(edge["node1"])
        for edge in edges
        if int(edge["edge_id"]) in selected
    })
    focus = next(edge for edge in edges if int(edge["edge_id"]) == focus_edge_id)
    focus_nodes = {int(focus["node0"]), int(focus["node1"])}
    volumes = sgdiag.node_volumes(nodes, triangles)
    sent_nodes = load_sentaurus_nodes(sentaurus)
    sent_fields = load_sentaurus_fields(sentaurus)
    electron_transport, hole_transport = compute_signed_terms(rows)
    srh = scalars.get("SRHRecombination", [0.0 for _ in nodes])

    for edge_id in selected:
        edge = next(item for item in edges if int(item["edge_id"]) == edge_id)
        row = source_by_edge[edge_id]
        i = int(edge["node0"])
        j = int(edge["node1"])
        si, _ = nearest_node(sent_nodes, nodes[i]["x_um"], nodes[i]["y_um"])
        sj, _ = nearest_node(sent_nodes, nodes[j]["x_um"], nodes[j]["y_um"])
        edge_area = float(row["edge_area_m2"])
        sent_generation = 0.5 * (
            scalar_at(sent_fields, "ImpactIonization", si) +
            scalar_at(sent_fields, "ImpactIonization", sj)) * 1.0e6
        sent_e_current = 0.5 * (
            magnitude_at(sent_fields, "eCurrentDensity", si) +
            magnitude_at(sent_fields, "eCurrentDensity", sj))
        sent_h_current = 0.5 * (
            magnitude_at(sent_fields, "hCurrentDensity", si) +
            magnitude_at(sent_fields, "hCurrentDensity", sj))
        n_vela = 0.5 * (scalars["Electrons"][i] + scalars["Electrons"][j]) * 1.0e-6
        p_vela = 0.5 * (scalars["Holes"][i] + scalars["Holes"][j]) * 1.0e-6
        n_sent = 0.5 * (
            scalar_at(sent_fields, "eDensity", si) +
            scalar_at(sent_fields, "eDensity", sj))
        p_sent = 0.5 * (
            scalar_at(sent_fields, "hDensity", si) +
            scalar_at(sent_fields, "hDensity", sj))
        relation = "focus" if edge_id == focus_edge_id else "incident_to_focus"
        edge_rows.append({
            "bias_V": bias,
            "focus_edge_id": focus_edge_id,
            "edge_id": edge_id,
            "source_rank": row["rank"],
            "relation": relation,
            "node0": i,
            "node1": j,
            "sentaurus_node0": si,
            "sentaurus_node1": sj,
            "edge_mid_x_um": 0.5 * (float(row["x0_um"]) + float(row["x1_um"])),
            "edge_mid_y_um": 0.5 * (float(row["y0_um"]) + float(row["y1_um"])),
            "edge_class": row["edge_class"],
            "edge_area_m2": edge_area,
            "vela_dpsi_V": scalars["Potential"][j] - scalars["Potential"][i],
            "sentaurus_dpsi_V": (
                scalar_at(sent_fields, "ElectrostaticPotential", sj) -
                scalar_at(sent_fields, "ElectrostaticPotential", si)),
            "vela_dphin_V": scalars["ElectronQuasiFermi"][j] - scalars["ElectronQuasiFermi"][i],
            "sentaurus_dphin_V": (
                scalar_at(sent_fields, "eQuasiFermiPotential", sj) -
                scalar_at(sent_fields, "eQuasiFermiPotential", si)),
            "vela_dphip_V": scalars["HoleQuasiFermi"][j] - scalars["HoleQuasiFermi"][i],
            "sentaurus_dphip_V": (
                scalar_at(sent_fields, "hQuasiFermiPotential", sj) -
                scalar_at(sent_fields, "hQuasiFermiPotential", si)),
            "vela_electron_qf_field_V_m": row["electron_field_V_m"],
            "sentaurus_electron_qf_field_V_m": abs(
                scalar_at(sent_fields, "eQuasiFermiPotential", sj) -
                scalar_at(sent_fields, "eQuasiFermiPotential", si)) / float(edge["length_m"]),
            "vela_hole_qf_field_V_m": row["hole_field_V_m"],
            "sentaurus_hole_qf_field_V_m": abs(
                scalar_at(sent_fields, "hQuasiFermiPotential", sj) -
                scalar_at(sent_fields, "hQuasiFermiPotential", si)) / float(edge["length_m"]),
            "vela_electron_flux_abs_m2_s": row["electron_flux_abs"],
            "sentaurus_electron_flux_abs_m2_s": sent_e_current * 1.0e4 / Q,
            "vela_hole_flux_abs_m2_s": row["hole_flux_abs"],
            "sentaurus_hole_flux_abs_m2_s": sent_h_current * 1.0e4 / Q,
            "vela_source_density_m3_s": (
                float(row["source_integral"]) / edge_area if edge_area > 0.0 else 0.0),
            "sentaurus_generation_avg_m3_s": sent_generation,
            "log10_vela_over_sentaurus_generation": log10_ratio(
                float(row["source_integral"]) / edge_area if edge_area > 0.0 else 0.0,
                sent_generation),
            "vela_electron_density_avg_cm3": n_vela,
            "sentaurus_electron_density_avg_cm3": n_sent,
            "log10_vela_over_sentaurus_electron_density": log10_ratio(n_vela, n_sent),
            "vela_hole_density_avg_cm3": p_vela,
            "sentaurus_hole_density_avg_cm3": p_sent,
            "log10_vela_over_sentaurus_hole_density": log10_ratio(p_vela, p_sent),
        })

    for node_id in selected_nodes:
        node = nodes[node_id]
        sent_id, distance_um = nearest_node(sent_nodes, node["x_um"], node["y_um"])
        volume = volumes[node_id]
        vela_n = scalars["Electrons"][node_id] * 1.0e-6
        vela_p = scalars["Holes"][node_id] * 1.0e-6
        sent_n = scalar_at(sent_fields, "eDensity", sent_id)
        sent_p = scalar_at(sent_fields, "hDensity", sent_id)
        vela_psi = scalars["Potential"][node_id]
        sent_psi = scalar_at(sent_fields, "ElectrostaticPotential", sent_id)
        vela_phin = scalars["ElectronQuasiFermi"][node_id]
        sent_phin = scalar_at(sent_fields, "eQuasiFermiPotential", sent_id)
        vela_phip = scalars["HoleQuasiFermi"][node_id]
        sent_phip = scalar_at(sent_fields, "hQuasiFermiPotential", sent_id)
        vela_psi_minus_phin = vela_psi - vela_phin
        sent_psi_minus_phin = sent_psi - sent_phin
        vela_phip_minus_psi = vela_phip - vela_psi
        sent_phip_minus_psi = sent_phip - sent_psi
        vela_srh = srh[node_id] * volume
        vela_avalanche = node_avalanche[node_id]
        sent_gen = scalar_at(sent_fields, "ImpactIonization", sent_id) * 1.0e6 * volume
        sent_srh = scalar_at(sent_fields, "srhRecombination", sent_id) * 1.0e6 * volume
        relation = "focus_endpoint" if node_id in focus_nodes else "incident_neighbor"
        node_rows.append({
            "bias_V": bias,
            "focus_edge_id": focus_edge_id,
            "node_id": node_id,
            "relation": relation,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "sentaurus_node_id": sent_id,
            "sentaurus_distance_um": distance_um,
            "node_volume_m2": volume,
            "vela_psi_V": vela_psi,
            "sentaurus_psi_V": sent_psi,
            "vela_phin_V": vela_phin,
            "sentaurus_phin_V": sent_phin,
            "vela_phip_V": vela_phip,
            "sentaurus_phip_V": sent_phip,
            "vela_psi_minus_phin_V": vela_psi_minus_phin,
            "sentaurus_psi_minus_phin_V": sent_psi_minus_phin,
            "delta_psi_minus_phin_V": vela_psi_minus_phin - sent_psi_minus_phin,
            "vela_phip_minus_psi_V": vela_phip_minus_psi,
            "sentaurus_phip_minus_psi_V": sent_phip_minus_psi,
            "delta_phip_minus_psi_V": vela_phip_minus_psi - sent_phip_minus_psi,
            "vela_electron_density_cm3": vela_n,
            "sentaurus_electron_density_cm3": sent_n,
            "log10_vela_over_sentaurus_electron_density": log10_ratio(vela_n, sent_n),
            "vela_hole_density_cm3": vela_p,
            "sentaurus_hole_density_cm3": sent_p,
            "log10_vela_over_sentaurus_hole_density": log10_ratio(vela_p, sent_p),
            "vela_ni_eff_cm3": ni[node_id] * 1.0e-6,
            "sentaurus_ni_eff_from_electron_cm3": inferred_ni_cm3(
                sent_n,
                sent_psi,
                sent_phin,
                vt,
                "electron"),
            "sentaurus_ni_eff_from_hole_cm3": inferred_ni_cm3(
                sent_p,
                sent_psi,
                sent_phip,
                vt,
                "hole"),
            "vela_electron_density_reconstructed_cm3": density_reconstruction(
                ni[node_id],
                vela_psi,
                vela_phin,
                vt,
                "electron") * 1.0e-6,
            "vela_hole_density_reconstructed_cm3": density_reconstruction(
                ni[node_id],
                vela_psi,
                vela_phip,
                vt,
                "hole") * 1.0e-6,
            "vela_electron_density_from_sentaurus_qf_cm3": density_reconstruction(
                ni[node_id],
                sent_psi,
                sent_phin,
                vt,
                "electron") * 1.0e-6,
            "vela_hole_density_from_sentaurus_qf_cm3": density_reconstruction(
                ni[node_id],
                sent_psi,
                sent_phip,
                vt,
                "hole") * 1.0e-6,
            "vela_electron_transport_integral_s_inv": electron_transport[node_id],
            "vela_hole_transport_integral_s_inv": hole_transport[node_id],
            "vela_avalanche_node_integral_s_inv": vela_avalanche,
            "vela_srh_node_integral_s_inv": vela_srh,
            "vela_electron_residual_estimate_s_inv": (
                electron_transport[node_id] + vela_srh - vela_avalanche),
            "vela_hole_residual_estimate_s_inv": (
                hole_transport[node_id] + vela_srh - vela_avalanche),
            "sentaurus_generation_node_integral_s_inv": sent_gen,
            "sentaurus_srh_node_integral_s_inv": sent_srh,
        })

    focus_row = source_by_edge[focus_edge_id]
    return {
        "bias_V": bias,
        "focus_edge_id": focus_edge_id,
        "focus_source_rank": focus_row["rank"],
        "selected_edge_count": len(selected),
        "selected_node_count": len(selected_nodes),
    }


def main() -> int:
    args = parse_args()
    vt = sgdiag.K_B_OVER_Q * args.temperature_k
    vtks = discover_vela_vtks(args.vela_vtk_root)
    edge_rows: list[dict[str, Any]] = []
    node_rows: list[dict[str, Any]] = []
    summaries = []
    for bias in parse_biases(args.biases):
        vtk = vtks.get(bias_key(bias))
        if vtk is None:
            raise FileNotFoundError(f"missing Vela VTK for {bias:g} V")
        summaries.append(append_bias_rows(
            bias,
            args.edge_id,
            args.mesh,
            args.doping_csv,
            vtk,
            sentaurus_dir(args.sentaurus_root, bias),
            vt,
            edge_rows,
            node_rows,
        ))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "continuity_feedback_edges.csv", edge_rows, EDGE_FIELDS)
    write_csv(args.out_dir / "continuity_feedback_nodes.csv", node_rows, NODE_FIELDS)
    (args.out_dir / "continuity_feedback_summary.json").write_text(
        json.dumps(clean_json({
            "edge_id": args.edge_id,
            "biases": summaries,
            "edge_rows": len(edge_rows),
            "node_rows": len(node_rows),
        }), indent=2) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
