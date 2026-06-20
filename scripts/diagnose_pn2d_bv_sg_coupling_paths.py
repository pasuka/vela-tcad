#!/usr/bin/env python3
"""Trace strong SG carrier-coupling paths from active BV support to contacts."""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_active_edge_flux_factors as fluxfac
import diagnose_pn2d_bv_continuity_feedback as cont
import diagnose_pn2d_bv_sg_flux_forms as fluxforms


EDGE_FIELDS = [
    "bias_V",
    "coupling_state",
    "carrier",
    "start_node",
    "target_contact",
    "path_found",
    "step",
    "edge_id",
    "node_from",
    "node_to",
    "x_from_um",
    "y_from_um",
    "x_to_um",
    "y_to_um",
    "edge_class",
    "contact_names",
    "length_m",
    "couple_m",
    "coupling_psi_from_V",
    "coupling_psi_to_V",
    "coupling_qf_from_V",
    "coupling_qf_to_V",
    "coupling_ni_from_m3",
    "coupling_ni_to_m3",
    "coupling_density_exp_from",
    "coupling_density_exp_to",
    "coupling_eta",
    "coupling_forward_flux_from_m2_s",
    "coupling_backward_flux_from_m2_s",
    "coupling_derivative_from_s_inv_per_V",
    "coupling_forward_flux_to_m2_s",
    "coupling_backward_flux_to_m2_s",
    "coupling_derivative_to_s_inv_per_V",
    "vela_psi_drop_V",
    "vela_qf_drop_V",
    "sentaurus_node_from",
    "sentaurus_node_to",
    "sentaurus_distance_from_um",
    "sentaurus_distance_to_um",
    "sentaurus_psi_drop_V",
    "sentaurus_qf_drop_V",
    "delta_qf_drop_V",
    "vela_flux_signed_m2_s",
    "vela_flux_integral_s_inv",
    "coupling_abs_s_inv_per_V",
    "path_cumulative_vela_qf_drop_V",
    "path_cumulative_sentaurus_qf_drop_V",
]


SUMMARY_FIELDS = [
    "bias_V",
    "coupling_state",
    "carrier",
    "start_node",
    "target_contact",
    "target_node",
    "path_found",
    "edge_count",
    "total_cost",
    "vela_total_qf_drop_V",
    "sentaurus_total_qf_drop_V",
    "delta_total_qf_drop_V",
    "vela_total_psi_drop_V",
    "sentaurus_total_psi_drop_V",
    "min_coupling_abs_s_inv_per_V",
    "median_coupling_abs_s_inv_per_V",
    "max_coupling_abs_s_inv_per_V",
    "path_edge_ids",
    "path_node_ids",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--support-csv", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--nodes", default="", help="Comma-separated node ids; default uses active support nodes.")
    parser.add_argument("--target-contact", default="", help="Optional contact name; default accepts any contact.")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--material-ni-m3", type=float, default=1.4638914958767616e16)
    parser.add_argument("--bandgap-narrowing", choices=["none", "old_slotboom"], default="old_slotboom")
    parser.add_argument("--min-coupling-s-inv-per-v", type=float, default=0.0)
    parser.add_argument("--delta-v", type=float, default=1.0e-4)
    parser.add_argument("--coupling-state", choices=["vela", "sentaurus"], default="vela")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_node_filter(raw: str) -> set[int] | None:
    if not raw.strip():
        return None
    return {int(item.strip()) for item in raw.split(",") if item.strip()}


def support_nodes(path: Path, node_filter: set[int] | None) -> list[int]:
    nodes: list[int] = []
    for row in read_csv_rows(path):
        if row.get("support_class") in (None, "", "inactive"):
            continue
        node_id = int(row["node_id"])
        if node_filter is not None and node_id not in node_filter:
            continue
        nodes.append(node_id)
    return sorted(set(nodes))


def carrier_qf_name(carrier: str) -> str:
    return "ElectronQuasiFermi" if carrier == "electron" else "HoleQuasiFermi"


def sentaurus_qf_name(carrier: str) -> str:
    return "eQuasiFermiPotential" if carrier == "electron" else "hQuasiFermiPotential"


def signed_edge_flux(row: dict[str, Any], carrier: str) -> float:
    return float(row[f"{carrier}_flux_signed_m2_s"])


def edge_flux_with_qf(
    row: dict[str, Any],
    carrier: str,
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    qf0: float,
    qf1: float,
) -> float:
    i = int(row["node0"])
    j = int(row["node1"])
    h = float(row["length_m"])
    if h <= 1.0e-30:
        return 0.0
    psi = scalars["Potential"]
    if carrier == "electron":
        mu = float(row["electron_mobility_m2_V_s"])
        return cont.sgdiag.sg_electron_flux_qf_variable_ni(
            ni[i], ni[j], psi[i], psi[j], qf0, qf1, vt, mu * vt / h)
    mu = float(row["hole_mobility_m2_V_s"])
    return cont.sgdiag.sg_hole_flux_qf_variable_ni(
        ni[i], ni[j], psi[i], psi[j], qf0, qf1, vt, mu * vt / h)


def edge_differential_coupling(
    row: dict[str, Any],
    carrier: str,
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    delta_v: float,
) -> float:
    i = int(row["node0"])
    j = int(row["node1"])
    qf = scalars[carrier_qf_name(carrier)]
    qf0 = qf[i]
    qf1 = qf[j]
    couple = float(row["couple_m"])
    derivatives: list[float] = []
    for perturb0, perturb1 in ((delta_v, 0.0), (0.0, delta_v)):
        forward = edge_flux_with_qf(row, carrier, scalars, ni, vt, qf0 + perturb0, qf1 + perturb1)
        backward = edge_flux_with_qf(row, carrier, scalars, ni, vt, qf0 - perturb0, qf1 - perturb1)
        derivatives.append(abs((forward - backward) * couple / (2.0 * delta_v)))
    return max(derivatives)


def edge_flux_with_node_perturbation(
    row: dict[str, Any],
    carrier: str,
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    node_id: int,
    delta_v: float,
) -> tuple[float, float, float]:
    i = int(row["node0"])
    j = int(row["node1"])
    qf = scalars[carrier_qf_name(carrier)]
    qf0 = qf[i] + (delta_v if node_id == i else 0.0)
    qf1 = qf[j] + (delta_v if node_id == j else 0.0)
    forward = edge_flux_with_qf(row, carrier, scalars, ni, vt, qf0, qf1)
    qf0 = qf[i] - (delta_v if node_id == i else 0.0)
    qf1 = qf[j] - (delta_v if node_id == j else 0.0)
    backward = edge_flux_with_qf(row, carrier, scalars, ni, vt, qf0, qf1)
    derivative = (forward - backward) * float(row["couple_m"]) / (2.0 * delta_v)
    return forward, backward, derivative


def density_exponent(carrier: str, psi: float, qf: float, vt: float) -> float:
    if carrier == "electron":
        return (psi - qf) / vt
    return (qf - psi) / vt


def coupling_eta(carrier: str, psi_from: float, psi_to: float, ni_from: float, ni_to: float, vt: float) -> float:
    if ni_from <= 0.0 or ni_to <= 0.0:
        return math.nan
    if carrier == "electron":
        return (psi_to - psi_from) / vt + math.log(ni_to / ni_from)
    return (psi_to - psi_from) / vt + math.log(ni_from / ni_to)


def edge_coupling_diagnostics(
    row: dict[str, Any],
    carrier: str,
    node_from: int,
    node_to: int,
    coupling_scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    delta_v: float,
) -> dict[str, Any]:
    qf = coupling_scalars[carrier_qf_name(carrier)]
    psi = coupling_scalars["Potential"]
    forward_from, backward_from, derivative_from = edge_flux_with_node_perturbation(
        row, carrier, coupling_scalars, ni, vt, node_from, delta_v)
    forward_to, backward_to, derivative_to = edge_flux_with_node_perturbation(
        row, carrier, coupling_scalars, ni, vt, node_to, delta_v)
    return {
        "coupling_psi_from_V": psi[node_from],
        "coupling_psi_to_V": psi[node_to],
        "coupling_qf_from_V": qf[node_from],
        "coupling_qf_to_V": qf[node_to],
        "coupling_ni_from_m3": ni[node_from],
        "coupling_ni_to_m3": ni[node_to],
        "coupling_density_exp_from": density_exponent(carrier, psi[node_from], qf[node_from], vt),
        "coupling_density_exp_to": density_exponent(carrier, psi[node_to], qf[node_to], vt),
        "coupling_eta": coupling_eta(carrier, psi[node_from], psi[node_to], ni[node_from], ni[node_to], vt),
        "coupling_forward_flux_from_m2_s": forward_from,
        "coupling_backward_flux_from_m2_s": backward_from,
        "coupling_derivative_from_s_inv_per_V": derivative_from,
        "coupling_forward_flux_to_m2_s": forward_to,
        "coupling_backward_flux_to_m2_s": backward_to,
        "coupling_derivative_to_s_inv_per_V": derivative_to,
    }


def build_adjacency(
    rows: list[dict[str, Any]],
    carrier: str,
    scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    delta_v: float,
    min_coupling: float,
) -> dict[int, list[dict[str, Any]]]:
    adjacency: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        coupling = edge_differential_coupling(row, carrier, scalars, ni, vt, delta_v)
        if coupling <= min_coupling:
            continue
        i = int(row["node0"])
        j = int(row["node1"])
        cost = 1.0 / max(coupling, 1.0e-300)
        adjacency[i].append({"from": i, "to": j, "row": row, "cost": cost, "coupling": coupling})
        adjacency[j].append({"from": j, "to": i, "row": row, "cost": cost, "coupling": coupling})
    return adjacency


def contact_matches(node_id: int, contact_by_node: dict[int, set[str]], target_contact: str) -> bool:
    contacts = contact_by_node.get(node_id, set())
    if not contacts:
        return False
    return not target_contact or target_contact in contacts


def find_path(
    start: int,
    adjacency: dict[int, list[dict[str, Any]]],
    contact_by_node: dict[int, set[str]],
    target_contact: str,
) -> tuple[int | None, float, list[dict[str, Any]]]:
    heap: list[tuple[float, int]] = [(0.0, start)]
    best: dict[int, float] = {start: 0.0}
    previous: dict[int, tuple[int, dict[str, Any]]] = {}
    target: int | None = None
    while heap:
        cost, node_id = heapq.heappop(heap)
        if cost > best[node_id]:
            continue
        if node_id != start and contact_matches(node_id, contact_by_node, target_contact):
            target = node_id
            break
        for edge in adjacency.get(node_id, []):
            next_id = int(edge["to"])
            next_cost = cost + float(edge["cost"])
            if next_cost < best.get(next_id, math.inf):
                best[next_id] = next_cost
                previous[next_id] = (node_id, edge)
                heapq.heappush(heap, (next_cost, next_id))
    if target is None:
        return None, math.inf, []
    path: list[dict[str, Any]] = []
    node = target
    while node != start:
        prev, edge = previous[node]
        path.append(edge)
        node = prev
    path.reverse()
    return target, best[target], path


def nearest_cache(sentaurus_nodes: dict[int, dict[str, float]], nodes: dict[int, dict[str, float]]) -> dict[int, tuple[int, float]]:
    return {
        node_id: cont.nearest_node(sentaurus_nodes, node["x_um"], node["y_um"])
        for node_id, node in nodes.items()
    }


def load_sentaurus_path_fields(root: Path) -> dict[str, dict[int, list[float]]]:
    return {
        name: cont.load_components(root, name)
        for name in ("ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential")
    }


def sentaurus_coupling_scalars(
    nodes: dict[int, dict[str, float]],
    sentaurus_fields: dict[str, dict[int, list[float]]],
    sentaurus_nearest: dict[int, tuple[int, float]],
) -> dict[str, list[float]]:
    count = max(nodes) + 1
    replay = {
        "Potential": [0.0 for _ in range(count)],
        "ElectronQuasiFermi": [0.0 for _ in range(count)],
        "HoleQuasiFermi": [0.0 for _ in range(count)],
    }
    for node_id in nodes:
        sent_id, _ = sentaurus_nearest[node_id]
        replay["Potential"][node_id] = cont.scalar_at(sentaurus_fields, "ElectrostaticPotential", sent_id)
        replay["ElectronQuasiFermi"][node_id] = cont.scalar_at(sentaurus_fields, "eQuasiFermiPotential", sent_id)
        replay["HoleQuasiFermi"][node_id] = cont.scalar_at(sentaurus_fields, "hQuasiFermiPotential", sent_id)
    return replay


def finite_or_none(value: float) -> float | None:
    return value if math.isfinite(value) else None


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def path_target_contact(target_node: int | None, contact_by_node: dict[int, set[str]], requested: str) -> str:
    if requested:
        return requested
    if target_node is None:
        return ""
    return ";".join(sorted(contact_by_node.get(target_node, set())))


def make_path_rows(
    *,
    bias: float,
    coupling_state: str,
    carrier: str,
    start: int,
    target_node: int | None,
    target_contact: str,
    path: list[dict[str, Any]],
    nodes: dict[int, dict[str, float]],
    contact_by_node: dict[int, set[str]],
    scalars: dict[str, list[float]],
    coupling_scalars: dict[str, list[float]],
    ni: list[float],
    vt: float,
    delta_v: float,
    sentaurus_fields: dict[str, dict[int, list[float]]],
    sentaurus_nearest: dict[int, tuple[int, float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    qf_name = carrier_qf_name(carrier)
    sent_qf_name = sentaurus_qf_name(carrier)
    cumulative_vela = 0.0
    cumulative_sentaurus = 0.0
    for step, edge in enumerate(path, start=1):
        item = edge["row"]
        node_from = int(edge["from"])
        node_to = int(edge["to"])
        sent_from, sent_from_distance = sentaurus_nearest[node_from]
        sent_to, sent_to_distance = sentaurus_nearest[node_to]
        vela_psi_drop = scalars["Potential"][node_to] - scalars["Potential"][node_from]
        vela_qf_drop = scalars[qf_name][node_to] - scalars[qf_name][node_from]
        sent_psi_drop = (
            cont.scalar_at(sentaurus_fields, "ElectrostaticPotential", sent_to)
            - cont.scalar_at(sentaurus_fields, "ElectrostaticPotential", sent_from)
        )
        sent_qf_drop = (
            cont.scalar_at(sentaurus_fields, sent_qf_name, sent_to)
            - cont.scalar_at(sentaurus_fields, sent_qf_name, sent_from)
        )
        cumulative_vela += vela_qf_drop
        cumulative_sentaurus += sent_qf_drop
        contacts = cont.sgdiag.contact_names_for_edge(item, contact_by_node)
        coupling_diag = edge_coupling_diagnostics(
            item, carrier, node_from, node_to, coupling_scalars, ni, vt, delta_v)
        rows.append({
            "bias_V": bias,
            "coupling_state": coupling_state,
            "carrier": carrier,
            "start_node": start,
            "target_contact": path_target_contact(target_node, contact_by_node, target_contact),
            "path_found": True,
            "step": step,
            "edge_id": int(item["edge_id"]),
            "node_from": node_from,
            "node_to": node_to,
            "x_from_um": nodes[node_from]["x_um"],
            "y_from_um": nodes[node_from]["y_um"],
            "x_to_um": nodes[node_to]["x_um"],
            "y_to_um": nodes[node_to]["y_um"],
            "edge_class": item["edge_class"],
            "contact_names": ";".join(sorted(contacts)),
            "length_m": item["length_m"],
            "couple_m": item["couple_m"],
            **coupling_diag,
            "vela_psi_drop_V": vela_psi_drop,
            "vela_qf_drop_V": vela_qf_drop,
            "sentaurus_node_from": sent_from,
            "sentaurus_node_to": sent_to,
            "sentaurus_distance_from_um": sent_from_distance,
            "sentaurus_distance_to_um": sent_to_distance,
            "sentaurus_psi_drop_V": sent_psi_drop,
            "sentaurus_qf_drop_V": sent_qf_drop,
            "delta_qf_drop_V": vela_qf_drop - sent_qf_drop,
            "vela_flux_signed_m2_s": signed_edge_flux(item, carrier),
            "vela_flux_integral_s_inv": signed_edge_flux(item, carrier) * float(item["couple_m"]),
            "coupling_abs_s_inv_per_V": edge["coupling"],
            "path_cumulative_vela_qf_drop_V": cumulative_vela,
            "path_cumulative_sentaurus_qf_drop_V": cumulative_sentaurus,
        })
    return rows


def make_summary_row(
    *,
    bias: float,
    coupling_state: str,
    carrier: str,
    start: int,
    target_node: int | None,
    target_contact: str,
    total_cost: float,
    path: list[dict[str, Any]],
    edge_rows: list[dict[str, Any]],
    contact_by_node: dict[int, set[str]],
) -> dict[str, Any]:
    couplings = [float(row["coupling_abs_s_inv_per_V"]) for row in edge_rows]
    vela_qf = sum(float(row["vela_qf_drop_V"]) for row in edge_rows)
    sent_qf = sum(float(row["sentaurus_qf_drop_V"]) for row in edge_rows)
    vela_psi = sum(float(row["vela_psi_drop_V"]) for row in edge_rows)
    sent_psi = sum(float(row["sentaurus_psi_drop_V"]) for row in edge_rows)
    path_nodes = [str(start)]
    for item in edge_rows:
        path_nodes.append(str(item["node_to"]))
    return {
        "bias_V": bias,
        "coupling_state": coupling_state,
        "carrier": carrier,
        "start_node": start,
        "target_contact": path_target_contact(target_node, contact_by_node, target_contact),
        "target_node": target_node if target_node is not None else "",
        "path_found": bool(path),
        "edge_count": len(path),
        "total_cost": finite_or_none(total_cost),
        "vela_total_qf_drop_V": vela_qf,
        "sentaurus_total_qf_drop_V": sent_qf,
        "delta_total_qf_drop_V": vela_qf - sent_qf,
        "vela_total_psi_drop_V": vela_psi,
        "sentaurus_total_psi_drop_V": sent_psi,
        "min_coupling_abs_s_inv_per_V": min(couplings) if couplings else None,
        "median_coupling_abs_s_inv_per_V": statistics.median(couplings) if couplings else None,
        "max_coupling_abs_s_inv_per_V": max(couplings) if couplings else None,
        "path_edge_ids": ";".join(str(edge["row"]["edge_id"]) for edge in path),
        "path_node_ids": ";".join(path_nodes),
    }


def main() -> int:
    args = parse_args()
    node_filter = parse_node_filter(args.nodes)
    selected_nodes = support_nodes(args.support_csv, node_filter)
    if not selected_nodes:
        raise SystemExit("no active support nodes selected")

    nodes, triangles, contact_by_node = cont.sgdiag.read_mesh(args.mesh)
    edges = cont.sgdiag.build_edges(nodes, triangles)
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
    source_rows, _ = cont.source_rows_with_signed_flux(nodes, edges, contact_by_node, scalars, ni, vt)
    sentaurus_nodes = cont.load_sentaurus_nodes(args.sentaurus_dir)
    sentaurus_fields = load_sentaurus_path_fields(args.sentaurus_dir)
    sentaurus_nearest = nearest_cache(sentaurus_nodes, nodes)
    coupling_scalars = (
        sentaurus_coupling_scalars(nodes, sentaurus_fields, sentaurus_nearest)
        if args.coupling_state == "sentaurus"
        else scalars
    )

    edge_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for carrier in ("electron", "hole"):
        adjacency = build_adjacency(
            source_rows,
            carrier,
            coupling_scalars,
            ni,
            vt,
            args.delta_v,
            args.min_coupling_s_inv_per_v,
        )
        for start in selected_nodes:
            target_node, total_cost, path = find_path(
                start,
                adjacency,
                contact_by_node,
                args.target_contact,
            )
            rows = make_path_rows(
                bias=args.bias,
                coupling_state=args.coupling_state,
                carrier=carrier,
                start=start,
                target_node=target_node,
                target_contact=args.target_contact,
                path=path,
                nodes=nodes,
                contact_by_node=contact_by_node,
                scalars=scalars,
                coupling_scalars=coupling_scalars,
                ni=ni,
                vt=vt,
                delta_v=args.delta_v,
                sentaurus_fields=sentaurus_fields,
                sentaurus_nearest=sentaurus_nearest,
            )
            edge_rows.extend(rows)
            summary_rows.append(make_summary_row(
                bias=args.bias,
                coupling_state=args.coupling_state,
                carrier=carrier,
                start=start,
                target_node=target_node,
                target_contact=args.target_contact,
                total_cost=total_cost,
                path=path,
                edge_rows=rows,
                contact_by_node=contact_by_node,
            ))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sg_coupling_path_edges.csv", edge_rows, EDGE_FIELDS)
    write_csv(args.out_dir / "sg_coupling_path_summary.csv", summary_rows, SUMMARY_FIELDS)
    summary = {
        "bias_V": args.bias,
        "support_csv": str(args.support_csv),
        "mesh": str(args.mesh),
        "doping_csv": str(args.doping_csv),
        "sentaurus_dir": str(args.sentaurus_dir),
        "vela_vtk": str(vtk),
        "selected_nodes": selected_nodes,
        "target_contact": args.target_contact,
        "coupling_state": args.coupling_state,
        "path_count": len(summary_rows),
        "edge_row_count": len(edge_rows),
        "paths_found": sum(1 for row in summary_rows if row["path_found"]),
        "delta_v": args.delta_v,
        "coupling_definition": "max_abs_central_difference_of_edge_flux_integral_wrt_endpoint_qf",
    }
    (args.out_dir / "sg_coupling_path_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json(summary), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
