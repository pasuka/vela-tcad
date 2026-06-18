#!/usr/bin/env python3
"""Decompose PN2D SG avalanche source by edge class and boundary proximity."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


K_B_OVER_Q = 8.617333262145e-5


EDGE_FIELDS = [
    "rank",
    "edge_id",
    "node0",
    "node1",
    "x0_um",
    "y0_um",
    "x1_um",
    "y1_um",
    "length_m",
    "couple_m",
    "edge_area_m2",
    "adjacent_cell_count",
    "edge_class",
    "contact_proximity",
    "contact_names",
    "junction_crossing",
    "electron_source_integral",
    "hole_source_integral",
    "source_integral",
    "source_fraction",
    "electron_flux_abs",
    "hole_flux_abs",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "electron_field_V_m",
    "hole_field_V_m",
    "electric_field_V_m",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vtk", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float)
    parser.add_argument("--top", type=int, default=80)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--doping-csv", type=Path)
    parser.add_argument("--material-ni-m3", type=float, default=1.0e16)
    parser.add_argument(
        "--bandgap-narrowing",
        choices=["none", "old_slotboom"],
        default="old_slotboom",
    )
    parser.add_argument(
        "--driving-force",
        choices=["quasi_fermi_gradient", "electric_field"],
        default="quasi_fermi_gradient",
    )
    parser.add_argument(
        "--flux-form",
        choices=["quasi_fermi_variable_ni", "density"],
        default="quasi_fermi_variable_ni",
    )
    return parser.parse_args()


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    scalars: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS",
                    "VECTORS",
                    "FIELD",
                    "CELL_DATA",
                    "POINT_DATA",
                }:
                    break
                values.extend(float(value) for value in next_parts)
                i += 1
            scalars[name] = values
            continue
        i += 1
    return scalars


def read_mesh(path: Path) -> tuple[dict[int, dict[str, float]],
                                   list[dict[str, Any]],
                                   dict[int, set[str]]]:
    data = json.loads(path.read_text())
    nodes = {
        int(node["id"]): {"x_um": float(node["x"]), "y_um": float(node["y"])}
        for node in data["nodes"]
    }
    triangles = [
        {
            "id": int(triangle["id"]),
            "node_ids": [int(node_id) for node_id in triangle["node_ids"]],
        }
        for triangle in data["triangles"]
    ]
    contact_by_node: dict[int, set[str]] = defaultdict(set)
    for contact in data.get("contacts", []):
        name = str(contact.get("name", "contact"))
        for node_id in contact.get("node_ids", []):
            contact_by_node[int(node_id)].add(name)
    return nodes, triangles, contact_by_node


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a <= b else (b, a)


def point_m(nodes: dict[int, dict[str, float]], node_id: int) -> tuple[float, float]:
    node = nodes[node_id]
    return node["x_um"] * 1.0e-6, node["y_um"] * 1.0e-6


def distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(b[0] - a[0], b[1] - a[1])


def triangle_area(points: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = points
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))


def cotangent_at_opposite(a: tuple[float, float],
                          b: tuple[float, float],
                          opp: tuple[float, float]) -> float:
    ux = a[0] - opp[0]
    uy = a[1] - opp[1]
    vx = b[0] - opp[0]
    vy = b[1] - opp[1]
    cross = ux * vy - uy * vx
    if abs(cross) < 1.0e-30:
        return 0.0
    dot = ux * vx + uy * vy
    return dot / abs(cross)


def build_edges(nodes: dict[int, dict[str, float]],
                triangles: list[dict[str, Any]]) -> list[dict[str, Any]]:
    edge_map: dict[tuple[int, int], int] = {}
    edges: list[dict[str, Any]] = []
    for triangle in triangles:
        ids = triangle["node_ids"]
        for k in range(3):
            a, b = edge_key(ids[k], ids[(k + 1) % 3])
            key = (a, b)
            if key not in edge_map:
                pa = point_m(nodes, a)
                pb = point_m(nodes, b)
                edge_map[key] = len(edges)
                edges.append({
                    "edge_id": len(edges),
                    "node0": a,
                    "node1": b,
                    "length_m": distance(pa, pb),
                    "couple_m": 0.0,
                    "cell_ids": [],
                })

    for triangle in triangles:
        ids = triangle["node_ids"]
        points = [point_m(nodes, node_id) for node_id in ids]
        area = triangle_area(points)
        if area <= 1.0e-30:
            continue
        for k in range(3):
            a = ids[k]
            b = ids[(k + 1) % 3]
            opp = ids[(k + 2) % 3]
            key = edge_key(a, b)
            edge = edges[edge_map[key]]
            edge["cell_ids"].append(triangle["id"])
            h = float(edge["length_m"])
            if h <= 1.0e-30:
                continue
            cot = cotangent_at_opposite(point_m(nodes, a), point_m(nodes, b), point_m(nodes, opp))
            local_couple = 0.5 * cot * h
            if cot < 0.0:
                local_couple = area / (3.0 * h)
            edge["couple_m"] += max(local_couple, 0.0)
    return edges


def bernoulli(u: float) -> float:
    if abs(u) < 1.0e-8:
        return 1.0 - u / 2.0 + u * u / 12.0
    if u > 700.0:
        return 0.0
    if u < -700.0:
        return -u
    return u / math.expm1(u)


def limited_exp(value: float) -> float:
    return math.exp(max(-500.0, min(500.0, value)))


def sg_electron_flux(n0: float, n1: float, dpsi: float, vt: float, coef: float) -> float:
    u = dpsi / vt
    return coef * (bernoulli(-u) * n0 - bernoulli(u) * n1)


def sg_hole_flux(p0: float, p1: float, dpsi: float, vt: float, coef: float) -> float:
    u = dpsi / vt
    return coef * (bernoulli(u) * p0 - bernoulli(-u) * p1)


def sg_electron_flux_qf_variable_ni(ni0: float,
                                    ni1: float,
                                    psi0: float,
                                    psi1: float,
                                    phin0: float,
                                    phin1: float,
                                    vt: float,
                                    coef: float) -> float:
    if ni0 <= 0.0 or ni1 <= 0.0:
        return sg_electron_flux(
            ni0 * limited_exp((psi0 - phin0) / vt),
            ni1 * limited_exp((psi1 - phin1) / vt),
            psi1 - psi0,
            vt,
            coef,
        )
    eta = (psi1 - psi0) / vt + math.log(ni1 / ni0)
    prefactor = ni1 * limited_exp(psi1 / vt)
    return coef * bernoulli(eta) * prefactor * (
        limited_exp(-phin0 / vt) - limited_exp(-phin1 / vt))


def sg_hole_flux_qf_variable_ni(ni0: float,
                                ni1: float,
                                psi0: float,
                                psi1: float,
                                phip0: float,
                                phip1: float,
                                vt: float,
                                coef: float) -> float:
    if ni0 <= 0.0 or ni1 <= 0.0:
        return sg_hole_flux(
            ni0 * limited_exp((phip0 - psi0) / vt),
            ni1 * limited_exp((phip1 - psi1) / vt),
            psi1 - psi0,
            vt,
            coef,
        )
    eta = (psi1 - psi0) / vt + math.log(ni0 / ni1)
    prefactor = ni0 * limited_exp(-psi0 / vt)
    return coef * bernoulli(eta) * prefactor * (
        limited_exp(phip0 / vt) - limited_exp(phip1 / vt))


def estimate_effective_ni(electrons: list[float],
                          holes: list[float],
                          psi: list[float],
                          phin: list[float],
                          phip: list[float],
                          vt: float) -> list[float]:
    ni: list[float] = []
    for n, p, potential, qfn, qfp in zip(electrons, holes, psi, phin, phip):
        estimates = []
        if n > 0.0:
            estimates.append(n * limited_exp((qfn - potential) / vt))
        if p > 0.0:
            estimates.append(p * limited_exp((potential - qfp) / vt))
        positive = [value for value in estimates if value > 0.0 and math.isfinite(value)]
        if len(positive) == 2:
            ni.append(math.sqrt(positive[0] * positive[1]))
        elif positive:
            ni.append(positive[0])
        else:
            ni.append(0.0)
    return ni


def read_total_impurity_m3(path: Path, node_count: int) -> list[float]:
    impurity = [0.0 for _ in range(node_count)]
    seen = [False for _ in range(node_count)]
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            donors = float(row.get("donors_cm3", 0.0) or 0.0) * 1.0e6
            acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0) * 1.0e6
            impurity[node_id] = abs(donors) + abs(acceptors)
            seen[node_id] = True
    missing = [index for index, present in enumerate(seen) if not present]
    if missing:
        raise RuntimeError(f"doping CSV is missing node ids: {missing[:8]}")
    return impurity


def old_slotboom_delta_eg_eV(impurity: float) -> float:
    if impurity <= 0.0:
        return 0.0
    reference = 1.0e23
    coefficient = 9.0e-3
    smoothing = 0.5
    offset = -1.595e-2
    x = math.log(impurity / reference)
    delta = offset + coefficient * (x + math.sqrt(x * x + smoothing))
    return max(delta, 0.0)


def effective_ni_from_doping(doping_csv: Path,
                             node_count: int,
                             material_ni_m3: float,
                             vt: float,
                             bandgap_narrowing: str) -> list[float]:
    impurity = read_total_impurity_m3(doping_csv, node_count)
    ni = []
    for total in impurity:
        delta = old_slotboom_delta_eg_eV(total) if bandgap_narrowing == "old_slotboom" else 0.0
        ni.append(material_ni_m3 * math.exp(delta / (2.0 * vt)) if delta > 0.0 else material_ni_m3)
    return ni


def van_overstraeten_alpha(field: float,
                           low_a: float,
                           high_a: float,
                           low_b: float,
                           high_b: float,
                           switch_field: float = 4.0e7) -> float:
    field = abs(field)
    if field <= 0.0:
        return 0.0
    prefactor = low_a if field < switch_field else high_a
    critical = low_b if field < switch_field else high_b
    return prefactor * math.exp(max(-700.0, min(0.0, -critical / field)))


def required_scalars(scalars: dict[str, list[float]]) -> None:
    required = [
        "Potential",
        "ElectronQuasiFermi",
        "HoleQuasiFermi",
        "Electrons",
        "Holes",
        "ElectronMobility",
        "HoleMobility",
        "AvalancheGeneration",
    ]
    missing = [name for name in required if name not in scalars]
    if missing:
        raise RuntimeError(f"VTK is missing required scalars: {', '.join(missing)}")


def field_between(values: list[float], i: int, j: int, h: float) -> float:
    if h <= 1.0e-30:
        return 0.0
    return abs(values[j] - values[i]) / h


def avg(values: list[float], i: int, j: int) -> float:
    return 0.5 * (values[i] + values[j])


def contact_names_for_edge(edge: dict[str, Any],
                           contact_by_node: dict[int, set[str]]) -> set[str]:
    return set(contact_by_node.get(edge["node0"], set())) | set(
        contact_by_node.get(edge["node1"], set()))


def classify_edge(edge: dict[str, Any],
                  contact_names: set[str],
                  net_doping: list[float] | None) -> tuple[str, bool]:
    is_boundary = len(edge["cell_ids"]) == 1
    crosses_junction = False
    if net_doping is not None:
        d0 = net_doping[edge["node0"]]
        d1 = net_doping[edge["node1"]]
        crosses_junction = (d0 < 0.0 < d1) or (d1 < 0.0 < d0)
    if is_boundary and contact_names:
        edge_class = "contact_boundary"
    elif is_boundary:
        edge_class = "boundary_noncontact"
    elif contact_names:
        edge_class = "contact_adjacent_interior"
    elif crosses_junction:
        edge_class = "interior_junction"
    else:
        edge_class = "interior_bulk"
    return edge_class, crosses_junction


def source_rows(nodes: dict[int, dict[str, float]],
                edges: list[dict[str, Any]],
                contact_by_node: dict[int, set[str]],
                scalars: dict[str, list[float]],
                vt: float,
                driving_force: str,
                flux_form: str,
                ni: list[float] | None = None) -> tuple[list[dict[str, Any]], list[float]]:
    required_scalars(scalars)
    psi = scalars["Potential"]
    phin = scalars["ElectronQuasiFermi"]
    phip = scalars["HoleQuasiFermi"]
    electrons = scalars["Electrons"]
    holes = scalars["Holes"]
    electron_mobility = scalars["ElectronMobility"]
    hole_mobility = scalars["HoleMobility"]
    net_doping = scalars.get("NetDoping")
    if ni is None:
        ni = estimate_effective_ni(electrons, holes, psi, phin, phip, vt)
    node_source = [0.0 for _ in nodes]
    rows: list[dict[str, Any]] = []
    for edge in edges:
        i = int(edge["node0"])
        j = int(edge["node1"])
        h = float(edge["length_m"])
        couple = float(edge["couple_m"])
        if h <= 1.0e-30 or couple <= 0.0:
            continue
        edge_area = 0.5 * h * couple
        electric_field = field_between(psi, i, j, h)
        electron_qf_field = field_between(phin, i, j, h)
        hole_qf_field = field_between(phip, i, j, h)
        electron_field = electron_qf_field if driving_force == "quasi_fermi_gradient" else electric_field
        hole_field = hole_qf_field if driving_force == "quasi_fermi_gradient" else electric_field
        alpha_n = van_overstraeten_alpha(electron_field, 7.03e7, 7.03e7, 1.231e8, 1.231e8)
        alpha_p = van_overstraeten_alpha(hole_field, 1.582e8, 6.71e7, 2.036e8, 1.693e8)
        mu_n = max(0.0, avg(electron_mobility, i, j))
        mu_p = max(0.0, avg(hole_mobility, i, j))
        dpsi = psi[j] - psi[i]
        if flux_form == "quasi_fermi_variable_ni":
            electron_flux = sg_electron_flux_qf_variable_ni(
                ni[i], ni[j], psi[i], psi[j], phin[i], phin[j], vt, mu_n * vt / h)
            hole_flux = sg_hole_flux_qf_variable_ni(
                ni[i], ni[j], psi[i], psi[j], phip[i], phip[j], vt, mu_p * vt / h)
        else:
            electron_flux = sg_electron_flux(
                electrons[i], electrons[j], dpsi, vt, mu_n * vt / h)
            hole_flux = sg_hole_flux(
                holes[i], holes[j], dpsi, vt, mu_p * vt / h)
        electron_source = alpha_n * abs(electron_flux) * edge_area
        hole_source = alpha_p * abs(hole_flux) * edge_area
        source = electron_source + hole_source
        half_source = 0.5 * source
        node_source[i] += half_source
        node_source[j] += half_source
        contacts = contact_names_for_edge(edge, contact_by_node)
        edge_class, crosses_junction = classify_edge(edge, contacts, net_doping)
        n0 = nodes[i]
        n1 = nodes[j]
        rows.append({
            "edge_id": int(edge["edge_id"]),
            "node0": i,
            "node1": j,
            "x0_um": n0["x_um"],
            "y0_um": n0["y_um"],
            "x1_um": n1["x_um"],
            "y1_um": n1["y_um"],
            "length_m": h,
            "couple_m": couple,
            "edge_area_m2": edge_area,
            "adjacent_cell_count": len(edge["cell_ids"]),
            "edge_class": edge_class,
            "contact_proximity": "contact_endpoint" if contacts else (
                "boundary_edge" if len(edge["cell_ids"]) == 1 else "interior"),
            "contact_names": ";".join(sorted(contacts)),
            "junction_crossing": crosses_junction,
            "electron_source_integral": electron_source,
            "hole_source_integral": hole_source,
            "source_integral": source,
            "electron_flux_abs": abs(electron_flux),
            "hole_flux_abs": abs(hole_flux),
            "electron_alpha_m_inv": alpha_n,
            "hole_alpha_m_inv": alpha_p,
            "electron_field_V_m": electron_field,
            "hole_field_V_m": hole_field,
            "electric_field_V_m": electric_field,
            "electron_mobility_m2_V_s": mu_n,
            "hole_mobility_m2_V_s": mu_p,
        })
    total = sum(row["source_integral"] for row in rows)
    for row in rows:
        row["source_fraction"] = row["source_integral"] / total if total > 0.0 else 0.0
    rows.sort(key=lambda item: abs(float(item["source_integral"])), reverse=True)
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank
    return rows, node_source


def node_volumes(nodes: dict[int, dict[str, float]],
                 triangles: list[dict[str, Any]]) -> list[float]:
    volumes = [0.0 for _ in nodes]
    for triangle in triangles:
        ids = triangle["node_ids"]
        area = triangle_area([point_m(nodes, node_id) for node_id in ids])
        for node_id in ids:
            volumes[node_id] += area / 3.0
    return volumes


def aggregate(rows: list[dict[str, Any]], key: str) -> list[dict[str, Any]]:
    totals: dict[str, dict[str, Any]] = {}
    total_source = sum(float(row["source_integral"]) for row in rows)
    for row in rows:
        group = str(row[key])
        item = totals.setdefault(group, {
            key: group,
            "edge_count": 0,
            "source_integral": 0.0,
            "electron_source_integral": 0.0,
            "hole_source_integral": 0.0,
            "max_edge_source_integral": 0.0,
        })
        source = float(row["source_integral"])
        item["edge_count"] += 1
        item["source_integral"] += source
        item["electron_source_integral"] += float(row["electron_source_integral"])
        item["hole_source_integral"] += float(row["hole_source_integral"])
        item["max_edge_source_integral"] = max(item["max_edge_source_integral"], source)
    result = []
    for item in totals.values():
        item["source_fraction"] = (
            item["source_integral"] / total_source if total_source > 0.0 else 0.0)
        result.append(item)
    result.sort(key=lambda item: item["source_integral"], reverse=True)
    return result


def node_classes(nodes: dict[int, dict[str, float]],
                 edges: list[dict[str, Any]],
                 contact_by_node: dict[int, set[str]]) -> list[str]:
    boundary_nodes: set[int] = set()
    for edge in edges:
        if len(edge["cell_ids"]) == 1:
            boundary_nodes.add(int(edge["node0"]))
            boundary_nodes.add(int(edge["node1"]))
    classes = []
    for node_id in range(len(nodes)):
        if contact_by_node.get(node_id):
            classes.append("contact")
        elif node_id in boundary_nodes:
            classes.append("boundary")
        else:
            classes.append("interior")
    return classes


def aggregate_nodes(classes: list[str],
                    vtk_node_integral: list[float],
                    reconstructed_node_source: list[float]) -> list[dict[str, Any]]:
    totals: dict[str, dict[str, Any]] = {}
    vtk_total = sum(vtk_node_integral)
    reconstructed_total = sum(reconstructed_node_source)
    for node_id, node_class in enumerate(classes):
        item = totals.setdefault(node_class, {
            "node_class": node_class,
            "node_count": 0,
            "vtk_source_integral": 0.0,
            "reconstructed_source_integral": 0.0,
        })
        item["node_count"] += 1
        item["vtk_source_integral"] += vtk_node_integral[node_id]
        item["reconstructed_source_integral"] += reconstructed_node_source[node_id]
    result = []
    for item in totals.values():
        item["vtk_source_fraction"] = (
            item["vtk_source_integral"] / vtk_total if vtk_total > 0.0 else 0.0)
        item["reconstructed_source_fraction"] = (
            item["reconstructed_source_integral"] / reconstructed_total
            if reconstructed_total > 0.0 else 0.0)
        result.append(item)
    result.sort(key=lambda item: item["vtk_source_integral"], reverse=True)
    return result


def log10_p95(reference: list[float], candidate: list[float]) -> float | None:
    errors = []
    eps = 1.0e-300
    for ref, cand in zip(reference, candidate):
        if ref > 0.0 or cand > 0.0:
            errors.append(abs(math.log10(max(abs(cand), eps)) - math.log10(max(abs(ref), eps))))
    if not errors:
        return None
    errors.sort()
    index = min(len(errors) - 1, int(math.ceil(0.95 * len(errors))) - 1)
    return errors[index]


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


def make_summary(args: argparse.Namespace,
                 rows: list[dict[str, Any]],
                 node_source: list[float],
                 vtk_node_integral: list[float],
                 node_class_values: list[str]) -> dict[str, Any]:
    total_source = sum(float(row["source_integral"]) for row in rows)
    total_vtk = sum(vtk_node_integral)
    top_source = [float(row["source_integral"]) for row in rows]
    def concentration(n: int) -> float:
        return sum(top_source[:min(n, len(top_source))]) / total_source if total_source > 0.0 else 0.0
    top = rows[0] if rows else {}
    return {
        "bias_V": args.bias,
        "vtk": str(args.vtk),
        "mesh": str(args.mesh),
        "driving_force": args.driving_force,
        "flux_form": args.flux_form,
        "ni_source": (
            f"doping_csv:{args.doping_csv};bandgap:{args.bandgap_narrowing}"
            if args.doping_csv else "vtk_state_backsolve"),
        "mobility_source": "VTK endpoint average",
        "edge_area_policy": "0.5 * edge.length * edge.couple",
        "edge_count": len(rows),
        "total_source_integral_reconstructed": total_source,
        "total_source_integral_vtk": total_vtk,
        "total_source_relative_error": (
            (total_source - total_vtk) / total_vtk if abs(total_vtk) > 0.0 else None),
        "node_source_log10_p95_error": log10_p95(vtk_node_integral, node_source),
        "top_edge": top,
        "concentration": {
            "top1_fraction": concentration(1),
            "top5_fraction": concentration(5),
            "top20_fraction": concentration(20),
        },
        "by_edge_class": aggregate(rows, "edge_class"),
        "by_contact_proximity": aggregate(rows, "contact_proximity"),
        "by_junction_crossing": aggregate(rows, "junction_crossing"),
        "by_node_class": aggregate_nodes(node_class_values, vtk_node_integral, node_source),
    }


def main() -> None:
    args = parse_args()
    if args.top <= 0:
        raise SystemExit("--top must be positive")
    nodes, triangles, contact_by_node = read_mesh(args.mesh)
    edges = build_edges(nodes, triangles)
    scalars = parse_vtk_scalars(args.vtk)
    vt = K_B_OVER_Q * args.temperature_k
    ni = (
        effective_ni_from_doping(
            args.doping_csv, len(nodes), args.material_ni_m3, vt, args.bandgap_narrowing)
        if args.doping_csv else None
    )
    rows, reconstructed_node_source = source_rows(
        nodes, edges, contact_by_node, scalars, vt, args.driving_force, args.flux_form, ni)
    volumes = node_volumes(nodes, triangles)
    vtk_node_integral = [
        scalars["AvalancheGeneration"][node_id] * volumes[node_id]
        for node_id in range(len(volumes))
    ]
    node_class_values = node_classes(nodes, edges, contact_by_node)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sg_avalanche_edges.csv", rows[:args.top], EDGE_FIELDS)
    summary = make_summary(
        args, rows, reconstructed_node_source, vtk_node_integral, node_class_values)
    (args.out_dir / "sg_avalanche_edge_summary.json").write_text(
        json.dumps(clean_json(summary), indent=2) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
