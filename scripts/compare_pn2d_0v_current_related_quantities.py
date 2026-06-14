#!/usr/bin/env python3
"""Extract and compare PN2D 0V current-related Sentaurus and Vela quantities."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


REPORT_NAME = "pn2d_0v_current_related_comparison.json"
MARKDOWN_NAME = "pn2d_0v_current_related_comparison.md"
COMPONENT_ZERO_FLOOR_A_PER_UM = 1.0e-24


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in fieldnames})


def as_float(value: Any, default: float = 0.0) -> float:
    if value in (None, ""):
        return default
    return float(value)


def field_files(root: Path, field: str) -> list[Path]:
    exact = root / f"{field}_region0.csv"
    if exact.is_file():
        return [exact]
    return sorted(root.glob(f"{field}_region*.csv"))


def read_sentaurus_scalar(fields_root: Path, field: str) -> dict[int, float]:
    values: dict[int, float] = {}
    paths = field_files(fields_root, field)
    if not paths:
        raise FileNotFoundError(f"missing Sentaurus field export: {field}")
    for path in paths:
        for row in read_csv_rows(path):
            values[int(float(row["node_id"]))] = float(row["component0"])
    return values


def read_nodes(path: Path) -> dict[int, dict[str, float]]:
    nodes: dict[int, dict[str, float]] = {}
    for row in read_csv_rows(path):
        node_id = int(float(row["id"]))
        nodes[node_id] = {
            "id": node_id,
            "x_um": float(row.get("x_um", row.get("x", 0.0))),
            "y_um": float(row.get("y_um", row.get("y", 0.0))),
        }
    return nodes


def read_contacts(path: Path) -> dict[str, list[int]]:
    contacts: dict[str, list[int]] = {}
    for row in read_csv_rows(path):
        contacts[row["name"]] = [
            int(item) for item in row.get("node_ids", "").split(";") if item.strip()
        ]
    return contacts


def read_elements(path: Path) -> list[list[int]]:
    elements: list[list[int]] = []
    for row in read_csv_rows(path):
        elements.append([int(row["node0"]), int(row["node1"]), int(row["node2"])])
    return elements


def parse_vtk(path: Path) -> tuple[list[tuple[float, float, float]], dict[str, list[float]]]:
    lines = path.read_text().splitlines()
    points: list[tuple[float, float, float]] = []
    fields: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "POINTS":
            npoints = int(parts[1])
            i += 1
            raw: list[float] = []
            while i < len(lines) and len(raw) < 3 * npoints:
                raw.extend(float(item) for item in lines[i].split())
                i += 1
            points = [
                (raw[3 * idx], raw[3 * idx + 1], raw[3 * idx + 2])
                for idx in range(npoints)
            ]
            continue
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < len(points):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                head = line.split()[0]
                if head in {"SCALARS", "CELL_DATA", "POINT_DATA", "VECTORS"}:
                    break
                values.extend(float(item) for item in line.split())
                i += 1
            fields[name] = values[:len(points)]
            continue
        i += 1
    return points, fields


def triangle_gradient(points: list[tuple[float, float, float]],
                      values: list[float],
                      ids: list[int]) -> float | None:
    i0, i1, i2 = ids
    if max(ids) >= len(points) or max(ids) >= len(values):
        return None
    x0, y0, _ = points[i0]
    x1, y1, _ = points[i1]
    x2, y2, _ = points[i2]
    v0, v1, v2 = values[i0], values[i1], values[i2]
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if abs(det) <= 1.0e-300:
        return None
    dv_dx = ((v1 - v0) * (y2 - y0) - (v2 - v0) * (y1 - y0)) / det
    dv_dy = ((x1 - x0) * (v2 - v0) - (x2 - x0) * (v1 - v0)) / det
    return math.hypot(dv_dx, dv_dy)


def derive_nodal_gradient(points: list[tuple[float, float, float]],
                          values: list[float],
                          elements: list[list[int]]) -> dict[int, float]:
    accum: dict[int, list[float]] = {}
    for ids in elements:
        magnitude = triangle_gradient(points, values, ids)
        if magnitude is None:
            continue
        for node_id in ids:
            accum.setdefault(node_id, []).append(magnitude)
    return {node_id: sum(vals) / len(vals) for node_id, vals in accum.items() if vals}


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    index = int(round((len(ordered) - 1) * p))
    return ordered[index]


def stats(values: list[float]) -> dict[str, Any]:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return {"points": 0, "min": None, "mean": None, "median": None, "p95": None, "max": None}
    return {
        "points": len(finite),
        "min": min(finite),
        "mean": sum(finite) / len(finite),
        "median": statistics.median(finite),
        "p95": quantile(finite, 0.95),
        "max": max(finite),
    }


def compare_maps(reference: dict[int, float],
                 candidate: dict[int, float],
                 groups: dict[str, set[int]]) -> dict[str, Any]:
    common = sorted(set(reference) & set(candidate))
    abs_diffs = [abs(candidate[node] - reference[node]) for node in common]
    rel_diffs = [
        abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
        for node in common
    ]
    group_stats = {}
    for name, ids in groups.items():
        selected = [node for node in common if node in ids]
        selected_abs = [abs(candidate[node] - reference[node]) for node in selected]
        selected_rel = [
            abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
            for node in selected
        ]
        group_stats[name] = {
            "points": len(selected),
            "mean_abs_diff": sum(selected_abs) / len(selected_abs) if selected_abs else None,
            "max_abs_diff": max(selected_abs) if selected_abs else None,
            "mean_rel_diff": sum(selected_rel) / len(selected_rel) if selected_rel else None,
            "max_rel_diff": max(selected_rel) if selected_rel else None,
        }
    worst = None
    if common:
        worst_node = max(common, key=lambda node: abs(candidate[node] - reference[node]))
        worst = {
            "node_id": worst_node,
            "sentaurus": reference[worst_node],
            "vela": candidate[worst_node],
            "abs_diff": abs(candidate[worst_node] - reference[worst_node]),
            "rel_diff": abs(candidate[worst_node] - reference[worst_node])
            / max(abs(reference[worst_node]), 1.0e-300),
        }
    return {
        "points_compared": len(common),
        "sentaurus_points": len(reference),
        "vela_points": len(candidate),
        "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
        "median_abs_diff": statistics.median(abs_diffs) if abs_diffs else None,
        "p95_abs_diff": quantile(abs_diffs, 0.95),
        "max_abs_diff": max(abs_diffs) if abs_diffs else None,
        "mean_rel_diff": sum(rel_diffs) / len(rel_diffs) if rel_diffs else None,
        "median_rel_diff": statistics.median(rel_diffs) if rel_diffs else None,
        "p95_rel_diff": quantile(rel_diffs, 0.95),
        "max_rel_diff": max(rel_diffs) if rel_diffs else None,
        "worst": worst,
        "groups": group_stats,
    }


def contact_flux_by_name(fields_root: Path, manifest_path: Path) -> dict[str, float]:
    manifest = read_json(manifest_path)
    out: dict[str, float] = {}
    for entry in manifest.get("fields", []):
        if entry.get("name") != "ContactCurrentFlux":
            continue
        region = int(entry["region"])
        region_name = str(entry["region_name"])
        path = fields_root / f"ContactCurrentFlux_region{region}.csv"
        rows = read_csv_rows(path)
        if rows:
            out[region_name] = float(rows[0]["component0"])
    return out


def contact_boundary_edges(elements: list[list[int]], contact_nodes: set[int]) -> list[tuple[int, int]]:
    counts: dict[tuple[int, int], int] = {}
    for a, b, c in elements:
        for u, v in ((a, b), (b, c), (c, a)):
            edge = tuple(sorted((u, v)))
            counts[edge] = counts.get(edge, 0) + 1
    return [
        edge for edge, count in counts.items()
        if count == 1 and edge[0] in contact_nodes and edge[1] in contact_nodes
    ]


def sentaurus_current_density_boundary_integrals(
    fields_root: Path,
    nodes: dict[int, dict[str, float]],
    elements: list[list[int]],
    contacts: dict[str, list[int]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    fields = {
        "electron": read_sentaurus_scalar(fields_root, "eCurrentDensity"),
        "hole": read_sentaurus_scalar(fields_root, "hCurrentDensity"),
        "total": read_sentaurus_scalar(fields_root, "TotalCurrentDensity"),
    }
    rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    for contact, ids in contacts.items():
        edges = contact_boundary_edges(elements, set(ids))
        by_quantity: dict[str, dict[str, float]] = {}
        for quantity, values in fields.items():
            signed_sum = 0.0
            abs_sum = 0.0
            edge_values: list[float] = []
            for node0, node1 in edges:
                n0 = nodes[node0]
                n1 = nodes[node1]
                length_um = math.hypot(n1["x_um"] - n0["x_um"], n1["y_um"] - n0["y_um"])
                length_cm = length_um * 1.0e-4
                j0 = values.get(node0)
                j1 = values.get(node1)
                if j0 is None or j1 is None:
                    continue
                avg_j = 0.5 * (j0 + j1)
                edge_current = avg_j * length_cm
                signed_sum += edge_current
                abs_sum += abs(edge_current)
                edge_values.append(edge_current)
                rows.append({
                    "contact": contact,
                    "quantity": quantity,
                    "node0": node0,
                    "node1": node1,
                    "edge_length_um": length_um,
                    "edge_length_cm": length_cm,
                    "j0_A_per_cm2": j0,
                    "j1_A_per_cm2": j1,
                    "avg_j_A_per_cm2": avg_j,
                    "edge_current_A_per_cm_width": edge_current,
                    "direction_status": "direction_unresolved",
                })
            by_quantity[quantity] = {
                "edges": len(edge_values),
                "signed_sum_A_per_cm_width": signed_sum,
                "abs_sum_A_per_cm_width": abs_sum,
                "mean_edge_current_A_per_cm_width": (
                    sum(edge_values) / len(edge_values) if edge_values else 0.0
                ),
                "max_abs_edge_current_A_per_cm_width": (
                    max((abs(value) for value in edge_values), default=0.0)
                ),
                "direction_status": "direction_unresolved",
            }
        summary[contact] = {
            "boundary_edges": len(edges),
            "by_quantity": by_quantity,
        }
    return rows, summary


def make_groups(nodes: dict[int, dict[str, float]], contacts: dict[str, list[int]]) -> dict[str, set[int]]:
    all_ids = set(nodes)
    groups: dict[str, set[int]] = {"all": all_ids}
    for name, ids in contacts.items():
        groups[f"contact:{name}"] = set(ids)
    if nodes:
        xs = [node["x_um"] for node in nodes.values()]
        x_mid = 0.5 * (min(xs) + max(xs))
        halfwidth = max((max(xs) - min(xs)) * 0.02, 0.02)
        groups["junction_near"] = {
            node_id for node_id, node in nodes.items()
            if abs(node["x_um"] - x_mid) <= halfwidth
        }
    return groups


def terminal_rows(balance: dict[str, Any], terminal_csv: Path) -> list[dict[str, Any]]:
    sentaurus = balance["sentaurus_current_reference"]["final_coupled"]["contacts"]
    vela_rows = {row["contact"]: row for row in read_csv_rows(terminal_csv)}
    component_values = {
        contact: {
            "electron": float(row["current_electron_A_per_um"]),
            "hole": float(row["current_hole_A_per_um"]),
        }
        for contact, row in vela_rows.items()
    }

    def pair_values(component: str) -> tuple[dict[str, float], float, str]:
        raw = {
            contact: values[component]
            for contact, values in component_values.items()
            if contact in {"Anode", "Cathode"}
        }
        clamped = {
            contact: 0.0 if abs(value) <= COMPONENT_ZERO_FLOOR_A_PER_UM else value
            for contact, value in raw.items()
        }
        raw_sum = sum(raw.values())
        clamped_sum = sum(clamped.values())
        max_abs = max((abs(value) for value in clamped.values()), default=0.0)
        if abs(raw_sum) <= COMPONENT_ZERO_FLOOR_A_PER_UM:
            status = "pass_absolute_floor"
        elif max_abs > 0.0 and abs(clamped_sum) / max_abs <= 0.05:
            status = "pass_relative"
        else:
            status = "fail"
        return clamped, clamped_sum, status

    electron_pair, electron_pair_sum, electron_pair_status = pair_values("electron")
    hole_pair, hole_pair_sum, hole_pair_status = pair_values("hole")
    rows = []
    for contact in ["Anode", "Cathode"]:
        vrow = vela_rows[contact]
        srow = sentaurus[contact]
        sentaurus_total = float(srow["total_current"])
        vela_total_a_per_um = float(vrow["current_total_A_per_um"])
        vela_total_a_per_m = float(vrow["current_total"])
        same_sign = (
            (sentaurus_total == 0.0 and vela_total_a_per_um == 0.0)
            or (sentaurus_total > 0.0 and vela_total_a_per_um > 0.0)
            or (sentaurus_total < 0.0 and vela_total_a_per_um < 0.0)
        )
        conversion_candidates = {
            "sentaurus_raw_vs_vela_A_per_um": sentaurus_total / vela_total_a_per_um,
            "sentaurus_A_per_um_vs_vela_A_per_um": sentaurus_total / vela_total_a_per_um,
            "sentaurus_A_per_cm_width_to_A_per_um": (sentaurus_total / 1.0e4) / vela_total_a_per_um,
            "sentaurus_A_per_m_width_to_A_per_um": (sentaurus_total / 1.0e6) / vela_total_a_per_um,
            "vela_A_per_m_vs_sentaurus": vela_total_a_per_m / sentaurus_total,
        }
        rows.append({
            "contact": contact,
            "vela_electron_A_per_um": vrow["current_electron_A_per_um"],
            "vela_hole_A_per_um": vrow["current_hole_A_per_um"],
            "vela_total_A_per_um": vrow["current_total_A_per_um"],
            "vela_total_A_per_m": vrow["current_total"],
            "component_zero_floor_A_per_um": COMPONENT_ZERO_FLOOR_A_PER_UM,
            "vela_electron_pair_A_per_um": electron_pair.get(contact, 0.0),
            "vela_hole_pair_A_per_um": hole_pair.get(contact, 0.0),
            "vela_electron_pair_sum_A_per_um": electron_pair_sum,
            "vela_hole_pair_sum_A_per_um": hole_pair_sum,
            "vela_electron_pair_status": electron_pair_status,
            "vela_hole_pair_status": hole_pair_status,
            "sentaurus_electron_current": srow["electron_current"],
            "sentaurus_hole_current": srow["hole_current"],
            "sentaurus_total_current": srow["total_current"],
            "abs_ratio_sentaurus_to_vela": abs(sentaurus_total) / max(abs(vela_total_a_per_um), 1.0e-300),
            "sign_relation": "same_sign" if same_sign else "opposite_sign",
            **conversion_candidates,
        })
    return rows


def edge_summary(edge_csv: Path) -> dict[str, dict[str, Any]]:
    rows = read_csv_rows(edge_csv)
    by_contact: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_contact.setdefault(row["current_contact"], []).append(row)
    summary = {}
    for contact, items in by_contact.items():
        total_a_per_m = sum(as_float(row["current_total"]) for row in items)
        electron_a_per_m = sum(as_float(row["current_electron"]) for row in items)
        hole_a_per_m = sum(as_float(row["current_hole"]) for row in items)
        total_length_m = sum(as_float(row["edge_length_m"]) for row in items)
        def density_values(column: str) -> list[float]:
            return [
                as_float(row[column]) / max(as_float(row["edge_length_m"]), 1.0e-300) / 1.0e4
                for row in items
            ]
        summary[contact] = {
            "edges": len(items),
            "total_edge_length_m": total_length_m,
            "electron_current_A_per_m": electron_a_per_m,
            "hole_current_A_per_m": hole_a_per_m,
            "total_current_A_per_m": total_a_per_m,
            "total_current_A_per_um": total_a_per_m / 1.0e6,
            "electron_normal_current_density_A_per_cm2": stats(density_values("current_electron")),
            "hole_normal_current_density_A_per_cm2": stats(density_values("current_hole")),
            "total_normal_current_density_A_per_cm2": stats(density_values("current_total")),
            "electron_drift_current_A_per_m": sum(as_float(row["current_electron_drift"]) for row in items),
            "electron_diffusion_current_A_per_m": sum(as_float(row["current_electron_diffusion"]) for row in items),
            "hole_drift_current_A_per_m": sum(as_float(row["current_hole_drift"]) for row in items),
            "hole_diffusion_current_A_per_m": sum(as_float(row["current_hole_diffusion"]) for row in items),
        }
    return summary


def current_density_rows(fields_root: Path,
                         contacts: dict[str, list[int]],
                         edge_stats: dict[str, dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    sentaurus_fields = {
        "electron": read_sentaurus_scalar(fields_root, "eCurrentDensity"),
        "hole": read_sentaurus_scalar(fields_root, "hCurrentDensity"),
        "total": read_sentaurus_scalar(fields_root, "TotalCurrentDensity"),
    }
    rows = []
    summary: dict[str, Any] = {}
    for contact, ids in contacts.items():
        contact_summary: dict[str, Any] = {}
        for carrier, field in sentaurus_fields.items():
            values = [field[node_id] for node_id in ids if node_id in field]
            abs_values = [abs(value) for value in values]
            contact_summary[f"sentaurus_{carrier}_current_density_A_per_cm2"] = stats(values)
            contact_summary[f"sentaurus_abs_{carrier}_current_density_A_per_cm2"] = stats(abs_values)
            for node_id in ids:
                if node_id in field:
                    rows.append({
                        "contact": contact,
                        "node_id": node_id,
                        "quantity": f"sentaurus_{carrier}_current_density_A_per_cm2",
                        "value": field[node_id],
                    })
        contact_summary["vela_edge_current_density"] = edge_stats.get(contact, {})
        summary[contact] = contact_summary
    return rows, summary


def safe_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or abs(denominator) <= 1.0e-300:
        return None
    return numerator / denominator


def contact_current_driver_rows(
    fields_root: Path,
    nodes: dict[int, dict[str, float]],
    contacts: dict[str, list[int]],
    vtk_fields: dict[str, list[float]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    sentaurus = {
        "potential": read_sentaurus_scalar(fields_root, "ElectrostaticPotential"),
        "electron_qf": read_sentaurus_scalar(fields_root, "eQuasiFermiPotential"),
        "hole_qf": read_sentaurus_scalar(fields_root, "hQuasiFermiPotential"),
        "electron_density": read_sentaurus_scalar(fields_root, "eDensity"),
        "hole_density": read_sentaurus_scalar(fields_root, "hDensity"),
        "electron_current_density": read_sentaurus_scalar(fields_root, "eCurrentDensity"),
        "hole_current_density": read_sentaurus_scalar(fields_root, "hCurrentDensity"),
        "total_current_density": read_sentaurus_scalar(fields_root, "TotalCurrentDensity"),
    }
    vela = {
        "potential": {idx: value for idx, value in enumerate(vtk_fields["Potential"])},
        "electron_qf": {idx: value for idx, value in enumerate(vtk_fields["ElectronQuasiFermi"])},
        "hole_qf": {idx: value for idx, value in enumerate(vtk_fields["HoleQuasiFermi"])},
        "electron_density": {idx: value * 1.0e-6 for idx, value in enumerate(vtk_fields["Electrons"])},
        "hole_density": {idx: value * 1.0e-6 for idx, value in enumerate(vtk_fields["Holes"])},
    }
    rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {}
    for contact, ids in contacts.items():
        contact_rows: list[dict[str, Any]] = []
        for node_id in ids:
            row = {
                "contact": contact,
                "node_id": node_id,
                "x_um": nodes[node_id]["x_um"],
                "y_um": nodes[node_id]["y_um"],
                "sentaurus_potential_V": sentaurus["potential"].get(node_id),
                "vela_potential_V": vela["potential"].get(node_id),
                "potential_diff_V": (
                    vela["potential"].get(node_id, 0.0) - sentaurus["potential"].get(node_id, 0.0)
                ),
                "sentaurus_eQF_V": sentaurus["electron_qf"].get(node_id),
                "vela_eQF_V": vela["electron_qf"].get(node_id),
                "eQF_diff_V": (
                    vela["electron_qf"].get(node_id, 0.0) - sentaurus["electron_qf"].get(node_id, 0.0)
                ),
                "sentaurus_hQF_V": sentaurus["hole_qf"].get(node_id),
                "vela_hQF_V": vela["hole_qf"].get(node_id),
                "hQF_diff_V": (
                    vela["hole_qf"].get(node_id, 0.0) - sentaurus["hole_qf"].get(node_id, 0.0)
                ),
                "sentaurus_eDensity_cm3": sentaurus["electron_density"].get(node_id),
                "vela_eDensity_cm3": vela["electron_density"].get(node_id),
                "eDensity_ratio_vela_to_sentaurus": safe_ratio(
                    vela["electron_density"].get(node_id), sentaurus["electron_density"].get(node_id)
                ),
                "sentaurus_hDensity_cm3": sentaurus["hole_density"].get(node_id),
                "vela_hDensity_cm3": vela["hole_density"].get(node_id),
                "hDensity_ratio_vela_to_sentaurus": safe_ratio(
                    vela["hole_density"].get(node_id), sentaurus["hole_density"].get(node_id)
                ),
                "sentaurus_eCurrentDensity_A_cm2": sentaurus["electron_current_density"].get(node_id),
                "sentaurus_hCurrentDensity_A_cm2": sentaurus["hole_current_density"].get(node_id),
                "sentaurus_totalCurrentDensity_A_cm2": sentaurus["total_current_density"].get(node_id),
            }
            contact_rows.append(row)
            rows.append(row)
        summary[contact] = {
            "nodes": len(contact_rows),
            "max_abs_eQF_diff_V": max((abs(row["eQF_diff_V"]) for row in contact_rows), default=0.0),
            "max_abs_hQF_diff_V": max((abs(row["hQF_diff_V"]) for row in contact_rows), default=0.0),
            "mean_eDensity_ratio_vela_to_sentaurus": statistics.mean([
                row["eDensity_ratio_vela_to_sentaurus"] for row in contact_rows
                if row["eDensity_ratio_vela_to_sentaurus"] is not None
            ]),
            "mean_hDensity_ratio_vela_to_sentaurus": statistics.mean([
                row["hDensity_ratio_vela_to_sentaurus"] for row in contact_rows
                if row["hDensity_ratio_vela_to_sentaurus"] is not None
            ]),
        }
    return rows, summary


def contact_edge_driver_rows(
    edge_csv: Path,
    elements: list[list[int]],
    contacts: dict[str, list[int]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    triangle_lookup: dict[tuple[int, int], list[int]] = {}
    for tri in elements:
        a, b, c = tri
        for u, v, w in ((a, b, c), (b, c, a), (c, a, b)):
            triangle_lookup.setdefault(tuple(sorted((u, v))), []).append(w)

    rows: list[dict[str, Any]] = []
    by_contact: dict[str, list[dict[str, Any]]] = {}
    for row in read_csv_rows(edge_csv):
        contact = row["current_contact"]
        contact_nodes = set(contacts.get(contact, []))
        node0 = int(row["node0"])
        node1 = int(row["node1"])
        on_contact = [node for node in (node0, node1) if node in contact_nodes]
        off_contact = [node for node in (node0, node1) if node not in contact_nodes]
        inferred_interior = ""
        if len(off_contact) == 1:
            inferred_interior = off_contact[0]
        elif len(on_contact) == 2:
            candidates = [
                node for node in triangle_lookup.get(tuple(sorted((node0, node1))), [])
                if node not in contact_nodes
            ]
            inferred_interior = candidates[0] if candidates else ""
        contact_node0 = on_contact[0] if len(on_contact) >= 1 else ""
        contact_node1 = on_contact[1] if len(on_contact) >= 2 else ""

        def avg_for(column0: str, column1: str) -> float:
            values = []
            if node0 in contact_nodes:
                values.append(as_float(row[column0]))
            if node1 in contact_nodes:
                values.append(as_float(row[column1]))
            return sum(values) / len(values) if values else 0.5 * (as_float(row[column0]) + as_float(row[column1]))

        def interior_for(column0: str, column1: str) -> float | None:
            if inferred_interior == node0:
                return as_float(row[column0])
            if inferred_interior == node1:
                return as_float(row[column1])
            return None

        out = {
            "contact": contact,
            "edge_id": row["edge_id"],
            "contact_node0": contact_node0,
            "contact_node1": contact_node1,
            "interior_node": inferred_interior,
            "psi_contact_avg": avg_for("psi0", "psi1"),
            "psi_interior": interior_for("psi0", "psi1"),
            "phin_contact_avg": avg_for("phin0", "phin1"),
            "phin_interior": interior_for("phin0", "phin1"),
            "phip_contact_avg": avg_for("phip0", "phip1"),
            "phip_interior": interior_for("phip0", "phip1"),
            "n_contact_avg": avg_for("n0", "n1"),
            "n_interior": interior_for("n0", "n1"),
            "p_contact_avg": avg_for("p0", "p1"),
            "p_interior": interior_for("p0", "p1"),
            "vela_current_electron_A_per_m": row["current_electron"],
            "vela_current_hole_A_per_m": row["current_hole"],
            "vela_current_total_A_per_m": row["current_total"],
        }
        out["abs_phin_delta_contact_to_interior_V"] = (
            abs(as_float(out["phin_interior"]) - as_float(out["phin_contact_avg"]))
            if out["phin_interior"] not in (None, "") else ""
        )
        out["abs_phip_delta_contact_to_interior_V"] = (
            abs(as_float(out["phip_interior"]) - as_float(out["phip_contact_avg"]))
            if out["phip_interior"] not in (None, "") else ""
        )
        rows.append(out)
        by_contact.setdefault(contact, []).append(out)

    summary = {}
    for contact, items in by_contact.items():
        phin = [
            as_float(row["abs_phin_delta_contact_to_interior_V"])
            for row in items if row["abs_phin_delta_contact_to_interior_V"] not in ("", None)
        ]
        phip = [
            as_float(row["abs_phip_delta_contact_to_interior_V"])
            for row in items if row["abs_phip_delta_contact_to_interior_V"] not in ("", None)
        ]
        currents = [abs(as_float(row["vela_current_total_A_per_m"])) for row in items]
        summary[contact] = {
            "edges": len(items),
            "max_abs_phin_delta_contact_to_interior_V": max(phin, default=0.0),
            "max_abs_phip_delta_contact_to_interior_V": max(phip, default=0.0),
            "mean_abs_total_current_A_per_m": sum(currents) / len(currents) if currents else 0.0,
            "max_abs_total_current_A_per_m": max(currents, default=0.0),
        }
    return rows, summary


def state_comparisons(fields_root: Path,
                      nodes: dict[int, dict[str, float]],
                      groups: dict[str, set[int]],
                      vtk_fields: dict[str, list[float]],
                      vela_field_v_per_cm: dict[int, float]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    specs = [
        ("potential_V", "ElectrostaticPotential", "Potential", 1.0),
        ("electron_qf_V", "eQuasiFermiPotential", "ElectronQuasiFermi", 1.0),
        ("hole_qf_V", "hQuasiFermiPotential", "HoleQuasiFermi", 1.0),
        ("electron_density_cm3", "eDensity", "Electrons", 1.0e-6),
        ("hole_density_cm3", "hDensity", "Holes", 1.0e-6),
        ("net_doping_cm3", "DopingConcentration", "NetDoping", 1.0e-6),
    ]
    comparisons: dict[str, Any] = {}
    rows_by_node: dict[int, dict[str, Any]] = {
        node_id: {"node_id": node_id, "x_um": node["x_um"], "y_um": node["y_um"]}
        for node_id, node in nodes.items()
    }
    for name, sentaurus_field, vela_field, vela_scale in specs:
        sentaurus = read_sentaurus_scalar(fields_root, sentaurus_field)
        vela = {idx: value * vela_scale for idx, value in enumerate(vtk_fields[vela_field])}
        comparisons[name] = compare_maps(sentaurus, vela, groups)
        for node_id in set(sentaurus) & set(vela):
            row = rows_by_node[node_id]
            row[f"{name}_sentaurus"] = sentaurus[node_id]
            row[f"{name}_vela"] = vela[node_id]
            row[f"{name}_diff"] = vela[node_id] - sentaurus[node_id]
    sentaurus_field_v_per_cm = read_sentaurus_scalar(fields_root, "ElectricField")
    comparisons["electric_field_V_per_cm"] = compare_maps(
        sentaurus_field_v_per_cm, vela_field_v_per_cm, groups)
    for node_id in set(sentaurus_field_v_per_cm) & set(vela_field_v_per_cm):
        row = rows_by_node[node_id]
        row["electric_field_V_per_cm_sentaurus"] = sentaurus_field_v_per_cm[node_id]
        row["electric_field_V_per_cm_vela"] = vela_field_v_per_cm[node_id]
        row["electric_field_V_per_cm_diff"] = vela_field_v_per_cm[node_id] - sentaurus_field_v_per_cm[node_id]
    return comparisons, [rows_by_node[node_id] for node_id in sorted(rows_by_node)]


def qf_profile_rows(nodes: dict[int, dict[str, float]],
                    contacts: dict[str, list[int]],
                    fields_root: Path,
                    vtk_fields: dict[str, list[float]]) -> list[dict[str, Any]]:
    sent_e = read_sentaurus_scalar(fields_root, "eQuasiFermiPotential")
    sent_h = read_sentaurus_scalar(fields_root, "hQuasiFermiPotential")
    vela_e = {idx: value for idx, value in enumerate(vtk_fields["ElectronQuasiFermi"])}
    vela_h = {idx: value for idx, value in enumerate(vtk_fields["HoleQuasiFermi"])}
    rows = []
    for contact, ids in contacts.items():
        for node_id in sorted(ids, key=lambda nid: (nodes[nid]["y_um"], nodes[nid]["x_um"])):
            rows.append({
                "contact": contact,
                "node_id": node_id,
                "x_um": nodes[node_id]["x_um"],
                "y_um": nodes[node_id]["y_um"],
                "sentaurus_electron_qf_V": sent_e.get(node_id),
                "vela_electron_qf_V": vela_e.get(node_id),
                "electron_qf_diff_V": vela_e.get(node_id, 0.0) - sent_e.get(node_id, 0.0),
                "sentaurus_hole_qf_V": sent_h.get(node_id),
                "vela_hole_qf_V": vela_h.get(node_id),
                "hole_qf_diff_V": vela_h.get(node_id, 0.0) - sent_h.get(node_id, 0.0),
            })
    return rows


def try_plot(output_dir: Path,
             terminal: list[dict[str, Any]],
             flux_rows: list[dict[str, Any]],
             current_density_summary: dict[str, Any],
             qf_rows: list[dict[str, Any]],
             state_summary: dict[str, Any]) -> dict[str, str]:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return {}

    plots: dict[str, str] = {}
    font = ImageFont.load_default()

    def canvas(title: str) -> tuple[Any, Any]:
        image = Image.new("RGB", (1200, 760), "white")
        draw = ImageDraw.Draw(image)
        draw.text((40, 24), title, fill=(20, 24, 32), font=font)
        return image, draw

    def save(image: Any, name: str) -> None:
        path = output_dir / name
        image.save(path)
        plots[name[:-4]] = str(path.resolve())

    def log_height(value: float, scale_min: float, scale_max: float, height: int) -> int:
        value = max(abs(value), scale_min)
        lo = math.log10(scale_min)
        hi = math.log10(scale_max)
        return int(height * (math.log10(value) - lo) / max(hi - lo, 1.0e-12))

    image, draw = canvas("Terminal and contact-flux current magnitude")
    items = []
    for row in terminal:
        items.append((f"{row['contact']} Vela", abs(float(row["vela_total_A_per_um"])), "#2563eb"))
        items.append((f"{row['contact']} Sentaurus", abs(float(row["sentaurus_total_current"])), "#dc2626"))
    for row in flux_rows:
        items.append((f"{row['contact']} TDR flux", abs(float(row["sentaurus_contact_current_flux_A"])), "#059669"))
    max_val = max(value for _, value, _ in items)
    min_val = max(min(value for _, value, _ in items if value > 0.0), 1.0e-30)
    colors = {"#2563eb": (37, 99, 235), "#dc2626": (220, 38, 38), "#059669": (5, 150, 105)}
    x0, y0, w, h = 80, 100, 1040, 480
    draw.line((x0, y0 + h, x0 + w, y0 + h), fill=(30, 41, 59))
    bar_w = max(24, int(w / max(len(items), 1) * 0.58))
    step = w / max(len(items), 1)
    for idx, (label, value, color) in enumerate(items):
        bh = log_height(value, min_val, max_val, h)
        x = int(x0 + idx * step + 0.2 * step)
        draw.rectangle((x, y0 + h - bh, x + bar_w, y0 + h), fill=colors[color])
        draw.text((x - 8, y0 + h + 12), label[:16], fill=(15, 23, 42), font=font)
        draw.text((x - 8, y0 + h - bh - 16), f"{value:.2e}", fill=(15, 23, 42), font=font)
    draw.text((80, 650), "Log-scale magnitudes; Vela terminal values are A/um, Sentaurus/TDR values use exported units.", fill=(71, 85, 105), font=font)
    save(image, "terminal_contact_currents.png")

    image, draw = canvas("Contact current-density statistics")
    cd_items = []
    for contact, summary in current_density_summary.items():
        sent_total = summary["sentaurus_abs_total_current_density_A_per_cm2"]["mean"] or 0.0
        vela_total = summary["vela_edge_current_density"]["total_normal_current_density_A_per_cm2"]["mean"] or 0.0
        cd_items.append((f"{contact} Sentaurus |J|", sent_total, (220, 38, 38)))
        cd_items.append((f"{contact} Vela edge Jn", abs(vela_total), (37, 99, 235)))
    max_val = max((value for _, value, _ in cd_items), default=1.0)
    min_val = max(min((value for _, value, _ in cd_items if value > 0.0), default=1.0e-30), 1.0e-30)
    for idx, (label, value, color) in enumerate(cd_items):
        bh = log_height(value, min_val, max_val, 480)
        x = 120 + idx * 220
        draw.rectangle((x, 580 - bh, x + 95, 580), fill=color)
        draw.text((x - 20, 598), label[:24], fill=(15, 23, 42), font=font)
        draw.text((x - 8, 560 - bh), f"{value:.2e}", fill=(15, 23, 42), font=font)
    draw.text((80, 650), "Sentaurus uses contact-node |TotalCurrentDensity|; Vela uses edge normal current / edge length.", fill=(71, 85, 105), font=font)
    save(image, "contact_current_density_stats.png")

    image, draw = canvas("Contact quasi-Fermi potential profiles")
    panel_w = 520
    colors_qf = {
        "sentaurus_electron_qf_V": (220, 38, 38),
        "vela_electron_qf_V": (37, 99, 235),
        "sentaurus_hole_qf_V": (245, 158, 11),
        "vela_hole_qf_V": (5, 150, 105),
    }
    all_qf = [
        float(row[key]) for row in qf_rows
        for key in colors_qf
        if row.get(key) not in (None, "")
    ]
    qmin, qmax = min(all_qf), max(all_qf)
    for panel_idx, contact in enumerate(["Anode", "Cathode"]):
        rows = [row for row in qf_rows if row["contact"] == contact]
        ys = [float(row["y_um"]) for row in rows]
        ymin, ymax = min(ys), max(ys)
        left = 60 + panel_idx * 580
        top = 100
        draw.rectangle((left, top, left + panel_w, top + 500), outline=(203, 213, 225))
        draw.text((left, top - 22), contact, fill=(15, 23, 42), font=font)
        for row in rows:
            x_base = left + int((float(row["y_um"]) - ymin) / max(ymax - ymin, 1.0e-12) * panel_w)
            for key, color in colors_qf.items():
                y_plot = top + 500 - int((float(row[key]) - qmin) / max(qmax - qmin, 1.0e-12) * 500)
                draw.ellipse((x_base - 2, y_plot - 2, x_base + 2, y_plot + 2), fill=color)
    draw.text((80, 650), "Red/blue: electron QF Sentaurus/Vela. Orange/green: hole QF Sentaurus/Vela.", fill=(71, 85, 105), font=font)
    save(image, "contact_qf_profiles.png")

    image, draw = canvas("Node-state mean relative error")
    state_items = [
        (name, summary.get("mean_rel_diff") or 0.0)
        for name, summary in state_summary.items()
    ]
    max_val = max((value for _, value in state_items), default=1.0)
    min_val = max(min((value for _, value in state_items if value > 0.0), default=1.0e-18), 1.0e-18)
    for idx, (label, value) in enumerate(state_items):
        bh = log_height(value, min_val, max_val, 460)
        x = 80 + idx * 145
        draw.rectangle((x, 570 - bh, x + 70, 570), fill=(79, 70, 229))
        draw.text((x - 18, 592), label[:18], fill=(15, 23, 42), font=font)
        draw.text((x - 4, 550 - bh), f"{value:.2e}", fill=(15, 23, 42), font=font)
    draw.text((80, 650), "Log-scale mean relative error over common nodes.", fill=(71, 85, 105), font=font)
    save(image, "state_error_summary.png")
    return plots


def write_markdown(path: Path, report: dict[str, Any]) -> None:
    lines = [
        "# PN2D 0V current-related quantity comparison",
        "",
        "## Outputs",
    ]
    for name, value in report["outputs"].items():
        lines.append(f"- {name}: `{value}`")
    lines.extend(["", "## Key terminal currents", ""])
    lines.append("| contact | Vela total A/um | Sentaurus total | ratio S/V | sign |")
    lines.append("|---|---:|---:|---:|---|")
    for row in report["terminal_currents"]["rows"]:
        lines.append(
            f"| {row['contact']} | {float(row['vela_total_A_per_um']):.6e} | "
            f"{float(row['sentaurus_total_current']):.6e} | "
            f"{float(row['abs_ratio_sentaurus_to_vela']):.6e} | {row['sign_relation']} |"
        )
    lines.extend([
        "",
        "Vela carrier-component rows keep the raw terminal values, but values below "
        f"`{COMPONENT_ZERO_FLOOR_A_PER_UM:.1e} A/um` are also reported with a zero-floor "
        "view for pair-balance interpretation.",
        "",
        "| contact | raw Vela electron A/um | raw Vela hole A/um | pair electron A/um | pair hole A/um | electron pair status | hole pair status |",
        "|---|---:|---:|---:|---:|---|---|",
    ])
    for row in report["terminal_currents"]["rows"]:
        lines.append(
            f"| {row['contact']} | {float(row['vela_electron_A_per_um']):.6e} | "
            f"{float(row['vela_hole_A_per_um']):.6e} | "
            f"{float(row['vela_electron_pair_A_per_um']):.6e} | "
            f"{float(row['vela_hole_pair_A_per_um']):.6e} | "
            f"{row['vela_electron_pair_status']} | {row['vela_hole_pair_status']} |"
        )
    lines.extend(["", "## Node-state comparison", ""])
    lines.append("| quantity | points | mean abs diff | mean rel diff | max abs diff |")
    lines.append("|---|---:|---:|---:|---:|")
    for name, summary in report["node_state_comparison"].items():
        lines.append(
            f"| {name} | {summary['points_compared']} | "
            f"{summary['mean_abs_diff']:.6e} | {summary['mean_rel_diff']:.6e} | "
            f"{summary['max_abs_diff']:.6e} |"
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", type=Path,
                        default=Path("build/reference_tcad/pn2d_sentaurus2018"))
    parser.add_argument("--vtk", type=Path,
                        default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/0v_qf_drivers/baseline/baseline_0000_0V.vtk"))
    parser.add_argument("--current-balance", type=Path,
                        default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/0v_current_balance/pn2d_0v_current_balance.json"))
    parser.add_argument("--terminal-csv", type=Path,
                        default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/0v_current_balance/probe/terminal_balance.csv"))
    parser.add_argument("--edge-csv", type=Path,
                        default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/0v_current_balance/probe/contact_edges.csv"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("build/reference_tcad/pn2d_sentaurus2018/reports/0v_current_related"))
    args = parser.parse_args()

    reference_root = args.reference_root
    sim_root = reference_root / "sim_fields" / "0v"
    fields_root = sim_root / "fields"
    manifest_path = sim_root / "field_manifest.json"
    nodes = read_nodes(sim_root / "nodes.csv")
    contacts = read_contacts(sim_root / "contacts.csv")
    groups = make_groups(nodes, contacts)
    points, vtk_fields = parse_vtk(args.vtk)
    elements = read_elements(sim_root / "elements.csv")
    vela_field_v_per_cm = {
        node_id: value * 1.0e-2
        for node_id, value in derive_nodal_gradient(points, vtk_fields["Potential"], elements).items()
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    balance = read_json(args.current_balance)
    terminal = terminal_rows(balance, args.terminal_csv)
    terminal_csv_out = args.output_dir / "terminal_currents_compare.csv"
    write_csv(terminal_csv_out, terminal, list(terminal[0]))

    flux = contact_flux_by_name(fields_root, manifest_path)
    boundary_rows, boundary_summary = sentaurus_current_density_boundary_integrals(
        fields_root, nodes, elements, contacts)
    flux_rows = []
    for row in terminal:
        contact = row["contact"]
        total_integral = boundary_summary.get(contact, {}).get("by_quantity", {}).get("total", {})
        flux_rows.append({
            "contact": contact,
            "sentaurus_contact_current_flux_A": flux.get(contact),
            "sentaurus_plt_total_current": row["sentaurus_total_current"],
            "sentaurus_total_current_density_signed_A_per_cm_width": total_integral.get("signed_sum_A_per_cm_width"),
            "sentaurus_total_current_density_abs_A_per_cm_width": total_integral.get("abs_sum_A_per_cm_width"),
            "vela_total_A_per_um": row["vela_total_A_per_um"],
        })
    flux_csv_out = args.output_dir / "contact_flux_compare.csv"
    write_csv(flux_csv_out, flux_rows, list(flux_rows[0]))
    boundary_csv_out = args.output_dir / "sentaurus_boundary_current_density_integral.csv"
    write_csv(
        boundary_csv_out,
        boundary_rows,
        [
            "contact",
            "quantity",
            "node0",
            "node1",
            "edge_length_um",
            "edge_length_cm",
            "j0_A_per_cm2",
            "j1_A_per_cm2",
            "avg_j_A_per_cm2",
            "edge_current_A_per_cm_width",
            "direction_status",
        ],
    )

    edge_stats = edge_summary(args.edge_csv)
    cd_rows, cd_summary = current_density_rows(fields_root, contacts, edge_stats)
    current_density_csv_out = args.output_dir / "contact_current_density_compare.csv"
    write_csv(current_density_csv_out, cd_rows, ["contact", "node_id", "quantity", "value"])

    driver_rows, driver_summary = contact_current_driver_rows(
        fields_root, nodes, contacts, vtk_fields)
    driver_csv_out = args.output_dir / "contact_current_driver_nodes.csv"
    write_csv(driver_csv_out, driver_rows, list(driver_rows[0]))

    edge_driver_rows, edge_driver_summary = contact_edge_driver_rows(
        args.edge_csv, elements, contacts)
    edge_driver_csv_out = args.output_dir / "contact_edge_driver_nodes.csv"
    write_csv(edge_driver_csv_out, edge_driver_rows, list(edge_driver_rows[0]))

    qf_rows = qf_profile_rows(nodes, contacts, fields_root, vtk_fields)
    qf_csv_out = args.output_dir / "qf_contact_profiles.csv"
    write_csv(qf_csv_out, qf_rows, list(qf_rows[0]))

    state_summary, state_rows = state_comparisons(fields_root, nodes, groups, vtk_fields, vela_field_v_per_cm)
    state_csv_out = args.output_dir / "state_node_comparison.csv"
    state_fieldnames = sorted({key for row in state_rows for key in row})
    state_fieldnames = ["node_id", "x_um", "y_um"] + [
        name for name in state_fieldnames if name not in {"node_id", "x_um", "y_um"}
    ]
    write_csv(state_csv_out, state_rows, state_fieldnames)

    plots = try_plot(args.output_dir, terminal, flux_rows, cd_summary, qf_rows, state_summary)
    report = {
        "status": "generated",
        "inputs": {
            "reference_root": str(reference_root),
            "vtk": str(args.vtk),
            "current_balance": str(args.current_balance),
            "terminal_csv": str(args.terminal_csv),
            "edge_csv": str(args.edge_csv),
        },
        "outputs": {
            "terminal_currents_csv": str(terminal_csv_out.resolve()),
            "contact_flux_csv": str(flux_csv_out.resolve()),
            "sentaurus_boundary_current_density_integral_csv": str(boundary_csv_out.resolve()),
            "contact_current_density_csv": str(current_density_csv_out.resolve()),
            "contact_current_driver_nodes_csv": str(driver_csv_out.resolve()),
            "contact_edge_driver_nodes_csv": str(edge_driver_csv_out.resolve()),
            "qf_contact_profiles_csv": str(qf_csv_out.resolve()),
            "state_node_comparison_csv": str(state_csv_out.resolve()),
            **{f"{name}_png": path for name, path in plots.items()},
        },
        "terminal_currents": {
            "rows": terminal,
            "conservation_summary": balance.get("conservation_summary", {}),
        },
        "contact_flux_comparison": {
            "rows": flux_rows,
            "note": "Sentaurus ContactCurrentFlux is exported from TDR contact regions and is kept in its native A unit.",
        },
        "sentaurus_boundary_current_density_integral": boundary_summary,
        "contact_current_density_comparison": cd_summary,
        "contact_current_driver_nodes": driver_summary,
        "contact_edge_driver_nodes": edge_driver_summary,
        "qf_contact_profiles": {
            "rows": len(qf_rows),
            "contacts": {name: len(ids) for name, ids in contacts.items()},
        },
        "node_state_comparison": state_summary,
    }
    report_path = args.output_dir / REPORT_NAME
    write_json(report_path, report)
    write_markdown(args.output_dir / MARKDOWN_NAME, report)
    print(report_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
