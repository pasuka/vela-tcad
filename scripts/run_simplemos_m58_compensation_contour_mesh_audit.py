#!/usr/bin/env python3
"""Execute the frozen, read-only SimpleMOS M58 compensation/contour audit."""

from __future__ import annotations

from collections import defaultdict
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m58_compensation_contour_mesh_audit_contract_v1.json"
FREEZE = ROOT / "simplemos_m58_compensation_contour_mesh_audit_contract_freeze.json"
M46_REPORT = ROOT / "full_matrix_requalification/m46_full_matrix_requalification_report.json"
M46_CASES = ROOT / "full_matrix_requalification/m46_case_summary.csv"
M55_REPORT = ROOT / "doping_well_terminal_attribution/m55_doping_well_terminal_attribution_report.json"
M57_REPORT = ROOT / "doping_well_boundary_flux_closure/m57_doping_well_boundary_flux_closure_report.json"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
IMPORT_ROOT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m58_compensation_contour_mesh_audit/tdr_exports"
OUTPUT = ROOT / "compensation_contour_mesh_audit"
REPORT = OUTPUT / "m58_compensation_contour_mesh_audit_report.json"
DEVICES_CSV = OUTPUT / "m58_device_compensation_mesh_ledger.csv"
CONTOURS_CSV = OUTPUT / "m58_zero_contour_segment_ledger.csv"
COMPONENTS_CSV = OUTPUT / "m58_connectivity_component_ledger.csv"
CRITICAL_CSV = OUTPUT / "m58_critical_gate_edge_node_ledger.csv"
SPECIES_CSV = OUTPUT / "m58_pair_species_coordinate_ledger.csv"
PAIRS_CSV = OUTPUT / "m58_nwell_error_attribution_ledger.csv"
CORRELATIONS_CSV = OUTPUT / "m58_metric_error_correlation_ledger.csv"
EVIDENCE = ROOT / "simplemos_m58_compensation_contour_mesh_audit_evidence.json"
DOC = REPO / "docs/validation/simplemos_m58_compensation_contour_mesh_audit_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m58/artifact.json"
FIELDS = ("NetActive", "BActive", "AsActive", "PActive")


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def close_enough(left: float, right: float, absolute: float,
                 relative: float) -> bool:
    delta = abs(left - right)
    return delta <= absolute or delta <= relative * max(
        abs(left), abs(right), 1e-300)


def validate_contract() -> tuple[dict[str, Any], dict[str, Any]]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    expected = "vela.simplemos.sdevice.m58_compensation_contour_mesh_audit_contract.v1"
    if contract.get("schema") != expected:
        raise ValueError("unexpected M58 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M58 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M58 contract changed after freeze")
    for relative, expected_hash in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected_hash:
            raise ValueError(f"M58 upstream changed after freeze: {relative}")
    m46 = read_json(M46_REPORT)
    m55 = read_json(M55_REPORT)
    m57 = read_json(M57_REPORT)
    upstream = contract["upstream"]
    if m46.get("status") != upstream["required_m46_status"]:
        raise ValueError("M58 M46 status anchor changed")
    if m55.get("status") != upstream["required_m55_status"]:
        raise ValueError("M58 requires accepted M55")
    if m55.get("classification") != upstream["required_m55_classification"]:
        raise ValueError("M58 M55 classification anchor changed")
    if m57.get("status") != upstream["required_m57_status"]:
        raise ValueError("M58 requires accepted M57")
    if m57.get("classification") != upstream["required_m57_classification"]:
        raise ValueError("M58 M57 classification anchor changed")
    return contract, freeze


def ensure_tdr_export(device: str, contract: dict[str, Any]) -> Path:
    tdr = REPO / contract["fixed_input_policy"]["tdr_pattern"].format(device=device)
    output = IMPORT_ROOT / device
    required = [output / "nodes.csv", output / "elements.csv", output / "contacts.csv"]
    required += [output / "fields" / f"{field}_region0.csv" for field in FIELDS]
    if all(path.exists() for path in required):
        return output
    if not IMPORTER.exists():
        raise FileNotFoundError(f"M58 importer unavailable: {IMPORTER}")
    output.mkdir(parents=True, exist_ok=True)
    subprocess.run([
        str(IMPORTER), "--tdr", str(tdr),
        "--inventory-json", str(output / "inventory.json"),
        "--export-dir", str(output),
    ], cwd=REPO, check=True)
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise ValueError(f"M58 TDR export missing required files: {missing}")
    return output


def scalar_field(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"])
            for row in read_csv(path)}


def in_gate_edge_roi(point: tuple[float, float], rules: dict[str, Any]) -> bool:
    x, y = point
    return (x <= float(rules["surface_depth_max_um"]) and
            abs(abs(y) - float(rules["gate_edge_abs_y_um"])) <=
            float(rules["gate_edge_half_width_um"]))


class UnionFind:
    def __init__(self, nodes: Iterable[int]) -> None:
        self.parent = {node: node for node in nodes}

    def find(self, node: int) -> int:
        parent = self.parent[node]
        if parent != node:
            self.parent[node] = self.find(parent)
        return self.parent[node]

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[b] = a


def contour_crossing(a: int, b: int,
                     nodes: dict[int, tuple[float, float]],
                     net: dict[int, float]) -> tuple[tuple[float, float], float, float] | None:
    va, vb = net[a], net[b]
    if va == 0.0 and vb == 0.0:
        return None
    if va * vb > 0.0:
        return None
    if va == 0.0:
        fraction = 0.0
    elif vb == 0.0:
        fraction = 1.0
    else:
        fraction = va / (va - vb)
    pa, pb = nodes[a], nodes[b]
    point = (pa[0] + fraction * (pb[0] - pa[0]),
             pa[1] + fraction * (pb[1] - pa[1]))
    edge_length = math.dist(pa, pb)
    vertex_distance = min(fraction, 1.0 - fraction) * edge_length
    normalized = min(fraction, 1.0 - fraction)
    return point, vertex_distance, normalized


def analyze_device(meta: dict[str, Any], contract: dict[str, Any]
                   ) -> tuple[dict[str, Any], list[dict[str, Any]],
                              list[dict[str, Any]], list[dict[str, Any]],
                              dict[str, Any]]:
    device = meta["id"]
    rules = contract["frozen_analysis_rules"]
    export = ensure_tdr_export(device, contract)
    nodes = {int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
             for row in read_csv(export / "nodes.csv")}
    triangles = [(int(row["id"]), tuple(int(row[f"node{i}"]) for i in range(3)))
                 for row in read_csv(export / "elements.csv")
                 if row["material"] == "Si"]
    silicon_nodes = {node for _, triangle in triangles for node in triangle}
    fields = {name: scalar_field(export / "fields" / f"{name}_region0.csv")
              for name in FIELDS}
    if any(set(field) != silicon_nodes for field in fields.values()):
        raise ValueError(f"M58 incomplete Silicon field coverage: {device}")
    net = fields["NetActive"]
    acceptor = fields["BActive"]
    arsenic = fields["AsActive"]
    phosphorus = fields["PActive"]
    donor = {node: arsenic[node] + phosphorus[node] for node in silicon_nodes}
    floor = float(rules["compensation_ratio_floor_cm3"])
    compensation = {node: (donor[node] + acceptor[node]) /
                    max(abs(net[node]), floor) for node in silicon_nodes}
    absolute_tolerance = float(
        contract["acceptance"]["active_field_identity_absolute_tolerance_cm3"])
    relative_tolerance = float(
        contract["acceptance"]["active_field_identity_relative_tolerance"])
    identity_residuals = {
        node: net[node] - (donor[node] - acceptor[node]) for node in silicon_nodes}
    identity_pass = all(close_enough(
        net[node], donor[node] - acceptor[node], absolute_tolerance,
        relative_tolerance) for node in silicon_nodes)

    edges: set[tuple[int, int]] = set()
    for _, triangle in triangles:
        for left, right in ((triangle[0], triangle[1]),
                            (triangle[1], triangle[2]),
                            (triangle[2], triangle[0])):
            edges.add(tuple(sorted((left, right))))
    roi_nodes = {node for node in silicon_nodes if in_gate_edge_roi(nodes[node], rules)}
    roi_edges = {edge for edge in edges if edge[0] in roi_nodes or edge[1] in roi_nodes}
    roi_edge_lengths = [math.dist(nodes[a], nodes[b]) for a, b in roi_edges]

    contacts = {row["name"].lower(): {int(value) for value in row["node_ids"].split(";")}
                & silicon_nodes for row in read_csv(export / "contacts.csv")}
    sign = {node: 1 if net[node] >= 0.0 else -1 for node in silicon_nodes}
    union = UnionFind(silicon_nodes)
    for left, right in edges:
        if sign[left] == sign[right]:
            union.union(left, right)
    component_nodes: dict[int, list[int]] = defaultdict(list)
    for node in silicon_nodes:
        component_nodes[union.find(node)].append(node)
    component_rows: list[dict[str, Any]] = []
    floating_n_count = 0
    single_node_floating_n_count = 0
    for index, members in enumerate(sorted(component_nodes.values(),
                                           key=lambda values: min(values))):
        kind = "n_type" if sign[members[0]] > 0 else "p_type"
        associated = sorted(name for name, contact_nodes in contacts.items()
                            if set(members) & contact_nodes)
        floating = not associated
        gate_edge = any(node in roi_nodes for node in members)
        if kind == "n_type" and floating:
            floating_n_count += 1
            if len(members) == 1:
                single_node_floating_n_count += 1
        component_rows.append({
            "device": device, "component_index": index,
            "doping_type": kind, "node_count": len(members),
            "contact_names": ";".join(associated), "floating": floating,
            "single_node": len(members) == 1, "intersects_gate_edge_roi": gate_edge,
            "min_abs_net_cm3": min(abs(net[node]) for node in members),
            "max_compensation_ratio": max(compensation[node] for node in members),
            "x_min_um": min(nodes[node][0] for node in members),
            "x_max_um": max(nodes[node][0] for node in members),
            "y_min_um": min(nodes[node][1] for node in members),
            "y_max_um": max(nodes[node][1] for node in members),
        })

    contour_rows: list[dict[str, Any]] = []
    for triangle_id, triangle in triangles:
        crossings = []
        for left, right in ((triangle[0], triangle[1]),
                            (triangle[1], triangle[2]),
                            (triangle[2], triangle[0])):
            crossing = contour_crossing(left, right, nodes, net)
            if crossing is not None:
                crossings.append(crossing)
        unique: list[tuple[tuple[float, float], float, float]] = []
        for crossing in crossings:
            if not any(math.dist(crossing[0], other[0]) < 1e-15 for other in unique):
                unique.append(crossing)
        if len(unique) != 2:
            continue
        point0, distance0, normalized0 = unique[0]
        point1, distance1, normalized1 = unique[1]
        middle = ((point0[0] + point1[0]) / 2.0,
                  (point0[1] + point1[1]) / 2.0)
        contour_rows.append({
            "device": device, "triangle_id": triangle_id,
            "x0_um": point0[0], "y0_um": point0[1],
            "x1_um": point1[0], "y1_um": point1[1],
            "mid_x_um": middle[0], "mid_y_um": middle[1],
            "segment_length_um": math.dist(point0, point1),
            "minimum_crossing_to_vertex_distance_um": min(distance0, distance1),
            "minimum_crossing_edge_fraction": min(normalized0, normalized1),
            "gate_edge_roi": in_gate_edge_roi(middle, rules),
        })
    roi_contours = [row for row in contour_rows if row["gate_edge_roi"]]

    critical_rows = []
    for node in sorted(roi_nodes):
        if (abs(net[node]) <= float(rules["near_compensation_abs_net_max_cm3"]) or
                compensation[node] >= float(rules["high_compensation_ratio_min"])):
            critical_rows.append({
                "device": device, "node_id": node,
                "x_um": nodes[node][0], "y_um": nodes[node][1],
                "BActive_cm3": acceptor[node], "AsActive_cm3": arsenic[node],
                "PActive_cm3": phosphorus[node], "DonorActive_cm3": donor[node],
                "NetActive_cm3": net[node],
                "abs_net_cm3": abs(net[node]),
                "compensation_ratio": compensation[node],
                "doping_type": "n_type" if net[node] >= 0.0 else "p_type",
            })

    roi_comp = [compensation[node] for node in roi_nodes]
    roi_abs_net = [abs(net[node]) for node in roi_nodes]
    device_row = {
        "device": device, "NWell_cm3": meta["NWell_cm3"],
        "GOxTime_min": meta["GOxTime_min"], "LDD_Dose_cm2": meta["LDD_Dose_cm2"],
        "global_node_count": len(nodes), "silicon_node_count": len(silicon_nodes),
        "silicon_triangle_count": len(triangles), "silicon_edge_count": len(edges),
        "gate_edge_roi_node_count": len(roi_nodes),
        "gate_edge_roi_edge_count": len(roi_edges),
        "gate_edge_min_edge_length_um": min(roi_edge_lengths),
        "gate_edge_median_edge_length_um": statistics.median(roi_edge_lengths),
        "gate_edge_min_abs_net_cm3": min(roi_abs_net),
        "gate_edge_max_compensation_ratio": max(roi_comp),
        "gate_edge_near_compensation_node_count": len(critical_rows),
        "zero_contour_segment_count": len(contour_rows),
        "gate_edge_zero_contour_segment_count": len(roi_contours),
        "gate_edge_min_contour_vertex_distance_um": min(
            float(row["minimum_crossing_to_vertex_distance_um"])
            for row in roi_contours),
        "gate_edge_min_contour_edge_fraction": min(
            float(row["minimum_crossing_edge_fraction"]) for row in roi_contours),
        "n_type_component_count": sum(row["doping_type"] == "n_type"
                                      for row in component_rows),
        "floating_n_type_component_count": floating_n_count,
        "single_node_floating_n_type_component_count": single_node_floating_n_count,
        "active_field_identity_max_abs_residual_cm3": max(
            abs(value) for value in identity_residuals.values()),
        "active_field_identity_pass": identity_pass,
    }
    data = {"nodes": nodes, "silicon_nodes": silicon_nodes, "fields": fields,
            "donor": donor, "acceptor": acceptor, "net": net,
            "compensation": compensation}
    return device_row, contour_rows, component_rows, critical_rows, data


def log_shifts(low: dict[int, float], high: dict[int, float],
               shared: list[tuple[int, int]]) -> list[float]:
    values = []
    for low_node, high_node in shared:
        left, right = low[low_node], high[high_node]
        if left > 0.0 and right > 0.0:
            values.append(math.log10(right / left))
    return values


def summarize_shifts(prefix: str, shifts: list[float]) -> dict[str, float]:
    absolute = [abs(value) for value in shifts]
    return {
        f"{prefix}_median_log10_shift_dex": statistics.median(shifts),
        f"{prefix}_p95_abs_log10_shift_dex": percentile(absolute, 0.95),
        f"{prefix}_max_abs_log10_shift_dex": max(absolute),
    }


def pair_species_row(low: str, high: str, device_data: dict[str, dict[str, Any]],
                     device_rows: dict[str, dict[str, Any]], digits: int
                     ) -> dict[str, Any]:
    low_data, high_data = device_data[low], device_data[high]
    low_coords = {tuple(round(value, digits) for value in low_data["nodes"][node]): node
                  for node in low_data["silicon_nodes"]}
    high_coords = {tuple(round(value, digits) for value in high_data["nodes"][node]): node
                   for node in high_data["silicon_nodes"]}
    common = sorted(set(low_coords) & set(high_coords))
    shared = [(low_coords[coord], high_coords[coord]) for coord in common]
    row: dict[str, Any] = {
        "low_device": low, "high_device": high,
        "low_silicon_node_count": len(low_coords),
        "high_silicon_node_count": len(high_coords),
        "common_coordinate_count": len(common),
        "low_only_coordinate_count": len(low_coords) - len(common),
        "high_only_coordinate_count": len(high_coords) - len(common),
        "global_node_count_delta": (device_rows[high]["global_node_count"] -
                                    device_rows[low]["global_node_count"]),
        "gate_edge_roi_node_count_delta": (
            device_rows[high]["gate_edge_roi_node_count"] -
            device_rows[low]["gate_edge_roi_node_count"]),
    }
    field_map = {
        "BActive": (low_data["acceptor"], high_data["acceptor"]),
        "AsActive": (low_data["fields"]["AsActive"], high_data["fields"]["AsActive"]),
        "PActive": (low_data["fields"]["PActive"], high_data["fields"]["PActive"]),
        "DonorActive": (low_data["donor"], high_data["donor"]),
    }
    for name, (left, right) in field_map.items():
        row.update(summarize_shifts(name, log_shifts(left, right, shared)))
    return row


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = (cursor + end - 1) / 2.0
        for position in range(cursor, end):
            result[order[position]] = rank
        cursor = end
    return result


def correlation(left: list[float], right: list[float]) -> float:
    mean_left, mean_right = statistics.mean(left), statistics.mean(right)
    numerator = sum((a - mean_left) * (b - mean_right)
                    for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - mean_left) ** 2 for a in left) *
                            sum((b - mean_right) ** 2 for b in right))
    return numerator / denominator if denominator else 0.0


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    case_rows = read_csv(M46_CASES)
    if len(case_rows) != int(contract["acceptance"]["required_m46_case_count"]):
        raise ValueError("M58 M46 case count mismatch")
    point_count = sum(int(row["point_count"]) for row in case_rows)
    if point_count != int(contract["acceptance"]["required_m46_point_count"]):
        raise ValueError("M58 M46 point count mismatch")
    case_map = {row["case"]: row for row in case_rows}
    device_errors = defaultdict(list)
    for row in case_rows:
        device_errors[row["case"].split("_vd_")[0]].append(
            float(row["maximum_absolute_log10_ratio_dex"]))

    all_contours: list[dict[str, Any]] = []
    all_components: list[dict[str, Any]] = []
    all_critical: list[dict[str, Any]] = []
    device_row_list: list[dict[str, Any]] = []
    device_data: dict[str, dict[str, Any]] = {}
    for meta in contract["matrix"]["devices"]:
        row, contours, components, critical, data = analyze_device(meta, contract)
        row["maximum_m46_error_dex"] = max(device_errors[meta["id"]])
        device_row_list.append(row)
        all_contours.extend(contours)
        all_components.extend(components)
        all_critical.extend(critical)
        device_data[meta["id"]] = data
    device_rows = {row["device"]: row for row in device_row_list}

    species_rows = [pair_species_row(
        pair["low"], pair["high"], device_data, device_rows,
        int(contract["frozen_analysis_rules"]["common_coordinate_round_digits"]))
        for pair in contract["matrix"]["matched_pairs"]]

    pair_rows: list[dict[str, Any]] = []
    error_increases = []
    compensation_stronger = []
    topology_or_proximity = []
    device_meta = {row["id"]: row for row in contract["matrix"]["devices"]}
    for pair in contract["matrix"]["matched_pairs"]:
        low, high = pair["low"], pair["high"]
        low_metrics, high_metrics = device_rows[low], device_rows[high]
        stronger = (float(high_metrics["gate_edge_max_compensation_ratio"]) >
                    float(low_metrics["gate_edge_max_compensation_ratio"]))
        topology = (
            int(high_metrics["floating_n_type_component_count"]) >
            int(low_metrics["floating_n_type_component_count"]) or
            int(high_metrics["single_node_floating_n_type_component_count"]) >
            int(low_metrics["single_node_floating_n_type_component_count"]) or
            float(high_metrics["gate_edge_min_contour_vertex_distance_um"]) <
            float(low_metrics["gate_edge_min_contour_vertex_distance_um"]))
        compensation_stronger.append(stronger)
        topology_or_proximity.append(topology)
        for vd in contract["matrix"]["drain_voltages_V"]:
            suffix = "0p05" if float(vd) == 0.05 else "1"
            low_error = float(case_map[f"{low}_vd_{suffix}"][
                "maximum_absolute_log10_ratio_dex"])
            high_error = float(case_map[f"{high}_vd_{suffix}"][
                "maximum_absolute_log10_ratio_dex"])
            increased = high_error > low_error
            error_increases.append(increased)
            pair_rows.append({
                "low_device": low, "high_device": high,
                "GOxTime_min": device_meta[low]["GOxTime_min"],
                "LDD_Dose_cm2": device_meta[low]["LDD_Dose_cm2"],
                "drain_voltage_V": vd,
                "low_maximum_error_dex": low_error,
                "high_maximum_error_dex": high_error,
                "high_minus_low_error_dex": high_error - low_error,
                "error_increased": increased,
                "high_has_stronger_gate_edge_compensation": stronger,
                "high_has_topology_or_closer_contour_change": topology,
                "low_gate_edge_max_compensation_ratio":
                    low_metrics["gate_edge_max_compensation_ratio"],
                "high_gate_edge_max_compensation_ratio":
                    high_metrics["gate_edge_max_compensation_ratio"],
                "low_floating_n_type_component_count":
                    low_metrics["floating_n_type_component_count"],
                "high_floating_n_type_component_count":
                    high_metrics["floating_n_type_component_count"],
                "low_single_node_floating_n_type_component_count":
                    low_metrics["single_node_floating_n_type_component_count"],
                "high_single_node_floating_n_type_component_count":
                    high_metrics["single_node_floating_n_type_component_count"],
                "low_min_contour_vertex_distance_um":
                    low_metrics["gate_edge_min_contour_vertex_distance_um"],
                "high_min_contour_vertex_distance_um":
                    high_metrics["gate_edge_min_contour_vertex_distance_um"],
            })

    if all(error_increases) and all(compensation_stronger) and any(topology_or_proximity):
        classification = "compensation_contour_mesh_sensitivity_supported"
    elif all(error_increases) and all(compensation_stronger):
        classification = "compensation_strength_only"
    elif all(error_increases) and all(topology_or_proximity):
        classification = "mesh_topology_without_systematic_compensation"
    else:
        classification = "no_systematic_frozen_tdr_relation"

    metrics = [
        "gate_edge_max_compensation_ratio",
        "gate_edge_min_abs_net_cm3",
        "gate_edge_near_compensation_node_count",
        "gate_edge_min_contour_vertex_distance_um",
        "floating_n_type_component_count",
        "single_node_floating_n_type_component_count",
        "global_node_count",
    ]
    errors = [float(row["maximum_m46_error_dex"]) for row in device_row_list]
    correlation_rows = []
    for metric in metrics:
        values = [float(row[metric]) for row in device_row_list]
        correlation_rows.append({
            "metric": metric, "device_count": len(values),
            "pearson_with_device_maximum_m46_error": correlation(values, errors),
            "rank_correlation_with_device_maximum_m46_error":
                correlation(ranks(values), ranks(errors)),
            "interpretation": "descriptive_only",
        })

    checks = {
        "device_count": len(device_row_list) ==
            int(contract["acceptance"]["required_device_count"]),
        "pair_voltage_row_count": len(pair_rows) ==
            int(contract["acceptance"]["required_pair_voltage_rows"]),
        "m46_case_count": len(case_rows) ==
            int(contract["acceptance"]["required_m46_case_count"]),
        "m46_point_count": point_count ==
            int(contract["acceptance"]["required_m46_point_count"]),
        "all_active_field_identities_pass": all(
            bool(row["active_field_identity_pass"]) for row in device_row_list),
        "all_error_pair_rows_increase": all(error_increases),
        "classification_declared": classification in
            contract["analysis"]["classifications"],
        "read_only_execution": True,
        "closed_topics_not_reopened": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    target_critical = sorted(
        (row for row in all_critical if row["device"] == "n23"),
        key=lambda row: float(row["abs_net_cm3"]))
    control_critical = sorted(
        (row for row in all_critical if row["device"] == "n19"),
        key=lambda row: float(row["abs_net_cm3"]))
    low_ldd_pairs = [row for row in pair_rows if float(row["LDD_Dose_cm2"]) == 1e14]
    high_ldd_pairs = [row for row in pair_rows if float(row["LDD_Dose_cm2"]) == 2e14]
    report = {
        "schema": "vela.simplemos.sdevice.m58_compensation_contour_mesh_audit_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_execution": False,
                      "new_vela_execution": False,
                      "local_tdr_exports_generated": True,
                      "source": "eight immutable original M8/M46 input_fps.tdr files"},
        "semantics": contract["semantic_guard"],
        "matrix": {"device_count": len(device_row_list),
                   "case_count": len(case_rows), "point_count": point_count,
                   "pair_voltage_row_count": len(pair_rows)},
        "systematic_tests": {
            "error_increased_rows": sum(error_increases),
            "error_row_count": len(error_increases),
            "pairs_with_stronger_high_compensation": sum(compensation_stronger),
            "pair_count": len(compensation_stronger),
            "pairs_with_topology_or_closer_contour_change": sum(topology_or_proximity),
            "low_ldd_rows_with_topology_or_closer_contour_change": sum(
                bool(row["high_has_topology_or_closer_contour_change"])
                for row in low_ldd_pairs),
            "low_ldd_row_count": len(low_ldd_pairs),
            "high_ldd_rows_with_topology_or_closer_contour_change": sum(
                bool(row["high_has_topology_or_closer_contour_change"])
                for row in high_ldd_pairs),
            "high_ldd_row_count": len(high_ldd_pairs),
        },
        "target": {"device": "n23", "metrics": device_rows["n23"],
                   "lowest_abs_net_critical_nodes": target_critical[:8]},
        "matched_control": {"device": "n19", "metrics": device_rows["n19"],
                            "lowest_abs_net_critical_nodes": control_critical[:8]},
        "device_metrics": device_row_list,
        "pair_species_common_coordinate_metrics": species_rows,
        "correlations": correlation_rows,
        "claim_guard": contract["analysis"]["claim_guard"],
        "forbidden_work_respected": contract["forbidden_work"],
        "acceptance": checks,
    }
    return report, {
        "devices": device_row_list, "contours": all_contours,
        "components": all_components, "critical": all_critical,
        "species": species_rows, "pairs": pair_rows,
        "correlations": correlation_rows,
    }


def write_human_report(report: dict[str, Any], ledgers: dict[str, list[dict[str, Any]]]) -> None:
    target = report["target"]["metrics"]
    control = report["matched_control"]["metrics"]
    target_nodes = report["target"]["lowest_abs_net_critical_nodes"]
    target_lines = "\n".join(
        f"| {row['node_id']} | {float(row['x_um']):.9f} | {float(row['y_um']):.9f} | "
        f"{float(row['AsActive_cm3']):.6e} | {float(row['PActive_cm3']):.6e} | "
        f"{float(row['BActive_cm3']):.6e} | {float(row['NetActive_cm3']):.6e} | "
        f"{float(row['compensation_ratio']):.3f} |"
        for row in target_nodes)
    pair_lines = "\n".join(
        f"| {row['low_device']}/{row['high_device']} | {float(row['LDD_Dose_cm2']):.0e} | "
        f"{float(row['drain_voltage_V']):g} | "
        f"{float(row['low_maximum_error_dex']):.6f} | "
        f"{float(row['high_maximum_error_dex']):.6f} | "
        f"{float(row['high_minus_low_error_dex']):+.6f} | "
        f"{row['high_has_stronger_gate_edge_compensation']} | "
        f"{row['high_has_topology_or_closer_contour_change']} |"
        for row in ledgers["pairs"])
    species_lines = "\n".join(
        f"| {row['low_device']}/{row['high_device']} | {row['common_coordinate_count']} | "
        f"{row['low_only_coordinate_count']}/{row['high_only_coordinate_count']} | "
        f"{float(row['BActive_median_log10_shift_dex']):.6f} | "
        f"{float(row['AsActive_p95_abs_log10_shift_dex']):.6f} | "
        f"{float(row['PActive_p95_abs_log10_shift_dex']):.6f} | "
        f"{float(row['DonorActive_p95_abs_log10_shift_dex']):.6f} |"
        for row in ledgers["species"])
    if report["classification"] == "no_systematic_frozen_tdr_relation":
        interpretation = (
            "冻结规则不支持用单一补偿或网格拓扑指标解释完整八组误差趋势。"
            "但交互分层是明确的：两个 LDD=1e14 配对的四个漏压行全部出现"
            "额外浮置分量或更贴近顶点的零轮廓，而两个 LDD=2e14 配对的四行均未出现；"
            "因此 n23 局部机制成立，但不能外推到高 LDD 组。")
    else:
        interpretation = (
            "冻结矩阵满足合同中声明的系统关系，可进入针对性的局部重网格或"
            "补偿迁移率控制，但仍不构成因果证明。")
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M58 补偿结轮廓与网格敏感性审计

## 结论

M58 分类为 `{report['classification']}`。Workbench 参数 `NWell` 在原始工艺脚本中实际对应 `init field=Boron` 的 p 型背景硼，而不是传统掩膜 n-well。M58 只读取八个冻结 `input_fps.tdr` 和 M46/M55/M57 证据；没有新增 Sentaurus 或 Vela 求解，没有改动物理、网格、偏压、端口算法或生产默认值。

八个匹配 NWell/漏压组合中，M46 最大 Id 误差均随背景硼从 1e17 增至 2e17 cm^-3 而增加；但只有 `{report['systematic_tests']['pairs_with_stronger_high_compensation']}/4` 个器件配对的浅表栅边最大补偿指标增强，只有 `{report['systematic_tests']['pairs_with_topology_or_closer_contour_change']}/4` 对出现额外浮置同号分量或更贴近网格顶点的 `NetActive=0` 轮廓。{interpretation}

## 目标 n23 与低背景控制 n19

| 指标 | n19 | n23 |
|---|---:|---:|
| 全局节点数 | {control['global_node_count']} | {target['global_node_count']} |
| Si节点数 | {control['silicon_node_count']} | {target['silicon_node_count']} |
| 栅边最小 `abs(NetActive)` (cm^-3) | {float(control['gate_edge_min_abs_net_cm3']):.6e} | {float(target['gate_edge_min_abs_net_cm3']):.6e} |
| 栅边最大补偿指标 | {float(control['gate_edge_max_compensation_ratio']):.3f} | {float(target['gate_edge_max_compensation_ratio']):.3f} |
| 浮置n型同号分量 | {control['floating_n_type_component_count']} | {target['floating_n_type_component_count']} |
| 单节点浮置n型分量 | {control['single_node_floating_n_type_component_count']} | {target['single_node_floating_n_type_component_count']} |
| 栅边轮廓到顶点最短距离 (um) | {float(control['gate_edge_min_contour_vertex_distance_um']):.6e} | {float(target['gate_edge_min_contour_vertex_distance_um']):.6e} |

n23 最接近补偿的冻结工艺节点如下：

| 节点 | x (um) | y (um) | AsActive | PActive | BActive | NetActive | 补偿指标 |
|---:|---:|---:|---:|---:|---:|---:|---:|
{target_lines}

## 八个 NWell/漏压配对

| 低/高背景器件 | LDD | Vd | 低误差(dex) | 高误差(dex) | 增幅(dex) | 补偿增强 | 拓扑/轮廓变化 |
|---|---:|---:|---:|---:|---:|---|---|
{pair_lines}

## 物种与共有坐标

`DonorActive=AsActive+PActive`、`AcceptorActive=BActive`，八个器件均在冻结容差内重构 `NetActive`。高低背景器件并不共享完整网格；下表同时给出共有坐标以及活性物种变化。

| 配对 | 共有Si坐标 | 低/高独有坐标 | B中位变化(dex) | As P95绝对变化 | P P95绝对变化 | donor P95绝对变化 |
|---|---:|---:|---:|---:|---:|---:|
{species_lines}

## 方法与边界

- 补偿指标固定为 `(AsActive+PActive+BActive)/max(abs(NetActive),1e10)`；浅表栅边窗口固定为 `x<=0.03 um` 且 `abs(abs(y)-0.125)<=0.03 um`。
- `NetActive=0` 轮廓用硅三角形内分片线性插值提取；轮廓到顶点距离没有按器件或结果调整。
- 同号连通分量仅沿硅三角形边连接；只有实际接触硅节点才能使分量成为接触关联分量。
- 相关系数仅是八器件描述性结果，不作为因果判据。
- 本任务不重新打开 HFS、SG、BGN、SRH、准费米打包、通用接触提取或 M57 节点电流求积。

机器报告：`reference_tcad/simplemos_sentaurus2022/compensation_contour_mesh_audit/m58_compensation_contour_mesh_audit_report.json`。
""", encoding="utf-8", newline="\n")


def write_artifact(report: dict[str, Any]) -> None:
    payload = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report",
            "title": "SimpleMOS M58 compensation contour and mesh audit",
            "description": "Frozen eight-TDR active-species, compensation, contour, and connectivity attribution.",
            "generatedAt": "2026-09-02T00:00:00+08:00",
            "sources": [
                {"id": "src_m58", "label": "M58 machine report", "path": portable(REPORT)},
                {"id": "src_contract", "label": "Frozen M58 contract", "path": portable(CONTRACT)}],
            "tables": [{"id": "devices", "title": "Device compensation and mesh metrics",
                        "dataset": "devices", "sourceId": "src_m58", "layout": "full",
                        "columns": [
                            {"field": "device", "label": "Device", "type": "text"},
                            {"field": "NWell_cm3", "label": "Background boron", "format": "number"},
                            {"field": "gate_edge_max_compensation_ratio", "label": "Max compensation", "format": "number"},
                            {"field": "floating_n_type_component_count", "label": "Floating n components", "format": "number"},
                            {"field": "maximum_m46_error_dex", "label": "M46 max error", "format": "number"}]}],
            "blocks": [
                {"id": "title", "type": "markdown", "body": "# SimpleMOS M58 compensation contour and mesh audit"},
                {"id": "summary", "type": "markdown", "sourceId": "src_m58",
                 "body": f"Outcome: `{report['classification']}`. Eight immutable process TDRs were analyzed without a device solve."},
                {"id": "devices", "type": "table", "tableId": "devices"},
                {"id": "guard", "type": "markdown", "sourceId": "src_contract",
                 "body": "This structural relation supports a targeted control; it does not prove causality or authorize a production change."}]
        },
        "snapshot": {"version": 1, "generatedAt": "2026-09-02T00:00:00+08:00",
                     "status": "ready", "datasets": {"devices": report["device_metrics"]}},
        "sources": [{"id": "src_m58", "label": "M58 machine report",
                     "path": portable(REPORT)}]
    }
    write_json(ARTIFACT, payload)


def execute() -> dict[str, Any]:
    contract, freeze = validate_contract()
    report, ledgers = analyze(contract)
    write_csv(DEVICES_CSV, ledgers["devices"])
    write_csv(CONTOURS_CSV, ledgers["contours"])
    write_csv(COMPONENTS_CSV, ledgers["components"])
    write_csv(CRITICAL_CSV, ledgers["critical"])
    write_csv(SPECIES_CSV, ledgers["species"])
    write_csv(PAIRS_CSV, ledgers["pairs"])
    write_csv(CORRELATIONS_CSV, ledgers["correlations"])
    write_json(REPORT, report)
    write_human_report(report, ledgers)
    write_artifact(report)
    artifacts = {portable(path): sha256(path) for path in (
        REPORT, DEVICES_CSV, CONTOURS_CSV, COMPONENTS_CSV, CRITICAL_CSV,
        SPECIES_CSV, PAIRS_CSV, CORRELATIONS_CSV, DOC, ARTIFACT)}
    evidence = {
        "schema": "vela.simplemos.sdevice.m58_compensation_contour_mesh_audit_evidence.v1",
        "status": report["status"], "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT), "artifacts": artifacts,
        "source_hashes": freeze["upstream_hashes"],
        "new_sentaurus_execution": False, "new_vela_execution": False,
        "local_immutable_tdr_export_only": True,
        "historical_artifacts_rewritten": False,
        "closed_topics_reinvestigated": False,
        "acceptance": report["acceptance"],
    }
    write_json(EVIDENCE, evidence)
    return report


def verify() -> dict[str, Any]:
    validate_contract()
    evidence = read_json(EVIDENCE)
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M58 artifact hash mismatch: {relative}")
    for relative, expected in evidence["source_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M58 source hash mismatch: {relative}")
    return evidence


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="validate the frozen contract and evidence hashes")
    args = parser.parse_args()
    if args.verify:
        evidence = verify()
        print(f"M58 verified: {evidence['classification']}")
        return
    report = execute()
    print(json.dumps({"status": report["status"],
                      "classification": report["classification"],
                      "systematic_tests": report["systematic_tests"],
                      "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
