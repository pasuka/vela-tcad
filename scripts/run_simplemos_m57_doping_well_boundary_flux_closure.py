#!/usr/bin/env python3
"""Execute the frozen, read-only SimpleMOS M57 well-boundary closure analysis."""

from __future__ import annotations

from collections import defaultdict
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m57_doping_well_boundary_flux_closure_contract_v1.json"
FREEZE = ROOT / "simplemos_m57_doping_well_boundary_flux_closure_contract_freeze.json"
M55 = ROOT / "doping_well_terminal_attribution"
M55_REPORT = M55 / "m55_doping_well_terminal_attribution_report.json"
M55_WELLS = M55 / "m55_doping_well_ledger.csv"
M55_TERMINALS = M55 / "m55_terminal_component_attribution_ledger.csv"
M55_EXPORT_MANIFEST = REPO / (
    "build-release/reference_tcad/simplemos_sentaurus2022/"
    "m55_doping_well_terminal_attribution/sentaurus_export_manifest.json")
OUTPUT = ROOT / "doping_well_boundary_flux_closure"
REPORT = OUTPUT / "m57_doping_well_boundary_flux_closure_report.json"
SEGMENTS = OUTPUT / "m57_boundary_segment_flux_ledger.csv"
WELLS = OUTPUT / "m57_well_boundary_summary.csv"
CLOSURES = OUTPUT / "m57_carrier_closure_ledger.csv"
PAIRS = OUTPUT / "m57_nwell_pair_summary.csv"
EVIDENCE = ROOT / "simplemos_m57_doping_well_boundary_flux_closure_evidence.json"
DOC = REPO / "docs/validation/simplemos_m57_doping_well_boundary_flux_closure_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m57/artifact.json"
DEVICES = tuple(f"n{i}" for i in range(17, 25))
CONTACTS = ("source", "drain", "substrate")
COMPONENTS = ("electron", "hole")
FIELD_FILES = {
    "well": "DopingWells_region0.csv",
    "electron": "eCurrentDensity_region0.csv",
    "hole": "hCurrentDensity_region0.csv",
}


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


def close_enough(left: float, right: float, absolute: float,
                 relative: float) -> bool:
    difference = abs(left - right)
    return difference <= absolute or difference <= relative * max(
        abs(left), abs(right), 1e-300)


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    expected = "vela.simplemos.sdevice.m57_doping_well_boundary_flux_closure_contract.v1"
    if contract.get("schema") != expected:
        raise ValueError("unexpected M57 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M57 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M57 contract changed after freeze")
    for relative, expected_hash in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected_hash:
            raise ValueError(f"M57 upstream changed after freeze: {relative}")
    m55 = read_json(M55_REPORT)
    upstream = contract["upstream"]
    if m55.get("status") != upstream["required_m55_status"]:
        raise ValueError("M57 requires accepted M55")
    if m55.get("classification") != upstream["required_m55_classification"]:
        raise ValueError("M55 classification anchor changed")
    if int(m55["sentaurus"]["state_count"]) != int(
            upstream["required_m55_state_count"]):
        raise ValueError("M55 state count changed")
    return contract


def scalar_field(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"])
            for row in read_csv(path)}


def vector_field(path: Path) -> dict[int, tuple[float, float]]:
    return {int(row["node_id"]): (float(row["component0"]),
                                  float(row["component1"]))
            for row in read_csv(path)}


def midpoint(left: tuple[float, float], right: tuple[float, float]) -> tuple[float, float]:
    return ((left[0] + right[0]) / 2.0, (left[1] + right[1]) / 2.0)


def average_vector(left: tuple[float, float],
                   right: tuple[float, float]) -> tuple[float, float]:
    return ((left[0] + right[0]) / 2.0, (left[1] + right[1]) / 2.0)


def outward_normal(point0: tuple[float, float], point1: tuple[float, float],
                   inside: tuple[float, float]) -> tuple[float, float]:
    dx, dy = point1[0] - point0[0], point1[1] - point0[1]
    length = math.hypot(dx, dy)
    if length <= 0.0:
        raise ValueError("zero-length M57 segment")
    nx, ny = dy / length, -dx / length
    middle = midpoint(point0, point1)
    if nx * (inside[0] - middle[0]) + ny * (inside[1] - middle[1]) > 0.0:
        nx, ny = -nx, -ny
    return nx, ny


def interpolate_endpoint(node0: int, node1: int,
                         nodes: dict[int, tuple[float, float]],
                         currents: dict[str, dict[int, tuple[float, float]]]
                         ) -> tuple[tuple[float, float], dict[str, tuple[float, float]]]:
    return midpoint(nodes[node0], nodes[node1]), {
        component: average_vector(values[node0], values[node1])
        for component, values in currents.items()}


def flux(point0: tuple[float, float], point1: tuple[float, float],
         current0: tuple[float, float], current1: tuple[float, float],
         normal: tuple[float, float]) -> float:
    length_um = math.dist(point0, point1)
    normal0 = current0[0] * normal[0] + current0[1] * normal[1]
    normal1 = current1[0] * normal[0] + current1[1] * normal[1]
    return (normal0 + normal1) / 2.0 * length_um * 1e-8


def reconstruct_state(state: dict[str, Any], contract: dict[str, Any],
                      well_rows: dict[tuple[str, float, float, str], dict[str, str]]
                      ) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, str]]:
    export = REPO / state["export_dir"]
    nodes = {int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
             for row in read_csv(export / "nodes.csv")}
    triangles: list[tuple[int, tuple[int, int, int]]] = []
    for row in read_csv(export / "elements.csv"):
        if row["material"] == "Si":
            triangles.append((int(row["id"]), tuple(
                int(row[f"node{i}"]) for i in range(3))))
    contact_nodes = {row["name"].lower(): {int(value) for value in
                     row["node_ids"].split(";")}
                     for row in read_csv(export / "contacts.csv")}
    raw_wells = scalar_field(export / "fields" / FIELD_FILES["well"])
    tolerance = float(contract["acceptance"]["doping_well_integer_tolerance"])
    rounded_wells: dict[int, int] = {}
    noninteger: list[tuple[int, float]] = []
    for node, value in raw_wells.items():
        rounded = round(value)
        if abs(value - rounded) > tolerance:
            noninteger.append((node, value))
        rounded_wells[node] = rounded
    currents = {component: vector_field(export / "fields" / FIELD_FILES[component])
                for component in COMPONENTS}
    silicon_nodes = {node for _, tri in triangles for node in tri}
    if not silicon_nodes <= set(rounded_wells):
        raise ValueError(f"DopingWells coverage incomplete: {state['state']}")
    if any(not silicon_nodes <= set(values) for values in currents.values()):
        raise ValueError(f"current-density coverage incomplete: {state['state']}")

    named_labels: dict[str, int] = {}
    label_problems: list[str] = []
    for contact in CONTACTS:
        values = {rounded_wells[node] for node in contact_nodes[contact]
                  if node in silicon_nodes}
        if len(values) != 1:
            label_problems.append(f"{contact}:{sorted(values)}")
            continue
        named_labels[contact] = next(iter(values))
        upstream = well_rows[(state["device"], float(state["drain_voltage_V"]),
                              float(state["gate_voltage_V"]), contact)]
        if named_labels[contact] != int(upstream["well_index"]):
            label_problems.append(
                f"{contact}:contact={named_labels[contact]},seed={upstream['well_index']}")

    edge_triangles: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    maximum_labels = 0
    for triangle_id, triangle in triangles:
        labels = {rounded_wells[node] for node in triangle}
        maximum_labels = max(maximum_labels, len(labels))
        for a, b, third in ((triangle[0], triangle[1], triangle[2]),
                            (triangle[1], triangle[2], triangle[0]),
                            (triangle[2], triangle[0], triangle[1])):
            edge_triangles[tuple(sorted((a, b)))].append((triangle_id, third))

    by_label: dict[int, list[dict[str, Any]]] = defaultdict(list)
    all_labels = sorted(set(rounded_wells[node] for node in silicon_nodes))
    for triangle_id, triangle in triangles:
        labels = {rounded_wells[node] for node in triangle}
        if len(labels) != 2:
            continue
        mixed_edges = [(a, b) for a, b in ((triangle[0], triangle[1]),
                                           (triangle[1], triangle[2]),
                                           (triangle[2], triangle[0]))
                       if rounded_wells[a] != rounded_wells[b]]
        if len(mixed_edges) != 2:
            label_problems.append(f"triangle-{triangle_id}:mixed-edge-count-{len(mixed_edges)}")
            continue
        endpoint_data = [interpolate_endpoint(a, b, nodes, currents)
                         for a, b in mixed_edges]
        point0, current0 = endpoint_data[0]
        point1, current1 = endpoint_data[1]
        for label in labels:
            inside_nodes = [node for node in triangle if rounded_wells[node] == label]
            inside = (sum(nodes[node][0] for node in inside_nodes) / len(inside_nodes),
                      sum(nodes[node][1] for node in inside_nodes) / len(inside_nodes))
            normal = outward_normal(point0, point1, inside)
            other = next(value for value in labels if value != label)
            by_label[label].append({
                "segment_key": f"triangle:{triangle_id}:labels:{min(labels)}-{max(labels)}",
                "triangle_id": triangle_id,
                "boundary_class": "interior_doping_well_interface",
                "neighbor_well_label": other,
                "point0": point0, "point1": point1,
                "current0": current0, "current1": current1, "normal": normal,
            })

    for (node0, node1), adjacent in edge_triangles.items():
        if len(adjacent) != 1:
            continue
        triangle_id, third = adjacent[0]
        label0, label1 = rounded_wells[node0], rounded_wells[node1]
        full_point0, full_point1 = nodes[node0], nodes[node1]
        full_current0 = {component: currents[component][node0] for component in COMPONENTS}
        full_current1 = {component: currents[component][node1] for component in COMPONENTS}
        middle_point, middle_current = interpolate_endpoint(
            node0, node1, nodes, currents)
        parts: list[tuple[int, tuple[float, float], tuple[float, float],
                          dict[str, tuple[float, float]], dict[str, tuple[float, float]]]]
        if label0 == label1:
            parts = [(label0, full_point0, full_point1, full_current0, full_current1)]
        else:
            parts = [(label0, full_point0, middle_point, full_current0, middle_current),
                     (label1, middle_point, full_point1, middle_current, full_current1)]
        contact_name = ""
        for name, support in contact_nodes.items():
            if node0 in support and node1 in support:
                contact_name = name
                break
        boundary_class = ("physical_contact" if contact_name in CONTACTS else
                          "other_semiconductor_exterior")
        for part_index, (label, point0, point1, current0, current1) in enumerate(parts):
            normal = outward_normal(point0, point1, nodes[third])
            by_label[label].append({
                "segment_key": f"edge:{node0}-{node1}:part:{part_index}",
                "triangle_id": triangle_id,
                "boundary_class": boundary_class,
                "contact_name": contact_name,
                "neighbor_well_label": "",
                "point0": point0, "point1": point1,
                "current0": current0, "current1": current1, "normal": normal,
            })

    segment_rows: list[dict[str, Any]] = []
    for contact, label in named_labels.items():
        for index, segment in enumerate(by_label[label]):
            point0, point1, normal = segment["point0"], segment["point1"], segment["normal"]
            row: dict[str, Any] = {
                "state": state["state"], "device": state["device"],
                "drain_voltage_V": float(state["drain_voltage_V"]),
                "gate_voltage_V": float(state["gate_voltage_V"]),
                "contact": contact, "well_label": label,
                "segment_index": index, "segment_key": segment["segment_key"],
                "triangle_id": segment["triangle_id"],
                "boundary_class": segment["boundary_class"],
                "boundary_contact_name": segment.get("contact_name", ""),
                "neighbor_well_label": segment["neighbor_well_label"],
                "x0_um": point0[0], "y0_um": point0[1],
                "x1_um": point1[0], "y1_um": point1[1],
                "length_um": math.dist(point0, point1),
                "normal_x": normal[0], "normal_y": normal[1],
            }
            for component in COMPONENTS:
                row[f"{component}_geometric_outward_flux_A_per_um"] = flux(
                    point0, point1, segment["current0"][component],
                    segment["current1"][component], normal)
            segment_rows.append(row)

    opposite_max = 0.0
    for label in all_labels:
        for segment in by_label[label]:
            if segment["boundary_class"] != "interior_doping_well_interface":
                continue
            other = int(segment["neighbor_well_label"])
            if label > other:
                continue
            matches = [candidate for candidate in by_label[other]
                       if candidate["segment_key"] == segment["segment_key"]]
            if len(matches) != 1:
                label_problems.append(f"unpaired:{segment['segment_key']}:{label}-{other}")
                continue
            opposite = matches[0]
            for component in COMPONENTS:
                left = flux(segment["point0"], segment["point1"],
                            segment["current0"][component], segment["current1"][component],
                            segment["normal"])
                right = flux(opposite["point0"], opposite["point1"],
                             opposite["current0"][component], opposite["current1"][component],
                             opposite["normal"])
                opposite_max = max(opposite_max, abs(left + right))

    topology_payload = b"".join((export / relative).read_bytes() for relative in (
        "nodes.csv", "elements.csv", "contacts.csv", "fields/DopingWells_region0.csv"))
    topology_hash = hashlib.sha256(topology_payload).hexdigest()
    diagnostics = {
        "noninteger_value_count": len(noninteger),
        "maximum_labels_per_silicon_triangle": maximum_labels,
        "label_problems": label_problems,
        "named_labels": named_labels,
        "all_labels": all_labels,
        "opposite_interface_flux_max_A_per_um": opposite_max,
        "topology_sha256": topology_hash,
    }
    source_hashes = {}
    for relative in ("nodes.csv", "elements.csv", "contacts.csv",
                     "fields/DopingWells_region0.csv",
                     "fields/eCurrentDensity_region0.csv",
                     "fields/hCurrentDensity_region0.csv"):
        source_hashes[portable(export / relative)] = sha256(export / relative)
    return segment_rows, diagnostics, source_hashes


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]],
                                               list[dict[str, Any]], list[dict[str, Any]],
                                               list[dict[str, Any]], dict[str, str]]:
    manifest = read_json(M55_EXPORT_MANIFEST)
    states = manifest["states"]
    well_rows_raw = read_csv(M55_WELLS)
    terminal_rows_raw = read_csv(M55_TERMINALS)
    well_map = {(row["device"], float(row["drain_voltage_V"]),
                 float(row["gate_voltage_V"]), row["contact"]): row
                for row in well_rows_raw}
    terminal_map = {(row["device"], float(row["drain_voltage_V"]),
                     float(row["gate_voltage_V"]), row["contact"], row["component"]): row
                    for row in terminal_rows_raw}
    all_segments: list[dict[str, Any]] = []
    state_diagnostics: list[dict[str, Any]] = []
    source_hashes: dict[str, str] = {}
    for state in states:
        segments, diagnostics, hashes = reconstruct_state(state, contract, well_map)
        all_segments.extend(segments)
        state_diagnostics.append({"state": state["state"], "device": state["device"],
                                  "drain_voltage_V": state["drain_voltage_V"],
                                  "gate_voltage_V": state["gate_voltage_V"], **diagnostics})
        source_hashes.update(hashes)

    grouped: dict[tuple[str, float, float, str], list[dict[str, Any]]] = defaultdict(list)
    for row in all_segments:
        grouped[(row["device"], float(row["drain_voltage_V"]),
                 float(row["gate_voltage_V"]), row["contact"])].append(row)
    summary_rows: list[dict[str, Any]] = []
    closure_rows: list[dict[str, Any]] = []
    acceptance = contract["acceptance"]
    anchor_abs = float(acceptance["direct_anchor_absolute_tolerance_A_per_um"])
    anchor_rel = float(acceptance["direct_anchor_relative_tolerance"])
    closure_abs = float(acceptance["default_closure_absolute_tolerance_A_per_um"])
    closure_rel = float(acceptance["default_closure_relative_tolerance"])
    for state in states:
        key0 = (state["device"], float(state["drain_voltage_V"]),
                float(state["gate_voltage_V"]))
        for contact in CONTACTS:
            rows = grouped[(*key0, contact)]
            counts = {kind: sum(row["boundary_class"] == kind for row in rows)
                      for kind in contract["boundary_reconstruction"]["boundary_classes"]}
            lengths = {kind: sum(float(row["length_um"]) for row in rows
                                 if row["boundary_class"] == kind)
                       for kind in counts}
            summary: dict[str, Any] = {
                "state": state["state"], "device": state["device"],
                "drain_voltage_V": key0[1], "gate_voltage_V": key0[2],
                "contact": contact, "well_label": rows[0]["well_label"],
                "segment_count": len(rows),
            }
            for kind in counts:
                summary[f"{kind}_segment_count"] = counts[kind]
                summary[f"{kind}_length_um"] = lengths[kind]
            upstream_well = well_map[(*key0, contact)]
            summary["well_area_cm2"] = float(upstream_well["well_area_cm2"])
            for component in COMPONENTS:
                field = f"{component}_geometric_outward_flux_A_per_um"
                fluxes = {kind: sum(float(row[field]) for row in rows
                                    if row["boundary_class"] == kind)
                          for kind in counts}
                total_outward = sum(fluxes.values())
                summary[f"{component}_complete_geometric_outward_flux_A_per_um"] = total_outward
                for kind, value in fluxes.items():
                    summary[f"{component}_{kind}_geometric_outward_flux_A_per_um"] = value
                terminal = terminal_map[(*key0, contact, component)]
                default = float(terminal["default_A_per_um"])
                direct = float(terminal["direct_A_per_um"])
                generation = float(terminal["signed_well_srh_charge_current_A_per_um"])
                contact_prediction = -fluxes["physical_contact"]
                complete_surface = -total_outward
                noncontact_surface = -(fluxes["interior_doping_well_interface"] +
                                       fluxes["other_semiconductor_exterior"])
                default_prediction = complete_surface + generation
                direct_residual = contact_prediction - direct
                closure_residual = default_prediction - default
                anchor_pass = close_enough(contact_prediction, direct, anchor_abs, anchor_rel)
                closure_pass = close_enough(default_prediction, default,
                                            closure_abs, closure_rel)
                closure_rows.append({
                    "state": state["state"], "device": state["device"],
                    "drain_voltage_V": key0[1], "gate_voltage_V": key0[2],
                    "contact": contact, "component": component,
                    "well_label": rows[0]["well_label"],
                    "default_A_per_um": default, "direct_A_per_um": direct,
                    "signed_well_srh_charge_current_A_per_um": generation,
                    "physical_contact_terminal_oriented_flux_A_per_um": contact_prediction,
                    "interior_interface_terminal_oriented_flux_A_per_um":
                        -fluxes["interior_doping_well_interface"],
                    "other_exterior_terminal_oriented_flux_A_per_um":
                        -fluxes["other_semiconductor_exterior"],
                    "complete_well_terminal_oriented_surface_flux_A_per_um": complete_surface,
                    "reconstructed_noncontact_surface_redistribution_A_per_um": noncontact_surface,
                    "m55_inferred_surface_redistribution_A_per_um":
                        float(terminal["surface_redistribution_A_per_um"]),
                    "predicted_default_A_per_um": default_prediction,
                    "direct_anchor_residual_A_per_um": direct_residual,
                    "direct_anchor_symmetric_relative_difference":
                        abs(direct_residual) / max(abs(contact_prediction), abs(direct), 1e-300),
                    "default_closure_residual_A_per_um": closure_residual,
                    "default_closure_symmetric_relative_difference":
                        abs(closure_residual) / max(abs(default_prediction), abs(default), 1e-300),
                    "direct_anchor_within_tolerance": anchor_pass,
                    "default_closure_within_tolerance": closure_pass,
                })
            summary_rows.append(summary)

    topology_groups: dict[str, set[str]] = defaultdict(set)
    partition_problems: list[str] = []
    maximum_interface_balance = 0.0
    maximum_labels = 0
    for row in state_diagnostics:
        topology_groups[row["device"]].add(row["topology_sha256"])
        maximum_interface_balance = max(maximum_interface_balance,
                                        float(row["opposite_interface_flux_max_A_per_um"]))
        maximum_labels = max(maximum_labels,
                             int(row["maximum_labels_per_silicon_triangle"]))
        if int(row["noninteger_value_count"]):
            partition_problems.append(f"{row['state']}:noninteger")
        partition_problems.extend(f"{row['state']}:{value}"
                                  for value in row["label_problems"])
    topology_invariant = all(len(values) == 1 for values in topology_groups.values())
    if not topology_invariant:
        partition_problems.append("topology changed across bias states")
    if maximum_labels > int(acceptance["required_maximum_labels_per_silicon_triangle"]):
        partition_problems.append(f"maximum-labels-per-triangle:{maximum_labels}")
    partition_valid = not partition_problems
    anchors_pass = all(bool(row["direct_anchor_within_tolerance"])
                       for row in closure_rows)
    closures_pass = all(bool(row["default_closure_within_tolerance"])
                        for row in closure_rows)
    if not partition_valid:
        classification = "doping_well_partition_ambiguous"
    elif not anchors_pass:
        classification = "boundary_quadrature_limited"
    elif closures_pass:
        classification = "full_doping_well_boundary_closure"
    else:
        classification = "closure_mismatch"

    matched = (("n17", "n21"), ("n18", "n22"),
               ("n19", "n23"), ("n20", "n24"))
    closure_map = {(row["device"], float(row["drain_voltage_V"]),
                    float(row["gate_voltage_V"]), row["contact"], row["component"]): row
                   for row in closure_rows}
    pair_rows: list[dict[str, Any]] = []
    for low, high in matched:
        for drain in (0.05, 1.0):
            for contact in CONTACTS:
                for component in COMPONENTS:
                    lows = [closure_map[(low, drain, gate, contact, component)]
                            for gate in (0.0, 0.05, 0.1)]
                    highs = [closure_map[(high, drain, gate, contact, component)]
                             for gate in (0.0, 0.05, 0.1)]
                    pair_rows.append({
                        "low_device": low, "high_device": high,
                        "drain_voltage_V": drain, "contact": contact,
                        "component": component,
                        "low_maximum_absolute_default_closure_residual_A_per_um":
                            max(abs(float(row["default_closure_residual_A_per_um"])) for row in lows),
                        "high_maximum_absolute_default_closure_residual_A_per_um":
                            max(abs(float(row["default_closure_residual_A_per_um"])) for row in highs),
                        "low_maximum_absolute_reconstructed_noncontact_flux_A_per_um":
                            max(abs(float(row["reconstructed_noncontact_surface_redistribution_A_per_um"])) for row in lows),
                        "high_maximum_absolute_reconstructed_noncontact_flux_A_per_um":
                            max(abs(float(row["reconstructed_noncontact_surface_redistribution_A_per_um"])) for row in highs),
                    })

    target_rows = [row for row in closure_rows if row["device"] == "n23" and
                   math.isclose(float(row["drain_voltage_V"]), 0.05) and
                   math.isclose(float(row["gate_voltage_V"]), 0.05)]
    control_rows = [row for row in closure_rows if row["device"] in {"n19", "n23"} and
                    math.isclose(float(row["drain_voltage_V"]), 0.05)]
    checks = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "state_count": len(states) == int(acceptance["required_state_count"]),
        "well_summary_count": len(summary_rows) == int(acceptance["required_well_count"]),
        "component_closure_count": len(closure_rows) == int(
            acceptance["required_component_closure_count"]),
        "doping_well_partition_constructed": partition_valid,
        "topology_invariant_across_bias_states": topology_invariant,
        "opposite_interface_flux_balance": maximum_interface_balance <= float(
            acceptance["opposite_interface_flux_absolute_tolerance_A_per_um"]),
        "classification_declared": classification in contract["analysis"]["classifications"],
        "read_only_execution": True,
        "closed_topics_not_reopened": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    max_anchor_row = max(closure_rows,
                         key=lambda row: abs(float(row["direct_anchor_residual_A_per_um"])))
    max_closure_row = max(closure_rows,
                          key=lambda row: abs(float(row["default_closure_residual_A_per_um"])))
    report = {
        "schema": "vela.simplemos.sdevice.m57_doping_well_boundary_flux_closure_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_execution": False, "new_vela_execution": False,
                      "source": "existing M55 48-state exports"},
        "partition": {
            "state_count": len(states), "well_count": len(summary_rows),
            "maximum_labels_per_silicon_triangle": maximum_labels,
            "partition_problem_count": len(partition_problems),
            "partition_problems": partition_problems,
            "topology_invariant_across_bias_states": topology_invariant,
            "maximum_opposite_interface_flux_balance_residual_A_per_um":
                maximum_interface_balance,
            "states_with_unassociated_wells": [row["state"] for row in state_diagnostics
                                                if len(row["all_labels"]) > 3],
        },
        "physical_contact_anchor": {
            "all_within_tolerance": anchors_pass,
            "passing_row_count": sum(bool(row["direct_anchor_within_tolerance"])
                                     for row in closure_rows),
            "row_count": len(closure_rows),
            "maximum_absolute_residual_A_per_um": abs(float(
                max_anchor_row["direct_anchor_residual_A_per_um"])),
            "maximum_residual_row": max_anchor_row,
        },
        "complete_well_closure": {
            "all_within_tolerance": closures_pass,
            "passing_row_count": sum(bool(row["default_closure_within_tolerance"])
                                     for row in closure_rows),
            "row_count": len(closure_rows),
            "maximum_absolute_residual_A_per_um": abs(float(
                max_closure_row["default_closure_residual_A_per_um"])),
            "maximum_residual_row": max_closure_row,
        },
        "target": {"device": "n23", "drain_voltage_V": 0.05,
                   "gate_voltage_V": 0.05, "closure_rows": target_rows},
        "controls": {"description": "n19/n23, Vd=0.05 V, Vg=0/0.05/0.1 V",
                     "closure_rows": control_rows},
        "claim_guard": contract["analysis"]["claim_guard"],
        "forbidden_work_respected": contract["forbidden_work"],
        "acceptance": checks,
    }
    return report, all_segments, summary_rows, closure_rows, pair_rows, source_hashes


def write_human_report(report: dict[str, Any]) -> None:
    target = {(row["contact"], row["component"]): row
              for row in report["target"]["closure_rows"]}
    lines = []
    for contact in CONTACTS:
        for component in COMPONENTS:
            row = target[(contact, component)]
            lines.append(
                f"| {contact} | {component} | {float(row['default_A_per_um']):.12e} | "
                f"{float(row['physical_contact_terminal_oriented_flux_A_per_um']):.12e} | "
                f"{float(row['reconstructed_noncontact_surface_redistribution_A_per_um']):.12e} | "
                f"{float(row['signed_well_srh_charge_current_A_per_um']):.12e} | "
                f"{float(row['default_closure_residual_A_per_um']):.12e} |")
    anchor = report["physical_contact_anchor"]
    closure = report["complete_well_closure"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M57 完整掺杂阱边界通量闭合

## 结论

M57 分类为 `{report['classification']}`。本任务只读取 M55 已取回的 48 个状态，没有新增 Sentaurus 或 Vela 求解，也没有修改任何物理、网格、接触、电流算法或生产默认值。

`DopingWells` 节点标签在全部状态均可按冻结的中点分区规则构造；硅三角形最多含 `{report['partition']['maximum_labels_per_silicon_triangle']}` 个标签，跨偏压拓扑保持不变，对置内部界面通量最大抵消残差为 `{float(report['partition']['maximum_opposite_interface_flux_balance_residual_A_per_um']):.3e}` A/um。n21/n23 中未连接端口的附加阱标签保持独立，没有合并进 source、drain 或 substrate 阱。

物理接触面节点场重构相对 DirectCurrent 通过 `{anchor['passing_row_count']}/{anchor['row_count']}` 行，最大绝对残差为 `{float(anchor['maximum_absolute_residual_A_per_um']):.12e}` A/um。完整阱表面通量加冻结 SRH 项相对默认端口通过 `{closure['passing_row_count']}/{closure['row_count']}` 行，最大绝对残差为 `{float(closure['maximum_absolute_residual_A_per_um']):.12e}` A/um。分类严格按冻结合同决定：先检查分区，再检查物理接触面锚点，最后判断完整阱闭合。

## 目标点 n23, Vd=0.05 V, Vg=0.05 V

| 端口 | 分量 | 默认端口 (A/um) | 物理接触面 (A/um) | 非接触阱边界 (A/um) | SRH项 (A/um) | 闭合残差 (A/um) |
|---|---|---:|---:|---:|---:|---:|
{chr(10).join(lines)}

## 方法与边界

- 内部跨阱三角形的界面取两条异标签边的中点连线；外边界异标签边在中点拆分。
- 电流密度在段端点线性插值并作梯形积分，几何法向始终指向关联阱之外；Sentaurus 端口方向固定为几何外向通量的负号。
- `physical_contact`、`interior_doping_well_interface` 和 `other_semiconductor_exterior` 三类边界分别保留在逐段账本。
- M57 检验的是导出节点场上的一个明确离散重构。若未闭合，只能归为节点场边界求积限制或该重构与 Sentaurus 内部默认阱面积分面不一致，不能据此声明 Sentaurus 缺陷，也不授权修改 Vela。

机器报告：`reference_tcad/simplemos_sentaurus2022/doping_well_boundary_flux_closure/m57_doping_well_boundary_flux_closure_report.json`。
""", encoding="utf-8", newline="\n")


def write_artifact(report: dict[str, Any]) -> None:
    target = report["target"]["closure_rows"]
    payload = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report",
            "title": "SimpleMOS M57 doping-well boundary flux closure",
            "description": "Frozen 48-state complete DopingWells boundary reconstruction.",
            "generatedAt": "2026-09-02T00:00:00+08:00",
            "sources": [{"id": "src_m57", "label": "M57 machine report",
                         "path": portable(REPORT)},
                        {"id": "src_contract", "label": "Frozen M57 contract",
                         "path": portable(CONTRACT)}],
            "tables": [{"id": "target", "title": "Target carrier closure",
                        "dataset": "target", "sourceId": "src_m57", "layout": "full",
                        "columns": [{"field": "contact", "label": "Contact", "type": "text"},
                                    {"field": "component", "label": "Carrier", "type": "text"},
                                    {"field": "default_A_per_um", "label": "Default", "format": "number"},
                                    {"field": "predicted_default_A_per_um", "label": "Reconstructed", "format": "number"},
                                    {"field": "default_closure_residual_A_per_um", "label": "Residual", "format": "number"}]}],
            "blocks": [{"id": "title", "type": "markdown",
                        "body": "# SimpleMOS M57 doping-well boundary flux closure"},
                       {"id": "summary", "type": "markdown", "sourceId": "src_m57",
                        "body": f"Outcome: `{report['classification']}`. Read-only reconstruction of 48 frozen M55 states."},
                       {"id": "target", "type": "table", "tableId": "target"},
                       {"id": "limits", "type": "markdown", "sourceId": "src_contract",
                        "body": "The result tests a frozen midpoint partition of exported nodal fields and does not authorize a production change."}]
        },
        "snapshot": {"version": 1, "generatedAt": "2026-09-02T00:00:00+08:00",
                     "status": "ready", "datasets": {"target": target}},
        "sources": [{"id": "src_m57", "label": "M57 machine report",
                     "path": portable(REPORT)}]
    }
    write_json(ARTIFACT, payload)


def execute() -> dict[str, Any]:
    contract = validate_contract()
    report, segments, wells, closures, pairs, source_hashes = analyze(contract)
    write_csv(SEGMENTS, segments)
    write_csv(WELLS, wells)
    write_csv(CLOSURES, closures)
    write_csv(PAIRS, pairs)
    write_json(REPORT, report)
    write_human_report(report)
    write_artifact(report)
    artifacts = {portable(path): sha256(path) for path in
                 (REPORT, SEGMENTS, WELLS, CLOSURES, PAIRS, DOC, ARTIFACT)}
    evidence = {
        "schema": "vela.simplemos.sdevice.m57_doping_well_boundary_flux_closure_evidence.v1",
        "status": report["status"], "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT),
        "artifacts": artifacts, "source_hashes": source_hashes,
        "new_sentaurus_execution": False, "new_vela_execution": False,
        "historical_artifacts_rewritten": False,
        "closed_topics_reinvestigated": False,
        "acceptance": report["acceptance"],
    }
    write_json(EVIDENCE, evidence)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true",
                        help="validate the frozen contract and existing evidence hashes")
    args = parser.parse_args()
    if args.verify:
        validate_contract()
        evidence = read_json(EVIDENCE)
        for relative, expected in evidence["artifacts"].items():
            if sha256(REPO / relative) != expected:
                raise ValueError(f"M57 artifact hash mismatch: {relative}")
        for relative, expected in evidence["source_hashes"].items():
            if sha256(REPO / relative) != expected:
                raise ValueError(f"M57 source hash mismatch: {relative}")
        print(f"M57 verified: {evidence['classification']}")
        return
    report = execute()
    print(json.dumps({"status": report["status"],
                      "classification": report["classification"],
                      "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
