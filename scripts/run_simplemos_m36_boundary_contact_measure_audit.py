#!/usr/bin/env python3
"""Audit n23 Sentaurus AverageBox measures by boundary/contact class."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import run_simplemos_m34_sentaurus_interface_box_probe as m34


CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m36_boundary_contact_measure_audit_contract_v1.json")
BUILD = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
         / "m34_interface_box_probe")
IMPORTED = BUILD / "n23_import"
DEBUG = BUILD / "sentaurus_raw/sentaurus_bundle/n23/MeasureCoefficients.debug"
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "boundary_contact_measure_audit")
MEASURE_PERMUTATION = (0, 2, 1)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def relative_error(candidate: float, reference: float) -> float:
    return abs(candidate - reference) / max(abs(reference), 1e-300)


def summary(values: list[float]) -> dict[str, float | int]:
    ordered = sorted(values)
    return {
        "count": len(ordered),
        "minimum": min(ordered),
        "median": statistics.median(ordered),
        "p95": ordered[min(len(ordered) - 1, math.ceil(0.95 * len(ordered)) - 1)],
        "maximum": max(ordered),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--build-root", type=Path, default=BUILD)
    args = parser.parse_args()
    contract = read_json(args.contract)
    if contract["schema"] != "vela.simplemos.sdevice.m36_boundary_contact_measure_audit.v1":
        raise ValueError("unexpected M36 contract schema")

    imported = args.build_root / "n23_import"
    debug_path = args.build_root / "sentaurus_raw/sentaurus_bundle/n23/MeasureCoefficients.debug"
    nodes = {int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
             for row in read_csv(imported / "nodes.csv")}
    elements = {
        int(row["id"]): {
            "nodes": [int(row[f"node{local}"]) for local in range(3)],
            "region": row["region"],
            "material": row["material"],
        }
        for row in read_csv(imported / "elements.csv")
    }
    contacts = {
        row["name"]: {int(value) for value in row["node_ids"].split(";") if value}
        for row in read_csv(imported / "contacts.csv")
    }
    contact_by_node: dict[int, set[str]] = defaultdict(set)
    for name, support in contacts.items():
        for node in support:
            contact_by_node[node].add(name)

    edge_cells: dict[tuple[int, int], list[int]] = defaultdict(list)
    node_neighbors: dict[int, set[int]] = defaultdict(set)
    node_elements: dict[int, list[int]] = defaultdict(list)
    for element, item in elements.items():
        ids = item["nodes"]
        for node in ids:
            node_elements[node].append(element)
        for local in range(3):
            a, b = ids[local], ids[(local + 1) % 3]
            edge_cells[edge_key(a, b)].append(element)
            node_neighbors[a].add(b)
            node_neighbors[b].add(a)

    external_edges = {edge for edge, adjacent in edge_cells.items()
                      if len(adjacent) == 1}
    material_edges = {
        edge for edge, adjacent in edge_cells.items()
        if len(adjacent) == 2
        and len({elements[element]["material"] for element in adjacent}) > 1
    }
    contact_edges: dict[str, set[tuple[int, int]]] = {}
    for name, support in contacts.items():
        contact_edges[name] = {
            edge for edge in external_edges
            if edge[0] in support and edge[1] in support
        }
    contact_first_layer = {
        name: {
            neighbour for node in support for neighbour in node_neighbors[node]
            if neighbour not in support
        }
        for name, support in contacts.items()
    }
    external_nodes = {node for edge in external_edges for node in edge}
    material_interface_nodes = {node for edge in material_edges for node in edge}

    text = debug_path.read_text(encoding="utf-8", errors="replace")
    measures = m34.parse_debug_block(text, "Measure")
    coefficients = m34.parse_debug_block(text, "Coefficients")

    local_rows: list[dict[str, Any]] = []
    node_accumulator: dict[int, dict[str, float]] = defaultdict(
        lambda: {"sentaurus_si": 0.0, "raw_si": 0.0,
                 "sentaurus_all": 0.0, "raw_all": 0.0})
    for element, item in elements.items():
        ids = item["nodes"]
        points = [nodes[node] for node in ids]
        raw = m34.signed_average_box_measures(points)
        sent = [float(measures[element]["values"][MEASURE_PERMUTATION[local]])
                for local in range(3)]
        tri_edges = {edge_key(ids[local], ids[(local + 1) % 3])
                     for local in range(3)}
        direct_names = sorted({name for name, support in contacts.items()
                               if any(node in support for node in ids)})
        edge_names = sorted({name for name, support in contact_edges.items()
                             if tri_edges & support})
        touches_external = bool(tri_edges & external_edges)
        touches_material = bool(tri_edges & material_edges)
        if edge_names:
            element_class = "direct_contact"
        elif direct_names:
            element_class = "contact_adjacent"
        elif touches_external:
            element_class = "external_boundary"
        elif touches_material:
            element_class = "material_interface"
        else:
            element_class = "interior"
        sent_sum = sum(sent)
        raw_sum = sum(raw)
        area = m34.area(points)
        for local, node in enumerate(ids):
            acc = node_accumulator[node]
            acc["sentaurus_all"] += sent[local]
            acc["raw_all"] += raw[local]
            if item["material"] == "Si":
                acc["sentaurus_si"] += sent[local]
                acc["raw_si"] += raw[local]
            local_rows.append({
                "element": element,
                "local_vertex": local,
                "node": node,
                "region": item["region"],
                "material": item["material"],
                "element_class": element_class,
                "direct_contact_names": ";".join(direct_names),
                "contact_edge_names": ";".join(edge_names),
                "touches_external_boundary": touches_external,
                "touches_material_interface": touches_material,
                "sentaurus_measure_um2": sent[local],
                "raw_signed_measure_um2": raw[local],
                "absolute_error_um2": abs(raw[local] - sent[local]),
                "relative_error": relative_error(raw[local], sent[local]),
                "sentaurus_element_sum_um2": sent_sum,
                "raw_element_sum_um2": raw_sum,
                "triangle_area_um2": area,
                "sentaurus_area_closure_fraction": sent_sum / area,
                "raw_area_closure_fraction": raw_sum / area,
                "minimum_debug_coefficient": min(coefficients[element]["values"]),
            })

    node_rows: list[dict[str, Any]] = []
    for node in sorted(nodes):
        direct = sorted(contact_by_node[node])
        first_layer = sorted(name for name, support in contact_first_layer.items()
                             if node in support)
        if direct:
            node_class = "direct_contact"
        elif first_layer:
            node_class = "contact_first_layer"
        elif node in external_nodes:
            node_class = "external_boundary"
        elif node in material_interface_nodes:
            node_class = "material_interface"
        else:
            node_class = "interior"
        acc = node_accumulator[node]
        has_si = any(elements[element]["material"] == "Si"
                     for element in node_elements[node])
        node_rows.append({
            "node": node,
            "x_um": nodes[node][0],
            "y_um": nodes[node][1],
            "node_class": node_class,
            "direct_contact_names": ";".join(direct),
            "first_layer_contact_names": ";".join(first_layer),
            "has_silicon": has_si,
            **acc,
            "si_absolute_error_um2": abs(acc["raw_si"] - acc["sentaurus_si"]),
            "si_relative_error": relative_error(acc["raw_si"], acc["sentaurus_si"])
                if has_si else 0.0,
            "all_absolute_error_um2": abs(acc["raw_all"] - acc["sentaurus_all"]),
            "all_relative_error": relative_error(acc["raw_all"], acc["sentaurus_all"]),
        })

    class_rows: list[dict[str, Any]] = []
    for node_class in sorted({row["node_class"] for row in node_rows}):
        subset = [row for row in node_rows
                  if row["node_class"] == node_class and row["has_silicon"]]
        if not subset:
            continue
        errors = [float(row["si_relative_error"]) for row in subset]
        class_rows.append({
            "node_class": node_class,
            **{f"relative_error_{key}": value
               for key, value in summary(errors).items()},
            "sentaurus_si_measure_sum_um2": sum(
                float(row["sentaurus_si"]) for row in subset),
            "raw_si_measure_sum_um2": sum(float(row["raw_si"]) for row in subset),
        })

    contact_rows: list[dict[str, Any]] = []
    for name in sorted(contacts):
        for support_class, support in (
            ("direct_contact", contacts[name]),
            ("first_layer", contact_first_layer[name]),
        ):
            subset = [row for row in node_rows
                      if int(row["node"]) in support and row["has_silicon"]]
            if not subset:
                continue
            errors = [float(row["si_relative_error"]) for row in subset]
            contact_rows.append({
                "contact": name,
                "support_class": support_class,
                **{f"relative_error_{key}": value
                   for key, value in summary(errors).items()},
                "sentaurus_si_measure_sum_um2": sum(
                    float(row["sentaurus_si"]) for row in subset),
                "raw_si_measure_sum_um2": sum(
                    float(row["raw_si"]) for row in subset),
            })

    interface_subset = [row for row in node_rows
                        if row["node_class"] == "material_interface"
                        and row["has_silicon"]]
    limits = contract["acceptance"]
    present_classes = {row["node_class"] for row in node_rows
                       if row["has_silicon"]}
    maximum_si_node_relative_error = max(
        float(row["si_relative_error"]) for row in node_rows if row["has_silicon"])
    checks = {
        "triangles": len(elements) == limits["triangle_count"],
        "source_nodes": len(contacts.get("source", set())) >= limits["minimum_source_nodes"],
        "drain_nodes": len(contacts.get("drain", set())) >= limits["minimum_drain_nodes"],
        "substrate_nodes": len(contacts.get("substrate", set())) >= limits["minimum_substrate_nodes"],
        "required_classes": set(contract["required_classes"]) <= present_classes,
        "assembled_si_node_measure": maximum_si_node_relative_error <= limits[
            "maximum_assembled_si_node_relative_error"],
        "finite": all(math.isfinite(float(row["relative_error"]))
                      for row in local_rows),
    }
    report = {
        "schema": "vela.simplemos.sdevice.m36_boundary_contact_measure_audit_report.v1",
        "status": "complete" if all(checks.values()) else "failed_acceptance",
        "execution": {
            "new_sentaurus_execution": False,
            "raw_oracle_reused": portable(debug_path),
            "cpp_changed": False,
            "default_model_changed": False,
        },
        "mesh": {
            "node_count": len(nodes),
            "triangle_count": len(elements),
            "external_edge_count": len(external_edges),
            "material_interface_edge_count": len(material_edges),
            "contacts": {name: len(support) for name, support in contacts.items()},
            "contact_edges": {name: len(support) for name, support in contact_edges.items()},
            "contact_first_layer_nodes": {
                name: len(support) for name, support in contact_first_layer.items()
            },
        },
        "acceptance": {"checks": checks, "all_checks_pass": all(checks.values())},
        "node_class_summary": class_rows,
        "contact_support_summary": contact_rows,
        "maximum_local_relative_error": max(
            float(row["relative_error"]) for row in local_rows),
        "maximum_si_node_relative_error": maximum_si_node_relative_error,
        "claim_policy": contract["claim_policy"],
        "artifacts": {
            "local_measure_ledger": portable(PORTABLE / "m36_local_measure_ledger.csv"),
            "node_measure_ledger": portable(PORTABLE / "m36_node_measure_ledger.csv"),
            "node_class_summary": portable(PORTABLE / "m36_node_class_summary.csv"),
            "contact_support_summary": portable(PORTABLE / "m36_contact_support_summary.csv"),
        },
    }
    PORTABLE.mkdir(parents=True, exist_ok=True)
    write_csv(PORTABLE / "m36_local_measure_ledger.csv", local_rows)
    write_csv(PORTABLE / "m36_node_measure_ledger.csv", node_rows)
    write_csv(PORTABLE / "m36_node_class_summary.csv", class_rows)
    write_csv(PORTABLE / "m36_contact_support_summary.csv", contact_rows)
    write_json(PORTABLE / "m36_boundary_contact_measure_audit_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "checks": checks,
        "node_class_summary": class_rows,
        "contact_support_summary": contact_rows,
    }, indent=2))
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
