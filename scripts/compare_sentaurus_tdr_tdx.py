#!/usr/bin/env python3
"""Compare Sentaurus TDR neutral exports with tdx-generated DF-ISE text files."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def read_json_if_exists(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {}
    return json.loads(path.read_text())


def parse_count(text: str, name: str) -> int | None:
    match = re.search(rf"\b{name}\s*=\s*(\d+)", text)
    return int(match.group(1)) if match else None


def parse_ise_list(text: str, name: str) -> list[str]:
    match = re.search(rf"\b{name}\s*=\s*\[(.*?)\]", text, re.S)
    if not match:
        return []
    raw = match.group(1)
    return [token[0] or token[1] for token in re.findall(r'"([^"]+)"|([A-Za-z0-9_.:+-]+)', raw)]


def parse_grd(path: Path) -> dict[str, Any]:
    text = path.read_text(errors="replace")
    regions = parse_ise_list(text, "regions")
    materials = parse_ise_list(text, "materials")
    region_materials = {
        region: materials[index]
        for index, region in enumerate(regions)
        if index < len(materials)
    }
    region_blocks: dict[str, dict[str, Any]] = {}
    matches = list(re.finditer(r'Region\s*\("([^"]+)"\)\s*\{', text))
    for index, match in enumerate(matches):
        name = match.group(1)
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        body = text[match.end():end]
        material_match = re.search(r"\bmaterial\s*=\s*([A-Za-z0-9_.:+-]+)", body)
        elements_match = re.search(r"\bElements\s*\((\d+)\)", body)
        material = material_match.group(1) if material_match else region_materials.get(name)
        region_blocks[name] = {
            "material": material,
            "elements": int(elements_match.group(1)) if elements_match else None,
        }
    contact_counts = {
        name: int(info["elements"])
        for name, info in region_blocks.items()
        if str(info.get("material", "")).lower() == "contact" and info.get("elements") is not None
    }
    vertices = parse_numeric_rows(text, "Vertices")
    edges = [[int(value) for value in row] for row in parse_numeric_rows(text, "Edges")]
    elements = [[int(value) for value in row] for row in parse_numeric_rows(text, "Elements")]
    region_elements = parse_region_elements(text)
    return {
        "path": str(path),
        "vertex_count": parse_count(text, "nb_vertices"),
        "element_count": parse_count(text, "nb_elements"),
        "region_count": parse_count(text, "nb_regions"),
        "regions": regions,
        "materials": materials,
        "region_blocks": region_blocks,
        "contact_element_counts": contact_counts,
        "vertices": vertices,
        "edges": edges,
        "elements": elements,
        "region_elements": region_elements,
    }


def parse_numeric_rows(text: str, block_name: str) -> list[list[float]]:
    match = re.search(rf"\b{block_name}\s*\((\d+)\)\s*\{{(?P<body>.*?)\n\s*\}}", text, re.S)
    if not match:
        return []
    rows: list[list[float]] = []
    for line in match.group("body").splitlines():
        values = [float(item) for item in NUMBER_RE.findall(line)]
        if values:
            rows.append(values)
    return rows


def parse_region_elements(text: str) -> dict[str, list[int]]:
    elements_by_region: dict[str, list[int]] = {}
    matches = list(re.finditer(r'Region\s*\("([^"]+)"\)\s*\{', text))
    for index, match in enumerate(matches):
        name = match.group(1)
        end = matches[index + 1].start() if index + 1 < len(matches) else len(text)
        body = text[match.end():end]
        elements_match = re.search(r"\bElements\s*\((\d+)\)\s*\{(?P<values>.*?)\n\s*\}", body, re.S)
        if not elements_match:
            elements_by_region[name] = []
            continue
        elements_by_region[name] = [int(item) for item in NUMBER_RE.findall(elements_match.group("values"))]
    return elements_by_region


def parse_dat(path: Path) -> dict[str, dict[str, Any]]:
    text = path.read_text(errors="replace")
    datasets: dict[str, dict[str, Any]] = {}
    pattern = re.compile(
        r'Dataset\s*\("([^"]+)"\)\s*\{(?P<header>.*?)\bValues\s*\((\d+)\)\s*\{(?P<values>.*?)\n\s*\}',
        re.S,
    )
    for match in pattern.finditer(text):
        name = match.group(1)
        header = match.group("header")
        dimension_match = re.search(r"\bdimension\s*=\s*(\d+)", header)
        dimension = int(dimension_match.group(1)) if dimension_match else 1
        values = [float(item) for item in NUMBER_RE.findall(match.group("values"))]
        count = int(match.group(3))
        datasets[name] = {
            "count": count,
            "dimension": dimension,
            "values": values,
        }
    return datasets


def count_csv_rows(path: Path) -> int | None:
    if not path.is_file():
        return None
    with path.open(newline="") as handle:
        return max(sum(1 for _ in handle) - 1, 0)


def load_field_csv(path: Path) -> list[list[float]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        component_columns = sorted(
            (name for name in reader.fieldnames if name.startswith("component")),
            key=lambda name: int(name.removeprefix("component") or "0"),
        )
        if not component_columns:
            component_columns = [name for name in reader.fieldnames if name != "node_id"]
        rows = []
        for row in reader:
            rows.append([float(row[name]) for name in component_columns])
        return rows


def find_tdr_field(root: Path, field: str) -> Path | None:
    fields_root = root / "fields"
    exact = fields_root / f"{field}_region0.csv"
    if exact.is_file():
        return exact
    matches = sorted(fields_root.glob(f"{field}_region*.csv"))
    return matches[0] if matches else None


def compare_flat_values(tdr_rows: list[list[float]], tdx_dataset: dict[str, Any]) -> dict[str, Any]:
    tdr_values = [component for row in tdr_rows for component in row]
    tdx_values = list(tdx_dataset.get("values", []))
    points = min(len(tdr_values), len(tdx_values))
    max_abs = 0.0
    max_rel = 0.0
    worst: dict[str, Any] | None = None
    for index in range(points):
        abs_diff = abs(tdr_values[index] - tdx_values[index])
        rel_diff = abs_diff / max(abs(tdx_values[index]), 1.0e-300)
        if abs_diff > max_abs:
            max_abs = abs_diff
            worst = {
                "index": index,
                "tdr": tdr_values[index],
                "tdx": tdx_values[index],
                "abs_diff": abs_diff,
                "rel_diff": rel_diff,
            }
        max_rel = max(max_rel, rel_diff)
    return {
        "points_compared": points,
        "tdr_values": len(tdr_values),
        "tdx_values": len(tdx_values),
        "max_abs_diff": max_abs,
        "max_rel_diff": max_rel,
        "worst": worst,
    }


def tdr_contact_counts(root: Path) -> dict[str, int]:
    metadata = read_json_if_exists(root / "metadata.json")
    counts: dict[str, int] = {}
    for region in metadata.get("regions", []):
        edges = region.get("edges")
        material = str(region.get("material", ""))
        if edges is not None and (material.lower() == "contact" or int(region.get("triangles", 0)) == 0):
            counts[str(region.get("name"))] = int(edges)
    return counts


def status_entry(tdr: int | None, tdx: int | None) -> dict[str, Any]:
    return {"tdr": tdr, "tdx": tdx, "match": tdr == tdx if tdr is not None and tdx is not None else None}


def non_contact_element_count(grd: dict[str, Any]) -> int | None:
    total = 0
    found = False
    for info in grd.get("region_blocks", {}).values():
        if str(info.get("material", "")).lower() == "contact":
            continue
        elements = info.get("elements")
        if elements is None:
            continue
        total += int(elements)
        found = True
    return total if found else grd.get("element_count")


def tdr_total_element_count(root: Path, contact_counts: dict[str, int]) -> int | None:
    bulk_count = count_csv_rows(root / "elements.csv")
    if bulk_count is None:
        return None
    return bulk_count + sum(contact_counts.values())


def read_tdr_nodes(root: Path) -> list[dict[str, Any]]:
    inventory = read_tdr_geometry_inventory(root)
    vertices = inventory.get("vertices", [])
    if vertices:
        return [
            {"id": index, "x_um": float(vertex[0]), "y_um": float(vertex[1])}
            for index, vertex in enumerate(vertices)
        ]
    with (root / "nodes.csv").open(newline="") as handle:
        rows = [
            {"id": int(row["id"]), "x_um": float(row["x_um"]), "y_um": float(row["y_um"])}
            for row in csv.DictReader(handle)
        ]
    return sorted(rows, key=lambda row: row["id"])


def read_tdr_bulk_elements(root: Path) -> list[dict[str, Any]]:
    inventory = read_tdr_geometry_inventory(root)
    regions = inventory.get("regions", [])
    if regions:
        rows = []
        cell_id = 0
        for region in sorted(regions, key=lambda item: int(item.get("index", 0))):
            if int(region.get("type", 99)) != 0:
                continue
            for triangle in region.get("triangles", []):
                rows.append({
                    "id": cell_id,
                    "nodes": sorted(int(node) for node in triangle),
                    "region": str(region.get("name", "")),
                    "material": str(region.get("material", "")),
                })
                cell_id += 1
        return rows
    with (root / "elements.csv").open(newline="") as handle:
        rows = []
        for row in csv.DictReader(handle):
            rows.append({
                "id": int(row["id"]),
                "nodes": sorted([int(row["node0"]), int(row["node1"]), int(row["node2"])]),
                "region": row.get("region", ""),
                "material": row.get("material", ""),
            })
    return sorted(rows, key=lambda row: row["id"])


def read_tdr_contacts(root: Path) -> dict[str, dict[str, Any]]:
    inventory = read_tdr_geometry_inventory(root)
    regions = inventory.get("regions", [])
    if regions:
        contacts: dict[str, dict[str, Any]] = {}
        for region in regions:
            if int(region.get("type", 99)) != 1:
                continue
            node_ids = sorted({int(node) for edge in region.get("edges", []) for node in edge})
            contacts[str(region.get("name", ""))] = {
                "node_ids": node_ids,
                "region": "",
                "edge_count": len(region.get("edges", [])),
            }
        return contacts
    contacts: dict[str, dict[str, Any]] = {}
    contact_path = root / "contacts.csv"
    if not contact_path.is_file():
        return contacts
    with contact_path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            name = str(row["name"])
            raw_nodes = str(row.get("node_ids", ""))
            contacts[name] = {
                "node_ids": sorted(int(item) for item in raw_nodes.split(";") if item != ""),
                "region": row.get("region", ""),
            }
    return contacts


def read_tdr_geometry_inventory(root: Path) -> dict[str, Any]:
    candidates = [
        root / "tdr_inventory" / "mesh.json",
        root / "tdr_inventory.json",
    ]
    for path in candidates:
        data = read_json_if_exists(path)
        geometry = data.get("geometry")
        if isinstance(geometry, dict):
            return geometry
    return {}


def decode_grd_element_nodes(record: list[int], edges: list[list[int]]) -> list[int]:
    if not record:
        return []
    kind = record[0]
    if kind == 1:
        if len(record) == 3:
            return sorted(record[1:])
        if len(record) == 2:
            return sorted(edges[abs(record[1])])
    if kind == 2:
        nodes: set[int] = set()
        for edge_ref in record[1:]:
            edge_index = edge_ref if edge_ref >= 0 else -edge_ref - 1
            edge = edges[edge_index]
            nodes.update(edge)
        return sorted(nodes)
    return sorted(record[1:])


def non_contact_regions(grd: dict[str, Any]) -> list[str]:
    names = []
    for name, info in grd.get("region_blocks", {}).items():
        if str(info.get("material", "")).lower() != "contact":
            names.append(name)
    return sorted(names)


def contact_regions(grd: dict[str, Any]) -> list[str]:
    names = []
    for name, info in grd.get("region_blocks", {}).items():
        if str(info.get("material", "")).lower() == "contact":
            names.append(name)
    return sorted(names)


def compare_geometry(tdr_export: Path, grd: dict[str, Any], coordinate_epsilon: float) -> dict[str, Any]:
    tdr_nodes = read_tdr_nodes(tdr_export)
    grd_vertices = grd.get("vertices", [])
    nodes_report = compare_geometry_nodes(tdr_nodes, grd_vertices, coordinate_epsilon)

    tdr_elements = read_tdr_bulk_elements(tdr_export)
    elements_report = compare_geometry_elements(tdr_elements, grd)

    tdr_contacts = read_tdr_contacts(tdr_export)
    boundaries_report = compare_geometry_boundaries(tdr_contacts, grd)

    failed = (
        not nodes_report["count_match"]
        or nodes_report["status"] != "pass"
        or not elements_report["count_match"]
        or elements_report["status"] != "pass"
        or boundaries_report["status"] != "pass"
    )
    return {
        "status": "fail" if failed else "pass",
        "nodes": nodes_report,
        "elements": elements_report,
        "boundaries": boundaries_report,
    }


def compare_geometry_nodes(tdr_nodes: list[dict[str, Any]],
                           grd_vertices: list[list[float]],
                           coordinate_epsilon: float) -> dict[str, Any]:
    points = min(len(tdr_nodes), len(grd_vertices))
    max_abs = 0.0
    max_rel = 0.0
    worst: dict[str, Any] | None = None
    failed = len(tdr_nodes) != len(grd_vertices)
    for index in range(points):
        tdr_node = tdr_nodes[index]
        if int(tdr_node["id"]) != index:
            failed = True
            worst = worst or {
                "node_id": int(tdr_node["id"]),
                "component": "id",
                "tdr": int(tdr_node["id"]),
                "grd": index,
                "abs_diff": abs(int(tdr_node["id"]) - index),
                "rel_diff": None,
            }
            continue
        for component_index, component in enumerate(("x_um", "y_um")):
            tdr_value = float(tdr_node[component])
            grd_value = float(grd_vertices[index][component_index])
            abs_diff = abs(tdr_value - grd_value)
            rel_diff = abs_diff / abs(grd_value) if grd_value != 0.0 else 0.0
            threshold = coordinate_epsilon
            component_failed = abs_diff > threshold if grd_value == 0.0 else rel_diff > threshold
            if component_failed:
                failed = True
            if abs_diff > max_abs or (abs_diff == max_abs and rel_diff > max_rel):
                max_abs = abs_diff
                max_rel = rel_diff
                worst = {
                    "node_id": index,
                    "component": component,
                    "tdr": tdr_value,
                    "grd": grd_value,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_diff,
                }
            max_rel = max(max_rel, rel_diff)
    return {
        "status": "fail" if failed else "pass",
        "tdr_count": len(tdr_nodes),
        "grd_count": len(grd_vertices),
        "count_match": len(tdr_nodes) == len(grd_vertices),
        "epsilon": coordinate_epsilon,
        "max_abs_diff": max_abs,
        "max_rel_diff": max_rel,
        "worst": worst,
    }


def compare_geometry_elements(tdr_elements: list[dict[str, Any]], grd: dict[str, Any]) -> dict[str, Any]:
    elements = grd.get("elements", [])
    edges = grd.get("edges", [])
    bulk_ids: list[int] = []
    for region in non_contact_regions(grd):
        bulk_ids.extend(grd.get("region_elements", {}).get(region, []))
    bulk_ids = sorted(bulk_ids)
    first_mismatch: dict[str, Any] | None = None
    points = min(len(tdr_elements), len(bulk_ids))
    for index in range(points):
        tdr = tdr_elements[index]
        grd_element_id = bulk_ids[index]
        grd_nodes = decode_grd_element_nodes(elements[grd_element_id], edges)
        if int(tdr["id"]) != grd_element_id or tdr["nodes"] != grd_nodes:
            first_mismatch = {
                "index": index,
                "tdr_id": int(tdr["id"]),
                "grd_id": grd_element_id,
                "tdr_nodes": tdr["nodes"],
                "grd_nodes": grd_nodes,
            }
            break
    region_membership = {}
    for region in non_contact_regions(grd):
        grd_ids = sorted(grd.get("region_elements", {}).get(region, []))
        tdr_ids = sorted(int(row["id"]) for row in tdr_elements if row["region"] == region)
        region_membership[region] = {
            "tdr": tdr_ids,
            "grd": grd_ids,
            "match": tdr_ids == grd_ids,
        }
    count_match = len(tdr_elements) == len(bulk_ids)
    membership_match = all(entry["match"] for entry in region_membership.values())
    return {
        "status": "pass" if count_match and first_mismatch is None and membership_match else "fail",
        "tdr_count": len(tdr_elements),
        "grd_count": len(bulk_ids),
        "count_match": count_match,
        "first_mismatch": first_mismatch,
        "region_membership": region_membership,
    }


def compare_geometry_boundaries(tdr_contacts: dict[str, dict[str, Any]], grd: dict[str, Any]) -> dict[str, Any]:
    elements = grd.get("elements", [])
    edges = grd.get("edges", [])
    grd_contacts = contact_regions(grd)
    contact_names = {
        "tdr": sorted(tdr_contacts),
        "grd": sorted(grd_contacts),
        "match": sorted(tdr_contacts) == sorted(grd_contacts),
    }
    contacts: dict[str, Any] = {}
    failed = not contact_names["match"]
    for name in sorted(set(tdr_contacts) | set(grd_contacts)):
        element_ids = sorted(grd.get("region_elements", {}).get(name, []))
        grd_nodes: set[int] = set()
        for element_id in element_ids:
            grd_nodes.update(decode_grd_element_nodes(elements[element_id], edges))
        tdr_nodes = tdr_contacts.get(name, {}).get("node_ids", [])
        tdr_edge_count = tdr_contacts.get(name, {}).get("edge_count")
        node_set_match = sorted(tdr_nodes) == sorted(grd_nodes)
        expected_edge_count = int(tdr_edge_count) if tdr_edge_count is not None else int(
            grd.get("region_blocks", {}).get(name, {}).get("elements", len(element_ids)))
        element_count_match = len(element_ids) == expected_edge_count
        if not node_set_match or not element_count_match:
            failed = True
        contacts[name] = {
            "tdr_node_ids": sorted(tdr_nodes),
            "grd_node_ids": sorted(grd_nodes),
            "node_set_match": node_set_match,
            "tdr_edge_count": tdr_edge_count,
            "grd_element_ids": element_ids,
            "grd_element_count": len(element_ids),
            "element_count_match": element_count_match,
        }
    return {
        "status": "fail" if failed else "pass",
        "contact_names": contact_names,
        "contacts": contacts,
    }


def choose_first(root: Path, pattern: str, preferred: list[str] | None = None) -> Path:
    if preferred:
        for name in preferred:
            candidate = root / name
            if candidate.is_file():
                return candidate
    matches = sorted(root.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"no {pattern} found in {root}")
    return matches[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdr-export", required=True, type=Path)
    parser.add_argument("--tdx-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--geometry-only", action="store_true",
                        help="Compare only TDR/GRD nodes, elements, and contact boundaries.")
    parser.add_argument("--coordinate-epsilon", type=float, default=sys.float_info.epsilon,
                        help="Absolute/relative coordinate tolerance for --geometry-only checks.")
    parser.add_argument("--max-abs-diff", type=float, default=0.0)
    parser.add_argument("--max-rel-diff", type=float, default=1.0e-12)
    args = parser.parse_args()

    grd = parse_grd(choose_first(args.tdx_dir, "*.grd", ["pn2d_msh.grd"]))
    if args.geometry_only:
        report = compare_geometry(args.tdr_export, grd, args.coordinate_epsilon)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        (args.output_dir / "tdr_tdx_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
        return 1 if report["status"] == "fail" else 0

    dat = parse_dat(choose_first(args.tdx_dir, "*.dat", ["pn2d_msh.dat"]))
    metadata = read_json_if_exists(args.tdr_export / "metadata.json")
    tdr_regions = [
        str(region.get("name"))
        for region in metadata.get("regions", [])
        if region.get("name") is not None
    ]
    if not tdr_regions and (args.tdr_export / "elements.csv").is_file():
        with (args.tdr_export / "elements.csv").open(newline="") as handle:
            reader = csv.DictReader(handle)
            tdr_regions = sorted({row.get("region", "") for row in reader if row.get("region")})
    if (args.tdr_export / "contacts.csv").is_file():
        with (args.tdr_export / "contacts.csv").open(newline="") as handle:
            reader = csv.DictReader(handle)
            contact_names = {row.get("name", "") for row in reader if row.get("name")}
        tdr_regions = sorted(set(tdr_regions) | contact_names)
    contact_counts = tdr_contact_counts(args.tdr_export)

    mesh = {
        "vertex_count": status_entry(
            int(metadata.get("vertex_count")) if metadata.get("vertex_count") is not None else count_csv_rows(args.tdr_export / "nodes.csv"),
            grd.get("vertex_count"),
        ),
        "element_count": status_entry(tdr_total_element_count(args.tdr_export, contact_counts), grd.get("element_count")),
        "bulk_element_count": status_entry(
            count_csv_rows(args.tdr_export / "elements.csv"),
            non_contact_element_count(grd),
        ),
        "region_names": {
            "tdr": sorted(tdr_regions),
            "tdx": sorted(grd.get("regions", [])),
            "match": sorted(tdr_regions) == sorted(grd.get("regions", [])) if tdr_regions else None,
        },
        "contact_element_counts": {},
    }
    for name in sorted(set(contact_counts) | set(grd.get("contact_element_counts", {}))):
        mesh["contact_element_counts"][name] = {
            "tdr": contact_counts.get(name),
            "tdx": grd.get("contact_element_counts", {}).get(name),
            "match": (
                contact_counts.get(name) == grd.get("contact_element_counts", {}).get(name)
                if name in contact_counts
                else None
            ),
        }

    fields: dict[str, Any] = {}
    for field, dataset in sorted(dat.items()):
        field_path = find_tdr_field(args.tdr_export, field)
        if field_path is None:
            continue
        report = compare_flat_values(load_field_csv(field_path), dataset)
        report["status"] = "pass"
        if (
            report["tdr_values"] != report["tdx_values"]
            or (
                report["max_abs_diff"] > args.max_abs_diff
                and report["max_rel_diff"] > args.max_rel_diff
            )
        ):
            report["status"] = "fail"
        fields[field] = report

    failed = any(
        entry.get("match") is False
        for key, entry in mesh.items()
        if isinstance(entry, dict) and key != "contact_element_counts"
    )
    failed = failed or any(report.get("status") == "fail" for report in fields.values())
    report = {
        "status": "fail" if failed else "pass",
        "mesh": mesh,
        "fields": fields,
        "tdx": {
            "grd": grd["path"],
            "datasets": sorted(dat),
        },
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "tdr_tdx_comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
