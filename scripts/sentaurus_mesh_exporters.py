#!/usr/bin/env python3
"""Export every M1 format from one normalized mesh/field bundle."""

from __future__ import annotations

import csv
import json
import math
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


POINT_FIELD_NAMES = (
    "DonorConcentration",
    "AcceptorConcentration",
    "DopingConcentration",
    "AbsDopingConcentration",
    "Log10AbsDopingConcentration",
)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _contact_lines(bundle: dict[str, Any], contact: dict[str, Any]) -> list[dict[str, Any]]:
    return [line for line in bundle["lines"]
            if contact["physical_tag"] in line["physical_curve_tags"]]


def bundle_to_vela_mesh(bundle: dict[str, Any]) -> dict[str, Any]:
    nodes = [
        {"id": item["vela_id"], "x": item["x"], "y": item["y"]}
        for item in bundle["nodes"]
    ]
    triangles = [
        {"id": item["vela_id"], "region_id": item["region_id"],
         "node_ids": item["node_ids"]}
        for item in bundle["triangles"]
    ]
    material_aliases = bundle.get("provenance", {}).get("material_aliases", {})
    material_name = {
        item["id"]: material_aliases.get(item["name"], item["name"])
        for item in bundle["materials"]
    }
    regions = []
    for region in bundle["regions"]:
        cells = [item for item in bundle["triangles"] if item["region_id"] == region["id"]]
        material_ids = {item["material_id"] for item in cells}
        if len(material_ids) != 1:
            raise ValueError(f"region {region['name']!r} has ambiguous material ownership")
        regions.append({
            "id": region["id"], "name": region["name"],
            "material": material_name[next(iter(material_ids))],
            "cell_ids": [item["vela_id"] for item in cells],
        })
    contacts = []
    for contact in bundle["contacts"]:
        lines = _contact_lines(bundle, contact)
        node_ids = sorted({node for line in lines for node in line["node_ids"]})
        touching = {
            triangle["region_id"] for triangle in bundle["triangles"]
            if set(triangle["node_ids"]) & set(node_ids)
        }
        if len(touching) != 1:
            raise ValueError(f"contact {contact['name']!r} has ambiguous region ownership")
        contacts.append({
            "id": contact["id"], "name": contact["name"],
            "region_id": next(iter(touching)), "node_ids": node_ids,
        })
    return {
        "_comment": "Derived from vela.mesh_field_bundle.v1; coordinates are micrometres.",
        "nodes": nodes, "triangles": triangles, "regions": regions, "contacts": contacts,
    }


def write_doping_csv(path: Path, bundle: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = bundle["fields"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "node_id", "gmsh_tag", "x_um", "y_um",
            "donors_cm3", "acceptors_cm3", "net_doping_cm3",
        ])
        for node in bundle["nodes"]:
            index = node["vela_id"]
            writer.writerow([
                index, node["gmsh_tag"], repr(float(node["x"])), repr(float(node["y"])),
                repr(float(fields["DonorConcentration"]["values"][index])),
                repr(float(fields["AcceptorConcentration"]["values"][index])),
                repr(float(fields["DopingConcentration"]["values"][index])),
            ])


def _numbers(values: Iterable[Any]) -> str:
    return " ".join(str(value) if isinstance(value, int) else f"{float(value):.17g}"
                    for value in values)


def _vtk_cells(bundle: dict[str, Any]) -> tuple[list[list[int]], list[int], list[dict[str, int]]]:
    cells: list[list[int]] = []
    types: list[int] = []
    metadata: list[dict[str, int]] = []
    for triangle in bundle["triangles"]:
        cells.append(triangle["node_ids"])
        types.append(5)
        metadata.append({
            "RegionId": triangle["region_id"],
            "MaterialId": triangle["material_id"],
            "PhysicalSurfaceTag": triangle["physical_surface_tags"][0],
            "ElementId": triangle["vela_id"],
            "ContactId": -1,
            "PhysicalCurveTag": -1,
        })
    contact_id_by_tag = {item["physical_tag"]: item["id"] for item in bundle["contacts"]}
    for line in bundle["lines"]:
        cells.append(line["node_ids"])
        types.append(3)
        physical_tag = line["physical_curve_tags"][0]
        metadata.append({
            "RegionId": -1, "MaterialId": -1, "PhysicalSurfaceTag": -1,
            "ElementId": line["vela_id"],
            "ContactId": contact_id_by_tag.get(physical_tag, -1),
            "PhysicalCurveTag": physical_tag,
        })
    return cells, types, metadata


def write_vtu(path: Path, bundle: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cells, types, metadata = _vtk_cells(bundle)
    offsets: list[int] = []
    running = 0
    for cell in cells:
        running += len(cell)
        offsets.append(running)
    root = ET.Element("VTKFile", type="UnstructuredGrid", version="0.1", byte_order="LittleEndian")
    grid = ET.SubElement(root, "UnstructuredGrid")
    piece = ET.SubElement(grid, "Piece", NumberOfPoints=str(len(bundle["nodes"])),
                          NumberOfCells=str(len(cells)))
    point_data = ET.SubElement(piece, "PointData")
    for name in POINT_FIELD_NAMES:
        array = ET.SubElement(point_data, "DataArray", type="Float64", Name=name, format="ascii")
        array.text = _numbers(bundle["fields"][name]["values"])
    for name, field in sorted(bundle["fields"]["species"].items()):
        array = ET.SubElement(point_data, "DataArray", type="Float64", Name=name, format="ascii")
        array.text = _numbers(field["values"])
    cell_data = ET.SubElement(piece, "CellData")
    for name in ("RegionId", "MaterialId", "PhysicalSurfaceTag", "ElementId",
                 "ContactId", "PhysicalCurveTag"):
        array = ET.SubElement(cell_data, "DataArray", type="Int64", Name=name, format="ascii")
        array.text = _numbers(item[name] for item in metadata)
    points = ET.SubElement(piece, "Points")
    point_array = ET.SubElement(points, "DataArray", type="Float64",
                                NumberOfComponents="3", format="ascii")
    point_array.text = _numbers(
        value for node in bundle["nodes"] for value in (node["x"], node["y"], node["z"]))
    cells_node = ET.SubElement(piece, "Cells")
    for name, values, value_type in (
        ("connectivity", (value for cell in cells for value in cell), "Int64"),
        ("offsets", offsets, "Int64"), ("types", types, "UInt8"),
    ):
        array = ET.SubElement(cells_node, "DataArray", type=value_type, Name=name, format="ascii")
        array.text = _numbers(values)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def write_legacy_vtk(path: Path, bundle: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cells, types, metadata = _vtk_cells(bundle)
    lines = [
        "# vtk DataFile Version 3.0", "Vela Sentaurus Gmsh M1", "ASCII",
        "DATASET UNSTRUCTURED_GRID", f"POINTS {len(bundle['nodes'])} double",
    ]
    lines.extend(
        f"{node['x']:.17g} {node['y']:.17g} {node['z']:.17g}"
        for node in bundle["nodes"])
    lines.append(f"CELLS {len(cells)} {sum(len(cell) + 1 for cell in cells)}")
    lines.extend(f"{len(cell)} {_numbers(cell)}" for cell in cells)
    lines.append(f"CELL_TYPES {len(types)}")
    lines.extend(str(value) for value in types)
    lines.append(f"POINT_DATA {len(bundle['nodes'])}")
    for name in POINT_FIELD_NAMES:
        lines.extend([f"SCALARS {name} double 1", "LOOKUP_TABLE default"])
        lines.extend(f"{float(value):.17g}" for value in bundle["fields"][name]["values"])
    for name, field in sorted(bundle["fields"]["species"].items()):
        lines.extend([f"SCALARS {name} double 1", "LOOKUP_TABLE default"])
        lines.extend(f"{float(value):.17g}" for value in field["values"])
    lines.append(f"CELL_DATA {len(cells)}")
    for name in ("RegionId", "MaterialId", "PhysicalSurfaceTag", "ElementId",
                 "ContactId", "PhysicalCurveTag"):
        lines.extend([f"SCALARS {name} int 1", "LOOKUP_TABLE default"])
        lines.extend(str(item[name]) for item in metadata)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _canonical_edges(bundle: dict[str, Any]) -> list[tuple[int, int]]:
    return sorted({tuple(sorted(edge)) for triangle in bundle["triangles"] for edge in (
        (triangle["node_ids"][0], triangle["node_ids"][1]),
        (triangle["node_ids"][1], triangle["node_ids"][2]),
        (triangle["node_ids"][2], triangle["node_ids"][0]),
    )})


def _edge_reference(edge: tuple[int, int], edge_ids: dict[tuple[int, int], int]) -> int:
    canonical = tuple(sorted(edge))
    index = edge_ids[canonical]
    return index if edge == canonical else -(index + 1)


def write_dfise(path_grd: Path, path_dat: Path, bundle: dict[str, Any]) -> None:
    path_grd.parent.mkdir(parents=True, exist_ok=True)
    edges = _canonical_edges(bundle)
    edge_ids = {edge: index for index, edge in enumerate(edges)}
    incidence = Counter(tuple(sorted(edge)) for triangle in bundle["triangles"] for edge in (
        (triangle["node_ids"][0], triangle["node_ids"][1]),
        (triangle["node_ids"][1], triangle["node_ids"][2]),
        (triangle["node_ids"][2], triangle["node_ids"][0]),
    ))
    locations = "".join("e" if incidence[edge] == 1 else "i" for edge in edges)
    triangle_rows = []
    for triangle in bundle["triangles"]:
        a, b, c = triangle["node_ids"]
        references = [_edge_reference(edge, edge_ids) for edge in ((a, b), (b, c), (c, a))]
        triangle_rows.append(" 2 " + " ".join(map(str, references)))
    contact_elements: dict[str, list[int]] = defaultdict(list)
    contact_name = {item["physical_tag"]: item["name"] for item in bundle["contacts"]}
    contact_rows = []
    for line in bundle["lines"]:
        tag = line["physical_curve_tags"][0]
        if tag not in contact_name:
            continue
        element_id = len(triangle_rows) + len(contact_rows)
        contact_elements[contact_name[tag]].append(element_id)
        contact_rows.append(f" 1 {line['node_ids'][0]} {line['node_ids'][1]}")
    device_region_names = [item["name"] for item in bundle["regions"]]
    contact_region_names = [item["name"] for item in bundle["contacts"]]
    all_regions = device_region_names + contact_region_names
    material_name = {item["id"]: item["name"] for item in bundle["materials"]}
    region_material: dict[int, str] = {}
    region_cells: dict[int, list[int]] = defaultdict(list)
    for triangle in bundle["triangles"]:
        region_cells[triangle["region_id"]].append(triangle["dfise_id"])
        region_material[triangle["region_id"]] = material_name[triangle["material_id"]]
    materials = [region_material[item["id"]] for item in bundle["regions"]] + ["Contact"] * len(contact_region_names)
    vertices = "\n".join(f" {item['x']:.15e} {item['y']:.15e}" for item in bundle["nodes"])
    edge_rows = "\n".join(f" {a} {b}" for a, b in edges)
    element_rows = "\n".join(triangle_rows + contact_rows)
    region_blocks = []
    for region in bundle["regions"]:
        values = " ".join(map(str, region_cells[region["id"]]))
        region_blocks.append(
            f'  Region ("{region["name"]}") {{\n'
            f'    material = {region_material[region["id"]]}\n'
            f'    Elements ({len(region_cells[region["id"]])}) {{\n {values}\n    }}\n  }}')
    for contact in bundle["contacts"]:
        values = " ".join(map(str, contact_elements[contact["name"]]))
        region_blocks.append(
            f'  Region ("{contact["name"]}") {{\n'
            f'    material = Contact\n'
            f'    Elements ({len(contact_elements[contact["name"]])}) {{\n {values}\n    }}\n  }}')
    quoted_regions = " ".join(f'"{name}"' for name in all_regions)
    grd = f"""DF-ISE text

Info {{
  version = 1.1
  type = grid
  dimension = 2
  nb_vertices = {len(bundle['nodes'])}
  nb_edges = {len(edges)}
  nb_faces = 0
  nb_elements = {len(triangle_rows) + len(contact_rows)}
  nb_regions = {len(all_regions)}
  regions = [ {quoted_regions} ]
  materials = [ {' '.join(materials)} ]
}}

Data {{
  CoordSystem {{
    translate = [ 0 0 0 ]
    transform = [ 1 0 0 0 1 0 0 0 1 ]
  }}
  Vertices ({len(bundle['nodes'])}) {{
{vertices}
  }}
  Edges ({len(edges)}) {{
{edge_rows}
  }}
  Locations ({len(edges)}) {{
{locations}
  }}
  Elements ({len(triangle_rows) + len(contact_rows)}) {{
{element_rows}
  }}
{chr(10).join(region_blocks)}
}}
"""
    path_grd.write_text(grd, encoding="utf-8")

    datasets = {name: bundle["fields"][name]["values"] for name in POINT_FIELD_NAMES}
    datasets.update({name: field["values"] for name, field in bundle["fields"]["species"].items()})
    blocks = []
    validity = " ".join(f'"{name}"' for name in device_region_names)
    for name, values in datasets.items():
        value_rows = "\n".join(f" {float(value):.15e}" for value in values)
        blocks.append(
            f'  Dataset ("{name}") {{\n    function = {name}\n    type = scalar\n'
            f'    dimension = 1\n    location = vertex\n    validity = [ {validity} ]\n'
            f'    Values ({len(values)}) {{\n{value_rows}\n    }}\n  }}')
    names = " ".join(f'"{name}"' for name in datasets)
    dat = f"""DF-ISE text

Info {{
  version = 1.0
  type = dataset
  dimension = 2
  nb_vertices = {len(bundle['nodes'])}
  nb_edges = {len(edges)}
  nb_faces = 0
  nb_elements = {len(triangle_rows) + len(contact_rows)}
  nb_regions = {len(all_regions)}
  datasets = [ {names} ]
  functions = [ {' '.join(datasets)} ]
}}

Data {{
{chr(10).join(blocks)}
}}
"""
    path_dat.write_text(dat, encoding="utf-8")


def read_vtu(path: Path) -> dict[str, Any]:
    root = ET.parse(path).getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("VTU has no Piece")
    point_fields = {
        item.attrib["Name"]: [float(value) for value in (item.text or "").split()]
        for item in piece.findall("./PointData/DataArray")
    }
    cell_fields = {
        item.attrib["Name"]: [int(value) for value in (item.text or "").split()]
        for item in piece.findall("./CellData/DataArray")
    }
    cell_arrays = {item.attrib["Name"]: (item.text or "").split()
                   for item in piece.findall("./Cells/DataArray")}
    connectivity = [int(value) for value in cell_arrays["connectivity"]]
    offsets = [int(value) for value in cell_arrays["offsets"]]
    cells: list[list[int]] = []
    begin = 0
    for end in offsets:
        cells.append(connectivity[begin:end])
        begin = end
    point_array = piece.find("./Points/DataArray")
    if point_array is None:
        raise ValueError("VTU has no point coordinates")
    coordinate_values = [float(value) for value in (point_array.text or "").split()]
    if len(coordinate_values) % 3:
        raise ValueError("VTU point coordinates are not 3-component tuples")
    points = [coordinate_values[index:index + 3]
              for index in range(0, len(coordinate_values), 3)]
    return {
        "point_count": int(piece.attrib["NumberOfPoints"]),
        "cell_count": int(piece.attrib["NumberOfCells"]),
        "point_fields": point_fields, "cell_fields": cell_fields,
        "cell_types": [int(value) for value in cell_arrays["types"]],
        "points": points, "cells": cells,
    }


def read_legacy_vtk(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points_index = next(index for index, line in enumerate(lines) if line.startswith("POINTS "))
    point_count = int(lines[points_index].split()[1])
    points = [[float(value) for value in line.split()]
              for line in lines[points_index + 1:points_index + 1 + point_count]]
    cells_index = next(index for index, line in enumerate(lines) if line.startswith("CELLS "))
    cell_count = int(lines[cells_index].split()[1])
    cells = []
    for line in lines[cells_index + 1:cells_index + 1 + cell_count]:
        values = [int(value) for value in line.split()]
        if not values or values[0] != len(values) - 1:
            raise ValueError("invalid legacy VTK cell row")
        cells.append(values[1:])
    types_index = next(index for index, line in enumerate(lines)
                       if line.startswith("CELL_TYPES "))
    type_count = int(lines[types_index].split()[1])
    cell_types = [int(value) for value in lines[types_index + 1:types_index + 1 + type_count]]
    point_index = next(index for index, line in enumerate(lines) if line.startswith("POINT_DATA "))
    count = int(lines[point_index].split()[1])
    point_fields: dict[str, list[float]] = {}
    index = point_index + 1
    while index < len(lines) and not lines[index].startswith("CELL_DATA "):
        if not lines[index].startswith("SCALARS "):
            index += 1
            continue
        name = lines[index].split()[1]
        if lines[index + 1] != "LOOKUP_TABLE default":
            raise ValueError("invalid VTK scalar lookup table")
        point_fields[name] = [float(value) for value in lines[index + 2:index + 2 + count]]
        index += count + 2
    if index >= len(lines):
        raise ValueError("legacy VTK has no CELL_DATA")
    declared_cells = int(lines[index].split()[1])
    cell_fields: dict[str, list[int]] = {}
    index += 1
    while index < len(lines):
        if not lines[index].startswith("SCALARS "):
            index += 1
            continue
        name = lines[index].split()[1]
        if index + 1 >= len(lines) or lines[index + 1] != "LOOKUP_TABLE default":
            raise ValueError("invalid VTK cell scalar lookup table")
        cell_fields[name] = [
            int(value) for value in lines[index + 2:index + 2 + declared_cells]
        ]
        index += declared_cells + 2
    return {
        "point_count": point_count, "cell_count": cell_count,
        "points": points, "cells": cells, "cell_types": cell_types,
        "point_fields": point_fields, "cell_fields": cell_fields,
    }


def read_legacy_vtk_point_fields(path: Path) -> dict[str, list[float]]:
    """Compatibility wrapper for callers that only need nodal fields."""
    return read_legacy_vtk(path)["point_fields"]


def read_doping_csv(path: Path) -> dict[str, list[float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {
        "DonorConcentration": [float(row["donors_cm3"]) for row in rows],
        "AcceptorConcentration": [float(row["acceptors_cm3"]) for row in rows],
        "DopingConcentration": [float(row["net_doping_cm3"]) for row in rows],
    }


def compare_nodal_fields(reference: dict[str, Any], candidates: dict[str, dict[str, list[float]]],
                         tolerance: float = 1.0e-12) -> dict[str, Any]:
    failures = []
    maxima: dict[str, dict[str, float]] = {}
    for source, fields in candidates.items():
        maxima[source] = {}
        for name in POINT_FIELD_NAMES:
            if name not in fields:
                failures.append(f"{source} missing {name}")
                continue
            expected = reference["fields"][name]["values"]
            actual = fields[name]
            if len(actual) != len(expected):
                failures.append(f"{source} {name} length mismatch")
                continue
            maximum = max((abs(a - b) for a, b in zip(actual, expected)), default=0.0)
            maxima[source][name] = maximum
            scale = max((abs(value) for value in expected), default=1.0)
            if maximum > tolerance * max(scale, 1.0):
                failures.append(f"{source} {name} differs by {maximum}")
    return {"status": "passed" if not failures else "failed",
            "max_abs_difference": maxima, "failures": failures}
