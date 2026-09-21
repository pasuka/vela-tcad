#!/usr/bin/env python3
"""Prepare and audit the exact-topology Templates/LDMOS SProcess export.

This script intentionally stops at the structural contract.  It does not add
materials, transport physics, avalanche, thermal, or continuation settings.
The proprietary TDR and neutral exports belong under ignored
``reference_staging/``; only this reproducible procedure is tracked.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_IMPORTER = REPO.parent.parent / "build-release" / "sentaurus_import.exe"
DEFAULT_CONVERTER = REPO / "scripts" / "convert_tcad_export.py"
SEMICONDUCTOR_MATERIALS = {"Silicon", "Germanium", "GaAs"}
POLY_MATERIALS = {"PolySilicon"}
DIELECTRIC_MATERIALS = {"Oxide", "Nitride", "SiO2"}


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def run_checked(command: list[str], cwd: Path | None = None) -> dict[str, Any]:
    started = time.perf_counter()
    result = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return {
        "command": command,
        "return_code": result.returncode,
        "wall_clock_seconds": elapsed,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def import_tdr(importer: Path,
               tdr: Path,
               output: Path,
               policy: str,
               coordinate_unit: str) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=False)
    return run_checked([
        str(importer),
        "--tdr", str(tdr),
        "--inventory-json", str(output / "inventory.json"),
        "--field-values-json", str(output / "field_values.json"),
        "--export-dir", str(output),
        "--coordinate-unit", coordinate_unit,
        "--compensated-doping-policy", policy,
    ])


def triangle_area(points: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = points
    return abs((x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)) * 0.5


def triangle_angles(points: list[tuple[float, float]]) -> tuple[float, float, float]:
    values: list[float] = []
    for i in range(3):
        origin = points[i]
        a = (points[(i + 1) % 3][0] - origin[0], points[(i + 1) % 3][1] - origin[1])
        b = (points[(i + 2) % 3][0] - origin[0], points[(i + 2) % 3][1] - origin[1])
        denominator = math.hypot(*a) * math.hypot(*b)
        if denominator == 0.0:
            values.append(float("nan"))
            continue
        cosine = max(-1.0, min(1.0, (a[0] * b[0] + a[1] * b[1]) / denominator))
        values.append(math.degrees(math.acos(cosine)))
    return values[0], values[1], values[2]


def percentile(values: list[float], probability: float) -> float | None:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return None
    index = probability * (len(finite) - 1)
    lo = math.floor(index)
    hi = math.ceil(index)
    if lo == hi:
        return finite[lo]
    return finite[lo] * (hi - index) + finite[hi] * (index - lo)


def material_class(material: str) -> str:
    if material in SEMICONDUCTOR_MATERIALS:
        return "transport_semiconductor"
    if material in POLY_MATERIALS:
        return "electrostatic_only_semiconductor_poly"
    if material in DIELECTRIC_MATERIALS:
        return "dielectric"
    if material in {"Gas", "Air"}:
        return "ignored_gas"
    return "unclassified"


def edge_key(a: int, b: int) -> tuple[int, int]:
    return (a, b) if a < b else (b, a)


def region_nodes(region: dict[str, Any]) -> list[int]:
    nodes: set[int] = set()
    for tri in region.get("triangles", []):
        nodes.update(int(value) for value in tri)
    for edge in region.get("edges", []):
        nodes.update(int(value) for value in edge)
    return sorted(nodes)


def supported_field_rows(inventory: dict[str, Any], field: dict[str, Any]) -> list[int]:
    region = next(
        item for item in inventory["geometry"]["regions"]
        if int(item["index"]) == int(field["region"])
    )
    if int(field.get("location_type", 0)) == 3 and len(field.get("raw_values", [])) == len(region.get("triangles", [])):
        return []
    count = len(field.get("raw_values", [])) // max(int(field.get("components", 1)), 1)
    if int(region["type"]) != 1 and count == len(inventory["geometry"]["vertices"]):
        return list(range(count))
    return region_nodes(region)


def audit_structure(neutral: Path, alternate: Path) -> dict[str, Any]:
    inventory = json.loads((neutral / "field_values.json").read_text(encoding="utf-8"))
    metadata = json.loads((neutral / "metadata.json").read_text(encoding="utf-8"))
    coordinate_scale = float(metadata.get("coordinate_to_um_scale", 1.0))
    raw_vertices = [tuple(map(float, point)) for point in inventory["geometry"]["vertices"]]
    vertices = [
        (point[0] * coordinate_scale, point[1] * coordinate_scale)
        for point in raw_vertices
    ]
    material_regions = [
        region for region in inventory["geometry"]["regions"] if int(region["type"]) == 0
    ]
    contact_regions = [
        region for region in inventory["geometry"]["regions"] if int(region["type"]) == 1
    ]

    nodes_rows = read_csv(neutral / "nodes.csv")
    exported_nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"])) for row in nodes_rows
    }
    coordinate_errors = [
        max(abs(exported_nodes[index][axis] - point[axis]) for axis in (0, 1))
        for index, point in enumerate(vertices)
    ]

    element_rows = read_csv(neutral / "elements.csv")
    exported_triangles: dict[str, set[tuple[int, int, int]]] = defaultdict(set)
    exported_areas: dict[str, float] = defaultdict(float)
    for row in element_rows:
        tri = (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        exported_triangles[row["region"]].add(tuple(sorted(tri)))
        exported_areas[row["region"]] += triangle_area([exported_nodes[node] for node in tri])
    source_triangles = {
        region["name"]: {tuple(sorted(map(int, tri))) for tri in region["triangles"]}
        for region in material_regions
    }

    areas: dict[str, float] = {}
    geometry_rows: list[dict[str, Any]] = []
    edge_opposites: dict[tuple[int, int], list[float]] = defaultdict(list)
    region_node_areas: dict[str, dict[int, float]] = {}
    degenerate = 0
    obtuse = 0
    negative_circumcentric_weights: list[dict[str, Any]] = []
    all_angles: list[float] = []
    edge_to_regions: dict[tuple[int, int], set[str]] = defaultdict(set)
    for region in material_regions:
        area = 0.0
        local_node_area: dict[int, float] = defaultdict(float)
        for local_cell, tri in enumerate(region["triangles"]):
            tri_ids = [int(value) for value in tri]
            points = [vertices[node] for node in tri_ids]
            cell_area = triangle_area(points)
            angles = triangle_angles(points)
            area += cell_area
            if cell_area <= 0.0 or not all(math.isfinite(value) for value in angles):
                degenerate += 1
            else:
                all_angles.extend(angles)
                if max(angles) > 90.0 + 1.0e-10:
                    obtuse += 1
                for node in tri_ids:
                    local_node_area[node] += cell_area / 3.0
                for local, pair in enumerate(((1, 2), (2, 0), (0, 1))):
                    key = edge_key(tri_ids[pair[0]], tri_ids[pair[1]])
                    edge_opposites[key].append(angles[local])
                    edge_to_regions[key].add(region["name"])
                    cotangent_half = 0.5 / math.tan(math.radians(angles[local]))
                    if cotangent_half < -1.0e-14:
                        negative_circumcentric_weights.append({
                            "region": region["name"],
                            "region_local_cell": local_cell,
                            "opposite_node": tri_ids[local],
                            "edge": list(key),
                            "opposite_angle_deg": angles[local],
                            "half_cotangent_weight": cotangent_half,
                        })
        areas[region["name"]] = area
        region_node_areas[region["name"]] = dict(local_node_area)
        exported_area = exported_areas[region["name"]]
        relative_area_error = abs(exported_area - area) / max(abs(area), 1.0e-300)
        geometry_rows.append({
            "name": region["name"],
            "material": region["material"],
            "classification": material_class(region["material"]),
            "triangle_count": len(region["triangles"]),
            "node_occurrence_count": len(region_nodes(region)),
            "area_um2": area,
            "exported_area_um2": exported_area,
            "relative_area_error": relative_area_error,
        })

    non_delaunay = [
        {"edge": list(edge), "opposite_angle_sum_deg": sum(opposites)}
        for edge, opposites in edge_opposites.items()
        if len(opposites) == 2 and sum(opposites) > 180.0 + 1.0e-10
    ]

    contact_rows = {row["name"]: row for row in read_csv(neutral / "contacts.csv")}
    material_by_region = {region["name"]: region["material"] for region in material_regions}
    contacts: list[dict[str, Any]] = []
    for contact in contact_regions:
        source_edges = {edge_key(int(edge[0]), int(edge[1])) for edge in contact["edges"]}
        exported_contact_nodes = {
            int(value) for value in contact_rows[contact["name"]]["node_ids"].split(";") if value
        }
        serialized_edges = contact_rows[contact["name"]].get("edge_node_ids", "")
        if serialized_edges:
            candidate_edges = {
                edge_key(*(int(value) for value in item.split("-")))
                for item in serialized_edges.split(";") if item
            }
            edge_mapping = "explicit_edge_node_ids"
        else:
            candidate_edges = {
                edge for edge in edge_to_regions
                if edge[0] in exported_contact_nodes and edge[1] in exported_contact_nodes
            }
            edge_mapping = "legacy_node_set_inference"
        adjacent_regions = sorted({name for edge in source_edges for name in edge_to_regions.get(edge, set())})
        contacts.append({
            "name": contact["name"],
            "electrical": contact["name"] in {"gate", "source", "drain", "substrate"},
            "thermal_only": contact["name"] == "th_lat",
            "owner_region": contact_rows[contact["name"]]["region"],
            "adjacent_regions": adjacent_regions,
            "adjacent_materials": sorted({material_by_region[name] for name in adjacent_regions}),
            "edge_count": len(source_edges),
            "edge_mapping": edge_mapping,
            "length_um": sum(math.dist(vertices[a], vertices[b]) for a, b in source_edges),
            "edge_set_exact": source_edges == candidate_edges,
            "missing_edges": [list(edge) for edge in sorted(source_edges - candidate_edges)],
            "extra_edges": [list(edge) for edge in sorted(candidate_edges - source_edges)],
        })

    doping_rows = {
        int(row["node_id"]): (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
        for row in read_csv(neutral / "doping.csv")
    }
    alternate_rows = {
        int(row["node_id"]): (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
        for row in read_csv(alternate / "doping.csv")
    }
    invalid_doping = [
        node for node, values in doping_rows.items()
        if any(not math.isfinite(value) or value < 0.0 for value in values)
    ]
    changed_policy_nodes = [node for node in doping_rows if doping_rows[node] != alternate_rows[node]]
    region_doping_integrals = []
    for region in material_regions:
        node_areas = region_node_areas[region["name"]]
        region_doping_integrals.append({
            "region": region["name"],
            "material": region["material"],
            "donor_integral_per_cm_out_of_plane": sum(
                doping_rows[node][0] * area * 1.0e-8 for node, area in node_areas.items()
            ),
            "acceptor_integral_per_cm_out_of_plane": sum(
                doping_rows[node][1] * area * 1.0e-8 for node, area in node_areas.items()
            ),
        })
    donor_integral = sum(item["donor_integral_per_cm_out_of_plane"] for item in region_doping_integrals)
    acceptor_integral = sum(item["acceptor_integral_per_cm_out_of_plane"] for item in region_doping_integrals)

    dex_errors: list[float] = []
    sign_mismatches = 0
    occurrence_count = 0
    material_by_index = {
        int(region["index"]): region.get("material", "") for region in material_regions
    }
    for field in inventory["fields"]:
        if field["name"] not in {"DopingConcentration", "NetActive"} or "raw_values" not in field:
            continue
        if material_by_index.get(int(field["region"]), "") not in (
            SEMICONDUCTOR_MATERIALS | POLY_MATERIALS
        ):
            continue
        components = max(int(field.get("components", 1)), 1)
        nodes = supported_field_rows(inventory, field)
        values = field["raw_values"]
        for row, node in enumerate(nodes[: len(values) // components]):
            reference = float(values[row * components])
            exported = doping_rows[node][0] - doping_rows[node][1]
            occurrence_count += 1
            if reference * exported < 0.0:
                sign_mismatches += 1
            floor = 1.0
            dex_errors.append(abs(
                math.log10(max(abs(reference), floor)) -
                math.log10(max(abs(exported), floor))
            ))

    source_bounds = {
        "x_min_um": min(point[0] for point in vertices),
        "x_max_um": max(point[0] for point in vertices),
        "y_min_um": min(point[1] for point in vertices),
        "y_max_um": max(point[1] for point in vertices),
    }
    triangle_exact = source_triangles == dict(exported_triangles)
    max_area_relative_error = max(
        (region["relative_area_error"] for region in geometry_rows), default=0.0
    )
    return {
        "schema": "vela.templates_ldmos.structure_audit.v1-draft",
        "coordinate_contract": {
            "source_unit": metadata.get("source_coordinate_unit", "um"),
            "exported_unit": metadata.get("exported_coordinate_unit", "um"),
            "coordinate_to_um_scale": coordinate_scale,
            "evidence": "explicit sentaurus_import --coordinate-unit contract and SProcess deck dimensions",
            "bounds": source_bounds,
            "source_vertex_count": len(vertices),
            "exported_node_count": len(exported_nodes),
            "mapping": "global_vertex_id_identity; region-side occurrence order is ascending global vertex id",
            "max_abs_coordinate_error_um": max(coordinate_errors, default=0.0),
        },
        "topology": {
            "source_triangle_count": sum(len(region["triangles"]) for region in material_regions),
            "exported_triangle_count": len(element_rows),
            "triangle_sets_exact_by_region": triangle_exact,
            "max_region_area_relative_error": max_area_relative_error,
            "regions": geometry_rows,
        },
        "contacts": contacts,
        "doping": {
            "unit": "cm^-3",
            "invalid_node_count": len(invalid_doping),
            "reported_vs_dominant_signed_region_changed_nodes": len(changed_policy_nodes),
            "changed_node_ids": changed_policy_nodes,
            "donor_integral_per_cm_out_of_plane": donor_integral,
            "acceptor_integral_per_cm_out_of_plane": acceptor_integral,
            "region_integrals": region_doping_integrals,
            "netactive_same_support_occurrences": occurrence_count,
            "netactive_p95_dex": percentile(dex_errors, 0.95),
            "netactive_max_dex": max(dex_errors, default=0.0),
            "netactive_sign_mismatches": sign_mismatches,
            "note": "Global Vela doping is compared to each Sentaurus region-side DopingConcentration occurrence; interface conflicts remain explicit.",
        },
        "mesh_quality": {
            "degenerate_triangle_count": degenerate,
            "obtuse_triangle_count": obtuse,
            "minimum_angle_deg": min(all_angles, default=None),
            "maximum_angle_deg": max(all_angles, default=None),
            "non_delaunay_interior_edge_count": len(non_delaunay),
            "non_delaunay_interior_edges": non_delaunay,
            "negative_raw_circumcentric_weight_count": len(negative_circumcentric_weights),
            "negative_raw_circumcentric_weights": negative_circumcentric_weights,
            "negative_box_weight_note": (
                "These are raw half-cotangent geometric weights, not the frozen barycentric "
                "candidate control volumes; profile-specific truncation remains explicit."
            ),
        },
        "gates": {
            "coordinate_max_abs_le_1e_10_um": max(coordinate_errors, default=0.0) <= 1.0e-10,
            "triangles_exact": triangle_exact,
            "region_area_relative_error_le_1e_10": max_area_relative_error <= 1.0e-10,
            "contacts_exact": all(item["edge_set_exact"] for item in contacts),
            "doping_finite_nonnegative": not invalid_doping,
            "netactive_p95_le_1e_4_dex": (percentile(dex_errors, 0.95) or 0.0) <= 1.0e-4,
            "mesh_anomalies_explicit": True,
        },
    }


def render_markdown(audit: dict[str, Any], repeatability: dict[str, Any]) -> str:
    coordinate = audit["coordinate_contract"]
    topology = audit["topology"]
    doping = audit["doping"]
    mesh = audit["mesh_quality"]
    p95_text = "not_available" if doping["netactive_p95_dex"] is None else f"{doping['netactive_p95_dex']:.3e}"
    lines = [
        "# Templates/LDMOS stage 1 structural audit",
        "",
        "This report covers exact-topology import only. It is not a physics-validation result.",
        "",
        "## Summary",
        "",
        f"- Vertices / triangles: {coordinate['source_vertex_count']} / {topology['source_triangle_count']}",
        f"- Coordinate max error: {coordinate['max_abs_coordinate_error_um']:.3e} um",
        f"- Exact triangle sets: {topology['triangle_sets_exact_by_region']}",
        f"- Exact contact edge sets: {all(item['edge_set_exact'] for item in audit['contacts'])}",
        f"- Obtuse / non-Delaunay / degenerate: {mesh['obtuse_triangle_count']} / {mesh['non_delaunay_interior_edge_count']} / {mesh['degenerate_triangle_count']}",
        f"- NetActive P95 / max: {p95_text} / {doping['netactive_max_dex']:.3e} dex",
        f"- Doping policy changed nodes: {doping['reported_vs_dominant_signed_region_changed_nodes']}",
        f"- Byte-for-byte repeatable neutral export: {repeatability['identical']}",
        "",
        "## Regions",
        "",
        "| Region | Material | Classification | Triangles | Area (um^2) |",
        "| --- | --- | --- | ---: | ---: |",
    ]
    for region in topology["regions"]:
        lines.append(
            f"| {region['name']} | {region['material']} | {region['classification']} | "
            f"{region['triangle_count']} | {region['area_um2']:.12g} |"
        )
    lines.extend([
        "",
        "## Contacts",
        "",
        "`th_lat`, when present, is thermal-only and must not become an electrical electrode.",
        "",
        "| Contact | Owner | Adjacent regions | Edges | Length (um) | Exact |",
        "| --- | --- | --- | ---: | ---: | --- |",
    ])
    for contact in audit["contacts"]:
        lines.append(
            f"| {contact['name']} | {contact['owner_region']} | "
            f"{', '.join(contact['adjacent_regions'])} | {contact['edge_count']} | "
            f"{contact['length_um']:.12g} | {contact['edge_set_exact']} |"
        )
    lines.extend([
        "",
        "## Qualification boundary",
        "",
        "The geometry audit does not authorize the PN2D-specific `element_edge_sg_gss_laux` bundle. "
        "The discretization draft remains fail-closed until WP1.75 supplies its schema and the exact "
        "mesh policy is independently approved.",
        "",
    ])
    return "\n".join(lines)


def discretization_draft(audit: dict[str, Any]) -> dict[str, Any]:
    has_obtuse = audit["mesh_quality"]["obtuse_triangle_count"] > 0
    return {
        "schema": "vela.templates_ldmos.discretization_contract.v1-draft-unvalidated",
        "status": "draft_pending_wp1_75_schema_and_independent_approval",
        "applicable_mesh": "exact T-2022.03-SP2 Templates/LDMOS SProcess topology",
        "profile_name": "templates_ldmos_exact_topology_legacy_candidate",
        "current_support": "scharfetter_gummel_edge_flux",
        "control_volume": "barycentric",
        "field_recovery": "cell_reconstructed",
        "volume_source_mapping": "cell_reconstructed",
        "contact_edge_integration": "exact_imported_contact_boundary_edges",
        "obtuse_policy": "accepted_and_listed_with_barycentric_positive_control_volumes" if has_obtuse else "not_present",
        "non_delaunay_policy": "accepted_only_for_structural_import; physics parity unresolved",
        "require_non_obtuse": False,
        "avalanche_profile": "legacy_cell_reconstructed",
        "forbidden_inference": "Do not apply PN2D element_edge_sg_gss_laux as a global or device default.",
        "physics_use_authorized": False,
    }


def separate_thermal_contacts(converted: Path, audit: dict[str, Any]) -> None:
    thermal = [item for item in audit["contacts"] if item["thermal_only"]]
    write_json(converted / "thermal_contacts.json", {
        "schema": "vela.templates_ldmos.thermal_contacts.v1-draft",
        "contacts": thermal,
        "electrical_use_authorized": False,
    })
    thermal_names = {item["name"] for item in thermal}
    contact_audit = {item["name"]: item for item in audit["contacts"]}
    for path in sorted(converted.glob("simulation_*.json")):
        deck = json.loads(path.read_text(encoding="utf-8"))
        deck["contacts"] = [
            contact for contact in deck.get("contacts", [])
            if contact.get("name") not in thermal_names
        ]
        for contact in deck["contacts"]:
            audit_entry = contact_audit.get(contact.get("name"), {})
            owner = str(audit_entry.get("owner_region", "")).lower()
            if contact.get("name") == "gate" and "oxide" in owner:
                # The imported gate boundary is on the oxide and must not be
                # interpreted as an Ohmic carrier contact.  Zero flat-band is
                # only the stage-1.5 solver-qualification reference; the
                # Sentaurus PolySi work-function mapping is a stage-2 physics
                # control and remains deliberately unresolved here.
                contact["type"] = "metal_gate"
                contact["flatband_voltage"] = 0.0
        deck["_comment"] = (
            str(deck.get("_comment", "")) + " Thermal-only contacts are retained in "
            "mesh.json/thermal_contacts.json but deliberately omitted from electrical biases. "
            "The oxide-owned gate is a metal_gate with a provisional 0 V flat-band value for "
            "stage-1.5 solver qualification only; stage 2 owns the PolySi work-function mapping."
        ).strip()
        write_json(path, deck)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tdr", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--python", default="python")
    parser.add_argument("--converter", type=Path, default=DEFAULT_CONVERTER)
    parser.add_argument(
        "--source-coordinate-unit", choices=["cm", "um"], default="cm",
        help="Explicit raw TDR geometry unit; Templates/LDMOS SProcess uses cm.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    tdr = args.tdr.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite stage-1 directory: {output}")
    if not tdr.is_file():
        raise FileNotFoundError(tdr)
    if not args.importer.is_file():
        raise FileNotFoundError(args.importer)
    output.mkdir(parents=True)
    primary = output / "neutral_reported"
    alternate = output / "neutral_dominant_signed_region"
    repeated = output / "neutral_reported_repeat"
    commands = [
        import_tdr(args.importer, tdr, primary, "reported", args.source_coordinate_unit),
        import_tdr(args.importer, tdr, alternate, "dominant_signed_region", args.source_coordinate_unit),
        import_tdr(args.importer, tdr, repeated, "reported", args.source_coordinate_unit),
    ]
    converted = output / "vela_exact_topology"
    commands.append(run_checked([
        args.python,
        str(args.converter),
        "--input-dir", str(primary),
        "--output-dir", str(converted),
        "--device", "ldmos2d",
        "--simulation-types", "iv,bv",
    ]))

    primary_hashes = tree_hashes(primary)
    repeated_hashes = tree_hashes(repeated)
    repeatability = {
        "identical": primary_hashes == repeated_hashes,
        "primary_hashes": primary_hashes,
        "repeat_hashes": repeated_hashes,
        "differing_files": sorted(
            name for name in set(primary_hashes) | set(repeated_hashes)
            if primary_hashes.get(name) != repeated_hashes.get(name)
        ),
    }
    audit = audit_structure(primary, alternate)
    separate_thermal_contacts(converted, audit)
    write_json(output / "reports" / "structure_audit.json", audit)
    write_json(output / "reports" / "repeatability.json", repeatability)
    write_json(output / "contracts" / "discretization_contract.json", discretization_draft(audit))
    write_json(output / "manifest.json", {
        "schema": "vela.templates_ldmos.stage1_manifest.v1-draft",
        "source_tdr": {"name": tdr.name, "sha256": sha256(tdr)},
        "importer": {"path": str(args.importer), "sha256": sha256(args.importer)},
        "converter": {"path": str(args.converter), "sha256": sha256(args.converter)},
        "commands": commands,
        "audit_gates": audit["gates"],
        "repeatable": repeatability["identical"],
    })
    (output / "reports" / "structure_audit.md").write_text(
        render_markdown(audit, repeatability), encoding="utf-8"
    )
    print(json.dumps({
        "output": str(output),
        "gates": audit["gates"],
        "repeatable": repeatability["identical"],
    }, indent=2))
    return 0 if all(audit["gates"].values()) and repeatability["identical"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
