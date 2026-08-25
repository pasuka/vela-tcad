#!/usr/bin/env python3
"""Versioned Gmsh OCC backend and normalized mesh/field bundle construction."""

from __future__ import annotations

import math
from collections import defaultdict, deque
from pathlib import Path
from typing import Any, Iterable

from gmsh_runtime import EXPECTED_GMSH_VERSION
from sentaurus_doping import DopingEvaluationError, evaluate_m2_fields


BUNDLE_SCHEMA = "vela.mesh_field_bundle.v1"
GEOM_TOL = 1.0e-9
REGION_TAG_BASE = 1000
MATERIAL_TAG_BASE = 2000
CONTACT_TAG_BASE = 3000
INTERFACE_TAG_BASE = 4000
LOG_DOPING_FLOOR_CM3 = 1.0


class GmshBackendError(ValueError):
    """A fail-closed capability, geometry, mesh, or field error."""

    def __init__(self, message: str, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.details = details or {}


def _rectangle(body: dict[str, Any]) -> tuple[float, float, float, float]:
    geometry = body["geometry"]
    x0, y0 = map(float, geometry["lower_left"])
    x1, y1 = map(float, geometry["upper_right"])
    if x1 <= x0 + GEOM_TOL or y1 <= y0 + GEOM_TOL:
        raise GmshBackendError(
            f"body {body['id']!r} is not a positive-area rectangle",
            {"stage": "geometry", "body": body["id"]},
        )
    return x0, y0, x1, y1


def validate_m1_device_ir(device_ir: dict[str, Any]) -> None:
    """Accept only the deliberately narrow M1 semantic subset."""
    if device_ir.get("schema") != "vela.sentaurus_device_ir.v2":
        raise GmshBackendError(
            "Gmsh M1 requires vela.sentaurus_device_ir.v2",
            {"stage": "input", "schema": device_ir.get("schema")},
        )
    bodies = device_ir["geometry"]["bodies"]
    if not bodies:
        raise GmshBackendError("Device IR contains no bodies", {"stage": "geometry"})
    if device_ir["geometry"]["inserted_vertices"]:
        raise GmshBackendError(
            "M1 does not support insert-vertex; use M2",
            {"stage": "capability", "capability": "insert-vertex"},
        )
    rectangles: list[tuple[dict[str, Any], tuple[float, float, float, float]]] = []
    for body in bodies:
        if body.get("shape") != "rectangle":
            raise GmshBackendError(
                f"M1 supports Rectangle bodies only, got {body.get('shape')!r}",
                {"stage": "capability", "body": body["id"], "shape": body.get("shape")},
            )
        rectangles.append((body, _rectangle(body)))
    for index, (first, a) in enumerate(rectangles):
        for second, b in rectangles[index + 1:]:
            overlap_x = min(a[2], b[2]) - max(a[0], b[0])
            overlap_y = min(a[3], b[3]) - max(a[1], b[1])
            if overlap_x > GEOM_TOL and overlap_y > GEOM_TOL:
                raise GmshBackendError(
                    "M1 refuses overlapping bodies; ABA/BAB ownership is M2",
                    {"stage": "capability", "capability": "boolean-ownership",
                     "bodies": [first["id"], second["id"]]},
                )
            if ((abs(overlap_x) <= GEOM_TOL and overlap_y > GEOM_TOL) or
                    (abs(overlap_y) <= GEOM_TOL and overlap_x > GEOM_TOL)):
                raise GmshBackendError(
                    "M1 requires strictly separated rectangles; conformal interfaces are M2",
                    {"stage": "capability", "capability": "conformal-interface",
                     "bodies": [first["id"], second["id"]]},
                )

    mesh_control = device_ir["mesh_control"]
    if mesh_control["refinement_functions"] or mesh_control["refinement_placements"]:
        raise GmshBackendError(
            "M1 does not support SDE refinement placements or functions; use M3",
            {"stage": "capability", "capability": "refinement"},
        )
    profiles = {item["name"]: item for item in device_ir["doping"]["profiles"]}
    for profile in profiles.values():
        if profile.get("type") != "constant":
            raise GmshBackendError(
                "M1 supports constant doping profiles only; Gaussian is M2",
                {"stage": "capability", "capability": "gaussian-doping",
                 "profile": profile.get("name")},
            )
    region_names = {body["region"] for body in bodies}
    material_names = {body["material"]["sde"] for body in bodies}
    for placement in device_ir["doping"]["placements"]:
        if placement["profile"] not in profiles:
            raise GmshBackendError(
                f"placement references unknown profile {placement['profile']!r}",
                {"stage": "doping", "placement": placement["name"]},
            )
        if placement.get("replacement") != "NoReplace":
            raise GmshBackendError(
                "M1 supports NoReplace constant placement only; ordered replacement is M2",
                {"stage": "capability", "capability": "doping-replacement",
                 "placement": placement["name"]},
            )
        kind, target = placement["target_kind"], placement["target"]
        if kind == "region" and target not in region_names:
            raise GmshBackendError(f"unknown placement region {target!r}")
        if kind == "material" and target not in material_names:
            raise GmshBackendError(f"unknown placement material {target!r}")
        if kind not in {"region", "material"}:
            raise GmshBackendError(
                "M1 supports region/material constant placement only",
                {"stage": "capability", "capability": "window-doping",
                 "placement": placement["name"]},
            )
    for contact in device_ir["contacts"]:
        if not contact["assignments"]:
            raise GmshBackendError(
                f"contact {contact['name']!r} has no assignment",
                {"stage": "contact", "contact": contact["name"]},
            )
        if any(item["kind"] != "pick_point" for item in contact["assignments"]):
            raise GmshBackendError(
                "M1 supports pick-point contact assignments only",
                {"stage": "capability", "capability": "contact-entity-list",
                 "contact": contact["name"]},
            )


def validate_m2_device_ir(device_ir: dict[str, Any]) -> None:
    """Validate the M2 rectangle/boolean/contact/doping semantic subset."""
    if device_ir.get("schema") != "vela.sentaurus_device_ir.v2":
        raise GmshBackendError(
            "Gmsh M2 requires vela.sentaurus_device_ir.v2",
            {"stage": "input", "schema": device_ir.get("schema")},
        )
    bodies = device_ir["geometry"]["bodies"]
    if not bodies:
        raise GmshBackendError("Device IR contains no bodies", {"stage": "geometry"})
    orders = [int(body["creation_order"]) for body in bodies]
    if sorted(orders) != list(range(len(bodies))):
        raise GmshBackendError(
            "M2 body creation_order must be unique and continuous",
            {"stage": "geometry", "creation_orders": orders},
        )
    for body in bodies:
        if body.get("shape") != "rectangle":
            raise GmshBackendError(
                "M2 supports Rectangle bodies only; Polygon is not implemented",
                {"stage": "capability", "capability": "polygon-body",
                 "body": body["id"], "shape": body.get("shape")},
            )
        _rectangle(body)
        if body.get("boolean_mode") not in {"ABA", "BAB"}:
            raise GmshBackendError(
                f"boolean mode {body.get('boolean_mode')!r} is unsupported in M2",
                {"stage": "capability", "capability": "boolean-mode",
                 "body": body["id"], "mode": body.get("boolean_mode")},
            )

    mesh_control = device_ir["mesh_control"]
    if mesh_control["refinement_functions"] or mesh_control["refinement_placements"]:
        raise GmshBackendError(
            "M2 does not support refinement placements or functions; use M3",
            {"stage": "capability", "capability": "refinement"},
        )
    windows = {item["name"]: item for item in device_ir["doping"]["windows"]}
    profiles = {item["name"]: item for item in device_ir["doping"]["profiles"]}
    for profile in profiles.values():
        if profile.get("type") not in {"constant", "gaussian"}:
            raise GmshBackendError(
                f"unsupported doping profile type {profile.get('type')!r}",
                {"stage": "capability", "capability": "doping-profile",
                 "profile": profile.get("name")},
            )
    region_names = {body["region"] for body in bodies}
    material_names = {body["material"]["sde"] for body in bodies}
    for placement in device_ir["doping"]["placements"]:
        profile = profiles.get(placement.get("profile"))
        if profile is None:
            raise GmshBackendError(
                f"placement references unknown profile {placement.get('profile')!r}",
                {"stage": "doping", "placement": placement.get("name")},
            )
        if placement.get("replacement", "NoReplace") not in {
            "Replace", "LocalReplace", "NoReplace"
        }:
            raise GmshBackendError(
                f"unsupported replacement mode {placement.get('replacement')!r}",
                {"stage": "doping", "placement": placement.get("name")},
            )
        if profile["type"] == "constant" and (
                abs(float(placement.get("decay_length", 0.0))) > GEOM_TOL or
                placement.get("decay", "none") != "none"):
            raise GmshBackendError(
                "M2 does not implement constant-profile decay; use M3",
                {"stage": "capability", "capability": "constant-profile-decay",
                 "placement": placement.get("name"),
                 "source_location": placement.get("provenance", {}).get(
                     "source_location")},
            )
        kind, target = placement.get("target_kind"), placement.get("target")
        valid = ((kind == "region" and target in region_names) or
                 (kind == "material" and target in material_names) or
                 (kind == "window" and target in windows))
        if not valid:
            raise GmshBackendError(
                f"invalid placement target {kind!r}:{target!r}",
                {"stage": "doping", "placement": placement.get("name")},
            )
        if profile["type"] == "gaussian" and (
                kind != "window" or windows[target].get("shape") != "line"):
            raise GmshBackendError(
                "M2 Gaussian placement requires a Line baseline window",
                {"stage": "doping", "placement": placement.get("name")},
            )
        if profile["type"] == "gaussian" and (
                abs(float(placement.get("restriction_decay_length", 0.0))) > GEOM_TOL or
                placement.get("restriction_decay", "none") != "none"):
            raise GmshBackendError(
                "M2 does not implement analytical restriction decay; use M3",
                {"stage": "capability", "capability": "analytical-restriction-decay",
                 "placement": placement.get("name"),
                 "source_location": placement.get("provenance", {}).get(
                     "source_location")},
            )
    for contact in device_ir["contacts"]:
        if not contact["assignments"]:
            raise GmshBackendError(
                f"contact {contact['name']!r} has no assignment",
                {"stage": "contact", "contact": contact["name"]},
            )
        if any(item["kind"] != "pick_point" for item in contact["assignments"]):
            raise GmshBackendError(
                "M2 supports pick-point contact assignments only",
                {"stage": "capability", "capability": "contact-entity-list",
                 "contact": contact["name"]},
            )


def _groups(names: Iterable[str], namespace: str, base: int) -> list[dict[str, Any]]:
    return [
        {"id": index, "name": name, "physical_name": f"{namespace}:{name}",
         "physical_tag": base + index + 1}
        for index, name in enumerate(sorted(set(names)))
    ]


def _curve_contains_pick(gmsh, curve: int, x: float, y: float) -> tuple[bool, bool]:
    xmin, ymin, _zmin, xmax, ymax, _zmax = gmsh.model.getBoundingBox(1, curve)
    scale = max(xmax - xmin, ymax - ymin, 1.0)
    # OpenCASCADE bounding boxes are padded by its modelling tolerance (about
    # 1e-7 in this build), which is intentionally larger than the authored SDE
    # coordinate comparison tolerance.
    tol = max(GEOM_TOL * scale, 5.0e-7 * scale)
    if not (xmin - tol <= x <= xmax + tol and ymin - tol <= y <= ymax + tol):
        return False, False
    vertical = xmax - xmin <= tol
    horizontal = ymax - ymin <= tol
    if not (vertical or horizontal):
        raise GmshBackendError(
            "M1 contact selection encountered a non-axis-aligned curve",
            {"stage": "contact", "curve": curve},
        )
    distance = abs(x - 0.5 * (xmin + xmax)) if vertical else abs(y - 0.5 * (ymin + ymax))
    if distance > tol:
        return False, False
    free, lo, hi = (y, ymin, ymax) if vertical else (x, xmin, xmax)
    endpoint = abs(free - lo) <= tol or abs(free - hi) <= tol
    return True, endpoint


def _boundary_curves(gmsh, surfaces: Iterable[int]) -> list[int]:
    curves = set()
    for surface in surfaces:
        for dimension, tag in gmsh.model.getBoundary(
                [(2, int(surface))], combined=False, oriented=False):
            if int(dimension) == 1:
                curves.add(int(tag))
    return sorted(curves)


def _surface_vertices(gmsh, surface: int) -> list[tuple[float, float]]:
    vertices = set()
    for curve in _boundary_curves(gmsh, [surface]):
        for dimension, tag in gmsh.model.getBoundary(
                [(1, curve)], combined=False, oriented=False):
            if int(dimension) != 0:
                continue
            values = gmsh.model.getValue(0, int(tag), [])
            vertices.add((float(values[0]), float(values[1])))
    if len(vertices) < 3:
        raise GmshBackendError(
            "OCC fragment has fewer than three boundary vertices",
            {"stage": "geometry", "surface": surface, "vertices": sorted(vertices)},
        )
    return sorted(vertices)


def _point_in_rectangle(bounds: tuple[float, float, float, float],
                        point: tuple[float, float]) -> bool:
    x0, y0, x1, y1 = bounds
    scale = max(x1 - x0, y1 - y0, 1.0)
    tolerance = 1.0e-7 * scale
    return (x0 - tolerance <= point[0] <= x1 + tolerance and
            y0 - tolerance <= point[1] <= y1 + tolerance)


def _fragment_m2_geometry(gmsh, bodies: list[dict[str, Any]],
                          inserted_vertices: list[dict[str, Any]]) -> tuple[
        dict[int, dict[str, Any]], dict[str, Any]]:
    """Fragment all source rectangles, audit the map, then replay ABA/BAB."""
    source_tags = []
    bounds_by_id = {}
    for body in bodies:
        x0, y0, x1, y1 = _rectangle(body)
        bounds_by_id[body["id"]] = (x0, y0, x1, y1)
        source_tags.append(gmsh.model.occ.addRectangle(x0, y0, 0.0, x1 - x0, y1 - y0))
    point_tags = []
    for record in inserted_vertices:
        x, y = map(float, record["position"])
        authored_segments = []
        for body in bodies:
            x0, y0, x1, y1 = bounds_by_id[body["id"]]
            if ((abs(x - x0) <= GEOM_TOL or abs(x - x1) <= GEOM_TOL) and
                    y0 + GEOM_TOL < y < y1 - GEOM_TOL):
                authored_segments.append(body["id"])
            if ((abs(y - y0) <= GEOM_TOL or abs(y - y1) <= GEOM_TOL) and
                    x0 + GEOM_TOL < x < x1 - GEOM_TOL):
                authored_segments.append(body["id"])
        if len(authored_segments) != 1:
            raise GmshBackendError(
                "insert-vertex must uniquely select one authored edge interior",
                {"stage": "contact", "position": [x, y],
                 "bodies": authored_segments,
                 "source_location": record["provenance"].get("source_location")},
            )
        point_tags.append(gmsh.model.occ.addPoint(x, y, 0.0))
    if len(source_tags) == 1 and not point_tags:
        mapped = [[(2, source_tags[0])]]
        output = [(2, source_tags[0])]
    else:
        output, mapped = gmsh.model.occ.fragment(
            [(2, source_tags[0])],
            ([(2, tag) for tag in source_tags[1:]] +
             [(0, tag) for tag in point_tags]),
            removeObject=True, removeTool=True,
        )
    gmsh.model.occ.synchronize()
    surfaces = sorted(int(tag) for dimension, tag in output if int(dimension) == 2)
    if not surfaces:
        raise GmshBackendError("OCC fragment produced no positive-area surface")
    mapped_surfaces = [
        {int(tag) for dimension, tag in entries if int(dimension) == 2}
        for entries in mapped[:len(bodies)]
    ]
    if len(mapped_surfaces) != len(bodies):
        raise GmshBackendError(
            "OCC fragment source map length is inconsistent",
            {"stage": "geometry", "sources": len(bodies), "maps": len(mapped_surfaces)},
        )
    vertices = {surface: _surface_vertices(gmsh, surface) for surface in surfaces}
    areas = {surface: float(gmsh.model.occ.getMass(2, surface)) for surface in surfaces}
    for surface, area in areas.items():
        if not math.isfinite(area) or area <= GEOM_TOL * GEOM_TOL:
            raise GmshBackendError(
                "OCC fragment has zero or non-finite area",
                {"stage": "geometry", "surface": surface, "area": area},
            )
    memberships: dict[int, set[int]] = {surface: set() for surface in surfaces}
    source_audit = []
    for index, body in enumerate(bodies):
        x0, y0, x1, y1 = bounds_by_id[body["id"]]
        source_area = (x1 - x0) * (y1 - y0)
        mapped_area = sum(areas[surface] for surface in mapped_surfaces[index])
        tolerance = 1.0e-8 * max(1.0, source_area)
        if abs(mapped_area - source_area) > tolerance:
            raise GmshBackendError(
                "OCC source-to-fragment area audit failed",
                {"stage": "geometry", "body": body["id"],
                 "source_area": source_area, "mapped_area": mapped_area},
            )
        for surface in surfaces:
            contained = all(_point_in_rectangle(bounds_by_id[body["id"]], point)
                            for point in vertices[surface])
            mapped_candidate = surface in mapped_surfaces[index]
            if contained != mapped_candidate:
                raise GmshBackendError(
                    "OCC fragment map disagrees with geometric containment audit",
                    {"stage": "geometry", "body": body["id"], "surface": surface,
                     "mapped": mapped_candidate, "contained": contained},
                )
            if mapped_candidate:
                memberships[surface].add(index)
        source_audit.append({
            "body": body["id"], "source_area": source_area,
            "mapped_area": mapped_area,
            "fragments": sorted(mapped_surfaces[index]),
        })
    surface_owner: dict[int, dict[str, Any]] = {}
    for surface, member_indices in memberships.items():
        if not member_indices:
            raise GmshBackendError(
                "OCC fragment has no source body membership",
                {"stage": "geometry", "surface": surface},
            )
        owner = None
        for index in sorted(member_indices):
            candidate = bodies[index]
            if owner is None or candidate["boolean_mode"] == "ABA":
                owner = candidate
            elif candidate["boolean_mode"] != "BAB":
                raise GmshBackendError(
                    f"unsupported boolean mode {candidate['boolean_mode']!r}",
                    {"stage": "geometry", "body": candidate["id"]},
                )
        assert owner is not None
        surface_owner[surface] = owner
    for record in inserted_vertices:
        x, y = map(float, record["position"])
        endpoint_curves = []
        for curve in _boundary_curves(gmsh, surface_owner):
            contains, endpoint = _curve_contains_pick(gmsh, curve, x, y)
            if contains and endpoint:
                endpoint_curves.append(curve)
        if len(set(endpoint_curves)) < 2:
            raise GmshBackendError(
                "OCC fragment did not retain insert-vertex as a split boundary point",
                {"stage": "contact", "position": [x, y],
                 "curves": sorted(set(endpoint_curves))},
            )
    region_areas: dict[str, float] = defaultdict(float)
    for surface, owner in surface_owner.items():
        region_areas[owner["region"]] += areas[surface]
    return surface_owner, {
        "algorithm": "occ.fragment+source-map+vertex-containment+boolean-replay",
        "source_audit": source_audit,
        "surface_memberships": {
            str(surface): [bodies[index]["id"] for index in sorted(members)]
            for surface, members in memberships.items()
        },
        "surface_owners": {
            str(surface): owner["id"] for surface, owner in surface_owner.items()
        },
        "region_areas_um2": dict(sorted(region_areas.items())),
    }


def _interface_curve_groups(gmsh, surface_owner: dict[int, dict[str, Any]]) -> dict[str, set[int]]:
    result: dict[str, set[int]] = defaultdict(set)
    for curve in _boundary_curves(gmsh, surface_owner):
        upward, _downward = gmsh.model.getAdjacencies(1, curve)
        adjacent = [int(surface) for surface in upward if int(surface) in surface_owner]
        if len(adjacent) > 2:
            raise GmshBackendError(
                "non-manifold curve has more than two adjacent device surfaces",
                {"stage": "geometry", "curve": curve, "surfaces": adjacent},
            )
        owners = {
            surface_owner[surface]["id"]: surface_owner[surface]
            for surface in adjacent
        }
        ordered = sorted(owners.values(), key=lambda body: body["creation_order"])
        regions = [body["region"] for body in ordered]
        if len(set(regions)) == 2:
            result["+".join(regions)].add(curve)
    return result


def configure_gmsh(gmsh, mesh_size: float) -> dict[str, float | int]:
    options: dict[str, float | int] = {
        "General.NumThreads": 1,
        "Mesh.Algorithm": 6,
        "Mesh.RandomFactor": 0,
        "Mesh.RecombineAll": 0,
        "Mesh.ElementOrder": 1,
        "Mesh.MshFileVersion": 4.1,
        "Mesh.Binary": 0,
        "Mesh.Optimize": 0,
        "Mesh.OptimizeNetgen": 0,
        "Mesh.SaveAll": 0,
        "Mesh.CharacteristicLengthMin": mesh_size,
        "Mesh.CharacteristicLengthMax": mesh_size,
    }
    for name, value in options.items():
        gmsh.option.setNumber(name, value)
    return options


def _triangles_for_surface(gmsh, surface: int) -> list[tuple[int, list[int]]]:
    kinds, tag_sets, node_sets = gmsh.model.mesh.getElements(2, surface)
    result: list[tuple[int, list[int]]] = []
    for kind, tags, connectivity in zip(kinds, tag_sets, node_sets):
        if int(kind) != 2:
            raise GmshBackendError(
                "M1 requires first-order Tri3 elements",
                {"stage": "mesh", "element_type": int(kind)},
            )
        if len(connectivity) != 3 * len(tags):
            raise GmshBackendError("invalid Tri3 connectivity length")
        for index, tag in enumerate(tags):
            result.append((int(tag), [int(value) for value in connectivity[3 * index:3 * index + 3]]))
    return result


def _lines_for_curve(gmsh, curve: int) -> list[tuple[int, list[int]]]:
    kinds, tag_sets, node_sets = gmsh.model.mesh.getElements(1, curve)
    result: list[tuple[int, list[int]]] = []
    for kind, tags, connectivity in zip(kinds, tag_sets, node_sets):
        if int(kind) != 1:
            raise GmshBackendError(
                "M1 requires first-order Line2 boundary elements",
                {"stage": "mesh", "element_type": int(kind)},
            )
        for index, tag in enumerate(tags):
            result.append((int(tag), [int(value) for value in connectivity[2 * index:2 * index + 2]]))
    return result


def _signed_area2(nodes: dict[int, tuple[float, float, float]], tags: list[int]) -> float:
    a, b, c = (nodes[tag] for tag in tags)
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _field(values: list[float], unit: str = "cm^-3") -> dict[str, Any]:
    return {"association": "node", "unit": unit, "values": values}


def _evaluate_constant_fields(
    device_ir: dict[str, Any], nodes: list[dict[str, Any]],
    triangles: list[dict[str, Any]], regions: list[dict[str, Any]],
    materials: list[dict[str, Any]],
) -> dict[str, Any]:
    profiles = {item["name"]: item for item in device_ir["doping"]["profiles"]}
    region_name = {item["id"]: item["name"] for item in regions}
    material_name = {item["id"]: item["name"] for item in materials}
    incident_regions: dict[int, set[int]] = defaultdict(set)
    incident_materials: dict[int, set[int]] = defaultdict(set)
    for triangle in triangles:
        for node_id in triangle["node_ids"]:
            incident_regions[node_id].add(triangle["region_id"])
            incident_materials[node_id].add(triangle["material_id"])
    placements = sorted(device_ir["doping"]["placements"], key=lambda item: item["priority"])
    species_names = sorted({profile["species"] for profile in profiles.values()})
    species = {name: [0.0] * len(nodes) for name in species_names}
    carrier = {profile["species"]: profile["carrier"] for profile in profiles.values()}
    for node in nodes:
        node_id = node["vela_id"]
        rnames = {region_name[value] for value in incident_regions[node_id]}
        mnames = {material_name[value] for value in incident_materials[node_id]}
        for placement in placements:
            active = (placement["target"] in rnames if placement["target_kind"] == "region"
                      else placement["target"] in mnames)
            if active:
                profile = profiles[placement["profile"]]
                species[profile["species"]][node_id] += float(profile["value_cm3"])
    donors = [0.0] * len(nodes)
    acceptors = [0.0] * len(nodes)
    for name, values in species.items():
        target = donors if carrier[name] == "donor" else acceptors
        for index, value in enumerate(values):
            target[index] += value
    net = [donor - acceptor for donor, acceptor in zip(donors, acceptors)]
    absolute = [abs(value) for value in net]
    logarithm = [math.log10(max(value, LOG_DOPING_FLOOR_CM3)) for value in absolute]
    return {
        "DonorConcentration": _field(donors),
        "AcceptorConcentration": _field(acceptors),
        "DopingConcentration": _field(net),
        "AbsDopingConcentration": _field(absolute),
        "Log10AbsDopingConcentration": _field(logarithm, "log10(cm^-3)"),
        "species": {name: _field(values) for name, values in species.items()},
    }


def _build_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float,
                  capability_level: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build OCC geometry, physical groups, Tri3 mesh, and nodal fields."""
    if capability_level == "M1":
        validate_m1_device_ir(device_ir)
    elif capability_level == "M2":
        validate_m2_device_ir(device_ir)
    else:
        raise GmshBackendError(f"unknown capability level {capability_level!r}")
    if not math.isfinite(mesh_size) or mesh_size <= 0.0:
        raise GmshBackendError("mesh_size must be positive and finite")
    options = configure_gmsh(gmsh, mesh_size)
    gmsh.model.add(f"vela_sentaurus_gmsh_{capability_level.lower()}")

    bodies = sorted(device_ir["geometry"]["bodies"], key=lambda item: item["creation_order"])
    if capability_level == "M1":
        surface_owner = {}
        for body in bodies:
            x0, y0, x1, y1 = _rectangle(body)
            surface = gmsh.model.occ.addRectangle(x0, y0, 0.0, x1 - x0, y1 - y0)
            surface_owner[surface] = body
        gmsh.model.occ.synchronize()
        geometry_audit = {
            "algorithm": "direct-strictly-separated-rectangles",
            "surface_owners": {
                str(surface): body["id"] for surface, body in surface_owner.items()
            },
            "region_areas_um2": {
                body["region"]: (_rectangle(body)[2] - _rectangle(body)[0]) *
                (_rectangle(body)[3] - _rectangle(body)[1]) for body in bodies
            },
        }
    else:
        surface_owner, geometry_audit = _fragment_m2_geometry(
            gmsh, bodies, device_ir["geometry"]["inserted_vertices"])

    surviving_bodies = {body["id"]: body for body in surface_owner.values()}
    regions = _groups(
        (body["region"] for body in surviving_bodies.values()), "region", REGION_TAG_BASE)
    materials = _groups(
        (body["material"]["sde"] for body in surviving_bodies.values()),
        "material", MATERIAL_TAG_BASE)
    contacts = _groups((item["name"] for item in device_ir["contacts"]), "contact", CONTACT_TAG_BASE)
    interface_curves = _interface_curve_groups(gmsh, surface_owner) if capability_level == "M2" else {}
    interfaces = _groups(interface_curves, "interface", INTERFACE_TAG_BASE)
    region_by_name = {item["name"]: item for item in regions}
    material_by_name = {item["name"]: item for item in materials}
    contact_by_name = {item["name"]: item for item in contacts}
    interface_by_name = {item["name"]: item for item in interfaces}

    for region in regions:
        surfaces = sorted(
            surface for surface, owner in surface_owner.items()
            if owner["region"] == region["name"])
        gmsh.model.addPhysicalGroup(2, surfaces, region["physical_tag"])
        gmsh.model.setPhysicalName(2, region["physical_tag"], region["physical_name"])
    for material in materials:
        surfaces = sorted(
            surface for surface, owner in surface_owner.items()
            if owner["material"]["sde"] == material["name"])
        gmsh.model.addPhysicalGroup(2, surfaces, material["physical_tag"])
        gmsh.model.setPhysicalName(2, material["physical_tag"], material["physical_name"])
    for interface in interfaces:
        curves = sorted(interface_curves[interface["name"]])
        gmsh.model.addPhysicalGroup(1, curves, interface["physical_tag"])
        gmsh.model.setPhysicalName(1, interface["physical_tag"], interface["physical_name"])

    contact_curves: dict[str, set[int]] = defaultdict(set)
    curve_contact: dict[int, str] = {}
    contact_owner: dict[str, set[str]] = defaultdict(set)
    for contact in device_ir["contacts"]:
        for assignment in contact["assignments"]:
            x, y = map(float, assignment["position"])
            candidates: list[tuple[int, bool]] = []
            for curve in _boundary_curves(gmsh, surface_owner):
                contains, endpoint = _curve_contains_pick(gmsh, curve, x, y)
                if contains:
                    candidates.append((curve, endpoint))
            candidates = sorted(set(candidates))
            interior = sorted({curve for curve, endpoint in candidates if not endpoint})
            if len(interior) != 1:
                raise GmshBackendError(
                    f"contact {contact['name']!r} pick point does not select one curve interior",
                    {"stage": "contact", "contact": contact["name"],
                     "pick_point": [x, y], "candidates": candidates,
                     "source_location": assignment["provenance"].get(
                         "source_location")},
                )
            curve = interior[0]
            previous = curve_contact.get(curve)
            if previous is not None and previous != contact["name"]:
                raise GmshBackendError(
                    "two contacts select the same complete curve",
                    {"stage": "contact", "contacts": [previous, contact["name"]],
                     "curve": curve},
                )
            curve_contact[curve] = contact["name"]
            contact_curves[contact["name"]].add(curve)
            upward, _downward = gmsh.model.getAdjacencies(1, curve)
            adjacent = sorted({int(surface) for surface in upward
                               if int(surface) in surface_owner})
            if len(adjacent) != 1:
                raise GmshBackendError(
                    "contact curve must have exactly one adjacent device surface",
                    {"stage": "contact", "contact": contact["name"],
                     "curve": curve, "surfaces": adjacent,
                     "source_location": assignment["provenance"].get("source_location")},
                )
            contact_owner[contact["name"]].add(surface_owner[adjacent[0]]["region"])
    for contact in contacts:
        curves = sorted(contact_curves[contact["name"]])
        gmsh.model.addPhysicalGroup(1, curves, contact["physical_tag"])
        gmsh.model.setPhysicalName(1, contact["physical_tag"], contact["physical_name"])

    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    gmsh.model.mesh.generate(2)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    coords_by_tag = {
        int(tag): (float(coordinates[3 * index]), float(coordinates[3 * index + 1]),
                   float(coordinates[3 * index + 2]))
        for index, tag in enumerate(node_tags)
    }
    ordered_node_tags = sorted(coords_by_tag)
    node_id = {tag: index for index, tag in enumerate(ordered_node_tags)}
    nodes = [
        {"gmsh_tag": tag, "vela_id": index, "dfise_id": index,
         "x": coords_by_tag[tag][0], "y": coords_by_tag[tag][1], "z": coords_by_tag[tag][2]}
        for index, tag in enumerate(ordered_node_tags)
    ]

    region_id = {item["name"]: item["id"] for item in regions}
    material_id = {item["name"]: item["id"] for item in materials}
    raw_triangles: list[dict[str, Any]] = []
    for surface, body in sorted(surface_owner.items()):
        region = region_by_name[body["region"]]
        material = material_by_name[body["material"]["sde"]]
        for tag, connectivity in _triangles_for_surface(gmsh, surface):
            if _signed_area2(coords_by_tag, connectivity) < 0.0:
                connectivity[1], connectivity[2] = connectivity[2], connectivity[1]
            raw_triangles.append({
                "gmsh_tag": tag,
                "node_tags": connectivity,
                "region_id": region_id[body["region"]],
                "material_id": material_id[body["material"]["sde"]],
                "physical_surface_tags": sorted([region["physical_tag"], material["physical_tag"]]),
            })
    triangles = []
    for index, item in enumerate(sorted(raw_triangles, key=lambda entry: entry["gmsh_tag"])):
        triangles.append({
            "gmsh_tag": item["gmsh_tag"], "vela_id": index, "dfise_id": index,
            "node_ids": [node_id[tag] for tag in item["node_tags"]],
            "region_id": item["region_id"], "material_id": item["material_id"],
            "physical_surface_tags": item["physical_surface_tags"],
        })

    raw_lines: list[dict[str, Any]] = []
    for curve, name in sorted(curve_contact.items()):
        physical_tag = contact_by_name[name]["physical_tag"]
        for tag, connectivity in _lines_for_curve(gmsh, curve):
            raw_lines.append({"gmsh_tag": tag, "node_tags": connectivity,
                              "kind": "contact", "physical_curve_tags": [physical_tag]})
    for name, curves in sorted(interface_curves.items()):
        physical_tag = interface_by_name[name]["physical_tag"]
        for curve in sorted(curves):
            for tag, connectivity in _lines_for_curve(gmsh, curve):
                raw_lines.append({
                    "gmsh_tag": tag, "node_tags": connectivity,
                    "kind": "interface", "physical_curve_tags": [physical_tag],
                })
    lines = []
    for index, item in enumerate(sorted(raw_lines, key=lambda entry: entry["gmsh_tag"])):
        lines.append({
            "gmsh_tag": item["gmsh_tag"], "vela_id": index, "dfise_id": index,
            "node_ids": [node_id[tag] for tag in item["node_tags"]],
            "kind": item["kind"], "physical_curve_tags": item["physical_curve_tags"],
        })

    try:
        fields = (_evaluate_constant_fields(device_ir, nodes, triangles, regions, materials)
                  if capability_level == "M1" else evaluate_m2_fields(
                      device_ir, nodes, triangles, regions))
    except DopingEvaluationError as error:
        raise GmshBackendError(str(error), error.details) from error
    bundle = {
        "schema": BUNDLE_SCHEMA,
        "coordinate_unit": "um",
        "field_units": {"concentration": "cm^-3"},
        "provenance": {
            "device_ir_schema": device_ir["schema"], "mesh_backend": "gmsh",
            "source": device_ir["source"]["path"],
            "capability_level": capability_level,
            "capability": ("strictly-separated-rectangles-constant-doping"
                           if capability_level == "M1" else
                           "aba-bab-insert-contact-ordered-replacement-gaussian"),
            "log10_abs_floor_cm3": LOG_DOPING_FLOOR_CM3,
            "geometry_audit": geometry_audit,
            "contact_owners": {
                name: sorted(owners) for name, owners in sorted(contact_owner.items())
            },
            "node_boundary_policy": (
                "boolean-replay-last-ABA-first-BAB" if capability_level == "M1" else
                "sentaurus2022-region-local-sum-at-shared-node"),
            "material_aliases": {
                body["material"]["sde"]: body["material"]["vela"] for body in bodies
            },
        },
        "gmsh": {"version": EXPECTED_GMSH_VERSION, "msh_version": "4.1", "options": options},
        "nodes": nodes, "triangles": triangles, "lines": lines,
        "regions": regions, "materials": materials, "contacts": contacts,
        "interfaces": interfaces,
        "fields": fields,
        "id_mapping": {
            "gmsh_node_to_vela": {str(item["gmsh_tag"]): item["vela_id"] for item in nodes},
            "gmsh_triangle_to_vela": {str(item["gmsh_tag"]): item["vela_id"] for item in triangles},
            "gmsh_line_to_vela": {str(item["gmsh_tag"]): item["vela_id"] for item in lines},
            "vela_node_to_dfise": {str(item["vela_id"]): item["dfise_id"] for item in nodes},
            "vela_triangle_to_dfise": {str(item["vela_id"]): item["dfise_id"] for item in triangles},
            "vela_line_to_dfise": {str(item["vela_id"]): item["dfise_id"] for item in lines},
        },
    }
    qualification = validate_bundle(bundle)
    qualification["geometry"] = geometry_audit
    return bundle, qualification


def build_m1_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float
                    ) -> tuple[dict[str, Any], dict[str, Any]]:
    return _build_bundle(gmsh, device_ir, mesh_size, "M1")


def build_m2_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float
                    ) -> tuple[dict[str, Any], dict[str, Any]]:
    return _build_bundle(gmsh, device_ir, mesh_size, "M2")


def _angles(points: list[tuple[float, float]]) -> list[float]:
    result = []
    for index in range(3):
        p, q, r = points[index], points[(index + 1) % 3], points[(index + 2) % 3]
        u, v = (q[0] - p[0], q[1] - p[1]), (r[0] - p[0], r[1] - p[1])
        cosine = (u[0] * v[0] + u[1] * v[1]) / (math.hypot(*u) * math.hypot(*v))
        result.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
    return result


def validate_bundle(bundle: dict[str, Any]) -> dict[str, Any]:
    """Enforce invariants JSON Schema cannot express."""
    if bundle.get("schema") != BUNDLE_SCHEMA:
        raise GmshBackendError("unexpected mesh bundle schema")
    nodes, triangles, lines = bundle["nodes"], bundle["triangles"], bundle["lines"]
    for name, records in (("node", nodes), ("triangle", triangles), ("line", lines)):
        ids = [item["vela_id"] for item in records]
        if ids != list(range(len(records))):
            raise GmshBackendError(f"{name} IDs are not contiguous from zero")
        tags = [item["gmsh_tag"] for item in records]
        if len(tags) != len(set(tags)) or any(tag <= 0 for tag in tags):
            raise GmshBackendError(f"{name} Gmsh tags are not unique positive integers")
    node_count = len(nodes)
    coords = {item["vela_id"]: (item["x"], item["y"]) for item in nodes}
    all_angles: list[float] = []
    adjacency: dict[int, set[int]] = defaultdict(set)
    node_owner: dict[int, list[int]] = defaultdict(list)
    for triangle in triangles:
        ids = triangle["node_ids"]
        if len(set(ids)) != 3 or any(value < 0 or value >= node_count for value in ids):
            raise GmshBackendError("triangle contains duplicate or invalid node IDs")
        area2 = ((coords[ids[1]][0] - coords[ids[0]][0]) *
                 (coords[ids[2]][1] - coords[ids[0]][1]) -
                 (coords[ids[1]][1] - coords[ids[0]][1]) *
                 (coords[ids[2]][0] - coords[ids[0]][0]))
        if area2 <= 0.0:
            raise GmshBackendError("triangle is degenerate or not CCW")
        all_angles.extend(_angles([coords[value] for value in ids]))
        for value in ids:
            node_owner[value].append(triangle["vela_id"])
    for owners in node_owner.values():
        for first in owners:
            adjacency[first].update(second for second in owners if second != first)
    components = 0
    unseen = set(range(len(triangles)))
    while unseen:
        components += 1
        queue = deque([unseen.pop()])
        while queue:
            for neighbour in adjacency[queue.popleft()] & unseen:
                unseen.remove(neighbour)
                queue.append(neighbour)
    for line in lines:
        if len(set(line["node_ids"])) != 2 or any(
                value < 0 or value >= node_count for value in line["node_ids"]):
            raise GmshBackendError("line contains duplicate or invalid node IDs")
    for name, field in {
        **{key: value for key, value in bundle["fields"].items() if key != "species"},
        **bundle["fields"]["species"],
    }.items():
        if len(field["values"]) != node_count or not all(
                math.isfinite(float(value)) for value in field["values"]):
            raise GmshBackendError(f"nodal field {name!r} has invalid length or values")
    expected_maps = {
        "gmsh_node_to_vela": {str(item["gmsh_tag"]): item["vela_id"] for item in nodes},
        "gmsh_triangle_to_vela": {str(item["gmsh_tag"]): item["vela_id"] for item in triangles},
        "gmsh_line_to_vela": {str(item["gmsh_tag"]): item["vela_id"] for item in lines},
    }
    for name, expected in expected_maps.items():
        if bundle["id_mapping"][name] != expected:
            raise GmshBackendError(f"ID mapping {name!r} is not reversible")
    contact_counts = []
    for contact in bundle["contacts"]:
        nodes_for_contact = {
            value for line in lines
            if contact["physical_tag"] in line["physical_curve_tags"]
            for value in line["node_ids"]
        }
        if len(nodes_for_contact) < 2:
            raise GmshBackendError(f"contact {contact['name']!r} has fewer than two nodes")
        contact_counts.append({"name": contact["name"], "node_count": len(nodes_for_contact)})
    return {
        "node_count": len(nodes), "cell_count": len(triangles),
        "line_count": len(lines), "connected_components": components,
        "min_angle_degrees": min(all_angles), "max_angle_degrees": max(all_angles),
        "all_triangles_ccw": True, "contacts": contact_counts,
        "node_volume_policy": "barycentric",
    }


def write_msh_with_node_data(gmsh, path: Path, bundle: dict[str, Any]) -> None:
    """Write MSH 4.1 ASCII and append five authoritative NodeData views."""
    path.parent.mkdir(parents=True, exist_ok=True)
    gmsh.write(str(path))
    # ``view.write(..., append=True)`` otherwise emits the complete mesh for
    # every view.  The resulting file can be opened by Gmsh, but contains
    # duplicate $Nodes/$Elements sections and duplicate element tags.  Keep
    # the single authoritative mesh written above and append NodeData only.
    gmsh.option.setNumber("PostProcessing.SaveMesh", 0)
    node_tags = [item["gmsh_tag"] for item in bundle["nodes"]]
    for name in (
        "DonorConcentration", "AcceptorConcentration", "DopingConcentration",
        "AbsDopingConcentration", "Log10AbsDopingConcentration",
    ):
        view = gmsh.view.add(name)
        data = [[float(value)] for value in bundle["fields"][name]["values"]]
        gmsh.view.addModelData(
            view, 0, gmsh.model.getCurrent(), "NodeData", node_tags, data, 0.0, 1)
        gmsh.view.write(view, str(path), True)


def read_msh_roundtrip(gmsh, path: Path) -> dict[str, Any]:
    """Reopen MSH through Gmsh and recover topology, groups and NodeData."""
    gmsh.clear()
    gmsh.open(str(path))
    physical_names = sorted(
        gmsh.model.getPhysicalName(dim, tag)
        for dim, tag in gmsh.model.getPhysicalGroups()
    )
    fields: dict[str, dict[int, float]] = {}
    for tag in gmsh.view.getTags():
        index = gmsh.view.getIndex(tag)
        name = gmsh.option.getString(f"View[{index}].Name")
        data_type, node_tags, data, _time, components = gmsh.view.getModelData(tag, 0)
        if data_type != "NodeData" or int(components) != 1:
            raise GmshBackendError(
                "MSH view is not scalar NodeData",
                {"stage": "msh_roundtrip", "view": name, "type": data_type},
            )
        fields[name] = {
            int(node_tag): float(values[0])
            for node_tag, values in zip(node_tags, data)
        }
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    node_coordinates = {
        int(tag): [float(value) for value in coordinates[index:index + 3]]
        for index, tag in zip(range(0, len(coordinates), 3), node_tags)
    }

    elements: dict[str, dict[int, list[int]]] = {"triangles": {}, "lines": {}}
    for dimension, expected_type, width, key in (
        (2, 2, 3, "triangles"), (1, 1, 2, "lines")
    ):
        element_types, element_tags, connectivity = gmsh.model.mesh.getElements(dimension)
        for element_type, tags, nodes in zip(element_types, element_tags, connectivity):
            if int(element_type) != expected_type:
                raise GmshBackendError(
                    "MSH round-trip contains a non-M1 element type",
                    {"stage": "msh_roundtrip", "dimension": dimension,
                     "element_type": int(element_type)},
                )
            for offset, tag in enumerate(tags):
                begin = offset * width
                elements[key][int(tag)] = [
                    int(value) for value in nodes[begin:begin + width]
                ]
    return {
        "physical_names": physical_names,
        "fields": fields,
        "node_coordinates": node_coordinates,
        "elements": elements,
    }
