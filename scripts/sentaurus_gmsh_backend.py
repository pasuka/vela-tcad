#!/usr/bin/env python3
"""Versioned Gmsh OCC backend and normalized mesh/field bundle construction."""

from __future__ import annotations

import copy
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
M5_MESH_STRATEGIES = {"auto", "transfinite-rect", "unstructured"}
M5_TRANSFINITE_DIAGONALS = {"Left", "Right", "AlternateLeft", "AlternateRight"}


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
        if body['geometry'].get('fillets') or body.get('removed'):
            raise GmshBackendError('M1 refuses fillet/lifecycle geometry; use M3')
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
    if (mesh_control["refinement_functions"] or mesh_control["refinement_placements"]
            or mesh_control.get('multibox_placements')):
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
        if body.get('removed') and body['geometry'].get('fillets'):
            raise GmshBackendError('body/remove supports unfilleted Rectangle bodies only')
        if body.get("boolean_mode") not in {"ABA", "BAB"}:
            raise GmshBackendError(
                f"boolean mode {body.get('boolean_mode')!r} is unsupported in M2",
                {"stage": "capability", "capability": "boolean-mode",
                 "body": body["id"], "mode": body.get("boolean_mode")},
            )

    mesh_control = device_ir["mesh_control"]
    if (mesh_control["refinement_functions"] or mesh_control["refinement_placements"]
            or mesh_control.get('multibox_placements')):
        raise GmshBackendError(
            "M2 does not support refinement placements or functions; use M3",
            {"stage": "capability", "capability": "refinement"},
        )
    windows = {item["name"]: item for item in device_ir["doping"]["windows"]}
    profiles = {item["name"]: item for item in device_ir["doping"]["profiles"]}
    for profile in profiles.values():
        if profile.get("type") not in {"constant", "gaussian", "general"}:
            raise GmshBackendError(
                f"unsupported doping profile type {profile.get('type')!r}",
                {"stage": "capability", "capability": "doping-profile",
                 "profile": profile.get("name")},
            )
    for profile in profiles.values():
        if profile["type"] == "general":
            from sentaurus_general import validate_profile
            try:
                validate_profile(profile)
            except (ValueError, KeyError, TypeError) as error:
                raise GmshBackendError(f"invalid General profile: {error}", {"stage": "doping"}) from error
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
        if profile["type"] == "general" and (
                placement.get("direction") != "Both" or placement.get("evaluation") != "Eval" or
                placement.get("restriction_decay", "none") != "none" or
                abs(float(placement.get("restriction_decay_length", 0.0))) > GEOM_TOL):
            source = placement.get("provenance", {}).get("source_location", "<unknown>")
            raise GmshBackendError(
                f"{source}: General requires Both/Eval without restriction decay",
                {"stage": "capability", "capability": "general-placement", "source_location": source})
        if profile["type"] in {"gaussian", "general"} and (
                kind != "window" or windows[target].get("shape") != "line"):
            raise GmshBackendError(
                "M2 analytical placement requires a Line baseline window",
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
        for assignment in contact['assignments']:
            if assignment['kind']=='body_boundary':
                matches=[b for b in bodies if assignment['entity_ids']==[b['creation_order']]]
                if len(matches)!=1 or not matches[0].get('removed') or matches[0]['removed']['contact']!=contact['name']:
                    raise GmshBackendError('body_boundary requires one matching removed Rectangle')
        if any(item["kind"] not in {"pick_point","body_boundary"} for item in contact["assignments"]):
            raise GmshBackendError(
                "M2 supports pick points and bounded Rectangle body/remove contacts only",
                {"stage": "capability", "capability": "contact-entity-list",
                 "contact": contact["name"]},
            )


def validate_m3_device_ir(device_ir: dict[str, Any]) -> None:
    """Validate the M2 core plus the phase-4A refinement subset."""
    core = copy.deepcopy(device_ir)
    core["mesh_control"]["refinement_placements"] = []
    core["mesh_control"]["refinement_functions"] = []
    core["mesh_control"]["multibox_placements"] = []
    validate_m2_device_ir(core)

    control = device_ir["mesh_control"]
    sizes = {item["name"]: item for item in control["refinement_sizes"]}
    windows = {item["name"]: item for item in control["windows"]}
    regions = {body["region"] for body in device_ir["geometry"]["bodies"]}
    materials = {
        body["material"]["sde"] for body in device_ir["geometry"]["bodies"]
    }
    for size in sizes.values():
        if size.get("provenance", {}).get("status") != "known":
            raise GmshBackendError(
                "M3 refinement size requires known source provenance",
                {"stage": "capability", "capability": "refinement-provenance",
                 "refinement": size["name"]},
            )
        for maximum, minimum in (("max_x", "min_x"), ("max_y", "min_y")):
            if float(size[maximum]) <= 0.0:
                raise GmshBackendError("refinement maximum spacing must be positive")
            if minimum in size and not (0.0 < float(size[minimum]) <= float(size[maximum])):
                raise GmshBackendError(
                    "refinement minimum spacing must be positive and no larger than maximum",
                    {"stage": "refinement", "refinement": size["name"],
                     "axis": maximum[-1]},
                )
    from sentaurus_multibox import validate
    validate(device_ir)
    for placement in control["refinement_placements"]:
        if placement.get("provenance", {}).get("status") != "known":
            raise GmshBackendError(
                "M3 refinement placement requires known source provenance",
                {"stage": "capability", "capability": "refinement-provenance",
                 "placement": placement["name"]},
            )
        if placement["refinement"] not in sizes:
            raise GmshBackendError(
                "refinement placement references an unknown definition",
                {"stage": "refinement", "placement": placement["name"]},
            )
        kind, target = placement["target_kind"], placement["target"]
        valid = ((kind == "region" and target in regions) or
                 (kind == "material" and target in materials) or
                 (kind == "window" and target in windows and
                  windows[target].get("shape") in {"line", "rectangle", "polygon"}))
        if not valid:
            raise GmshBackendError(
                f"invalid M3 refinement target {kind!r}:{target!r}",
                {"stage": "refinement", "placement": placement["name"]},
            )
    for function in control["refinement_functions"]:
        if function.get("provenance", {}).get("status") != "known":
            raise GmshBackendError(
                "M3 refinement function requires known source provenance",
                {"stage": "capability", "capability": "refinement-provenance",
                 "refinement": function["refinement"]},
            )
        if function["function"] == "MaxTransDiff":
            from sentaurus_adaptive import validate
            validate(device_ir)
            continue
        if function["function"] == "MaxGradient":
            raise GmshBackendError(
                f"{function['function']} belongs to M3 phase 4B and has no "
                "declared exact Gmsh mapping",
                {"stage": "capability", "capability": "refinement-4b",
                 "function": function["function"],
                 "source_location": function["provenance"]["source_location"]},
            )
        if function["function"] != "MaxLenInt":
            raise GmshBackendError(
                f"unsupported refinement function {function['function']!r}",
                {"stage": "capability", "capability": "refinement-function"},
            )
        first, second = function["interface"]
        namespace = regions if function["use_region_names"] else materials
        if first != "All" and first not in namespace:
            raise GmshBackendError(
                f"MaxLenInt first interface name {first!r} is unknown",
                {"stage": "refinement", "refinement": function["refinement"]},
            )
        if second not in {"", "All", "Contact"} and second not in namespace:
            raise GmshBackendError(
                f"MaxLenInt second interface name {second!r} is unknown",
                {"stage": "refinement", "refinement": function["refinement"]},
            )
        if second == "Contact":
            raise GmshBackendError(
                "M3 MaxLenInt contact interfaces are not implemented",
                {"stage": "capability", "capability": "maxlenint-contact",
                 "source_location": function["provenance"]["source_location"]},
            )


def validate_m5_device_ir(device_ir: dict[str, Any]) -> None:
    """M5 preserves the complete M3 semantic contract and adds mesh policy."""
    validate_m3_device_ir(device_ir)


def _m5_transfinite_eligibility(device_ir: dict[str, Any]) -> dict[str, Any]:
    """Return an auditable, deliberately narrow structured-mesh decision."""
    reasons: list[dict[str, Any]] = []
    bodies = device_ir["geometry"]["bodies"]
    if len(bodies) != 1:
        reasons.append({
            "code": "body_count", "message": "requires exactly one authored body",
            "actual": len(bodies),
        })
    elif bodies[0].get("shape") != "rectangle":
        reasons.append({
            "code": "body_shape", "message": "requires an axis-aligned Rectangle body",
            "actual": bodies[0].get("shape"),
        })
    if device_ir['mesh_control'].get('multibox_placements'):
        reasons.append({'code':'multibox','message':'graded metrics require unstructured mesh'})
    if any(b['geometry'].get('fillets') for b in bodies):
        reasons.append({'code':'fillets','message':'curved outline requires unstructured mesh'})
    if device_ir["geometry"].get("inserted_vertices"):
        reasons.append({
            "code": "inserted_vertices",
            "message": "authored insert-vertex topology is not yet mapped to a structured grid",
        })
    control = device_ir["mesh_control"]
    if control.get("refinement_functions"):
        reasons.append({
            "code": "refinement_functions",
            "message": "Transfinite M5 supports spacing placements, not refinement functions",
            "functions": [item.get("function") for item in control["refinement_functions"]],
        })
    windows = {item["name"]: item for item in control.get("windows", [])}
    for placement in control.get("refinement_placements", []):
        if placement.get("target_kind") != "window":
            continue
        window = windows.get(placement.get("target"))
        if window is None or window.get("shape") != "rectangle":
            reasons.append({
                "code": "window_shape",
                "message": "Transfinite M5 requires Rectangle refinement windows",
                "placement": placement.get("name"),
                "target": placement.get("target"),
                "shape": None if window is None else window.get("shape"),
            })
    return {
        "eligible": not reasons,
        "requirements": [
            "one axis-aligned Rectangle body", "no authored insert-vertex",
            "Rectangle refinement windows only", "no refinement functions",
        ],
        "reasons": reasons,
    }


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
    clip_tags=[]
    bounds_by_id = {}
    for body in bodies:
        x0,y0,x1,y1=_rectangle(body)
        bounds_by_id[body['id']]=(x0,y0,x1,y1)
        source_tags.append(gmsh.model.occ.addRectangle(x0,y0,0,x1-x0,y1-y0))
        if body['geometry'].get('fillets'):
            from sentaurus_curved import polygon
            points=[gmsh.model.occ.addPoint(x,y,0) for x,y in polygon(body)]
            lines=[gmsh.model.occ.addLine(a,b) for a,b in zip(points,points[1:]+points[:1])]
            clip_tags.append(gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addCurveLoop(lines)]))
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
    if len(source_tags) == 1 and not point_tags and not clip_tags:
        mapped = [[(2, source_tags[0])]]
        output = [(2, source_tags[0])]
    else:
        output, mapped = gmsh.model.occ.fragment(
            [(2, source_tags[0])],
            ([(2, tag) for tag in source_tags[1:]] +
             [(0, tag) for tag in point_tags]+[(2,tag) for tag in clip_tags]),
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
        from sentaurus_curved import area as body_area
        source_area = (x1-x0)*(y1-y0)
        mapped_area = sum(areas[surface] for surface in mapped_surfaces[index])
        tolerance = 1.0e-8 * max(1.0, source_area)
        if abs(mapped_area - source_area) > tolerance:
            raise GmshBackendError(
                "OCC source-to-fragment area audit failed",
                {"stage": "geometry", "body": body["id"],
                 "source_area": source_area, "mapped_area": mapped_area},
            )
        for surface in surfaces:
            from sentaurus_curved import contains
            contained = all(_point_in_rectangle(bounds_by_id[body['id']],point) for point in vertices[surface])
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
    discarded=[]
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
                from sentaurus_curved import contains
                owner = None if candidate.get('removed') or not all(
                    contains(candidate,point,1e-7) for point in vertices[surface]) else candidate
            elif candidate["boolean_mode"] != "BAB":
                raise GmshBackendError(
                    f"unsupported boolean mode {candidate['boolean_mode']!r}",
                    {"stage": "geometry", "body": candidate["id"]},
                )
        if owner is None:
            discarded.append(surface)
        else:surface_owner[surface] = owner
    if discarded:
        gmsh.model.occ.remove([(2,t) for t in discarded],recursive=True)
        gmsh.model.occ.synchronize()
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


def configure_gmsh(gmsh, mesh_size: float,
                   anisotropic: bool = False) -> dict[str, float | int]:
    options: dict[str, float | int] = {
        "General.NumThreads": 1,
        # Gmsh's official t17 example requires BAMG for a 2-D anisotropic
        # metric background field.  M1/M2 retain Frontal-Delaunay.
        "Mesh.Algorithm": 7 if anisotropic else 6,
        "Mesh.RandomFactor": 0,
        "Mesh.RecombineAll": 0,
        "Mesh.ElementOrder": 1,
        "Mesh.MshFileVersion": 4.1,
        "Mesh.Binary": 0,
        "Mesh.Optimize": 0,
        "Mesh.OptimizeNetgen": 0,
        "Mesh.SaveAll": 0,
        "Mesh.CharacteristicLengthMin": 0.0 if anisotropic else mesh_size,
        "Mesh.CharacteristicLengthMax": mesh_size,
    }
    if anisotropic:
        options.update({
            "Mesh.AnisoMax": 1000,
            "Mesh.MeshSizeExtendFromBoundary": 0,
            "Mesh.MeshSizeFromPoints": 0,
            "Mesh.MeshSizeFromCurvature": 0,
            "Mesh.SmoothRatio": 3,
        })
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


def _number(value: float) -> str:
    return format(float(value), ".17g")


def _union_conditions(conditions: Iterable[str]) -> str:
    values = [f"({value})" for value in conditions if value and value != "0"]
    if not values:
        return "0"
    if len(values) == 1:
        return values[0]
    return f"step(({'+'.join(values)})-0.5)"


def _rectangle_condition(x0: float, y0: float, x1: float, y1: float) -> str:
    tolerance = GEOM_TOL * max(abs(x0), abs(y0), abs(x1), abs(y1), 1.0)
    return "*".join((
        f"step(x-({_number(x0 - tolerance)}))",
        f"step(({_number(x1 + tolerance)})-x)",
        f"step(y-({_number(y0 - tolerance)}))",
        f"step(({_number(y1 + tolerance)})-y)",
    ))


def _cross(first: tuple[float, float], second: tuple[float, float],
           third: tuple[float, float]) -> float:
    return ((second[0] - first[0]) * (third[1] - first[1]) -
            (second[1] - first[1]) * (third[0] - first[0]))


def _triangulate_polygon(points: list[list[float]]) -> list[list[tuple[float, float]]]:
    vertices = [tuple(map(float, point)) for point in points]
    if len(vertices) < 3 or sum(
            _cross((0.0, 0.0), vertices[index], vertices[(index + 1) % len(vertices)])
            for index in range(len(vertices))) <= 0.0:
        raise GmshBackendError(
            "Polygon Ref/Eval window must be nondegenerate and counterclockwise",
            {"stage": "refinement", "capability": "polygon-window"},
        )
    remaining = list(range(len(vertices)))
    triangles: list[list[tuple[float, float]]] = []
    while len(remaining) > 3:
        clipped = False
        for offset, current in enumerate(remaining):
            previous = remaining[offset - 1]
            following = remaining[(offset + 1) % len(remaining)]
            a, b, c = vertices[previous], vertices[current], vertices[following]
            if _cross(a, b, c) <= GEOM_TOL:
                continue
            contains = False
            for candidate in remaining:
                if candidate in {previous, current, following}:
                    continue
                point = vertices[candidate]
                if (_cross(a, b, point) >= -GEOM_TOL and
                        _cross(b, c, point) >= -GEOM_TOL and
                        _cross(c, a, point) >= -GEOM_TOL):
                    contains = True
                    break
            if contains:
                continue
            triangles.append([a, b, c])
            del remaining[offset]
            clipped = True
            break
        if not clipped:
            raise GmshBackendError(
                "Polygon Ref/Eval window is self-intersecting or cannot be triangulated",
                {"stage": "refinement", "capability": "polygon-window"},
            )
    triangles.append([vertices[index] for index in remaining])
    return triangles


def _triangle_condition(points: list[tuple[float, float]]) -> str:
    tolerance = GEOM_TOL * max(
        1.0, *(abs(value) for point in points for value in point))
    terms = []
    for first, second in zip(points, points[1:] + points[:1]):
        dx, dy = second[0] - first[0], second[1] - first[1]
        expression = (
            f"({_number(dx)})*(y-({_number(first[1])}))-"
            f"({_number(dy)})*(x-({_number(first[0])}))")
        terms.append(f"step(({expression})+{_number(tolerance)})")
    return "*".join(terms)


def _line_distance_expression(window: dict[str, Any]) -> str:
    (x0, y0), (x1, y1) = map(tuple, (window["start"], window["end"]))
    dx, dy = float(x1) - float(x0), float(y1) - float(y0)
    length2 = dx * dx + dy * dy
    if length2 <= GEOM_TOL * GEOM_TOL:
        raise GmshBackendError("Line Ref/Eval window has zero length")
    projection = (
        f"(((x-({_number(x0)}))*({_number(dx)})+"
        f"(y-({_number(y0)}))*({_number(dy)}))/({_number(length2)}))")
    start_distance2 = f"((x-({_number(x0)}))^2+(y-({_number(y0)}))^2)"
    end_distance2 = f"((x-({_number(x1)}))^2+(y-({_number(y1)}))^2)"
    line_distance2 = (
        f"((({_number(dx)})*(y-({_number(y0)}))-"
        f"({_number(dy)})*(x-({_number(x0)})))^2/({_number(length2)}))")
    epsilon = "1e-12"
    return (
        "sqrt("
        f"step(-({projection})-{epsilon})*{start_distance2}+"
        f"step(({projection})+{epsilon})*step(1+{epsilon}-({projection}))*"
        f"{line_distance2}+step(({projection})-1-{epsilon})*{end_distance2})")


def _window_condition(window: dict[str, Any], line_width: float | None = None) -> str:
    shape = window["shape"]
    if shape == "rectangle":
        x0, y0 = map(float, window["lower_left"])
        x1, y1 = map(float, window["upper_right"])
        return _rectangle_condition(x0, y0, x1, y1)
    if shape == "polygon":
        return _union_conditions(
            _triangle_condition(triangle)
            for triangle in _triangulate_polygon(window["points"])
        )
    if shape == "line" and line_width is not None:
        return f"step({_number(line_width)}-({_line_distance_expression(window)}))"
    raise GmshBackendError(
        f"unsupported refinement window shape {shape!r}",
        {"stage": "capability", "capability": "refinement-window"},
    )


def _ownership_conditions(device_ir: dict[str, Any]) -> tuple[
        dict[str, str], dict[str, str]]:
    regions: dict[str, str] = {}
    occupied = "0"
    material_for_region: dict[str, str] = {}
    for body in sorted(device_ir["geometry"]["bodies"],
                       key=lambda item: item["creation_order"]):
        x0, y0, x1, y1 = _rectangle(body)
        condition = _rectangle_condition(x0, y0, x1, y1)
        region = body["region"]
        material_for_region[region] = body["material"]["sde"]
        if body["boolean_mode"] == "ABA":
            outside = f"step(0.5-({condition}))"
            regions = {name: f"({value})*({outside})"
                       for name, value in regions.items()}
            available = condition
        else:
            available = f"({condition})*step(0.5-({occupied}))"
        if not body.get('removed'):
            regions[region] = _union_conditions((regions.get(region, "0"), available))
        occupied = _union_conditions(regions.values())
    materials: dict[str, str] = defaultdict(str)
    for region, condition in regions.items():
        material = material_for_region[region]
        materials[material] = _union_conditions((materials.get(material, "0"), condition))
    return regions, dict(materials)


def _add_metric_field(gmsh, m11: str, m22: str) -> int:
    tag = int(gmsh.model.mesh.field.add("MathEvalAniso"))
    for option, expression in (
        ("M11", m11), ("M12", "0"), ("M13", "0"),
        ("M22", m22), ("M23", "0"), ("M33", "1"),
    ):
        gmsh.model.mesh.field.setString(tag, option, expression)
    return tag


def _target_condition(placement: dict[str, Any], size: dict[str, Any],
                      windows: dict[str, dict[str, Any]],
                      regions: dict[str, str], materials: dict[str, str]) -> str:
    kind, target = placement["target_kind"], placement["target"]
    if kind == "region":
        return regions.get(target,"0")
    if kind == "material":
        return materials.get(target,"0")
    line_width = max(float(size["max_x"]), float(size["max_y"]))
    if windows[target]["shape"] == "line":
        # A zero-measure field does not constrain adjacent 2-D triangles.
        # Extend the metric over three requested edge lengths; qualification
        # below still samples only the one-length core band.
        line_width *= 3.0
    return _window_condition(windows[target], line_width)


def _matching_interface_curves(gmsh, surface_owner: dict[int, dict[str, Any]],
                               function: dict[str, Any]) -> list[int]:
    first, second = function["interface"]
    use_regions = bool(function["use_region_names"])
    matches = []
    for curve in _boundary_curves(gmsh, surface_owner):
        upward, _downward = gmsh.model.getAdjacencies(1, curve)
        owners = [surface_owner[int(surface)] for surface in upward
                  if int(surface) in surface_owner]
        names = [owner["region"] if use_regions else owner["material"]["sde"]
                 for owner in owners]
        first_present = first == "All" or first in names
        if second == "":
            selected = len(owners) == 1 and first_present
        elif second == "All":
            selected = len(owners) == 2 and first_present
        else:
            selected = len(owners) == 2 and first_present and second in names and first != second
        if selected:
            matches.append(curve)
    if not matches:
        raise GmshBackendError(
            "MaxLenInt did not match any device interface curve",
            {"stage": "refinement", "refinement": function["refinement"],
             "interface": function["interface"],
             "source_location": function["provenance"]["source_location"]},
        )
    return sorted(set(matches))


def _build_m3_refinement_fields(
    gmsh, device_ir: dict[str, Any], surface_owner: dict[int, dict[str, Any]],
    mesh_size: float,
) -> dict[str, Any]:
    control = device_ir["mesh_control"]
    sizes = {item["name"]: item for item in control["refinement_sizes"]}
    windows = {item["name"]: item for item in control["windows"]}
    regions, materials = _ownership_conditions(device_ir)
    coarse_metric = 1.0 / (mesh_size * mesh_size)
    base = _add_metric_field(gmsh, _number(coarse_metric), _number(coarse_metric))
    tags = [base]
    nodes: list[dict[str, Any]] = [{
        "id": "base", "kind": "MathEvalAniso", "gmsh_field_tag": base,
        "max_spacing": [mesh_size, mesh_size],
    }]
    placements_by_refinement: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for placement in control["refinement_placements"]:
        size = sizes[placement["refinement"]]
        condition = _target_condition(placement, size, windows, regions, materials)
        mx = 1.0 / float(size["max_x"]) ** 2
        my = 1.0 / float(size["max_y"]) ** 2
        m11 = f"{_number(coarse_metric)}+({_number(mx - coarse_metric)})*({condition})"
        m22 = f"{_number(coarse_metric)}+({_number(my - coarse_metric)})*({condition})"
        tag = _add_metric_field(gmsh, m11, m22)
        tags.append(tag)
        placements_by_refinement[placement["refinement"]].append(placement)
        nodes.append({
            "id": f"placement:{placement['name']}", "kind": "MathEvalAniso",
            "gmsh_field_tag": tag, "refinement": placement["refinement"],
            "target": {"kind": placement["target_kind"], "name": placement["target"]},
            "max_spacing": [size["max_x"], size["max_y"]],
            "min_spacing_stop_target": [size.get("min_x"), size.get("min_y")],
            "condition": condition,
            "source_location": placement["provenance"]["source_location"],
        })
    for index, function in enumerate(control["refinement_functions"]):
        if function["function"] == "MaxTransDiff":
            continue
        curves = _matching_interface_curves(gmsh, surface_owner, function)
        distance = int(gmsh.model.mesh.field.add("Distance"))
        gmsh.model.mesh.field.setNumbers(distance, "CurvesList", curves)
        gmsh.model.mesh.field.setNumber(distance, "Sampling", 100)
        nodes.append({
            "id": f"maxlenint-distance:{index}", "kind": "Distance",
            "gmsh_field_tag": distance, "curves": curves,
        })
        value = float(function["value"])
        factor = function.get("factor")
        if factor is None:
            h = (f"{_number(value)}+({_number(mesh_size - value)})*"
                 f"step(F{distance}-{_number(value)})")
            growth = "one-layer-default"
        elif float(factor) == 1.0:
            h = _number(value)
            growth = "constant"
        else:
            h = (f"min({_number(mesh_size)},"
                 f"{_number(value)}*{_number(float(factor))}^(F{distance}/{_number(value)}))")
            growth = "geometric"
        scope = _union_conditions(
            _target_condition(placement, sizes[placement["refinement"]],
                              windows, regions, materials)
            for placement in placements_by_refinement.get(function["refinement"], [])
        )
        if scope == "0":
            nodes.append({
                "id": f"maxlenint:{index}", "kind": "inactive-definition",
                "refinement": function["refinement"], "reason": "no placement",
            })
            continue
        if not function["double_side"]:
            first = function["interface"][0]
            if first == "All":
                raise GmshBackendError(
                    "single-sided MaxLenInt cannot use 'All' as its first interface name",
                    {"stage": "capability", "capability": "maxlenint-all-side"},
                )
            side = regions[first] if function["use_region_names"] else materials[first]
            scope = f"({scope})*({side})"
        fine_metric = f"1/(({h})^2)"
        m11 = f"{_number(coarse_metric)}+(({fine_metric})-{_number(coarse_metric)})*({scope})"
        tag = _add_metric_field(gmsh, m11, m11)
        tags.append(tag)
        nodes.append({
            "id": f"maxlenint:{index}", "kind": "MathEvalAniso",
            "gmsh_field_tag": tag, "distance_field": distance,
            "curves": curves,
            "refinement": function["refinement"], "interface": function["interface"],
            "value": value, "factor": factor, "growth": growth,
            "double_side": function["double_side"],
            "use_region_names": function["use_region_names"],
            "condition": scope,
            "source_location": function["provenance"]["source_location"],
        })
    for placement in control.get('multibox_placements',[]):
        size=next(s for s in control['multibox_sizes'] if s['name']==placement['multibox'])
        window=windows[placement['window']];x0,y0=window['lower_left'];x1,y1=window['upper_right']
        active=False
        for surface in surface_owner:
            points=_surface_vertices(gmsh,surface)
            xs=[p[0] for p in points];ys=[p[1] for p in points]
            if min(x1,max(xs))-max(x0,min(xs))>1e-10 and min(y1,max(ys))-max(y0,min(ys))>1e-10:active=True
        if not active:
            nodes.append({'id':'multibox:'+placement['name'],'kind':'inactive-multibox','reason':'window has no surviving volume','definition':size,'window':window})
            continue
        condition=_window_condition(window)
        expressions=[]
        for axis,lo,hi in [('x',x0,x1),('y',y0,y1)]:
            ratio=size['ratio_'+axis];minimum=size['min_'+axis];maximum=size['max_'+axis]
            if ratio==0:h=_number(mesh_size)
            elif abs(ratio)==1:h=_number(maximum)
            else:
                distance=f'({axis}-({_number(lo)}))' if ratio>0 else f'(({_number(hi)})-{axis})'
                h=f'min({_number(maximum)},{_number(minimum)}+{_number(abs(ratio)-1)}*max(0,{distance}))'
            metric=f'({_number(coarse_metric)})+(1/(({h})^2)-({_number(coarse_metric)}))*({condition})'
            expressions.append(metric)
        tag=_add_metric_field(gmsh,*expressions);tags.append(tag)
        nodes.append({'id':'multibox:'+placement['name'],'kind':'MathEvalAniso','gmsh_field_tag':tag,'definition':size,'window':window,'spacing_expressions':expressions})
    root = int(gmsh.model.mesh.field.add("MinAniso"))
    gmsh.model.mesh.field.setNumbers(root, "FieldsList", tags)
    gmsh.model.mesh.field.setAsBackgroundMesh(root)
    nodes.append({
        "id": "root", "kind": "MinAniso", "gmsh_field_tag": root,
        "inputs": tags, "combination": "finest directional metric intersection",
    })
    return {
        "schema": "vela.sentaurus_gmsh_refinement_field_dag.v1",
        "engine": "Gmsh MathEvalAniso/Distance/MinAniso with BAMG",
        "nodes": nodes, "root_field_tag": root,
        "min_spacing_semantics": "refinement stop target, not a hard lower clamp",
        "unsupported_fail_closed": ["MaxGradient", "MaxInterval", "offsetting"],
    }


def _point_owner(device_ir: dict[str, Any], x: float, y: float) -> dict[str, Any] | None:
    owner = None
    for body in sorted(device_ir["geometry"]["bodies"],
                       key=lambda item: item["creation_order"]):
        x0, y0, x1, y1 = _rectangle(body)
        inside = (x0 - GEOM_TOL <= x <= x1 + GEOM_TOL and
                  y0 - GEOM_TOL <= y <= y1 + GEOM_TOL)
        if not inside:
            continue
        if owner is None or body["boolean_mode"] == "ABA":
            from sentaurus_curved import contains
            owner = None if body.get('removed') or not contains(body,(x,y)) else body
    return owner


def _point_in_polygon(points: list[list[float]], x: float, y: float) -> bool:
    inside = False
    previous = points[-1]
    for current in points:
        x0, y0 = map(float, previous)
        x1, y1 = map(float, current)
        if abs(_cross((x0, y0), (x1, y1), (x, y))) <= GEOM_TOL and \
                min(x0, x1) - GEOM_TOL <= x <= max(x0, x1) + GEOM_TOL and \
                min(y0, y1) - GEOM_TOL <= y <= max(y0, y1) + GEOM_TOL:
            return True
        if (y0 > y) != (y1 > y):
            intersection = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if x <= intersection + GEOM_TOL:
                inside = not inside
        previous = current
    return inside


def _point_in_window(window: dict[str, Any], x: float, y: float,
                     line_width: float) -> bool:
    if window["shape"] == "rectangle":
        x0, y0 = map(float, window["lower_left"])
        x1, y1 = map(float, window["upper_right"])
        return x0 - GEOM_TOL <= x <= x1 + GEOM_TOL and \
            y0 - GEOM_TOL <= y <= y1 + GEOM_TOL
    if window["shape"] == "polygon":
        return _point_in_polygon(window["points"], x, y)
    if window["shape"] == "line":
        (x0, y0), (x1, y1) = window["start"], window["end"]
        dx, dy = float(x1) - float(x0), float(y1) - float(y0)
        length2 = dx * dx + dy * dy
        fraction = max(0.0, min(1.0, (
            (x - float(x0)) * dx + (y - float(y0)) * dy) / length2))
        return math.hypot(
            x - (float(x0) + fraction * dx),
            y - (float(y0) + fraction * dy),
        ) <= line_width + GEOM_TOL
    return False


def _placement_active_at(device_ir: dict[str, Any], placement: dict[str, Any],
                         size: dict[str, Any], windows: dict[str, dict[str, Any]],
                         x: float, y: float) -> bool:
    owner = _point_owner(device_ir, x, y)
    if placement["target_kind"] == "region":
        return owner is not None and owner["region"] == placement["target"]
    if placement["target_kind"] == "material":
        return owner is not None and owner["material"]["sde"] == placement["target"]
    return _point_in_window(
        windows[placement["target"]], x, y,
        max(float(size["max_x"]), float(size["max_y"])),
    )


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[min(len(ordered) - 1, int(math.ceil(fraction * len(ordered))) - 1)]


def _qualify_m3_refinement(
    gmsh, device_ir: dict[str, Any], nodes: list[dict[str, Any]],
    triangles: list[dict[str, Any]], field_dag: dict[str, Any],
) -> dict[str, Any]:
    coordinates = {node["vela_id"]: (float(node["x"]), float(node["y"]))
                   for node in nodes}
    edges = {
        tuple(sorted((triangle["node_ids"][index],
                      triangle["node_ids"][(index + 1) % 3])))
        for triangle in triangles for index in range(3)
    }
    control = device_ir["mesh_control"]
    sizes = {item["name"]: item for item in control["refinement_sizes"]}
    windows = {item["name"]: item for item in control["windows"]}
    reports = []
    failures = []
    for placement in control["refinement_placements"]:
        size = sizes[placement["refinement"]]
        selected = []
        for first, second in edges:
            a, b = coordinates[first], coordinates[second]
            midpoint = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))
            if _placement_active_at(
                    device_ir, placement, size, windows, midpoint[0], midpoint[1]):
                selected.append((abs(a[0] - b[0]), abs(a[1] - b[1])))
        if not selected:
            if placement['target_kind']=='region' and any(b['region']==placement['target'] and b.get('removed') for b in device_ir['geometry']['bodies']):
                reports.append({'name':placement['name'],'active':False,'reason':'region removed before meshing'})
                continue
            failures.append(f"refinement placement {placement['name']} selected no mesh edges")
            continue
        p95_x = _percentile([value[0] for value in selected], 0.95)
        p95_y = _percentile([value[1] for value in selected], 0.95)
        limit_x, limit_y = float(size["max_x"]), float(size["max_y"])
        passed = p95_x <= 1.35 * limit_x and p95_y <= 1.35 * limit_y
        if not passed:
            failures.append(
                f"refinement placement {placement['name']} P95 edge components "
                "exceed the 1.35x BAMG conformance envelope")
        reports.append({
            "name": placement["name"], "target_kind": placement["target_kind"],
            "target": placement["target"], "edge_count": len(selected),
            "requested_max_spacing": [limit_x, limit_y],
            "p95_edge_component": [p95_x, p95_y],
            "max_edge_component": [max(value[0] for value in selected),
                                   max(value[1] for value in selected)],
            "min_spacing_stop_target": [size.get("min_x"), size.get("min_y")],
            "passed": passed,
        })
    interface_reports = []
    for node in field_dag["nodes"]:
        if not node["id"].startswith("maxlenint:") or node["kind"] != "MathEvalAniso":
            continue
        lengths = []
        for curve in node["curves"]:
            for _tag, connectivity in _lines_for_curve(gmsh, int(curve)):
                first = gmsh.model.mesh.getNode(int(connectivity[0]))[0]
                second = gmsh.model.mesh.getNode(int(connectivity[1]))[0]
                lengths.append(math.hypot(
                    float(first[0]) - float(second[0]),
                    float(first[1]) - float(second[1])))
        maximum = max(lengths) if lengths else math.inf
        passed = maximum <= 1.35 * float(node["value"])
        if not passed:
            failures.append(
                f"MaxLenInt {node['refinement']} interface edges exceed target envelope")
        interface_reports.append({
            "refinement": node["refinement"], "interface": node["interface"],
            "line_element_count": len(lengths), "requested_length": node["value"],
            "max_interface_edge_length": maximum, "passed": passed,
        })
    from sentaurus_multibox import qualify
    multibox = qualify(device_ir,coordinates,triangles,field_dag)
    failures.extend(multibox['failures'])
    return {
        "multibox": multibox['placements'],
        "status": "passed" if not failures else "failed",
        "field_dag_schema": field_dag["schema"],
        "root_field_tag": field_dag["root_field_tag"],
        "anisotropic_algorithm": "BAMG",
        "edge_conformance_envelope": 1.35,
        "placements": reports, "max_len_int": interface_reports,
        "failures": failures,
    }


def _m5_rectangle_window(window: dict[str, Any]) -> tuple[float, float, float, float]:
    lower = window.get("lower_left", window.get("start"))
    upper = window.get("upper_right", window.get("end"))
    if lower is None or upper is None:
        raise GmshBackendError(
            "Rectangle refinement window has no bounds",
            {"stage": "mesh_strategy", "window": window.get("name")},
        )
    x0, y0 = map(float, lower)
    x1, y1 = map(float, upper)
    return min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1)


def _m5_segment_spacing(device_ir: dict[str, Any], body: dict[str, Any],
                        axis: str, midpoint: tuple[float, float],
                        mesh_size: float) -> float:
    control = device_ir["mesh_control"]
    sizes = {item["name"]: item for item in control["refinement_sizes"]}
    windows = {item["name"]: item for item in control["windows"]}
    spacing = mesh_size
    for placement in control["refinement_placements"]:
        applies = False
        kind, target = placement["target_kind"], placement["target"]
        if kind == "region":
            applies = target == body["region"]
        elif kind == "material":
            applies = target == body["material"]["sde"]
        elif kind == "window":
            x0, y0, x1, y1 = _m5_rectangle_window(windows[target])
            # A Transfinite tensor grid cannot localize a spacing change to
            # one tile without producing incompatible opposite-side counts.
            # Conservatively propagate each Rectangle window along the other
            # axis instead of missing a window that does not cross the device
            # centre.
            applies = (
                x0 - GEOM_TOL <= midpoint[0] <= x1 + GEOM_TOL
                if axis == "x" else
                y0 - GEOM_TOL <= midpoint[1] <= y1 + GEOM_TOL
            )
        if applies:
            spacing = min(spacing, float(sizes[placement["refinement"]][f"max_{axis}"]))
    return spacing


def _m5_transfinite_geometry(gmsh, device_ir: dict[str, Any], mesh_size: float,
                             diagonal: str) -> tuple[
                                 dict[int, dict[str, Any]], dict[str, Any], dict[str, Any]]:
    """Partition one rectangle on refinement-window cuts and transfinite each tile."""
    body = device_ir["geometry"]["bodies"][0]
    x0, y0, x1, y1 = _rectangle(body)
    control = device_ir["mesh_control"]
    windows = {item["name"]: item for item in control["windows"]}
    xcuts, ycuts = {x0, x1}, {y0, y1}
    for placement in control["refinement_placements"]:
        if placement["target_kind"] != "window":
            continue
        wx0, wy0, wx1, wy1 = _m5_rectangle_window(windows[placement["target"]])
        if wx1 < x0 - GEOM_TOL or wx0 > x1 + GEOM_TOL or \
                wy1 < y0 - GEOM_TOL or wy0 > y1 + GEOM_TOL:
            continue
        for value in (wx0, wx1):
            if x0 + GEOM_TOL < value < x1 - GEOM_TOL:
                xcuts.add(value)
        for value in (wy0, wy1):
            if y0 + GEOM_TOL < value < y1 - GEOM_TOL:
                ycuts.add(value)
    xs, ys = sorted(xcuts), sorted(ycuts)
    surfaces: list[int] = []
    for xa, xb in zip(xs, xs[1:]):
        for ya, yb in zip(ys, ys[1:]):
            surfaces.append(gmsh.model.occ.addRectangle(xa, ya, 0.0, xb - xa, yb - ya))
    gmsh.model.occ.fragment(
        [(2, surfaces[0])], [(2, tag) for tag in surfaces[1:]],
        removeObject=True, removeTool=True,
    ) if len(surfaces) > 1 else None
    gmsh.model.occ.synchronize()
    final_surfaces = sorted(tag for _dim, tag in gmsh.model.getEntities(2))
    surface_owner = {int(tag): body for tag in final_surfaces}

    x_intervals = []
    for xa, xb in zip(xs, xs[1:]):
        spacing = _m5_segment_spacing(
            device_ir, body, "x", ((xa + xb) * 0.5, (y0 + y1) * 0.5), mesh_size)
        x_intervals.append(max(1, int(math.ceil((xb - xa) / spacing - 1.0e-12))))
    y_intervals = []
    for ya, yb in zip(ys, ys[1:]):
        spacing = _m5_segment_spacing(
            device_ir, body, "y", ((x0 + x1) * 0.5, (ya + yb) * 0.5), mesh_size)
        y_intervals.append(max(1, int(math.ceil((yb - ya) / spacing - 1.0e-12))))

    curve_nodes: dict[int, int] = {}
    for surface in final_surfaces:
        curves = _boundary_curves(gmsh, [surface])
        if len(curves) != 4:
            raise GmshBackendError(
                "Transfinite rectangle tile does not have four boundary curves",
                {"stage": "mesh_strategy", "surface": surface, "curves": curves},
            )
        for curve in curves:
            xmin, ymin, _zmin, xmax, ymax, _zmax = gmsh.model.getBoundingBox(1, curve)
            horizontal = (ymax - ymin) <= 1.0e-6 * max(x1 - x0, y1 - y0, 1.0)
            if horizontal:
                midpoint = 0.5 * (xmin + xmax)
                index = next((i for i, (a, b) in enumerate(zip(xs, xs[1:]))
                              if a - GEOM_TOL <= midpoint <= b + GEOM_TOL), None)
                count = x_intervals[index] + 1 if index is not None else None
            else:
                midpoint = 0.5 * (ymin + ymax)
                index = next((i for i, (a, b) in enumerate(zip(ys, ys[1:]))
                              if a - GEOM_TOL <= midpoint <= b + GEOM_TOL), None)
                count = y_intervals[index] + 1 if index is not None else None
            if count is None:
                raise GmshBackendError(
                    "cannot map a Transfinite curve to an authored axis segment",
                    {"stage": "mesh_strategy", "curve": curve},
                )
            previous = curve_nodes.setdefault(curve, count)
            if previous != count:
                raise GmshBackendError(
                    "shared Transfinite curve received inconsistent node counts",
                    {"stage": "mesh_strategy", "curve": curve,
                     "counts": [previous, count]},
                )
        gmsh.model.mesh.setTransfiniteSurface(int(surface), diagonal)
    for curve, count in sorted(curve_nodes.items()):
        gmsh.model.mesh.setTransfiniteCurve(curve, count)

    expected_nodes = (sum(x_intervals) + 1) * (sum(y_intervals) + 1)
    expected_cells = 2 * sum(x_intervals) * sum(y_intervals)
    contract = {
        "schema": "vela.gmsh_transfinite_mesh_contract.v1",
        "requested_strategy": "transfinite-rect", "actual_strategy": "transfinite-rect",
        "diagonal": diagonal, "recombine": False, "element_order": 1,
        "x_coordinates": xs, "y_coordinates": ys,
        "x_intervals": x_intervals, "y_intervals": y_intervals,
        "expected_node_count": expected_nodes, "expected_cell_count": expected_cells,
        "surface_count": len(final_surfaces),
    }
    geometry_audit = {
        "algorithm": "partitioned-axis-aligned-transfinite-rectangle",
        "surface_owners": {str(tag): body["id"] for tag in final_surfaces},
        "region_areas_um2": {body["region"]: (x1 - x0) * (y1 - y0)},
        "partition": {"x_coordinates": xs, "y_coordinates": ys},
    }
    return surface_owner, geometry_audit, contract


def _qualify_m5_transfinite(bundle: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    coords = {item["vela_id"]: (item["x"], item["y"]) for item in bundle["nodes"]}
    right_errors = []
    maximum_angle = 0.0
    for triangle in bundle["triangles"]:
        angles = _angles([coords[node] for node in triangle["node_ids"]])
        maximum_angle = max(maximum_angle, max(angles))
        right_errors.append(min(abs(value - 90.0) for value in angles))
    failures = []
    if len(bundle["nodes"]) != contract["expected_node_count"]:
        failures.append("node count differs from the Transfinite contract")
    if len(bundle["triangles"]) != contract["expected_cell_count"]:
        failures.append("cell count differs from the Transfinite contract")
    if maximum_angle > 90.000001:
        failures.append("mesh contains an obtuse triangle")
    if right_errors and max(right_errors) > 1.0e-6:
        failures.append("not every Tri3 element has a right angle")
    return {
        "status": "passed" if not failures else "failed",
        "strategy": "transfinite-rect", "diagonal": contract["diagonal"],
        "expected_node_count": contract["expected_node_count"],
        "actual_node_count": len(bundle["nodes"]),
        "expected_cell_count": contract["expected_cell_count"],
        "actual_cell_count": len(bundle["triangles"]),
        "max_angle_degrees": maximum_angle,
        "max_right_angle_error_degrees": max(right_errors, default=0.0),
        "all_triangles_right": not right_errors or max(right_errors) <= 1.0e-6,
        "non_obtuse": maximum_angle <= 90.000001,
        "failures": failures,
    }


def _m5_outer_side(gmsh, curve: int,
                   bounds: tuple[float, float, float, float]) -> str | None:
    x0, y0, x1, y1 = bounds
    xmin, ymin, _zmin, xmax, ymax, _zmax = gmsh.model.getBoundingBox(1, curve)
    tolerance = 1.0e-6 * max(x1 - x0, y1 - y0, 1.0)
    if xmax - xmin <= tolerance:
        coordinate = 0.5 * (xmin + xmax)
        if abs(coordinate - x0) <= tolerance:
            return "left"
        if abs(coordinate - x1) <= tolerance:
            return "right"
    if ymax - ymin <= tolerance:
        coordinate = 0.5 * (ymin + ymax)
        if abs(coordinate - y0) <= tolerance:
            return "bottom"
        if abs(coordinate - y1) <= tolerance:
            return "top"
    return None


def _build_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float,
                  capability_level: str, mesh_strategy: str = "unstructured",
                  transfinite_diagonal: str = "AlternateRight",
                  maxtransdiff_max_iterations: int = 12, maxtransdiff_max_nodes: int = 200000) -> tuple[
                      dict[str, Any], dict[str, Any]]:
    """Build OCC geometry, physical groups, Tri3 mesh, and nodal fields."""
    if capability_level == "M1":
        validate_m1_device_ir(device_ir)
    elif capability_level == "M2":
        validate_m2_device_ir(device_ir)
    elif capability_level == "M3":
        validate_m3_device_ir(device_ir)
    elif capability_level == "M5":
        validate_m5_device_ir(device_ir)
    else:
        raise GmshBackendError(f"unknown capability level {capability_level!r}")
    if not math.isfinite(mesh_size) or mesh_size <= 0.0:
        raise GmshBackendError("mesh_size must be positive and finite")
    if mesh_strategy not in M5_MESH_STRATEGIES:
        raise GmshBackendError(f"unknown M5 mesh strategy {mesh_strategy!r}")
    if transfinite_diagonal not in M5_TRANSFINITE_DIAGONALS:
        raise GmshBackendError(
            f"unknown Transfinite diagonal arrangement {transfinite_diagonal!r}")
    eligibility = _m5_transfinite_eligibility(device_ir) if capability_level == "M5" else None
    actual_strategy = mesh_strategy
    if capability_level == "M5":
        if mesh_strategy == "auto":
            actual_strategy = "transfinite-rect" if eligibility["eligible"] else "unstructured"
        elif mesh_strategy == "transfinite-rect" and not eligibility["eligible"]:
            raise GmshBackendError(
                "explicit transfinite-rect strategy is not eligible for this Device IR",
                {"stage": "mesh_strategy", "strategy": mesh_strategy,
                 "eligibility": eligibility},
            )
    else:
        actual_strategy = "unstructured"
    anisotropic = capability_level == "M3" or (
        capability_level == "M5" and actual_strategy == "unstructured")
    options = configure_gmsh(gmsh, mesh_size, anisotropic=anisotropic)
    gmsh.model.add(f"vela_sentaurus_gmsh_{capability_level.lower()}")

    bodies = sorted(device_ir["geometry"]["bodies"], key=lambda item: item["creation_order"])
    transfinite_contract = None
    if capability_level == "M5" and actual_strategy == "transfinite-rect":
        surface_owner, geometry_audit, transfinite_contract = _m5_transfinite_geometry(
            gmsh, device_ir, mesh_size, transfinite_diagonal)
    elif capability_level == "M1":
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
    interface_curves = (_interface_curve_groups(gmsh, surface_owner)
                        if capability_level in {"M2", "M3", "M5"} else {})
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
            if assignment['kind']=='body_boundary':
                selected_id=assignment['entity_ids'][0]
                removed_body=next(b for b in bodies if b['creation_order']==selected_id)
                x0,y0,x1,y1=_rectangle(removed_body)
                for curve in _boundary_curves(gmsh,surface_owner):
                    lo,hi=gmsh.model.getParametrizationBounds(1,curve)
                    value=gmsh.model.getValue(1,curve,[(float(lo[0])+float(hi[0]))/2])
                    x,y=map(float,value[:2])
                    if not (x0-1e-8<=x<=x1+1e-8 and y0-1e-8<=y<=y1+1e-8 and
                            min(abs(x-x0),abs(x-x1),abs(y-y0),abs(y-y1))<=1e-8):continue
                    up,_=gmsh.model.getAdjacencies(1,curve)
                    adjacent=[int(t) for t in up if int(t) in surface_owner]
                    if len(adjacent)!=1:continue
                    old=curve_contact.get(curve)
                    if old is not None and old!=contact['name']:raise GmshBackendError('conflicting body boundary contact')
                    curve_contact[curve]=contact['name'];contact_curves[contact['name']].add(curve)
                    contact_owner[contact['name']].add(surface_owner[adjacent[0]]['region'])
                continue
            x, y = map(float, assignment["position"])
            candidates: list[tuple[int, bool]] = []
            for curve in _boundary_curves(gmsh, surface_owner):
                contains, endpoint = _curve_contains_pick(gmsh, curve, x, y)
                if contains:
                    candidates.append((curve, endpoint))
            candidates = sorted(set(candidates))
            interior = sorted({curve for curve, endpoint in candidates if not endpoint})
            selected_curves: list[int] = []
            if transfinite_contract is not None:
                bounds = _rectangle(bodies[0])
                candidate_sides = {
                    side for curve, _endpoint in candidates
                    if (side := _m5_outer_side(gmsh, curve, bounds)) is not None
                }
                if len(candidate_sides) == 1:
                    selected_side = next(iter(candidate_sides))
                    selected_curves = [
                        curve for curve in _boundary_curves(gmsh, surface_owner)
                        if _m5_outer_side(gmsh, curve, bounds) == selected_side
                    ]
            elif len(interior) == 1:
                selected_curves = [interior[0]]
            if not selected_curves:
                raise GmshBackendError(
                    f"contact {contact['name']!r} pick point does not select one curve interior",
                    {"stage": "contact", "contact": contact["name"],
                     "pick_point": [x, y], "candidates": candidates,
                     "source_location": assignment["provenance"].get(
                         "source_location")},
                )
            for curve in selected_curves:
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

    refinement_field_dag = None
    if capability_level == "M3" or (
            capability_level == "M5" and actual_strategy == "unstructured"):
        refinement_field_dag = _build_m3_refinement_fields(
            gmsh, device_ir, surface_owner, mesh_size)
    else:
        gmsh.model.mesh.setSize(gmsh.model.getEntities(0), mesh_size)
    gmsh.model.mesh.generate(2)
    adaptive_report = None
    if refinement_field_dag is not None:
        from sentaurus_adaptive import remesh
        adaptive_report = remesh(gmsh, device_ir, refinement_field_dag, mesh_size,
                                 maxtransdiff_max_iterations, maxtransdiff_max_nodes, surface_owner)
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

    general_diagnostics = {}
    try:
        fields = (_evaluate_constant_fields(device_ir, nodes, triangles, regions, materials)
                  if capability_level == "M1" else evaluate_m2_fields(
                      device_ir, nodes, triangles, regions, general_diagnostics))
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
            "capability": ({
                "M1": "strictly-separated-rectangles-constant-doping",
                "M2": "aba-bab-insert-contact-ordered-replacement-gaussian",
                "M3": "m2-core-anisotropic-refinement-maxlenint-pn2d",
                "M5": "m3-core-selectable-transfinite-rectangle-or-unstructured",
            }[capability_level]),
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
            **({'inactive_region_references': [
                {'kind': kind, 'name': p['name'], 'target': p['target'],
                 'reason': 'authored region removed before meshing'}
                for kind, placements in [('doping',device_ir['doping']['placements']),
                                         ('refinement',device_ir['mesh_control']['refinement_placements'])]
                for p in placements if p['target_kind']=='region' and
                any(b['region']==p['target'] and b.get('removed') for b in bodies)]}
               if any(b.get('removed') for b in bodies) else {}),
            "refinement_field_dag": refinement_field_dag,
            "maxtransdiff": adaptive_report,
            "general_evaluation": general_diagnostics,
            "mesh_strategy": {
                "requested": mesh_strategy, "actual": actual_strategy,
                "transfinite_diagonal": transfinite_diagonal,
                "eligibility": eligibility,
            },
            "transfinite_contract": transfinite_contract,
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
    qualification["maxtransdiff"] = adaptive_report
    qualification["geometry"] = geometry_audit
    qualification["mesh_strategy"] = actual_strategy
    if transfinite_contract is not None:
        qualification["node_volume_policy"] = "mixed_voronoi"
        qualification["require_non_obtuse"] = True
        qualification["non_obtuse"] = qualification["max_angle_degrees"] <= 90.000001
        qualification["transfinite"] = _qualify_m5_transfinite(bundle, transfinite_contract)
        if qualification["transfinite"]["status"] != "passed":
            raise GmshBackendError(
                "M5 Transfinite right-triangle mesh contract failed",
                {"stage": "mesh_quality", "transfinite": qualification["transfinite"]},
            )
    if refinement_field_dag is not None:
        qualification["refinement"] = _qualify_m3_refinement(
            gmsh, device_ir, nodes, triangles, refinement_field_dag)
        if qualification["refinement"]["status"] != "passed":
            raise GmshBackendError(
                "M3 refinement mesh quality contract failed",
                {"stage": "mesh_quality",
                 "refinement": qualification["refinement"]},
            )
    return bundle, qualification


def build_m1_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float
                    ) -> tuple[dict[str, Any], dict[str, Any]]:
    return _build_bundle(gmsh, device_ir, mesh_size, "M1")


def build_m2_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float
                    ) -> tuple[dict[str, Any], dict[str, Any]]:
    return _build_bundle(gmsh, device_ir, mesh_size, "M2")


def build_m3_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float,
                    maxtransdiff_max_iterations: int = 12, maxtransdiff_max_nodes: int = 200000
                    ) -> tuple[dict[str, Any], dict[str, Any]]:
    return _build_bundle(gmsh, device_ir, mesh_size, "M3",
                         maxtransdiff_max_iterations=maxtransdiff_max_iterations,
                         maxtransdiff_max_nodes=maxtransdiff_max_nodes)


def build_m5_bundle(gmsh, device_ir: dict[str, Any], mesh_size: float,
                    mesh_strategy: str = "auto",
                    transfinite_diagonal: str = "AlternateRight",
                    maxtransdiff_max_iterations: int = 12, maxtransdiff_max_nodes: int = 200000
                    ) -> tuple[dict[str, Any], dict[str, Any]]:
    return _build_bundle(
        gmsh, device_ir, mesh_size, "M5", mesh_strategy, transfinite_diagonal,
        maxtransdiff_max_iterations, maxtransdiff_max_nodes)


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
