#!/usr/bin/env python3
"""Device IR -> Vela ``mesh.json`` + ``doping.csv``.

The generator deliberately targets the file contract that ``JsonMeshReader``
and ``DCSweep`` already implement. No third mesh format is introduced.

Geometry model
--------------
Regions are axis-aligned rectangles evaluated in creation order, matching SDE's
default boolean behaviour where a later body claims the overlap. A tensor-product
grid is built from the union of all region/window edges refined by the applicable
``sdedr:define-refinement-size`` extents, and each grid rectangle is split into
two right triangles. Right triangles are non-obtuse by construction, which is a
precondition for the ``mixed_voronoi`` node-volume policy used by BV decks.

Contacts
--------
Each ``(find-edge-id (position ...))`` pick point is snapped to the outer
boundary segment that contains it, and then expanded to every mesh node lying on
that segment -- including nodes introduced by refinement. Pick points that fall
on a corner, or on a T-junction between two collinear boundary segments, are
resolved by explicit owner rules rather than by floating-point luck.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence


# Geometric tolerance in micrometres. SDE coordinates are authored in um and the
# fixtures use 1e-3 um features, so 1e-9 um is far below any meaningful length
# while still absorbing binary rounding from the expression evaluator.
GEOM_TOL = 1.0e-9

# Right triangles have a 90-degree angle; BoxGeometryBuilder applies its own
# tolerance on top of this when require_non_obtuse is set.
MAX_ANGLE_DEGREES = 90.0 + 1.0e-6


class MeshQualificationError(ValueError):
    """Raised when a generated mesh fails an acceptance gate.

    Carries a machine-readable report so callers can persist it verbatim.
    """

    def __init__(self, message: str, report: dict[str, Any]):
        super().__init__(message)
        self.report = report


@dataclass(frozen=True)
class Rect:
    x0: float
    y0: float
    x1: float
    y1: float

    @staticmethod
    def from_ir(entry: dict[str, Any]) -> "Rect":
        lower = entry["lower_left"]
        upper = entry["upper_right"]
        return Rect(float(lower[0]), float(lower[1]), float(upper[0]), float(upper[1]))

    def contains_point(self, x: float, y: float, tol: float = GEOM_TOL) -> bool:
        return (self.x0 - tol <= x <= self.x1 + tol
                and self.y0 - tol <= y <= self.y1 + tol)

    def contains_center(self, x: float, y: float) -> bool:
        return self.x0 <= x <= self.x1 and self.y0 <= y <= self.y1


@dataclass(frozen=True)
class BoundarySegment:
    """An outer boundary segment of the device outline."""

    axis: str          # "x" for a vertical segment (constant x), "y" for horizontal
    coordinate: float  # the constant coordinate
    start: float       # inclusive lower bound along the free axis
    end: float         # inclusive upper bound along the free axis

    def contains(self, x: float, y: float, tol: float = GEOM_TOL) -> bool:
        if self.axis == "x":
            return abs(x - self.coordinate) <= tol and self.start - tol <= y <= self.end + tol
        return abs(y - self.coordinate) <= tol and self.start - tol <= x <= self.end + tol

    def free_value(self, x: float, y: float) -> float:
        return y if self.axis == "x" else x

    def length(self) -> float:
        return self.end - self.start

    def label(self) -> str:
        return (f"{self.axis}={self.coordinate:.12g} "
                f"[{self.start:.12g}, {self.end:.12g}]")


def _dedupe_sorted(values: Iterable[float], tol: float = GEOM_TOL) -> list[float]:
    ordered = sorted(values)
    result: list[float] = []
    for value in ordered:
        if not result or abs(value - result[-1]) > tol:
            result.append(value)
        else:
            # Keep the first representative so that authored coordinates are not
            # perturbed by nearly-identical duplicates.
            continue
    return result


def _subdivide(lo: float, hi: float, max_step: float) -> list[float]:
    """Return interior cut positions splitting [lo, hi] into <= max_step pieces."""
    span = hi - lo
    if max_step <= 0.0 or span <= max_step + GEOM_TOL:
        return []
    count = max(2, int(math.ceil(span / max_step - GEOM_TOL)))
    step = span / count
    return [lo + step * index for index in range(1, count)]


@dataclass
class GridSpec:
    xs: list[float]
    ys: list[float]


def _region_rects(device_ir: dict[str, Any]) -> list[tuple[str, Rect, str]]:
    material_by_region = {
        item["region"]: item["vela_material"] for item in device_ir["materials"]
    }
    regions: list[tuple[str, Rect, str]] = []
    for entry in device_ir["geometry"]["regions"]:
        name = entry["name"]
        if entry.get("shape") != "rectangle":
            raise MeshQualificationError(
                f"region '{name}' has unsupported shape '{entry.get('shape')}'",
                {"stage": "geometry", "region": name, "shape": entry.get("shape")},
            )
        regions.append((name, Rect.from_ir(entry), material_by_region[name]))
    return regions


def _windows_by_name(device_ir: dict[str, Any]) -> dict[str, Rect]:
    return {
        window["name"]: Rect.from_ir(window)
        for window in device_ir["mesh_control"]["windows"]
    }


def _refinement_steps(device_ir: dict[str, Any],
                      regions: dict[str, Rect],
                      windows: dict[str, Rect]) -> list[tuple[Rect, float, float]]:
    """Resolve refinement placements to (area, max_dx, max_dy) triples."""
    sizes = {entry["name"]: entry for entry in device_ir["mesh_control"]["refinement_sizes"]}
    steps: list[tuple[Rect, float, float]] = []
    for placement in device_ir["mesh_control"]["refinement_placements"]:
        size = sizes[placement["refinement"]]
        if placement["target_kind"] == "region":
            area = regions[placement["target"]]
        else:
            area = windows[placement["target"]]
        steps.append((area, float(size["max_x"]), float(size["max_y"])))
    return steps


def build_grid(regions: Sequence[tuple[str, Rect, str]],
               windows: dict[str, Rect],
               steps: Sequence[tuple[Rect, float, float]],
               fallback_divisions: int = 8) -> GridSpec:
    """Build the tensor grid honouring every region, window, and refinement."""
    xs: set[float] = set()
    ys: set[float] = set()
    for _name, rect, _material in regions:
        xs.update((rect.x0, rect.x1))
        ys.update((rect.y0, rect.y1))
    for rect in windows.values():
        xs.update((rect.x0, rect.x1))
        ys.update((rect.y0, rect.y1))
    for rect, _dx, _dy in steps:
        xs.update((rect.x0, rect.x1))
        ys.update((rect.y0, rect.y1))

    domain_x = (min(xs), max(xs))
    domain_y = (min(ys), max(ys))

    if not steps:
        # Without any refinement control, fall back to a uniform grid so the
        # generator still produces a usable mesh instead of a 1-cell device.
        span_x = domain_x[1] - domain_x[0]
        span_y = domain_y[1] - domain_y[0]
        steps = [(Rect(*domain_x[:1], domain_y[0], domain_x[1], domain_y[1]),
                  span_x / fallback_divisions,
                  span_y / fallback_divisions)]

    # Apply the finest applicable step to every base interval it covers.
    for _pass in range(2):
        base_x = _dedupe_sorted(xs)
        base_y = _dedupe_sorted(ys)
        added = False
        for rect, max_dx, max_dy in steps:
            for lo, hi in zip(base_x, base_x[1:]):
                if hi <= rect.x0 + GEOM_TOL or lo >= rect.x1 - GEOM_TOL:
                    continue
                cuts = _subdivide(max(lo, rect.x0), min(hi, rect.x1), max_dx)
                # Also honour the part of the interval outside the window so the
                # grid stays a proper tensor product.
                cuts.extend(_subdivide(lo, hi, max_dx) if rect.x0 <= lo and hi <= rect.x1 else [])
                for cut in cuts:
                    if lo + GEOM_TOL < cut < hi - GEOM_TOL:
                        xs.add(cut)
                        added = True
            for lo, hi in zip(base_y, base_y[1:]):
                if hi <= rect.y0 + GEOM_TOL or lo >= rect.y1 - GEOM_TOL:
                    continue
                cuts = _subdivide(max(lo, rect.y0), min(hi, rect.y1), max_dy)
                cuts.extend(_subdivide(lo, hi, max_dy) if rect.y0 <= lo and hi <= rect.y1 else [])
                for cut in cuts:
                    if lo + GEOM_TOL < cut < hi - GEOM_TOL:
                        ys.add(cut)
                        added = True
        if not added:
            break

    return GridSpec(xs=_dedupe_sorted(xs), ys=_dedupe_sorted(ys))


def outer_boundary_segments(regions: Sequence[tuple[str, Rect, str]],
                            grid: GridSpec) -> list[BoundarySegment]:
    """Derive the outer boundary of the occupied domain from the grid.

    A grid cell face is on the outer boundary when exactly one of the two cells
    sharing it is occupied by a region. Collinear faces are merged into maximal
    segments so that a pick point anywhere along an SDE edge selects the whole
    edge.
    """
    occupied: set[tuple[int, int]] = set()
    for ix, (x0, x1) in enumerate(zip(grid.xs, grid.xs[1:])):
        for iy, (y0, y1) in enumerate(zip(grid.ys, grid.ys[1:])):
            cx = 0.5 * (x0 + x1)
            cy = 0.5 * (y0 + y1)
            if any(rect.contains_center(cx, cy) for _n, rect, _m in regions):
                occupied.add((ix, iy))

    vertical: dict[float, list[tuple[float, float]]] = {}
    horizontal: dict[float, list[tuple[float, float]]] = {}
    for ix, iy in sorted(occupied):
        x0, x1 = grid.xs[ix], grid.xs[ix + 1]
        y0, y1 = grid.ys[iy], grid.ys[iy + 1]
        if (ix - 1, iy) not in occupied:
            vertical.setdefault(x0, []).append((y0, y1))
        if (ix + 1, iy) not in occupied:
            vertical.setdefault(x1, []).append((y0, y1))
        if (ix, iy - 1) not in occupied:
            horizontal.setdefault(y0, []).append((x0, x1))
        if (ix, iy + 1) not in occupied:
            horizontal.setdefault(y1, []).append((x0, x1))

    segments: list[BoundarySegment] = []
    for axis, table in (("x", vertical), ("y", horizontal)):
        for coordinate, spans in table.items():
            for start, end in _merge_spans(spans):
                segments.append(BoundarySegment(axis, coordinate, start, end))
    segments.sort(key=lambda seg: (seg.axis, seg.coordinate, seg.start))
    return segments


def _merge_spans(spans: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    merged: list[tuple[float, float]] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1] + GEOM_TOL:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def select_contact_segment(pick_x: float,
                           pick_y: float,
                           segments: Sequence[BoundarySegment]) -> BoundarySegment:
    """Resolve an SDE pick point to exactly one outer boundary segment.

    Owner rules, applied in order:

    1. Only segments that geometrically contain the pick point are candidates.
    2. A pick point strictly inside a segment (a T-junction against a
       perpendicular segment that merely touches its endpoint) selects the
       segment that contains it in its interior. This is the case the previous
       implementation got wrong: an endpoint-touching segment used to be able to
       win by ordering.
    3. A pick point on a shared endpoint (a corner, or a T-junction between two
       collinear segments) is ambiguous and is rejected, because either choice
       would silently produce a different contact.
    """
    candidates = [segment for segment in segments if segment.contains(pick_x, pick_y)]
    if not candidates:
        raise MeshQualificationError(
            f"contact pick point ({pick_x:.12g}, {pick_y:.12g}) is not on the "
            f"outer boundary",
            {
                "stage": "contact",
                "pick_point": [pick_x, pick_y],
                "reason": "no outer boundary segment contains the pick point",
                "boundary": [segment.label() for segment in segments],
            },
        )

    interior = [
        segment for segment in candidates
        if segment.start + GEOM_TOL < segment.free_value(pick_x, pick_y) < segment.end - GEOM_TOL
    ]
    if len(interior) == 1:
        return interior[0]
    if len(interior) > 1:
        raise MeshQualificationError(
            f"contact pick point ({pick_x:.12g}, {pick_y:.12g}) lies in the "
            f"interior of {len(interior)} boundary segments (non-manifold "
            f"boundary)",
            {
                "stage": "contact",
                "pick_point": [pick_x, pick_y],
                "reason": "non-manifold outer boundary",
                "candidates": [segment.label() for segment in interior],
            },
        )
    if len(candidates) == 1:
        return candidates[0]
    raise MeshQualificationError(
        f"contact pick point ({pick_x:.12g}, {pick_y:.12g}) sits on the shared "
        f"endpoint of {len(candidates)} boundary segments; move the pick point "
        f"into the interior of the intended edge",
        {
            "stage": "contact",
            "pick_point": [pick_x, pick_y],
            "reason": "ambiguous pick point on a segment endpoint",
            "candidates": [segment.label() for segment in candidates],
        },
    )


@dataclass
class GeneratedMesh:
    mesh: dict[str, Any]
    doping_rows: list[dict[str, Any]]
    qualification: dict[str, Any]
    node_coords: list[tuple[float, float]] = field(default_factory=list)


def _triangle_angles(a: tuple[float, float],
                     b: tuple[float, float],
                     c: tuple[float, float]) -> list[float]:
    def angle_at(p, q, r):
        ux, uy = q[0] - p[0], q[1] - p[1]
        vx, vy = r[0] - p[0], r[1] - p[1]
        nu = math.hypot(ux, uy)
        nv = math.hypot(vx, vy)
        if nu == 0.0 or nv == 0.0:
            return 180.0
        cosine = max(-1.0, min(1.0, (ux * vx + uy * vy) / (nu * nv)))
        return math.degrees(math.acos(cosine))

    return [angle_at(a, b, c), angle_at(b, c, a), angle_at(c, a, b)]


def evaluate_doping(device_ir: dict[str, Any],
                    regions: dict[str, Rect],
                    windows: dict[str, Rect],
                    x: float,
                    y: float) -> tuple[float, float]:
    """Evaluate donor/acceptor concentration at a node.

    Constant profiles are applied in placement order; a later placement replaces
    an earlier one where their targets overlap, matching SDE's
    last-definition-wins semantics.

    Only constant profiles are supported today; any analytic profile type is
    rejected fail-closed.
    """
    profiles = {profile["name"]: profile for profile in device_ir["doping"]["profiles"]}
    donors = 0.0
    acceptors = 0.0
    for placement in sorted(device_ir["doping"]["placements"],
                            key=lambda item: item["priority"]):
        profile = profiles[placement["profile"]]
        area = (regions[placement["target"]] if placement["target_kind"] == "region"
                else windows[placement["target"]])
        if not area.contains_point(x, y):
            continue
        if profile["type"] != "constant":
            raise MeshQualificationError(
                f"doping profile '{profile['name']}' has unsupported type "
                f"'{profile['type']}'",
                {"stage": "doping", "profile": profile["name"]},
            )
        value = float(profile["value_cm3"])
        if profile["carrier"] == "donor":
            donors = value
        else:
            acceptors = value
    return donors, acceptors


def build_mesh_and_doping(device_ir: dict[str, Any],
                          require_non_obtuse: bool = True,
                          node_volume_policy: str = "mixed_voronoi") -> GeneratedMesh:
    """Generate ``mesh.json`` content and node doping rows from a device IR."""
    if device_ir.get("schema") != "vela.sentaurus_device_ir.v1":
        raise MeshQualificationError(
            f"unsupported device IR schema {device_ir.get('schema')!r}",
            {"stage": "input", "schema": device_ir.get("schema")},
        )

    region_list = _region_rects(device_ir)
    if not region_list:
        raise MeshQualificationError(
            f"{device_ir.get('source', '<device ir>')}: no region was created by "
            f"this SDE file, so no mesh can be generated",
            {"stage": "input", "regions": []},
        )
    region_rects = {name: rect for name, rect, _material in region_list}
    windows = _windows_by_name(device_ir)
    steps = _refinement_steps(device_ir, region_rects, windows)
    grid = build_grid(region_list, windows, steps)

    node_index: dict[tuple[int, int], int] = {}
    nodes: list[dict[str, Any]] = []
    coords: list[tuple[float, float]] = []

    def node_id(ix: int, iy: int) -> int:
        key = (ix, iy)
        if key not in node_index:
            node_index[key] = len(nodes)
            nodes.append({"id": len(nodes), "x": grid.xs[ix], "y": grid.ys[iy]})
            coords.append((grid.xs[ix], grid.ys[iy]))
        return node_index[key]

    triangles: list[dict[str, Any]] = []
    region_cells: dict[str, list[int]] = {name: [] for name, _r, _m in region_list}

    for ix, (x0, x1) in enumerate(zip(grid.xs, grid.xs[1:])):
        for iy, (y0, y1) in enumerate(zip(grid.ys, grid.ys[1:])):
            cx = 0.5 * (x0 + x1)
            cy = 0.5 * (y0 + y1)
            owner: str | None = None
            for name, rect, _material in region_list:
                if rect.contains_center(cx, cy):
                    owner = name  # later regions win, matching creation order
            if owner is None:
                continue
            n00 = node_id(ix, iy)
            n10 = node_id(ix + 1, iy)
            n11 = node_id(ix + 1, iy + 1)
            n01 = node_id(ix, iy + 1)
            # Split along the diagonal so both halves are right triangles.
            for node_ids in ((n00, n10, n11), (n00, n11, n01)):
                cell_id = len(triangles)
                triangles.append({
                    "id": cell_id,
                    "region_id": 0,  # patched below once region ids are assigned
                    "node_ids": list(node_ids),
                })
                region_cells[owner].append(cell_id)

    if not triangles:
        raise MeshQualificationError(
            "generated mesh has no cells",
            {"stage": "mesh", "reason": "no grid cell is covered by a region"},
        )

    used_regions = [name for name, _rect, _material in region_list if region_cells[name]]
    empty_regions = [name for name, _rect, _material in region_list if not region_cells[name]]
    if empty_regions:
        raise MeshQualificationError(
            f"region(s) {', '.join(empty_regions)} received no cells; the "
            f"geometry is fully overlapped by a later region",
            {"stage": "geometry", "empty_regions": empty_regions},
        )

    material_by_region = {name: material for name, _rect, material in region_list}
    regions_json: list[dict[str, Any]] = []
    for region_id, name in enumerate(used_regions):
        for cell_id in region_cells[name]:
            triangles[cell_id]["region_id"] = region_id
        regions_json.append({
            "id": region_id,
            "name": name,
            "material": material_by_region[name],
            "cell_ids": sorted(region_cells[name]),
        })
    region_id_by_name = {entry["name"]: entry["id"] for entry in regions_json}

    qualification = _qualify(triangles, coords, require_non_obtuse, node_volume_policy)

    segments = outer_boundary_segments(region_list, grid)
    contacts_json: list[dict[str, Any]] = []
    for contact_id, contact in enumerate(device_ir["contacts"]):
        node_ids: list[int] = []
        chosen: list[BoundarySegment] = []
        if not contact["pick_points"]:
            raise MeshQualificationError(
                f"contact '{contact['name']}' has no "
                f"'sdegeo:define-2d-contact' pick point; Vela cannot derive its "
                f"boundary nodes",
                {"stage": "contact", "contact": contact["name"], "segments": []},
            )
        for point in contact["pick_points"]:
            segment = select_contact_segment(float(point["x"]), float(point["y"]), segments)
            chosen.append(segment)
            for (ix, iy), nid in node_index.items():
                nx, ny = grid.xs[ix], grid.ys[iy]
                if segment.contains(nx, ny):
                    node_ids.append(nid)
        node_ids = sorted(dict.fromkeys(node_ids))
        if len(node_ids) < 2:
            raise MeshQualificationError(
                f"contact '{contact['name']}' resolved to {len(node_ids)} node(s); "
                f"at least two boundary nodes are required",
                {
                    "stage": "contact",
                    "contact": contact["name"],
                    "segments": [segment.label() for segment in chosen],
                },
            )
        owner_region = _contact_region(node_ids, triangles, regions_json)
        contacts_json.append({
            "id": contact_id,
            "name": contact["name"],
            "region_id": owner_region,
            "node_ids": node_ids,
        })
        qualification.setdefault("contacts", []).append({
            "name": contact["name"],
            "segments": [segment.label() for segment in chosen],
            "node_count": len(node_ids),
            "region_id": owner_region,
        })

    if not contacts_json:
        raise MeshQualificationError(
            "device IR defines no contacts",
            {"stage": "contact", "reason": "no contact set with a pick point"},
        )

    overlapping = _overlapping_contacts(contacts_json)
    if overlapping:
        raise MeshQualificationError(
            f"contacts share boundary nodes: {overlapping}",
            {"stage": "contact", "reason": "overlapping contacts", "pairs": overlapping},
        )

    mesh = {
        "_comment": (
            "Generated by scripts/sentaurus_mesh_builder.py from "
            f"{device_ir.get('source', 'an SDE command file')}. Coordinates are "
            "in micrometres and require scaling.mode = unit_scaling."
        ),
        "nodes": nodes,
        "triangles": triangles,
        "regions": regions_json,
        "contacts": contacts_json,
    }

    doping_rows = []
    for nid, (x, y) in enumerate(coords):
        donors, acceptors = evaluate_doping(device_ir, region_rects, windows, x, y)
        doping_rows.append({
            "node_id": nid,
            "donors_cm3": donors,
            "acceptors_cm3": acceptors,
        })

    qualification["node_count"] = len(nodes)
    qualification["cell_count"] = len(triangles)
    qualification["region_ids"] = region_id_by_name
    return GeneratedMesh(
        mesh=mesh,
        doping_rows=doping_rows,
        qualification=qualification,
        node_coords=coords,
    )


def _contact_region(node_ids: Sequence[int],
                    triangles: Sequence[dict[str, Any]],
                    regions_json: Sequence[dict[str, Any]]) -> int:
    """Return the region that owns every cell touching the contact nodes."""
    wanted = set(node_ids)
    touching = {
        cell["region_id"] for cell in triangles
        if wanted & set(cell["node_ids"])
    }
    if len(touching) == 1:
        return next(iter(touching))
    # A contact spanning several regions is attached to the region with the most
    # incident cells, which keeps mesh validation happy while remaining
    # deterministic.
    counts: dict[int, int] = {}
    for cell in triangles:
        if wanted & set(cell["node_ids"]):
            counts[cell["region_id"]] = counts.get(cell["region_id"], 0) + 1
    return max(sorted(counts), key=lambda region_id: (counts[region_id], -region_id))


def _overlapping_contacts(contacts: Sequence[dict[str, Any]]) -> list[list[str]]:
    pairs: list[list[str]] = []
    for index, first in enumerate(contacts):
        for second in contacts[index + 1:]:
            if set(first["node_ids"]) & set(second["node_ids"]):
                pairs.append([first["name"], second["name"]])
    return pairs


def _qualify(triangles: Sequence[dict[str, Any]],
             coords: Sequence[tuple[float, float]],
             require_non_obtuse: bool,
             node_volume_policy: str) -> dict[str, Any]:
    if node_volume_policy not in {"mixed_voronoi", "barycentric"}:
        raise MeshQualificationError(
            f"unsupported node_volume_policy '{node_volume_policy}'",
            {"stage": "qualification", "node_volume_policy": node_volume_policy},
        )

    min_angle = 180.0
    max_angle = 0.0
    obtuse: list[int] = []
    for cell in triangles:
        a, b, c = (coords[nid] for nid in cell["node_ids"])
        angles = _triangle_angles(a, b, c)
        min_angle = min(min_angle, min(angles))
        max_angle = max(max_angle, max(angles))
        if max(angles) > MAX_ANGLE_DEGREES:
            obtuse.append(cell["id"])

    report = {
        "min_angle_degrees": min_angle,
        "max_angle_degrees": max_angle,
        "obtuse_cells": obtuse,
        "require_non_obtuse": require_non_obtuse,
        "node_volume_policy": node_volume_policy,
    }
    if require_non_obtuse and obtuse:
        raise MeshQualificationError(
            f"{len(obtuse)} generated cell(s) exceed the non-obtuse limit "
            f"(max angle {max_angle:.6g} degrees)",
            {"stage": "qualification", **report},
        )
    return report


def write_mesh_json(path: Path, mesh: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(mesh, indent=2) + "\n", encoding="utf-8")


def write_doping_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "donors_cm3", "acceptors_cm3"])
        for row in rows:
            writer.writerow([
                row["node_id"],
                repr(float(row["donors_cm3"])),
                repr(float(row["acceptors_cm3"])),
            ])
