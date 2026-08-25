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
class WindowGeometry:
    name: str
    shape: str
    points: tuple[tuple[float, float], ...]

    @property
    def bounds(self) -> Rect:
        xs = [point[0] for point in self.points]
        ys = [point[1] for point in self.points]
        return Rect(min(xs), min(ys), max(xs), max(ys))

    def contains_point(self, x: float, y: float, tol: float = GEOM_TOL) -> bool:
        if self.shape == "rectangle":
            return self.bounds.contains_point(x, y, tol)
        if self.shape == "line":
            (x0, y0), (x1, y1) = self.points
            cross = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0)
            return abs(cross) <= tol * max(1.0, math.hypot(x1 - x0, y1 - y0)) and \
                min(x0, x1) - tol <= x <= max(x0, x1) + tol and \
                min(y0, y1) - tol <= y <= max(y0, y1) + tol
        # Boundary-inclusive ray casting for polygon Ref/Eval windows.
        inside = False
        for first, second in zip(self.points, self.points[1:] + self.points[:1]):
            x0, y0 = first
            x1, y1 = second
            cross = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0)
            if abs(cross) <= tol and min(x0, x1) - tol <= x <= max(x0, x1) + tol \
                    and min(y0, y1) - tol <= y <= max(y0, y1) + tol:
                return True
            if (y0 > y) != (y1 > y):
                hit_x = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
                if hit_x >= x - tol:
                    inside = not inside
        return inside


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


def _windows_by_name(device_ir: dict[str, Any]) -> dict[str, WindowGeometry]:
    windows: dict[str, WindowGeometry] = {}
    for entry in device_ir["mesh_control"]["windows"]:
        shape = entry["shape"]
        if shape == "rectangle":
            lower = entry["lower_left"]
            upper = entry["upper_right"]
            points = ((float(lower[0]), float(lower[1])),
                      (float(upper[0]), float(upper[1])))
        elif shape == "line":
            points = (tuple(map(float, entry["start"])), tuple(map(float, entry["end"])))
        elif shape == "polygon":
            points = tuple(tuple(map(float, point)) for point in entry["points"])
        else:
            raise MeshQualificationError(
                f"unsupported window shape '{shape}'",
                {"stage": "mesh_control", "window": entry["name"], "shape": shape})
        windows[entry["name"]] = WindowGeometry(entry["name"], shape, points)
    return windows


def _refinement_steps(device_ir: dict[str, Any],
                      regions: dict[str, Rect],
                      windows: dict[str, WindowGeometry]) -> list[tuple[Rect, float, float]]:
    """Resolve refinement placements to (area, max_dx, max_dy) triples."""
    sizes = {entry["name"]: entry for entry in device_ir["mesh_control"]["refinement_sizes"]}
    material_by_region = {
        entry["region"]: entry["sde_material"] for entry in device_ir["materials"]
    }
    functions: dict[str, list[dict[str, Any]]] = {}
    for function in device_ir["mesh_control"].get("refinement_functions", []):
        functions.setdefault(function["refinement"], []).append(function)
    steps: list[tuple[Rect, float, float]] = []
    for placement in device_ir["mesh_control"]["refinement_placements"]:
        size = sizes[placement["refinement"]]
        max_dx = float(size["max_x"])
        max_dy = float(size["max_y"])
        for function in functions.get(placement["refinement"], []):
            if function["function"] in {"MaxTransDiff", "MaxGradient"}:
                # Conservative static realization: refine the full placement to
                # the documented stopping target. This never treats min spacing
                # as a hard lower bound and cannot miss a requested transition.
                max_dx = min(max_dx, float(size.get("min_x", max_dx)))
                max_dy = min(max_dy, float(size.get("min_y", max_dy)))
            elif function["function"] == "MaxLenInt":
                interface_length = float(function["arguments"][2])
                max_dx = min(max_dx, interface_length)
                max_dy = min(max_dy, interface_length)
        if placement["target_kind"] == "region":
            areas = [regions[placement["target"]]]
        elif placement["target_kind"] == "material":
            areas = [rect for name, rect in regions.items()
                     if material_by_region[name] == placement["target"]]
        else:
            areas = [windows[placement["target"]].bounds]
        steps.extend((area, max_dx, max_dy) for area in areas)
    return steps


def build_grid(regions: Sequence[tuple[str, Rect, str]],
               windows: dict[str, WindowGeometry],
               steps: Sequence[tuple[Rect, float, float]],
               inserted_vertices: Sequence[dict[str, Any]] = (),
               fallback_divisions: int = 8) -> GridSpec:
    """Build the tensor grid honouring every region, window, and refinement."""
    xs: set[float] = set()
    ys: set[float] = set()
    for _name, rect, _material in regions:
        xs.update((rect.x0, rect.x1))
        ys.update((rect.y0, rect.y1))
    for window in windows.values():
        rect = window.bounds
        xs.update((rect.x0, rect.x1))
        ys.update((rect.y0, rect.y1))
        for x, y in window.points:
            xs.add(x)
            ys.add(y)
    for vertex in inserted_vertices:
        xs.add(float(vertex["x"]))
        ys.add(float(vertex["y"]))
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
                            grid: GridSpec,
                            inserted_vertices: Sequence[dict[str, Any]] = ()) \
        -> list[BoundarySegment]:
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
    split_segments: list[BoundarySegment] = []
    unmatched = set(range(len(inserted_vertices)))
    for segment in segments:
        cuts = [segment.start, segment.end]
        for index, vertex in enumerate(inserted_vertices):
            x, y = float(vertex["x"]), float(vertex["y"])
            if segment.contains(x, y):
                free = segment.free_value(x, y)
                if segment.start + GEOM_TOL < free < segment.end - GEOM_TOL:
                    cuts.append(free)
                unmatched.discard(index)
        ordered = _dedupe_sorted(cuts)
        split_segments.extend(
            BoundarySegment(segment.axis, segment.coordinate, lo, hi)
            for lo, hi in zip(ordered, ordered[1:])
        )
    if unmatched:
        points = [inserted_vertices[index] for index in sorted(unmatched)]
        raise MeshQualificationError(
            f"inserted vertex is not on a unique outer 2-D straight edge: {points}",
            {"stage": "geometry", "reason": "invalid inserted vertex", "vertices": points})
    split_segments.sort(key=lambda seg: (seg.axis, seg.coordinate, seg.start))
    return split_segments


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
                    windows: dict[str, WindowGeometry],
                    x: float,
                    y: float) -> tuple[float, float]:
    """Evaluate donor/acceptor concentration at a node.

    Profiles are applied in source order.  NoReplace adds, LocalReplace clears
    only the same species, and Replace clears every previously active species.
    Donor and acceptor species remain separate until the final aggregation.
    """
    profiles = {profile["name"]: profile for profile in device_ir["doping"]["profiles"]}
    material_by_region = {
        item["region"]: item["sde_material"] for item in device_ir["materials"]
    }
    owner = _owner_at(device_ir, regions, x, y)
    species_values: dict[str, tuple[str, float]] = {}
    for placement in sorted(device_ir["doping"]["placements"],
                            key=lambda item: item["priority"]):
        profile = profiles[placement["profile"]]
        target_kind = placement["target_kind"]
        active = False
        constant_weight = 1.0
        if profile["type"] == "gaussian":
            active = _analytic_restriction_active(
                placement, owner, material_by_region, windows, x, y)
        else:
            constant_weight = _constant_placement_weight(
                placement, owner, material_by_region, regions, windows, x, y)
            active = constant_weight > 0.0
        if not active:
            continue
        if profile["type"] == "constant":
            value = float(profile["value_cm3"]) * constant_weight
        elif profile["type"] == "gaussian":
            value = _evaluate_gaussian(profile, placement, windows, x, y)
        else:
            raise MeshQualificationError(
                f"doping profile '{profile['name']}' has unsupported type "
                f"'{profile['type']}'",
                {"stage": "doping", "profile": profile["name"]},
            )
        replacement = placement.get("replacement", "NoReplace")
        species = profile["species"]
        if replacement == "Replace":
            species_values.clear()
        elif replacement == "LocalReplace":
            species_values.pop(species, None)
        elif replacement != "NoReplace":
            raise MeshQualificationError(
                f"unsupported replacement mode '{replacement}'",
                {"stage": "doping", "placement": placement["name"]})
        previous = species_values.get(species, (profile["carrier"], 0.0))[1]
        species_values[species] = (profile["carrier"], previous + value)
    donors = sum(value for carrier, value in species_values.values() if carrier == "donor")
    acceptors = sum(value for carrier, value in species_values.values() if carrier == "acceptor")
    return donors, acceptors


def _distance_to_rect(rect: Rect, x: float, y: float) -> float:
    dx = max(rect.x0 - x, 0.0, x - rect.x1)
    dy = max(rect.y0 - y, 0.0, y - rect.y1)
    return math.hypot(dx, dy)


def _distance_to_segment(points: tuple[tuple[float, float], ...],
                         x: float, y: float) -> float:
    (x0, y0), (x1, y1) = points[0], points[1]
    dx, dy = x1 - x0, y1 - y0
    scale = dx * dx + dy * dy
    fraction = 0.0 if scale == 0.0 else max(
        0.0, min(1.0, ((x - x0) * dx + (y - y0) * dy) / scale))
    return math.hypot(x - (x0 + fraction * dx), y - (y0 + fraction * dy))


def _constant_placement_weight(placement: dict[str, Any], owner: str | None,
                               materials: dict[str, str], regions: dict[str, Rect],
                               windows: dict[str, WindowGeometry],
                               x: float, y: float) -> float:
    kind = placement["target_kind"]
    target = placement["target"]
    inside = False
    distances: list[float] = []
    if kind == "region":
        inside = owner == target
        distances = [_distance_to_rect(regions[target], x, y)]
    elif kind == "material":
        inside = owner is not None and materials[owner] == target
        distances = [_distance_to_rect(regions[name], x, y)
                     for name, material in materials.items() if material == target]
    else:
        window = windows[target]
        inside = window.contains_point(x, y)
        if window.shape == "line":
            distances = [_distance_to_segment(window.points, x, y)]
        elif window.shape == "polygon":
            distances = [
                _distance_to_segment((first, second), x, y)
                for first, second in zip(
                    window.points, window.points[1:] + window.points[:1])
            ]
        else:
            distances = [_distance_to_rect(window.bounds, x, y)]
    if inside:
        return 1.0
    length = float(placement.get("decay_length", 0.0))
    if length <= 0.0 or not distances:
        return 0.0
    distance = min(distances)
    if placement.get("decay") == "gaussian":
        return math.exp(-(distance / length) ** 2)
    if placement.get("decay") == "erf":
        return 0.5 * math.erfc(distance / length)
    return 0.0


def _owner_at(device_ir: dict[str, Any], regions: dict[str, Rect],
              x: float, y: float) -> str | None:
    owner: str | None = None
    for entry in sorted(device_ir["geometry"]["regions"], key=lambda item: item["order"]):
        if not regions[entry["name"]].contains_point(x, y):
            continue
        mode = entry.get("boolean_mode", "ABA")
        if mode == "ABA" or owner is None:
            owner = entry["name"]
        elif mode != "BAB":
            raise MeshQualificationError(
                f"unsupported boolean mode '{mode}'",
                {"stage": "geometry", "region": entry["name"], "mode": mode})
    return owner


def _analytic_restriction_active(placement: dict[str, Any], owner: str | None,
                                 materials: dict[str, str],
                                 windows: dict[str, WindowGeometry],
                                 x: float, y: float) -> bool:
    kind = placement.get("restriction_kind")
    target = placement.get("restriction")
    if kind is None:
        return owner is not None
    if kind == "region":
        return owner == target
    if kind == "material":
        return owner is not None and materials[owner] == target
    return windows[target].contains_point(x, y)


def _evaluate_gaussian(profile: dict[str, Any], placement: dict[str, Any],
                       windows: dict[str, WindowGeometry], x: float, y: float) -> float:
    baseline = windows[placement["target"]]
    (x0, y0), (x1, y1) = baseline.points
    dx, dy = x1 - x0, y1 - y0
    length = math.hypot(dx, dy)
    tx, ty = dx / length, dy / length
    # T-2022.03-SP2 TDR oracle: Positive is the left normal from start to end.
    nx, ny = -ty, tx
    along = (x - x0) * tx + (y - y0) * ty
    normal = (x - x0) * nx + (y - y0) * ny
    direction = placement.get("direction", "Both")
    if direction == "Positive" and normal < -GEOM_TOL:
        return 0.0
    if direction == "Negative" and normal > GEOM_TOL:
        return 0.0
    if placement.get("evaluation") == "NoEval" and abs(normal) <= GEOM_TOL:
        return 0.0
    coordinate = abs(normal) if direction == "Both" else \
        normal if direction == "Positive" else -normal
    params = profile["parameters"]
    peak_pos = float(params["PeakPos"])
    if "StdDev" in params:
        sigma = abs(float(params["StdDev"]))
    elif "Length" in params:
        sigma = abs(float(params["Length"])) / math.sqrt(2.0)
    elif {"PeakVal", "ValueAtDepth", "Depth"} <= params.keys():
        peak = float(params["PeakVal"])
        depth_value = float(params["ValueAtDepth"])
        if not (0.0 < depth_value < peak):
            raise MeshQualificationError(
                "Gaussian ValueAtDepth must be between zero and PeakVal",
                {"stage": "doping", "profile": profile["name"]})
        sigma = abs(float(params["Depth"]) - peak_pos) / \
            math.sqrt(2.0 * math.log(peak / depth_value))
    else:
        raise MeshQualificationError(
            "Dose Gaussian requires StdDev or Length",
            {"stage": "doping", "profile": profile["name"]})
    if sigma <= 0.0:
        raise MeshQualificationError(
            "Gaussian standard deviation is zero",
            {"stage": "doping", "profile": profile["name"]})
    if "PeakVal" in params:
        peak_value = float(params["PeakVal"])
    else:
        peak_value = float(params["Dose"]) / (
            math.sqrt(2.0 * math.pi) * sigma * 1.0e-4)
    value = peak_value * math.exp(-0.5 * ((coordinate - peak_pos) / sigma) ** 2)
    lateral_distance = max(0.0, -along, along - length)
    if lateral_distance > 0.0:
        lateral_scale = params.get("StdDev", params.get("Length", sigma * math.sqrt(2.0)))
        lateral_scale = abs(float(params.get("Factor", 1.0))) * abs(float(lateral_scale))
        if lateral_scale <= 0.0:
            return 0.0
        if profile.get("lateral_type") == "erf":
            value *= 0.5 * math.erfc(lateral_distance / lateral_scale)
        else:
            value *= math.exp(-(lateral_distance / lateral_scale) ** 2)
    return value


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
    inserted_vertices = device_ir["geometry"].get("inserted_vertices", [])
    grid = build_grid(region_list, windows, steps, inserted_vertices)

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
            owner = _owner_at(device_ir, region_rects, cx, cy)
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
    qualification["removed_zero_area_regions"] = empty_regions

    segments = outer_boundary_segments(region_list, grid, inserted_vertices)
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
    names = [entry["name"] for entry in regions_json if entry["id"] in touching]
    raise MeshQualificationError(
        f"contact boundary touches multiple region owners: {names}",
        {"stage": "contact", "reason": "ambiguous region owner", "regions": names})


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


def qualify_mesh(mesh: dict[str, Any],
                 require_non_obtuse: bool = True,
                 node_volume_policy: str = "mixed_voronoi") -> dict[str, Any]:
    """Qualify an existing Vela mesh, including one imported from a TDR.

    The offline generator and the authoritative TDR importer must feed the
    same acceptance report shape to the E2E contract.  Keep this public helper
    limited to geometry/contact checks; doping provenance is recorded by the
    orchestration manifest.
    """
    nodes = list(mesh.get("nodes", []))
    triangles = list(mesh.get("triangles", []))
    ordered = sorted(nodes, key=lambda node: int(node["id"]))
    if [int(node["id"]) for node in ordered] != list(range(len(ordered))):
        raise MeshQualificationError(
            "existing mesh node IDs must be contiguous from zero",
            {"stage": "qualification", "reason": "non-contiguous node IDs"},
        )
    coords = [(float(node["x"]), float(node["y"])) for node in ordered]
    report = _qualify(
        triangles, coords, require_non_obtuse, node_volume_policy)
    contacts = list(mesh.get("contacts", []))
    empty = [str(contact.get("name", "")) for contact in contacts
             if len(contact.get("node_ids", [])) < 2]
    overlaps = _overlapping_contacts(contacts)
    if empty or overlaps:
        raise MeshQualificationError(
            "existing mesh has invalid contact topology",
            {
                "stage": "contact",
                "empty_or_singleton_contacts": empty,
                "overlapping_contacts": overlaps,
            },
        )
    report["contacts"] = [
        {
            "name": str(contact["name"]),
            "node_count": len(contact.get("node_ids", [])),
        }
        for contact in contacts
    ]
    report["node_count"] = len(nodes)
    report["cell_count"] = len(triangles)
    report["region_ids"] = {
        str(region["name"]): int(region["id"])
        for region in mesh.get("regions", [])
    }
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
