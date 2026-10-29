#!/usr/bin/env python3
"""Backend-independent M2 Sentaurus doping evaluation.

The evaluator consumes Device IR v2 and authored coordinates.  It deliberately
does not query Gmsh: mesh backends provide only the node coordinates, while the
ordered placement semantics remain shared and independently testable.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any


GEOM_TOL = 1.0e-9
LOG_DOPING_FLOOR_CM3 = 1.0
SEMICONDUCTOR_MATERIALS = {"Si", "PolySilicon"}


class DopingEvaluationError(ValueError):
    def __init__(self, message: str, details: dict[str, Any] | None = None):
        super().__init__(message)
        self.details = details or {}


@dataclass(frozen=True)
class Rect:
    x0: float
    y0: float
    x1: float
    y1: float

    def contains(self, x: float, y: float, tolerance: float = GEOM_TOL) -> bool:
        return (self.x0 - tolerance <= x <= self.x1 + tolerance and
                self.y0 - tolerance <= y <= self.y1 + tolerance)


@dataclass(frozen=True)
class Window:
    name: str
    shape: str
    points: tuple[tuple[float, float], ...]

    @property
    def bounds(self) -> Rect:
        xs = [point[0] for point in self.points]
        ys = [point[1] for point in self.points]
        return Rect(min(xs), min(ys), max(xs), max(ys))

    def contains(self, x: float, y: float) -> bool:
        if self.shape == "rectangle":
            return self.bounds.contains(x, y)
        if self.shape == "line":
            return _distance_to_segment(self.points[0], self.points[1], x, y) <= GEOM_TOL
        if self.shape != "polygon":
            raise DopingEvaluationError(
                f"unsupported window shape {self.shape!r}",
                {"stage": "doping", "window": self.name},
            )
        inside = False
        for first, second in zip(self.points, self.points[1:] + self.points[:1]):
            if _distance_to_segment(first, second, x, y) <= GEOM_TOL:
                return True
            x0, y0 = first
            x1, y1 = second
            if (y0 > y) != (y1 > y):
                hit_x = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
                if hit_x >= x - GEOM_TOL:
                    inside = not inside
        return inside


def _rectangles(device_ir: dict[str, Any]) -> dict[str, Rect]:
    result = {}
    for body in device_ir["geometry"]["bodies"]:
        geometry = body["geometry"]
        lower = geometry["lower_left"]
        upper = geometry["upper_right"]
        result[body["region"]] = Rect(
            float(lower[0]), float(lower[1]), float(upper[0]), float(upper[1]))
    return result


def _windows(device_ir: dict[str, Any]) -> dict[str, Window]:
    result = {}
    for item in device_ir["doping"]["windows"]:
        shape = item["shape"]
        if shape == "rectangle":
            lower, upper = item["lower_left"], item["upper_right"]
            points = ((float(lower[0]), float(lower[1])),
                      (float(upper[0]), float(upper[1])))
        elif shape == "line":
            points = (tuple(map(float, item["start"])), tuple(map(float, item["end"])))
        elif shape == "polygon":
            points = tuple(tuple(map(float, point)) for point in item["points"])
        else:
            raise DopingEvaluationError(
                f"unsupported window shape {shape!r}",
                {"stage": "doping", "window": item["name"]},
            )
        result[item["name"]] = Window(item["name"], shape, points)
    return result


def _owner_body(device_ir: dict[str, Any], rectangles: dict[str, Rect],
                x: float, y: float) -> dict[str, Any] | None:
    """Replay authored boolean ownership, including boundary points.

    This is the explicit M2 global-node boundary policy.  A later ABA body owns
    a point shared with an earlier body; BAB preserves the existing owner.
    """
    owner = None
    for body in sorted(device_ir["geometry"]["bodies"],
                       key=lambda item: item["creation_order"]):
        if not rectangles[body["region"]].contains(x, y):
            continue
        if owner is None or body["boolean_mode"] == "ABA":
            owner = body
        elif body["boolean_mode"] != "BAB":
            raise DopingEvaluationError(
                f"unsupported boolean mode {body['boolean_mode']!r}",
                {"stage": "geometry", "body": body["id"]},
            )
    return owner


def _distance_to_segment(first: tuple[float, float], second: tuple[float, float],
                         x: float, y: float) -> float:
    x0, y0 = first
    dx, dy = second[0] - x0, second[1] - y0
    scale = dx * dx + dy * dy
    fraction = 0.0 if scale == 0.0 else max(
        0.0, min(1.0, ((x - x0) * dx + (y - y0) * dy) / scale))
    return math.hypot(x - (x0 + fraction * dx), y - (y0 + fraction * dy))


def _constant_active(placement: dict[str, Any], owner: dict[str, Any] | None,
                     windows: dict[str, Window], x: float, y: float) -> bool:
    kind, target = placement["target_kind"], placement["target"]
    if kind == "region":
        return owner is not None and owner["region"] == target
    if kind == "material":
        return owner is not None and owner["material"]["sde"] == target
    return windows[target].contains(x, y)


def _restriction_active(placement: dict[str, Any], owner: dict[str, Any] | None,
                        windows: dict[str, Window], x: float, y: float) -> bool:
    kind, target = placement.get("restriction_kind"), placement.get("restriction")
    if kind is None:
        return owner is not None
    if kind == "region":
        return owner is not None and owner["region"] == target
    if kind == "material":
        return owner is not None and owner["material"]["sde"] == target
    return windows[target].contains(x, y)


def _gaussian_value(profile: dict[str, Any], placement: dict[str, Any],
                    windows: dict[str, Window], x: float, y: float) -> tuple[bool, float]:
    baseline = windows[placement["target"]]
    (x0, y0), (x1, y1) = baseline.points
    dx, dy = x1 - x0, y1 - y0
    length = math.hypot(dx, dy)
    if length <= GEOM_TOL:
        raise DopingEvaluationError(
            "Gaussian baseline has zero length",
            {"stage": "doping", "placement": placement["name"]},
        )
    tx, ty = dx / length, dy / length
    # T-2022.03-SP2 oracle: Positive is the left normal of the directed
    # baseline (for an east-going line this is the positive-y side).
    nx, ny = -ty, tx
    along = (x - x0) * tx + (y - y0) * ty
    normal = (x - x0) * nx + (y - y0) * ny
    direction = placement.get("direction", "Both")
    if direction == "Positive" and normal < -GEOM_TOL:
        return False, 0.0
    if direction == "Negative" and normal > GEOM_TOL:
        return False, 0.0
    if placement.get("evaluation") == "NoEval" and abs(normal) <= GEOM_TOL:
        return False, 0.0
    coordinate = abs(normal) if direction == "Both" else (
        normal if direction == "Positive" else -normal)
    params = profile["parameters"]
    peak_position = float(params["PeakPos"])
    if "StdDev" in params:
        sigma = abs(float(params["StdDev"]))
    elif "Length" in params:
        sigma = abs(float(params["Length"])) / math.sqrt(2.0)
    elif {"PeakVal", "ValueAtDepth", "Depth"} <= params.keys():
        peak = float(params["PeakVal"])
        depth_value = float(params["ValueAtDepth"])
        if not (0.0 < depth_value < peak):
            raise DopingEvaluationError(
                "Gaussian ValueAtDepth must be between zero and PeakVal",
                {"stage": "doping", "profile": profile["name"]},
            )
        sigma = abs(float(params["Depth"]) - peak_position) / math.sqrt(
            2.0 * math.log(peak / depth_value))
    else:
        raise DopingEvaluationError(
            "Dose Gaussian requires StdDev or Length",
            {"stage": "doping", "profile": profile["name"]},
        )
    if sigma <= 0.0:
        raise DopingEvaluationError(
            "Gaussian standard deviation is zero",
            {"stage": "doping", "profile": profile["name"]},
        )
    peak_value = float(params["PeakVal"]) if "PeakVal" in params else (
        float(params["Dose"]) / (math.sqrt(2.0 * math.pi) * sigma * 1.0e-4))
    value = peak_value * math.exp(-0.5 * ((coordinate - peak_position) / sigma) ** 2)
    lateral_distance = max(0.0, -along, along - length)
    if lateral_distance > 0.0:
        base_scale = params.get("StdDev", params.get("Length", sigma * math.sqrt(2.0)))
        lateral_scale = abs(float(params.get("Factor", 1.0))) * abs(float(base_scale))
        if lateral_scale <= 0.0:
            return False, 0.0
        if profile.get("lateral_type") == "erf":
            value *= 0.5 * math.erfc(lateral_distance / lateral_scale)
        else:
            value *= math.exp(-(lateral_distance / lateral_scale) ** 2)
    return True, value


def _evaluate_m2_doping_for_owner(
    device_ir: dict[str, Any], owner: dict[str, Any] | None,
    rectangles: dict[str, Rect], windows: dict[str, Window], x: float, y: float, diagnostics: dict | None = None,
) -> tuple[float, float, dict[str, float]]:
    if owner is None or owner["material"]["vela"] not in SEMICONDUCTOR_MATERIALS:
        return 0.0, 0.0, {}
    profiles = {item["name"]: item for item in device_ir["doping"]["profiles"]}
    values: dict[str, tuple[str, float]] = {}
    for placement in sorted(device_ir["doping"]["placements"],
                            key=lambda item: item["priority"]):
        profile = profiles[placement["profile"]]
        if profile["type"] == "constant":
            active = _constant_active(placement, owner, windows, x, y)
            value = float(profile["value_cm3"])
        elif profile["type"] == "gaussian":
            active = _restriction_active(placement, owner, windows, x, y)
            gaussian_active, value = _gaussian_value(profile, placement, windows, x, y)
            active = active and gaussian_active
        elif profile["type"] == "general":
            from sentaurus_general import evaluate
            active = _restriction_active(placement, owner, windows, x, y)
            if not active:
                continue
            try:
                value, fallback = evaluate(profile, x, y)
            except ValueError as error:
                raise DopingEvaluationError(str(error), {"stage": "doping", "profile": profile["name"]}) from error
            if diagnostics is not None:
                key = profile["name"]
                item = diagnostics.setdefault(key, {"evaluations": 0, "default_evaluations": 0})
                item["evaluations"] += 1
                item["default_evaluations"] += int(fallback)
        else:
            raise DopingEvaluationError(
                f"unsupported profile type {profile['type']!r}",
                {"stage": "doping", "profile": profile["name"]},
            )
        if not active:
            continue
        species = profile["species"]
        replacement = placement.get("replacement", "NoReplace")
        if replacement == "Replace":
            values.clear()
        elif replacement == "LocalReplace":
            values.pop(species, None)
        elif replacement != "NoReplace":
            raise DopingEvaluationError(
                f"unsupported replacement mode {replacement!r}",
                {"stage": "doping", "placement": placement["name"]},
            )
        previous = values.get(species, (profile["carrier"], 0.0))[1]
        values[species] = (profile["carrier"], previous + value)
    donors = sum(value for carrier, value in values.values() if carrier == "donor")
    acceptors = sum(value for carrier, value in values.values() if carrier == "acceptor")
    return donors, acceptors, {name: value for name, (_carrier, value) in values.items()}


def evaluate_m2_doping_at(device_ir: dict[str, Any], x: float, y: float
                          ) -> tuple[float, float, dict[str, float]]:
    """Evaluate one authored point using the deterministic boolean owner.

    Gmsh mesh nodes on a conformal interface need the separate region-local
    policy implemented by :func:`evaluate_m2_fields` below.
    """
    rectangles = _rectangles(device_ir)
    windows = _windows(device_ir)
    owner = _owner_body(device_ir, rectangles, x, y)
    return _evaluate_m2_doping_for_owner(
        device_ir, owner, rectangles, windows, x, y)


def _field(values: list[float], unit: str = "cm^-3") -> dict[str, Any]:
    return {"association": "node", "unit": unit, "values": values}


def evaluate_m2_fields(
    device_ir: dict[str, Any], nodes: list[dict[str, Any]],
    triangles: list[dict[str, Any]], regions: list[dict[str, Any]],
    diagnostics: dict | None = None,
) -> dict[str, Any]:
    """Evaluate fields using T-2022.03 region-local shared-node semantics."""
    profiles = {item["name"]: item for item in device_ir["doping"]["profiles"]}
    species_names = sorted({item["species"] for item in profiles.values()})
    species = {name: [0.0] * len(nodes) for name in species_names}
    rectangles = _rectangles(device_ir)
    windows = _windows(device_ir)
    body_by_region = {body["region"]: body for body in device_ir["geometry"]["bodies"]}
    region_name = {item["id"]: item["name"] for item in regions}
    incident_regions: dict[int, set[str]] = {node["vela_id"]: set() for node in nodes}
    for triangle in triangles:
        name = region_name[triangle["region_id"]]
        for node_id in triangle["node_ids"]:
            incident_regions[node_id].add(name)
    donors: list[float] = []
    acceptors: list[float] = []
    for node in nodes:
        donor = 0.0
        acceptor = 0.0
        active_species = {name: 0.0 for name in species_names}
        for region in sorted(incident_regions[node["vela_id"]]):
            owner = body_by_region.get(region)
            if owner is None:
                raise DopingEvaluationError(
                    f"mesh region {region!r} has no authored body",
                    {"stage": "doping", "region": region},
                )
            local_donor, local_acceptor, local_species = _evaluate_m2_doping_for_owner(
                device_ir, owner, rectangles, windows,
                float(node["x"]), float(node["y"]), diagnostics)
            donor += local_donor
            acceptor += local_acceptor
            for name, value in local_species.items():
                active_species[name] += value
        donors.append(donor)
        acceptors.append(acceptor)
        for name, value in active_species.items():
            species[name][node["vela_id"]] = value
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
