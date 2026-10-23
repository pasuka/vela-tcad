"""Read native GrdNumbering box evidence; never load it into a solver.

Measures are um^2 and coefficients are dimensionless for a 2-D mesh.  The
native local edge order is (0,1), (1,2), (2,0).  Bulk grid element IDs must be
matched to the original grid, not inferred from an internal element index.
"""
from __future__ import annotations

import math
import re


def read_box(text: str) -> dict:
    result = {}
    for name in ("Measure", "Coefficients"):
        match = re.search(r"\b" + name + r"\s*\{([^}]+)\}", text, re.I)
        if match is None:
            raise ValueError("missing native box section: " + name)
        rows = {}
        for line in match[1].splitlines():
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            if len(fields) < 3:
                raise ValueError("incomplete native box row")
            grid, internal, kind = map(int, fields[:3])
            if grid in rows:
                raise ValueError("duplicate native grid element")
            values = tuple(map(float, fields[3:]))
            if internal >= 0 and (kind != 2 or len(values) != 3):
                raise ValueError("only 2-D bulk triangles are supported")
            if any(not math.isfinite(x) or x < 0 for x in values):
                raise ValueError("nonfinite or negative native box value")
            rows[grid] = {"internal_id": internal, "kind": kind, "values": values}
        result[name] = rows
    if result["Measure"].keys() != result["Coefficients"].keys():
        raise ValueError("native measure/coefficient grid IDs differ")
    for grid, measure in result["Measure"].items():
        coeff = result["Coefficients"][grid]
        if (measure["internal_id"], measure["kind"]) != (coeff["internal_id"], coeff["kind"]):
            raise ValueError("native measure/coefficient topology differs")
    return result


def dielectric_entries(triangles: list, coefficients: dict, epsilon: dict) -> dict:
    """Assemble epsilon_r*kappa; absolute eps0 and voltage units are separate."""
    entries = {}
    seen = set()
    for cell in triangles:
        grid = cell["id"]
        if grid in seen:
            raise ValueError("duplicate mesh triangle")
        seen.add(grid)
        native = coefficients[grid]
        if native["internal_id"] < 0 or native["kind"] != 2:
            raise ValueError("mesh triangle is not native bulk support")
        nodes = cell["node_ids"]
        if len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError("invalid triangle nodes")
        eps = epsilon[cell["region_id"]]
        if not math.isfinite(eps) or eps <= 0:
            raise ValueError("invalid relative permittivity")
        for k, coefficient in enumerate(native["values"]):
            a, b = nodes[k], nodes[(k + 1) % 3]
            value = coefficient * eps
            for row, col, x in ((a, a, value), (b, b, value), (a, b, -value), (b, a, -value)):
                entries[row, col] = entries.get((row, col), 0.0) + x
    return entries


def node_measures(triangles: list, measures: dict, region_ids: set | None = None) -> dict:
    """Sum native element-vertex support in selected bulk material regions."""
    result = {}
    for cell in triangles:
        if region_ids is not None and cell["region_id"] not in region_ids:
            continue
        native = measures[cell["id"]]
        if native["internal_id"] < 0 or native["kind"] != 2:
            raise ValueError("invalid bulk measure topology")
        for node, value in zip(cell["node_ids"], native["values"], strict=True):
            result[node] = result.get(node, 0.0) + value
    return result


def box_measures_from_coefficients(points: list, coefficients: list) -> tuple:
    """Recover native 2-D element-vertex measures from nonnegative edge kappa.

    Coordinate units are squared in the result. Do not normalize the sum to
    the triangle area: native box transfer can redistribute measure between
    cells, including the non-obtuse neighbors of an obtuse cell.
    """
    if (len(points) != 3 or any(len(p) != 2 for p in points) or
            len(coefficients) != 3):
        raise ValueError("three 2-D points and three edge coefficients required")
    if (any(not math.isfinite(x) for p in points for x in p) or
            any(not math.isfinite(x) or x < 0 for x in coefficients)):
        raise ValueError("nonfinite geometry or negative native coefficient")
    a, b, c = points
    if (b[0]-a[0])*(c[1]-a[1]) == (b[1]-a[1])*(c[0]-a[0]):
        raise ValueError("degenerate triangle")
    products = [coefficients[k] * sum((points[(k+1) % 3][d]-points[k][d])**2
                                    for d in range(2)) for k in range(3)]
    return tuple((products[k] + products[(k-1) % 3]) / 4 for k in range(3))
