"""Shared Sentaurus avalanche log parsing and current reconstruction."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


TARGET_BIASES = (-1.0, -10.0, -20.0)


def parse_tokens(line: str) -> dict[str, str]:
    return dict(re.findall(r"(\w+)=([^\s]+)", line))


def typed_row(tokens: dict[str, str], integer_keys: set[str]) -> dict[str, Any]:
    row: dict[str, Any] = {}
    for key, value in tokens.items():
        if key in integer_keys:
            row[key] = int(value)
        else:
            row[key] = float(value)
    return row


def parse_log(
    path: Path,
    target_biases: tuple[float, ...] = TARGET_BIASES,
) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    prefixes = {
        "AVAL_PROBE_VERTEX ": ("vertices", {"bias_V", "vertex"}),
        "AVAL_PROBE_ELEMENT ": ("elements", {"bias_V", "element"}),
        "AVAL_PROBE_MEASURE ": (
            "measures",
            {"bias_V", "element", "local_vertex", "vertex"},
        ),
        "AVAL_PROBE_EDGE ": (
            "edges",
            {
                "bias_V",
                "element",
                "local_edge",
                "edge",
                "start",
                "end",
            },
        ),
        "AVAL_PROBE_INTEGRAL ": ("integrals", {"bias_V"}),
    }
    text = path.read_text(encoding="ascii", errors="strict")
    for line in text.splitlines():
        for prefix, (group, integer_keys) in prefixes.items():
            if line.startswith(prefix):
                groups[group].append(
                    typed_row(parse_tokens(line[len(prefix) :]), integer_keys)
                )
                break
    state_count = len(target_biases)
    expected = {
        "vertices": state_count * 10,
        "elements": state_count * 4,
        "measures": state_count * 12,
        "edges": state_count * 12,
        "integrals": state_count,
    }
    for group, count in expected.items():
        if len(groups[group]) != count:
            raise ValueError(
                f"{path}: expected {count} {group}, got {len(groups[group])}"
            )
        observed_biases = {float(row["bias_V"]) for row in groups[group]}
        if observed_biases != set(target_biases):
            raise ValueError(
                f"{path}: {group} bias matrix mismatch: "
                f"expected {sorted(target_biases)}, got {sorted(observed_biases)}"
            )
    return groups


def parse_plt(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    text = path.read_text(encoding="ascii", errors="strict")
    info, data = text.split("Data {", 1)
    match = re.search(r"datasets\s*=\s*\[(.*?)\]", info, re.S)
    if match is None:
        raise ValueError(f"{path}: missing datasets")
    names = re.findall(r'"([^"]+)"', match.group(1))
    values = [
        float(token)
        for token in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?", data
        )
    ]
    if len(values) % len(names):
        raise ValueError(f"{path}: non-rectangular CurrentPlot")
    rows = [
        dict(zip(names, values[index : index + len(names)], strict=True))
        for index in range(0, len(values), len(names))
    ]
    return names, rows


def currentplot_targets(
    path: Path,
    target_biases: tuple[float, ...] = TARGET_BIASES,
) -> list[dict[str, float]]:
    names, rows = parse_plt(path)
    voltage_name = next(
        name
        for name in names
        if name.endswith("Anode OuterVoltage") or name == "Anode OuterVoltage"
    )
    result = []
    for bias in target_biases:
        match = min(rows, key=lambda row: abs(row[voltage_name] - bias))
        if abs(match[voltage_name] - bias) > 1.0e-8:
            raise ValueError(f"{path}: missing CurrentPlot bias {bias:g}")
        result.append({"bias_V": bias, **match})
    return result


def solve_pair(
    tangent_a: tuple[float, float],
    value_a: float,
    tangent_b: tuple[float, float],
    value_b: float,
) -> tuple[float, float]:
    ax, ay = tangent_a
    bx, by = tangent_b
    determinant = ax * by - ay * bx
    if abs(determinant) <= 1.0e-14:
        raise ValueError("parallel edge tangents cannot reconstruct a vector")
    return (
        (value_a * by - ay * value_b) / determinant,
        (ax * value_b - value_a * bx) / determinant,
    )


def least_squares(
    rows: Iterable[tuple[tuple[float, float], float, float]]
) -> tuple[float, float]:
    a11 = a12 = a22 = b1 = b2 = 0.0
    for (tx, ty), value, weight in rows:
        a11 += weight * tx * tx
        a12 += weight * tx * ty
        a22 += weight * ty * ty
        b1 += weight * tx * value
        b2 += weight * ty * value
    determinant = a11 * a22 - a12 * a12
    if abs(determinant) <= 1.0e-30:
        raise ValueError("singular weighted edge reconstruction")
    return (
        (b1 * a22 - b2 * a12) / determinant,
        (a11 * b2 - a12 * b1) / determinant,
    )


def gss_laux_vector(
    edges: list[dict[str, Any]], current_key: str
) -> tuple[float, float]:
    pair_vectors: dict[tuple[int, int], tuple[float, float]] = {}
    for first in range(3):
        for second in range(first + 1, 3):
            pair_vectors[(first, second)] = solve_pair(
                (edges[first]["tangent_x"], edges[first]["tangent_y"]),
                edges[first][current_key],
                (edges[second]["tangent_x"], edges[second]["tangent_y"]),
                edges[second][current_key],
            )

    edge_vectors: list[tuple[tuple[float, float], float]] = []
    for target in range(3):
        others = [index for index in range(3) if index != target]
        first, second = others
        first_weight = edges[first]["kappa"] * edges[first]["length_um"]
        second_weight = edges[second]["kappa"] * edges[second]["length_um"]
        first_pair = pair_vectors[tuple(sorted((target, first)))]
        second_pair = pair_vectors[tuple(sorted((target, second)))]
        weight_sum = first_weight + second_weight
        if weight_sum <= 0.0:
            vector = (
                0.5 * (first_pair[0] + second_pair[0]),
                0.5 * (first_pair[1] + second_pair[1]),
            )
        else:
            vector = (
                (first_weight * first_pair[0] + second_weight * second_pair[0])
                / weight_sum,
                (first_weight * first_pair[1] + second_weight * second_pair[1])
                / weight_sum,
            )
        partial_area = (
            0.5
            * edges[target]["kappa"]
            * edges[target]["length_um"]
            * edges[target]["length_um"]
        )
        edge_vectors.append((vector, partial_area))
    total_area = sum(weight for _, weight in edge_vectors)
    if total_area <= 0.0:
        raise ValueError("triangle has no positive box partial area")
    return (
        sum(vector[0] * weight for vector, weight in edge_vectors) / total_area,
        sum(vector[1] * weight for vector, weight in edge_vectors) / total_area,
    )


def charon_whitney_vector(
    edges: list[dict[str, Any]],
    local_vertices: list[int],
    vertex_by_id: dict[int, dict[str, Any]],
    current_key: str,
) -> tuple[float, float]:
    if len(local_vertices) != 3:
        raise ValueError("triangle must have three local vertices")
    points = [
        (
            vertex_by_id[vertex]["x_um"] * 1.0e-4,
            vertex_by_id[vertex]["y_um"] * 1.0e-4,
        )
        for vertex in local_vertices
    ]
    reference_edges = (
        (local_vertices[0], local_vertices[1]),
        (local_vertices[1], local_vertices[2]),
        (local_vertices[2], local_vertices[0]),
    )
    edge_by_nodes = {
        tuple(sorted((int(edge["start"]), int(edge["end"])))): edge
        for edge in edges
    }
    dofs = []
    for start, end in reference_edges:
        edge = edge_by_nodes[tuple(sorted((start, end)))]
        orientation = (
            1.0
            if (int(edge["start"]), int(edge["end"])) == (start, end)
            else -1.0
        )
        length_cm = edge["length_um"] * 1.0e-4
        dofs.append(orientation * edge[current_key] * length_cm)

    # Lowest-order Whitney edge basis at the reference-triangle centroid.
    basis = ((2.0 / 3.0, 1.0 / 3.0), (-1.0 / 3.0, 1.0 / 3.0),
             (-1.0 / 3.0, -2.0 / 3.0))
    ref_x = sum(dof * item[0] for dof, item in zip(dofs, basis, strict=True))
    ref_y = sum(dof * item[1] for dof, item in zip(dofs, basis, strict=True))
    ax = points[1][0] - points[0][0]
    ay = points[1][1] - points[0][1]
    bx = points[2][0] - points[0][0]
    by = points[2][1] - points[0][1]
    determinant = ax * by - ay * bx
    if abs(determinant) <= 1.0e-30:
        raise ValueError("degenerate physical triangle")
    # Covariant Piola transform: physical = A^{-T} reference.
    return (
        (by * ref_x - ay * ref_y) / determinant,
        (-bx * ref_x + ax * ref_y) / determinant,
    )
