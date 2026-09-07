#!/usr/bin/env python3
"""Close general-Tri3 matching-support current vectors and box currents."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_general_tri3_element_edge_avalanche import (
    geometry_rows,
    parse_log,
)
from scripts.diagnose_pn2d_general_tri3_element_edge_current import (
    geometric_partial_volumes,
    gss_laux_vector,
)
from scripts.diagnose_pn2d_general_tri3_imported_state import (
    abs_dex,
    error_summary,
    vector_angle_deg,
)
from scripts.sentaurus_avalanche_replay import (
    currentplot_targets,
)
from scripts.pn2d_general_tri3_contract import EXACT_BIASES_V, SCHEMA_ID


CURRENT_SCHEMA = "pn2d_general_tri3_element_edge_current/v1"
OUTPUT_SCHEMA = "pn2d_general_tri3_element_edge_current_closure/v1"
CARRIERS = ("electron", "hole")
BOX_FLUX_SCALE = 1.0e-8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--edge-analysis-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def least_squares(
    rows: Iterable[tuple[tuple[float, float], float, float]],
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


def charon_whitney_vector(
    points_um: list[tuple[float, float]],
    signed_edge_current: list[float],
) -> tuple[float, float]:
    points = [(x * 1.0e-4, y * 1.0e-4) for x, y in points_um]
    dofs: list[float] = []
    for index in range(3):
        start = points[index]
        end = points[(index + 1) % 3]
        dofs.append(
            signed_edge_current[index] * math.dist(start, end)
        )
    basis = (
        (2.0 / 3.0, 1.0 / 3.0),
        (-1.0 / 3.0, 1.0 / 3.0),
        (-1.0 / 3.0, -2.0 / 3.0),
    )
    ref_x = sum(
        dof * item[0] for dof, item in zip(dofs, basis, strict=True)
    )
    ref_y = sum(
        dof * item[1] for dof, item in zip(dofs, basis, strict=True)
    )
    ax = points[1][0] - points[0][0]
    ay = points[1][1] - points[0][1]
    bx = points[2][0] - points[0][0]
    by = points[2][1] - points[0][1]
    determinant = ax * by - ay * bx
    if abs(determinant) <= 1.0e-30:
        raise ValueError("degenerate physical triangle")
    return (
        (by * ref_x - ay * ref_y) / determinant,
        (-bx * ref_x + ax * ref_y) / determinant,
    )


def candidate_vectors(
    points: list[tuple[float, float]],
    values: list[float],
    read_coefficients: list[float],
) -> dict[str, tuple[float, float]]:
    tangents: list[tuple[float, float]] = []
    lengths: list[float] = []
    for index in range(3):
        start = points[index]
        end = points[(index + 1) % 3]
        delta = (end[0] - start[0], end[1] - start[1])
        length = math.hypot(*delta)
        tangents.append((delta[0] / length, delta[1] / length))
        lengths.append(length)
    all_rows = [
        (tangent, value, 1.0)
        for tangent, value in zip(tangents, values, strict=True)
    ]
    active_rows = [
        (tangents[index], values[index], weight)
        for index in range(3)
        if (
            weight := max(
                read_coefficients[index] * lengths[index],
                0.0,
            )
        )
        > 0.0
    ]
    return {
        "gss_laux_truncated_support": gss_laux_vector(
            points,
            values,
            geometric_partial_volumes(points),
        ),
        "charon_whitney_hcurl_cell_average": charon_whitney_vector(
            points,
            values,
        ),
        "genius_least_squares_tangent": least_squares(all_rows),
        "box_active_edge_exact": least_squares(active_rows),
    }


def signal_status(value: float, state_maximum: float) -> str:
    if value == 0.0:
        return "zero"
    if value <= max(state_maximum * 1.0e-12, 1.0e-300):
        return "below_state_relative_floor"
    return "valid"


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    edge_root = args.edge_analysis_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    raw_manifest_path = raw_root / "manifest.json"
    raw_manifest = json.loads(raw_manifest_path.read_text(encoding="ascii"))
    if raw_manifest.get("schema") != SCHEMA_ID:
        raise ValueError("raw schema mismatch")
    if raw_manifest.get("status") != "passed":
        raise ValueError("raw root is incomplete")
    edge_manifest_path = edge_root / "analysis_manifest.json"
    edge_manifest = json.loads(edge_manifest_path.read_text(encoding="ascii"))
    if edge_manifest.get("schema") != CURRENT_SCHEMA:
        raise ValueError("edge-current schema mismatch")
    if edge_manifest.get("source_manifest_sha256") != digest(
        raw_manifest_path
    ):
        raise ValueError("edge-current/raw manifest mismatch")

    case_names = tuple(raw_manifest["cases"])
    if len(case_names) != 1:
        raise ValueError(f"expected one case, got {case_names}")
    case_name = case_names[0]
    variant_root = raw_root / case_name / "implicit_default"
    log_path = variant_root / "fetched" / "run_implicit_default.out"
    plt_path = (
        variant_root
        / "fetched"
        / "runtime_general_tri3_avalanche_probe_implicit_default.plt"
    )
    groups = parse_log(log_path)
    geometry = geometry_rows(groups, EXACT_BIASES_V[0])
    geometry_by_element = {
        int(row["element"]): row for row in geometry
    }
    vertices = {
        (float(row["bias_V"]), int(row["vertex"])): row
        for row in groups["vertices"]
    }
    local_vertices: dict[tuple[float, int], list[tuple[int, int]]] = (
        defaultdict(list)
    )
    for row in groups["measures"]:
        local_vertices[
            (float(row["bias_V"]), int(row["element"]))
        ].append((int(row["local_vertex"]), int(row["vertex"])))
    raw_edges: dict[tuple[float, int], list[dict[str, Any]]] = defaultdict(
        list
    )
    for row in groups["edges"]:
        raw_edges[
            (float(row["bias_V"]), int(row["element"]))
        ].append(row)

    edge_rows = read_csv(edge_root / "element_local_edges.csv")
    edge_by_key = {
        (
            float(row["bias_V"]),
            int(row["element"]),
            row["carrier"],
            int(row["local_edge"]),
        ): row
        for row in edge_rows
    }

    vector_rows: list[dict[str, Any]] = []
    geometry_rows_output: list[dict[str, Any]] = []
    state_reference_maximum: dict[tuple[float, str, str], float] = (
        defaultdict(float)
    )
    staged_vectors: list[dict[str, Any]] = []
    for key in sorted(local_vertices):
        bias, element_id = key
        ids = [vertex for _, vertex in sorted(local_vertices[key])]
        points = [
            (
                float(vertices[(bias, vertex)]["x_um"]),
                float(vertices[(bias, vertex)]["y_um"]),
            )
            for vertex in ids
        ]
        edge_metadata = [
            edge_by_key[(bias, element_id, "electron", local_edge)]
            for local_edge in range(3)
        ]
        coefficients = [
            float(row["read_coefficient"]) for row in edge_metadata
        ]
        lengths = [float(row["length_um"]) for row in edge_metadata]
        sentaurus_partials = [
            0.5 * coefficient * length * length
            for coefficient, length in zip(
                coefficients,
                lengths,
                strict=True,
            )
        ]
        vela_partials = geometric_partial_volumes(points)
        area = abs(float(geometry_by_element[element_id]["signed_area_um2"]))
        geometry_rows_output.append(
            {
                "case": case_name,
                "element": element_id,
                "element_class": (
                    "contact"
                    if geometry_by_element[element_id]["contact_adjacent"]
                    else "interior"
                ),
                "angle_class": geometry_by_element[element_id]["angle_class"],
                "area_um2": area,
                "sentaurus_signed_partial_sum_um2": sum(
                    sentaurus_partials
                ),
                "sentaurus_area_relative_error": abs(
                    sum(sentaurus_partials) - area
                )
                / area,
                "vela_truncated_partial_sum_um2": sum(vela_partials),
                "vela_truncated_area_relative_error": abs(
                    sum(vela_partials) - area
                )
                / area,
                "negative_read_coefficient_count": sum(
                    coefficient < 0.0 for coefficient in coefficients
                ),
                "zero_read_coefficient_count": sum(
                    coefficient == 0.0 for coefficient in coefficients
                ),
            }
        )
        for carrier in CARRIERS:
            carrier_edges = [
                edge_by_key[(bias, element_id, carrier, local_edge)]
                for local_edge in range(3)
            ]
            vela_values = [
                float(row["vela_sg_current_A_cm2"])
                for row in carrier_edges
            ]
            sentaurus_values = [
                float(row["sentaurus_box_operator_sg_current_A_cm2"])
                for row in carrier_edges
            ]
            vela_candidates = candidate_vectors(
                points,
                vela_values,
                coefficients,
            )
            sentaurus_candidates = candidate_vectors(
                points,
                sentaurus_values,
                coefficients,
            )
            for method in vela_candidates:
                vela_vector = vela_candidates[method]
                reference_vector = sentaurus_candidates[method]
                reference_magnitude = math.hypot(*reference_vector)
                state_reference_maximum[(bias, carrier, method)] = max(
                    state_reference_maximum[(bias, carrier, method)],
                    reference_magnitude,
                )
                staged_vectors.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "element": element_id,
                        "element_class": (
                            "contact"
                            if geometry_by_element[element_id][
                                "contact_adjacent"
                            ]
                            else "interior"
                        ),
                        "angle_class": geometry_by_element[element_id][
                            "angle_class"
                        ],
                        "carrier": carrier,
                        "method": method,
                        "vela_vector": vela_vector,
                        "reference_vector": reference_vector,
                    }
                )

    for row in staged_vectors:
        vela_vector = row.pop("vela_vector")
        reference_vector = row.pop("reference_vector")
        vela_magnitude = math.hypot(*vela_vector)
        reference_magnitude = math.hypot(*reference_vector)
        maximum = state_reference_maximum[
            (float(row["bias_V"]), row["carrier"], row["method"])
        ]
        status = signal_status(reference_magnitude, maximum)
        vector_rows.append(
            {
                **row,
                "vela_x_A_cm2": vela_vector[0],
                "vela_y_A_cm2": vela_vector[1],
                "vela_magnitude_A_cm2": vela_magnitude,
                "sentaurus_box_x_A_cm2": reference_vector[0],
                "sentaurus_box_y_A_cm2": reference_vector[1],
                "sentaurus_box_magnitude_A_cm2": reference_magnitude,
                "magnitude_absolute_error_dex": (
                    ""
                    if status != "valid"
                    else abs_dex(vela_magnitude, reference_magnitude)
                ),
                "vector_angle_error_deg": (
                    ""
                    if status != "valid"
                    else vector_angle_deg(vela_vector, reference_vector)
                ),
                "signal_status": status,
                "observation_label": "box_operator_reconstruction",
            }
        )

    currentplot = {
        float(row["bias_V"]): row
        for row in currentplot_targets(
            plt_path,
            tuple(float(value) for value in EXACT_BIASES_V),
        )
    }
    terminal_rows: list[dict[str, Any]] = []
    kcl_rows: list[dict[str, Any]] = []
    for bias in (float(value) for value in EXACT_BIASES_V):
        used_vertices = {
            vertex
            for (row_bias, _), values in local_vertices.items()
            if row_bias == bias
            for _, vertex in values
        }
        x_by_vertex = {
            vertex: float(vertices[(bias, vertex)]["x_um"])
            for vertex in used_vertices
        }
        min_x = min(x_by_vertex.values())
        max_x = max(x_by_vertex.values())
        contacts = {
            "Anode": {
                vertex for vertex, x in x_by_vertex.items()
                if abs(x - min_x) <= 1.0e-12
            },
            "Cathode": {
                vertex for vertex, x in x_by_vertex.items()
                if abs(x - max_x) <= 1.0e-12
            },
        }
        endpoints: dict[int, tuple[int, int]] = {}
        fluxes: dict[str, dict[str, defaultdict[int, float]]] = {
            method: {
                carrier: defaultdict(float) for carrier in CARRIERS
            }
            for method in ("sentaurus_box_operator", "vela_recomputed")
        }
        for (row_bias, element_id), raw_group in raw_edges.items():
            if row_bias != bias:
                continue
            ids = [
                vertex
                for _, vertex in sorted(
                    local_vertices[(row_bias, element_id)]
                )
            ]
            for local_edge in range(3):
                node0 = ids[local_edge]
                node1 = ids[(local_edge + 1) % 3]
                raw = next(
                    edge for edge in raw_group
                    if {
                        int(edge["start"]),
                        int(edge["end"]),
                    }
                    == {node0, node1}
                )
                edge_id = int(raw["edge"])
                endpoints[edge_id] = (
                    int(raw["start"]),
                    int(raw["end"]),
                )
                orientation = (
                    1.0
                    if endpoints[edge_id] == (node0, node1)
                    else -1.0
                )
                length_um = float(raw["length_um"])
                coefficient = float(raw["kappa"])
                for carrier, raw_key in (
                    ("electron", "box_flux_n_A_um"),
                    ("hole", "box_flux_p_A_um"),
                ):
                    fluxes["sentaurus_box_operator"][carrier][
                        edge_id
                    ] += float(raw[raw_key])
                    edge_row = edge_by_key[
                        (bias, element_id, carrier, local_edge)
                    ]
                    raw_vela_current = orientation * float(
                        edge_row["vela_sg_current_A_cm2"]
                    )
                    fluxes["vela_recomputed"][carrier][edge_id] += (
                        raw_vela_current
                        * coefficient
                        * length_um
                        * BOX_FLUX_SCALE
                    )

        balances: dict[
            str,
            dict[str, dict[int, float]],
        ] = defaultdict(dict)
        for method in fluxes:
            for carrier in CARRIERS:
                balance = {vertex: 0.0 for vertex in used_vertices}
                for edge_id, value in fluxes[method][carrier].items():
                    start, end = endpoints[edge_id]
                    balance[start] -= value
                    balance[end] += value
                balances[method][carrier] = balance
                label = "eCurrent" if carrier == "electron" else "hCurrent"
                for contact, nodes in contacts.items():
                    predicted = sum(balance[node] for node in nodes)
                    reference = float(currentplot[bias][f"{contact} {label}"])
                    terminal_rows.append(
                        {
                            "case": case_name,
                            "bias_V": bias,
                            "method": method,
                            "carrier": carrier,
                            "contact": contact,
                            "predicted_A_um": predicted,
                            "sentaurus_terminal_A_um": reference,
                            "absolute_error_A_um": abs(
                                predicted - reference
                            ),
                            "relative_error": abs(predicted - reference)
                            / max(abs(reference), 1.0e-300),
                        }
                    )
        for method in balances:
            for contact, nodes in contacts.items():
                predicted = sum(
                    balances[method][carrier][node]
                    for carrier in CARRIERS
                    for node in nodes
                )
                reference = float(
                    currentplot[bias][f"{contact} TotalCurrent"]
                )
                terminal_rows.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "method": method,
                        "carrier": "total",
                        "contact": contact,
                        "predicted_A_um": predicted,
                        "sentaurus_terminal_A_um": reference,
                        "absolute_error_A_um": abs(predicted - reference),
                        "relative_error": abs(predicted - reference)
                        / max(abs(reference), 1.0e-300),
                    }
                )
        terminal_scale = max(
            abs(float(currentplot[bias]["Anode TotalCurrent"])),
            abs(float(currentplot[bias]["Cathode TotalCurrent"])),
            1.0e-300,
        )
        internal = (
            used_vertices - contacts["Anode"] - contacts["Cathode"]
        )
        for method in balances:
            for vertex in sorted(internal):
                residual = sum(
                    balances[method][carrier][vertex]
                    for carrier in CARRIERS
                )
                kcl_rows.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "method": method,
                        "vertex": vertex,
                        "total_current_residual_A_um": residual,
                        "terminal_current_scale_A_um": terminal_scale,
                        "relative_residual": abs(residual) / terminal_scale,
                    }
                )

    outputs = {
        "matching_support_cell_vectors.csv": vector_rows,
        "geometry_support_closure.csv": geometry_rows_output,
        "terminal_current_closure.csv": terminal_rows,
        "internal_kcl.csv": kcl_rows,
    }
    for name, rows in outputs.items():
        write_csv(output_root / name, rows)

    vector_summary: dict[str, Any] = {}
    methods = sorted({row["method"] for row in vector_rows})
    for method in methods:
        vector_summary[method] = {}
        for carrier in CARRIERS:
            selected = [
                row for row in vector_rows
                if row["method"] == method
                and row["carrier"] == carrier
                and row["signal_status"] == "valid"
            ]
            vector_summary[method][carrier] = {
                "magnitude_error_dex": error_summary(
                    [
                        float(row["magnitude_absolute_error_dex"])
                        for row in selected
                    ]
                ),
                "angle_error_deg": error_summary(
                    [
                        float(row["vector_angle_error_deg"])
                        for row in selected
                        if row["vector_angle_error_deg"] != ""
                    ]
                ),
                "valid_count": len(selected),
            }
    closure_summary = {
        method: {
            "maximum_terminal_total_relative_error": max(
                float(row["relative_error"])
                for row in terminal_rows
                if row["method"] == method
                and row["carrier"] == "total"
            ),
            "maximum_terminal_carrier_relative_error": max(
                float(row["relative_error"])
                for row in terminal_rows
                if row["method"] == method
                and row["carrier"] != "total"
            ),
            "maximum_internal_total_kcl_relative_residual": max(
                float(row["relative_residual"])
                for row in kcl_rows
                if row["method"] == method
            ),
        }
        for method in ("sentaurus_box_operator", "vela_recomputed")
    }
    geometry_summary = {
        "maximum_sentaurus_area_relative_error": max(
            float(row["sentaurus_area_relative_error"])
            for row in geometry_rows_output
        ),
        "maximum_vela_truncated_area_relative_error": max(
            float(row["vela_truncated_area_relative_error"])
            for row in geometry_rows_output
        ),
        "negative_read_coefficient_element_count": sum(
            int(row["negative_read_coefficient_count"]) > 0
            for row in geometry_rows_output
        ),
    }
    manifest = {
        "schema": OUTPUT_SCHEMA,
        "status": "valid",
        "case_name": case_name,
        "source_raw_manifest_sha256": digest(raw_manifest_path),
        "source_edge_manifest_sha256": digest(edge_manifest_path),
        "source_log_sha256": digest(log_path),
        "source_currentplot_sha256": digest(plt_path),
        "exact_biases_V": list(EXACT_BIASES_V),
        "native_element_current_observation": (
            "insufficient_native_observation_undocumented_element_vector"
        ),
        "matching_support_cell_vector": vector_summary,
        "geometry_support": geometry_summary,
        "closure": closure_summary,
        "outputs": {
            name: digest(output_root / name) for name in outputs
        },
    }
    (output_root / "analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
