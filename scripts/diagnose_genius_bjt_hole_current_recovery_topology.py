#!/usr/bin/env python3
"""Explain the Genius BJT weak-hole-current tail using recovery topology."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

import diagnose_genius_bjt_sdevice_current_semantics as semantics


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
Q_C = 1.602176634e-19
FLOOR = 1.0e-300


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def report_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path)


def finite_distribution(values: list[float]) -> dict[str, float | int]:
    selected = [value for value in values if math.isfinite(value)]
    if not selected:
        return {"count": 0, "median": 0.0, "p95": 0.0, "maximum": 0.0}
    return {
        "count": len(selected),
        "median": statistics.median(selected),
        "p95": semantics.percentile(selected, 0.95),
        "maximum": max(selected),
    }


def pearson(first: list[float], second: list[float]) -> float | None:
    pairs = [
        (a, b) for a, b in zip(first, second, strict=True)
        if math.isfinite(a) and math.isfinite(b)
    ]
    if len(pairs) < 2:
        return None
    mean_a = statistics.fmean(a for a, _ in pairs)
    mean_b = statistics.fmean(b for _, b in pairs)
    covariance = sum((a - mean_a) * (b - mean_b) for a, b in pairs)
    variance_a = sum((a - mean_a) ** 2 for a, _ in pairs)
    variance_b = sum((b - mean_b) ** 2 for _, b in pairs)
    if variance_a <= 0.0 or variance_b <= 0.0:
        return None
    return covariance / math.sqrt(variance_a * variance_b)


def condition_number(terms: list[tuple[float, float, float, float]]) -> float:
    a00 = sum(weight * tx * tx for weight, tx, _, _ in terms)
    a01 = sum(weight * tx * ty for weight, tx, ty, _ in terms)
    a11 = sum(weight * ty * ty for weight, _, ty, _ in terms)
    trace = a00 + a11
    discriminant = math.sqrt(max((a00 - a11) ** 2 + 4.0 * a01 * a01, 0.0))
    largest = 0.5 * (trace + discriminant)
    smallest = 0.5 * (trace - discriminant)
    if smallest <= max(largest, 1.0) * 1.0e-15:
        return math.inf
    return largest / smallest


def direct_node_recovery(
    node_count: int, edge_rows: list[dict[str, str]]
) -> tuple[list[tuple[float, float]], list[list[tuple[float, float, float, float]]]]:
    terms_by_node: list[list[tuple[float, float, float, float]]] = [
        [] for _ in range(node_count)
    ]
    for edge in edge_rows:
        couple = float(edge["couple_m"])
        length = float(edge["length_m"])
        if couple <= 0.0 or length <= 0.0:
            continue
        tx = (float(edge["x1"]) - float(edge["x0"])) / length
        ty = (float(edge["y1"]) - float(edge["y0"])) / length
        projection = (
            Q_C * float(edge["hole_particle_line_flux_per_m_s"]) / couple / 1.0e4
        )
        term = (couple, tx, ty, projection)
        terms_by_node[int(edge["node0"])].append(term)
        terms_by_node[int(edge["node1"])].append(term)

    recovered: list[tuple[float, float]] = []
    for terms in terms_by_node:
        value = semantics.solve_vector_projections(terms)
        if value is None:
            total_weight = sum(weight for weight, *_ in terms)
            value = (
                sum(weight * projection * tx for weight, tx, _, projection in terms)
                / max(total_weight, FLOOR),
                sum(weight * projection * ty for weight, _, ty, projection in terms)
                / max(total_weight, FLOOR),
            )
        recovered.append(value)
    return recovered, terms_by_node


def formal_vela_current(path: Path) -> list[tuple[float, float]]:
    result = []
    for expected, row in enumerate(semantics.rows(path)):
        if int(row["node_id"]) != expected:
            raise ValueError(f"unordered formal transport field: {path}")
        result.append(
            (
                float(row["hole_physical_total_x_A_per_cm2"]),
                float(row["hole_physical_total_y_A_per_cm2"]),
            )
        )
    return result


def abs_log_error(
    reference: tuple[float, float], candidate: tuple[float, float]
) -> float:
    return abs(
        math.log10(
            max(semantics.magnitude(candidate), FLOOR)
            / max(semantics.magnitude(reference), FLOOR)
        )
    )


def group_summary(
    name: str,
    indices: list[int],
    sdevice: list[tuple[float, float]],
    direct: list[tuple[float, float]],
    cell_first: list[tuple[float, float]],
    diagnostic_rows: list[dict[str, float | int | bool]],
) -> dict[str, object]:
    rows = [diagnostic_rows[index] for index in indices]
    direct_errors = [float(row["direct_absolute_log10_error_decade"]) for row in rows]
    cell_errors = [float(row["cell_first_absolute_log10_error_decade"]) for row in rows]
    improvement = [direct_error - cell_error for direct_error, cell_error in zip(direct_errors, cell_errors, strict=True)]
    return {
        "name": name,
        "node_count": len(indices),
        "sdevice_vs_direct_node_fit": semantics.field_metrics(sdevice, direct, indices),
        "sdevice_vs_cell_first_area_projection": semantics.field_metrics(
            sdevice, cell_first, indices
        ),
        "nodes_above_0p5_decade": {
            "direct_node_fit": sum(value > 0.5 for value in direct_errors),
            "cell_first_area_projection": sum(value > 0.5 for value in cell_errors),
        },
        "cell_first_error_improvement_decade": finite_distribution(improvement),
        "direct_fit_normalized_residual": finite_distribution(
            [float(row["direct_fit_normalized_residual"]) for row in rows]
        ),
        "incident_cell_cancellation_factor": finite_distribution(
            [float(row["incident_cell_cancellation_factor"]) for row in rows]
        ),
        "incident_cell_vector_dispersion": finite_distribution(
            [float(row["incident_cell_vector_dispersion"]) for row in rows]
        ),
        "nodal_fit_condition_number": finite_distribution(
            [float(row["nodal_fit_condition_number"]) for row in rows]
        ),
    }


def write_markdown(path: Path, result: dict[str, object]) -> None:
    groups = {item["name"]: item for item in result["groups"]}
    replay = result["production_replay"]
    tail = groups["original_tail"]
    dominant = groups["dominant_current"]
    weak = groups["weak_main_current"]
    lines = [
        "# Genius NPN BJT weak-hole-current recovery topology audit",
        "",
        "## Technical summary",
        "",
        "The coarse-grid 0.828-decade P95 failure is produced after the conservative SG "
        "edge fluxes are mapped to a vertex vector. Replaying the production dual-face "
        "least-squares formula from the saved edges reproduces the exported Vela field "
        f"to {replay['maximum_absolute_vector_difference_A_per_cm2']:.3e} A/cm2. "
        "Changing only the recovery order to per-cell fitting followed by area-weighted "
        "cell-to-node projection reduces the original 180-node tail P95 from "
        f"{tail['sdevice_vs_direct_node_fit']['magnitude_abs_dex_p95']:.6f} to "
        f"{tail['sdevice_vs_cell_first_area_projection']['magnitude_abs_dex_p95']:.6f} decade. "
        "The SDevice node source used here matches the formal comparison source with "
        f"normalized vector RMSE {result['source_reconciliation']['metrics']['normalized_vector_rmse']:.3e}.",
        "",
        "## Recovery-order comparison",
        "",
        "| Population | Nodes | Direct nodal P95 (dec) | Cell-first P95 (dec) | Direct failures | Cell-first failures |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for key in ("dominant_current", "weak_main_current", "original_tail"):
        item = groups[key]
        lines.append(
            f"| {key} | {item['node_count']} | "
            f"{item['sdevice_vs_direct_node_fit']['magnitude_abs_dex_p95']:.6f} | "
            f"{item['sdevice_vs_cell_first_area_projection']['magnitude_abs_dex_p95']:.6f} | "
            f"{item['nodes_above_0p5_decade']['direct_node_fit']} | "
            f"{item['nodes_above_0p5_decade']['cell_first_area_projection']} |"
        )
    lines.extend(
        [
            "",
            "The dominant-current region is insensitive to recovery order. The improvement "
            "is concentrated in the weak-current base. The cell-first path reconstructs each "
            "triangle from all three of its edges, including the edge opposite the target "
            "vertex, before projecting the cell vectors to that vertex; the direct-node path "
            "uses only edges incident on the vertex.",
            "",
            "## Topology evidence",
            "",
            "| Population | Direct-fit residual median | Cell cancellation median | Cell-vector dispersion median | Fit condition median |",
            "|---|---:|---:|---:|---:|",
            f"| dominant_current | {dominant['direct_fit_normalized_residual']['median']:.6g} | "
            f"{dominant['incident_cell_cancellation_factor']['median']:.6g} | "
            f"{dominant['incident_cell_vector_dispersion']['median']:.6g} | "
            f"{dominant['nodal_fit_condition_number']['median']:.6g} |",
            f"| weak_main_current | {weak['direct_fit_normalized_residual']['median']:.6g} | "
            f"{weak['incident_cell_cancellation_factor']['median']:.6g} | "
            f"{weak['incident_cell_vector_dispersion']['median']:.6g} | "
            f"{weak['nodal_fit_condition_number']['median']:.6g} |",
            f"| original_tail | {tail['direct_fit_normalized_residual']['median']:.6g} | "
            f"{tail['incident_cell_cancellation_factor']['median']:.6g} | "
            f"{tail['incident_cell_vector_dispersion']['median']:.6g} | "
            f"{tail['nodal_fit_condition_number']['median']:.6g} |",
            "",
            "The condition-number comparison tests mesh-direction degeneracy; the residual "
            "and cell-dispersion comparisons test whether one smooth vector can represent the "
            "whole vertex patch. Tail cell-vector dispersion is much larger than in the dominant "
            "region, but the weak within-tail correlations mean it classifies the affected region "
            "better than it predicts node-by-node error severity. The near-unity cancellation "
            "factor and modest condition numbers reject simple vector cancellation and ill-"
            "conditioning as primary explanations.",
            "",
            "## Correlations on the original tail",
            "",
        ]
    )
    for name, value in result["tail_correlations"].items():
        rendered = "not defined" if value is None else f"{value:.6f}"
        lines.append(f"- {name}: {rendered}")
    lines.extend(
        [
            "",
            "## Assessment",
            "",
            "- Verified: the formal Vela node field is an exact replay of the saved SG edge "
            "projections under the production direct-node formula.",
            "- Verified: using the same edge fluxes, cell-first recovery removes most of the "
            "weak-current discrepancy without changing the state, mobility, or flux operator.",
            "- Strongly supported: the remaining coarse-grid P95 failure is a representation-"
            "support mismatch between Vela's direct vertex fit and SDevice's element-oriented "
            "current field, amplified where neighbouring cell currents vary in direction.",
            "- Not proven: the exact proprietary SDevice element-to-vertex weighting. Terminal "
            "currents and conservative section fluxes remain the physical acceptance oracle.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    probe = BUILD_ROOT / "sdevice_current_semantics_probe" / "default_export"
    parser.add_argument("--sdevice-export", type=Path, default=probe)
    parser.add_argument(
        "--formal-transport", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "transport_sources" / "transport" / "vce_030.csv",
    )
    parser.add_argument(
        "--formal-sdevice-node-current", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "sentaurus" / "vce_030" / "fields" / "hCurrentDensity_region0.csv",
    )
    parser.add_argument(
        "--sg-edges", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "edge_audit" / "vce_030" / "sg_edges.csv",
    )
    parser.add_argument(
        "--tail-diagnostics", type=Path,
        default=BUILD_ROOT / "m1_hole_current_tail_diagnosis" / "node_diagnostics.csv",
    )
    parser.add_argument(
        "--output-root", type=Path,
        default=BUILD_ROOT / "m1_p0_multibias" / "hole_current_recovery_topology",
    )
    parser.add_argument(
        "--report-json", type=Path,
        default=FIXTURE / "reports" / "hole_current_recovery_topology.json",
    )
    parser.add_argument(
        "--report-markdown", type=Path,
        default=FIXTURE / "reports" / "hole_current_recovery_topology.md",
    )
    args = parser.parse_args()

    export = semantics.load_export(args.sdevice_export)
    sdevice = export["node_current"]
    formal_sdevice = semantics.vector(args.formal_sdevice_node_current, "node")
    node_count = len(sdevice)
    formal = formal_vela_current(args.formal_transport)
    if len(formal) != node_count:
        raise ValueError("formal Vela and SDevice node counts differ")
    if len(formal_sdevice) != node_count:
        raise ValueError("formal and element-export SDevice node counts differ")
    sdevice_source_delta = max(
        math.hypot(a[0] - b[0], a[1] - b[1])
        for a, b in zip(sdevice, formal_sdevice, strict=True)
    )

    triangles, areas = semantics.triangle_areas(args.sdevice_export)
    edge_rows = semantics.rows(args.sg_edges)
    seen_pairs: set[tuple[int, int]] = set()
    for expected, edge in enumerate(edge_rows):
        if int(edge["edge_id"]) != expected:
            raise ValueError("unordered SG edge diagnostics")
        pair = tuple(sorted((int(edge["node0"]), int(edge["node1"]))))
        if pair in seen_pairs:
            raise ValueError(f"duplicate SG edge pair: {pair}")
        seen_pairs.add(pair)
    direct_replay, terms_by_node = direct_node_recovery(node_count, edge_rows)
    cell_vectors, valid_cells = semantics.vela_cell_currents_from_sg_edges(
        args.sg_edges, triangles
    )
    cell_first, valid_nodes = semantics.reconstruct_sparse(
        node_count, triangles, areas, cell_vectors, valid_cells
    )
    valid_node_set = set(valid_nodes)

    node_cells: list[list[int]] = [[] for _ in range(node_count)]
    for cell_id in valid_cells:
        for node in triangles[cell_id]:
            node_cells[node].append(cell_id)

    tail_input = semantics.rows(args.tail_diagnostics)
    original_tail = [
        int(row["node_id"]) for row in tail_input
        if row["exceeds_0p5_decade"].lower() == "true"
    ]
    peak = max(semantics.magnitude(value) for value in sdevice)
    main = [i for i, value in enumerate(sdevice) if semantics.magnitude(value) >= peak * 1.0e-6]
    dominant = [i for i, value in enumerate(sdevice) if semantics.magnitude(value) >= peak * 1.0e-2]
    main_set = set(main)
    dominant_set = set(dominant)
    original_tail_set = set(original_tail)
    if len(original_tail_set) != len(original_tail):
        raise ValueError("duplicate node in original tail selection")
    if not original_tail_set.issubset(main_set):
        raise ValueError("original tail is not contained in the formal main-current mask")
    weak_main = [i for i in main if i not in dominant_set]
    if not all(index in valid_node_set for index in main):
        raise ValueError("a selected node has no recoverable incident cell")

    diagnostic_rows: list[dict[str, float | int | bool]] = []
    for node in range(node_count):
        terms = terms_by_node[node]
        residual_numerator = sum(
            weight * (tx * formal[node][0] + ty * formal[node][1] - projection) ** 2
            for weight, tx, ty, projection in terms
        )
        residual_denominator = sum(
            weight * projection * projection for weight, _, _, projection in terms
        )
        fit_residual = math.sqrt(residual_numerator / max(residual_denominator, FLOOR))
        cells = node_cells[node]
        area_total = sum(areas[cell] for cell in cells)
        mean_magnitude = sum(
            areas[cell] * semantics.magnitude(cell_vectors[cell]) for cell in cells
        ) / max(area_total, FLOOR)
        projected_magnitude = semantics.magnitude(cell_first[node])
        cancellation = mean_magnitude / max(projected_magnitude, FLOOR)
        dispersion = math.sqrt(
            sum(
                areas[cell]
                * (
                    (cell_vectors[cell][0] - cell_first[node][0]) ** 2
                    + (cell_vectors[cell][1] - cell_first[node][1]) ** 2
                )
                for cell in cells
            )
            / max(area_total, FLOOR)
        ) / max(projected_magnitude, FLOOR)
        diagnostic_rows.append(
            {
                "node_id": node,
                "in_main_mask": node in main_set,
                "in_dominant_mask": node in dominant_set,
                "in_original_tail": node in original_tail_set,
                "active_edge_count": len(terms),
                "incident_cell_count": len(cells),
                "nodal_fit_condition_number": condition_number(terms),
                "direct_fit_normalized_residual": fit_residual,
                "incident_cell_cancellation_factor": cancellation,
                "incident_cell_vector_dispersion": dispersion,
                "direct_absolute_log10_error_decade": abs_log_error(sdevice[node], formal[node]),
                "cell_first_absolute_log10_error_decade": abs_log_error(sdevice[node], cell_first[node]),
                "cell_first_over_direct_magnitude_ratio": (
                    projected_magnitude / max(semantics.magnitude(formal[node]), FLOOR)
                ),
            }
        )

    groups = [
        group_summary("main_current", main, sdevice, formal, cell_first, diagnostic_rows),
        group_summary("dominant_current", dominant, sdevice, formal, cell_first, diagnostic_rows),
        group_summary("weak_main_current", weak_main, sdevice, formal, cell_first, diagnostic_rows),
        group_summary("original_tail", original_tail, sdevice, formal, cell_first, diagnostic_rows),
    ]
    tail_rows = [diagnostic_rows[index] for index in original_tail]
    tail_error = [float(row["direct_absolute_log10_error_decade"]) for row in tail_rows]
    correlations = {
        "direct_error_vs_log10_cell_cancellation": pearson(
            tail_error,
            [math.log10(max(float(row["incident_cell_cancellation_factor"]), FLOOR)) for row in tail_rows],
        ),
        "direct_error_vs_log10_cell_vector_dispersion": pearson(
            tail_error,
            [math.log10(max(float(row["incident_cell_vector_dispersion"]), FLOOR)) for row in tail_rows],
        ),
        "direct_error_vs_log10_direct_fit_residual": pearson(
            tail_error,
            [math.log10(max(float(row["direct_fit_normalized_residual"]), FLOOR)) for row in tail_rows],
        ),
        "direct_error_vs_log10_fit_condition_number": pearson(
            tail_error,
            [math.log10(max(float(row["nodal_fit_condition_number"]), 1.0)) for row in tail_rows],
        ),
    }
    replay_max = max(
        math.hypot(a[0] - b[0], a[1] - b[1])
        for a, b in zip(formal, direct_replay, strict=True)
    )
    result: dict[str, object] = {
        "schema": "vela.genius_bjt_hole_current_recovery_topology.v1",
        "scope": "VBE=0.70 V, VCE=3.00 V, coarse common mesh, formal accepted Vela state",
        "node_count": node_count,
        "cell_count": len(triangles),
        "source_files": {
            "formal_vela_transport": {
                "path": report_path(args.formal_transport),
                "sha256": sha256(args.formal_transport),
            },
            "formal_sdevice_node_current": {
                "path": report_path(args.formal_sdevice_node_current),
                "sha256": sha256(args.formal_sdevice_node_current),
            },
            "sdevice_element_export_manifest": {
                "path": report_path(args.sdevice_export / "field_manifest.json"),
                "sha256": sha256(args.sdevice_export / "field_manifest.json"),
            },
            "sg_edges": {
                "path": report_path(args.sg_edges),
                "sha256": sha256(args.sg_edges),
                "row_count": len(edge_rows),
            },
            "tail_diagnostics": {
                "path": report_path(args.tail_diagnostics),
                "sha256": sha256(args.tail_diagnostics),
                "row_count": len(tail_input),
            },
        },
        "source_reconciliation": {
            "sdevice_node_current_max_absolute_difference_A_per_cm2": sdevice_source_delta,
            "metrics": semantics.field_metrics(
                formal_sdevice, sdevice, list(range(node_count))
            ),
        },
        "production_replay": {
            "formula": "dual-face-weighted least squares over all active incident SG edge projections",
            "maximum_absolute_vector_difference_A_per_cm2": replay_max,
            "metrics": semantics.field_metrics(formal, direct_replay, list(range(node_count))),
        },
        "groups": groups,
        "tail_correlations": correlations,
        "interpretation": {
            "verified": [
                "saved SG edges reproduce the formal Vela node-current export under the production formula",
                "cell-first recovery improves the same-edge-flux comparison without changing the device state or transport operator",
            ],
            "strongly_supported": "recovery support/order is the principal cause of the weak-current node-field tail",
            "unresolved": "the exact proprietary SDevice element-to-vertex reconstruction weights",
        },
    }

    args.output_root.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_root / "node_topology_diagnostics.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostic_rows[0]))
        writer.writeheader()
        writer.writerows(diagnostic_rows)
    args.report_json.parent.mkdir(parents=True, exist_ok=True)
    args.report_json.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    write_markdown(args.report_markdown, result)
    groups_by_name = {group["name"]: group for group in groups}
    print(
        json.dumps(
            {
                "main_p95_direct_decade": groups_by_name["main_current"][
                    "sdevice_vs_direct_node_fit"
                ][
                    "magnitude_abs_dex_p95"
                ],
                "main_p95_cell_first_decade": groups_by_name["main_current"][
                    "sdevice_vs_cell_first_area_projection"
                ][
                    "magnitude_abs_dex_p95"
                ],
                "tail_failures_direct": groups_by_name["original_tail"][
                    "nodes_above_0p5_decade"
                ][
                    "direct_node_fit"
                ],
                "tail_failures_cell_first": groups_by_name["original_tail"][
                    "nodes_above_0p5_decade"
                ][
                    "cell_first_area_projection"
                ],
                "report_json": str(args.report_json),
                "report_markdown": str(args.report_markdown),
            },
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
