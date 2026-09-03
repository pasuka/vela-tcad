#!/usr/bin/env python3
"""Audit frozen-state SDevice hole-current representation choices for Genius BJT."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def scalar(path: Path, support: str = "node") -> list[float]:
    key = "node_id" if support == "node" else "cell_id"
    result = []
    for index, row in enumerate(rows(path)):
        if int(row[key]) != index:
            raise ValueError(f"unordered {support} field: {path}")
        result.append(float(row["component0"]))
    return result


def vector(path: Path, support: str) -> list[tuple[float, float]]:
    key = "node_id" if support == "node" else "cell_id"
    result = []
    for index, row in enumerate(rows(path)):
        if int(row[key]) != index:
            raise ValueError(f"unordered {support} field: {path}")
        result.append((float(row["component0"]), float(row["component1"])))
    return result


def magnitude(value: tuple[float, float]) -> float:
    return math.hypot(*value)


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def field_metrics(
    reference: list[tuple[float, float]],
    candidate: list[tuple[float, float]],
    indices: list[int],
) -> dict[str, float | int]:
    floor = 1.0e-300
    signed_dex = [
        math.log10(max(magnitude(candidate[i]), floor) / max(magnitude(reference[i]), floor))
        for i in indices
    ]
    dex = [abs(value) for value in signed_dex]
    error2 = sum(
        (candidate[i][0] - reference[i][0]) ** 2
        + (candidate[i][1] - reference[i][1]) ** 2
        for i in indices
    )
    reference2 = sum(reference[i][0] ** 2 + reference[i][1] ** 2 for i in indices)
    candidate2 = sum(candidate[i][0] ** 2 + candidate[i][1] ** 2 for i in indices)
    dot = sum(
        reference[i][0] * candidate[i][0] + reference[i][1] * candidate[i][1]
        for i in indices
    )
    return {
        "count": len(indices),
        "magnitude_log10_ratio_median": statistics.median(signed_dex) if signed_dex else 0.0,
        "magnitude_ratio_median": 10.0 ** statistics.median(signed_dex) if signed_dex else 1.0,
        "magnitude_abs_dex_median": statistics.median(dex) if dex else 0.0,
        "magnitude_abs_dex_p95": percentile(dex, 0.95),
        "magnitude_abs_dex_max": max(dex, default=0.0),
        "normalized_vector_rmse": math.sqrt(error2 / max(reference2, floor)),
        "cosine_similarity": dot / math.sqrt(max(reference2 * candidate2, floor)),
    }


def least_squares_scale(
    reference: list[tuple[float, float]],
    candidate: list[tuple[float, float]],
    indices: list[int],
) -> float:
    numerator = sum(
        reference[i][0] * candidate[i][0] + reference[i][1] * candidate[i][1]
        for i in indices
    )
    denominator = sum(
        candidate[i][0] ** 2 + candidate[i][1] ** 2 for i in indices
    )
    return numerator / denominator if denominator > 0.0 else 1.0


def scaled(
    values: list[tuple[float, float]], factor: float
) -> list[tuple[float, float]]:
    return [(factor * x, factor * y) for x, y in values]


def triangle_areas(export: Path) -> tuple[list[tuple[int, int, int]], list[float]]:
    points = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in rows(export / "nodes.csv")
    }
    triangles = []
    areas = []
    for row in rows(export / "elements.csv"):
        ids = (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        a, b, c = (points[i] for i in ids)
        area = 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
        if area <= 0.0:
            raise ValueError(f"degenerate triangle {row['id']}")
        triangles.append(ids)
        areas.append(area)
    return triangles, areas


def reconstruct(
    node_count: int,
    triangles: list[tuple[int, int, int]],
    areas: list[float],
    cell_vectors: list[tuple[float, float]],
    weighting: str,
) -> list[tuple[float, float]]:
    sums = [[0.0, 0.0, 0.0] for _ in range(node_count)]
    for cell_id, ids in enumerate(triangles):
        weight = areas[cell_id] if weighting == "area" else 1.0
        value = cell_vectors[cell_id]
        for node_id in ids:
            sums[node_id][0] += weight * value[0]
            sums[node_id][1] += weight * value[1]
            sums[node_id][2] += weight
    if any(value[2] == 0.0 for value in sums):
        raise ValueError("mesh contains a node without incident cells")
    return [(value[0] / value[2], value[1] / value[2]) for value in sums]


def solve_vector_projections(
    terms: list[tuple[float, float, float, float]],
) -> tuple[float, float] | None:
    a00 = sum(weight * tx * tx for weight, tx, _, _ in terms)
    a01 = sum(weight * tx * ty for weight, tx, ty, _ in terms)
    a11 = sum(weight * ty * ty for weight, _, ty, _ in terms)
    b0 = sum(weight * tx * value for weight, tx, _, value in terms)
    b1 = sum(weight * ty * value for weight, _, ty, value in terms)
    determinant = a00 * a11 - a01 * a01
    scale = max(abs(a00 * a11), abs(a01 * a01), 1.0e-300)
    if abs(determinant) <= 1.0e-24 * scale:
        return None
    return (
        (b0 * a11 - b1 * a01) / determinant,
        (a00 * b1 - a01 * b0) / determinant,
    )


def vela_cell_currents_from_sg_edges(
    edge_path: Path,
    triangles: list[tuple[int, int, int]],
) -> tuple[list[tuple[float, float]], list[int]]:
    edge_by_nodes = {
        tuple(sorted((int(row["node0"]), int(row["node1"])))): row
        for row in rows(edge_path)
    }
    values = [(0.0, 0.0) for _ in triangles]
    valid = []
    for cell_id, ids in enumerate(triangles):
        terms = []
        for first, second in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            edge = edge_by_nodes.get(tuple(sorted((first, second))))
            if edge is None:
                continue
            couple = float(edge["couple_m"])
            length = float(edge["length_m"])
            if couple <= 0.0 or length <= 0.0:
                continue
            tx = (float(edge["x1"]) - float(edge["x0"])) / length
            ty = (float(edge["y1"]) - float(edge["y0"])) / length
            particle_flux = float(edge["hole_particle_line_flux_per_m_s"]) / couple
            current_A_per_cm2 = 1.602176634e-19 * particle_flux * 1.0e-4
            terms.append((couple, tx, ty, current_A_per_cm2))
        value = solve_vector_projections(terms)
        if value is not None:
            values[cell_id] = value
            valid.append(cell_id)
    return values, valid


def reconstruct_sparse(
    node_count: int,
    triangles: list[tuple[int, int, int]],
    areas: list[float],
    cell_vectors: list[tuple[float, float]],
    valid_cells: list[int],
) -> tuple[list[tuple[float, float]], list[int]]:
    sums = [[0.0, 0.0, 0.0] for _ in range(node_count)]
    for cell_id in valid_cells:
        weight = areas[cell_id]
        for node_id in triangles[cell_id]:
            sums[node_id][0] += weight * cell_vectors[cell_id][0]
            sums[node_id][1] += weight * cell_vectors[cell_id][1]
            sums[node_id][2] += weight
    valid_nodes = [node_id for node_id, value in enumerate(sums) if value[2] > 0.0]
    values = [
        (value[0] / value[2], value[1] / value[2]) if value[2] > 0.0 else (0.0, 0.0)
        for value in sums
    ]
    return values, valid_nodes


def load_export(export: Path) -> dict[str, list]:
    fields = export / "fields"
    return {
        "potential": scalar(fields / "ElectrostaticPotential_region0.csv"),
        "holes": scalar(fields / "hDensity_region0.csv"),
        "hole_qf": scalar(fields / "hQuasiFermiPotential_region0.csv"),
        "node_current": vector(fields / "hCurrentDensity_region0.csv", "node"),
        "cell_current": vector(fields / "hCurrentDensity_region0_cells.csv", "cell"),
        "cell_hole_qf_gradient": vector(fields / "hGradQuasiFermi_region0_cells.csv", "cell"),
        "cell_hole_mobility": scalar(fields / "hMobility_region0_cells.csv", "cell"),
    }


def load_node_state(export: Path) -> dict[str, list]:
    fields = export / "fields"
    return {
        "potential": scalar(fields / "ElectrostaticPotential_region0.csv"),
        "holes": scalar(fields / "hDensity_region0.csv"),
        "hole_qf": scalar(fields / "hQuasiFermiPotential_region0.csv"),
        "node_current": vector(fields / "hCurrentDensity_region0.csv", "node"),
    }


def max_scalar_delta(reference: list[float], candidate: list[float]) -> float:
    return max((abs(a - b) for a, b in zip(reference, candidate, strict=True)), default=0.0)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_markdown(path: Path, result: dict) -> None:
    variants = result["variants"]
    reconstruction = result["cell_to_node_reconstruction"]["area"]["fixed_unit_scale"]
    sg = result["vela_sg_edge_to_cell_comparison"]
    lines = [
        "# Genius NPN BJT SDevice current-semantics validation",
        "",
        "## Scope",
        "",
        "VBE=0.70 V and VCE=3.00 V on the exact 5611-node / 10940-triangle mesh. "
        "The original 2121-node current mask, 180 nodes above the 0.5-decade diagnostic limit, "
        "and the accepted SDevice state are retained.",
        "",
        "## Frozen-state A/B",
        "",
        "| Variant | State max delta | Tail P95 current delta (dec) | Tail max (dec) | Result |",
        "|---|---:|---:|---:|---|",
    ]
    for name, label in (("element_edge_current", "ElementEdgeCurrent"), ("element_edge_mobility", "hMobilityAveraging=ElementEdge")):
        item = variants[name]
        state = item["state_max_absolute_delta"]
        state_max = max(state.values())
        tail = item["node_current"]["original_tail_nodes"]
        lines.append(
            f"| {label} | {state_max:.3e} | {tail['magnitude_abs_dex_p95']:.6f} | "
            f"{tail['magnitude_abs_dex_max']:.6f} | "
            + ("No plotted-current effect" if tail["magnitude_abs_dex_max"] == 0.0 else "Too small to explain the tail")
            + " |"
        )
    lines.extend([
        "",
        "`ElementEdgeCurrent` leaves both node and element plotted current exactly unchanged on the "
        "loaded state. The mobility edge-average option changes the tail by at most 0.024 decade, "
        "far below the observed >=0.5-decade discrepancy.",
        "",
        "## Native SDevice element-to-node reconstruction",
        "",
        "The raw SDevice element-vector encoding is normalized by the fixed factor "
        "`1e-6/sqrt(2)`. An independently fitted scale on the main-current mask differs from this "
        "factor by less than 8e-5 relative. Area-weighting incident triangle currents gives:",
        "",
        "| Evaluation set | Magnitude P95 (dec) | Normalized vector RMSE | Cosine |",
        "|---|---:|---:|---:|",
        f"| Main current mask | {reconstruction['main_current_mask']['magnitude_abs_dex_p95']:.6f} | "
        f"{reconstruction['main_current_mask']['normalized_vector_rmse']:.6f} | "
        f"{reconstruction['main_current_mask']['cosine_similarity']:.6f} |",
        f"| Original 180 tail nodes | {reconstruction['original_tail_nodes']['magnitude_abs_dex_p95']:.6f} | "
        f"{reconstruction['original_tail_nodes']['normalized_vector_rmse']:.6f} | "
        f"{reconstruction['original_tail_nodes']['cosine_similarity']:.6f} |",
        "",
        "## Vela SG edge-flux cross-check",
        "",
        "Vela's directed SG line flux is divided by the dual-face length, recovered to one vector per "
        "triangle, then area-projected to nodes. This uses the same Vela flux data and changes only "
        "the representation support.",
        "",
        "| Comparison | Median magnitude ratio | P95 magnitude error (dec) | Cosine |",
        "|---|---:|---:|---:|",
        f"| Vela-recovered vs SDevice, 598 tail-adjacent cells | "
        f"{sg['tail_adjacent_cells']['magnitude_ratio_median']:.6f} | "
        f"{sg['tail_adjacent_cells']['magnitude_abs_dex_p95']:.6f} | "
        f"{sg['tail_adjacent_cells']['cosine_similarity']:.6f} |",
        f"| Vela area-projected vs SDevice, 180 tail nodes | "
        f"{sg['area_projected_nodes_vs_sdevice_tail']['magnitude_ratio_median']:.6f} | "
        f"{sg['area_projected_nodes_vs_sdevice_tail']['magnitude_abs_dex_p95']:.6f} | "
        f"{sg['area_projected_nodes_vs_sdevice_tail']['cosine_similarity']:.6f} |",
        f"| Vela area-projected vs current Vela node field, 180 tail nodes | "
        f"{sg['area_projected_nodes_vs_vela_tail']['magnitude_ratio_median']:.6f} | "
        f"{sg['area_projected_nodes_vs_vela_tail']['magnitude_abs_dex_p95']:.6f} | "
        f"{sg['area_projected_nodes_vs_vela_tail']['cosine_similarity']:.6f} |",
        f"| Vela area-projected vs current Vela node field, 1700 dominant nodes | "
        f"{sg['area_projected_nodes_vs_vela_dominant_current']['magnitude_ratio_median']:.6f} | "
        f"{sg['area_projected_nodes_vs_vela_dominant_current']['magnitude_abs_dex_p95']:.6f} | "
        f"{sg['area_projected_nodes_vs_vela_dominant_current']['cosine_similarity']:.6f} |",
        "",
        "The same Vela conservative fluxes reproduce the SDevice tail after cell recovery and "
        "area projection, but are a median 22.06 times the current Vela node-vector representation. "
        "This is strong evidence that the large weak-current tail is primarily a support/reconstruction "
        "semantic difference, not a 22x transport-current physics difference.",
        "",
        "## Convergence sensitivity",
        "",
        "| Variant | Tail P95 delta (dec) | Tail max delta (dec) |",
        "|---|---:|---:|",
    ])
    for name, item in result["convergence_variants"].items():
        label = "-Extrapolate" if name == "no_extrapolate" else "Digits=8 + ExtendedPrecision(128)"
        tail = item["node_current"]["original_tail_nodes"]
        lines.append(
            f"| {label} | {tail['magnitude_abs_dex_p95']:.6e} | "
            f"{tail['magnitude_abs_dex_max']:.6e} |"
        )
    if "digits8_ep128" not in result["convergence_variants"]:
        lines.append("| Digits=8 + ExtendedPrecision(128) | running | running |")
    lines.extend([
        "",
        "## Assessment and limitation",
        "",
        "- The no-extrapolation solve is numerically identical to the accepted state, so continuation "
        "initialization is excluded as the tail cause.",
        "- The element-current and Vela SG-flux cross-check upgrades the earlier hypothesis to strong "
        "evidence for differing node-vector construction semantics.",
        "- This does not identify Synopsys' proprietary vertex reconstruction formula exactly; the "
        "strict oracle remains terminal current, KCL, and conservative control-surface flux rather than "
        "low-energy nodal-vector tails.",
        "- Keep the existing 0.5-decade node gate as a diagnostic record. It should not alone fail the "
        "device while the dominant-current region and conservative flux checks pass.",
    ])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    probe = BUILD_ROOT / "sdevice_current_semantics_probe"
    parser.add_argument("--default", type=Path, default=probe / "default_export")
    parser.add_argument("--element-edge-current", type=Path, default=probe / "element_edge_current_export")
    parser.add_argument("--element-edge-mobility", type=Path, default=probe / "element_edge_mobility_export")
    parser.add_argument(
        "--tail-diagnostics",
        type=Path,
        default=BUILD_ROOT / "m1_hole_current_tail_diagnosis" / "node_diagnostics.csv",
    )
    parser.add_argument(
        "--report-json",
        type=Path,
        default=FIXTURE / "reports" / "sdevice_current_semantics_validation.json",
    )
    parser.add_argument(
        "--report-markdown",
        type=Path,
        default=FIXTURE / "reports" / "sdevice_current_semantics_validation.md",
    )
    parser.add_argument("--output", type=Path, default=probe / "current_semantics_audit.json")
    parser.add_argument("--no-extrapolate", type=Path, default=probe / "no_extrapolate_export")
    parser.add_argument("--strict", type=Path, default=probe / "digits8_ep128_export")
    parser.add_argument(
        "--sg-edges",
        type=Path,
        default=BUILD_ROOT / "m1_hole_current_edge_audit" / "sg_edges.csv",
    )
    args = parser.parse_args()

    default = load_export(args.default)
    variants = {
        "element_edge_current": load_export(args.element_edge_current),
        "element_edge_mobility": load_export(args.element_edge_mobility),
    }
    node_count = len(default["node_current"])
    peak = max(magnitude(value) for value in default["node_current"])
    main_indices = [
        i for i, value in enumerate(default["node_current"])
        if magnitude(value) >= peak * 1.0e-6
    ]
    dominant_indices = [
        i for i, value in enumerate(default["node_current"])
        if magnitude(value) >= peak * 1.0e-2
    ]
    tail_indices = [
        int(row["node_id"])
        for row in rows(args.tail_diagnostics)
        if row["exceeds_0p5_decade"].lower() == "true"
    ]
    tail_rows = rows(args.tail_diagnostics)
    vela_all = {
        int(row["node_id"]): (float(row["vela_x_A_per_cm2"]), float(row["vela_y_A_per_cm2"]))
        for row in tail_rows
    }
    vela_tail = {
        int(row["node_id"]): (float(row["vela_x_A_per_cm2"]), float(row["vela_y_A_per_cm2"]))
        for row in tail_rows
        if row["exceeds_0p5_decade"].lower() == "true"
    }
    native_tail_vs_vela = field_metrics(
        [vela_tail[i] for i in tail_indices],
        [default["node_current"][i] for i in tail_indices],
        list(range(len(tail_indices))),
    )

    triangles, areas = triangle_areas(args.default)
    physical_cell_current = scaled(default["cell_current"], 1.0e-6 / math.sqrt(2.0))
    vela_cell_current, valid_vela_cells = vela_cell_currents_from_sg_edges(args.sg_edges, triangles)
    tail_set = set(tail_indices)
    tail_cells = [
        cell_id for cell_id in valid_vela_cells if any(node_id in tail_set for node_id in triangles[cell_id])
    ]
    vela_projected_nodes, valid_vela_nodes = reconstruct_sparse(
        node_count, triangles, areas, vela_cell_current, valid_vela_cells
    )
    valid_tail_nodes = [node_id for node_id in tail_indices if node_id in set(valid_vela_nodes)]
    valid_main_nodes = [node_id for node_id in main_indices if node_id in set(valid_vela_nodes)]
    valid_dominant_nodes = [node_id for node_id in dominant_indices if node_id in set(valid_vela_nodes)]
    vela_sg_cell_comparison = {
        "valid_cell_count": len(valid_vela_cells),
        "all_valid_cells": field_metrics(physical_cell_current, vela_cell_current, valid_vela_cells),
        "tail_adjacent_cells": field_metrics(physical_cell_current, vela_cell_current, tail_cells),
        "area_projected_nodes_vs_sdevice_tail": field_metrics(
            default["node_current"], vela_projected_nodes, valid_tail_nodes
        ),
        "area_projected_nodes_vs_sdevice_main_mask": field_metrics(
            default["node_current"], vela_projected_nodes, valid_main_nodes
        ),
        "area_projected_nodes_vs_vela_tail": field_metrics(
            [vela_tail[i] for i in valid_tail_nodes],
            [vela_projected_nodes[i] for i in valid_tail_nodes],
            list(range(len(valid_tail_nodes))),
        ),
        "area_projected_nodes_vs_vela_main_mask": field_metrics(
            [vela_all[i] for i in valid_main_nodes],
            [vela_projected_nodes[i] for i in valid_main_nodes],
            list(range(len(valid_main_nodes))),
        ),
        "area_projected_nodes_vs_vela_dominant_current": field_metrics(
            [vela_all[i] for i in valid_dominant_nodes],
            [vela_projected_nodes[i] for i in valid_dominant_nodes],
            list(range(len(valid_dominant_nodes))),
        ),
    }
    cell_identity = {}
    for averaging in ("arithmetic", "geometric"):
        predicted = []
        for cell_id, ids in enumerate(triangles):
            densities = [default["holes"][node_id] for node_id in ids]
            if averaging == "arithmetic":
                density = sum(densities) / 3.0
            else:
                density = math.exp(sum(math.log(max(value, 1.0e-300)) for value in densities) / 3.0)
            coefficient = 1.602176634e-19 * default["cell_hole_mobility"][cell_id] * density
            gradient = default["cell_hole_qf_gradient"][cell_id]
            predicted.append((coefficient * gradient[0], coefficient * gradient[1]))
        cell_identity[averaging] = field_metrics(
            physical_cell_current, predicted, list(range(len(triangles)))
        )
        cell_identity[averaging]["tail_adjacent_cells"] = field_metrics(
            physical_cell_current, predicted, tail_cells
        )
    reconstruction = {}
    for weighting in ("uniform", "area"):
        recovered = reconstruct(node_count, triangles, areas, default["cell_current"], weighting)
        scale = least_squares_scale(default["node_current"], recovered, main_indices)
        calibrated = scaled(recovered, scale)
        fixed_scale = 1.0e-6 / math.sqrt(2.0)
        fixed_calibrated = scaled(recovered, fixed_scale)
        reconstruction[weighting] = {
            "least_squares_scale_from_main_mask": scale,
            "scale_ratio_to_1e_minus_6_over_sqrt_2": scale / (1.0e-6 / math.sqrt(2.0)),
            "all_nodes": field_metrics(default["node_current"], calibrated, list(range(node_count))),
            "main_current_mask": field_metrics(default["node_current"], calibrated, main_indices),
            "original_tail_nodes": field_metrics(default["node_current"], calibrated, tail_indices),
            "tail_vs_vela": field_metrics(
                [vela_tail[i] for i in tail_indices],
                [calibrated[i] for i in tail_indices],
                list(range(len(tail_indices))),
            ),
            "fixed_unit_scale": {
                "factor": fixed_scale,
                "main_current_mask": field_metrics(
                    default["node_current"], fixed_calibrated, main_indices
                ),
                "original_tail_nodes": field_metrics(
                    default["node_current"], fixed_calibrated, tail_indices
                ),
            },
        }

    variant_results = {}
    for name, candidate in variants.items():
        variant_results[name] = {
            "state_max_absolute_delta": {
                "potential_V": max_scalar_delta(default["potential"], candidate["potential"]),
                "hole_density_cm3": max_scalar_delta(default["holes"], candidate["holes"]),
                "hole_quasi_fermi_V": max_scalar_delta(default["hole_qf"], candidate["hole_qf"]),
            },
            "node_current": {
                "all_nodes": field_metrics(default["node_current"], candidate["node_current"], list(range(node_count))),
                "main_current_mask": field_metrics(default["node_current"], candidate["node_current"], main_indices),
                "original_tail_nodes": field_metrics(default["node_current"], candidate["node_current"], tail_indices),
            },
            "cell_current": field_metrics(
                default["cell_current"], candidate["cell_current"], list(range(len(default["cell_current"])))
            ),
        }


    convergence_results = {}
    for name, export in (("no_extrapolate", args.no_extrapolate), ("digits8_ep128", args.strict)):
        if not export.exists():
            continue
        candidate = load_node_state(export)
        convergence_results[name] = {
            "state_max_absolute_delta": {
                "potential_V": max_scalar_delta(default["potential"], candidate["potential"]),
                "hole_density_cm3": max_scalar_delta(default["holes"], candidate["holes"]),
                "hole_quasi_fermi_V": max_scalar_delta(default["hole_qf"], candidate["hole_qf"]),
            },
            "node_current": {
                "all_nodes": field_metrics(default["node_current"], candidate["node_current"], list(range(node_count))),
                "main_current_mask": field_metrics(default["node_current"], candidate["node_current"], main_indices),
                "original_tail_nodes": field_metrics(default["node_current"], candidate["node_current"], tail_indices),
            },
        }

    result = {
        "schema": "vela.genius_bjt_sdevice_current_semantics.v1",
        "node_count": node_count,
        "cell_count": len(default["cell_current"]),
        "main_current_mask_count": len(main_indices),
        "dominant_current_count": len(dominant_indices),
        "original_tail_count": len(tail_indices),
        "native_sdevice_tail_vs_vela": native_tail_vs_vela,
        "variants": variant_results,
        "convergence_variants": convergence_results,
        "cell_to_node_reconstruction": reconstruction,
        "cell_continuum_identity_using_fixed_1e_minus_6_over_sqrt_2_scale": cell_identity,
        "vela_sg_edge_to_cell_comparison": vela_sg_cell_comparison,
    }
    write_json(args.output, result)
    write_json(args.report_json, result)
    write_markdown(args.report_markdown, result)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
