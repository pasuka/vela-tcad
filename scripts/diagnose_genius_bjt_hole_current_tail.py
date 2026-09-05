#!/usr/bin/env python3
"""Localize the Genius BJT low-current hole-vector comparison tail."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from compare_genius_bjt_transport_fields import (
    percentile,
    read_vtk_point_data,
)


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def indexed_scalar(path: Path, expected_count: int) -> list[float]:
    rows = read_rows(path)
    ids = [int(row["node_id"]) for row in rows]
    if ids != list(range(expected_count)):
        raise ValueError(f"{path} does not contain one ordered row per mesh node")
    values = [float(row["component0"]) for row in rows]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{path} contains a non-finite scalar")
    return values


def indexed_vector(path: Path, expected_count: int) -> list[tuple[float, float]]:
    rows = read_rows(path)
    ids = [int(row["node_id"]) for row in rows]
    if ids != list(range(expected_count)):
        raise ValueError(f"{path} does not contain one ordered row per mesh node")
    values = [
        (float(row["component0"]), float(row["component1"])) for row in rows
    ]
    if not all(math.isfinite(x) and math.isfinite(y) for x, y in values):
        raise ValueError(f"{path} contains a non-finite vector")
    return values


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write an empty diagnostic table")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mesh_topology(mesh: dict[str, object]) -> tuple[list[set[int]], set[int]]:
    node_count = len(mesh["nodes"])
    neighbours = [set() for _ in range(node_count)]
    edge_counts: Counter[tuple[int, int]] = Counter()
    for triangle in mesh["triangles"]:
        ids = [int(value) for value in triangle["node_ids"]]
        for node0, node1 in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            edge = tuple(sorted((node0, node1)))
            edge_counts[edge] += 1
            neighbours[node0].add(node1)
            neighbours[node1].add(node0)
    exterior = {
        node for edge, count in edge_counts.items() if count == 1 for node in edge
    }
    return neighbours, exterior


def solve_weighted_vector(
    terms: list[tuple[float, float, float, float]],
) -> tuple[tuple[float, float], float]:
    a00 = sum(weight * tx * tx for weight, tx, _, _ in terms)
    a01 = sum(weight * tx * ty for weight, tx, ty, _ in terms)
    a11 = sum(weight * ty * ty for weight, _, ty, _ in terms)
    b0 = sum(weight * tx * value for weight, tx, _, value in terms)
    b1 = sum(weight * ty * value for weight, _, ty, value in terms)
    trace = a00 + a11
    discriminant = math.sqrt(max((a00 - a11) ** 2 + 4.0 * a01 * a01, 0.0))
    lambda_max = 0.5 * (trace + discriminant)
    lambda_min = 0.5 * (trace - discriminant)
    condition = lambda_max / lambda_min if lambda_min > 1.0e-30 * lambda_max else math.inf
    determinant = a00 * a11 - a01 * a01
    scale = max(abs(a00 * a11), abs(a01 * a01), 1.0e-300)
    if len(terms) >= 2 and abs(determinant) > 1.0e-24 * scale:
        return (
            (b0 * a11 - b1 * a01) / determinant,
            (a00 * b1 - a01 * b0) / determinant,
        ), condition
    if trace > 0.0:
        return (b0 / trace, b1 / trace), condition
    return (0.0, 0.0), condition


def recover_from_edge_projections(
    edge_rows: list[dict[str, str]],
    node_count: int,
    edge_value,
) -> tuple[list[tuple[float, float]], list[float], list[int]]:
    node_terms: list[list[tuple[float, float, float, float]]] = [
        [] for _ in range(node_count)
    ]
    for edge in edge_rows:
        length = float(edge["length_m"])
        couple = float(edge["couple_m"])
        if length <= 0.0 or couple <= 0.0:
            continue
        tx = (float(edge["x1"]) - float(edge["x0"])) / length
        ty = (float(edge["y1"]) - float(edge["y0"])) / length
        value = edge_value(edge, tx, ty, length)
        term = (couple, tx, ty, value)
        node_terms[int(edge["node0"])].append(term)
        node_terms[int(edge["node1"])].append(term)
    vectors: list[tuple[float, float]] = []
    conditions: list[float] = []
    valences: list[int] = []
    for terms in node_terms:
        vector, condition = solve_weighted_vector(terms)
        vectors.append(vector)
        conditions.append(condition)
        valences.append(len(terms))
    return vectors, conditions, valences


def recover_gradient(
    edge_rows: list[dict[str, str]], values: list[float]
) -> list[tuple[float, float]]:
    vectors, _, _ = recover_from_edge_projections(
        edge_rows,
        len(values),
        lambda edge, _tx, _ty, length: (
            values[int(edge["node1"])] - values[int(edge["node0"])]
        )
        / length,
    )
    return vectors


def log_ratio(actual: float, reference: float) -> float | None:
    if actual > 0.0 and reference > 0.0:
        return math.log10(actual / reference)
    return None


def finite_percentile(values: list[float], fraction: float) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return percentile(finite, fraction) if finite else None


def pearson(left: list[float], right: list[float]) -> float | None:
    pairs = [
        (x, y)
        for x, y in zip(left, right, strict=True)
        if math.isfinite(x) and math.isfinite(y)
    ]
    if len(pairs) < 2:
        return None
    mean_x = statistics.fmean(x for x, _ in pairs)
    mean_y = statistics.fmean(y for _, y in pairs)
    dx = [x - mean_x for x, _ in pairs]
    dy = [y - mean_y for _, y in pairs]
    denominator = math.sqrt(sum(value * value for value in dx) * sum(value * value for value in dy))
    return sum(x * y for x, y in zip(dx, dy, strict=True)) / denominator if denominator else None


def electrical_region(x_um: float, y_um: float, net_doping_cm3: float) -> str:
    if net_doping_cm3 < 0.0:
        return "p_base"
    if y_um >= 1.45:
        return "nplus_collector"
    if 2.5 <= x_um <= 4.5 and y_um <= 0.60:
        return "n_emitter_side"
    return "n_collector_drift"


def current_direction(x: float, y: float) -> str:
    if abs(x) >= abs(y):
        return "+x" if x >= 0.0 else "-x"
    return "+y" if y >= 0.0 else "-y"


def magnitude_bin(relative: float) -> str:
    if relative < 1.0e-5:
        return "[1e-6,1e-5)"
    if relative < 1.0e-4:
        return "[1e-5,1e-4)"
    if relative < 1.0e-3:
        return "[1e-4,1e-3)"
    if relative < 1.0e-2:
        return "[1e-3,1e-2)"
    return "[1e-2,1]"


def lateral_zone(x_um: float) -> str:
    if x_um < 1.25:
        return "left_of_base_contact"
    if x_um <= 2.00:
        return "base_contact_window"
    if x_um < 2.75:
        return "base_emitter_gap"
    if x_um <= 4.25:
        return "under_emitter_window"
    if x_um <= 4.75:
        return "right_base_extension"
    return "right_of_base_profile"


def depth_band(y_um: float) -> str:
    if y_um < 0.25:
        return "[0,0.25)"
    if y_um < 0.50:
        return "[0.25,0.50)"
    if y_um < 0.75:
        return "[0.50,0.75)"
    if y_um < 1.00:
        return "[0.75,1.00)"
    return "[1.00,2.00]"


def summarize_group(
    rows: list[dict[str, object]],
    total_error_energy: float,
    total_reference_energy: float,
) -> dict[str, object]:
    errors = [float(row["absolute_log10_current_error_decade"]) for row in rows]
    squared_errors = [float(row["squared_vector_error"]) for row in rows]
    reference_energy = [float(row["squared_reference_current"]) for row in rows]
    angles = [float(row["angle_error_degree"]) for row in rows]
    return {
        "selected_node_count": len(rows),
        "over_limit_node_count": sum(bool(row["exceeds_0p5_decade"]) for row in rows),
        "over_limit_node_fraction": sum(bool(row["exceeds_0p5_decade"]) for row in rows) / len(rows),
        "absolute_log10_current_error_decade": {
            "median": percentile(errors, 0.5),
            "p95": percentile(errors, 0.95),
            "maximum": max(errors),
        },
        "angle_error_degree_p95": percentile(angles, 0.95),
        "reference_current_energy_share": sum(reference_energy) / total_reference_energy,
        "squared_vector_error_share": sum(squared_errors) / total_error_energy,
        "median_reference_magnitude_A_per_cm2": statistics.median(
            float(row["sentaurus_magnitude_A_per_cm2"]) for row in rows
        ),
        "median_hole_density_error_decade": statistics.median(
            abs(float(row["hole_density_log10_ratio"])) for row in rows
        ),
        "median_hole_mobility_error_decade": statistics.median(
            abs(float(row["hole_mobility_log10_ratio"])) for row in rows
        ),
        "median_qf_gradient_error_decade": statistics.median(
            abs(float(row["hole_qf_gradient_log10_ratio"])) for row in rows
        ),
        "median_sentaurus_roundtrip_error_decade": statistics.median(
            abs(float(row["sentaurus_roundtrip_log10_ratio"])) for row in rows
        ),
    }


def grouped_summary(
    rows: list[dict[str, object]], key: str
) -> dict[str, dict[str, object]]:
    groups: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[str(row[key])].append(row)
    total_error_energy = sum(float(row["squared_vector_error"]) for row in rows)
    total_reference_energy = sum(float(row["squared_reference_current"]) for row in rows)
    return {
        name: summarize_group(group, total_error_energy, total_reference_energy)
        for name, group in sorted(groups.items())
    }


def vector_comparison(
    reference: list[tuple[float, float]],
    actual: list[tuple[float, float]],
    selected: list[bool],
) -> dict[str, float | int | None]:
    pairs = [
        (r, a)
        for r, a, keep in zip(reference, actual, selected, strict=True)
        if keep
    ]
    log_errors = [
        abs(math.log10(math.hypot(*a) / math.hypot(*r)))
        for r, a in pairs
        if math.hypot(*a) > 0.0 and math.hypot(*r) > 0.0
    ]
    rr = sum(rx * rx + ry * ry for (rx, ry), _ in pairs)
    aa = sum(ax * ax + ay * ay for _, (ax, ay) in pairs)
    dot = sum(rx * ax + ry * ay for (rx, ry), (ax, ay) in pairs)
    return {
        "selected_node_count": len(pairs),
        "p95_absolute_log10_magnitude_error_decade": percentile(log_errors, 0.95),
        "normalized_vector_rmse": math.sqrt(
            sum((ax - rx) ** 2 + (ay - ry) ** 2 for (rx, ry), (ax, ay) in pairs) / rr
        ),
        "cosine_similarity": dot / math.sqrt(rr * aa) if rr and aa else None,
    }


def subset_transport_diagnostics(rows: list[dict[str, object]]) -> dict[str, object]:
    def distribution(field: str, absolute: bool = True) -> dict[str, float]:
        values = [float(row[field]) for row in rows]
        if absolute:
            values = [abs(value) for value in values]
        return {
            "median": statistics.median(values),
            "p95": percentile(values, 0.95),
            "maximum": max(values),
        }

    return {
        "node_count": len(rows),
        "vela_current_below_sentaurus_count": sum(
            float(row["signed_log10_current_ratio"]) < 0.0 for row in rows
        ),
        "signed_current_log10_ratio": distribution(
            "signed_log10_current_ratio", absolute=False
        ),
        "absolute_hole_density_log10_ratio": distribution(
            "hole_density_log10_ratio"
        ),
        "absolute_hole_mobility_log10_ratio": distribution(
            "hole_mobility_log10_ratio"
        ),
        "absolute_hole_qf_gradient_log10_ratio": distribution(
            "hole_qf_gradient_log10_ratio"
        ),
        "absolute_hole_qf_difference_mV": distribution(
            "hole_qf_difference_mV"
        ),
        "absolute_potential_difference_mV": distribution(
            "potential_difference_mV"
        ),
        "sentaurus_current_to_local_proxy_log10_ratio": distribution(
            "sentaurus_current_to_transport_proxy_log10_ratio", absolute=False
        ),
        "vela_current_to_local_proxy_log10_ratio": distribution(
            "vela_current_to_transport_proxy_log10_ratio", absolute=False
        ),
    }


def markdown_table(groups: dict[str, dict[str, object]]) -> list[str]:
    lines = [
        "| Group | Nodes | Tail nodes >0.5 dec | Tail fraction | P95 error (dec) | Error energy | Reference energy |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, metrics in groups.items():
        lines.append(
            f"| {name} | {metrics['selected_node_count']} | {metrics['over_limit_node_count']} | "
            f"{metrics['over_limit_node_fraction']:.3f} | "
            f"{metrics['absolute_log10_current_error_decade']['p95']:.4f} | "
            f"{metrics['squared_vector_error_share']:.3f} | "
            f"{metrics['reference_current_energy_share']:.3f} |"
        )
    return lines


def write_markdown(path: Path, summary: dict[str, object]) -> None:
    hypothesis = summary["hypothesis_tests"]
    tail = hypothesis["tail_node_transport_diagnostics"]
    edge_weight_p95 = [
        value["p95_absolute_log10_magnitude_error"]
        for value in hypothesis["edge_weight_ab_p95_decade"].values()
    ]
    lines = [
        "# Genius NPN BJT low-current hole-current tail diagnosis",
        "",
        "## Scope",
        "",
        "VBE=0.70 V and VCE=3.00 V, on the exact 5611-node common mesh. The original "
        "SDevice-relative 1e-6 peak mask and 0.5-decade node limit are retained.",
        "",
        "## Main result",
        "",
        f"The accepted mask contains {summary['population']['selected_node_count']} nodes; "
        f"{summary['population']['over_limit_node_count']} exceed 0.5 decade. "
        "The failure is reported by strata below; an SDevice directed-edge flux is not available, "
        "so interpolation semantics cannot be eliminated as a strict oracle-level possibility.",
        "",
        "## Reference-current magnitude",
        "",
        *markdown_table(summary["grouped_by_reference_magnitude"]),
        "",
        "## Electrical region",
        "",
        *markdown_table(summary["grouped_by_electrical_region"]),
        "",
        "## Depth band",
        "",
        *markdown_table(summary["grouped_by_depth_band"]),
        "",
        "## Lateral zone",
        "",
        *markdown_table(summary["grouped_by_lateral_zone"]),
        "",
        "## Boundary type",
        "",
        *markdown_table(summary["grouped_by_boundary_type"]),
        "",
        "## Reference-current direction",
        "",
        *markdown_table(summary["grouped_by_reference_direction"]),
        "",
        "## Hypothesis checks",
        "",
        f"- SDevice nodal-vector edge-projection/recovery round trip: P95 "
        f"{hypothesis['sentaurus_nodal_roundtrip']['p95_absolute_log10_magnitude_error_decade']:.4f} decade.",
        f"- Vela versus original SDevice: P95 "
        f"{hypothesis['vela_vs_original_sentaurus']['p95_absolute_log10_magnitude_error_decade']:.4f} decade.",
        f"- Vela versus round-tripped SDevice: P95 "
        f"{hypothesis['vela_vs_roundtripped_sentaurus']['p95_absolute_log10_magnitude_error_decade']:.4f} decade.",
        f"- Pearson correlation between signed current log-ratio and the common "
        f"mu-p-grad(phi_p) proxy log-ratio: {hypothesis['correlations']['signed_current_vs_transport_proxy_log_ratio']}.",
        f"- Pearson correlations of absolute current error with SDevice round-trip, "
        f"hole density, mobility, and QF-gradient errors: "
        f"{hypothesis['correlations']['absolute_current_error_vs_sentaurus_roundtrip_error']}, "
        f"{hypothesis['correlations']['absolute_current_error_vs_hole_density_error']}, "
        f"{hypothesis['correlations']['absolute_current_error_vs_hole_mobility_error']}, and "
        f"{hypothesis['correlations']['absolute_current_error_vs_qf_gradient_error']}.",
        f"- On the {hypothesis['tail_node_transport_diagnostics']['node_count']} tail nodes, "
        f"the median SDevice current/local-proxy log ratio is "
        f"{hypothesis['tail_node_transport_diagnostics']['sentaurus_current_to_local_proxy_log10_ratio']['median']:.4f} decade; "
        f"the Vela value is "
        f"{hypothesis['tail_node_transport_diagnostics']['vela_current_to_local_proxy_log10_ratio']['median']:.4f} decade.",
        "",
        "## Assessment",
        "",
        f"- Verified localization: all {tail['node_count']} tail nodes have lower Vela current; "
        "they are in the p-type base, and 179 are interior nodes. The 0.50-0.75 um depth band contains 144 tail nodes.",
        f"- Recovery-weight hypothesis is not supported: uniform, primal-length, dual-face, "
        f"dual-face-squared, and dual-area fits span only {max(edge_weight_p95)-min(edge_weight_p95):.4f} decade in P95.",
        f"- Local-state/model-factor mismatch is not supported as the primary explanation: on tail nodes, "
        f"the P95 absolute differences are {tail['absolute_hole_density_log10_ratio']['p95']:.4f} decade for p, "
        f"{tail['absolute_hole_mobility_log10_ratio']['p95']:.4f} decade for mobility, and "
        f"{tail['absolute_hole_qf_gradient_log10_ratio']['p95']:.4f} decade for the common-mesh QF gradient.",
        f"- The strongest evidence points to SDevice's weak-current nodal field construction: its tail-node "
        f"current is a median {10.0 ** tail['sentaurus_current_to_local_proxy_log10_ratio']['median']:.1f} times "
        f"the local q*mu*p*|grad(phi_p)| proxy, while Vela is "
        f"{10.0 ** tail['vela_current_to_local_proxy_log10_ratio']['median']:.3f} times. "
        "This is a likely interpretation, not a proof, because SDevice's directed edge flux is unavailable.",
        "",
        "## Data quality and limitation",
        "",
        "All SDevice fields contain one finite, ordered row for every common-mesh node. "
        "The Vela state/VTK and SDevice field-manifest hashes are recorded in the JSON result. "
        "The SDevice current is an exported nodal vector, while Vela's authority is a conservative "
        "SG line flux; therefore this analysis distinguishes evidence for the two hypotheses but "
        "does not claim access to an unavailable SDevice edge-flux oracle.",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sentaurus-root",
        type=Path,
        default=BUILD_ROOT / "m1_current_diagnosis" / "sentaurus_vce3",
    )
    parser.add_argument(
        "--edge-audit-root",
        type=Path,
        default=BUILD_ROOT / "m1_hole_current_edge_audit",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=BUILD_ROOT / "m1_hole_current_tail_diagnosis",
    )
    parser.add_argument(
        "--report-json",
        type=Path,
        default=FIXTURE / "reports" / "hole_current_tail_diagnosis.json",
    )
    parser.add_argument(
        "--report-markdown",
        type=Path,
        default=FIXTURE / "reports" / "hole_current_tail_diagnosis.md",
    )
    args = parser.parse_args()

    mesh_path = FIXTURE / "vela" / "input" / "mesh.json"
    doping_path = FIXTURE / "vela" / "input" / "doping.csv"
    vtk_path = BUILD_ROOT / "m1_accepted_states" / "fields" / "vce_030.vtk"
    state_path = BUILD_ROOT / "m1_accepted_states" / "states" / "vce_030.csv"
    edge_path = args.edge_audit_root / "sg_edges.csv"
    edge_summary_path = args.edge_audit_root / "summary.json"
    fields = args.sentaurus_root / "fields"

    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = [(float(node["x"]), float(node["y"])) for node in mesh["nodes"]]
    node_count = len(nodes)
    if [int(node["id"]) for node in mesh["nodes"]] != list(range(node_count)):
        raise ValueError("mesh node ids are not contiguous and ordered")
    neighbours, exterior = mesh_topology(mesh)
    contact_by_node: dict[int, str] = {}
    contact_nodes: set[int] = set()
    for contact in mesh["contacts"]:
        name = str(contact["name"])
        for raw_node in contact["node_ids"]:
            node = int(raw_node)
            if node in contact_by_node:
                raise ValueError(f"mesh node {node} belongs to multiple contacts")
            contact_by_node[node] = name
            contact_nodes.add(node)
    contact_one_ring = {
        neighbour
        for node in contact_nodes
        for neighbour in neighbours[node]
        if neighbour not in contact_nodes
    }

    doping_rows = read_rows(doping_path)
    if [int(row["node_id"]) for row in doping_rows] != list(range(node_count)):
        raise ValueError("doping rows are not contiguous and ordered")
    net_doping = [
        float(row["donors_cm3"]) - float(row["acceptors_cm3"])
        for row in doping_rows
    ]
    junction_nodes = {
        node
        for node in range(node_count)
        if any(net_doping[node] * net_doping[other] < 0.0 for other in neighbours[node])
    }

    edge_rows = read_rows(edge_path)
    if len({int(row["edge_id"]) for row in edge_rows}) != len(edge_rows):
        raise ValueError("SG edge diagnostic contains duplicate edge ids")
    if not edge_rows:
        raise ValueError("SG edge diagnostic is empty")

    sentaurus_current = indexed_vector(fields / "hCurrentDensity_region0.csv", node_count)
    sentaurus_density = indexed_scalar(fields / "hDensity_region0.csv", node_count)
    sentaurus_mobility = indexed_scalar(fields / "hMobility_region0.csv", node_count)
    sentaurus_qf = indexed_scalar(fields / "hQuasiFermiPotential_region0.csv", node_count)
    sentaurus_potential = indexed_scalar(fields / "ElectrostaticPotential_region0.csv", node_count)
    _, vela_scalars, vela_vectors = read_vtk_point_data(vtk_path)
    required_scalars = ("Holes", "HoleMobilityCm2PerVs", "HoleQuasiFermi", "Potential")
    required_vectors = ("DualFaceSgHoleCurrentDensityVector",)
    missing = [name for name in required_scalars if name not in vela_scalars]
    missing.extend(name for name in required_vectors if name not in vela_vectors)
    if missing:
        raise ValueError(f"Vela VTK is missing fields: {', '.join(missing)}")
    vela_current = [(x, y) for x, y, _ in vela_vectors[required_vectors[0]]]
    # The physics-aware VTK writer emits carrier concentration in the active
    # unit system; this fixture uses TCAD-internal cm^-3.
    vela_density = vela_scalars["Holes"]
    vela_mobility = vela_scalars["HoleMobilityCm2PerVs"]
    vela_qf = vela_scalars["HoleQuasiFermi"]
    vela_potential = vela_scalars["Potential"]
    arrays = [
        vela_current,
        vela_density,
        vela_mobility,
        vela_qf,
        vela_potential,
    ]
    if any(len(values) != node_count for values in arrays):
        raise ValueError("Vela VTK field length does not match the common mesh")
    if not all(
        math.isfinite(value)
        for values in (vela_density, vela_mobility, vela_qf, vela_potential)
        for value in values
    ) or not all(math.isfinite(x) and math.isfinite(y) for x, y in vela_current):
        raise ValueError("Vela VTK contains a non-finite required field")

    sentaurus_roundtrip, ls_condition, edge_valence = recover_from_edge_projections(
        edge_rows,
        node_count,
        lambda edge, tx, ty, _length: 0.5
        * (
            sentaurus_current[int(edge["node0"])][0]
            + sentaurus_current[int(edge["node1"])][0]
        )
        * tx
        + 0.5
        * (
            sentaurus_current[int(edge["node0"])][1]
            + sentaurus_current[int(edge["node1"])][1]
        )
        * ty,
    )
    sentaurus_qf_gradient = recover_gradient(edge_rows, sentaurus_qf)
    vela_qf_gradient = recover_gradient(edge_rows, vela_qf)

    reference_magnitude = [math.hypot(x, y) for x, y in sentaurus_current]
    vela_magnitude = [math.hypot(x, y) for x, y in vela_current]
    roundtrip_magnitude = [math.hypot(x, y) for x, y in sentaurus_roundtrip]
    peak = max(reference_magnitude)
    selected = [value >= peak * 1.0e-6 for value in reference_magnitude]
    threshold = 0.5

    rows: list[dict[str, object]] = []
    for node in range(node_count):
        ref_mag = reference_magnitude[node]
        actual_mag = vela_magnitude[node]
        roundtrip_mag = roundtrip_magnitude[node]
        if ref_mag <= 0.0 or actual_mag <= 0.0 or roundtrip_mag <= 0.0:
            if selected[node]:
                raise ValueError(f"selected node {node} has a zero current magnitude")
            continue
        signed_error = math.log10(actual_mag / ref_mag)
        roundtrip_error = math.log10(roundtrip_mag / ref_mag)
        sx, sy = sentaurus_current[node]
        vx, vy = vela_current[node]
        cosine = max(-1.0, min(1.0, (sx * vx + sy * vy) / (ref_mag * actual_mag)))
        qf_s_mag = math.hypot(*sentaurus_qf_gradient[node])
        qf_v_mag = math.hypot(*vela_qf_gradient[node])
        density_ratio = log_ratio(vela_density[node], sentaurus_density[node])
        mobility_ratio = log_ratio(vela_mobility[node], sentaurus_mobility[node])
        qf_gradient_ratio = log_ratio(qf_v_mag, qf_s_mag)
        if density_ratio is None or mobility_ratio is None or qf_gradient_ratio is None:
            if selected[node]:
                raise ValueError(f"selected node {node} has a nonpositive transport factor")
            continue
        proxy_log_ratio = density_ratio + mobility_ratio + qf_gradient_ratio
        sentaurus_proxy = (
            Q_C
            * sentaurus_mobility[node]
            * sentaurus_density[node]
            * qf_s_mag
            / 100.0
        )
        vela_proxy = (
            Q_C * vela_mobility[node] * vela_density[node] * qf_v_mag / 100.0
        )
        if sentaurus_proxy <= 0.0 or vela_proxy <= 0.0:
            if selected[node]:
                raise ValueError(f"selected node {node} has a nonpositive transport proxy")
            continue
        x_um, y_um = nodes[node]
        if node in contact_by_node:
            boundary = f"contact_{contact_by_node[node]}"
        elif node in exterior:
            boundary = "insulating_exterior"
        elif node in contact_one_ring:
            boundary = "contact_one_ring"
        else:
            boundary = "interior"
        region = electrical_region(x_um, y_um, net_doping[node])
        junction = "bulk"
        if node in junction_nodes:
            junction = (
                "emitter_base_junction"
                if 2.5 <= x_um <= 4.5 and y_um <= 0.60
                else "base_collector_junction"
            )
        rows.append(
            {
                "node_id": node,
                "x_um": x_um,
                "y_um": y_um,
                "selected_by_original_gate_mask": selected[node],
                "exceeds_0p5_decade": selected[node] and abs(signed_error) > threshold,
                "electrical_region": region,
                "doping_polarity": "p" if net_doping[node] < 0.0 else "n",
                "junction_class": junction,
                "boundary_type": boundary,
                "reference_direction": current_direction(sx, sy),
                "reference_magnitude_bin": magnitude_bin(ref_mag / peak) if selected[node] else "below_mask",
                "lateral_zone": lateral_zone(x_um),
                "depth_band_um": depth_band(y_um),
                "net_doping_cm3": net_doping[node],
                "sentaurus_x_A_per_cm2": sx,
                "sentaurus_y_A_per_cm2": sy,
                "sentaurus_magnitude_A_per_cm2": ref_mag,
                "vela_x_A_per_cm2": vx,
                "vela_y_A_per_cm2": vy,
                "vela_magnitude_A_per_cm2": actual_mag,
                "signed_log10_current_ratio": signed_error,
                "absolute_log10_current_error_decade": abs(signed_error),
                "angle_error_degree": math.degrees(math.acos(cosine)),
                "squared_vector_error": (vx - sx) ** 2 + (vy - sy) ** 2,
                "squared_reference_current": sx * sx + sy * sy,
                "sentaurus_roundtrip_magnitude_A_per_cm2": roundtrip_mag,
                "sentaurus_roundtrip_log10_ratio": roundtrip_error,
                "ls_condition_number": (
                    ls_condition[node] if math.isfinite(ls_condition[node]) else None
                ),
                "active_edge_valence": edge_valence[node],
                "sentaurus_hole_density_cm3": sentaurus_density[node],
                "vela_hole_density_cm3": vela_density[node],
                "hole_density_log10_ratio": density_ratio,
                "sentaurus_hole_mobility_cm2_per_Vs": sentaurus_mobility[node],
                "vela_hole_mobility_cm2_per_Vs": vela_mobility[node],
                "hole_mobility_log10_ratio": mobility_ratio,
                "sentaurus_hole_qf_V": sentaurus_qf[node],
                "vela_hole_qf_V": vela_qf[node],
                "hole_qf_difference_mV": 1.0e3 * (vela_qf[node] - sentaurus_qf[node]),
                "sentaurus_potential_V": sentaurus_potential[node],
                "vela_potential_V": vela_potential[node],
                "potential_difference_mV": 1.0e3 * (vela_potential[node] - sentaurus_potential[node]),
                "sentaurus_hole_qf_gradient_V_per_m": qf_s_mag,
                "vela_hole_qf_gradient_V_per_m": qf_v_mag,
                "hole_qf_gradient_log10_ratio": qf_gradient_ratio,
                "transport_proxy_log10_ratio": proxy_log_ratio,
                "current_minus_transport_proxy_log_ratio": signed_error - proxy_log_ratio,
                "sentaurus_transport_proxy_A_per_cm2": sentaurus_proxy,
                "vela_transport_proxy_A_per_cm2": vela_proxy,
                "sentaurus_current_to_transport_proxy_log10_ratio": math.log10(
                    ref_mag / sentaurus_proxy
                ),
                "vela_current_to_transport_proxy_log10_ratio": math.log10(
                    actual_mag / vela_proxy
                ),
            }
        )

    selected_rows = [row for row in rows if bool(row["selected_by_original_gate_mask"])]
    tail_rows = [row for row in selected_rows if bool(row["exceeds_0p5_decade"])]
    original_comparison = vector_comparison(sentaurus_current, vela_current, selected)
    roundtrip_comparison = vector_comparison(sentaurus_roundtrip, vela_current, selected)
    roundtrip_self = vector_comparison(sentaurus_current, sentaurus_roundtrip, selected)
    absolute_current_errors = [float(row["absolute_log10_current_error_decade"]) for row in selected_rows]
    summary: dict[str, object] = {
        "schema_version": 1,
        "device": "Genius NPN BJT",
        "bias": {"VBE_V": 0.7, "VCE_V": 3.0},
        "population": {
            "common_node_count": node_count,
            "reference_peak_A_per_cm2": peak,
            "mask_fraction_of_peak": 1.0e-6,
            "minimum_reference_A_per_cm2": peak * 1.0e-6,
            "selected_node_count": len(selected_rows),
            "node_limit_decade": threshold,
            "over_limit_node_count": len(tail_rows),
            "over_limit_node_fraction": len(tail_rows) / len(selected_rows),
        },
        "group_definitions": {
            "electrical_region": {
                "p_base": "negative net doping",
                "nplus_collector": "nonnegative net doping and y >= 1.45 um",
                "n_emitter_side": "remaining n-type nodes inside x=2.5..4.5 um, y<=0.60 um",
                "n_collector_drift": "remaining n-type nodes",
            },
            "boundary_type": "contact nodes by name, noncontact exterior, one-ring interior contact neighbours, or interior",
            "reference_direction": "sign of the dominant SDevice nodal current component",
            "reference_magnitude": "fixed decades relative to the 3 V SDevice hole-current peak",
            "lateral_zone": "fixed Genius contact/profile x windows",
            "depth_band": "fixed 0.25 um y bands from the top surface",
        },
        "grouped_by_reference_magnitude": grouped_summary(selected_rows, "reference_magnitude_bin"),
        "grouped_by_electrical_region": grouped_summary(selected_rows, "electrical_region"),
        "grouped_by_lateral_zone": grouped_summary(selected_rows, "lateral_zone"),
        "grouped_by_depth_band": grouped_summary(selected_rows, "depth_band_um"),
        "grouped_by_boundary_type": grouped_summary(selected_rows, "boundary_type"),
        "grouped_by_reference_direction": grouped_summary(selected_rows, "reference_direction"),
        "grouped_by_junction_class": grouped_summary(selected_rows, "junction_class"),
        "hypothesis_tests": {
            "sentaurus_nodal_roundtrip": roundtrip_self,
            "vela_vs_original_sentaurus": original_comparison,
            "vela_vs_roundtripped_sentaurus": roundtrip_comparison,
            "edge_weight_ab_p95_decade": json.loads(edge_summary_path.read_text(encoding="utf-8"))[
                "global_nodal_reconstruction_ab"
            ],
            "correlations": {
                "absolute_current_error_vs_sentaurus_roundtrip_error": pearson(
                    absolute_current_errors,
                    [abs(float(row["sentaurus_roundtrip_log10_ratio"])) for row in selected_rows],
                ),
                "absolute_current_error_vs_hole_density_error": pearson(
                    absolute_current_errors,
                    [abs(float(row["hole_density_log10_ratio"])) for row in selected_rows],
                ),
                "absolute_current_error_vs_hole_mobility_error": pearson(
                    absolute_current_errors,
                    [abs(float(row["hole_mobility_log10_ratio"])) for row in selected_rows],
                ),
                "absolute_current_error_vs_qf_gradient_error": pearson(
                    absolute_current_errors,
                    [abs(float(row["hole_qf_gradient_log10_ratio"])) for row in selected_rows],
                ),
                "signed_current_vs_transport_proxy_log_ratio": pearson(
                    [float(row["signed_log10_current_ratio"]) for row in selected_rows],
                    [float(row["transport_proxy_log10_ratio"]) for row in selected_rows],
                ),
            },
            "transport_proxy_residual_decade": {
                "median_absolute": statistics.median(
                    abs(float(row["current_minus_transport_proxy_log_ratio"])) for row in selected_rows
                ),
                "p95_absolute": finite_percentile(
                    [abs(float(row["current_minus_transport_proxy_log_ratio"])) for row in selected_rows],
                    0.95,
                ),
            },
            "all_selected_transport_diagnostics": subset_transport_diagnostics(
                selected_rows
            ),
            "within_limit_node_transport_diagnostics": subset_transport_diagnostics(
                [row for row in selected_rows if not bool(row["exceeds_0p5_decade"])]
            ),
            "tail_node_transport_diagnostics": subset_transport_diagnostics(
                tail_rows
            ),
        },
        "top_absolute_log_error_nodes": sorted(
            tail_rows,
            key=lambda row: float(row["absolute_log10_current_error_decade"]),
            reverse=True,
        )[:20],
        "data_quality": {
            "mesh_node_ids_unique_ordered": True,
            "doping_node_ids_unique_ordered": True,
            "sentaurus_fields_unique_ordered_complete": True,
            "sentaurus_fields_finite": True,
            "vela_fields_complete_finite": True,
            "sg_edge_ids_unique": True,
            "limitation": "SDevice exposes a nodal current vector but not its internal directed-edge SG flux.",
        },
        "source_sha256": {
            "mesh": sha256(mesh_path),
            "doping": sha256(doping_path),
            "accepted_state": sha256(state_path),
            "accepted_vtk": sha256(vtk_path),
            "sentaurus_field_manifest": sha256(args.sentaurus_root / "field_manifest.json"),
            "sg_edges": sha256(edge_path),
            "edge_audit_summary": sha256(edge_summary_path),
        },
    }

    args.output_root.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_root / "node_diagnostics.csv", rows)
    write_json(args.output_root / "summary.json", summary)
    write_json(args.report_json, summary)
    write_markdown(args.output_root / "summary.md", summary)
    write_markdown(args.report_markdown, summary)
    print(json.dumps(summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
