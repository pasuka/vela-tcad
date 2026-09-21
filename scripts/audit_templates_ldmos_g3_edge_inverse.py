#!/usr/bin/env python3
"""WP3/T2 conservative edge-correction inverse for Templates/LDMOS G3.

This tool consumes only the exact eight-state T1 replay package.  It builds the
R4 canonical incidence matrix on the frozen seven free-silicon continuity rows
and their incident active silicon transport edges, solves ``B delta = -r`` with
three weighted minimum-norm contracts, and evaluates one-parameter edge-family
models with deterministic label permutation tests.

``delta`` is a required feasible correction that would make the Sentaurus state
a zero of the Vela operator.  It is not a measured Sentaurus/Vela edge-flux
difference, and this script cannot by itself establish a root cause.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np

try:
    from scripts.audit_templates_ldmos_g3_residual_scaling import (
        DEFAULT_NODES,
        EXPECTED_BIASES,
        current_scale_from_sg,
    )
    from scripts.audit_templates_ldmos_interface_pair_box import region_local_geometry
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_residual_scaling import (  # type: ignore
        DEFAULT_NODES,
        EXPECTED_BIASES,
        current_scale_from_sg,
    )
    from audit_templates_ldmos_interface_pair_box import region_local_geometry  # type: ignore


T1_SUMMARY_SCHEMA = "vela.templates_ldmos.g3_wp3_t1_residual_scaling.v1"
T1_MANIFEST_SCHEMA = "vela.templates_ldmos.g3_wp3_t1_state_manifest.v1"
OUTPUT_SCHEMA = "vela.templates_ldmos.g3_wp3_t2_edge_inverse.v1"
REGULARIZERS = ("unweighted_l2", "abs_phi_weighted", "couple_over_length_weighted")
FAMILY_NAMES = (
    "bernoulli_abs_eta",
    "interface_normal_fraction",
    "interface_parallel_fraction",
    "obtuse_incident_fraction",
    "doping_gradient_abs",
    "hfs_drive_abs",
)
ENDPOINTS = (1.0, 7.0 / 6.0)
ELEMENTARY_CHARGE_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite(value: str, field: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite {field}")
    return result


def percentile(values: Iterable[float], fraction: float) -> float:
    data = sorted(values)
    if not data:
        raise ValueError("percentile of empty data")
    position = fraction * (len(data) - 1)
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return data[lower]
    weight = position - lower
    return data[lower] * (1.0 - weight) + data[upper] * weight


def norm(values: np.ndarray) -> float:
    return float(np.linalg.norm(values))


def inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def load_t1_contract(root: Path) -> tuple[dict[str, Any], dict[str, Any], Path]:
    root = root.resolve()
    summary_path = root / "summary.json"
    manifest_path = root / "state_manifest.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if summary.get("schema") != T1_SUMMARY_SCHEMA:
        raise ValueError("T1 summary schema mismatch")
    if manifest.get("schema") != T1_MANIFEST_SCHEMA:
        raise ValueError("T1 manifest schema mismatch")
    contract = summary.get("contract", {})
    if contract.get("state_count") != 8:
        raise ValueError("T1 contract must contain exactly eight states")
    biases = [float(value) for value in contract.get("bias_points_V", [])]
    if len(biases) != 8 or any(
        not math.isclose(actual, expected, rel_tol=0.0, abs_tol=5.0e-7)
        for actual, expected in zip(biases, EXPECTED_BIASES, strict=True)
    ):
        raise ValueError("T1 bias contract is not 0:1/6:7/6")
    if tuple(contract.get("selected_nodes", [])) != tuple(DEFAULT_NODES):
        raise ValueError("T1 selected-node contract mismatch")
    if contract.get("drain_bias_V") != 0.1:
        raise ValueError("T1 drain-bias contract mismatch")
    forbidden = (
        contract.get("production_defaults_changed") is not False
        or contract.get("idvd_bv_full_physics_states_used") is not False
        or contract.get("ialmob") != "disabled"
        or contract.get("predictor") != "disabled"
    )
    if forbidden:
        raise ValueError("T1 physical/freeze contract mismatch")
    states = manifest.get("states", [])
    if len(states) != 8:
        raise ValueError("T1 manifest must contain eight states")
    mesh = Path(manifest["mesh"]).resolve()
    if not mesh.is_file() or sha256(mesh) != manifest.get("mesh_sha256"):
        raise ValueError("T1 mesh hash mismatch")
    for state, expected in zip(states, EXPECTED_BIASES, strict=True):
        if not math.isclose(float(state["bias_V"]), expected, rel_tol=0.0, abs_tol=5e-7):
            raise ValueError("T1 manifest state order/bias mismatch")
        for field, hash_field in (
            ("state_path", "state_sha256"),
            ("carrier_config", "carrier_config_sha256"),
            ("sg_config", "sg_config_sha256"),
        ):
            path = Path(state[field]).resolve()
            if not path.is_file() or sha256(path) != state[hash_field]:
                raise ValueError(f"T1 manifest hash mismatch: {field} at {expected}")
        for field in ("carrier_csv", "sg_csv"):
            path = Path(state[field]).resolve()
            if not path.is_file() or not inside(path, root):
                raise ValueError(f"T1 generated table escapes package: {field}")
    return summary, manifest, mesh


def build_incidence(
    node_ids: tuple[int, ...], edge_pairs: list[tuple[int, int]]
) -> np.ndarray:
    if len(set(node_ids)) != len(node_ids):
        raise ValueError("duplicate graph nodes")
    node_index = {node: index for index, node in enumerate(node_ids)}
    matrix = np.zeros((len(node_ids), len(edge_pairs)), dtype=float)
    seen: set[tuple[int, int]] = set()
    for column, pair in enumerate(edge_pairs):
        tail, head = pair
        if tail >= head:
            raise ValueError("edges must use canonical low-ID to high-ID orientation")
        if pair in seen:
            raise ValueError(f"duplicate canonical edge {pair}")
        seen.add(pair)
        if tail not in node_index and head not in node_index:
            raise ValueError("edge is not incident on graph domain")
        if tail in node_index:
            matrix[node_index[tail], column] = 1.0
        if head in node_index:
            matrix[node_index[head], column] = -1.0
    return matrix


def regularization_scale(kind: str, phi: np.ndarray, couple_over_length: np.ndarray) -> tuple[np.ndarray, float]:
    if kind == "unweighted_l2":
        return np.ones_like(phi), 0.0
    raw = np.abs(phi) if kind == "abs_phi_weighted" else np.abs(couple_over_length)
    if kind not in REGULARIZERS:
        raise ValueError(f"unknown regularizer: {kind}")
    maximum = float(np.max(raw))
    if not math.isfinite(maximum) or maximum <= 0.0:
        raise ValueError(f"{kind} has no positive scale")
    # A declared relative floor avoids converting an exact-zero baseline flux
    # into an impossible edge while retaining nine decades of weight contrast.
    floor = maximum * 1.0e-9
    return np.maximum(raw, floor) / maximum, floor


def weighted_minimum_norm(
    incidence: np.ndarray,
    required: np.ndarray,
    scale: np.ndarray,
) -> dict[str, Any]:
    if incidence.shape[0] != required.size or incidence.shape[1] != scale.size:
        raise ValueError("weighted inverse dimension mismatch")
    if np.any(~np.isfinite(incidence)) or np.any(~np.isfinite(required)) or np.any(~np.isfinite(scale)):
        raise ValueError("weighted inverse received non-finite data")
    if np.any(scale <= 0.0):
        raise ValueError("regularization scales must be positive")
    transformed = incidence * scale[np.newaxis, :]
    z, _, rank, singular = np.linalg.lstsq(transformed, required, rcond=1.0e-13)
    if rank != incidence.shape[0]:
        raise ValueError(f"incidence inverse is row-rank deficient: {rank}")
    delta = scale * z
    closure = incidence @ delta - required
    relative = norm(closure) / max(norm(required), 1.0e-300)
    if relative > 1.0e-8:
        raise ValueError(f"inverse conservation closure failed: {relative}")
    return {
        "delta": delta,
        "closure": closure,
        "relative_closure_l2": relative,
        "rank": int(rank),
        "condition_number": float(singular[0] / singular[-1]),
        "penalty_l2": norm(delta / scale),
    }


def triangle_is_obtuse(points: list[tuple[float, float]]) -> bool:
    if len(points) != 3:
        raise ValueError("triangle must have three points")
    squared = []
    for index in range(3):
        x0, y0 = points[index]
        x1, y1 = points[(index + 1) % 3]
        squared.append((x1 - x0) ** 2 + (y1 - y0) ** 2)
    squared.sort()
    return squared[2] > squared[0] + squared[1] + 1.0e-12 * squared[2]


def mesh_edge_attributes(mesh: dict[str, Any]) -> dict[tuple[int, int], dict[str, float]]:
    nodes = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    silicon_regions = {
        int(region["id"])
        for region in mesh["regions"]
        if str(region.get("material", "")).strip().lower() in {"si", "silicon"}
    }
    incident: dict[tuple[int, int], list[bool]] = defaultdict(list)
    silicon_centroids: dict[int, list[tuple[float, float]]] = defaultdict(list)
    for triangle in mesh["triangles"]:
        if int(triangle["region_id"]) not in silicon_regions:
            continue
        ids = [int(value) for value in triangle["node_ids"]]
        points = [nodes[value] for value in ids]
        obtuse = triangle_is_obtuse(points)
        centroid = (
            sum(point[0] for point in points) / 3.0,
            sum(point[1] for point in points) / 3.0,
        )
        for node in ids:
            silicon_centroids[node].append(centroid)
        for index in range(3):
            pair = tuple(sorted((ids[index], ids[(index + 1) % 3])))
            incident[pair].append(obtuse)

    _, pairs = region_local_geometry(mesh, {"si", "silicon"})
    direct = {int(row["global_node_id"]) for row in pairs}
    normals: dict[int, np.ndarray] = {}
    for node in direct:
        centers = silicon_centroids.get(node, [])
        if not centers:
            raise ValueError(f"interface node lacks silicon triangle support: {node}")
        point = np.asarray(nodes[node], dtype=float)
        vector = np.mean(np.asarray(centers, dtype=float), axis=0) - point
        magnitude = norm(vector)
        if magnitude <= 0.0:
            raise ValueError(f"interface normal is undefined at node {node}")
        normals[node] = vector / magnitude

    result: dict[tuple[int, int], dict[str, float]] = {}
    for pair, flags in incident.items():
        tail, head = pair
        edge = np.asarray(nodes[head], dtype=float) - np.asarray(nodes[tail], dtype=float)
        edge_norm = norm(edge)
        anchors = [node for node in pair if node in normals]
        if not anchors:
            midpoint = (
                (nodes[tail][0] + nodes[head][0]) / 2.0,
                (nodes[tail][1] + nodes[head][1]) / 2.0,
            )
            anchors = [min(direct, key=lambda node: math.dist(nodes[node], midpoint))]
        normal = np.mean(np.asarray([normals[node] for node in anchors]), axis=0)
        normal /= norm(normal)
        normal_fraction = abs(float(np.dot(edge / edge_norm, normal)))
        result[pair] = {
            "obtuse_incident_fraction": sum(flags) / len(flags),
            "interface_normal_fraction": normal_fraction,
            "interface_parallel_fraction": math.sqrt(max(0.0, 1.0 - normal_fraction ** 2)),
        }
    return result


def one_parameter_fit(required: np.ndarray, predicted_direction: np.ndarray) -> dict[str, float]:
    denominator = float(np.dot(predicted_direction, predicted_direction))
    if denominator <= 1.0e-300:
        return {"alpha": 0.0, "explained_residual_energy": 0.0}
    alpha = float(np.dot(predicted_direction, required) / denominator)
    remainder = required - alpha * predicted_direction
    total = float(np.dot(required, required))
    explained = 1.0 - float(np.dot(remainder, remainder)) / max(total, 1.0e-300)
    return {"alpha": alpha, "explained_residual_energy": max(0.0, min(1.0, explained))}


def permutation_p_value(
    incidence: np.ndarray,
    required: np.ndarray,
    signed_phi: np.ndarray,
    labels: np.ndarray,
    observed: float,
    *,
    count: int,
    seed: int,
) -> float:
    if count < 99:
        raise ValueError("permutation count must be at least 99")
    if np.allclose(labels, labels[0], rtol=0.0, atol=0.0):
        return 1.0
    generator = np.random.default_rng(seed)
    exceed = 0
    for _ in range(count):
        direction = incidence @ (signed_phi * generator.permutation(labels))
        score = one_parameter_fit(required, direction)["explained_residual_energy"]
        exceed += score >= observed - 1.0e-15
    return (exceed + 1.0) / (count + 1.0)


def vector_similarity(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    denominator = norm(left) * norm(right)
    cosine = float(np.dot(left, right) / denominator) if denominator > 0.0 else 1.0
    threshold = max(float(np.max(np.abs(left))), float(np.max(np.abs(right)))) * 1.0e-9
    active = (np.abs(left) > threshold) | (np.abs(right) > threshold)
    sign_agreement = float(np.mean(np.sign(left[active]) == np.sign(right[active]))) if np.any(active) else 1.0
    k = min(5, left.size)
    top_left = set(np.argsort(np.abs(left))[-k:])
    top_right = set(np.argsort(np.abs(right))[-k:])
    return {
        "cosine": cosine,
        "active_sign_agreement": sign_agreement,
        "top5_jaccard": len(top_left & top_right) / len(top_left | top_right),
    }


def state_graph(
    state: dict[str, Any],
    node_ids: tuple[int, ...],
    static_attributes: dict[tuple[int, int], dict[str, float]],
) -> dict[str, Any]:
    edges_all = read_csv(Path(state["sg_csv"]))
    carrier_all = read_csv(Path(state["carrier_csv"]))
    carrier = {int(row["node_id"]): row for row in carrier_all}
    if set(node_ids) - set(carrier):
        raise ValueError("T1 carrier table lacks graph nodes")
    for node in node_ids:
        row = carrier[node]
        if finite(row["electron_flux_abs_sum"], "electron_flux_abs_sum") <= 0.0:
            raise ValueError(f"node {node} is not a transport row")
        if abs(finite(row["electron_gauge"], "electron_gauge")) > 1.0e-300 or abs(finite(row["electron_boundary"], "electron_boundary")) > 1.0e-300:
            raise ValueError(f"node {node} is not a free electron continuity row")
    selected = []
    node_set = set(node_ids)
    for row in edges_all:
        n0, n1 = int(row["node0"]), int(row["node1"])
        if n0 not in node_set and n1 not in node_set:
            continue
        if finite(row["couple_m"], "couple_m") <= 0.0 or finite(row["electron_mobility_m2_V_s"], "electron_mobility") <= 0.0:
            continue
        if finite(row["length_m"], "length_m") <= 0.0:
            raise ValueError("active transport edge has non-positive length")
        selected.append(row)
    selected.sort(key=lambda row: tuple(sorted((int(row["node0"]), int(row["node1"])))))
    pairs = [tuple(sorted((int(row["node0"]), int(row["node1"])))) for row in selected]
    if len(pairs) != 24:
        raise ValueError(f"R4 frozen edge domain must contain 24 active edges, got {len(pairs)}")
    incidence = build_incidence(node_ids, pairs)
    if np.linalg.matrix_rank(incidence) != len(node_ids):
        raise ValueError("R4 incidence matrix is not full row rank")
    scale_info = current_scale_from_sg(edges_all)
    current_scale = scale_info["particle_per_m_s_per_scaled"] * ELEMENTARY_CHARGE_C / 1.0e6
    residual = np.asarray([finite(carrier[node]["electron_residual"], "electron_residual") * current_scale for node in node_ids])
    required = -residual
    signed_phi = []
    couple_length = []
    attributes: dict[str, list[float]] = {name: [] for name in FAMILY_NAMES}
    edge_records = []
    for row, pair in zip(selected, pairs, strict=True):
        n0 = int(row["node0"])
        orientation = 1.0 if n0 == pair[0] else -1.0
        flux = orientation * finite(row["electron_flux"], "electron_flux") * current_scale
        signed_phi.append(flux)
        ratio = finite(row["couple_m"], "couple_m") / finite(row["length_m"], "length_m")
        couple_length.append(ratio)
        doping0 = finite(carrier[pair[0]]["net_doping_m3"], "net_doping_m3") * 1.0e6
        doping1 = finite(carrier[pair[1]]["net_doping_m3"], "net_doping_m3") * 1.0e6
        exported_doping_average = finite(row["net_doping_avg_m3"], "net_doping_avg_m3")
        reconstructed_doping_average = 0.5 * (doping0 + doping1)
        doping_average_error = abs(exported_doping_average - reconstructed_doping_average) / max(
            abs(exported_doping_average), abs(reconstructed_doping_average), 1.0
        )
        if doping_average_error > 1.0e-12:
            raise ValueError(f"carrier/SG doping unit contract mismatch on {pair}")
        mesh_attrs = static_attributes.get(pair)
        if mesh_attrs is None:
            raise ValueError(f"active edge absent from silicon mesh: {pair}")
        values = {
            "bernoulli_abs_eta": abs(finite(row["electron_bernoulli_argument"], "eta")),
            "interface_normal_fraction": mesh_attrs["interface_normal_fraction"],
            "interface_parallel_fraction": mesh_attrs["interface_parallel_fraction"],
            "obtuse_incident_fraction": mesh_attrs["obtuse_incident_fraction"],
            "doping_gradient_abs": abs(doping1 - doping0) / finite(row["length_m"], "length_m"),
            "hfs_drive_abs": abs(finite(row["electron_mobility_field_V_m"], "HFS drive")),
        }
        for name, value in values.items():
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"invalid edge attribute {name} on {pair}")
            attributes[name].append(value)
        edge_records.append({"edge_id": int(row["edge_id"]), "tail": pair[0], "head": pair[1], "cross_domain": int((pair[0] in node_set) ^ (pair[1] in node_set)), **values})
    signed_phi_array = np.asarray(signed_phi)
    reconstructed = incidence @ signed_phi_array
    transport_residual = np.asarray([finite(carrier[node]["electron_flux"], "electron_flux") * current_scale for node in node_ids])
    reconstruction_error = norm(reconstructed - transport_residual) / max(norm(transport_residual), 1.0e-300)
    if reconstruction_error > 1.0e-10:
        raise ValueError(f"canonical edge orientation fails T1 transport reconstruction: {reconstruction_error}")
    return {
        "bias_V": float(state["bias_V"]), "pairs": pairs, "incidence": incidence,
        "required": required, "residual": residual, "signed_phi": signed_phi_array,
        "couple_over_length": np.asarray(couple_length),
        "attributes": {name: np.asarray(values) for name, values in attributes.items()},
        "edge_records": edge_records, "carrier": carrier,
        "current_scale": current_scale, "reconstruction_error": reconstruction_error,
        "input_sha256": {"sg_csv": sha256(Path(state["sg_csv"])), "carrier_csv": sha256(Path(state["carrier_csv"]))},
    }


def analyze(root: Path, permutation_count: int) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    t1_summary, manifest, mesh_path = load_t1_contract(root)
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    static_attributes = mesh_edge_attributes(mesh)
    node_ids = tuple(int(value) for value in t1_summary["contract"]["selected_nodes"])
    graphs = [state_graph(state, node_ids, static_attributes) for state in manifest["states"]]
    expected_pairs = graphs[0]["pairs"]
    if any(graph["pairs"] != expected_pairs for graph in graphs[1:]):
        raise ValueError("T1 active edge topology changes across states")

    solution_rows: list[dict[str, Any]] = []
    node_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    state_summaries = []
    endpoint_stability = []
    for state_index, graph in enumerate(graphs):
        solutions: dict[str, np.ndarray] = {}
        regularization_reports = {}
        for kind in REGULARIZERS:
            scale, floor = regularization_scale(kind, graph["signed_phi"], graph["couple_over_length"])
            solved = weighted_minimum_norm(graph["incidence"], graph["required"], scale)
            solutions[kind] = solved["delta"]
            regularization_reports[kind] = {key: value for key, value in solved.items() if key not in {"delta", "closure"}}
            regularization_reports[kind]["raw_scale_floor"] = floor
            for edge, delta, local_scale in zip(graph["edge_records"], solved["delta"], scale, strict=True):
                solution_rows.append({
                    "bias_V": graph["bias_V"], "regularizer": kind,
                    **edge, "phi_A_per_um": graph["signed_phi"][expected_pairs.index((edge['tail'], edge['head']))],
                    "delta_required_A_per_um": delta, "regularization_scale": local_scale,
                })
        similarities = {}
        for left_index, left in enumerate(REGULARIZERS):
            for right in REGULARIZERS[left_index + 1:]:
                similarities[f"{left}__vs__{right}"] = vector_similarity(solutions[left], solutions[right])
        stable = all(
            report["cosine"] >= 0.9 and report["active_sign_agreement"] >= 0.85
            for report in similarities.values()
        )
        if any(math.isclose(graph["bias_V"], endpoint, abs_tol=5e-7) for endpoint in ENDPOINTS):
            endpoint_stability.append(stable)

        nodal_scales = []
        for row_index, node in enumerate(node_ids):
            incident = np.abs(graph["incidence"][row_index]) > 0.0
            denominator = float(np.sum(np.abs(graph["signed_phi"][incident])))
            nodal_scale = graph["required"][row_index] / max(denominator, 1.0e-300)
            nodal_scales.append(nodal_scale)
            node_rows.append({"bias_V": graph["bias_V"], "node_id": node, "required_A_per_um": graph["required"][row_index], "incident_abs_phi_A_per_um": denominator, "required_over_incident_abs_phi": nodal_scale})

        family_reports = {}
        for family_index, name in enumerate(FAMILY_NAMES):
            labels = graph["attributes"][name]
            label_scale = float(np.max(labels))
            normalized = labels / label_scale if label_scale > 0.0 else labels.copy()
            direction = graph["incidence"] @ (graph["signed_phi"] * normalized)
            fit = one_parameter_fit(graph["required"], direction)
            p_value = permutation_p_value(
                graph["incidence"], graph["required"], graph["signed_phi"], normalized,
                fit["explained_residual_energy"], count=permutation_count,
                seed=0x5EED + 1009 * state_index + 97 * family_index,
            )
            report = {**fit, "permutation_p_value": p_value, "label_max": label_scale}
            association = {}
            edge_basis = graph["signed_phi"] * normalized
            for kind in REGULARIZERS:
                solution = solutions[kind]
                fit_solution = one_parameter_fit(solution, edge_basis)
                association[kind] = fit_solution["explained_residual_energy"]
            report["inverse_solution_edge_energy_by_regularizer"] = association
            report["minimum_inverse_solution_edge_energy"] = min(association.values())
            family_reports[name] = report
            fit_rows.append({
                "bias_V": graph["bias_V"], "family": name,
                "alpha": report["alpha"],
                "explained_residual_energy": report["explained_residual_energy"],
                "permutation_p_value": report["permutation_p_value"],
                "label_max": report["label_max"],
                "minimum_inverse_solution_edge_energy": report["minimum_inverse_solution_edge_energy"],
                **{f"inverse_edge_energy_{kind}": association[kind] for kind in REGULARIZERS},
            })
        best = max(family_reports, key=lambda name: family_reports[name]["explained_residual_energy"])
        state_summaries.append({
            "bias_V": graph["bias_V"], "node_count": len(node_ids), "edge_count": len(expected_pairs),
            "internal_shared_edge_count": sum(pair[0] in node_ids and pair[1] in node_ids for pair in expected_pairs),
            "cross_domain_edge_count": sum((pair[0] in node_ids) ^ (pair[1] in node_ids) for pair in expected_pairs),
            "incidence_rank": int(np.linalg.matrix_rank(graph["incidence"])),
            "transport_reconstruction_relative_l2": graph["reconstruction_error"],
            "regularizers": regularization_reports, "regularization_similarity": similarities,
            "regularization_structure_stable": stable, "best_family": best,
            "family_fits": family_reports, "input_sha256": graph["input_sha256"],
            "nodal_required_scale_abs_min": min(abs(value) for value in nodal_scales),
            "nodal_required_scale_abs_max": max(abs(value) for value in nodal_scales),
        })

    endpoint_reports = [state for state in state_summaries if any(math.isclose(state["bias_V"], endpoint, abs_tol=5e-7) for endpoint in ENDPOINTS)]
    nominations = []
    for name in FAMILY_NAMES:
        if all(
            state["family_fits"][name]["explained_residual_energy"] >= 0.9
            and state["family_fits"][name]["permutation_p_value"] <= 0.05
            for state in endpoint_reports
        ) and all(endpoint_stability):
            nominations.append(name)

    endpoint_node_rows = [row for row in node_rows if any(math.isclose(row["bias_V"], endpoint, abs_tol=5e-7) for endpoint in ENDPOINTS)]
    nodal_values = [float(row["required_over_incident_abs_phi"]) for row in endpoint_node_rows]
    dominant_sign = 1.0 if np.median(nodal_values) >= 0.0 else -1.0
    sign_flip_fraction = sum(value * dominant_sign < 0.0 for value in nodal_values) / len(nodal_values)
    nonzero_abs = [abs(value) for value in nodal_values if abs(value) > 1.0e-300]
    order_span = math.log10(max(nonzero_abs) / min(nonzero_abs)) if nonzero_abs else math.inf
    max_endpoint_explanation = max(
        min(state["family_fits"][name]["explained_residual_energy"] for state in endpoint_reports)
        for name in FAMILY_NAMES
    )
    max_inverse_attribute_energy = max(
        state["family_fits"][name]["minimum_inverse_solution_edge_energy"]
        for state in endpoint_reports for name in FAMILY_NAMES
    )
    no_attribute_50_all_regularizers = (
        max_endpoint_explanation < 0.5 and max_inverse_attribute_energy < 0.5
    )
    stop_evidence = sign_flip_fraction > 0.15 and order_span >= 2.0 and no_attribute_50_all_regularizers
    summary = {
        "schema": OUTPUT_SCHEMA,
        "diagnostic_only": True,
        "interpretation": "required feasible Vela edge correction, not measured Sentaurus/Vela edge-flux difference",
        "contract": {
            "frozen_plan_commit": "01f20ace88dc2b68a7bc3788534244dadf0d18f2",
            "t1_root": str(root.resolve()), "t1_summary_sha256": sha256(root / "summary.json"),
            "t1_manifest_sha256": sha256(root / "state_manifest.json"), "mesh_sha256": sha256(mesh_path),
            "node_domain": list(node_ids), "orientation": "low_node_id_to_high_node_id; B(tail)=+1; B(head)=-1",
            "cross_domain_rule": "only the in-domain endpoint receives an incidence entry",
            "regularizers": list(REGULARIZERS), "permutation_count": permutation_count,
            "stop_sign_metric": "fraction of 14 endpoint node-state required/incident-abs-flux values opposite the pooled median sign",
            "stop_order_metric": "log10(max/min nonzero absolute endpoint nodal required/incident-abs-flux)",
            "production_defaults_changed": False, "vm_used": False, "reclose_run": False, "idvg_31_point_run": False,
        },
        "states": state_summaries,
        "gates": {
            "candidate_explanation_threshold": 0.9, "permutation_p_threshold": 0.05,
            "regularization_cosine_threshold": 0.9, "regularization_sign_agreement_threshold": 0.85,
            "nominated_candidate_families": nominations,
            "candidate_nomination_passed": bool(nominations),
            "candidate_nomination_is_root_cause_confirmation": False,
            "stop_sign_flip_fraction": sign_flip_fraction, "stop_sign_flip_threshold_strict": 0.15,
            "stop_magnitude_order_span": order_span, "stop_order_threshold": 2.0,
            "stop_max_single_attribute_endpoint_energy": max_endpoint_explanation,
            "stop_max_inverse_attribute_energy_across_all_regularizers": max_inverse_attribute_energy,
            "stop_no_attribute_explains_50_percent_all_regularizers": no_attribute_50_all_regularizers,
            "t2_stop_evidence_passed": stop_evidence,
        },
    }
    return summary, solution_rows, node_rows, fit_rows


def markdown_report(summary: dict[str, Any]) -> str:
    gates = summary["gates"]
    endpoints = [state for state in summary["states"] if any(math.isclose(state["bias_V"], endpoint, abs_tol=5e-7) for endpoint in ENDPOINTS)]
    lines = [
        "# Templates/LDMOS G3 WP3 T2 edge-correction inverse",
        "", "## Decision", "",
        "This is a feasibility/screening inverse only. The inferred edge correction is not a measured Sentaurus/Vela edge-flux difference and cannot confirm root cause.",
        "", f"Candidate families nominated: `{gates['nominated_candidate_families']}`. T2 stop-evidence gate: `{gates['t2_stop_evidence_passed']}`.",
        "", "## R4 graph and conservation", "",
        "The node domain is the frozen seven free-silicon electron-continuity rows. Every active silicon edge is oriented low-ID to high-ID; internal shared edges enter both endpoints with opposite signs and cross-domain edges enter only the in-domain row.",
        "", "| Vg | nodes | edges | shared | cross-domain | rank | max inverse closure | stable |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for state in summary["states"]:
        maximum = max(report["relative_closure_l2"] for report in state["regularizers"].values())
        lines.append(f"| {state['bias_V']:.6g} | {state['node_count']} | {state['edge_count']} | {state['internal_shared_edge_count']} | {state['cross_domain_edge_count']} | {state['incidence_rank']} | {maximum:.3e} | {state['regularization_structure_stable']} |")
    lines += ["", "## Endpoint single-parameter fits", "", "| family | Vg=1 residual energy | p | min inverse-edge energy | Vg=7/6 residual energy | p | min inverse-edge energy |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"]
    for name in FAMILY_NAMES:
        left = endpoints[0]["family_fits"][name]
        right = endpoints[1]["family_fits"][name]
        lines.append(f"| {name} | {left['explained_residual_energy']:.6g} | {left['permutation_p_value']:.4g} | {left['minimum_inverse_solution_edge_energy']:.6g} | {right['explained_residual_energy']:.6g} | {right['permutation_p_value']:.4g} | {right['minimum_inverse_solution_edge_energy']:.6g} |")
    lines += [
        "", "A family is only nominated when it explains at least 90% at both endpoints, passes p<=0.05 at both endpoints, and the inverse structure is stable across all three regularizers. A nomination would remain screening evidence, not root-cause confirmation.",
        "", "## Stop-evidence metrics", "",
        f"- Endpoint nodal sign-flip fraction: `{gates['stop_sign_flip_fraction']:.6g}` (strictly >0.15 required).",
        f"- Endpoint nodal magnitude span: `{gates['stop_magnitude_order_span']:.6g}` decades (>=2 required).",
        f"- Best single-attribute worst-endpoint energy: `{gates['stop_max_single_attribute_endpoint_energy']:.6g}` (<0.5 required, with regularization stability).",
        f"- Best attribute minimum inverse-edge energy across all regularizers: `{gates['stop_max_inverse_attribute_energy_across_all_regularizers']:.6g}` (<0.5 required).",
        f"- Combined T2 stop evidence: `{gates['t2_stop_evidence_passed']}`.",
        "", "## Files", "",
        "- `summary.json`: authoritative gates, hashes, state graphs, fits, and stability.",
        "- `edge_corrections.csv`: every state/regularizer canonical edge correction.",
        "- `nodal_required_scales.csv`: endpoint-consistency diagnostics.",
        "- `single_parameter_fits.csv`: all family fits and permutation p-values.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--t1-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=999)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    summary, corrections, nodes, fits = analyze(args.t1_root.resolve(), args.permutations)
    write_csv(output / "edge_corrections.csv", corrections)
    write_csv(output / "nodal_required_scales.csv", nodes)
    write_csv(output / "single_parameter_fits.csv", fits)
    summary["artifacts"] = {
        "edge_corrections": str((output / "edge_corrections.csv").resolve()),
        "nodal_required_scales": str((output / "nodal_required_scales.csv").resolve()),
        "single_parameter_fits": str((output / "single_parameter_fits.csv").resolve()),
        "report": str((output / "report.md").resolve()),
    }
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "report.md").write_text(markdown_report(summary), encoding="utf-8")
    print(json.dumps({"summary": str((output / 'summary.json').resolve()), "gates": summary["gates"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
