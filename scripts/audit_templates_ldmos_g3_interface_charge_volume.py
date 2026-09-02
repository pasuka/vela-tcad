#!/usr/bin/env python3
"""Audit the LDMOS Si/SiO2 interface semiconductor-charge volume contract.

This is a read-only P0 audit.  It reproduces Vela's legacy Poisson operator,
changes only the semiconductor charge volume from the global node volume to a
transport-cell barycentric volume, and reports three pre-registered channel
charge windows.  It never edits production inputs or runs a nonlinear solve.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


Q = 1.602176634e-19
EPS0 = 8.8541878128e-12
UM = 1.0e-6
EPS_R = {"Si": 11.7, "SiO2": 3.9}

DEFAULT_STATE_DIRS = {
    "0.166667": "reference_staging/templates_ldmos_g3_state_feedback_averagebox_node_local_20260831/vg_0p166667",
    "0.333333": "reference_staging/templates_ldmos_g3_state_feedback_averagebox_node_local_20260831/vg_0p333333",
    "0.500000": "reference_staging/templates_ldmos_g3_state_feedback_averagebox_node_local_20260831/vg_0p500000",
    "0.666667": "reference_staging/templates_ldmos_g3_state_feedback_averagebox_node_local_20260831/vg_0p666667",
    "1.000000": "reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/vg_1p000000",
    "1.166667": "reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/vg_1p166667",
}

ENDPOINT_SUMMARIES = {
    "1.000000": "reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/summary.json",
    "1.166667": "reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/summary.json",
}

WINDOWS = {
    "source_half": {"y_min_um": 2.06, "y_max_um": 2.50, "x_max_um": -9.85},
    "drain_half": {"y_min_um": 2.50, "y_max_um": 3.50, "x_max_um": -9.85},
    "whole_channel": {"y_min_um": 2.06, "y_max_um": 3.50, "x_max_um": -9.85},
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ids_hash(ids: np.ndarray) -> str:
    text = ",".join(str(int(value)) for value in ids)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def load_state(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"empty state file: {path}")
    fields = ("psi", "phin", "phip", "electrons_m3", "holes_m3")
    result = {field: np.zeros(len(rows), dtype=float) for field in fields}
    node_ids = np.zeros(len(rows), dtype=int)
    for row_index, row in enumerate(rows):
        node_ids[row_index] = int(row["node_id"])
        for field in fields:
            result[field][row_index] = float(row[field])
    order = np.argsort(node_ids)
    sorted_ids = node_ids[order]
    if not np.array_equal(sorted_ids, np.arange(len(rows), dtype=int)):
        raise ValueError(f"state node IDs are not the dense range [0,N): {path}")
    return {field: values[order] for field, values in result.items()}


def load_current_ratio(path: Path, bias: float) -> dict[str, float]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    point = min(payload["points"], key=lambda row: abs(float(row["bias_V"]) - bias))
    sentaurus = float(point["sentaurus_terminal_A_per_um"])
    vela = float(point["vela_curve_A_per_um"])
    return {
        "sentaurus_A_per_um": sentaurus,
        "vela_A_per_um": vela,
        "vela_to_sentaurus": vela / sentaurus,
    }


def build_geometry(mesh: dict[str, Any]) -> dict[str, Any]:
    nodes = mesh["nodes"]
    n_nodes = len(nodes)
    x_um = np.array([float(node["x"]) for node in nodes])
    y_um = np.array([float(node["y"]) for node in nodes])
    x = x_um * UM
    y = y_um * UM

    material_by_cell: dict[int, str] = {}
    for region in mesh["regions"]:
        for cell_id in region["cell_ids"]:
            material_by_cell[int(cell_id)] = str(region["material"])

    global_volume = np.zeros(n_nodes)
    transport_volume = np.zeros(n_nodes)
    silicon_nodes: set[int] = set()
    oxide_nodes: set[int] = set()
    edge_couple: dict[tuple[int, int], float] = {}
    edge_eps: dict[tuple[int, int], list[float]] = {}
    negative_cotangent_fallbacks = 0

    for triangle in mesh["triangles"]:
        cell_id = int(triangle["id"])
        ids = [int(value) for value in triangle["node_ids"]]
        if len(ids) != 3:
            raise ValueError(f"P0 audit requires Tri3 cells; cell {cell_id} has {len(ids)} nodes")
        material = material_by_cell[cell_id]
        if material not in EPS_R:
            raise ValueError(f"unsupported material in frozen LDMOS audit: {material}")
        points = np.array([[x[node], y[node]] for node in ids])
        area = 0.5 * abs(
            (points[1, 0] - points[0, 0]) * (points[2, 1] - points[0, 1])
            - (points[2, 0] - points[0, 0]) * (points[1, 1] - points[0, 1])
        )
        for node in ids:
            global_volume[node] += area / 3.0
            if material == "Si":
                transport_volume[node] += area / 3.0
                silicon_nodes.add(node)
            else:
                oxide_nodes.add(node)

        for local_index in range(3):
            node0 = ids[local_index]
            node1 = ids[(local_index + 1) % 3]
            opposite = ids[(local_index + 2) % 3]
            key = (min(node0, node1), max(node0, node1))
            length = math.hypot(x[node1] - x[node0], y[node1] - y[node0])
            u = np.array([x[node0] - x[opposite], y[node0] - y[opposite]])
            v = np.array([x[node1] - x[opposite], y[node1] - y[opposite]])
            cross = float(u[0] * v[1] - u[1] * v[0])
            cotangent = float(np.dot(u, v) / abs(cross)) if abs(cross) > 1.0e-60 else 0.0
            local_couple = 0.5 * cotangent * length
            if cotangent < 0.0:
                local_couple = area / (3.0 * length)
                negative_cotangent_fallbacks += 1
            edge_couple[key] = edge_couple.get(key, 0.0) + local_couple
            edge_eps.setdefault(key, []).append(EPS_R[material])

    edge_keys = list(edge_couple)
    edge_i = np.array([key[0] for key in edge_keys], dtype=int)
    edge_j = np.array([key[1] for key in edge_keys], dtype=int)
    edge_length = np.hypot(x[edge_i] - x[edge_j], y[edge_i] - y[edge_j])
    conductance = np.array(
        [np.mean(edge_eps[key]) * EPS0 * edge_couple[key] for key in edge_keys]
    ) / edge_length

    contact_nodes: set[int] = set()
    for contact in mesh["contacts"]:
        contact_nodes.update(int(value) for value in contact["node_ids"])

    silicon = np.array(sorted(silicon_nodes), dtype=int)
    interface = np.array(sorted((silicon_nodes & oxide_nodes) - contact_nodes), dtype=int)
    interior = np.array(sorted(silicon_nodes - oxide_nodes - contact_nodes), dtype=int)
    channel_interface = interface[
        (y_um[interface] >= 2.04) & (y_um[interface] <= 3.53)
    ]

    return {
        "x_um": x_um,
        "y_um": y_um,
        "global_volume": global_volume,
        "transport_volume": transport_volume,
        "silicon_nodes": silicon,
        "interface_nodes": interface,
        "interior_nodes": interior,
        "channel_interface_nodes": channel_interface,
        "edge_i": edge_i,
        "edge_j": edge_j,
        "conductance": conductance,
        "negative_cotangent_fallbacks": negative_cotangent_fallbacks,
    }


def poisson_flux(psi: np.ndarray, geometry: dict[str, Any]) -> np.ndarray:
    edge_i = geometry["edge_i"]
    edge_j = geometry["edge_j"]
    contribution = geometry["conductance"] * (psi[edge_j] - psi[edge_i])
    residual = np.zeros(len(psi))
    np.add.at(residual, edge_i, contribution)
    np.add.at(residual, edge_j, -contribution)
    return residual


def distribution(values: np.ndarray, indices: np.ndarray) -> dict[str, float]:
    selected = np.abs(values[indices])
    return {
        "median": float(np.median(selected)),
        "p95": float(np.percentile(selected, 95.0)),
        "maximum": float(np.max(selected)),
    }


def load_doping(path: Path, n_nodes: int) -> np.ndarray:
    doping = np.zeros(n_nodes)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            doping[int(row["node_id"])] = (
                float(row["donors_cm3"]) - float(row["acceptors_cm3"])
            ) * 1.0e6
    return doping


def load_measure(path: Path, n_nodes: int) -> np.ndarray:
    measure = np.zeros(n_nodes)
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            measure[int(row["node_id"])] = (
                float(row["averagebox_silicon_measure_um2"]) * UM * UM
            )
    return measure


def window_nodes(geometry: dict[str, Any], window: dict[str, float]) -> np.ndarray:
    silicon = geometry["silicon_nodes"]
    x_um = geometry["x_um"]
    y_um = geometry["y_um"]
    return silicon[
        (y_um[silicon] >= window["y_min_um"])
        & (y_um[silicon] <= window["y_max_um"])
        & (x_um[silicon] < window["x_max_um"])
    ]


def endpoint_audit(
    state_dir: Path,
    geometry: dict[str, Any],
    doping: np.ndarray,
    sentaurus_measure: np.ndarray,
    current: dict[str, float] | None,
) -> dict[str, Any]:
    sentaurus_path = state_dir / "sentaurus_state.csv"
    vela_path = state_dir / "VVV/hybrid_state.csv"
    sentaurus = load_state(sentaurus_path)
    vela = load_state(vela_path)
    global_volume = geometry["global_volume"]
    transport_volume = geometry["transport_volume"]

    rho_sentaurus = Q * (
        sentaurus["holes_m3"] - sentaurus["electrons_m3"] + doping
    )
    rho_vela = Q * (vela["holes_m3"] - vela["electrons_m3"] + doping)
    flux_sentaurus = poisson_flux(sentaurus["psi"], geometry)
    flux_vela = poisson_flux(vela["psi"], geometry)
    residual_vela = flux_vela + rho_vela * global_volume
    residual_global = flux_sentaurus + rho_sentaurus * global_volume
    residual_local = flux_sentaurus + rho_sentaurus * transport_volume

    scale_vela = (
        Q
        * (vela["electrons_m3"] + vela["holes_m3"] + np.abs(doping))
        * global_volume
        + 1.0e-300
    )
    scale_sentaurus = (
        Q
        * (sentaurus["electrons_m3"] + sentaurus["holes_m3"] + np.abs(doping))
        * sentaurus_measure
        + 1.0e-300
    )
    normalized_vela = residual_vela / scale_vela
    normalized_global = residual_global / scale_sentaurus
    normalized_local = residual_local / scale_sentaurus

    channel = geometry["channel_interface_nodes"]
    global_stats = distribution(normalized_global, channel)
    local_stats = distribution(normalized_local, channel)
    reduction = {
        key: 1.0 - local_stats[key] / max(global_stats[key], 1.0e-300)
        for key in ("median", "p95", "maximum")
    }

    charge_windows: dict[str, Any] = {}
    for name, definition in WINDOWS.items():
        selected = window_nodes(geometry, definition)
        sentaurus_charge = float(np.sum(sentaurus["electrons_m3"][selected] * sentaurus_measure[selected]))
        vela_global_charge = float(np.sum(vela["electrons_m3"][selected] * global_volume[selected]))
        vela_local_charge = float(np.sum(vela["electrons_m3"][selected] * transport_volume[selected]))
        row: dict[str, Any] = {
            "node_count": int(len(selected)),
            "node_ids_sha256": ids_hash(selected),
            "definition": definition,
            "vela_global_to_sentaurus_measure": vela_global_charge / sentaurus_charge,
            "vela_transport_to_sentaurus_measure": vela_local_charge / sentaurus_charge,
        }
        if current is not None:
            row["current_ratio"] = current["vela_to_sentaurus"]
            row["transport_charge_minus_current_abs"] = abs(
                row["vela_transport_to_sentaurus_measure"] - current["vela_to_sentaurus"]
            )
        charge_windows[name] = row

    return {
        "inputs": {
            "sentaurus_state": str(sentaurus_path),
            "sentaurus_state_sha256": sha256_file(sentaurus_path),
            "vela_state": str(vela_path),
            "vela_state_sha256": sha256_file(vela_path),
        },
        "operator_replica_on_vela_state": {
            "interior": distribution(normalized_vela, geometry["interior_nodes"]),
            "interface": distribution(normalized_vela, geometry["interface_nodes"]),
        },
        "sentaurus_state_volume_ab": {
            "global": global_stats,
            "material_local_barycentric": local_stats,
            "fractional_reduction": reduction,
        },
        "current": current,
        "charge_windows": charge_windows,
    }


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=repo)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo = args.repo.resolve()
    mesh_path = repo / "reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/stage1_v4/vela_exact_topology/mesh.json"
    doping_path = repo / "reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/stage1_v4/vela_exact_topology/doping.csv"
    measure_path = repo / "reference_staging/templates_ldmos_averagebox_full_mesh_20260831/node_residuals.csv"

    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    geometry = build_geometry(mesh)
    doping = load_doping(doping_path, len(mesh["nodes"]))
    measure = load_measure(measure_path, len(mesh["nodes"]))

    points: dict[str, Any] = {}
    for bias_text, relative in DEFAULT_STATE_DIRS.items():
        current = None
        if bias_text in ENDPOINT_SUMMARIES:
            current = load_current_ratio(repo / ENDPOINT_SUMMARIES[bias_text], float(bias_text))
        points[bias_text] = endpoint_audit(
            repo / relative, geometry, doping, measure, current
        )

    endpoint_rows = [points[key] for key in ("1.000000", "1.166667")]
    replica_pass = all(
        row["operator_replica_on_vela_state"][scope]["median"] <= 1.0e-10
        for row in points.values()
        for scope in ("interior", "interface")
    )
    volume_pass = all(
        row["sentaurus_state_volume_ab"]["fractional_reduction"][metric] >= 0.90
        for row in endpoint_rows
        for metric in ("median", "p95")
    )
    window_sensitivity = {
        name: {
            "accounting_ratio_min": min(
                row["charge_windows"][name]["vela_global_to_sentaurus_measure"]
                for row in endpoint_rows
            ),
            "accounting_ratio_max": max(
                row["charge_windows"][name]["vela_global_to_sentaurus_measure"]
                for row in endpoint_rows
            ),
            "max_transport_charge_minus_current_abs": max(
                row["charge_windows"][name]["transport_charge_minus_current_abs"]
                for row in endpoint_rows
            ),
        }
        for name in WINDOWS
    }
    whole = window_sensitivity["whole_channel"]
    whole_channel_pass = (
        whole["accounting_ratio_min"] >= 0.98
        and whole["accounting_ratio_max"] <= 1.05
        and whole["max_transport_charge_minus_current_abs"] <= 0.05
    )

    output = {
        "schema": "vela.templates_ldmos_g3_interface_charge_volume_audit.v1",
        "contract": {
            "read_only": True,
            "poisson_couples": "vela_legacy_edge_average_epsilon_cotangent_with_negative_fallback",
            "baseline_charge_volume": "global_barycentric",
            "candidate_charge_volume": "transport_cell_barycentric",
            "transport_material": "Si",
            "current_contract": {
                "predictor": "disabled",
                "ialmob": "disabled",
                "carrier_transport_couple_profile": "templates_ldmos_external_averagebox",
            },
        },
        "inputs": {
            "mesh": str(mesh_path),
            "mesh_sha256": sha256_file(mesh_path),
            "doping": str(doping_path),
            "doping_sha256": sha256_file(doping_path),
            "sentaurus_measure": str(measure_path),
            "sentaurus_measure_sha256": sha256_file(measure_path),
        },
        "masks": {
            "interface": {
                "count": int(len(geometry["interface_nodes"])),
                "node_ids_sha256": ids_hash(geometry["interface_nodes"]),
            },
            "channel_interface": {
                "count": int(len(geometry["channel_interface_nodes"])),
                "node_ids_sha256": ids_hash(geometry["channel_interface_nodes"]),
                "y_min_um": 2.04,
                "y_max_um": 3.53,
            },
            "charge_windows": WINDOWS,
        },
        "geometry": {
            "edge_count": int(len(geometry["edge_i"])),
            "negative_cotangent_fallback_count": int(geometry["negative_cotangent_fallbacks"]),
            "free_interface_transport_to_global_volume_median": float(
                np.median(
                    geometry["transport_volume"][geometry["interface_nodes"]]
                    / geometry["global_volume"][geometry["interface_nodes"]]
                )
            ),
        },
        "coverage": {
            "matched_state_count": len(points),
            "requested_state_count": 8,
            "missing_biases_V": [0.0, 0.833333],
            "coverage_complete": len(points) == 8,
        },
        "points": points,
        "window_sensitivity": window_sensitivity,
        "gates": {
            "operator_replica_median_le_1e-10": replica_pass,
            "endpoint_volume_ab_median_p95_reduction_ge_0p90": volume_pass,
            "whole_channel_accounting_and_current_tracking": whole_channel_pass,
            "coverage_eight_matched_states": len(points) == 8,
        },
        "verdict": (
            "P0_CORE_PASS_COVERAGE_INCOMPLETE"
            if replica_pass and volume_pass and whole_channel_pass
            else "P0_CORE_FAIL"
        ),
    }

    from templates_ldmos_contracts import validate_document
    validate_document(output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output), "verdict": output["verdict"], "gates": output["gates"]}, indent=2))
    return 0 if output["verdict"].startswith("P0_CORE_PASS") else 2


if __name__ == "__main__":
    raise SystemExit(main())
