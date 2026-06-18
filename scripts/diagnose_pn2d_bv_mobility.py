#!/usr/bin/env python3
"""Diagnose PN2D BV mobility decomposition against Sentaurus multibias exports."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="0,-0.1,-0.2")
    parser.add_argument("--electron-saturation-velocity-m-s", type=float, default=1.0e5)
    parser.add_argument("--hole-saturation-velocity-m-s", type=float, default=1.0e5)
    parser.add_argument("--electron-field-beta", type=float, default=2.0)
    parser.add_argument("--hole-field-beta", type=float, default=2.0)
    parser.add_argument("--top-error-neighborhood-radius-um", type=float, default=0.05)
    parser.add_argument("--top-error-neighborhood-max-nodes", type=int, default=25)
    parser.add_argument("--anchor-profile-tolerance-um", type=float, default=1.0e-5)
    parser.add_argument("--anchor-profile-max-nodes", type=int, default=256)
    return parser.parse_args()


def parse_biases(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def bias_key(value: float) -> float:
    return round(value, 12)


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def bias_token(bias: float) -> str:
    return f"{bias:g}".replace("-", "m").replace("+", "").replace(".", "p")


def discover_sentaurus_exports(root: Path, biases: list[float]) -> dict[float, Path]:
    result: dict[float, Path] = {}
    for bias in biases:
        candidates = [
            root / f"sentaurus_{signed_bias_token(bias)}v",
            root / f"sentaurus_{bias_token(bias)}v",
            root / f"sentaurus_{bias_token(abs(bias))}v",
        ]
        for candidate in candidates:
            if (candidate / "nodes.csv").exists():
                result[bias_key(bias)] = candidate
                break
    return result


def discover_vela_vtks(root: Path) -> dict[float, Path]:
    result: dict[float, Path] = {}
    mtimes: dict[float, float] = {}
    pattern = re.compile(r"_(?P<bias>[-+0-9.]+)V\.vtk$")
    for path in sorted(root.glob("*.vtk")):
        match = pattern.search(path.name)
        if not match:
            continue
        key = bias_key(float(match.group("bias")))
        mtime = path.stat().st_mtime
        if key not in result or mtime >= mtimes[key]:
            result[key] = path
            mtimes[key] = mtime
    return result


def load_sentaurus_mesh(root: Path) -> tuple[np.ndarray, np.ndarray, list[int]]:
    nodes: dict[int, tuple[float, float]] = {}
    with (root / "nodes.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            nodes[int(row["id"])] = (float(row["x_um"]), float(row["y_um"]))
    node_ids = sorted(nodes)
    xy = np.array([nodes[node_id] for node_id in node_ids], dtype=float)
    return xy[:, 0], xy[:, 1], node_ids


def load_sentaurus_triangles(root: Path, node_ids: list[int]) -> np.ndarray:
    index_by_node = {node_id: idx for idx, node_id in enumerate(node_ids)}
    triangles: list[list[int]] = []
    path = root / "elements.csv"
    if not path.exists():
        return np.empty((0, 3), dtype=int)
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            nodes = [int(row[name]) for name in ("node0", "node1", "node2") if row.get(name, "") != ""]
            if len(nodes) == 3 and all(node in index_by_node for node in nodes):
                triangles.append([index_by_node[node] for node in nodes])
    return np.array(triangles, dtype=int)


def load_sentaurus_scalar(root: Path, field: str) -> np.ndarray:
    path = root / "fields" / f"{field}_region0.csv"
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        component_cols = [name for name in reader.fieldnames or [] if name != "node_id"]
        for row in reader:
            comps = [float(row[name]) for name in component_cols if row.get(name, "") != ""]
            values[int(row["node_id"])] = comps[0] if len(comps) == 1 else math.sqrt(sum(c * c for c in comps))
    return np.array([values[node_id] for node_id in sorted(values)], dtype=float)


def load_optional_sentaurus_scalar(root: Path, field: str, size: int) -> np.ndarray:
    path = root / "fields" / f"{field}_region0.csv"
    if not path.exists():
        return np.full(size, np.nan, dtype=float)
    return load_sentaurus_scalar(root, field)


def parse_vtk(path: Path) -> dict[str, Any]:
    lines = path.read_text().splitlines()
    points: list[tuple[float, float]] = []
    triangles: list[list[int]] = []
    scalars: dict[str, np.ndarray] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if not parts:
            i += 1
            continue
        if parts[0] == "POINTS":
            count = int(parts[1])
            tokens: list[str] = []
            i += 1
            while i < len(lines) and len(tokens) < 3 * count:
                tokens.extend(lines[i].split())
                i += 1
            if len(tokens) < 3 * count:
                raise ValueError(f"VTK POINTS block in {path} ended early")
            for offset in range(0, 3 * count, 3):
                x, y = tokens[offset], tokens[offset + 1]
                points.append((float(x) * 1.0e6, float(y) * 1.0e6))
            continue
        if parts[0] == "CELLS":
            count = int(parts[1])
            tokens: list[str] = []
            i += 1
            parsed = 0
            cursor = 0
            while i < len(lines) and parsed < count:
                tokens.extend(lines[i].split())
                while cursor < len(tokens) and parsed < count:
                    width = int(tokens[cursor])
                    if cursor + 1 + width > len(tokens):
                        break
                    cell = [width] + [int(item) for item in tokens[cursor + 1:cursor + 1 + width]]
                    cursor += 1 + width
                    parsed += 1
                    if cell and cell[0] == 3:
                        triangles.append(cell[1:4])
                i += 1
            if parsed < count:
                raise ValueError(f"VTK CELLS block in {path} ended early")
            continue
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA",
                }:
                    break
                values.extend(float(value) for value in next_parts)
                i += 1
            scalars[name] = np.array(values, dtype=float)
            continue
        i += 1
    xy = np.array(points, dtype=float)
    return {
        "x": xy[:, 0],
        "y": xy[:, 1],
        "triangles": np.array(triangles, dtype=int),
        "scalars": scalars,
    }


def nearest_values(
    src_x: np.ndarray,
    src_y: np.ndarray,
    values: np.ndarray,
    dst_x: np.ndarray,
    dst_y: np.ndarray,
) -> np.ndarray:
    matched = np.empty_like(dst_x, dtype=float)
    src = np.column_stack([src_x, src_y])
    for start in range(0, len(dst_x), 512):
        end = min(start + 512, len(dst_x))
        dst = np.column_stack([dst_x[start:end], dst_y[start:end]])
        dist2 = np.sum((dst[:, None, :] - src[None, :, :]) ** 2, axis=2)
        matched[start:end] = values[np.argmin(dist2, axis=1)]
    return matched


def nearest_indices_and_distances(
    src_x: np.ndarray,
    src_y: np.ndarray,
    dst_x: np.ndarray,
    dst_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    indices = np.empty_like(dst_x, dtype=int)
    distances = np.empty_like(dst_x, dtype=float)
    src = np.column_stack([src_x, src_y])
    for start in range(0, len(dst_x), 512):
        end = min(start + 512, len(dst_x))
        dst = np.column_stack([dst_x[start:end], dst_y[start:end]])
        dist2 = np.sum((dst[:, None, :] - src[None, :, :]) ** 2, axis=2)
        local_indices = np.argmin(dist2, axis=1)
        indices[start:end] = local_indices
        distances[start:end] = np.sqrt(dist2[np.arange(end - start), local_indices])
    return indices, distances


def mapped_values(values: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return values[indices]


def scalar_gradient_node_magnitudes(
    x_um: np.ndarray,
    y_um: np.ndarray,
    triangles: np.ndarray,
    values: np.ndarray,
) -> np.ndarray:
    edges: set[tuple[int, int]] = set()
    for tri in triangles:
        ids = [int(item) for item in tri]
        for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            edges.add((a, b) if a < b else (b, a))
    gradients = np.zeros_like(values, dtype=float)
    for a, b in edges:
        dx_m = (x_um[b] - x_um[a]) * 1.0e-6
        dy_m = (y_um[b] - y_um[a]) * 1.0e-6
        length = math.hypot(float(dx_m), float(dy_m))
        if length <= 0.0:
            continue
        grad_v_cm = abs(float(values[b] - values[a])) / length / 100.0
        gradients[a] = max(gradients[a], grad_v_cm)
        gradients[b] = max(gradients[b], grad_v_cm)
    return gradients


def stats(prefix: str, values: np.ndarray) -> dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return {f"{prefix}_{name}": math.nan for name in ("min", "p5", "median", "p95", "max")}
    return {
        f"{prefix}_min": float(np.nanmin(finite)),
        f"{prefix}_p5": float(np.nanpercentile(finite, 5)),
        f"{prefix}_median": float(np.nanmedian(finite)),
        f"{prefix}_p95": float(np.nanpercentile(finite, 95)),
        f"{prefix}_max": float(np.nanmax(finite)),
    }


def log_error(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    eps = 1.0e-300
    return np.abs(
        np.log10(np.maximum(np.abs(candidate), eps))
        - np.log10(np.maximum(np.abs(reference), eps))
    )


def relative_error(reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    return np.abs(candidate - reference) / np.maximum(np.abs(reference), 1.0e-30)


def implied_high_field_drive_v_cm(
    low_field_mobility_cm2_v_s: np.ndarray,
    limited_mobility_cm2_v_s: np.ndarray,
    saturation_velocity_m_s: float,
    beta: float,
) -> np.ndarray:
    low = np.asarray(low_field_mobility_cm2_v_s, dtype=float)
    limited = np.asarray(limited_mobility_cm2_v_s, dtype=float)
    ratio = limited / np.maximum(low, 1.0e-300)
    drive = np.full_like(low, np.nan, dtype=float)
    mask = (
        np.isfinite(low)
        & np.isfinite(limited)
        & (low > 0.0)
        & (limited > 0.0)
        & (ratio < 1.0)
    )
    if not np.any(mask):
        return drive
    ratio = np.clip(ratio[mask], 1.0e-300, 1.0)
    low_m2_v_s = low[mask] * 1.0e-4
    term = np.maximum(np.power(ratio, -beta) - 1.0, 0.0)
    drive_v_m = saturation_velocity_m_s / low_m2_v_s * np.power(term, 1.0 / beta)
    drive[mask] = drive_v_m / 100.0
    return drive


def ratio_stats(prefix: str, numerator: np.ndarray, denominator: np.ndarray) -> dict[str, float]:
    ratio = np.asarray(numerator, dtype=float) / np.maximum(np.asarray(denominator, dtype=float), 1.0e-300)
    return stats(prefix, ratio)


def finite_percentile(values: np.ndarray, percentile: float) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return math.nan
    return float(np.nanpercentile(finite, percentile))


def finite_max(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return math.nan
    return float(np.nanmax(finite))


def aggregate_drive_ratios(prefix: str, implied: np.ndarray, vela_drive: np.ndarray) -> dict[str, float]:
    implied_p95 = finite_percentile(implied, 95)
    vela_p95 = finite_percentile(vela_drive, 95)
    implied_max = finite_max(implied)
    vela_max = finite_max(vela_drive)
    return {
        f"{prefix}_p95_over_vela_drive_p95": implied_p95 / max(vela_p95, 1.0e-300),
        f"{prefix}_max_over_vela_drive_max": implied_max / max(vela_max, 1.0e-300),
    }


def local_diagnostic_row(
    bias: float,
    anchor: str,
    anchor_idx: int,
    idx: int,
    distance_um: float,
    sx: np.ndarray,
    sy: np.ndarray,
    sent_e_mob: np.ndarray,
    sent_h_mob: np.ndarray,
    vela_e_mob: np.ndarray,
    vela_h_mob: np.ndarray,
    sent_e_den: np.ndarray,
    sent_h_den: np.ndarray,
    vela_e_den: np.ndarray,
    vela_h_den: np.ndarray,
    vela_e_limiter: np.ndarray,
    vela_h_limiter: np.ndarray,
    sent_e_limiter_inferred: np.ndarray,
    sent_h_limiter_inferred: np.ndarray,
    vela_e_drive: np.ndarray,
    vela_h_drive: np.ndarray,
    sent_e_implied_drive: np.ndarray,
    sent_h_implied_drive: np.ndarray,
    nearest_vela_distance_um: np.ndarray,
    sent_potential: np.ndarray,
    sent_e_qf: np.ndarray,
    sent_h_qf: np.ndarray,
    sent_electric_field: np.ndarray,
    sent_grad_potential: np.ndarray,
    sent_grad_e_qf: np.ndarray,
    sent_grad_h_qf: np.ndarray,
    vela_potential: np.ndarray,
    vela_e_qf: np.ndarray,
    vela_h_qf: np.ndarray,
    vela_electric_field: np.ndarray,
    vela_grad_psi: np.ndarray,
    vela_grad_phin: np.ndarray,
    vela_grad_phip: np.ndarray,
) -> dict[str, Any]:
    return {
        "bias_V": bias,
        "anchor": anchor,
        "anchor_node": anchor_idx,
        "node_id": idx,
        "x_um": float(sx[idx]),
        "y_um": float(sy[idx]),
        "distance_from_anchor_um": distance_um,
        "nearest_vela_node_distance_um": float(nearest_vela_distance_um[idx]),
        "sentaurus_potential": float(sent_potential[idx]),
        "vela_potential": float(vela_potential[idx]),
        "sentaurus_eQuasiFermiPotential": float(sent_e_qf[idx]),
        "vela_electron_quasi_fermi": float(vela_e_qf[idx]),
        "sentaurus_hQuasiFermiPotential": float(sent_h_qf[idx]),
        "vela_hole_quasi_fermi": float(vela_h_qf[idx]),
        "sentaurus_electric_field_V_per_cm": float(sent_electric_field[idx]),
        "vela_electric_field_V_per_cm": float(vela_electric_field[idx]),
        "sentaurus_grad_potential_V_per_cm": float(sent_grad_potential[idx]),
        "vela_grad_psi_V_per_cm": float(vela_grad_psi[idx]),
        "sentaurus_grad_eQuasiFermiPotential_V_per_cm": float(sent_grad_e_qf[idx]),
        "vela_grad_phin_V_per_cm": float(vela_grad_phin[idx]),
        "sentaurus_grad_hQuasiFermiPotential_V_per_cm": float(sent_grad_h_qf[idx]),
        "vela_grad_phip_V_per_cm": float(vela_grad_phip[idx]),
        "sentaurus_eMobility": float(sent_e_mob[idx]),
        "vela_electron_mobility": float(vela_e_mob[idx]),
        "electron_mobility_relative_error": float(relative_error(sent_e_mob, vela_e_mob)[idx]),
        "sentaurus_hMobility": float(sent_h_mob[idx]),
        "vela_hole_mobility": float(vela_h_mob[idx]),
        "hole_mobility_relative_error": float(relative_error(sent_h_mob, vela_h_mob)[idx]),
        "sentaurus_eDensity": float(sent_e_den[idx]),
        "vela_electron_density": float(vela_e_den[idx]),
        "electron_density_log10_error": float(log_error(sent_e_den, vela_e_den)[idx]),
        "sentaurus_hDensity": float(sent_h_den[idx]),
        "vela_hole_density": float(vela_h_den[idx]),
        "hole_density_log10_error": float(log_error(sent_h_den, vela_h_den)[idx]),
        "vela_electron_limiter": float(vela_e_limiter[idx]),
        "sentaurus_electron_limiter_inferred_from_vela_low_field": float(sent_e_limiter_inferred[idx]),
        "vela_hole_limiter": float(vela_h_limiter[idx]),
        "sentaurus_hole_limiter_inferred_from_vela_low_field": float(sent_h_limiter_inferred[idx]),
        "vela_electron_high_field_drive_V_per_cm": float(vela_e_drive[idx]),
        "sentaurus_electron_implied_high_field_drive_V_per_cm": float(sent_e_implied_drive[idx]),
        "vela_hole_high_field_drive_V_per_cm": float(vela_h_drive[idx]),
        "sentaurus_hole_implied_high_field_drive_V_per_cm": float(sent_h_implied_drive[idx]),
    }


def top_error_location(
    x_um: np.ndarray,
    y_um: np.ndarray,
    err: np.ndarray,
) -> tuple[int, float, float, float]:
    finite = np.where(np.isfinite(err), err, -math.inf)
    if not len(finite) or not np.any(np.isfinite(err)):
        return -1, math.nan, math.nan, math.nan
    idx = int(np.argmax(finite))
    return idx, float(x_um[idx]), float(y_um[idx]), float(err[idx])


def neighborhood_indices(
    x_um: np.ndarray,
    y_um: np.ndarray,
    anchor_idx: int,
    radius_um: float,
    max_nodes: int,
) -> list[tuple[int, float]]:
    if anchor_idx < 0 or anchor_idx >= len(x_um):
        return []
    distances = np.hypot(x_um - x_um[anchor_idx], y_um - y_um[anchor_idx])
    if radius_um > 0.0:
        candidates = np.where(distances <= radius_um)[0]
    else:
        candidates = np.arange(len(x_um))
    order = sorted((int(idx) for idx in candidates), key=lambda idx: (float(distances[idx]), idx))
    return [(idx, float(distances[idx])) for idx in order[:max(1, max_nodes)]]


def append_top_error_neighborhood(
    rows: list[dict[str, Any]],
    bias: float,
    anchor: str,
    anchor_idx: int,
    sx: np.ndarray,
    sy: np.ndarray,
    sent_e_mob: np.ndarray,
    sent_h_mob: np.ndarray,
    vela_e_mob: np.ndarray,
    vela_h_mob: np.ndarray,
    sent_e_den: np.ndarray,
    sent_h_den: np.ndarray,
    vela_e_den: np.ndarray,
    vela_h_den: np.ndarray,
    vela_e_limiter: np.ndarray,
    vela_h_limiter: np.ndarray,
    sent_e_limiter_inferred: np.ndarray,
    sent_h_limiter_inferred: np.ndarray,
    vela_e_drive: np.ndarray,
    vela_h_drive: np.ndarray,
    sent_e_implied_drive: np.ndarray,
    sent_h_implied_drive: np.ndarray,
    nearest_vela_distance_um: np.ndarray,
    sent_potential: np.ndarray,
    sent_e_qf: np.ndarray,
    sent_h_qf: np.ndarray,
    sent_electric_field: np.ndarray,
    sent_grad_potential: np.ndarray,
    sent_grad_e_qf: np.ndarray,
    sent_grad_h_qf: np.ndarray,
    vela_potential: np.ndarray,
    vela_e_qf: np.ndarray,
    vela_h_qf: np.ndarray,
    vela_electric_field: np.ndarray,
    vela_grad_psi: np.ndarray,
    vela_grad_phin: np.ndarray,
    vela_grad_phip: np.ndarray,
    args: argparse.Namespace,
) -> None:
    for idx, distance_um in neighborhood_indices(
        sx,
        sy,
        anchor_idx,
        args.top_error_neighborhood_radius_um,
        args.top_error_neighborhood_max_nodes,
    ):
        rows.append(local_diagnostic_row(
            bias, anchor, anchor_idx, idx, distance_um, sx, sy,
            sent_e_mob, sent_h_mob, vela_e_mob, vela_h_mob,
            sent_e_den, sent_h_den, vela_e_den, vela_h_den,
            vela_e_limiter, vela_h_limiter,
            sent_e_limiter_inferred, sent_h_limiter_inferred,
            vela_e_drive, vela_h_drive,
            sent_e_implied_drive, sent_h_implied_drive,
            nearest_vela_distance_um,
            sent_potential, sent_e_qf, sent_h_qf, sent_electric_field,
            sent_grad_potential, sent_grad_e_qf, sent_grad_h_qf,
            vela_potential, vela_e_qf, vela_h_qf, vela_electric_field,
            vela_grad_psi, vela_grad_phin, vela_grad_phip,
        ))


def append_anchor_profiles(
    rows: list[dict[str, Any]],
    bias: float,
    anchor: str,
    anchor_idx: int,
    sx: np.ndarray,
    sy: np.ndarray,
    sent_e_mob: np.ndarray,
    sent_h_mob: np.ndarray,
    vela_e_mob: np.ndarray,
    vela_h_mob: np.ndarray,
    sent_e_den: np.ndarray,
    sent_h_den: np.ndarray,
    vela_e_den: np.ndarray,
    vela_h_den: np.ndarray,
    vela_e_limiter: np.ndarray,
    vela_h_limiter: np.ndarray,
    sent_e_limiter_inferred: np.ndarray,
    sent_h_limiter_inferred: np.ndarray,
    vela_e_drive: np.ndarray,
    vela_h_drive: np.ndarray,
    sent_e_implied_drive: np.ndarray,
    sent_h_implied_drive: np.ndarray,
    nearest_vela_distance_um: np.ndarray,
    sent_potential: np.ndarray,
    sent_e_qf: np.ndarray,
    sent_h_qf: np.ndarray,
    sent_electric_field: np.ndarray,
    sent_grad_potential: np.ndarray,
    sent_grad_e_qf: np.ndarray,
    sent_grad_h_qf: np.ndarray,
    vela_potential: np.ndarray,
    vela_e_qf: np.ndarray,
    vela_h_qf: np.ndarray,
    vela_electric_field: np.ndarray,
    vela_grad_psi: np.ndarray,
    vela_grad_phin: np.ndarray,
    vela_grad_phip: np.ndarray,
    args: argparse.Namespace,
) -> None:
    if anchor_idx < 0 or anchor_idx >= len(sx):
        return
    axes = [
        ("horizontal_y", np.where(np.abs(sy - sy[anchor_idx]) <= args.anchor_profile_tolerance_um)[0], sx),
        ("vertical_x", np.where(np.abs(sx - sx[anchor_idx]) <= args.anchor_profile_tolerance_um)[0], sy),
    ]
    for axis, indices, coordinate in axes:
        ordered = sorted((int(idx) for idx in indices), key=lambda idx: (float(coordinate[idx]), idx))
        for idx in ordered[:max(1, args.anchor_profile_max_nodes)]:
            row = local_diagnostic_row(
                bias, anchor, anchor_idx, idx,
                float(math.hypot(sx[idx] - sx[anchor_idx], sy[idx] - sy[anchor_idx])),
                sx, sy,
                sent_e_mob, sent_h_mob, vela_e_mob, vela_h_mob,
                sent_e_den, sent_h_den, vela_e_den, vela_h_den,
                vela_e_limiter, vela_h_limiter,
                sent_e_limiter_inferred, sent_h_limiter_inferred,
                vela_e_drive, vela_h_drive,
                sent_e_implied_drive, sent_h_implied_drive,
                nearest_vela_distance_um,
                sent_potential, sent_e_qf, sent_h_qf, sent_electric_field,
                sent_grad_potential, sent_grad_e_qf, sent_grad_h_qf,
                vela_potential, vela_e_qf, vela_h_qf, vela_electric_field,
                vela_grad_psi, vela_grad_phin, vela_grad_phip,
            )
            row["profile_axis"] = axis
            row["profile_coordinate_um"] = float(coordinate[idx])
            rows.append(row)


def clean_for_csv(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def compare_bias(
    bias: float,
    sentaurus_dir: Path | None,
    vtk_path: Path | None,
    args: argparse.Namespace,
    neighborhood_rows: list[dict[str, Any]] | None = None,
    profile_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    row: dict[str, Any] = {
        "bias_V": bias,
        "status": "ok",
        "mobility_model_electron_saturation_velocity_m_s": args.electron_saturation_velocity_m_s,
        "mobility_model_hole_saturation_velocity_m_s": args.hole_saturation_velocity_m_s,
        "mobility_model_electron_field_beta": args.electron_field_beta,
        "mobility_model_hole_field_beta": args.hole_field_beta,
    }
    if sentaurus_dir is None or vtk_path is None:
        row["status"] = "missing_input"
        return row

    sx, sy, node_ids = load_sentaurus_mesh(sentaurus_dir)
    sent_triangles = load_sentaurus_triangles(sentaurus_dir, node_ids)
    vtk = parse_vtk(vtk_path)
    vx = np.asarray(vtk["x"], dtype=float)
    vy = np.asarray(vtk["y"], dtype=float)
    scalars = vtk["scalars"]
    nearest_idx, nearest_distance_um = nearest_indices_and_distances(vx, vy, sx, sy)

    sent_e_mob = load_sentaurus_scalar(sentaurus_dir, "eMobility")
    sent_h_mob = load_sentaurus_scalar(sentaurus_dir, "hMobility")
    sent_e_den = load_sentaurus_scalar(sentaurus_dir, "eDensity")
    sent_h_den = load_sentaurus_scalar(sentaurus_dir, "hDensity")
    sent_potential = load_optional_sentaurus_scalar(sentaurus_dir, "ElectrostaticPotential", len(sx))
    sent_e_qf = load_optional_sentaurus_scalar(sentaurus_dir, "eQuasiFermiPotential", len(sx))
    sent_h_qf = load_optional_sentaurus_scalar(sentaurus_dir, "hQuasiFermiPotential", len(sx))
    sent_electric_field = load_optional_sentaurus_scalar(sentaurus_dir, "ElectricField", len(sx))

    vela_e_mob = mapped_values(scalars["ElectronMobility"] * 1.0e4, nearest_idx)
    vela_h_mob = mapped_values(scalars["HoleMobility"] * 1.0e4, nearest_idx)
    vela_e_low = mapped_values(scalars.get("ElectronLowFieldMobility", scalars["ElectronMobility"]) * 1.0e4, nearest_idx)
    vela_h_low = mapped_values(scalars.get("HoleLowFieldMobility", scalars["HoleMobility"]) * 1.0e4, nearest_idx)
    vela_e_limiter = mapped_values(scalars.get("ElectronMobilityLimiter", np.ones_like(vx)), nearest_idx)
    vela_h_limiter = mapped_values(scalars.get("HoleMobilityLimiter", np.ones_like(vx)), nearest_idx)
    vela_e_drive = mapped_values(scalars.get("ElectronHighFieldDrive", np.zeros_like(vx)), nearest_idx)
    vela_h_drive = mapped_values(scalars.get("HoleHighFieldDrive", np.zeros_like(vx)), nearest_idx)
    vela_e_den = mapped_values(scalars["Electrons"] * 1.0e-6, nearest_idx)
    vela_h_den = mapped_values(scalars["Holes"] * 1.0e-6, nearest_idx)
    vela_potential = mapped_values(scalars.get("Potential", np.full_like(vx, np.nan)), nearest_idx)
    vela_e_qf = mapped_values(scalars.get("ElectronQuasiFermi", np.full_like(vx, np.nan)), nearest_idx)
    vela_h_qf = mapped_values(scalars.get("HoleQuasiFermi", np.full_like(vx, np.nan)), nearest_idx)
    vela_electric_field = mapped_values(scalars.get("ElectricField", np.full_like(vx, np.nan)), nearest_idx)
    sent_grad_potential = scalar_gradient_node_magnitudes(sx, sy, sent_triangles, sent_potential)
    sent_grad_e_qf = scalar_gradient_node_magnitudes(sx, sy, sent_triangles, sent_e_qf)
    sent_grad_h_qf = scalar_gradient_node_magnitudes(sx, sy, sent_triangles, sent_h_qf)
    vela_grad_psi = mapped_values(
        scalar_gradient_node_magnitudes(vx, vy, vtk["triangles"], scalars.get("Potential", np.full_like(vx, np.nan))),
        nearest_idx,
    )
    vela_grad_phin = mapped_values(
        scalar_gradient_node_magnitudes(vx, vy, vtk["triangles"], scalars.get("ElectronQuasiFermi", np.full_like(vx, np.nan))),
        nearest_idx,
    )
    vela_grad_phip = mapped_values(
        scalar_gradient_node_magnitudes(vx, vy, vtk["triangles"], scalars.get("HoleQuasiFermi", np.full_like(vx, np.nan))),
        nearest_idx,
    )

    row.update(stats("sentaurus_eMobility", sent_e_mob))
    row.update(stats("sentaurus_hMobility", sent_h_mob))
    row.update(stats("sentaurus_to_vela_nearest_distance_um", nearest_distance_um))
    row.update(stats("sentaurus_electric_field_V_per_cm", sent_electric_field))
    row.update(stats("sentaurus_grad_potential_V_per_cm", sent_grad_potential))
    row.update(stats("sentaurus_grad_eQuasiFermiPotential_V_per_cm", sent_grad_e_qf))
    row.update(stats("sentaurus_grad_hQuasiFermiPotential_V_per_cm", sent_grad_h_qf))
    row.update(stats("vela_electron_mobility", vela_e_mob))
    row.update(stats("vela_hole_mobility", vela_h_mob))
    row.update(stats("vela_electron_low_field_mobility", vela_e_low))
    row.update(stats("vela_hole_low_field_mobility", vela_h_low))
    row.update(stats("vela_electron_limiter", vela_e_limiter))
    row.update(stats("vela_hole_limiter", vela_h_limiter))
    row.update(stats("vela_electron_high_field_drive_V_per_cm", vela_e_drive))
    row.update(stats("vela_hole_high_field_drive_V_per_cm", vela_h_drive))

    row.update(stats("vela_grad_psi_V_per_cm", vela_grad_psi))
    row.update(stats("vela_grad_phin_V_per_cm", vela_grad_phin))
    row.update(stats("vela_grad_phip_V_per_cm", vela_grad_phip))

    e_mob_idx, e_mob_x, e_mob_y, e_mob_err = top_error_location(
        sx, sy, relative_error(sent_e_mob, vela_e_mob))
    h_mob_idx, h_mob_x, h_mob_y, h_mob_err = top_error_location(
        sx, sy, relative_error(sent_h_mob, vela_h_mob))
    e_den_idx, e_den_x, e_den_y, e_den_err = top_error_location(
        sx, sy, log_error(sent_e_den, vela_e_den))
    h_den_idx, h_den_x, h_den_y, h_den_err = top_error_location(
        sx, sy, log_error(sent_h_den, vela_h_den))
    e_lim_idx, e_lim_x, e_lim_y, e_lim_strength = top_error_location(
        sx, sy, np.abs(1.0 - vela_e_limiter))
    h_lim_idx, h_lim_x, h_lim_y, h_lim_strength = top_error_location(
        sx, sy, np.abs(1.0 - vela_h_limiter))

    row.update({
        "electron_mobility_top_error_node": e_mob_idx,
        "electron_mobility_top_error_x_um": e_mob_x,
        "electron_mobility_top_error_y_um": e_mob_y,
        "electron_mobility_top_relative_error": e_mob_err,
        "hole_mobility_top_error_node": h_mob_idx,
        "hole_mobility_top_error_x_um": h_mob_x,
        "hole_mobility_top_error_y_um": h_mob_y,
        "hole_mobility_top_relative_error": h_mob_err,
        "electron_density_top_error_node": e_den_idx,
        "electron_density_top_error_x_um": e_den_x,
        "electron_density_top_error_y_um": e_den_y,
        "electron_density_top_log10_error": e_den_err,
        "hole_density_top_error_node": h_den_idx,
        "hole_density_top_error_x_um": h_den_x,
        "hole_density_top_error_y_um": h_den_y,
        "hole_density_top_log10_error": h_den_err,
        "electron_limiter_top_node": e_lim_idx,
        "electron_limiter_top_x_um": e_lim_x,
        "electron_limiter_top_y_um": e_lim_y,
        "electron_limiter_top_strength": e_lim_strength,
        "hole_limiter_top_node": h_lim_idx,
        "hole_limiter_top_x_um": h_lim_x,
        "hole_limiter_top_y_um": h_lim_y,
        "hole_limiter_top_strength": h_lim_strength,
    })
    if e_mob_idx >= 0 and e_den_idx >= 0:
        row["electron_density_to_mobility_top_error_distance_um"] = float(
            math.hypot(e_mob_x - e_den_x, e_mob_y - e_den_y))
    if e_mob_idx >= 0 and e_lim_idx >= 0:
        row["electron_limiter_to_mobility_top_error_distance_um"] = float(
            math.hypot(e_mob_x - e_lim_x, e_mob_y - e_lim_y))
    if e_den_idx >= 0 and e_lim_idx >= 0:
        row["electron_limiter_to_density_top_error_distance_um"] = float(
            math.hypot(e_den_x - e_lim_x, e_den_y - e_lim_y))
    if h_mob_idx >= 0 and h_den_idx >= 0:
        row["hole_density_to_mobility_top_error_distance_um"] = float(
            math.hypot(h_mob_x - h_den_x, h_mob_y - h_den_y))
    if h_mob_idx >= 0 and h_lim_idx >= 0:
        row["hole_limiter_to_mobility_top_error_distance_um"] = float(
            math.hypot(h_mob_x - h_lim_x, h_mob_y - h_lim_y))
    if h_den_idx >= 0 and h_lim_idx >= 0:
        row["hole_limiter_to_density_top_error_distance_um"] = float(
            math.hypot(h_den_x - h_lim_x, h_den_y - h_lim_y))

    sent_e_limiter_inferred = sent_e_mob / np.maximum(vela_e_low, 1.0e-300)
    sent_h_limiter_inferred = sent_h_mob / np.maximum(vela_h_low, 1.0e-300)
    row.update(stats("sentaurus_electron_limiter_inferred_from_vela_low_field", sent_e_limiter_inferred))
    row.update(stats("sentaurus_hole_limiter_inferred_from_vela_low_field", sent_h_limiter_inferred))
    sent_e_implied_drive = implied_high_field_drive_v_cm(
        vela_e_low,
        sent_e_mob,
        saturation_velocity_m_s=args.electron_saturation_velocity_m_s,
        beta=args.electron_field_beta)
    sent_h_implied_drive = implied_high_field_drive_v_cm(
        vela_h_low,
        sent_h_mob,
        saturation_velocity_m_s=args.hole_saturation_velocity_m_s,
        beta=args.hole_field_beta)
    row.update(stats("sentaurus_electron_implied_high_field_drive_V_per_cm", sent_e_implied_drive))
    row.update(stats("sentaurus_hole_implied_high_field_drive_V_per_cm", sent_h_implied_drive))
    row.update(ratio_stats("electron_implied_to_vela_drive_ratio", sent_e_implied_drive, vela_e_drive))
    row.update(ratio_stats("hole_implied_to_vela_drive_ratio", sent_h_implied_drive, vela_h_drive))
    row.update(aggregate_drive_ratios(
        "electron_implied_drive", sent_e_implied_drive, vela_e_drive))
    row.update(aggregate_drive_ratios(
        "hole_implied_drive", sent_h_implied_drive, vela_h_drive))
    row["electron_limiter_weaker_than_sentaurus_median"] = (
        float(np.nanmedian(vela_e_limiter - sent_e_limiter_inferred)))
    row["hole_limiter_weaker_than_sentaurus_median"] = (
        float(np.nanmedian(vela_h_limiter - sent_h_limiter_inferred)))
    if neighborhood_rows is not None:
        for anchor, anchor_idx in (
            ("electron_mobility_top_error", e_mob_idx),
            ("electron_density_top_error", e_den_idx),
            ("electron_limiter_top", e_lim_idx),
            ("hole_mobility_top_error", h_mob_idx),
            ("hole_density_top_error", h_den_idx),
            ("hole_limiter_top", h_lim_idx),
        ):
            append_top_error_neighborhood(
                neighborhood_rows,
                bias,
                anchor,
                anchor_idx,
                sx,
                sy,
                sent_e_mob,
                sent_h_mob,
                vela_e_mob,
                vela_h_mob,
                sent_e_den,
                sent_h_den,
                vela_e_den,
                vela_h_den,
                vela_e_limiter,
                vela_h_limiter,
                sent_e_limiter_inferred,
                sent_h_limiter_inferred,
                vela_e_drive,
                vela_h_drive,
                sent_e_implied_drive,
                sent_h_implied_drive,
                nearest_distance_um,
                sent_potential,
                sent_e_qf,
                sent_h_qf,
                sent_electric_field,
                sent_grad_potential,
                sent_grad_e_qf,
                sent_grad_h_qf,
                vela_potential,
                vela_e_qf,
                vela_h_qf,
                vela_electric_field,
                vela_grad_psi,
                vela_grad_phin,
                vela_grad_phip,
                args,
            )
    if profile_rows is not None:
        for anchor, anchor_idx in (
            ("electron_mobility_top_error", e_mob_idx),
            ("electron_density_top_error", e_den_idx),
            ("electron_limiter_top", e_lim_idx),
            ("hole_mobility_top_error", h_mob_idx),
            ("hole_density_top_error", h_den_idx),
            ("hole_limiter_top", h_lim_idx),
        ):
            append_anchor_profiles(
                profile_rows,
                bias,
                anchor,
                anchor_idx,
                sx,
                sy,
                sent_e_mob,
                sent_h_mob,
                vela_e_mob,
                vela_h_mob,
                sent_e_den,
                sent_h_den,
                vela_e_den,
                vela_h_den,
                vela_e_limiter,
                vela_h_limiter,
                sent_e_limiter_inferred,
                sent_h_limiter_inferred,
                vela_e_drive,
                vela_h_drive,
                sent_e_implied_drive,
                sent_h_implied_drive,
                nearest_distance_um,
                sent_potential,
                sent_e_qf,
                sent_h_qf,
                sent_electric_field,
                sent_grad_potential,
                sent_grad_e_qf,
                sent_grad_h_qf,
                vela_potential,
                vela_e_qf,
                vela_h_qf,
                vela_electric_field,
                vela_grad_psi,
                vela_grad_phin,
                vela_grad_phip,
                args,
            )
    return row


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: clean_for_csv(row.get(key, "")) for key in fields})


def main() -> int:
    args = parse_args()
    biases = parse_biases(args.biases)
    sentaurus_exports = discover_sentaurus_exports(args.sentaurus_root, biases)
    vela_vtks = discover_vela_vtks(args.vela_vtk_root)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    neighborhood_rows: list[dict[str, Any]] = []
    profile_rows: list[dict[str, Any]] = []
    rows = [
        compare_bias(
            bias,
            sentaurus_exports.get(bias_key(bias)),
            vela_vtks.get(bias_key(bias)),
            args,
            neighborhood_rows,
            profile_rows,
        )
        for bias in biases
    ]
    write_csv(args.out_dir / "mobility_decomposition_summary.csv", rows)
    write_csv(args.out_dir / "top_error_neighborhood.csv", neighborhood_rows)
    write_csv(args.out_dir / "quasi_fermi_anchor_profiles.csv", profile_rows)
    (args.out_dir / "mobility_decomposition_summary.json").write_text(
        json.dumps(rows, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    (args.out_dir / "README.md").write_text(
        "# PN2D BV Mobility Debug\n\n"
        "Generated by `scripts/diagnose_pn2d_bv_mobility.py`.\n"
        "Mobility values are reported in cm^2/V/s and high-field drives in V/cm.\n"
        "`top_error_neighborhood.csv` lists local nodes around mobility, density, and limiter error anchors.\n"
        "`quasi_fermi_anchor_profiles.csv` lists horizontal and vertical quasi-Fermi profiles through those anchors.\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
