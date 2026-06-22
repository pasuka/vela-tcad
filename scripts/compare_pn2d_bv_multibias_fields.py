#!/usr/bin/env python3
"""Compare PN2D BV multibias Sentaurus exports against Vela VTK outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

import numpy as np


FIELD_SPECS: dict[str, dict[str, Any]] = {
    "potential": {
        "sentaurus": ["ElectrostaticPotential", "Potential"],
        "vela": "Potential",
        "scale": 1.0,
        "metric": "rms",
    },
    "electric_field": {
        "sentaurus": ["ElectricField"],
        "vela": "ElectricField",
        "scale": 1.0,
        "metric": "relative_p95",
    },
    "electron_density": {
        "sentaurus": ["eDensity"],
        "vela": "Electrons",
        "scale": 1.0e-6,
        "metric": "log10_p95",
    },
    "hole_density": {
        "sentaurus": ["hDensity"],
        "vela": "Holes",
        "scale": 1.0e-6,
        "metric": "log10_p95",
    },
    "srh_recombination": {
        "sentaurus": ["SRHRecombination"],
        "vela": "SRHRecombination",
        "scale": 1.0e-6,
        "metric": "log10_p95",
    },
    "avalanche_generation": {
        "sentaurus": ["AvalancheGeneration", "ImpactIonization"],
        "vela": "AvalancheGeneration",
        "scale": 1.0e-6,
        "metric": "log10_p95",
    },
    "electron_mobility": {
        "sentaurus": ["eMobility"],
        "vela": "ElectronMobility",
        "scale": 1.0e4,
        "metric": "relative_p95",
    },
    "hole_mobility": {
        "sentaurus": ["hMobility"],
        "vela": "HoleMobility",
        "scale": 1.0e4,
        "metric": "relative_p95",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--vela-vtk-root", type=Path, required=True)
    parser.add_argument("--curve-reference", type=Path, required=True)
    parser.add_argument("--curve-candidate", type=Path)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--current-floor",
        type=float,
        default=1.0e-25,
        help="Disable decade current error when either curve current is below this absolute value.",
    )
    parser.add_argument(
        "--avalanche-floor",
        type=float,
        default=1.0,
        help="Treat avalanche p99 values below this absolute field value as near-floor diagnostics.",
    )
    parser.add_argument(
        "--biases",
        default="0,-0.5,-2,-5,-10,-20",
        help="Comma-separated BV biases to compare in volts.",
    )
    parser.add_argument(
        "--quantities",
        default=",".join(FIELD_SPECS),
        help="Comma-separated quantities to compare.",
    )
    return parser.parse_args()


def parse_csv_list(raw: str) -> list[str]:
    return [item.strip() for item in raw.split(",") if item.strip()]


def parse_biases(raw: str) -> list[float]:
    return [float(item) for item in parse_csv_list(raw)]


def bias_key(value: float) -> float:
    return round(value, 12)


def bias_token(bias: float) -> str:
    token = f"{bias:g}".replace("-", "m").replace("+", "").replace(".", "p")
    return token


def signed_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("+", "").replace(".", "p")


def discover_sentaurus_exports(root: Path, biases: list[float]) -> dict[float, Path]:
    result: dict[float, Path] = {}
    for bias in biases:
        candidates = [
            root / f"sentaurus_{signed_bias_token(bias)}v",
            root / f"sentaurus_{bias_token(bias)}v",
            root / f"sentaurus_{bias_token(abs(bias))}v",
            root / f"sentaurus_{bias:g}v",
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
        if match:
            key = bias_key(float(match.group("bias")))
            mtime = path.stat().st_mtime
            if key not in result or mtime >= mtimes[key]:
                result[key] = path
                mtimes[key] = mtime
    return result


def load_sentaurus_mesh(root: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nodes: dict[int, tuple[float, float]] = {}
    with (root / "nodes.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            nodes[int(row["id"])] = (float(row["x_um"]), float(row["y_um"]))
    triangles: list[list[int]] = []
    with (root / "elements.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            triangles.append([int(row["node0"]), int(row["node1"]), int(row["node2"])])
    ids = np.array(sorted(nodes), dtype=int)
    xy = np.array([nodes[int(node_id)] for node_id in ids], dtype=float)
    return xy[:, 0], xy[:, 1], np.array(triangles, dtype=int)


def load_sentaurus_node_ids(root: Path) -> list[int]:
    with (root / "nodes.csv").open(newline="") as handle:
        return sorted(int(row["id"]) for row in csv.DictReader(handle))


def parse_node_list(raw: str) -> list[int]:
    return [int(item) for item in re.split(r"[;\s,]+", raw.strip()) if item]


def load_region_masks(root: Path, x_um: np.ndarray, y_um: np.ndarray) -> dict[str, np.ndarray]:
    node_ids = load_sentaurus_node_ids(root)
    node_index = {node_id: index for index, node_id in enumerate(node_ids)}
    count = len(node_ids)
    masks = {
        "p_contact": np.zeros(count, dtype=bool),
        "n_contact": np.zeros(count, dtype=bool),
        "junction": np.zeros(count, dtype=bool),
        "bulk": np.ones(count, dtype=bool),
    }

    contacts_path = root / "contacts.csv"
    doping_path = root / "doping.csv"
    net_doping: dict[int, float] = {}
    if doping_path.exists():
        with doping_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                net_doping[int(row["node_id"])] = (
                    float(row.get("donors_cm3", 0.0) or 0.0)
                    - float(row.get("acceptors_cm3", 0.0) or 0.0)
                )

    if contacts_path.exists():
        with contacts_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                nodes = [node for node in parse_node_list(row.get("node_ids", "")) if node in node_index]
                if not nodes:
                    continue
                name = row.get("name", "").lower()
                if "anode" in name or name.startswith("p"):
                    target = "p_contact"
                elif "cathode" in name or name.startswith("n"):
                    target = "n_contact"
                else:
                    mean_net = sum(net_doping.get(node, 0.0) for node in nodes) / len(nodes)
                    target = "n_contact" if mean_net >= 0.0 else "p_contact"
                for node in nodes:
                    masks[target][node_index[node]] = True

    elements_path = root / "elements.csv"
    if elements_path.exists() and net_doping:
        with elements_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                tri = [int(row[name]) for name in ("node0", "node1", "node2")]
                signs = [math.copysign(1.0, net_doping.get(node, 0.0))
                         if net_doping.get(node, 0.0) != 0.0 else 0.0
                         for node in tri]
                if min(signs) < 0.0 < max(signs):
                    for node in tri:
                        if node in node_index:
                            masks["junction"][node_index[node]] = True

    if not np.any(masks["junction"]) and count:
        x_mid = float(np.nanmedian(x_um))
        width = max(float(np.nanmax(x_um) - np.nanmin(x_um)), 1.0e-9)
        masks["junction"] = np.abs(x_um - x_mid) <= 0.1 * width

    occupied = masks["p_contact"] | masks["n_contact"] | masks["junction"]
    masks["bulk"] = ~occupied
    return masks


def load_sentaurus_scalar(root: Path, aliases: list[str]) -> tuple[str, np.ndarray] | None:
    for field in aliases:
        for path in (root / "fields" / f"{field}_region0.csv", root / f"{field}_region0.csv"):
            if not path.exists():
                continue
            values: dict[int, float] = {}
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                component_cols = [name for name in reader.fieldnames or [] if name != "node_id"]
                for row in reader:
                    comps = [float(row[name]) for name in component_cols if row.get(name, "") != ""]
                    value = comps[0] if len(comps) == 1 else math.sqrt(sum(c * c for c in comps))
                    values[int(row["node_id"])] = value
            return field, np.array([values[node_id] for node_id in sorted(values)], dtype=float)
    return None


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
            for raw in lines[i + 1:i + 1 + count]:
                x, y, *_ = raw.split()
                points.append((float(x) * 1.0e6, float(y) * 1.0e6))
            i += 1 + count
            continue
        if parts[0] == "CELLS":
            count = int(parts[1])
            for raw in lines[i + 1:i + 1 + count]:
                cell = [int(item) for item in raw.split()]
                if cell and cell[0] == 3:
                    triangles.append(cell[1:4])
            i += 1 + count
            continue
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            i += 2
            values: list[float] = []
            while i < len(lines):
                next_parts = lines[i].split()
                if not next_parts or next_parts[0] in {
                    "SCALARS",
                    "VECTORS",
                    "FIELD",
                    "CELL_DATA",
                    "POINT_DATA",
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


def log10_p95(reference: np.ndarray, candidate: np.ndarray) -> float:
    eps = 1.0e-300
    mask = np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(mask):
        return math.nan
    err = np.abs(np.log10(np.maximum(np.abs(candidate[mask]), eps))
                 - np.log10(np.maximum(np.abs(reference[mask]), eps)))
    return float(np.nanpercentile(err, 95))


def relative_p95(reference: np.ndarray, candidate: np.ndarray) -> float:
    mask = np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(mask):
        return math.nan
    denom = np.maximum(np.abs(reference[mask]), 1.0e-30)
    err = np.abs(candidate[mask] - reference[mask]) / denom
    return float(np.nanpercentile(err, 95))


def rms(reference: np.ndarray, candidate: np.ndarray) -> float:
    mask = np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(mask):
        return math.nan
    diff = candidate[mask] - reference[mask]
    return float(math.sqrt(float(np.mean(diff * diff))))


def centerline_indices(y: np.ndarray) -> np.ndarray:
    y_mid = float(np.nanmedian(y))
    span = max(float(np.nanmax(y) - np.nanmin(y)), 1.0e-9)
    tol = max(0.01 * span, 0.02)
    return np.where(np.abs(y - y_mid) <= tol)[0]


def metric_value(metric: str, reference: np.ndarray, candidate: np.ndarray) -> float:
    if metric == "rms":
        return rms(reference, candidate)
    if metric == "relative_p95":
        return relative_p95(reference, candidate)
    if metric == "log10_p95":
        return log10_p95(reference, candidate)
    raise ValueError(f"unknown metric {metric}")


def per_node_error(metric: str, reference: np.ndarray, candidate: np.ndarray) -> np.ndarray:
    err = np.full_like(reference, np.nan, dtype=float)
    mask = np.isfinite(reference) & np.isfinite(candidate)
    if not np.any(mask):
        return err
    if metric == "rms":
        err[mask] = np.abs(candidate[mask] - reference[mask])
    elif metric == "relative_p95":
        err[mask] = np.abs(candidate[mask] - reference[mask]) / np.maximum(np.abs(reference[mask]), 1.0e-30)
    elif metric == "log10_p95":
        eps = 1.0e-300
        err[mask] = np.abs(
            np.log10(np.maximum(np.abs(candidate[mask]), eps))
            - np.log10(np.maximum(np.abs(reference[mask]), eps))
        )
    else:
        raise ValueError(f"unknown metric {metric}")
    return err


def regional_metric(metric: str, reference: np.ndarray, candidate: np.ndarray, mask: np.ndarray) -> float:
    if not np.any(mask):
        return math.nan
    return metric_value(metric, reference[mask], candidate[mask])


def load_curve_points(path: Path) -> list[tuple[float, float]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    if not rows:
        return []
    fields = reader.fieldnames or []
    bias_col = next((name for name in fields if name in {"bias_V", "voltage", "V"}), fields[0])
    preferred_current_cols = ["current_total_A_per_um", "current_total", "current"]
    current_col = next((name for name in preferred_current_cols if name in fields), fields[-1])
    points: list[tuple[float, float]] = []
    for row in rows:
        if "converged" in row and str(row.get("converged", "")).strip() not in {"1", "true", "True"}:
            continue
        try:
            points.append((float(row[bias_col]), float(row[current_col])))
        except (TypeError, ValueError):
            continue
    return sorted(points, key=lambda item: item[0])


def interpolate_curve(points: list[tuple[float, float]], bias: float) -> float | None:
    if not points:
        return None
    for point_bias, current in points:
        if abs(point_bias - bias) <= 1.0e-12:
            return current
    for (b0, i0), (b1, i1) in zip(points, points[1:]):
        if b0 <= bias <= b1:
            weight = (bias - b0) / (b1 - b0)
            return i0 + weight * (i1 - i0)
    return None


def curve_compare_rows(
    biases: list[float],
    reference: list[tuple[float, float]],
    candidate: list[tuple[float, float]],
    current_floor: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias in biases:
        ref = interpolate_curve(reference, bias)
        cand = interpolate_curve(candidate, bias)
        status = "ok"
        log_ratio = math.nan
        if ref is None:
            status = "missing_reference"
        elif cand is None:
            status = "missing_candidate"
        elif abs(ref) < current_floor or abs(cand) < current_floor:
            status = "below_current_floor"
        elif ref == 0.0 or cand == 0.0:
            status = "zero_current"
        else:
            log_ratio = math.log10(abs(cand) / abs(ref))
        rows.append({
            "bias_V": bias,
            "sentaurus_current_A": ref,
            "vela_current_A_per_um": cand,
            "log10_abs_ratio": log_ratio,
            "abs_log10_error": abs(log_ratio) if math.isfinite(log_ratio) else math.nan,
            "status": status,
        })
    return rows


def evaluate_current_window(
    name: str,
    bias_min: float,
    bias_max: float,
    min_ratio: float,
    max_ratio: float,
    reference: list[tuple[float, float]],
    candidate: list[tuple[float, float]],
    current_floor: float,
) -> dict[str, Any]:
    sample_biases = sorted({bias_min, (bias_min + bias_max) / 2.0, bias_max})
    ratios: list[float] = []
    missing_biases: list[float] = []
    floor_biases: list[float] = []
    for bias in sample_biases:
        ref = interpolate_curve(reference, bias)
        cand = interpolate_curve(candidate, bias)
        if ref is None or cand is None:
            missing_biases.append(bias)
            continue
        if abs(ref) < current_floor or abs(cand) < current_floor or ref == 0.0:
            floor_biases.append(bias)
            continue
        ratios.append(abs(cand) / abs(ref))

    status = "pass"
    reason = ""
    if missing_biases:
        status = "not_evaluated"
        reason = "missing reference or candidate points in requested bias window"
    elif floor_biases:
        status = "not_evaluated"
        reason = "reference or candidate current below floor in requested bias window"
    elif not ratios:
        status = "not_evaluated"
        reason = "no current ratios available in requested bias window"
    elif min(ratios) < min_ratio or max(ratios) > max_ratio:
        status = "fail"
        reason = "current ratio outside configured band"

    return {
        "name": name,
        "bias_min": bias_min,
        "bias_max": bias_max,
        "min_ratio": min_ratio,
        "max_ratio": max_ratio,
        "sample_biases": sample_biases,
        "ratios": ratios,
        "status": status,
        "reason": reason,
    }


def build_bv_trend_summary(
    biases: list[float],
    curve_reference: list[tuple[float, float]],
    curve_candidate: list[tuple[float, float]],
    field_rows: list[dict[str, Any]],
    current_floor: float,
) -> dict[str, Any]:
    _ = field_rows
    current_windows = [
        evaluate_current_window(
            "mid_bias_current_band",
            -13.2,
            -13.0,
            0.6,
            1.4,
            curve_reference,
            curve_candidate,
            current_floor,
        )
    ]
    return {
        "biases_compared": biases,
        "max_field_monotonic": None,
        "max_field_monotonic_status": "not_evaluated",
        "current_windows": current_windows,
        "high_bias_knee_shape": {
            "status": "diagnostic",
            "reason": "current evidence does not support full -20 V absolute-current parity",
        },
    }


def compare_bias(
    bias: float,
    sentaurus_dir: Path | None,
    vtk_path: Path | None,
    quantities: list[str],
    avalanche_floor: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if sentaurus_dir is None or vtk_path is None:
        for quantity in quantities:
            rows.append({"bias_V": bias, "quantity": quantity, "status": "missing_input"})
        return rows

    sx, sy, _ = load_sentaurus_mesh(sentaurus_dir)
    node_ids = load_sentaurus_node_ids(sentaurus_dir)
    region_masks = load_region_masks(sentaurus_dir, sx, sy)
    vtk = parse_vtk(vtk_path)
    vx = np.asarray(vtk["x"], dtype=float)
    vy = np.asarray(vtk["y"], dtype=float)
    center_ids = centerline_indices(sy)

    for quantity in quantities:
        spec = FIELD_SPECS[quantity]
        loaded = load_sentaurus_scalar(sentaurus_dir, spec["sentaurus"])
        vraw = vtk["scalars"].get(spec["vela"])
        if loaded is None or vraw is None:
            rows.append({"bias_V": bias, "quantity": quantity, "status": "missing_field"})
            continue
        sentaurus_field, svals = loaded
        vvals = nearest_values(vx, vy, np.asarray(vraw, dtype=float) * spec["scale"], sx, sy)
        metric = metric_value(spec["metric"], svals, vvals)
        center_metric = (
            metric_value(spec["metric"], svals[center_ids], vvals[center_ids])
            if len(center_ids) > 0
            else math.nan
        )
        node_error = per_node_error(spec["metric"], svals, vvals)
        finite_error = np.where(np.isfinite(node_error), node_error, -math.inf)
        top_index = int(np.argmax(finite_error)) if len(finite_error) and np.any(np.isfinite(node_error)) else -1
        row = {
            "bias_V": bias,
            "quantity": quantity,
            "status": "ok",
            "sentaurus_field": sentaurus_field,
            "metric": spec["metric"],
            "field_error": metric,
            "centerline_error": center_metric,
            "top_error_node_id": node_ids[top_index] if top_index >= 0 and top_index < len(node_ids) else "",
            "top_error_x_um": float(sx[top_index]) if top_index >= 0 else math.nan,
            "top_error_y_um": float(sy[top_index]) if top_index >= 0 else math.nan,
            "top_error_value": float(node_error[top_index]) if top_index >= 0 else math.nan,
            "p_contact_error": regional_metric(spec["metric"], svals, vvals, region_masks["p_contact"]),
            "n_contact_error": regional_metric(spec["metric"], svals, vvals, region_masks["n_contact"]),
            "junction_error": regional_metric(spec["metric"], svals, vvals, region_masks["junction"]),
            "bulk_error": regional_metric(spec["metric"], svals, vvals, region_masks["bulk"]),
        }
        if quantity == "avalanche_generation":
            s_abs = np.abs(svals)
            v_abs = np.abs(vvals)
            s_peak = int(np.nanargmax(s_abs)) if len(svals) else -1
            v_peak = int(np.nanargmax(v_abs)) if len(vvals) else -1
            s_threshold = float(np.nanpercentile(s_abs, 99)) if len(s_abs) else math.nan
            v_threshold = float(np.nanpercentile(v_abs, 99)) if len(v_abs) else math.nan
            row.update({
                "sentaurus_peak_x_um": float(sx[s_peak]) if s_peak >= 0 else math.nan,
                "sentaurus_peak_y_um": float(sy[s_peak]) if s_peak >= 0 else math.nan,
                "vela_peak_x_um": float(sx[v_peak]) if v_peak >= 0 else math.nan,
                "vela_peak_y_um": float(sy[v_peak]) if v_peak >= 0 else math.nan,
                "avalanche_status": (
                    "near_floor"
                    if (not math.isfinite(s_threshold) or not math.isfinite(v_threshold)
                        or max(s_threshold, v_threshold) < avalanche_floor)
                    else "thresholded_peak"
                ),
                "sentaurus_p99": s_threshold,
                "vela_p99": v_threshold,
            })
        rows.append(row)
    return rows


def finite_or_blank(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float) and not math.isfinite(value):
        return ""
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "bias_V",
        "quantity",
        "status",
        "bv_current_log_error",
        "sentaurus_field",
        "metric",
        "field_error",
        "centerline_error",
        "sentaurus_peak_x_um",
        "sentaurus_peak_y_um",
        "vela_peak_x_um",
        "vela_peak_y_um",
        "avalanche_status",
        "sentaurus_p99",
        "vela_p99",
        "top_error_node_id",
        "top_error_x_um",
        "top_error_y_um",
        "top_error_value",
        "p_contact_error",
        "n_contact_error",
        "junction_error",
        "bulk_error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: finite_or_blank(row.get(key)) for key in fields})


def write_curve_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "bias_V",
        "sentaurus_current_A",
        "vela_current_A_per_um",
        "log10_abs_ratio",
        "abs_log10_error",
        "status",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: finite_or_blank(row.get(key)) for key in fields})


def clean_for_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_for_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_for_json(item) for item in value]
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def build_debug_ranking(field_rows: list[dict[str, Any]], curve_rows: list[dict[str, Any]]) -> dict[str, Any]:
    debug_priority = {
        "electron_density": 0,
        "electron_mobility": 1,
        "electric_field": 2,
        "hole_density": 3,
        "hole_mobility": 4,
        "potential": 5,
        "srh_recombination": 6,
        "avalanche_generation": 7,
    }
    ranked_fields = [
        row for row in field_rows
        if row.get("status") == "ok" and isinstance(row.get("field_error"), (int, float))
        and math.isfinite(float(row["field_error"]))
        and not (row.get("quantity") == "avalanche_generation" and row.get("avalanche_status") == "near_floor")
    ]
    ranked_fields.sort(
        key=lambda row: (
            debug_priority.get(str(row.get("quantity")), 99),
            -float(row["field_error"]),
        )
    )

    curve_errors = [
        float(row["abs_log10_error"]) for row in curve_rows
        if row.get("status") == "ok"
        and isinstance(row.get("abs_log10_error"), (int, float))
        and math.isfinite(float(row["abs_log10_error"]))
    ]
    return {
        "curve": {
            "max_abs_log10_error": max(curve_errors) if curve_errors else None,
            "current_floor_A": "see curve_compare.csv",
            "rows": curve_rows,
        },
        "field_rankings": ranked_fields[:20],
        "recommended_debug_order": [
            "electron_density",
            "electron_mobility",
            "electric_field",
            "avalanche_generation_thresholded",
        ],
    }


def write_readme(
    path: Path,
    biases: list[float],
    quantities: list[str],
    current_floor: float,
    avalanche_floor: float,
) -> None:
    text = f"""# PN2D BV Validation Summary

This directory is generated by `scripts/compare_pn2d_bv_multibias_fields.py`.

## Scope

- Biases: {", ".join(f"{bias:g} V" for bias in biases)}
- Quantities: {", ".join(quantities)}
- Curve current columns: Sentaurus `current_total` in A, Vela `current_total_A_per_um`.
- Decade current errors are disabled when either absolute current is below `{current_floor:g}`.
- Avalanche peak/rank diagnostics are marked near-floor when both p99 values are below `{avalanche_floor:g}`.

## Files

- `curve_compare.csv`: BV curve current ratio by bias.
- `field_compare.csv`: spatial field metrics with centerline, contact-local, junction-local, bulk, and top-error-node diagnostics.
- `debug_ranking.json`: ranked field discrepancies and recommended debug order.
- `pn2d_bv_multibias_field_compare.csv/json`: compatibility copies for existing consumers.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    biases = parse_biases(args.biases)
    quantities = parse_csv_list(args.quantities)
    unknown = [quantity for quantity in quantities if quantity not in FIELD_SPECS]
    if unknown:
        raise ValueError(f"unknown quantities: {unknown}; expected one of {sorted(FIELD_SPECS)}")

    sentaurus_exports = discover_sentaurus_exports(args.sentaurus_root, biases)
    vela_vtks = discover_vela_vtks(args.vela_vtk_root)
    curve_reference = load_curve_points(args.curve_reference)
    curve_candidate = load_curve_points(args.curve_candidate) if args.curve_candidate else []
    curve_rows = curve_compare_rows(biases, curve_reference, curve_candidate, args.current_floor)
    curve_error_by_bias = {
        bias_key(float(row["bias_V"])): row["abs_log10_error"]
        for row in curve_rows
        if row.get("status") == "ok"
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        bias_rows = compare_bias(
            bias,
            sentaurus_exports.get(bias_key(bias)),
            vela_vtks.get(bias_key(bias)),
            quantities,
            args.avalanche_floor,
        )
        curve_error = curve_error_by_bias.get(bias_key(bias), math.nan)
        for row in bias_rows:
            row["bv_current_log_error"] = curve_error
        rows.extend(bias_rows)

    write_curve_csv(args.out_dir / "curve_compare.csv", curve_rows)
    write_csv(args.out_dir / "field_compare.csv", rows)
    write_csv(args.out_dir / "pn2d_bv_multibias_field_compare.csv", rows)
    debug_ranking = build_debug_ranking(rows, curve_rows)
    (args.out_dir / "debug_ranking.json").write_text(
        json.dumps(clean_for_json(debug_ranking), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    trend_summary = build_bv_trend_summary(
        biases,
        curve_reference,
        curve_candidate,
        rows,
        args.current_floor,
    )
    (args.out_dir / "bv_trend_summary.json").write_text(
        json.dumps(clean_for_json(trend_summary), indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_readme(args.out_dir / "README.md", biases, quantities, args.current_floor, args.avalanche_floor)
    (args.out_dir / "pn2d_bv_multibias_field_compare.json").write_text(
        json.dumps(clean_for_json({"rows": rows, "curve_rows": curve_rows}), indent=2, sort_keys=True),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
