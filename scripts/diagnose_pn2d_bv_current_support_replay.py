#!/usr/bin/env python3
"""Replay PN2D avalanche source with alternate carrier-current supports.

The replay keeps Vela's impact-ionization coefficients fixed and changes only
the particle-flux support used by the current-density source term.  The named
variants are intentionally offline diagnostics; they are not solver modes.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import pn2d_fixed_state_manifest as fixed_manifest


Q_C = 1.602176634e-19

EDGE_FIELDS = [
    "bias_V",
    "edge_id",
    "node0",
    "node1",
    "sent_high_field_p95",
    "sent_high_h_current_p95",
    "sent_impact_m3_s",
    "edge_area_proxy_m2",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "sg_edge_e_flux_m2_s",
    "sg_edge_h_flux_m2_s",
    "sent_edge_e_flux_m2_s",
    "sent_edge_h_flux_m2_s",
    "sentaurus_node_e_flux_m2_s",
    "sentaurus_node_h_flux_m2_s",
    "charon_hcurl_e_flux_m2_s",
    "charon_hcurl_h_flux_m2_s",
    "genius_edge_gradqf_density_m3_s",
    "charon_hcurl_centroid_density_m3_s",
    "sentaurus_node_current_injected_density_m3_s",
    "sentaurus_edge_current_injected_density_m3_s",
    "genius_edge_gradqf_over_sentaurus",
    "charon_hcurl_centroid_over_sentaurus",
    "sentaurus_node_current_injected_over_sentaurus",
    "sentaurus_edge_current_injected_over_sentaurus",
    # Back-compatible aliases from the first replay script.
    "vela_baseline_density_m3_s",
    "sentaurus_hole_support_density_m3_s",
    "sentaurus_eh_support_density_m3_s",
    "genius_like_hole_support_density_m3_s",
    "charon_like_eh_support_density_m3_s",
    "vela_baseline_over_sentaurus",
    "sentaurus_hole_support_over_sentaurus",
    "sentaurus_eh_support_over_sentaurus",
    "genius_like_hole_support_over_sentaurus",
    "charon_like_eh_support_over_sentaurus",
]

VARIANT_DENSITY_KEYS = [
    "genius_edge_gradqf_density_m3_s",
    "charon_hcurl_centroid_density_m3_s",
    "sentaurus_node_current_injected_density_m3_s",
    "sentaurus_edge_current_injected_density_m3_s",
]

SUMMARY_FIELDS = [
    "bias_V",
    "edge_count",
    "sentaurus_source_integral",
    "genius_edge_gradqf_source_integral",
    "charon_hcurl_centroid_source_integral",
    "sentaurus_node_current_injected_source_integral",
    "sentaurus_edge_current_injected_source_integral",
    "genius_edge_gradqf_over_sentaurus",
    "charon_hcurl_centroid_over_sentaurus",
    "sentaurus_node_current_injected_over_sentaurus",
    "sentaurus_edge_current_injected_over_sentaurus",
    "charon_hcurl_over_genius_edge_gradqf",
]


def parse_args() -> argparse.Namespace:
    root = Path("build-release/reference_tcad/pn2d_sentaurus2018")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overlap-csv",
        type=Path,
        default=root / "reports/high_field_overlap_compare/high_field_overlap_edges.csv",
    )
    parser.add_argument(
        "--sg-edge-csv",
        type=Path,
        default=root / "reports/ionization_magnitude/sg_avalanche_edges.csv",
        help="Single sg_avalanche_edges CSV. Ignored when --sg-edge-manifest is set.",
    )
    parser.add_argument(
        "--sg-edge-manifest",
        type=Path,
        help="Fixed-state manifest whose entries contain sg_edges_csv paths.",
    )
    parser.add_argument(
        "--sentaurus-root",
        type=Path,
        default=root / "sentaurus_multibias_0p25v",
    )
    parser.add_argument(
        "--mesh",
        type=Path,
        default=Path("build/pn2d_bv_gradqf_full20/mesh.json"),
    )
    parser.add_argument("--bias-min", type=float, default=-20.0)
    parser.add_argument("--bias-max", type=float, default=-16.0)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def optional_float(raw: Any) -> float | None:
    if raw in (None, ""):
        return None
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def finite_float(raw: Any, default: float = 0.0) -> float:
    value = optional_float(raw)
    return value if value is not None else default


def bool_cell(raw: Any) -> bool:
    if isinstance(raw, bool):
        return raw
    return str(raw).strip().lower() in {"1", "true", "yes", "y"}


def ratio(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None or reference == 0.0:
        return None
    return candidate / reference


def in_bias_window(bias: float, lo: float, hi: float) -> bool:
    return min(lo, hi) - 1.0e-9 <= bias <= max(lo, hi) + 1.0e-9


def bias_key(bias: float) -> str:
    return f"{bias:.9f}"


def sentaurus_bias_dir(root: Path, bias: float) -> Path:
    text = f"{abs(bias):g}".replace(".", "p")
    return root / f"sentaurus_-{text}v"


def numeric_components(row: dict[str, str]) -> list[float]:
    values: list[float] = []
    for key, raw in row.items():
        if key == "node_id":
            continue
        value = optional_float(raw)
        if value is not None:
            values.append(value)
    return values


def load_field_magnitude(fields_dir: Path, name: str) -> dict[int, float]:
    result: dict[int, float] = {}
    for path in sorted(fields_dir.glob(f"{name}_region*.csv")):
        for row in read_csv(path):
            values = numeric_components(row)
            result[int(row["node_id"])] = math.sqrt(sum(v * v for v in values)) if values else 0.0
    return result


def current_to_particle_flux(current_a_cm2: float) -> float:
    return abs(current_a_cm2) * 1.0e4 / Q_C


def load_node_fluxes(sentaurus_root: Path, biases: Iterable[float]) -> dict[str, dict[str, dict[int, float]]]:
    result: dict[str, dict[str, dict[int, float]]] = {}
    for bias in biases:
        fields_dir = sentaurus_bias_dir(sentaurus_root, bias) / "fields"
        if not fields_dir.exists():
            continue
        e_current = load_field_magnitude(fields_dir, "eCurrentDensity")
        h_current = load_field_magnitude(fields_dir, "hCurrentDensity")
        result[bias_key(bias)] = {
            "electron": {node: current_to_particle_flux(value) for node, value in e_current.items()},
            "hole": {node: current_to_particle_flux(value) for node, value in h_current.items()},
        }
    return result


def load_overlap_rows(path: Path, bias_lo: float, bias_hi: float) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            bias = float(row["bias_V"])
            if in_bias_window(bias, bias_lo, bias_hi):
                rows.append(row)
    return rows


def sg_row_flux(row: dict[str, str], carrier: str) -> float:
    raw = finite_float(row.get(f"{carrier}_raw_flux_proxy"), math.nan)
    if math.isfinite(raw):
        return abs(raw)
    return abs(finite_float(row.get(f"{carrier}_flux_proxy")))


def sg_row_signed_flux(row: dict[str, str], carrier: str) -> float:
    signed = optional_float(row.get(f"{carrier}_raw_signed_flux_proxy"))
    if signed is not None:
        return signed
    # Older SG diagnostics only carried magnitudes. Use the row orientation as
    # positive so the HCurl replay still runs, while preserving the audit trail.
    return sg_row_flux(row, carrier)


def load_sg_edges_from_csv(path: Path, bias_lo: float, bias_hi: float) -> dict[tuple[str, int], dict[str, str]]:
    result: dict[tuple[str, int], dict[str, str]] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            bias = float(row["bias_V"])
            if in_bias_window(bias, bias_lo, bias_hi):
                result[(bias_key(bias), int(row["edge_id"]))] = row
    return result


def load_sg_edges_from_manifest(path: Path, bias_lo: float, bias_hi: float) -> dict[tuple[str, int], dict[str, str]]:
    manifest = json.loads(path.read_text(encoding="utf-8-sig"))
    result: dict[tuple[str, int], dict[str, str]] = {}
    for item in manifest:
        bias = float(item["bias_V"])
        if not in_bias_window(bias, bias_lo, bias_hi):
            continue
        fixed_manifest.require_materials_file(item, path)
        csv_path = Path(item["sg_edges_csv"])
        for row in read_csv(csv_path):
            result[(bias_key(float(row["bias_V"])), int(row["edge_id"]))] = row
    return result


def load_sg_edges(args: argparse.Namespace) -> dict[tuple[str, int], dict[str, str]]:
    if args.sg_edge_manifest:
        return load_sg_edges_from_manifest(args.sg_edge_manifest, args.bias_min, args.bias_max)
    return load_sg_edges_from_csv(args.sg_edge_csv, args.bias_min, args.bias_max)


def read_mesh(path: Path) -> tuple[dict[int, tuple[float, float]], list[list[int]], dict[tuple[int, int], int], dict[int, tuple[int, int]]]:
    data = json.loads(path.read_text(encoding="utf-8-sig"))
    # The Vela JSON mesh stores coordinates in microns for this reference case.
    nodes = {int(node["id"]): (float(node["x"]) * 1.0e-6, float(node["y"]) * 1.0e-6) for node in data["nodes"]}
    triangles = [[int(node_id) for node_id in tri["node_ids"]] for tri in data["triangles"]]
    edge_map: dict[tuple[int, int], int] = {}
    edge_nodes: dict[int, tuple[int, int]] = {}
    next_edge = 0
    for ids in triangles:
        for i in range(3):
            a = ids[i]
            b = ids[(i + 1) % 3]
            key = (a, b) if a <= b else (b, a)
            if key not in edge_map:
                edge_map[key] = next_edge
                edge_nodes[next_edge] = (a, b)
                next_edge += 1
    return nodes, triangles, edge_map, edge_nodes


def build_cell_maps(triangles: list[list[int]], edge_map: dict[tuple[int, int], int]) -> tuple[dict[int, list[int]], dict[int, list[int]]]:
    cell_edges: dict[int, list[int]] = {}
    edge_cells: dict[int, list[int]] = defaultdict(list)
    for cell_index, ids in enumerate(triangles):
        edges: list[int] = []
        for i in range(3):
            a = ids[i]
            b = ids[(i + 1) % 3]
            key = (a, b) if a <= b else (b, a)
            edge_id = edge_map[key]
            edges.append(edge_id)
            edge_cells[edge_id].append(cell_index)
        cell_edges[cell_index] = edges
    return cell_edges, dict(edge_cells)


def endpoint_average_flux(row: dict[str, str], node_flux: dict[int, float]) -> float | None:
    node0 = int(row["node0"])
    node1 = int(row["node1"])
    values = [node_flux[node] for node in (node0, node1) if node in node_flux]
    return statistics.fmean(values) if values else None


def tri_grad_lambdas(points: list[tuple[float, float]]) -> tuple[list[tuple[float, float]], float] | None:
    (x0, y0), (x1, y1), (x2, y2) = points
    area2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if abs(area2) <= 1.0e-40:
        return None
    grads = [
        ((y1 - y2) / area2, (x2 - x1) / area2),
        ((y2 - y0) / area2, (x0 - x2) / area2),
        ((y0 - y1) / area2, (x1 - x0) / area2),
    ]
    return grads, abs(area2) * 0.5


def local_signed_flux(row: dict[str, str], carrier: str, from_node: int, to_node: int) -> float:
    value = sg_row_signed_flux(row, carrier)
    row_from = int(row["node0"])
    row_to = int(row["node1"])
    if row_from == from_node and row_to == to_node:
        return value
    if row_from == to_node and row_to == from_node:
        return -value
    # This should not happen for a consistent mesh; fall back to the magnitude
    # so one malformed edge does not erase the whole replay.
    return abs(value)


def charon_cell_centroid_magnitude(
    cell_id: int,
    carrier: str,
    bias_sg_edges: dict[int, dict[str, str]],
    nodes: dict[int, tuple[float, float]],
    triangles: list[list[int]],
    cell_edges: dict[int, list[int]],
) -> float | None:
    ids = triangles[cell_id]
    points = [nodes[node_id] for node_id in ids]
    grad_data = tri_grad_lambdas(points)
    if grad_data is None:
        return None
    grads, _area = grad_data
    barycentric_points = [
        (7.0 / 12.0, 5.0 / 24.0, 5.0 / 24.0),
        (5.0 / 24.0, 7.0 / 12.0, 5.0 / 24.0),
        (5.0 / 24.0, 5.0 / 24.0, 7.0 / 12.0),
    ]
    edge_pairs = [(0, 1), (1, 2), (2, 0)]
    magnitudes: list[float] = []
    for lambdas in barycentric_points:
        jx = 0.0
        jy = 0.0
        used = 0
        for local_edge_index, (i, j) in enumerate(edge_pairs):
            edge_id = cell_edges[cell_id][local_edge_index]
            row = bias_sg_edges.get(edge_id)
            if row is None:
                continue
            from_node = ids[i]
            to_node = ids[j]
            x0, y0 = nodes[from_node]
            x1, y1 = nodes[to_node]
            length = math.hypot(x1 - x0, y1 - y0)
            if length <= 1.0e-30:
                continue
            dof = local_signed_flux(row, carrier, from_node, to_node) * length
            gx = lambdas[i] * grads[j][0] - lambdas[j] * grads[i][0]
            gy = lambdas[i] * grads[j][1] - lambdas[j] * grads[i][1]
            jx += dof * gx
            jy += dof * gy
            used += 1
        if used:
            magnitudes.append(math.hypot(jx, jy))
    return statistics.fmean(magnitudes) if magnitudes else None


def charon_hcurl_edge_fluxes(
    sg_edges: dict[tuple[str, int], dict[str, str]],
    nodes: dict[int, tuple[float, float]],
    triangles: list[list[int]],
    cell_edges: dict[int, list[int]],
    edge_cells: dict[int, list[int]],
    carrier: str,
) -> dict[tuple[str, int], float]:
    by_bias: dict[str, dict[int, dict[str, str]]] = defaultdict(dict)
    for (bkey, edge_id), row in sg_edges.items():
        by_bias[bkey][edge_id] = row

    result: dict[tuple[str, int], float] = {}
    for bkey, bias_edges in by_bias.items():
        cell_magnitudes: dict[int, float] = {}
        for cell_id in range(len(triangles)):
            value = charon_cell_centroid_magnitude(
                cell_id, carrier, bias_edges, nodes, triangles, cell_edges)
            if value is not None:
                cell_magnitudes[cell_id] = value
        for edge_id, cells in edge_cells.items():
            values = [cell_magnitudes[cell_id] for cell_id in cells if cell_id in cell_magnitudes]
            if values:
                result[(bkey, edge_id)] = statistics.fmean(values)
            elif (bkey, edge_id) in sg_edges:
                result[(bkey, edge_id)] = sg_row_flux(sg_edges[(bkey, edge_id)], carrier)
    return result


def finite_values(rows: list[dict[str, Any]], key: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = optional_float(row.get(key))
        if value is not None:
            values.append(value)
    return values


def mean(values: list[float]) -> float | None:
    return statistics.fmean(values) if values else None


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    pos = int(round((len(ordered) - 1) * q))
    return ordered[pos]


def clean_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def make_replay_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    overlap_rows = load_overlap_rows(args.overlap_csv, args.bias_min, args.bias_max)
    biases = sorted({float(row["bias_V"]) for row in overlap_rows})
    sg_edges = load_sg_edges(args)
    node_fluxes = load_node_fluxes(args.sentaurus_root, biases)
    nodes, triangles, edge_map, _edge_nodes = read_mesh(args.mesh)
    cell_edges, edge_cells = build_cell_maps(triangles, edge_map)
    charon_e_fluxes = charon_hcurl_edge_fluxes(sg_edges, nodes, triangles, cell_edges, edge_cells, "electron")
    charon_h_fluxes = charon_hcurl_edge_fluxes(sg_edges, nodes, triangles, cell_edges, edge_cells, "hole")

    rows: list[dict[str, Any]] = []
    for row in overlap_rows:
        bias = float(row["bias_V"])
        bkey = bias_key(bias)
        edge_id = int(row["edge_id"])
        sg = sg_edges.get((bkey, edge_id), {})

        edge_area = finite_float(row.get("edge_area_proxy_m2"))
        sent_impact = finite_float(row.get("sent_impact_m3_s"))
        alpha_e = finite_float(sg.get("electron_alpha_m_inv"))
        alpha_h = finite_float(sg.get("hole_alpha_m_inv"))
        sg_e_flux = sg_row_flux(sg, "electron") if sg else abs(finite_float(row.get("vela_e_flux_proxy_m2_s")))
        sg_h_flux = sg_row_flux(sg, "hole") if sg else abs(finite_float(row.get("vela_h_flux_proxy_m2_s")))
        sent_edge_e_flux = abs(finite_float(row.get("sent_e_particle_flux_m2_s")))
        sent_edge_h_flux = abs(finite_float(row.get("sent_h_particle_flux_m2_s")))

        node_flux = node_fluxes.get(bkey, {})
        sent_node_e_flux = endpoint_average_flux(row, node_flux.get("electron", {}))
        sent_node_h_flux = endpoint_average_flux(row, node_flux.get("hole", {}))
        charon_e_flux = charon_e_fluxes.get((bkey, edge_id))
        charon_h_flux = charon_h_fluxes.get((bkey, edge_id))

        genius_edge_gradqf = alpha_e * sg_e_flux + alpha_h * sg_h_flux
        charon_hcurl = alpha_e * (charon_e_flux or 0.0) + alpha_h * (charon_h_flux or 0.0)
        sentaurus_node = alpha_e * (sent_node_e_flux or 0.0) + alpha_h * (sent_node_h_flux or 0.0)
        sentaurus_edge = alpha_e * sent_edge_e_flux + alpha_h * sent_edge_h_flux
        sentaurus_hole = alpha_e * sg_e_flux + alpha_h * sent_edge_h_flux

        item: dict[str, Any] = {
            "bias_V": bias,
            "edge_id": edge_id,
            "node0": int(row["node0"]),
            "node1": int(row["node1"]),
            "sent_high_field_p95": bool_cell(row.get("sent_high_field_p95")),
            "sent_high_h_current_p95": bool_cell(row.get("sent_high_h_current_p95")),
            "sent_impact_m3_s": sent_impact,
            "edge_area_proxy_m2": edge_area,
            "electron_alpha_m_inv": alpha_e,
            "hole_alpha_m_inv": alpha_h,
            "sg_edge_e_flux_m2_s": sg_e_flux,
            "sg_edge_h_flux_m2_s": sg_h_flux,
            "sent_edge_e_flux_m2_s": sent_edge_e_flux,
            "sent_edge_h_flux_m2_s": sent_edge_h_flux,
            "sentaurus_node_e_flux_m2_s": sent_node_e_flux,
            "sentaurus_node_h_flux_m2_s": sent_node_h_flux,
            "charon_hcurl_e_flux_m2_s": charon_e_flux,
            "charon_hcurl_h_flux_m2_s": charon_h_flux,
            "genius_edge_gradqf_density_m3_s": genius_edge_gradqf,
            "charon_hcurl_centroid_density_m3_s": charon_hcurl,
            "sentaurus_node_current_injected_density_m3_s": sentaurus_node,
            "sentaurus_edge_current_injected_density_m3_s": sentaurus_edge,
            "vela_baseline_density_m3_s": genius_edge_gradqf,
            "sentaurus_hole_support_density_m3_s": sentaurus_hole,
            "sentaurus_eh_support_density_m3_s": sentaurus_edge,
            "genius_like_hole_support_density_m3_s": genius_edge_gradqf,
            "charon_like_eh_support_density_m3_s": charon_hcurl,
        }
        for key in VARIANT_DENSITY_KEYS:
            item[key.replace("_density_m3_s", "_over_sentaurus")] = ratio(item[key], sent_impact)
        item["vela_baseline_over_sentaurus"] = ratio(genius_edge_gradqf, sent_impact)
        item["sentaurus_hole_support_over_sentaurus"] = ratio(sentaurus_hole, sent_impact)
        item["sentaurus_eh_support_over_sentaurus"] = ratio(sentaurus_edge, sent_impact)
        item["genius_like_hole_support_over_sentaurus"] = ratio(genius_edge_gradqf, sent_impact)
        item["charon_like_eh_support_over_sentaurus"] = ratio(charon_hcurl, sent_impact)
        rows.append(item)
    return rows


def summarize_subset(rows: list[dict[str, Any]]) -> dict[str, Any]:
    item: dict[str, Any] = {"count": len(rows)}
    sent_integral = sum(float(row["sent_impact_m3_s"]) * float(row["edge_area_proxy_m2"]) for row in rows)
    item["sentaurus_source_integral"] = sent_integral
    for density_key in VARIANT_DENSITY_KEYS:
        ratio_key = density_key.replace("_density_m3_s", "_over_sentaurus")
        ratios = finite_values(rows, ratio_key)
        variant_integral = sum(float(row[density_key]) * float(row["edge_area_proxy_m2"]) for row in rows)
        prefix = density_key.replace("_density_m3_s", "")
        item[f"{prefix}_source_integral"] = variant_integral
        item[f"{prefix}_ratio_median"] = median(ratios)
        item[f"{prefix}_ratio_mean"] = mean(ratios)
        item[f"{prefix}_ratio_p10"] = quantile(ratios, 0.10)
        item[f"{prefix}_ratio_p90"] = quantile(ratios, 0.90)
        item[f"{prefix}_area_weighted_ratio"] = ratio(variant_integral, sent_integral)
    item["charon_hcurl_over_genius_edge_gradqf"] = ratio(
        item.get("charon_hcurl_centroid_source_integral"),
        item.get("genius_edge_gradqf_source_integral"),
    )
    return item


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_bias: dict[str, Any] = {}
    for bias in sorted({float(row["bias_V"]) for row in rows}):
        subset = [row for row in rows if abs(float(row["bias_V"]) - bias) <= 1.0e-9]
        high_field = [row for row in subset if row["sent_high_field_p95"]]
        high_h = [row for row in subset if row["sent_high_h_current_p95"]]
        overlap = [row for row in high_field if row["sent_high_h_current_p95"]]
        by_bias[f"{bias:g}"] = {
            "all_edges": summarize_subset(subset),
            "sentaurus_p95_high_field": summarize_subset(high_field),
            "sentaurus_p95_high_h_current": summarize_subset(high_h),
            "p95_high_field_and_high_h_current_count": len(overlap),
        }
    return {
        "schema": "pn2d.current_support_replay.v2",
        "bias_count": len(by_bias),
        "edge_row_count": len(rows),
        "variants": {
            "genius_edge_gradqf": "Vela alpha with Vela SG edge electron+hole particle flux support; Genius-like GradQF source support before endpoint directional distribution.",
            "charon_hcurl_centroid": "Vela alpha with signed Vela SG edge flux reconstructed through Tri3 HCurl/Nedelec basis at sub-CV centroid proxies, then edge-averaged.",
            "sentaurus_node_current_injected": "Vela alpha with endpoint-averaged Sentaurus node |Jn|/q and |Jp|/q support.",
            "sentaurus_edge_current_injected": "Vela alpha with Sentaurus edge |Jn|/q and |Jp|/q support from the overlap report.",
        },
        "by_bias": by_bias,
    }


def summary_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias_text, bias_payload in sorted(payload["by_bias"].items(), key=lambda item: float(item[0])):
        all_edges = bias_payload["all_edges"]
        row = {
            "bias_V": float(bias_text),
            "edge_count": all_edges["count"],
            "sentaurus_source_integral": all_edges["sentaurus_source_integral"],
            "genius_edge_gradqf_source_integral": all_edges["genius_edge_gradqf_source_integral"],
            "charon_hcurl_centroid_source_integral": all_edges["charon_hcurl_centroid_source_integral"],
            "sentaurus_node_current_injected_source_integral": all_edges["sentaurus_node_current_injected_source_integral"],
            "sentaurus_edge_current_injected_source_integral": all_edges["sentaurus_edge_current_injected_source_integral"],
            "genius_edge_gradqf_over_sentaurus": all_edges["genius_edge_gradqf_area_weighted_ratio"],
            "charon_hcurl_centroid_over_sentaurus": all_edges["charon_hcurl_centroid_area_weighted_ratio"],
            "sentaurus_node_current_injected_over_sentaurus": all_edges["sentaurus_node_current_injected_area_weighted_ratio"],
            "sentaurus_edge_current_injected_over_sentaurus": all_edges["sentaurus_edge_current_injected_area_weighted_ratio"],
            "charon_hcurl_over_genius_edge_gradqf": all_edges["charon_hcurl_over_genius_edge_gradqf"],
        }
        rows.append(row)
    return rows


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    rows = make_replay_rows(args)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "current_support_replay_edges.csv", EDGE_FIELDS, rows)
    payload = summarize(rows)
    write_csv(args.out_dir / "current_support_replay_summary.csv", SUMMARY_FIELDS, summary_rows(payload))
    (args.out_dir / "current_support_replay_summary.json").write_text(
        json.dumps(clean_json(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(clean_json({
        "bias_count": payload["bias_count"],
        "edge_row_count": payload["edge_row_count"],
        "out_dir": str(args.out_dir),
    }), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
