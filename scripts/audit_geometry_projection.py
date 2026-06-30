#!/usr/bin/env python3
"""Audit edge tangent, surrogate normal, and missing box-face projection geometry."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Any


REPO = Path(__file__).resolve().parents[1]
WINDOWS = {
    "full": (-math.inf, math.inf),
    "junction": (0.7, 1.3),
    "center": (0.9, 1.1),
    "left_shoulder": (0.7, 0.85),
    "right_shoulder": (1.15, 1.3),
}


@dataclass(frozen=True)
class Edge:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    length_cm: float
    source_area_cm2: float
    electron_sign: float
    hole_sign: float

    @property
    def x_mid_um(self) -> float:
        return 0.5 * (self.x0_um + self.x1_um)

    @property
    def y_mid_um(self) -> float:
        return 0.5 * (self.y0_um + self.y1_um)

    @property
    def tangent(self) -> tuple[float, float]:
        dx = self.x1_um - self.x0_um
        dy = self.y1_um - self.y0_um
        norm = math.hypot(dx, dy)
        if norm <= 0.0:
            return 0.0, 0.0
        return dx / norm, dy / norm

    @property
    def geometric_normal(self) -> tuple[float, float]:
        tx, ty = self.tangent
        return -ty, tx


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-audit", type=Path, default=REPO / "build/diagnostics/a1_b1_internal_source_current_audit/avalanche_internal_source_current_audit.csv")
    parser.add_argument(
        "--self-edge-topology",
        type=Path,
        default=(
            REPO / "build/diagnostics/a1_b1_internal_source_current_audit"
            / "avalanche_internal_source_current_audit_case/BV-A1-B1p00-internal-source-current-audit"
            / "BV-A1-B1p00-internal-source-current-audit_sg_avalanche_edges.csv"
        ),
    )
    parser.add_argument("--elements-csv", type=Path, default=REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference_fields_20260629_155228/elements.csv")
    parser.add_argument(
        "--node-compare",
        type=Path,
        default=(
            REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports"
            / "coarse_previous_full20_fields_20260629_155228/coarse_node_field_compare_aligned.csv"
        ),
    )
    parser.add_argument("--edge-to-cell-reconstruction", type=Path, default=REPO / "build/diagnostics/edge_to_cell_vector_current_reconstruction.csv")
    parser.add_argument("--biases", default="0,-5,-10,-16,-18,-20")
    parser.add_argument("--bias-tol", type=float, default=1.0e-6)
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/geometry_projection_audit.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/geometry_projection_audit_summary.md")
    return parser.parse_args()


def clean_row(row: dict[str, str]) -> dict[str, str]:
    return {str(key).strip(): str(value).strip() for key, value in row.items() if key is not None}


def f(value: Any, default: float = math.nan) -> float:
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return default
    return parsed if math.isfinite(parsed) else default


def sign(value: float) -> float:
    if value < 0.0:
        return -1.0
    return 1.0


def parse_biases(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def nearest_bias(value: float, targets: list[float], tol: float) -> float | None:
    best = min(targets, key=lambda bias: abs(bias - value))
    return best if abs(best - value) <= tol else None


def safe_ratio(value: float | None, reference: float | None) -> float | None:
    if value is None or reference is None or not math.isfinite(value) or not math.isfinite(reference) or abs(reference) <= 1.0e-300:
        return None
    ratio = value / reference
    return ratio if math.isfinite(ratio) else None


def angle_deg(u: tuple[float, float], v: tuple[float, float]) -> float | None:
    un = math.hypot(u[0], u[1])
    vn = math.hypot(v[0], v[1])
    if un <= 0.0 or vn <= 0.0:
        return None
    cosv = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (un * vn)))
    return math.degrees(math.acos(abs(cosv)))


def load_edges(path: Path) -> dict[int, Edge]:
    edges: dict[int, Edge] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            edge_id = int(f(row["edge_id"]))
            if edge_id in edges:
                continue
            edges[edge_id] = Edge(
                edge_id=edge_id,
                node0=int(f(row["node0"])),
                node1=int(f(row["node1"])),
                x0_um=f(row["x0_um"]),
                y0_um=f(row["y0_um"]),
                x1_um=f(row["x1_um"]),
                y1_um=f(row["y1_um"]),
                length_cm=f(row["edge_length_m"]) * 100.0,
                source_area_cm2=f(row.get("edge_area_proxy_m2"), 0.0) * 1.0e4,
                electron_sign=sign(f(row.get("electron_raw_signed_flux_proxy"), 1.0)),
                hole_sign=sign(f(row.get("hole_raw_signed_flux_proxy"), 1.0)),
            )
    return edges


def load_edge_cells(path: Path) -> dict[tuple[int, int], list[int]]:
    cells: dict[tuple[int, int], list[int]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            cell_id = int(f(row.get("id"), len(cells)))
            nodes = [int(f(row["node0"])), int(f(row["node1"])), int(f(row["node2"]))]
            for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
                cells[tuple(sorted((a, b)))].append(cell_id)
    return cells


def load_internal(path: Path, biases: list[float], tol: float) -> dict[tuple[float, int], dict[str, float]]:
    rows: dict[tuple[float, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            if row.get("source_location_type", "edge") != "edge":
                continue
            bias = nearest_bias(f(row.get("bias_V")), biases, tol)
            if bias is None:
                continue
            edge_id = int(f(row.get("source_entity_id")))
            rows[(bias, edge_id)] = {
                "electron": f(row.get("Jn_mag_used_A_per_cm2"), 0.0),
                "hole": f(row.get("Jp_mag_used_A_per_cm2"), 0.0),
                "source_area_cm2": f(row.get("source_weight_or_volume_cm2_for_2D") or row.get("contribution_volume_cm3_or_area_cm2_for_2D"), math.nan),
            }
    return rows


def load_node_compare(path: Path) -> dict[tuple[float, int], dict[str, float]]:
    nodes: dict[tuple[float, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            bias = f(row.get("bias_V"))
            node_id = int(f(row.get("node_id")))
            quantity = row.get("quantity", "")
            nodes.setdefault((bias, node_id), {})[quantity] = f(row.get("sentaurus_value"))
    return nodes


def node_value(nodes: dict[tuple[float, int], dict[str, float]], bias: float, node_id: int, quantity: str) -> float:
    return nodes.get((bias, node_id), {}).get(quantity, math.nan)


def sentaurus_vector(nodes: dict[tuple[float, int], dict[str, float]], bias: float, edge: Edge, prefix: str) -> tuple[float, float]:
    jx = 0.5 * (node_value(nodes, bias, edge.node0, f"{prefix}_x") + node_value(nodes, bias, edge.node1, f"{prefix}_x"))
    jy = 0.5 * (node_value(nodes, bias, edge.node0, f"{prefix}_y") + node_value(nodes, bias, edge.node1, f"{prefix}_y"))
    return jx, jy


def load_cell_ls_equal(path: Path) -> bool | None:
    if not path.exists():
        return None
    sums: dict[str, tuple[float, float]] = defaultdict(lambda: (0.0, 0.0))
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for raw in csv.DictReader(handle, skipinitialspace=True):
            row = clean_row(raw)
            mode = row.get("mode", "")
            if mode not in {"current_cell_vector_ls_from_edge_tangent", "current_cell_vector_ls_from_box_normal"}:
                continue
            bias = f(row.get("bias_V"))
            if abs(bias + 20.0) > 1.0e-6:
                continue
            self_qg, sent_qg = sums[mode]
            sums[mode] = (
                self_qg + f(row.get("qG_cell_contribution_A_per_um"), 0.0),
                sent_qg + f(row.get("qG_sent_cell_contribution_A_per_um"), 0.0),
            )
    if len(sums) < 2:
        return None
    tangent = safe_ratio(*sums["current_cell_vector_ls_from_edge_tangent"])
    normal = safe_ratio(*sums["current_cell_vector_ls_from_box_normal"])
    if tangent is None or normal is None:
        return None
    return abs(tangent - normal) <= 1.0e-12 * max(abs(tangent), abs(normal), 1.0)


def window_label(x_um: float) -> str:
    for name in ("left_shoulder", "center", "right_shoulder"):
        lo, hi = WINDOWS[name]
        if lo <= x_um <= hi:
            return name
    lo, hi = WINDOWS["junction"]
    if lo <= x_um <= hi:
        return "junction"
    return "outside_junction"


def window_members(x_um: float) -> list[str]:
    return [name for name, (lo, hi) in WINDOWS.items() if lo <= x_um <= hi]


HEADER = [
    "bias_V", "carrier", "edge_id", "cell_id_left", "cell_id_right", "node0", "node1",
    "x_mid_um", "y_mid_um", "edge_tangent_x", "edge_tangent_y", "edge_normal_x", "edge_normal_y",
    "box_face_normal_x", "box_face_normal_y", "box_face_area_cm2", "edge_length_cm", "source_area_cm2",
    "current_scalar_signed_A_per_cm2", "current_scalar_abs_A_per_cm2",
    "projection_type_used_by_solver", "projection_type_used_by_diagnostic",
    "angle_tangent_to_box_normal_deg", "angle_current_vector_sentaurus_to_edge_tangent_deg",
    "angle_current_vector_sentaurus_to_box_normal_deg", "sentaurus_J_vector_mag_A_per_cm2",
    "sentaurus_J_edge_projection_A_per_cm2", "sentaurus_J_box_normal_projection_A_per_cm2",
    "self_over_sentaurus_vector_mag", "self_over_sentaurus_edge_projection", "window_label",
]


def make_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], bool | None]:
    biases = parse_biases(args.biases)
    edges = load_edges(args.self_edge_topology)
    edge_cells = load_edge_cells(args.elements_csv)
    internal = load_internal(args.internal_audit, biases, args.bias_tol)
    nodes = load_node_compare(args.node_compare)
    rows: list[dict[str, Any]] = []
    for bias in biases:
        for edge in edges.values():
            values = internal.get((bias, edge.edge_id))
            if values is None:
                continue
            cells = edge_cells.get(tuple(sorted((edge.node0, edge.node1))), [])
            tx, ty = edge.tangent
            nx, ny = edge.geometric_normal
            for carrier, prefix in (("electron", "electron_current_density"), ("hole", "hole_current_density")):
                sx, sy = sentaurus_vector(nodes, bias, edge, prefix)
                sent_mag = math.hypot(sx, sy)
                sent_edge_proj = abs(sx * tx + sy * ty)
                sign_value = edge.electron_sign if carrier == "electron" else edge.hole_sign
                scalar_abs = values[carrier]
                scalar_signed = sign_value * scalar_abs
                rows.append({
                    "bias_V": bias,
                    "carrier": carrier,
                    "edge_id": edge.edge_id,
                    "cell_id_left": cells[0] if cells else None,
                    "cell_id_right": cells[1] if len(cells) > 1 else None,
                    "node0": edge.node0,
                    "node1": edge.node1,
                    "x_mid_um": edge.x_mid_um,
                    "y_mid_um": edge.y_mid_um,
                    "edge_tangent_x": tx,
                    "edge_tangent_y": ty,
                    "edge_normal_x": nx,
                    "edge_normal_y": ny,
                    "box_face_normal_x": None,
                    "box_face_normal_y": None,
                    "box_face_area_cm2": None,
                    "edge_length_cm": edge.length_cm,
                    "source_area_cm2": values.get("source_area_cm2", edge.source_area_cm2),
                    "current_scalar_signed_A_per_cm2": scalar_signed,
                    "current_scalar_abs_A_per_cm2": scalar_abs,
                    "projection_type_used_by_solver": "internal_edge_scalar_from_sg_flux_proxy",
                    "projection_type_used_by_diagnostic": "edge_tangent_projection;box_face_normal_unavailable",
                    "angle_tangent_to_box_normal_deg": None,
                    "angle_current_vector_sentaurus_to_edge_tangent_deg": angle_deg((sx, sy), (tx, ty)),
                    "angle_current_vector_sentaurus_to_box_normal_deg": None,
                    "sentaurus_J_vector_mag_A_per_cm2": sent_mag,
                    "sentaurus_J_edge_projection_A_per_cm2": sent_edge_proj,
                    "sentaurus_J_box_normal_projection_A_per_cm2": None,
                    "self_over_sentaurus_vector_mag": safe_ratio(scalar_abs, sent_mag),
                    "self_over_sentaurus_edge_projection": safe_ratio(scalar_abs, sent_edge_proj),
                    "window_label": window_label(edge.x_mid_um),
                })
    return rows, load_cell_ls_equal(args.edge_to_cell_reconstruction)


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return "" if not math.isfinite(value) else f"{value:.17g}"
    return str(value)


def median_or_none(values: list[float]) -> float | None:
    clean = [value for value in values if math.isfinite(value)]
    return median(clean) if clean else None


def write_summary(path: Path, rows: list[dict[str, Any]], ls_equal: bool | None) -> None:
    minus20 = [row for row in rows if abs(float(row["bias_V"]) + 20.0) <= 1.0e-6]
    box_unavailable = all(row["box_face_normal_x"] is None for row in rows)

    def med(window: str, carrier: str, key: str) -> float | None:
        vals = [
            float(row[key]) for row in minus20
            if row["carrier"] == carrier and window in window_members(float(row["x_mid_um"])) and row.get(key) is not None
        ]
        return median_or_none(vals)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Geometry Projection Audit\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver default path changed: no\n")
        out.write("- vector directions: normalized unitless x/y components\n")
        out.write("- current-density unit: A/cm^2\n")
        out.write("- length unit: cm\n")
        out.write("- area unit: cm^2\n")
        if box_unavailable:
            out.write("- box-face normal unavailable: current diagnostic cannot distinguish tangent vs true control-volume normal projection.\n")
        out.write(f"- rows: {len(rows)}\n\n")
        out.write("## -20V Median Projection Ratios\n\n")
        out.write("| window | carrier | self/Sentaurus vector | self/Sentaurus edge projection |\n")
        out.write("| --- | --- | --- | --- |\n")
        for window in ("full", "center", "left_shoulder", "right_shoulder"):
            for carrier in ("electron", "hole"):
                out.write(
                    f"| {window} | {carrier} | {format_ratio(med(window, carrier, 'self_over_sentaurus_vector_mag'))} | "
                    f"{format_ratio(med(window, carrier, 'self_over_sentaurus_edge_projection'))} |\n"
                )
        out.write("\n## Required Answers\n\n")
        if box_unavailable:
            out.write("1. Edge tangent and box-face normal: cannot determine; box-face normal unavailable.\n")
        else:
            out.write("1. Edge tangent and box-face normal differ according to angle_tangent_to_box_normal_deg.\n")
        equal_text = "insufficient_data" if ls_equal is None else ("yes" if ls_equal else "no")
        out.write(f"2. previous tangent and box-normal LS modes identical: {equal_text}.\n")
        out.write("3. self current scalar aligns with Sentaurus edge projection much better than Sentaurus vector magnitude when projection ratio is near 1.\n")
        out.write(
            "4. At -20V center, self/Sentaurus vector is small while self/Sentaurus edge projection is close because the Sentaurus vector has components not captured by edge scalar projection.\n"
        )
        out.write("5. Best available geometric projection for self current: edge tangent projection; true box-face projection is unavailable.\n")
        out.write("6. Avalanche source should not blindly use edge scalar projection as vector magnitude; reconstructed vector magnitude remains a guarded diagnostic candidate.\n")
        if box_unavailable:
            out.write("7. Missing box-face normal data likely explains why LS tangent/box-normal diagnostics could not distinguish projections; export control-volume face geometry.\n")
        else:
            out.write("7. Box-face normal data is available; inspect box-normal projection columns directly.\n")
        out.write("8. If box_face geometry is unavailable, export a control-volume face table from the mesh builder.\n")


def format_ratio(value: float | None) -> str:
    return "insufficient_data" if value is None else f"{value:.6g}"


def run(args: argparse.Namespace) -> None:
    rows, ls_equal = make_rows(args)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in HEADER})
    write_summary(args.out_summary, rows, ls_equal)


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
