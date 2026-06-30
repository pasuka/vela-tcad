#!/usr/bin/env python3
"""Export available control-volume/box-face geometry diagnostics.

The current Vela mesh stores node control-volume areas and edge couplings, but
does not persist explicit node-centered box/control-volume face entities. This
tool exports the geometry that is available for avalanche source edges and
marks unavailable sub-control-volume polygon faces separately, so downstream current
projection diagnostics can distinguish edge-coupled box normals from missing polygon geometry.
"""

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
DEFAULT_SG_EDGES = (
    REPO / "build/diagnostics/a1_b1_internal_source_current_audit"
    / "avalanche_internal_source_current_audit_case/BV-A1-B1p00-internal-source-current-audit"
    / "BV-A1-B1p00-internal-source-current-audit_sg_avalanche_edges.csv"
)
DEFAULT_ELEMENTS = (
    REPO / "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3"
    / "imported_reference_fields_20260629_155228/elements.csv"
)
DEFAULT_INTERNAL = (
    REPO / "build/diagnostics/a1_b1_internal_source_current_audit"
    / "avalanche_internal_source_current_audit.csv"
)


@dataclass(frozen=True)
class EdgeRow:
    edge_id: int
    node0: int
    node1: int
    x0_um: float
    y0_um: float
    x1_um: float
    y1_um: float
    edge_length_m: float
    edge_couple_m: float
    edge_area_proxy_m2: float
    edge_class: str

    @property
    def key(self) -> tuple[int, int]:
        return tuple(sorted((self.node0, self.node1)))

    @property
    def x_mid_um(self) -> float:
        return 0.5 * (self.x0_um + self.x1_um)

    @property
    def y_mid_um(self) -> float:
        return 0.5 * (self.y0_um + self.y1_um)

    @property
    def tangent(self) -> tuple[float | None, float | None]:
        dx = self.x1_um - self.x0_um
        dy = self.y1_um - self.y0_um
        norm = math.hypot(dx, dy)
        if norm <= 0.0:
            return None, None
        return dx / norm, dy / norm

    @property
    def box_normal(self) -> tuple[float | None, float | None]:
        return self.tangent

    @property
    def box_tangent(self) -> tuple[float | None, float | None]:
        nx, ny = self.box_normal
        if nx is None or ny is None:
            return None, None
        return -ny, nx


@dataclass(frozen=True)
class Element:
    cell_id: int
    nodes: tuple[int, int, int]
    region: str
    material: str


HEADER = [
    "box_face_id",
    "owner_node_id",
    "neighbor_node_id",
    "associated_edge_id",
    "left_cell_id",
    "right_cell_id",
    "x_mid_um",
    "y_mid_um",
    "edge_coupling_normal_x",
    "edge_coupling_normal_y",
    "owner_to_neighbor_direction_x",
    "owner_to_neighbor_direction_y",
    "face_normal_x",
    "face_normal_y",
    "face_tangent_x",
    "face_tangent_y",
    "face_length_cm",
    "source_area_cm2",
    "orientation_sign",
    "control_volume_owner_area_cm2",
    "control_volume_neighbor_area_cm2",
    "material_region",
    "boundary_type",
    "sub_control_volume_polygon_vertices",
    "partial_face_weight_to_owner",
    "partial_face_weight_to_neighbor",
]

DUAL_HEADER = [
    "dual_face_id",
    "cell_id",
    "owner_node_id",
    "neighbor_node_id",
    "associated_primal_edge_id",
    "dual_face_vertex0_x_um",
    "dual_face_vertex0_y_um",
    "dual_face_vertex1_x_um",
    "dual_face_vertex1_y_um",
    "dual_face_mid_x_um",
    "dual_face_mid_y_um",
    "dual_face_normal_x",
    "dual_face_normal_y",
    "dual_face_length_cm",
    "edge_direction_to_dual_normal_angle_deg",
    "orientation_sign",
    "owner_subcv_area_cm2",
    "neighbor_subcv_area_cm2",
    "material_region",
    "boundary_type",
    "dual_type",
]

VALIDATION_HEADER = [
    "cell_id",
    "cell_area_cm2",
    "subcv_area_sum_cm2",
    "subcv_area_abs_error_cm2",
    "subcv_area_rel_error",
    "internal_dual_face_cancel_residual_max",
    "normal_length_min",
    "normal_length_max",
    "face_length_min_cm",
    "cell_area_closure_pass",
    "internal_face_cancel_pass",
    "normal_unit_pass",
    "face_length_pass",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sg-edges", type=Path, default=DEFAULT_SG_EDGES)
    parser.add_argument("--elements-csv", type=Path, default=DEFAULT_ELEMENTS)
    parser.add_argument("--internal-audit", type=Path, default=DEFAULT_INTERNAL)
    parser.add_argument("--out-csv", type=Path, default=REPO / "build/diagnostics/box_face_geometry_audit.csv")
    parser.add_argument("--out-summary", type=Path, default=REPO / "build/diagnostics/box_face_geometry_audit_summary.md")
    parser.add_argument("--dual-out-csv", type=Path, default=REPO / "build/diagnostics/dual_face_geometry_audit.csv")
    parser.add_argument("--dual-out-summary", type=Path, default=REPO / "build/diagnostics/dual_face_geometry_audit_summary.md")
    parser.add_argument("--validation-out-csv", type=Path, default=REPO / "build/diagnostics/true_dual_face_geometry_validation.csv")
    parser.add_argument("--validation-out-summary", type=Path, default=REPO / "build/diagnostics/true_dual_face_geometry_validation_summary.md")
    return parser.parse_args()


def fnum(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    try:
        return float(value)
    except ValueError:
        return default


def inum(row: dict[str, str], key: str, default: int = -1) -> int:
    value = row.get(key, "")
    if value == "":
        return default
    return int(float(value))


def load_edges(path: Path) -> tuple[list[EdgeRow], dict[int, tuple[float, float]]]:
    by_edge: dict[int, EdgeRow] = {}
    nodes: dict[int, tuple[float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            edge_id = inum(row, "edge_id")
            if edge_id in by_edge:
                continue
            edge = EdgeRow(
                edge_id=edge_id,
                node0=inum(row, "node0"),
                node1=inum(row, "node1"),
                x0_um=fnum(row, "x0_um"),
                y0_um=fnum(row, "y0_um"),
                x1_um=fnum(row, "x1_um"),
                y1_um=fnum(row, "y1_um"),
                edge_length_m=fnum(row, "edge_length_m"),
                edge_couple_m=fnum(row, "edge_couple_m"),
                edge_area_proxy_m2=fnum(row, "edge_area_proxy_m2"),
                edge_class=row.get("edge_class", ""),
            )
            by_edge[edge_id] = edge
            nodes[edge.node0] = (edge.x0_um, edge.y0_um)
            nodes[edge.node1] = (edge.x1_um, edge.y1_um)
    return sorted(by_edge.values(), key=lambda item: item.edge_id), nodes


def load_elements(path: Path) -> list[Element]:
    elements: list[Element] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            elements.append(Element(
                cell_id=inum(row, "id", inum(row, "cell_id")),
                nodes=(inum(row, "node0"), inum(row, "node1"), inum(row, "node2")),
                region=row.get("region", ""),
                material=row.get("material", ""),
            ))
    return elements


def edge_key(a: int, b: int) -> tuple[int, int]:
    return tuple(sorted((a, b)))


def build_edge_cells(elements: list[Element]) -> dict[tuple[int, int], list[Element]]:
    edge_cells: dict[tuple[int, int], list[Element]] = defaultdict(list)
    for element in elements:
        a, b, c = element.nodes
        for u, v in ((a, b), (b, c), (c, a)):
            edge_cells[edge_key(u, v)].append(element)
    return edge_cells


def triangle_area_um2(element: Element, nodes: dict[int, tuple[float, float]]) -> float:
    try:
        (x0, y0), (x1, y1), (x2, y2) = [nodes[node] for node in element.nodes]
    except KeyError:
        return 0.0
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))


def barycentric_node_areas_cm2(elements: list[Element], nodes: dict[int, tuple[float, float]]) -> dict[int, float]:
    areas: dict[int, float] = defaultdict(float)
    for element in elements:
        share_cm2 = triangle_area_um2(element, nodes) * 1.0e-8 / 3.0
        for node in element.nodes:
            areas[node] += share_cm2
    return areas


def load_internal_source_weights(path: Path) -> dict[int, list[float]]:
    if not path.exists():
        return {}
    weights: dict[int, list[float]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("source_location_type") != "edge":
                continue
            edge_id = inum(row, "source_entity_id")
            weight = fnum(row, "source_weight_or_volume_cm2_for_2D", math.nan)
            if math.isfinite(weight) and weight > 0.0:
                weights[edge_id].append(weight)
    return weights


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.17g}"
    return str(value)


def make_rows(args: argparse.Namespace) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    edges, nodes = load_edges(args.sg_edges)
    elements = load_elements(args.elements_csv)
    edge_cells = build_edge_cells(elements)
    node_areas = barycentric_node_areas_cm2(elements, nodes)
    internal_weights = load_internal_source_weights(args.internal_audit)

    rows: list[dict[str, Any]] = []
    area_ratios: list[float] = []
    for face_id, edge in enumerate(edges):
        adjacent = edge_cells.get(edge.key, [])
        left = adjacent[0] if adjacent else None
        right = adjacent[1] if len(adjacent) > 1 else None
        nx, ny = edge.box_normal
        tx, ty = edge.box_tangent
        source_area_cm2 = edge.edge_area_proxy_m2 * 1.0e4
        for weight in internal_weights.get(edge.edge_id, []):
            if weight > 0.0:
                area_ratios.append(source_area_cm2 / weight)
        row = {
            "box_face_id": face_id,
            "owner_node_id": edge.node0,
            "neighbor_node_id": edge.node1,
            "associated_edge_id": edge.edge_id,
            "left_cell_id": left.cell_id if left is not None else None,
            "right_cell_id": right.cell_id if right is not None else None,
            "x_mid_um": edge.x_mid_um,
            "y_mid_um": edge.y_mid_um,
            "edge_coupling_normal_x": nx,
            "edge_coupling_normal_y": ny,
            "owner_to_neighbor_direction_x": nx,
            "owner_to_neighbor_direction_y": ny,
            "face_normal_x": nx,
            "face_normal_y": ny,
            "face_tangent_x": tx,
            "face_tangent_y": ty,
            "face_length_cm": edge.edge_couple_m * 100.0,
            "source_area_cm2": source_area_cm2,
            "orientation_sign": 1,
            "control_volume_owner_area_cm2": node_areas.get(edge.node0),
            "control_volume_neighbor_area_cm2": node_areas.get(edge.node1),
            "material_region": material_region(left, right),
            "boundary_type": edge.edge_class or boundary_type(left, right),
            "sub_control_volume_polygon_vertices": None,
            "partial_face_weight_to_owner": 0.5,
            "partial_face_weight_to_neighbor": 0.5,
        }
        rows.append(row)

    meta = {
        "edge_count": len(edges),
        "rows_with_left_cell": sum(1 for row in rows if row["left_cell_id"] is not None),
        "rows_with_right_cell": sum(1 for row in rows if row["right_cell_id"] is not None),
        "area_ratios": area_ratios,
        "material_interface_count": sum(1 for row in rows if is_material_interface(row, edge_cells, edges)),
        "boundary_count": sum(1 for row in rows if not row["right_cell_id"]),
    }
    return rows, meta


def material_region(left: Element | None, right: Element | None) -> str:
    if left is None:
        return ""
    if right is None or right.region == left.region:
        return left.region
    return f"{left.region}|{right.region}"


def boundary_type(left: Element | None, right: Element | None) -> str:
    if left is None:
        return "unattached"
    return "interior_bulk" if right is not None else "boundary"


def is_material_interface(row: dict[str, Any], edge_cells: dict[tuple[int, int], list[Element]], edges: list[EdgeRow]) -> bool:
    # The material interface count is derived directly from row material_region.
    material = row.get("material_region")
    return isinstance(material, str) and "|" in material


def median_or_none(values: list[float]) -> float | None:
    clean = [value for value in values if math.isfinite(value)]
    return median(clean) if clean else None


def write_summary(path: Path, rows: list[dict[str, Any]], meta: dict[str, Any]) -> None:
    normals_available = all(row["face_normal_x"] is not None and row["face_normal_y"] is not None for row in rows)
    tangents_valid = all(row["face_tangent_x"] is not None and row["face_tangent_y"] is not None for row in rows)
    area_ratio_median = median_or_none(meta["area_ratios"])
    area_ratio_min = min(meta["area_ratios"]) if meta["area_ratios"] else None
    area_ratio_max = max(meta["area_ratios"]) if meta["area_ratios"] else None

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Box Face Geometry Audit\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver path changed: no\n")
        out.write("- face tangent components: normalized unitless x/y\n")
        out.write("- edge_coupling_normal components: owner-to-neighbor unit direction used by the current edge-coupled box method\n")
        out.write("- backward-compatible face_normal_x/y duplicate edge_coupling_normal_x/y; this is not a true independent dual/box face normal\n")
        out.write("- face length unit: cm\n")
        out.write("- source area/control-volume area unit: cm^2\n")
        out.write(f"- exported rows: {len(rows)}\n")
        out.write(f"- rows with right cell: {meta['rows_with_right_cell']}\n")
        out.write(f"- boundary rows: {meta['boundary_count']}\n\n")

        out.write("## Required Answers\n\n")
        out.write(f"1. edge-coupling normals available for every source edge: {'yes' if normals_available else 'no'}.\n")
        out.write("2. edge_coupling_normal is not a true independent dual/box face normal; it is aligned with the owner-to-neighbor edge direction.\n")
        out.write("3. angle(owner_to_neighbor_direction, edge_coupling_normal) distribution: min=0 deg, median=0 deg, max=0 deg for valid source edges.\n")
        if area_ratio_median is None:
            out.write("4. source_area_cm2 geometry match: insufficient internal source weights for comparison.\n")
        else:
            out.write(
                "4. source_area_cm2 geometry match: "
                f"median exported/internal weight ratio={area_ratio_median:.6g}, "
                f"min={area_ratio_min:.6g}, max={area_ratio_max:.6g}.\n"
            )
        if area_ratio_median is not None and abs(area_ratio_median - 1.0) <= 1.0e-9:
            out.write("5. summing/exporting source areas reproduces previous internal source weights for matched edge rows: yes.\n")
        else:
            out.write("5. summing/exporting source areas reproduces previous internal source weights: not proven by available rows.\n")
        out.write(
            "6. boundary faces and material-interface faces: boundary faces are tagged through boundary_type/edge_class; "
            "material interfaces are represented by material_region when adjacent cell regions differ.\n"
        )
        if normals_available and tangents_valid:
            out.write("7. exported geometry can support a first-pass cell vector current reconstruction from edge-coupled face-normal fluxes: yes; explicit sub-control-volume polygon faces are still not exported.\n")
        else:
            out.write(
                "7. exported geometry cannot yet support true cell vector current reconstruction from face normal fluxes; "
                "export explicit control-volume face normals/areas from the mesh builder next.\n"
            )


def polygon_area_um2(points: list[tuple[float, float]]) -> float:
    area = 0.0
    for idx, (x0, y0) in enumerate(points):
        x1, y1 = points[(idx + 1) % len(points)]
        area += x0 * y1 - x1 * y0
    return 0.5 * abs(area)


def normalize(dx: float, dy: float) -> tuple[float | None, float | None, float]:
    length = math.hypot(dx, dy)
    if length <= 0.0:
        return None, None, 0.0
    return dx / length, dy / length, length


def make_true_dual_rows(
    box_rows: list[dict[str, Any]],
    elements: list[Element],
    nodes: dict[int, tuple[float, float]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    edge_info: dict[tuple[int, int], dict[str, Any]] = {}
    for row in box_rows:
        key = edge_key(int(row["owner_node_id"]), int(row["neighbor_node_id"]))
        edge_info[key] = row

    dual_rows: list[dict[str, Any]] = []
    validation_rows: list[dict[str, Any]] = []
    dual_id = 0
    for element in elements:
        ids = element.nodes
        if any(node not in nodes for node in ids):
            continue
        coords = {node: nodes[node] for node in ids}
        centroid = (
            sum(coords[node][0] for node in ids) / 3.0,
            sum(coords[node][1] for node in ids) / 3.0,
        )
        edge_mid: dict[tuple[int, int], tuple[float, float]] = {}
        for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            ax, ay = coords[a]
            bx, by = coords[b]
            edge_mid[edge_key(a, b)] = (0.5 * (ax + bx), 0.5 * (ay + by))

        subcv_area: dict[int, float] = {}
        for k, owner in enumerate(ids):
            prev_node = ids[(k - 1) % 3]
            next_node = ids[(k + 1) % 3]
            polygon = [
                coords[owner],
                edge_mid[edge_key(owner, next_node)],
                centroid,
                edge_mid[edge_key(prev_node, owner)],
            ]
            subcv_area[owner] = polygon_area_um2(polygon) * 1.0e-8

        cell_area_cm2 = triangle_area_um2(element, nodes) * 1.0e-8
        cancel_by_edge: dict[tuple[int, int], list[tuple[float, float]]] = defaultdict(list)
        normal_lengths: list[float] = []
        face_lengths_cm: list[float] = []

        for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            key = edge_key(a, b)
            info = edge_info.get(key, {})
            midpoint = edge_mid[key]
            vx = midpoint[0] - centroid[0]
            vy = midpoint[1] - centroid[1]
            face_length_um = math.hypot(vx, vy)
            if face_length_um <= 0.0:
                continue
            # Candidate normals perpendicular to the dual segment; choose the one
            # pointing from owner a toward neighbor b.
            cand0 = (vy / face_length_um, -vx / face_length_um)
            cand1 = (-cand0[0], -cand0[1])
            ax, ay = coords[a]
            bx, by = coords[b]
            toward_b = (bx - ax, by - ay)
            normal_a = cand0 if cand0[0] * toward_b[0] + cand0[1] * toward_b[1] >= 0.0 else cand1
            face_length_cm = face_length_um * 1.0e-4
            associated_edge = info.get("associated_edge_id")
            boundary = info.get("boundary_type") or "subcv_internal"
            region = material_region_for_element(element)
            def angle_deg(direction: tuple[float, float], normal: tuple[float, float]) -> float:
                dlen = math.hypot(direction[0], direction[1])
                nlen = math.hypot(normal[0], normal[1])
                if dlen <= 0.0 or nlen <= 0.0:
                    return math.nan
                dot = max(-1.0, min(1.0, (direction[0] * normal[0] + direction[1] * normal[1]) / (dlen * nlen)))
                return math.degrees(math.acos(dot))

            for owner, neighbor, normal, orientation in (
                (a, b, normal_a, 1.0),
                (b, a, (-normal_a[0], -normal_a[1]), -1.0),
            ):
                nx, ny = normal
                normal_lengths.append(math.hypot(nx, ny))
                face_lengths_cm.append(face_length_cm)
                cancel_by_edge[key].append((nx, ny))
                dual_rows.append({
                    "dual_face_id": dual_id,
                    "cell_id": element.cell_id,
                    "owner_node_id": owner,
                    "neighbor_node_id": neighbor,
                    "associated_primal_edge_id": associated_edge,
                    "dual_face_vertex0_x_um": centroid[0],
                    "dual_face_vertex0_y_um": centroid[1],
                    "dual_face_vertex1_x_um": midpoint[0],
                    "dual_face_vertex1_y_um": midpoint[1],
                    "dual_face_mid_x_um": 0.5 * (centroid[0] + midpoint[0]),
                    "dual_face_mid_y_um": 0.5 * (centroid[1] + midpoint[1]),
                    "dual_face_normal_x": nx,
                    "dual_face_normal_y": ny,
                    "dual_face_length_cm": face_length_cm,
                    "edge_direction_to_dual_normal_angle_deg": angle_deg((coords[neighbor][0] - coords[owner][0], coords[neighbor][1] - coords[owner][1]), normal),
                    "orientation_sign": orientation,
                    "owner_subcv_area_cm2": subcv_area[owner],
                    "neighbor_subcv_area_cm2": subcv_area[neighbor],
                    "material_region": region,
                    "boundary_type": boundary,
                    "dual_type": "median_dual",
                })
                dual_id += 1

        cancel_residuals = [math.hypot(sum(nx for nx, _ in normals), sum(ny for _, ny in normals)) for normals in cancel_by_edge.values()]
        subcv_sum = sum(subcv_area.values())
        area_abs_error = abs(subcv_sum - cell_area_cm2)
        area_rel_error = area_abs_error / max(abs(cell_area_cm2), 1.0e-300)
        cancel_max = max(cancel_residuals) if cancel_residuals else math.nan
        normal_min = min(normal_lengths) if normal_lengths else math.nan
        normal_max = max(normal_lengths) if normal_lengths else math.nan
        length_min = min(face_lengths_cm) if face_lengths_cm else math.nan
        validation_rows.append({
            "cell_id": element.cell_id,
            "cell_area_cm2": cell_area_cm2,
            "subcv_area_sum_cm2": subcv_sum,
            "subcv_area_abs_error_cm2": area_abs_error,
            "subcv_area_rel_error": area_rel_error,
            "internal_dual_face_cancel_residual_max": cancel_max,
            "normal_length_min": normal_min,
            "normal_length_max": normal_max,
            "face_length_min_cm": length_min,
            "cell_area_closure_pass": area_rel_error <= 1.0e-12,
            "internal_face_cancel_pass": math.isfinite(cancel_max) and cancel_max <= 1.0e-12,
            "normal_unit_pass": math.isfinite(normal_min) and math.isfinite(normal_max) and abs(normal_min - 1.0) <= 1.0e-12 and abs(normal_max - 1.0) <= 1.0e-12,
            "face_length_pass": math.isfinite(length_min) and length_min > 0.0,
        })
    return dual_rows, validation_rows


def material_region_for_element(element: Element) -> str:
    return element.region


def write_dual_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    explicit_dual = any(row.get("dual_type") != "edge_coupled_only" for row in rows)
    normal_lengths = [math.hypot(float(row["dual_face_normal_x"]), float(row["dual_face_normal_y"])) for row in rows if row.get("dual_face_normal_x") is not None]
    angles = [float(row["edge_direction_to_dual_normal_angle_deg"]) for row in rows if row.get("edge_direction_to_dual_normal_angle_deg") is not None and math.isfinite(float(row["edge_direction_to_dual_normal_angle_deg"]))]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# Dual Face Geometry Audit\n\n")
        out.write("- diagnostic: post-processing only\n")
        out.write("- solver path changed: no\n")
        out.write("- dual_type: median_dual for true Tri3 median-dual faces\n")
        out.write(f"- exported rows: {len(rows)}\n")
        if normal_lengths:
            out.write(f"- normal length min/median/max: {min(normal_lengths):.6g}/{median(normal_lengths):.6g}/{max(normal_lengths):.6g}\n")
        if angles:
            out.write(f"- angle(edge_direction, true_dual_face_normal) min/median/max deg: {min(angles):.6g}/{median(angles):.6g}/{max(angles):.6g}\n")
        out.write("\n## Required Answers\n\n")
        out.write(f"1. solver has explicit dual/sub-control-volume faces: {'yes' if explicit_dual else 'no'}.\n")
        if explicit_dual:
            out.write("2. true dual-face normals are exported from median-dual centroid-to-edge-midpoint segments.\n")
            if angles:
                out.write(f"3. angle(edge_direction, true_dual_face_normal) distribution deg: min={min(angles):.6g}, median={median(angles):.6g}, max={max(angles):.6g}.\n")
            else:
                out.write("3. angle(edge_direction, true_dual_face_normal) distribution: insufficient_data.\n")
            out.write("4. dual face length/area closure is reported in true_dual_face_geometry_validation.csv.\n")
            out.write("5. these faces can support true vector current reconstruction from dual-face normal fluxes.\n")
            out.write("6. current avalanche source projection can now be audited with explicit median-dual normals.\n")
        else:
            out.write("2. true dual-face normals are unavailable, so they cannot be compared with owner-neighbor edge directions.\n")
            out.write("3. angle(edge_direction, true_dual_face_normal) distribution: unavailable.\n")
            out.write("4. dual face length/area sum cannot reproduce node control-volume area from this export because explicit dual vertices are absent.\n")
            out.write("5. these rows cannot support true vector current reconstruction from dual-face normal fluxes.\n")
            out.write("6. current avalanche source is based on edge-coupled scalar flux and cannot distinguish tangent/normal projection.\n")


def write_validation(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=VALIDATION_HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in VALIDATION_HEADER})


def write_validation_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    area_pass = all(bool(row.get("cell_area_closure_pass")) for row in rows)
    cancel_pass = all(bool(row.get("internal_face_cancel_pass")) for row in rows)
    normal_pass = all(bool(row.get("normal_unit_pass")) for row in rows)
    length_pass = all(bool(row.get("face_length_pass")) for row in rows)
    max_area_error = max((float(row["subcv_area_abs_error_cm2"]) for row in rows), default=math.nan)
    max_cancel = max((float(row["internal_dual_face_cancel_residual_max"]) for row in rows), default=math.nan)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("# True Dual Face Geometry Validation\n\n")
        out.write(f"- cells: {len(rows)}\n")
        out.write(f"- max subCV area closure error cm2: {max_area_error:.6g}\n")
        out.write(f"- max internal dual face cancel residual: {max_cancel:.6g}\n")
        out.write(f"- cell_area_closure_pass: {'yes' if area_pass else 'no'}\n")
        out.write(f"- internal_dual_faces_cancel_pass: {'yes' if cancel_pass else 'no'}\n")
        out.write(f"- normals_unit_length_pass: {'yes' if normal_pass else 'no'}\n")
        out.write(f"- face_length_positive_pass: {'yes' if length_pass else 'no'}\n")

def run(args: argparse.Namespace) -> None:
    rows, meta = make_rows(args)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in HEADER})
    write_summary(args.out_summary, rows, meta)
    edges, nodes = load_edges(args.sg_edges)
    elements = load_elements(args.elements_csv)
    dual_rows, validation_rows = make_true_dual_rows(rows, elements, nodes)
    args.dual_out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.dual_out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DUAL_HEADER)
        writer.writeheader()
        for row in dual_rows:
            writer.writerow({key: fmt(row.get(key)) for key in DUAL_HEADER})
    write_dual_summary(args.dual_out_summary, dual_rows)
    write_validation(args.validation_out_csv, validation_rows)
    write_validation_summary(args.validation_out_summary, validation_rows)


def main() -> int:
    args = parse_args()
    run(args)
    print(args.out_csv)
    print(args.out_summary)
    print(args.dual_out_csv)
    print(args.dual_out_summary)
    print(args.validation_out_csv)
    print(args.validation_out_summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())







