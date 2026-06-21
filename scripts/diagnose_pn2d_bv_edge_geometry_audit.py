#!/usr/bin/env python3
"""Step 2: per-edge Scharfetter-Gummel avalanche geometry audit.

Step 1 (``diagnose_pn2d_bv_sentaurus_volume_source.py``) proved that the nodal
control-volume convention is not the cause of the ~0.373x Vela/Sentaurus
avalanche-source deficit: the barycentric (``area / 3``) and box/Voronoi
(``0.25 * h * couple``) lumped volumes agree (median ratio 1.0, total domain
area ratio 0.9993), and the source ratio stays at ~0.373x under both.

This script answers the remaining geometry question at the *edge* level instead
of inferring a factor through the terminal current. For every Scharfetter-Gummel
edge it compares:

* the C++ ``edge_area_proxy`` actually used by the avalanche source
  (``0.5 * h * couple``) against an independently recomputed box area from the
  shared mesh, to confirm the C++ geometry is faithful;
* the box area against the element-partition area (sum over adjacent triangles
  of the edge's circumcenter share), to check for a hidden boundary/obtuse
  factor;
* the Vela per-edge avalanche source against a Sentaurus per-edge proxy built
  from the shared-mesh nodal ``ImpactIonization`` field and the same box area.

All quantities are then broken down by ``edge_class`` (interior_bulk, boundary,
contact_edge) and weighted by the actual avalanche source, so the dominant
junction-hotspot edges drive the conclusion. The goal is to decide whether the
real deficit is geometric (concentrated on boundary/contact/obtuse edges) or in
the integrand (flux x alpha), which would redirect the fix away from any global
source-geometry multiplier.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag  # noqa: E402
import diagnose_pn2d_bv_sentaurus_volume_source as step1  # noqa: E402

Q = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-mesh", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=-13.2)
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
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


def load_cxx_edges(path: Path, bias: float) -> dict[tuple[int, int], dict[str, Any]]:
    """Authoritative C++ SG edge records at the requested bias, keyed by node
    pair (sorted)."""
    records: dict[tuple[int, int], dict[str, Any]] = {}
    for row in read_csv_rows(path):
        row_bias = optional_float(row.get("bias_V"))
        if row_bias is None or abs(row_bias - bias) > 1.0e-6:
            continue
        n0 = int(row["node0"])
        n1 = int(row["node1"])
        key = (n0, n1) if n0 <= n1 else (n1, n0)
        records[key] = {
            "edge_length_m": optional_float(row.get("edge_length_m")) or 0.0,
            "edge_couple_m": optional_float(row.get("edge_couple_m")) or 0.0,
            "edge_area_proxy_m2": optional_float(row.get("edge_area_proxy_m2")) or 0.0,
            "electric_field_V_per_m": optional_float(row.get("electric_field_V_per_m")) or 0.0,
            "electron_alpha_m_inv": optional_float(row.get("electron_alpha_m_inv")) or 0.0,
            "hole_alpha_m_inv": optional_float(row.get("hole_alpha_m_inv")) or 0.0,
            "electron_flux_proxy": optional_float(row.get("electron_flux_proxy")) or 0.0,
            "hole_flux_proxy": optional_float(row.get("hole_flux_proxy")) or 0.0,
            "edge_source_integral": optional_float(row.get("edge_source_integral")) or 0.0,
            "edge_class": row.get("edge_class") or "unknown",
        }
    return records


def edge_geometry(nodes: dict[int, dict[str, float]],
                  triangles: list[dict[str, Any]]) -> dict[tuple[int, int], dict[str, Any]]:
    """Independent per-edge geometry from the shared mesh.

    Reports, for each edge: box couple (Voronoi, with obtuse circumcenter
    fallback), the raw cotangent couple (no fallback), the element-partition
    area share, adjacent triangle count, and obtuse-adjacent count.
    """
    geom: dict[tuple[int, int], dict[str, Any]] = {}
    for triangle in triangles:
        ids = triangle["node_ids"]
        points = [sgdiag.point_m(nodes, nid) for nid in ids]
        area = sgdiag.triangle_area(points)
        if area <= 1.0e-30:
            continue
        for k in range(3):
            a = ids[k]
            b = ids[(k + 1) % 3]
            opp = ids[(k + 2) % 3]
            key = (a, b) if a <= b else (b, a)
            entry = geom.setdefault(key, {
                "length_m": sgdiag.distance(sgdiag.point_m(nodes, a),
                                            sgdiag.point_m(nodes, b)),
                "couple_box_m": 0.0,
                "couple_raw_m": 0.0,
                "element_area_share_m2": 0.0,
                "adjacent_cell_count": 0,
                "obtuse_adjacent_count": 0,
            })
            h = entry["length_m"]
            if h <= 1.0e-30:
                continue
            cot = sgdiag.cotangent_at_opposite(
                sgdiag.point_m(nodes, a), sgdiag.point_m(nodes, b),
                sgdiag.point_m(nodes, opp))
            raw_couple = 0.5 * cot * h
            box_couple = raw_couple
            if cot < 0.0:
                box_couple = area / (3.0 * h)
                entry["obtuse_adjacent_count"] += 1
            entry["couple_raw_m"] += raw_couple
            entry["couple_box_m"] += max(box_couple, 0.0)
            # Box-method area share of this triangle assigned to this edge:
            # 0.5 * h * (circumcenter perpendicular share) == 0.5 * h * box_couple_local.
            entry["element_area_share_m2"] += 0.5 * h * max(box_couple, 0.0)
            entry["adjacent_cell_count"] += 1
    for entry in geom.values():
        entry["box_area_m2"] = 0.5 * entry["length_m"] * entry["couple_box_m"]
    return geom


def safe_ratio(num: float, den: float) -> float | None:
    if den == 0.0:
        return None
    return num / den


def stat_block(values: list[float]) -> dict[str, float] | None:
    clean = [v for v in values if v is not None and math.isfinite(v)]
    if not clean:
        return None
    return {
        "count": len(clean),
        "median": statistics.median(clean),
        "min": min(clean),
        "max": max(clean),
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes, triangles, _ = sgdiag.read_mesh(args.vela_mesh)
    num_nodes = len(nodes)
    geom = edge_geometry(nodes, triangles)
    cxx = load_cxx_edges(args.sg_edge_csv, args.bias)

    # Sentaurus nodal generation density on the shared mesh (m^-3 s^-1).
    impact = step1.load_sentaurus_impact(args.sentaurus_dir)
    sent_gen = [impact.get(nid, 0.0) * 1.0e6 for nid in range(num_nodes)]

    fields = [
        "node0", "node1", "x_mid_um", "y_mid_um", "edge_class",
        "adjacent_cell_count", "obtuse_adjacent_count",
        "length_m", "couple_box_m", "couple_raw_m",
        "box_area_m2", "cxx_area_proxy_m2", "cxx_over_py_box_area",
        "element_area_share_m2", "box_over_element_share",
        "electric_field_V_per_m", "edge_source_integral_s_inv",
        "sentaurus_edge_source_s_inv", "vela_over_sentaurus_edge_source",
    ]

    rows_out: list[dict[str, Any]] = []
    by_class: dict[str, dict[str, Any]] = defaultdict(lambda: {
        "edge_count": 0,
        "vela_source_sum": 0.0,
        "sentaurus_source_sum": 0.0,
        "box_area_sum": 0.0,
        "element_area_share_sum": 0.0,
        "cxx_over_py_area": [],
        "box_over_element": [],
        "obtuse_edge_count": 0,
    })

    for key, record in cxx.items():
        g = geom.get(key)
        if g is None:
            continue
        n0, n1 = key
        mid_x = 0.5 * (nodes[n0]["x_um"] + nodes[n1]["x_um"])
        mid_y = 0.5 * (nodes[n0]["y_um"] + nodes[n1]["y_um"])
        box_area = g["box_area_m2"]
        cxx_area = record["edge_area_proxy_m2"]
        element_share = g["element_area_share_m2"]
        vela_src = record["edge_source_integral"]
        sent_src = 0.5 * (sent_gen[n0] + sent_gen[n1]) * box_area
        edge_class = record["edge_class"]

        agg = by_class[edge_class]
        agg["edge_count"] += 1
        agg["vela_source_sum"] += vela_src
        agg["sentaurus_source_sum"] += sent_src
        agg["box_area_sum"] += box_area
        agg["element_area_share_sum"] += element_share
        agg["cxx_over_py_area"].append(safe_ratio(cxx_area, box_area))
        agg["box_over_element"].append(safe_ratio(box_area, element_share))
        if g["obtuse_adjacent_count"] > 0:
            agg["obtuse_edge_count"] += 1

        rows_out.append({
            "node0": n0,
            "node1": n1,
            "x_mid_um": mid_x,
            "y_mid_um": mid_y,
            "edge_class": edge_class,
            "adjacent_cell_count": g["adjacent_cell_count"],
            "obtuse_adjacent_count": g["obtuse_adjacent_count"],
            "length_m": g["length_m"],
            "couple_box_m": g["couple_box_m"],
            "couple_raw_m": g["couple_raw_m"],
            "box_area_m2": box_area,
            "cxx_area_proxy_m2": cxx_area,
            "cxx_over_py_box_area": safe_ratio(cxx_area, box_area),
            "element_area_share_m2": element_share,
            "box_over_element_share": safe_ratio(box_area, element_share),
            "electric_field_V_per_m": record["electric_field_V_per_m"],
            "edge_source_integral_s_inv": vela_src,
            "sentaurus_edge_source_s_inv": sent_src,
            "vela_over_sentaurus_edge_source": safe_ratio(vela_src, sent_src),
        })

    rows_out.sort(key=lambda r: r["edge_source_integral_s_inv"], reverse=True)
    detail_path = args.out_dir / "edge_geometry_audit_edges.csv"
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_out)

    total_vela = sum(r["edge_source_integral_s_inv"] for r in rows_out)
    total_sent = sum(r["sentaurus_edge_source_s_inv"] for r in rows_out)

    class_summary: dict[str, Any] = {}
    for edge_class, agg in by_class.items():
        class_summary[edge_class] = {
            "edge_count": agg["edge_count"],
            "obtuse_edge_count": agg["obtuse_edge_count"],
            "vela_source_sum_s_inv": agg["vela_source_sum"],
            "sentaurus_source_sum_s_inv": agg["sentaurus_source_sum"],
            "vela_source_fraction": safe_ratio(agg["vela_source_sum"], total_vela),
            "sentaurus_source_fraction": safe_ratio(agg["sentaurus_source_sum"], total_sent),
            "vela_over_sentaurus_source": safe_ratio(
                agg["vela_source_sum"], agg["sentaurus_source_sum"]),
            "box_area_sum_m2": agg["box_area_sum"],
            "element_area_share_sum_m2": agg["element_area_share_sum"],
            "box_over_element_area_sum": safe_ratio(
                agg["box_area_sum"], agg["element_area_share_sum"]),
            "cxx_over_py_box_area": stat_block(agg["cxx_over_py_area"]),
            "box_over_element_share": stat_block(agg["box_over_element"]),
        }

    # Source-weighted geometry consistency over all edges: if the C++ geometry
    # faithfully reproduces the independent box area, this is ~1.0 everywhere.
    cxx_over_py_all = stat_block([r["cxx_over_py_box_area"] for r in rows_out])
    box_over_element_all = stat_block([r["box_over_element_share"] for r in rows_out])

    # Where does the avalanche source live? Top-source edges and their class.
    top_edges = rows_out[:25]
    top_boundary = sum(
        1 for r in top_edges if r["edge_class"] in ("boundary", "contact_edge"))

    report: dict[str, Any] = {
        "bias_V": args.bias,
        "num_nodes": num_nodes,
        "matched_edges": len(rows_out),
        "total_vela_source_s_inv": total_vela,
        "total_sentaurus_edge_proxy_s_inv": total_sent,
        "vela_over_sentaurus_total": safe_ratio(total_vela, total_sent),
        "cxx_over_py_box_area_all": cxx_over_py_all,
        "box_over_element_share_all": box_over_element_all,
        "top25_source_edges_on_boundary_or_contact": top_boundary,
        "by_edge_class": class_summary,
        "interpretation": {
            "geometry_faithful": (
                cxx_over_py_all is not None
                and abs(cxx_over_py_all["median"] - 1.0) < 1.0e-3
                and abs(cxx_over_py_all["min"] - 1.0) < 1.0e-3
                and abs(cxx_over_py_all["max"] - 1.0) < 1.0e-3),
            "box_equals_element_share": (
                box_over_element_all is not None
                and abs(box_over_element_all["median"] - 1.0) < 1.0e-3),
        },
    }

    summary_path = args.out_dir / "edge_geometry_audit_summary.json"
    summary_path.write_text(json.dumps(report, indent=2))

    print(json.dumps(report, indent=2))
    print(f"\nWrote {summary_path}")
    print(f"Wrote {detail_path}")


if __name__ == "__main__":
    main()
