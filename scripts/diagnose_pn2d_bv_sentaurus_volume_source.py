#!/usr/bin/env python3
"""Step 1: reconstruct the PN2D BV avalanche source integral with each side's
own consistent control volume, instead of multiplying the Sentaurus nodal
generation density by a single borrowed volume.

Vela and Sentaurus share the same mesh for this fixture, so the remaining
ambiguity is the *area convention*: the C++ Scharfetter-Gummel avalanche source
integrates each edge with the box-method area ``0.5 * h * couple`` (so each node
sees ``0.25 * h * couple`` per incident edge, i.e. the Voronoi/box area), while
the parity comparison multiplies the Sentaurus nodal ``ImpactIonization`` field
by the barycentric lumped area ``area / 3``.

This script integrates the total avalanche generation three independent ways and
reports how much the Vela/Sentaurus ratio moves with the area convention, then
cross-checks the generation-rate current proxy ``q * integral(G dV)`` against the
known terminal-current ratio (Vela/Sentaurus ~= 0.698 at -13.2 V).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

import diagnose_pn2d_bv_sg_avalanche_edges as sgdiag  # noqa: E402


Q = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-mesh", type=Path, required=True)
    parser.add_argument("--vela-vtk", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--sentaurus-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=-13.2)
    parser.add_argument(
        "--active-fraction",
        type=float,
        default=1.0e-3,
        help="Active support = nodes whose generation density exceeds this "
        "fraction of the per-side maximum.",
    )
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


def triangle_area(points: list[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = points
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))


def barycentric_volumes(nodes: dict[int, dict[str, float]],
                        triangles: list[dict[str, Any]]) -> list[float]:
    """Lumped area / 3 control volume (what the parity comparison uses)."""
    volumes = [0.0 for _ in nodes]
    for triangle in triangles:
        ids = triangle["node_ids"]
        area = triangle_area([sgdiag.point_m(nodes, nid) for nid in ids])
        for nid in ids:
            volumes[nid] += area / 3.0
    return volumes


def voronoi_volumes(num_nodes: int, edges: list[dict[str, Any]]) -> list[float]:
    """Box-method control volume sum_incident 0.25 * h * couple (what the C++
    SG avalanche source effectively assigns to each node)."""
    volumes = [0.0 for _ in range(num_nodes)]
    for edge in edges:
        h = float(edge["length_m"])
        couple = float(edge["couple_m"])
        if h <= 1.0e-30 or couple <= 0.0:
            continue
        half_box = 0.25 * h * couple
        volumes[int(edge["node0"])] += half_box
        volumes[int(edge["node1"])] += half_box
    return volumes


def load_sentaurus_impact(sentaurus_dir: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    field_dir = sentaurus_dir / "fields"
    for path in sorted(field_dir.glob("ImpactIonization_region*.csv")):
        for row in read_csv_rows(path):
            value = optional_float(row.get("component0"))
            if value is None:
                # fall back to first non-id numeric column
                for key, raw in row.items():
                    if key == "node_id":
                        continue
                    value = optional_float(raw)
                    if value is not None:
                        break
            if value is not None:
                result[int(row["node_id"])] = value
    return result


def boundary_nodes(num_nodes: int, edges: list[dict[str, Any]]) -> set[int]:
    boundary: set[int] = set()
    for edge in edges:
        if len(edge["cell_ids"]) == 1:
            boundary.add(int(edge["node0"]))
            boundary.add(int(edge["node1"]))
    return boundary


def load_sg_edge_totals(path: Path, bias: float) -> dict[str, float]:
    edge_total = 0.0
    node_total = 0.0
    matched_rows = 0
    for row in read_csv_rows(path):
        row_bias = optional_float(row.get("bias_V"))
        if row_bias is None or abs(row_bias - bias) > 1.0e-6:
            continue
        matched_rows += 1
        edge_total += optional_float(row.get("edge_source_integral")) or 0.0
        node_total += optional_float(row.get("node0_source_integral")) or 0.0
        node_total += optional_float(row.get("node1_source_integral")) or 0.0
    return {
        "edge_rows": matched_rows,
        "cxx_edge_source_total_s_inv": edge_total,
        "cxx_node_source_total_s_inv": node_total,
    }


def integrate(values: list[float], volumes: list[float],
              node_filter: set[int] | None) -> float:
    total = 0.0
    for nid, (val, vol) in enumerate(zip(values, volumes)):
        if node_filter is not None and nid not in node_filter:
            continue
        total += val * vol
    return total


def support_set(density: list[float], fraction: float) -> set[int]:
    peak = max(density) if density else 0.0
    if peak <= 0.0:
        return set()
    threshold = peak * fraction
    return {nid for nid, val in enumerate(density) if val > threshold}


def safe_ratio(num: float, den: float) -> float | None:
    if den == 0.0:
        return None
    return num / den


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes, triangles, _ = sgdiag.read_mesh(args.vela_mesh)
    edges = sgdiag.build_edges(nodes, triangles)
    num_nodes = len(nodes)

    vol_lump = barycentric_volumes(nodes, triangles)
    vol_voronoi = voronoi_volumes(num_nodes, edges)
    total_area_lump = sum(vol_lump)
    total_area_voronoi = sum(vol_voronoi)

    scalars = sgdiag.parse_vtk_scalars(args.vela_vtk)
    vela_gen = scalars.get("AvalancheGeneration")
    if vela_gen is None or len(vela_gen) != num_nodes:
        raise RuntimeError("Vela VTK missing AvalancheGeneration for all nodes")

    impact = load_sentaurus_impact(args.sentaurus_dir)
    sent_gen = [impact.get(nid, 0.0) * 1.0e6 for nid in range(num_nodes)]

    sg_totals = load_sg_edge_totals(args.sg_edge_csv, args.bias)

    boundary = boundary_nodes(num_nodes, edges)
    vela_support = support_set(vela_gen, args.active_fraction)
    sent_support = support_set(sent_gen, args.active_fraction)
    union_support = vela_support | sent_support

    def integrals(node_filter: set[int] | None) -> dict[str, Any]:
        vela_lump = integrate(vela_gen, vol_lump, node_filter)
        vela_vor = integrate(vela_gen, vol_voronoi, node_filter)
        sent_lump = integrate(sent_gen, vol_lump, node_filter)
        sent_vor = integrate(sent_gen, vol_voronoi, node_filter)
        return {
            "node_count": (num_nodes if node_filter is None else len(node_filter)),
            "vela_vtk_total_lump_s_inv": vela_lump,
            "vela_vtk_total_voronoi_s_inv": vela_vor,
            "sentaurus_total_lump_s_inv": sent_lump,
            "sentaurus_total_voronoi_s_inv": sent_vor,
            "vela_over_sentaurus_lump": safe_ratio(vela_lump, sent_lump),
            "vela_over_sentaurus_voronoi": safe_ratio(vela_vor, sent_vor),
            "vela_lump_over_voronoi": safe_ratio(vela_lump, vela_vor),
            "sentaurus_lump_over_voronoi": safe_ratio(sent_lump, sent_vor),
        }

    report: dict[str, Any] = {
        "bias_V": args.bias,
        "num_nodes": num_nodes,
        "num_edges": len(edges),
        "boundary_node_count": len(boundary),
        "total_area_lump_m2": total_area_lump,
        "total_area_voronoi_m2": total_area_voronoi,
        "total_area_lump_over_voronoi": safe_ratio(total_area_lump, total_area_voronoi),
        "active_fraction": args.active_fraction,
        "vela_support_node_count": len(vela_support),
        "sentaurus_support_node_count": len(sent_support),
        "union_support_node_count": len(union_support),
        "cxx_sg_edge_records": sg_totals,
        "integrals_all_nodes": integrals(None),
        "integrals_union_support": integrals(union_support),
    }

    # Consistency: does the C++ assembled SG edge source equal the VTK field
    # integrated with the box/Voronoi volume? It should, because the VTK field
    # is AvalancheGeneration = node_source_integral / node_box_volume.
    cxx_edge_total = sg_totals["cxx_edge_source_total_s_inv"]
    vela_vor_all = report["integrals_all_nodes"]["vela_vtk_total_voronoi_s_inv"]
    vela_lump_all = report["integrals_all_nodes"]["vela_vtk_total_lump_s_inv"]
    report["cxx_edge_total_over_vtk_voronoi"] = safe_ratio(cxx_edge_total, vela_vor_all)
    report["cxx_edge_total_over_vtk_lump"] = safe_ratio(cxx_edge_total, vela_lump_all)

    # Generation-rate current proxy q * integral(G dV) (per micron depth, A/um).
    # The Vela/Sentaurus ratio of this proxy should be compared to the terminal
    # current ratio. Geometry must be consistent on both sides for the ratio to
    # be meaningful, so we report both conventions.
    union = report["integrals_union_support"]
    report["generation_current_proxy_A_per_um"] = {
        "vela_lump": Q * union["vela_vtk_total_lump_s_inv"] * 1.0e-6,
        "vela_voronoi": Q * union["vela_vtk_total_voronoi_s_inv"] * 1.0e-6,
        "sentaurus_lump": Q * union["sentaurus_total_lump_s_inv"] * 1.0e-6,
        "sentaurus_voronoi": Q * union["sentaurus_total_voronoi_s_inv"] * 1.0e-6,
        "vela_over_sentaurus_lump": union["vela_over_sentaurus_lump"],
        "vela_over_sentaurus_voronoi": union["vela_over_sentaurus_voronoi"],
    }

    # Per-node detail for the union support, to localize where the convention
    # changes the integrand most.
    detail_path = args.out_dir / "sentaurus_volume_source_nodes.csv"
    fields = [
        "node_id", "x_um", "y_um", "is_boundary",
        "vol_lump_m2", "vol_voronoi_m2", "vol_lump_over_voronoi",
        "vela_generation_m3_s", "sentaurus_generation_m3_s",
        "vela_source_lump_s_inv", "vela_source_voronoi_s_inv",
        "sentaurus_source_lump_s_inv", "sentaurus_source_voronoi_s_inv",
        "vela_over_sentaurus_lump", "vela_over_sentaurus_voronoi",
    ]
    rows_out: list[dict[str, Any]] = []
    for nid in sorted(union_support,
                      key=lambda i: sent_gen[i] + vela_gen[i], reverse=True):
        node = nodes[nid]
        vl = vol_lump[nid]
        vv = vol_voronoi[nid]
        vela_l = vela_gen[nid] * vl
        vela_v = vela_gen[nid] * vv
        sent_l = sent_gen[nid] * vl
        sent_v = sent_gen[nid] * vv
        rows_out.append({
            "node_id": nid,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "is_boundary": int(nid in boundary),
            "vol_lump_m2": vl,
            "vol_voronoi_m2": vv,
            "vol_lump_over_voronoi": safe_ratio(vl, vv),
            "vela_generation_m3_s": vela_gen[nid],
            "sentaurus_generation_m3_s": sent_gen[nid],
            "vela_source_lump_s_inv": vela_l,
            "vela_source_voronoi_s_inv": vela_v,
            "sentaurus_source_lump_s_inv": sent_l,
            "sentaurus_source_voronoi_s_inv": sent_v,
            "vela_over_sentaurus_lump": safe_ratio(vela_l, sent_l),
            "vela_over_sentaurus_voronoi": safe_ratio(vela_v, sent_v),
        })
    with detail_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows_out)

    # Volume-convention statistics on the union support.
    lump_over_vor = [
        r["vol_lump_over_voronoi"] for r in rows_out
        if r["vol_lump_over_voronoi"] is not None
    ]
    if lump_over_vor:
        report["union_support_vol_lump_over_voronoi"] = {
            "median": statistics.median(lump_over_vor),
            "min": min(lump_over_vor),
            "max": max(lump_over_vor),
        }

    summary_path = args.out_dir / "sentaurus_volume_source_summary.json"
    summary_path.write_text(json.dumps(report, indent=2))

    print(json.dumps(report, indent=2))
    print(f"\nWrote {summary_path}")
    print(f"Wrote {detail_path}")


if __name__ == "__main__":
    main()
