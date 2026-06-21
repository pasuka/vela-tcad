#!/usr/bin/env python3
"""Localize the BV depletion-region SG flux-divergence delta to individual edges.

Step D established (per node) that when Sentaurus' converged quasi-Fermi levels
are substituted into Vela's discretized equations, the residual disagreement is
carried almost entirely by the Scharfetter-Gummel carrier-continuity flux
divergence (median d(flux) ~ +2.8e-16 electron / +3.1e-16 hole), while SRH and
avalanche stay matched.  Step D aggregated the flux PER NODE.

This Direction-3 probe re-evaluates the SAME two states with the new C++
``sg_edge_flux_probe`` (faithful to the residual edge loop) and emits the SG
continuity flux PER EDGE, so the +2.8e-16 per-node divergence delta can be
localized to specific edges / mesh locations / field regimes.

It reuses the prepared external states and probe configs that the Step D driver
already wrote under ``<stepd-dir>/configs`` and ``<stepd-dir>/states``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from pathlib import Path
from typing import Any

BASELINE = "vela_baseline"
SENTQF = "vela_psi_sentaurus_qf"
VARIANTS = [BASELINE, SENTQF]

FLOAT_COLUMNS = [
    "x0", "y0", "x1", "y1", "length_m", "couple_m", "net_doping_avg_m3",
    "ni0_m3", "ni1_m3", "psi0_V", "psi1_V", "phin0_V", "phin1_V",
    "phip0_V", "phip1_V", "electric_field_V_m",
    "electron_mobility_m2_V_s", "hole_mobility_m2_V_s",
    "electron_flux", "hole_flux",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stepd-dir", type=Path, required=True,
        help="Existing Step D output dir with configs/ and states/ subfolders.")
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--depletion-half-width-um", type=float, default=0.25,
        help="Half width (um) of the depletion window centered on the junction.")
    parser.add_argument(
        "--active-field-threshold", type=float, default=1.0e7,
        help="Edge |E| (V/m) above which an edge is counted as high-field active.")
    parser.add_argument("--top-n", type=int, default=25)
    return parser.parse_args()


def run_probe(runner: Path, stepd_dir: Path, out_dir: Path, variant: str) -> Path:
    src_cfg = stepd_dir / "configs" / f"{variant}_newton_carrier_term_probe.json"
    if not src_cfg.is_file():
        raise SystemExit(f"missing Step D config: {src_cfg}")
    cfg = json.loads(src_cfg.read_text())
    edge_csv = (out_dir / "edge_flux" / f"{variant}_edge_flux.csv").resolve()
    edge_csv.parent.mkdir(parents=True, exist_ok=True)
    cfg["simulation_type"] = "sg_edge_flux_probe"
    cfg["output_csv"] = str(edge_csv)
    cfg.pop("carrier_term_probe", None)
    dst_cfg = out_dir / "configs" / f"{variant}_sg_edge_flux_probe.json"
    dst_cfg.parent.mkdir(parents=True, exist_ok=True)
    dst_cfg.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(dst_cfg.resolve())],
        cwd=dst_cfg.parent, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {variant}: "
            f"{result.stderr.strip() or result.stdout.strip()}")
    return edge_csv


def read_edges(path: Path) -> dict[int, dict[str, Any]]:
    edges: dict[int, dict[str, Any]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            record: dict[str, Any] = {
                "edge_id": int(row["edge_id"]),
                "node0": int(row["node0"]),
                "node1": int(row["node1"]),
            }
            for key in FLOAT_COLUMNS:
                record[key] = float(row[key])
            edges[record["edge_id"]] = record
    return edges


def detect_junction_x(edges: dict[int, dict[str, Any]]) -> float:
    """Junction x (m): midpoint of the edge whose endpoints straddle net=0."""
    crossings: list[float] = []
    for e in edges.values():
        mid_x = 0.5 * (e["x0"] + e["x1"])
        crossings.append((abs(e["net_doping_avg_m3"]), mid_x))
    crossings.sort()
    # Average the few midpoints closest to net-doping zero for robustness.
    closest = [mx for _, mx in crossings[:5]]
    return statistics.fmean(closest)


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def summarize(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"count": 0, "median": None, "mean": None, "sum": None,
                "min": None, "max": None}
    return {
        "count": len(values),
        "median": statistics.median(values),
        "mean": statistics.fmean(values),
        "sum": math.fsum(values),
        "min": min(values),
        "max": max(values),
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    edges_by_variant: dict[str, dict[int, dict[str, Any]]] = {}
    for variant in VARIANTS:
        csv_path = run_probe(args.runner, args.stepd_dir, args.out_dir, variant)
        edges_by_variant[variant] = read_edges(csv_path)

    base = edges_by_variant[BASELINE]
    sent = edges_by_variant[SENTQF]
    common_ids = sorted(set(base) & set(sent))
    junction_x = detect_junction_x(base)
    half_width_m = args.depletion_half_width_um * 1.0e-6

    rows: list[dict[str, Any]] = []
    for eid in common_ids:
        b = base[eid]
        s = sent[eid]
        mid_x = 0.5 * (b["x0"] + b["x1"])
        d_e = s["electron_flux"] - b["electron_flux"]
        d_h = s["hole_flux"] - b["hole_flux"]
        rows.append({
            "edge_id": eid,
            "node0": b["node0"],
            "node1": b["node1"],
            "mid_x_um": mid_x * 1.0e6,
            "dist_to_junction_um": (mid_x - junction_x) * 1.0e6,
            "length_um": b["length_m"] * 1.0e6,
            "couple_um": b["couple_m"] * 1.0e6,
            "net_doping_avg_m3": b["net_doping_avg_m3"],
            "electric_field_V_m": b["electric_field_V_m"],
            "baseline_electron_flux": b["electron_flux"],
            "sentqf_electron_flux": s["electron_flux"],
            "delta_electron_flux": d_e,
            "baseline_hole_flux": b["hole_flux"],
            "sentqf_hole_flux": s["hole_flux"],
            "delta_hole_flux": d_h,
            "abs_delta_total": abs(d_e) + abs(d_h),
            "in_depletion": abs(mid_x - junction_x) <= half_width_m,
            "active_high_field": b["electric_field_V_m"] >= args.active_field_threshold,
        })

    # Per-node divergence delta reconstruction (cross-check vs Step D per-node).
    node_div_e: dict[int, float] = {}
    node_div_h: dict[int, float] = {}
    for r in rows:
        node_div_e[r["node0"]] = node_div_e.get(r["node0"], 0.0) + r["delta_electron_flux"]
        node_div_e[r["node1"]] = node_div_e.get(r["node1"], 0.0) - r["delta_electron_flux"]
        node_div_h[r["node0"]] = node_div_h.get(r["node0"], 0.0) + r["delta_hole_flux"]
        node_div_h[r["node1"]] = node_div_h.get(r["node1"], 0.0) - r["delta_hole_flux"]

    depletion_rows = [r for r in rows if r["in_depletion"]]
    active_rows = [r for r in rows if r["active_high_field"]]
    depletion_node_ids = {r["node0"] for r in depletion_rows} | {r["node1"] for r in depletion_rows}

    # Write the per-edge CSV.
    edge_out = args.out_dir / "edge_flux_divergence_edges.csv"
    fieldnames = list(rows[0].keys()) if rows else []
    with edge_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Distribution binned by distance-to-junction (0.05 um bins).
    bins: dict[int, dict[str, float]] = {}
    for r in rows:
        key = int(math.floor(r["dist_to_junction_um"] / 0.05))
        slot = bins.setdefault(key, {"count": 0.0, "sum_de": 0.0, "sum_dh": 0.0,
                                     "sum_abs": 0.0})
        slot["count"] += 1
        slot["sum_de"] += r["delta_electron_flux"]
        slot["sum_dh"] += r["delta_hole_flux"]
        slot["sum_abs"] += r["abs_delta_total"]
    bin_rows = []
    for key in sorted(bins):
        slot = bins[key]
        bin_rows.append({
            "dist_bin_um_low": round(key * 0.05, 4),
            "dist_bin_um_high": round((key + 1) * 0.05, 4),
            "edge_count": int(slot["count"]),
            "sum_delta_electron_flux": slot["sum_de"],
            "sum_delta_hole_flux": slot["sum_dh"],
            "sum_abs_delta_total": slot["sum_abs"],
        })
    bin_out = args.out_dir / "edge_flux_divergence_distance_bins.csv"
    with bin_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(bin_rows[0].keys()))
        writer.writeheader()
        writer.writerows(bin_rows)

    top = sorted(rows, key=lambda r: r["abs_delta_total"], reverse=True)[: args.top_n]

    summary: dict[str, Any] = {
        "junction_x_um": junction_x * 1.0e6,
        "depletion_half_width_um": args.depletion_half_width_um,
        "active_field_threshold_V_m": args.active_field_threshold,
        "edge_count_total": len(rows),
        "edge_count_depletion": len(depletion_rows),
        "edge_count_active_high_field": len(active_rows),
        "all_edges": {
            "delta_electron_flux": summarize([r["delta_electron_flux"] for r in rows]),
            "delta_hole_flux": summarize([r["delta_hole_flux"] for r in rows]),
        },
        "depletion_edges": {
            "delta_electron_flux": summarize([r["delta_electron_flux"] for r in depletion_rows]),
            "delta_hole_flux": summarize([r["delta_hole_flux"] for r in depletion_rows]),
        },
        "active_high_field_edges": {
            "delta_electron_flux": summarize([r["delta_electron_flux"] for r in active_rows]),
            "delta_hole_flux": summarize([r["delta_hole_flux"] for r in active_rows]),
        },
        "per_node_divergence_delta_cross_check": {
            "all_nodes": {
                "electron": summarize(list(node_div_e.values())),
                "hole": summarize(list(node_div_h.values())),
            },
            "depletion_nodes": {
                "electron": summarize([node_div_e[n] for n in depletion_node_ids if n in node_div_e]),
                "hole": summarize([node_div_h[n] for n in depletion_node_ids if n in node_div_h]),
            },
        },
        "concentration": {
            "abs_delta_total_sum_all": math.fsum(r["abs_delta_total"] for r in rows),
            "abs_delta_total_sum_depletion": math.fsum(r["abs_delta_total"] for r in depletion_rows),
            "abs_delta_total_sum_active": math.fsum(r["abs_delta_total"] for r in active_rows),
        },
        "top_edges": [
            {
                "edge_id": r["edge_id"],
                "mid_x_um": round(r["mid_x_um"], 5),
                "dist_to_junction_um": round(r["dist_to_junction_um"], 5),
                "electric_field_V_m": r["electric_field_V_m"],
                "net_doping_avg_m3": r["net_doping_avg_m3"],
                "delta_electron_flux": r["delta_electron_flux"],
                "delta_hole_flux": r["delta_hole_flux"],
                "abs_delta_total": r["abs_delta_total"],
            }
            for r in top
        ],
    }
    depl_abs = summary["concentration"]["abs_delta_total_sum_depletion"]
    all_abs = summary["concentration"]["abs_delta_total_sum_all"]
    summary["concentration"]["depletion_fraction_of_abs_delta"] = (
        depl_abs / all_abs if all_abs else None)
    act_abs = summary["concentration"]["abs_delta_total_sum_active"]
    summary["concentration"]["active_fraction_of_abs_delta"] = (
        act_abs / all_abs if all_abs else None)

    summary_out = args.out_dir / "edge_flux_divergence_summary.json"
    summary_out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({
        "junction_x_um": summary["junction_x_um"],
        "edge_count_total": summary["edge_count_total"],
        "edge_count_depletion": summary["edge_count_depletion"],
        "depletion_fraction_of_abs_delta":
            summary["concentration"]["depletion_fraction_of_abs_delta"],
        "active_fraction_of_abs_delta":
            summary["concentration"]["active_fraction_of_abs_delta"],
        "depletion_node_div_electron_median":
            summary["per_node_divergence_delta_cross_check"]["depletion_nodes"]["electron"]["median"],
        "depletion_node_div_hole_median":
            summary["per_node_divergence_delta_cross_check"]["depletion_nodes"]["hole"]["median"],
        "outputs": {
            "edges_csv": str(edge_out),
            "distance_bins_csv": str(bin_out),
            "summary_json": str(summary_out),
        },
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
