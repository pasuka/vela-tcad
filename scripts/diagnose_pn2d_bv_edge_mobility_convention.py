#!/usr/bin/env python3
"""Direction 1a audit: is the edge-mobility AVERAGING CONVENTION a lever for BV?

Step4 reported Vela/Sentaurus active-edge mobility ~= 0.983 (e)/0.984 (h) and the
proposal's Direction 1a asked whether detail::edgeMobility's averaging
(arithmetic / harmonic / cell-based) is inconsistent with Sentaurus.

Two code facts settle the convention question:
  * The 0.983 metric (scripts/diagnose_pn2d_bv_sg_flux_forms.py:edge_row ->
    avg(state["mun"], i, j)) is the ENDPOINT-AVERAGE of each engine's own
    exported NODAL mobility, applied IDENTICALLY to Vela and Sentaurus. Same
    averaging on both sides => the 2% is a converged-STATE difference (field /
    density feeding the same Masetti+field model), not a convention artifact.
  * Vela's solver edgeMobility evaluates the model at the arithmetic edge-
    averaged net doping and arithmetic-averages over adjacent transport cells.
    For a UNIFORM-doping single-material edge, evaluate-at-avg-doping,
    endpoint-average-of-nodal, and harmonic all COINCIDE. The conventions can
    only diverge where N0 and N1 differ strongly (junction-straddling) or across
    a material interface.

This read-only check joins the Step D per-node net doping with the Step F
per-edge SG flux-delta budget and asks: do the edges that actually carry the
flux-delta sit in uniform doping (=> averaging convention is irrelevant to BV)?
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stepd-nodes", type=Path, required=True,
                   help="Step D exact_carrier_term_state_nodes.csv")
    p.add_argument("--stepf-edges", type=Path, required=True,
                   help="Step F edge_flux_divergence_edges.csv")
    p.add_argument("--variant", default="vela_baseline")
    p.add_argument("--spread-threshold", type=float, default=0.10,
                   help="Relative |N0-N1| above which conventions diverge.")
    p.add_argument("--top-n", type=int, default=25)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    node_doping: dict[int, float] = {}
    with args.stepd_nodes.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("variant") != args.variant:
                continue
            node_doping[int(row["node_id"])] = float(row["net_doping_m3"])

    edges: list[dict[str, float]] = []
    with args.stepf_edges.open(newline="") as fh:
        for row in csv.DictReader(fh):
            n0 = int(row["node0"])
            n1 = int(row["node1"])
            if n0 not in node_doping or n1 not in node_doping:
                continue
            d0 = node_doping[n0]
            d1 = node_doping[n1]
            denom = max(abs(d0), abs(d1), 1.0)
            spread = abs(d0 - d1) / denom
            edges.append({
                "edge_id": int(row["edge_id"]),
                "mid_x_um": float(row["mid_x_um"]),
                "field": float(row["electric_field_V_m"]),
                "n0_doping": d0,
                "n1_doping": d1,
                "doping_rel_spread": spread,
                "abs_delta_total": float(row["abs_delta_total"]),
            })

    total_abs = math.fsum(e["abs_delta_total"] for e in edges)
    straddling = [e for e in edges if e["doping_rel_spread"] > args.spread_threshold]
    uniform = [e for e in edges if e["doping_rel_spread"] <= args.spread_threshold]
    straddle_abs = math.fsum(e["abs_delta_total"] for e in straddling)
    uniform_abs = math.fsum(e["abs_delta_total"] for e in uniform)

    top = sorted(edges, key=lambda e: e["abs_delta_total"], reverse=True)[: args.top_n]
    top_spread = [e["doping_rel_spread"] for e in top]

    print("edges_joined:", len(edges))
    print("total_abs_delta:", f"{total_abs:.6e}")
    print(f"straddling_edges (rel_spread>{args.spread_threshold}):", len(straddling))
    print("  abs_delta carried by straddling edges:",
          f"{straddle_abs:.6e}",
          f"({(straddle_abs/total_abs if total_abs else 0):.4%})")
    print("uniform_edges:", len(uniform),
          " abs_delta fraction:",
          f"{(uniform_abs/total_abs if total_abs else 0):.4%}")
    print()
    print(f"top-{args.top_n} flux-delta edges doping_rel_spread:")
    print("  max:", f"{max(top_spread):.3e}", " mean:", f"{sum(top_spread)/len(top_spread):.3e}")
    print()
    print("top edges (x_um | field | rel_spread | abs_delta):")
    for e in top[:12]:
        print(f"  x={e['mid_x_um']:.4f}  E={e['field']:.2e}  "
              f"spread={e['doping_rel_spread']:.2e}  abs_delta={e['abs_delta_total']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
