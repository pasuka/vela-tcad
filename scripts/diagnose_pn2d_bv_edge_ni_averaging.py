#!/usr/bin/env python3
"""Direction 1b audit: is the VariableNi edge ni handling a lever for BV?

The SG VariableNi continuity flux (src/discretization/ScharfetterGummel.cpp
sgElectronContinuityFluxFromQuasiFermiVariableNi) uses each node's OWN ni:
    eta = (psi1-psi0)/Vt + log(ni1/ni0)            # ni-gradient drift
    n0  = ni0 * exp((psi0-phin0)/Vt)
    n1  = ni1 * exp((psi1-phin1)/Vt)
    flux = coef*(B(-eta)*n0 - B(eta)*n1)
Vt is a single isothermal scalar => edge-averaging Vt is a no-op. So Direction
1b reduces to: replace nodal ni0/ni1 with an edge-averaged ni (geometric mean).

On a UNIFORM-doping edge ni0 == ni1, so log(ni1/ni0)=0 and geomean == ni0 == ni1
=> the nodal-ni and edge-averaged-ni forms are IDENTICAL. ni0 and ni1 differ
only where BGN ni_eff (hence doping) varies between the two nodes.

This read-only check uses the per-edge ni0_m3/ni1_m3 already emitted by the
sg_edge_flux_probe (Step F baseline CSV), joined with the Step F per-edge SG
flux-delta budget, and asks: do the edges that carry the flux-delta have
ni0 != ni1 (so that edge-averaging ni could change anything)?
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--baseline-edges", type=Path, required=True,
                   help="Step F <variant>_edge_flux.csv (has ni0_m3, ni1_m3).")
    p.add_argument("--divergence-edges", type=Path, required=True,
                   help="Step F edge_flux_divergence_edges.csv (has abs_delta_total).")
    p.add_argument("--spread-threshold", type=float, default=1.0e-3,
                   help="Relative |ni0-ni1| above which the ni form matters.")
    p.add_argument("--top-n", type=int, default=25)
    return p.parse_args()


def main() -> int:
    args = parse_args()

    ni_by_edge: dict[int, tuple[float, float]] = {}
    with args.baseline_edges.open(newline="") as fh:
        for row in csv.DictReader(fh):
            ni_by_edge[int(row["edge_id"])] = (
                float(row["ni0_m3"]), float(row["ni1_m3"]))

    edges: list[dict[str, float]] = []
    with args.divergence_edges.open(newline="") as fh:
        for row in csv.DictReader(fh):
            eid = int(row["edge_id"])
            if eid not in ni_by_edge:
                continue
            ni0, ni1 = ni_by_edge[eid]
            denom = max(ni0, ni1, 1.0)
            spread = abs(ni0 - ni1) / denom
            edges.append({
                "edge_id": eid,
                "mid_x_um": float(row["mid_x_um"]),
                "field": float(row["electric_field_V_m"]),
                "ni0": ni0,
                "ni1": ni1,
                "ni_rel_spread": spread,
                "abs_delta_total": float(row["abs_delta_total"]),
            })

    total_abs = math.fsum(e["abs_delta_total"] for e in edges)
    graded = [e for e in edges if e["ni_rel_spread"] > args.spread_threshold]
    uniform = [e for e in edges if e["ni_rel_spread"] <= args.spread_threshold]
    graded_abs = math.fsum(e["abs_delta_total"] for e in graded)
    uniform_abs = math.fsum(e["abs_delta_total"] for e in uniform)

    top = sorted(edges, key=lambda e: e["abs_delta_total"], reverse=True)[: args.top_n]
    top_spread = [e["ni_rel_spread"] for e in top]

    ni_values = sorted({round(e["ni0"], 6) for e in edges} | {round(e["ni1"], 6) for e in edges})

    print("edges_joined:", len(edges))
    print("distinct ni values across all edge endpoints:", len(ni_values))
    print("  ni range:", f"{min(ni_values):.4e}", "..", f"{max(ni_values):.4e}")
    print("total_abs_delta:", f"{total_abs:.6e}")
    print(f"graded-ni edges (rel_spread>{args.spread_threshold}):", len(graded))
    print("  abs_delta carried by graded-ni edges:",
          f"{graded_abs:.6e}",
          f"({(graded_abs/total_abs if total_abs else 0):.4%})")
    print("uniform-ni edges:", len(uniform),
          " abs_delta fraction:",
          f"{(uniform_abs/total_abs if total_abs else 0):.4%}")
    print()
    print(f"top-{args.top_n} flux-delta edges ni_rel_spread:")
    print("  max:", f"{max(top_spread):.3e}", " mean:", f"{sum(top_spread)/len(top_spread):.3e}")
    print()
    print("top edges (x_um | field | ni_rel_spread | abs_delta):")
    for e in top[:12]:
        print(f"  x={e['mid_x_um']:.4f}  E={e['field']:.2e}  "
              f"ni_spread={e['ni_rel_spread']:.2e}  abs_delta={e['abs_delta_total']:.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
