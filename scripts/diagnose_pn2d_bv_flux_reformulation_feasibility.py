#!/usr/bin/env python3
"""Direction 2 FEASIBILITY: can an SG flux-divergence reformulation close 0.73x?

The only way a continuity-flux reformulation closes the BV current deficit is if,
AT THE SENTAURUS STATE, it drives Vela's per-node continuity residual to ~0
(then Vela's converged fixed point == Sentaurus, raising the minority QF split
and the current). Step D established the residual at the Sentaurus state is
carried by the SG flux-divergence term.

A reformulation that multiplies the per-node flux-divergence term by (1+s) makes
the residual:  (1+s)*flux_i + recomb_i + impact_i + gauge_i + boundary_i.
Setting that to zero gives the REQUIRED per-node scale:
    s_i = -(flux_i + recomb_i + impact_i + ...) / flux_i = -residual_i / flux_i.

FEASIBILITY VERDICT:
  * If s_i is tightly clustered around a single value across the depletion nodes,
    a principled GLOBAL flux-coefficient reformulation could close the gap
    (low risk -> worth implementing behind a flag).
  * If s_i varies over orders of magnitude AND/OR flips sign, then NO consistent
    (global or smooth-local) reformulation can zero the residual without per-node
    tuning -> Direction 2 is ill-conditioned (high risk / low payoff).

Also reports the cancellation/conditioning factor: per-node incident |edge flux|
sum divided by the net divergence, i.e. how large the fluxes being perturbed are
relative to the net the reformulation must control.

Read-only. No solver changes.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stepd-nodes", type=Path, required=True,
                   help="Step D exact_carrier_term_state_nodes.csv")
    p.add_argument("--sentqf-edges", type=Path, required=True,
                   help="Step F vela_psi_sentaurus_qf_edge_flux.csv")
    p.add_argument("--variant", default="vela_psi_sentaurus_qf")
    p.add_argument("--support-classes", default="active",
                   help="Comma list of support_class values to include.")
    p.add_argument("--min-abs-flux", type=float, default=1.0e-30,
                   help="Skip nodes whose |flux term| is below this (avoid 0/0).")
    return p.parse_args()


def quantiles(values: list[float]) -> dict[str, float]:
    if not values:
        return {}
    s = sorted(values)
    n = len(s)

    def q(p: float) -> float:
        if n == 1:
            return s[0]
        idx = p * (n - 1)
        lo = int(math.floor(idx))
        hi = min(lo + 1, n - 1)
        return s[lo] + (s[hi] - s[lo]) * (idx - lo)

    return {
        "count": n,
        "min": s[0], "p05": q(0.05), "p25": q(0.25), "median": q(0.50),
        "p75": q(0.75), "p95": q(0.95), "max": s[-1],
        "mean": statistics.fmean(s),
    }


def main() -> int:
    args = parse_args()
    keep = {c.strip() for c in args.support_classes.split(",") if c.strip()}

    # --- per-node residual + flux-term (Step D) for the Sentaurus-QF state ----
    nodes: dict[int, dict[str, float]] = {}
    with args.stepd_nodes.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("variant") != args.variant:
                continue
            if keep and row.get("support_class") not in keep:
                continue
            nid = int(row["node_id"])
            nodes[nid] = {
                "x": float(row["x"]),
                "e_flux": float(row["electron_flux"]),
                "e_resid": float(row["electron_residual"]),
                "h_flux": float(row["hole_flux"]),
                "h_resid": float(row["hole_residual"]),
            }

    # --- per-edge fluxes (Step F sentaurus-qf) for the conditioning factor -----
    incident_abs_e: dict[int, float] = {}
    incident_abs_h: dict[int, float] = {}
    with args.sentqf_edges.open(newline="") as fh:
        for row in csv.DictReader(fh):
            n0 = int(row["node0"])
            n1 = int(row["node1"])
            fe = abs(float(row["electron_flux"]))
            fh_ = abs(float(row["hole_flux"]))
            for n in (n0, n1):
                incident_abs_e[n] = incident_abs_e.get(n, 0.0) + fe
                incident_abs_h[n] = incident_abs_h.get(n, 0.0) + fh_

    req_s_e: list[float] = []
    req_s_h: list[float] = []
    cond_e: list[float] = []
    cond_h: list[float] = []
    sign_flip_e = 0
    sign_flip_h = 0
    for nid, rec in nodes.items():
        if abs(rec["e_flux"]) > args.min_abs_flux:
            s = -rec["e_resid"] / rec["e_flux"]
            req_s_e.append(s)
            if s < 0:
                sign_flip_e += 1
            inc = incident_abs_e.get(nid)
            if inc and abs(rec["e_flux"]) > 0:
                cond_e.append(inc / abs(rec["e_flux"]))
        if abs(rec["h_flux"]) > args.min_abs_flux:
            s = -rec["h_resid"] / rec["h_flux"]
            req_s_h.append(s)
            if s < 0:
                sign_flip_h += 1
            inc = incident_abs_h.get(nid)
            if inc and abs(rec["h_flux"]) > 0:
                cond_h.append(inc / abs(rec["h_flux"]))

    def show(label: str, vals: list[float]) -> None:
        q = quantiles(vals)
        if not q:
            print(f"{label}: (no data)")
            return
        print(f"{label}: n={q['count']}")
        print(f"   median={q['median']:.4e}  mean={q['mean']:.4e}")
        print(f"   p05={q['p05']:.4e}  p25={q['p25']:.4e}  p75={q['p75']:.4e}  p95={q['p95']:.4e}")
        print(f"   min={q['min']:.4e}  max={q['max']:.4e}")

    print("=== Required per-node flux-divergence scale s_i = -residual_i/flux_i ===")
    print("(A consistent reformulation needs s_i clustered & single-signed.)")
    show("electron s_i", req_s_e)
    print(f"   sign-flips (s_i<0): {sign_flip_e}/{len(req_s_e)}")
    show("hole s_i", req_s_h)
    print(f"   sign-flips (s_i<0): {sign_flip_h}/{len(req_s_h)}")
    # spread metric: p95/p05 magnitude ratio and IQR/median
    for label, vals in (("electron", req_s_e), ("hole", req_s_h)):
        q = quantiles(vals)
        if not q:
            continue
        med = q["median"]
        iqr = q["p75"] - q["p25"]
        rng = q["max"] - q["min"]
        print(f"   {label} spread: IQR/|median|={iqr/abs(med) if med else float('inf'):.2f}  "
              f"range/|median|={rng/abs(med) if med else float('inf'):.2f}")

    print()
    print("=== Conditioning factor (incident |edge flux| sum)/|node flux term| ===")
    print("(How large the perturbed fluxes are vs the net the reformulation controls.)")
    show("electron conditioning", cond_e)
    show("hole conditioning", cond_h)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
