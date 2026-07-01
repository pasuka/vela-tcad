#!/usr/bin/env python3
"""Residual reconciliation at node7: does the electron continuity actually balance?

Paradox: node7 electron continuity receives a huge avalanche source (~2.5e9 @-20V,
mostly hole-driven pair generation) yet converged n7~=0 (probe B). This checks the
converged electron-continuity balance directly from the sg_avalanche_edges dump.

CoupledDDAssembler electron residual convention (src/equation/CoupledDDAssembler.cpp):
  per edge (node0,node1): nFlux = signedFlux01 * couple   (coef carries the couple)
    r(phin, node0) += nFlux ;  r(phin, node1) -= nFlux
  plus  r(phin, i) -= source_i   (avalanche, L377, = combined node0/1_source_integral)
At convergence r=0  =>  sum_edges(+/- nFlux) == source(node7)  (ignoring SRH, added sep.)

The dump's electron_raw_signed_flux_proxy = signedFlux01 computed with coef=mun*Vt/h
(NO couple), so the residual flux = electron_raw_signed_flux_proxy * edge_couple_m.
Sum the oriented electron outflow at the node and compare to the applied source.
Std lib only.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict


def load(path):
    out = defaultdict(list)
    with open(path, newline="") as h:
        r = csv.DictReader(h)
        for row in r:
            out[round(float(row["bias_V"]), 4)].append(row)
    return out


def nearest(biasmap, target, tol=0.03):
    best = None
    for b in biasmap:
        if abs(b - target) <= tol and (best is None or abs(b - target) < abs(best - target)):
            best = b
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", default="build/qf_cap_exp/baseline/sg_avalanche_edges.csv")
    ap.add_argument("--biases", type=float, nargs="*", default=[-18.0, -20.0])
    ap.add_argument("--nodes", type=int, nargs="*", default=[7, 10])
    args = ap.parse_args()
    data = load(args.csv)

    for tb in args.biases:
        b = nearest(data, tb)
        if b is None:
            continue
        print(f"\n### bias {b} V ###")
        for node in args.nodes:
            e_out = 0.0     # signed electron outflow (sum of residual nFlux contributions)
            h_out = 0.0     # signed hole outflow
            src_combined = 0.0
            edges = []
            for e in data[b]:
                n0 = int(e["node0"]); n1 = int(e["node1"])
                if node != n0 and node != n1:
                    continue
                couple = float(e["edge_couple_m"])
                e_signed = float(e["electron_raw_signed_flux_proxy"]) * couple
                h_signed = float(e["hole_raw_signed_flux_proxy"]) * couple
                n0src = float(e["node0_source_integral"])
                n1src = float(e["node1_source_integral"])
                if node == n0:
                    e_out += e_signed
                    h_out += h_signed
                    src_combined += n0src
                    other = n1
                else:
                    e_out += -e_signed
                    h_out += -h_signed
                    src_combined += n1src
                    other = n0
                edges.append((node, other,
                              e_signed if node == n0 else -e_signed,
                              n0src if node == n0 else n1src))
            print(f"  node {node}:")
            print(f"    {'->neighbor':>10} {'e_outflow(resid)':>18} {'src_to_node':>14}")
            for (a, o, ef, s) in edges:
                print(f"    {o:>10} {ef:>18.4g} {s:>14.4g}")
            print(f"    SUM electron outflow   = {e_out:.6g}")
            print(f"    applied avalanche src  = {src_combined:.6g}   (combined e+h, -> electron continuity)")
            resid = e_out - src_combined
            rel = resid / src_combined if src_combined != 0 else float("nan")
            print(f"    electron residual (outflow - src) = {resid:.4g}   (rel {rel:.3g})")
            print(f"    [hole outflow = {h_out:.4g}; hole residual = {h_out - src_combined:.4g}]")


if __name__ == "__main__":
    main()
