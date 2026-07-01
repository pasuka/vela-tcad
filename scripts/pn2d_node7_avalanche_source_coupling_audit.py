#!/usr/bin/env python3
"""Audit avalanche electron-source coupling at the upstream p-side node7.

Candidate (3): does the SG edge avalanche generation actually SEED the electron
continuity at node7 (x=0.75, p-side depletion edge), or is it routed away so
node7 is starved -> minority-electron extinction (per probe B)?

Reads a Vela sg_avalanche_edges.csv dump (build/qf_cap_exp/baseline). For the
requested biases, prints every edge incident to the audited nodes with the
per-edge electron alpha, impact field, electron flux proxy, electron source
integral, and how the edge source is distributed to its two end nodes.
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
    ap.add_argument("--biases", type=float, nargs="*", default=[-16.0, -18.0, -20.0])
    ap.add_argument("--nodes", type=int, nargs="*", default=[7, 10])
    args = ap.parse_args()
    data = load(args.csv)

    audit = set(args.nodes)
    for tb in args.biases:
        b = nearest(data, tb)
        if b is None:
            print(f"\n### bias {tb}: not found")
            continue
        print(f"\n### bias {b} V  (edges incident to nodes {sorted(audit)}) ###")
        print(f"{'edge(n0,n1)':>12} {'x0':>5}{'x1':>5} {'edge_src':>11}"
              f" | {'node0_src':>11}{'node1_src':>11}  (directionally weighted -> each end)")
        node_recv = defaultdict(float)   # directionally-weighted source actually assigned to node
        node_edge_count = defaultdict(int)
        for e in data[b]:
            n0 = int(e["node0"]); n1 = int(e["node1"])
            if n0 not in audit and n1 not in audit:
                continue
            tsrc = float(e["edge_source_integral"])
            n0src = float(e["node0_source_integral"])
            n1src = float(e["node1_source_integral"])
            x0 = float(e["x0_um"]); x1 = float(e["x1_um"])
            if n0 in audit:
                node_recv[n0] += n0src
                node_edge_count[n0] += 1
            if n1 in audit:
                node_recv[n1] += n1src
                node_edge_count[n1] += 1
            mark0 = "*" if n0 in audit else " "
            mark1 = "*" if n1 in audit else " "
            print(f"  ({n0:>2},{n1:>2})  {x0:>5.2f}{x1:>5.2f} {tsrc:>11.3g}"
                  f" | {mark0}{n0src:>10.3g}{mark1}{n1src:>10.3g}")
        print("  -- per-node RECEIVED source (directionally-weighted node0/node1_source_integral):")
        for n in sorted(audit):
            print(f"     node {n}: received = {node_recv[n]:.4g}  ({node_edge_count[n]} edges)")


if __name__ == "__main__":
    main()
