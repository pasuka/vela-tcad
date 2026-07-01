#!/usr/bin/env python3
"""Compare junction-cutline ionization integral between two sg_avalanche_edges dumps.

For each requested bias, integrate the electron avalanche coefficient along the
junction cutline (horizontal interior-row edges at y=Y_ROW, x in [X_LO,X_HI]):
    I = sum_edges  electron_alpha_m_inv * |dx|_m
and also report the peak electron_alpha and the total assembled electron avalanche
source  sum edge_source_integral (electron part) over ALL edges.

Lets us test whether junction MESH REFINEMENT increases int(alpha dl) at a fixed,
commonly-converged bias (mesh-robust, independent of breakdown continuation).
Std lib only.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict


def load(path):
    """Return {bias: list of edge dicts}."""
    out = defaultdict(list)
    with open(path, newline="") as h:
        r = csv.DictReader(h)
        for row in r:
            out[round(float(row["bias_V"]), 4)].append(row)
    return out


def cutline_integral(edges, y_row, x_lo, x_hi, tol=1e-3):
    total = 0.0
    peak = 0.0
    for e in edges:
        y0 = float(e["y0_um"]); y1 = float(e["y1_um"])
        x0 = float(e["x0_um"]); x1 = float(e["x1_um"])
        a = float(e["electron_alpha_m_inv"])
        peak = max(peak, a)
        if abs(y0 - y_row) > tol or abs(y1 - y_row) > tol:
            continue
        if abs(y1 - y0) > tol:  # must be horizontal
            continue
        xmid = 0.5 * (x0 + x1)
        if not (x_lo - tol <= xmid <= x_hi + tol):
            continue
        dx_m = abs(x1 - x0) * 1e-6
        total += a * dx_m
    return total, peak


def total_source(edges):
    tot = 0.0
    for e in edges:
        v = e.get("electron_source_integral", "")
        try:
            tot += abs(float(v))
        except ValueError:
            pass
    return tot


def nearest_bias(biasmap, target, tol=0.03):
    best = None
    for b in biasmap:
        if abs(b - target) <= tol and (best is None or abs(b - target) < abs(best - target)):
            best = b
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--refined", required=True)
    ap.add_argument("--biases", type=float, nargs="*", default=[-3.0, -4.0, -5.0])
    ap.add_argument("--y-row", type=float, default=0.25)
    ap.add_argument("--x-lo", type=float, default=0.5)
    ap.add_argument("--x-hi", type=float, default=1.5)
    args = ap.parse_args()

    base = load(args.baseline)
    ref = load(args.refined)

    print(f"# junction-cutline electron ionization integral  (y={args.y_row}, "
          f"x in [{args.x_lo},{args.x_hi}] um)")
    print(f"# {'bias':>6} | {'I_base':>11} {'peak_base':>11} {'src_base':>11}"
          f" | {'I_ref':>11} {'peak_ref':>11} {'src_ref':>11} | {'I_ref/I_base':>12}")
    for tb in args.biases:
        bb = nearest_bias(base, tb)
        rb = nearest_bias(ref, tb)
        if bb is None or rb is None:
            print(f"  {tb:>6.2f} | missing (base={bb}, ref={rb})")
            continue
        Ib, pkb = cutline_integral(base[bb], args.y_row, args.x_lo, args.x_hi)
        Ir, pkr = cutline_integral(ref[rb], args.y_row, args.x_lo, args.x_hi)
        sb = total_source(base[bb])
        sr = total_source(ref[rb])
        ratio = Ir / Ib if Ib > 0 else float("nan")
        print(f"  {tb:>6.2f} | {Ib:>11.4g} {pkb:>11.4g} {sb:>11.4g}"
              f" | {Ir:>11.4g} {pkr:>11.4g} {sr:>11.4g} | {ratio:>12.3f}"
              f"   (base@{bb} ref@{rb})")


if __name__ == "__main__":
    main()
