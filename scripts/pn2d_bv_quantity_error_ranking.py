#!/usr/bin/env python3
"""Rank exported physical quantities by their impact on the pn2d BV calculation.

Breakdown voltage (BV) is governed by avalanche multiplication in the junction:
the terminal current runs away when the ionization integral int(alpha dl) -> 1.
This script scans the coarse7x3 aligned comparison CSV and, per quantity per
bias, reports the Vela-vs-Sentaurus magnitude error (in decades = |log10(ratio)|)
restricted to nodes where the quantity is physically significant, so that
near-zero noise does not dominate the ranking.

For each (bias, quantity):
  * consider only nodes where |sentaurus_value| >= sig_frac * max|sentaurus_value|
  * ratio = |vela| / |sentaurus|  (magnitude; BV cares about magnitude)
  * report max/median |log10(ratio)| and the worst node.

The BV-relevant quantities (avalanche generation, ionization coefficients,
ionization integral, field, seed currents) are highlighted first.

Pure standard library.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

DEFAULT_CSV = os.path.join(
    "build-release",
    "reference_tcad",
    "pn2d_sentaurus2018_coarse7x3",
    "reports",
    "coarse_previous_full20_vector_current_20260630",
    "coarse_node_field_compare_aligned.csv",
)

# Ordered by expected BV leverage (exponential/multiplicative first).
BV_QUANTITIES = [
    "electron_alpha_avalanche",
    "hole_alpha_avalanche",
    "avalanche",
    "electron_ion_integral",
    "hole_ion_integral",
    "mean_ion_integral",
    "electric_field_mag",
    "electric_field_x",
    "electron_current_total",
    "hole_current_total",
    "total_current",
    "electron_density",
    "hole_density",
    "srh",
    "electron_qf",
    "hole_qf",
    "potential",
]


def _fnum(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def load(path, biases):
    """Return {bias: {quantity: [(node, x, sen, vela), ...]}}."""
    out = defaultdict(lambda: defaultdict(list))
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = [c.strip() for c in next(reader)]
        idx = {name: i for i, name in enumerate(header)}
        for row in reader:
            if len(row) < len(header):
                continue
            b = _fnum(row[idx["bias_V"]])
            if biases and all(abs(b - t) > 1e-9 for t in biases):
                continue
            q = row[idx["quantity"]].strip()
            node = int(_fnum(row[idx["node_id"]]))
            x = _fnum(row[idx["x_um"]])
            sen = _fnum(row[idx["sentaurus_value"]])
            vela = _fnum(row[idx["vela_value_scaled_to_sentaurus_units"]])
            out[b][q].append((node, x, sen, vela))
    return out


def analyze_quantity(rows, sig_frac, x_lo, x_hi):
    """Return (max_dex, median_dex, worst_node, n_used) over significant nodes."""
    in_band = [r for r in rows if x_lo - 1e-9 <= r[1] <= x_hi + 1e-9]
    if not in_band:
        in_band = rows
    max_abs_sen = max((abs(r[2]) for r in in_band if math.isfinite(r[2])), default=0.0)
    if max_abs_sen <= 0.0:
        return None
    floor = sig_frac * max_abs_sen
    dex = []
    worst = (0.0, None)
    for node, x, sen, vela in in_band:
        if not (math.isfinite(sen) and math.isfinite(vela)):
            continue
        if abs(sen) < floor or sen == 0.0 or vela == 0.0:
            continue
        d = abs(math.log10(abs(vela) / abs(sen)))
        dex.append(d)
        if d > worst[0]:
            worst = (d, node)
    if not dex:
        return None
    dex.sort()
    median = dex[len(dex) // 2]
    return max(dex), median, worst[1], len(dex)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument(
        "--biases",
        type=float,
        nargs="*",
        default=[-5.0, -10.0, -16.0, -18.0, -20.0],
        help="Biases to analyze (default reverse sweep).",
    )
    parser.add_argument("--sig-frac", type=float, default=0.05,
                        help="Significance floor as fraction of peak |sentaurus| per quantity.")
    parser.add_argument("--x-lo", type=float, default=0.75, help="Junction band low x [um].")
    parser.add_argument("--x-hi", type=float, default=1.25, help="Junction band high x [um].")
    args = parser.parse_args()

    data = load(args.csv, args.biases)

    print(f"# BV-impact quantity error ranking  (junction band x in [{args.x_lo},{args.x_hi}] um)")
    print(f"# metric = |log10(|vela|/|sentaurus|)| in DECADES over nodes with "
          f"|sen| >= {args.sig_frac:g}*peak; sorted by max-dex at the most negative bias.")
    print()

    biases = sorted(args.biases)  # most negative first when reversed below
    for b in reversed(biases):
        qmap = data.get(b, {})
        results = []
        for q in BV_QUANTITIES:
            if q not in qmap:
                continue
            res = analyze_quantity(qmap[q], args.sig_frac, args.x_lo, args.x_hi)
            if res is None:
                continue
            results.append((q, *res))
        results.sort(key=lambda r: r[1], reverse=True)

        print("=" * 74)
        print(f"bias = {b:g} V")
        print(f"  {'quantity':<26}{'max_dex':>9}{'med_dex':>9}{'x_ratio':>12}{'worst_node':>11}{'n':>4}")
        print("  " + "-" * 68)
        for q, max_dex, med_dex, node, n in results:
            xr = 10.0 ** max_dex
            print(f"  {q:<26}{max_dex:>9.3f}{med_dex:>9.3f}{xr:>12.3g}{node:>11}{n:>4}")
        print()


if __name__ == "__main__":
    main()
