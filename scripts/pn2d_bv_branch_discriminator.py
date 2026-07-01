#!/usr/bin/env python3
"""Branch discrimination: Vela vs Sentaurus junction profile at -20V.

Vela lands on the NON-multiplying branch (p-side node7 dead, n-side floods, flat
terminal current); Sentaurus on the MULTIPLYING branch (n7=8870). This pulls the
interior-row (y=0.25) profile x=0.5..1.5 (nodes 4,7,10,13,16) of electron/hole
density, electron/hole current x-component (signed, direction), and avalanche
generation for BOTH codes, to locate where the branches diverge (e.g. is the
avalanche current fed to BOTH sides or only the n-side?).

Steady-state current continuity: Jn_x + Jp_x = const (= terminal J). In avalanche
Jn_x grows toward cathode (electrons), Jp_x grows toward anode (holes). Reads the
aligned comparison CSV. Std lib only.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict

DEFAULT_CSV = os.path.join(
    "build-release", "reference_tcad", "pn2d_sentaurus2018_coarse7x3",
    "reports", "coarse_previous_full20_vector_current_20260630",
    "coarse_node_field_compare_aligned.csv",
)
# interior row y=0.25 : x = 0.5,0.75,1.0,1.25,1.5  (also anode-side 1 / cathode-side 19,22)
ROW = [(1, 0.25), (4, 0.5), (7, 0.75), (10, 1.0), (13, 1.25), (16, 1.5), (19, 1.75), (22, 2.0)]
QUANTS = ("electron_density", "hole_density",
          "electron_current_density_x", "hole_current_density_x", "avalanche", "srh")


def _f(t):
    try:
        return float(t)
    except (TypeError, ValueError):
        return float("nan")


def parse(path, bias):
    data = defaultdict(dict)  # node -> quant -> (sen, vela)
    with open(path, newline="") as h:
        r = csv.reader(h)
        header = [c.strip() for c in next(r)]
        idx = {n: i for i, n in enumerate(header)}
        for row in r:
            if len(row) < len(header):
                continue
            if abs(_f(row[idx["bias_V"]]) - bias) > 0.03:
                continue
            q = row[idx["quantity"]].strip()
            if q not in QUANTS:
                continue
            node = int(_f(row[idx["node_id"]]))
            data[node][q] = (_f(row[idx["sentaurus_value"]]),
                             _f(row[idx["vela_value_scaled_to_sentaurus_units"]]))
    return data


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--bias", type=float, default=-20.0)
    args = ap.parse_args()
    data = parse(args.csv, args.bias)

    def col(q, who):
        return f"{q}_{who}"

    print(f"# Junction interior-row profile @ {args.bias}V  (Sen | Vela)")
    for q in QUANTS:
        print(f"\n## {q}")
        print(f"  {'x_um':>5} {'node':>4} | {'Sentaurus':>13} {'Vela':>13} {'Vela/Sen':>10}")
        for node, x in ROW:
            e = data.get(node, {}).get(q)
            if not e:
                continue
            sen, vela = e
            ratio = vela / sen if sen not in (0.0,) and sen == sen and sen != 0 else float("nan")
            print(f"  {x:>5.2f} {node:>4} | {sen:>13.4g} {vela:>13.4g} {ratio:>10.3g}")


if __name__ == "__main__":
    main()
