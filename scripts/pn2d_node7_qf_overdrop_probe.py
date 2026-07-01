#!/usr/bin/env python3
"""Node-7 minority quasi-Fermi over-drop mechanism probe (coarse7x3 BV).

(i)/edge-operator analysis localized the -20V ionization-integral deficit to the
upstream p-side junction edge (7,10): Vela's phin gradient there is ~7% low.
Entry 74 shows node7 (x=0.75, p-side deep-depletion edge) has Vela phin pinned
HIGH (dphin ~ +0.76 V @-20V -> electrons -> 0). This probe tests WHY: it tracks,
across bias, for the upstream path nodes (4 x=0.5, 7 x=0.75, 10 x=1.0), the
electron quasi-Fermi level, electron density, and electron current magnitude
(Vela vs Sentaurus). Hypothesis: node7 phin floats to the anode value because the
minority electron current/generation seed is far weaker in Vela than Sentaurus
(few/no avalanche electrons transported through node7), leaving phin weakly
determined and pinned high -> flat phin gradient on edge(7,10).

Reads the aligned comparison CSV. Std lib only. phin/psi in V, Vt=25.852 mV.
"""

from __future__ import annotations

import argparse
import csv
import math
import os

K_B_OVER_Q = 8.617333262145e-5
DEFAULT_CSV = os.path.join(
    "build-release", "reference_tcad", "pn2d_sentaurus2018_coarse7x3",
    "reports", "coarse_previous_full20_vector_current_20260630",
    "coarse_node_field_compare_aligned.csv",
)
NODES = [4, 7, 10]
QUANTS = ("electron_qf", "electron_density", "electron_current",
          "electron_alpha_avalanche", "potential", "hole_qf")


def _f(t):
    try:
        return float(t)
    except (TypeError, ValueError):
        return float("nan")


def parse(path):
    # data[bias][node][quant] = (sen, vela)
    data = {}
    with open(path, newline="") as h:
        r = csv.reader(h)
        header = [c.strip() for c in next(r)]
        idx = {n: i for i, n in enumerate(header)}
        for row in r:
            if len(row) < len(header):
                continue
            q = row[idx["quantity"]].strip()
            if q not in QUANTS:
                continue
            node = int(_f(row[idx["node_id"]]))
            if node not in NODES:
                continue
            bias = round(_f(row[idx["bias_V"]]), 3)
            sen = _f(row[idx["sentaurus_value"]])
            vela = _f(row[idx["vela_value_scaled_to_sentaurus_units"]])
            data.setdefault(bias, {}).setdefault(node, {})[q] = (sen, vela)
    return data


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv", default=DEFAULT_CSV)
    ap.add_argument("--biases", type=float, nargs="*",
                    default=[-5.0, -10.0, -16.0, -20.0])
    args = ap.parse_args()
    data = parse(args.csv)

    for node in NODES:
        print(f"\n=== node {node} (x={0.25 + 0.125 * node:.3f} approx) upstream path node ===")
        print(f"# {'bias':>6} | {'phin_sen':>10}{'phin_vela':>10}{'dphin_mV':>10}"
              f" | {'n_sen':>10}{'n_vela':>10}{'n_v/n_s':>10}"
              f" | {'Jn_sen':>11}{'Jn_vela':>11}{'Jn_v/s':>9}")
        for tb in args.biases:
            b = round(tb, 3)
            nd = data.get(b, {}).get(node)
            if not nd:
                continue
            qs, qv = nd.get("electron_qf", (float("nan"),) * 2)
            ns, nv = nd.get("electron_density", (float("nan"),) * 2)
            js, jv = nd.get("electron_current", (float("nan"),) * 2)
            dphin_mV = (qv - qs) * 1e3
            nratio = nv / ns if ns not in (0.0, float("nan")) and math.isfinite(ns) and ns != 0 else float("nan")
            jratio = jv / js if math.isfinite(js) and js != 0 else float("nan")
            print(f"  {b:>6.1f} | {qs:>10.4f}{qv:>10.4f}{dphin_mV:>10.2f}"
                  f" | {ns:>10.3g}{nv:>10.3g}{nratio:>10.3g}"
                  f" | {js:>11.3g}{jv:>11.3g}{jratio:>9.3g}")

    # edge(7,10) phin drop comparison: Sentaurus vs Vela
    print("\n=== edge(7,10) phin drop  phin(7)-phin(10)  (drives upstream field) ===")
    print(f"# {'bias':>6} | {'drop_sen_V':>11}{'drop_vela_V':>12}{'vela/sen':>10}")
    for tb in args.biases:
        b = round(tb, 3)
        n7 = data.get(b, {}).get(7, {}).get("electron_qf")
        n10 = data.get(b, {}).get(10, {}).get("electron_qf")
        if not n7 or not n10:
            continue
        drop_s = n7[0] - n10[0]
        drop_v = n7[1] - n10[1]
        ratio = drop_v / drop_s if drop_s != 0 else float("nan")
        print(f"  {b:>6.1f} | {drop_s:>11.4f}{drop_v:>12.4f}{ratio:>10.3f}")


if __name__ == "__main__":
    main()
