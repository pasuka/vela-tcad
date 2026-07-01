#!/usr/bin/env python3
"""Test whether finer junction-field sampling grows Vela's local ionization integral.

The multiplication is governed by the ionization integral  I = int alpha_n dx
along the current path. On the coarse7x3 mesh Vela evaluates alpha pointwise from
the local |grad phi| (== |grad psi| off-peak, within 2%), which is too peaked.

This experiment reconstructs the field along the junction axis (interior row
y=0.25) from the coarse potential samples at TWO resolutions and applies the
SAME shared van Overstraeten LOCAL model, then compares the resulting ionization
integral against Sentaurus's exported-alpha integral (the target):

  * COARSE   : nodal trapezoid of the local-model alpha at the coarse nodes.
  * FINE-LIN : field = |d psi/dx| from LINEAR psi interpolation (piecewise
               constant per segment) -> lower bound; adds no new peak.
  * FINE-PCHIP: field from a monotone cubic (Fritsch-Carlson) psi(x) -> sharpens
               the junction peak -> optimistic UPPER bound for a LOCAL model.
  * SENTAURUS: nodal trapezoid of Sentaurus's exported alpha (target).

If even FINE-PCHIP (peak-sharpened local model) falls short of SENTAURUS, then
NO amount of local mesh refinement can close the gap -> the deficit is nonlocal.

van Overstraeten electron single branch (== Sentaurus .par):
  alpha[1/m] = 7.03e7 * exp(-1.231e6 / E[V/cm]),  gamma=1.
Std lib only. x [um], psi [V], alpha [1/m], integral dimensionless.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

A_E_PER_M = 7.03e7
B_E_V_PER_CM = 1.231e6
GRAD_V_PER_UM_TO_V_PER_CM = 1.0e4
UM_TO_M = 1.0e-6

DEFAULT_CSV = os.path.join(
    "build-release", "reference_tcad", "pn2d_sentaurus2018_coarse7x3",
    "reports", "coarse_previous_full20_vector_current_20260630",
    "coarse_node_field_compare_aligned.csv",
)


def _fnum(t):
    try:
        return float(t)
    except (TypeError, ValueError):
        return float("nan")


def parse(path, want_row_y):
    """Return {bias: sorted list of (x, psi_s, alpha_s, alpha_v)} for interior row."""
    rows = defaultdict(dict)
    with open(path, newline="") as h:
        reader = csv.reader(h)
        header = [c.strip() for c in next(reader)]
        idx = {n: i for i, n in enumerate(header)}
        for row in reader:
            if len(row) < len(header):
                continue
            q = row[idx["quantity"]].strip()
            if q not in ("potential", "electron_alpha_avalanche"):
                continue
            y = _fnum(row[idx["y_um"]])
            if abs(y - want_row_y) > 1e-6:
                continue
            bias = _fnum(row[idx["bias_V"]])
            x = _fnum(row[idx["x_um"]])
            e = rows[bias].setdefault(x, {})
            if q == "potential":
                e["psi_s"] = _fnum(row[idx["sentaurus_value"]])
            else:
                e["alpha_s"] = _fnum(row[idx["sentaurus_value"]])
                e["alpha_v"] = _fnum(row[idx["vela_value_scaled_to_sentaurus_units"]])
    out = {}
    for bias, xmap in rows.items():
        pts = sorted((x, d.get("psi_s", float("nan")),
                      d.get("alpha_s", float("nan")), d.get("alpha_v", float("nan")))
                     for x, d in xmap.items())
        out[bias] = pts
    return out


def local_alpha(field_v_per_cm):
    if field_v_per_cm <= 0.0:
        return 0.0
    return A_E_PER_M * math.exp(min(0.0, -B_E_V_PER_CM / field_v_per_cm))


def trapz(xs_um, vals):
    total = 0.0
    for i in range(len(xs_um) - 1):
        dx = (xs_um[i + 1] - xs_um[i]) * UM_TO_M
        total += 0.5 * (vals[i] + vals[i + 1]) * dx
    return total


def pchip_slopes(xs, ys):
    n = len(xs)
    h = [xs[i + 1] - xs[i] for i in range(n - 1)]
    delta = [(ys[i + 1] - ys[i]) / h[i] for i in range(n - 1)]
    d = [0.0] * n
    d[0] = delta[0]
    d[-1] = delta[-1]
    for i in range(1, n - 1):
        if delta[i - 1] * delta[i] <= 0.0:
            d[i] = 0.0
        else:
            w1 = 2 * h[i] + h[i - 1]
            w2 = h[i] + 2 * h[i - 1]
            d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
    return d


def pchip_eval_deriv(xs, ys, d, x):
    # locate interval
    i = 0
    while i < len(xs) - 2 and x > xs[i + 1]:
        i += 1
    h = xs[i + 1] - xs[i]
    t = (x - xs[i]) / h
    # derivative of Hermite cubic
    dydx = (
        ys[i] * (6 * t * t - 6 * t) / h
        + d[i] * (3 * t * t - 4 * t + 1)
        + ys[i + 1] * (-6 * t * t + 6 * t) / h
        + d[i + 1] * (3 * t * t - 2 * t)
    )
    return dydx


def fine_integral(xs_um, psi, use_pchip, n_sub=200):
    if use_pchip:
        d = pchip_slopes(xs_um, psi)
    total = 0.0
    x0, x1 = xs_um[0], xs_um[-1]
    dx_um = (x1 - x0) / n_sub
    prev_a = None
    prev_x = None
    for k in range(n_sub + 1):
        x = x0 + k * dx_um
        if use_pchip:
            slope = pchip_eval_deriv(xs_um, psi, d, x)
        else:
            # piecewise-linear psi -> segment slope
            i = 0
            while i < len(xs_um) - 2 and x > xs_um[i + 1]:
                i += 1
            slope = (psi[i + 1] - psi[i]) / (xs_um[i + 1] - xs_um[i])
        field = abs(slope) * GRAD_V_PER_UM_TO_V_PER_CM
        a = local_alpha(field)
        if prev_a is not None:
            total += 0.5 * (prev_a + a) * dx_um * UM_TO_M
        prev_a = a
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--biases", type=float, nargs="*", default=[-5.0, -10.0, -16.0, -20.0])
    parser.add_argument("--row-y", type=float, default=0.25, help="Interior row y [um].")
    parser.add_argument("--x-lo", type=float, default=0.5)
    parser.add_argument("--x-hi", type=float, default=1.5)
    args = parser.parse_args()

    data = parse(args.csv, args.row_y)
    print(f"# Ionization integral I = int alpha_n dx  (interior row y={args.row_y}, "
          f"x in [{args.x_lo},{args.x_hi}] um, van Overstraeten LOCAL model)")
    print(f"# {'bias':>6}{'I_sentaurus':>14}{'I_vela_exp':>13}{'I_coarse_loc':>14}"
          f"{'I_fineLIN':>12}{'I_finePCHIP':>13}{'PCHIP/Sen':>11}{'Vela/Sen':>10}")
    for bias in args.biases:
        pts = data.get(bias)
        if not pts:
            continue
        band = [(x, ps, as_, av) for (x, ps, as_, av) in pts
                if args.x_lo - 1e-9 <= x <= args.x_hi + 1e-9]
        xs = [p[0] for p in band]
        psi = [p[1] for p in band]
        alpha_s = [p[2] if math.isfinite(p[2]) else 0.0 for p in band]
        alpha_v = [p[3] if math.isfinite(p[3]) else 0.0 for p in band]

        # local-model alpha at coarse nodes from edge-averaged field
        # (nodal field via central/one-sided difference of psi)
        loc_alpha_nodes = []
        for i in range(len(xs)):
            if i == 0:
                slope = (psi[1] - psi[0]) / (xs[1] - xs[0])
            elif i == len(xs) - 1:
                slope = (psi[-1] - psi[-2]) / (xs[-1] - xs[-2])
            else:
                slope = (psi[i + 1] - psi[i - 1]) / (xs[i + 1] - xs[i - 1])
            loc_alpha_nodes.append(local_alpha(abs(slope) * GRAD_V_PER_UM_TO_V_PER_CM))

        I_sen = trapz(xs, alpha_s)
        I_vela = trapz(xs, alpha_v)
        I_coarse_loc = trapz(xs, loc_alpha_nodes)
        I_fine_lin = fine_integral(xs, psi, use_pchip=False)
        I_fine_pchip = fine_integral(xs, psi, use_pchip=True)

        pchip_ratio = I_fine_pchip / I_sen if I_sen > 0 else float("nan")
        vela_ratio = I_vela / I_sen if I_sen > 0 else float("nan")
        print(f"  {bias:>6.0f}{I_sen:>14.4g}{I_vela:>13.4g}{I_coarse_loc:>14.4g}"
              f"{I_fine_lin:>12.4g}{I_fine_pchip:>13.4g}{pchip_ratio:>11.3f}{vela_ratio:>10.3f}")


if __name__ == "__main__":
    main()
