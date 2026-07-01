#!/usr/bin/env python3
"""Isolate the error in the impact-ionization driving field |grad phi_n|.

Vela drives avalanche with the quasi-Fermi-potential gradient magnitude,
reconstructed by an area-weighted P1 cell-gradient operator (box method) and
fed into van Overstraeten:  alpha = gamma*A * exp(-gamma*B / F),  F = |grad phi|.

The exported CSV carries the *nodal* quasi-Fermi potential (phi_n, phi_p) and
the *nodal* ionization coefficients (electron/hole_alpha_avalanche) for both
Vela and Sentaurus, but NOT the reconstructed driving field F itself. This
script re-derives F with the same cell-gradient operator so we can separate:

  (A) SOLUTION error   -- feed each code's OWN phi_n through the SAME operator
                          and compare F_vela vs F_sen (operator cancels in the
                          ratio, so this is purely the phi_n solution effect);
  (B) OPERATOR / MODEL -- regress ln(alpha) against 1/F for each code over the
                          avalanche band. If the reconstruction matches the
                          code's internal driving field, the fit is a clean
                          line (R^2 ~ 1) with slope = -gamma*B and
                          intercept = ln(gamma*A). Per-node residuals from that
                          line expose where alpha is NOT consistent with
                          cell-gradient(phi_n) -- i.e. boundary operator / branch
                          artifacts beyond the phi_n solution.

Interior vs boundary rows are reported separately because the box cell-gradient
is exact on interior P1 patches but aliases the strong x-gradient into the
transverse direction on the structured-mesh boundary rows.

Pure standard library. Units: phi [V], coords [um], F [V/cm], alpha [1/cm].
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

GRAD_TO_V_PER_CM = 1.0e4  # V/um -> V/cm
COORD_TOL = 1.0e-9

DEFAULT_DIR = os.path.join(
    "build-release", "reference_tcad", "pn2d_sentaurus2018_coarse7x3",
    "reports", "coarse_previous_full20_vector_current_20260630",
)


def _fnum(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def parse_elements(path):
    tris = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            tris.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
    return tris


WANT = {
    "electron_qf": ("phin_s", "phin_v"),
    "hole_qf": ("phip_s", "phip_v"),
    "electron_alpha_avalanche": ("alpha_e_s", "alpha_e_v"),
    "hole_alpha_avalanche": ("alpha_p_s", "alpha_p_v"),
}


def parse_aligned(path):
    data = defaultdict(lambda: defaultdict(dict))
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = [c.strip() for c in next(reader)]
        idx = {name: i for i, name in enumerate(header)}
        for row in reader:
            if len(row) < len(header):
                continue
            q = row[idx["quantity"]].strip()
            if q not in WANT:
                continue
            bias = _fnum(row[idx["bias_V"]])
            nid = int(_fnum(row[idx["node_id"]]))
            entry = data[bias][nid]
            entry["x"] = _fnum(row[idx["x_um"]])
            entry["y"] = _fnum(row[idx["y_um"]])
            s_key, v_key = WANT[q]
            entry[s_key] = _fnum(row[idx["sentaurus_value"]])
            entry[v_key] = _fnum(row[idx["vela_value_scaled_to_sentaurus_units"]])
    return data


def cell_gradient(coords, values, tri):
    """P1 gradient (dv/dx, dv/dy) and area for one triangle; None if degenerate."""
    i, j, k = tri
    x0, y0 = coords[i]
    x1, y1 = coords[j]
    x2, y2 = coords[k]
    det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
    area = 0.5 * abs(det)
    if area <= 0.0 or not math.isfinite(det):
        return None
    dv1 = values[j] - values[i]
    dv2 = values[k] - values[i]
    gx = (dv1 * (y2 - y0) - (y1 - y0) * dv2) / det
    gy = ((x1 - x0) * dv2 - dv1 * (x2 - x0)) / det
    return gx, gy, area


def node_cell_gradient_magnitude(coords, values, node_cells, tris):
    """Area-weighted average of P1 cell gradients, then magnitude (V/cm).

    Mirrors Vela's computeNodeCellGradientMagnitudes.
    """
    fields = {}
    for node, cells in node_cells.items():
        gx = gy = tot = 0.0
        for c in cells:
            cg = cell_gradient(coords, values, tris[c])
            if cg is None:
                continue
            cgx, cgy, area = cg
            gx += area * cgx
            gy += area * cgy
            tot += area
        if tot > 0.0:
            fields[node] = math.hypot(gx / tot, gy / tot) * GRAD_TO_V_PER_CM
        else:
            fields[node] = 0.0
    return fields


def linear_fit(xs, ys):
    """Least-squares y = a + b*x; return (a, b, r2, n)."""
    n = len(xs)
    if n < 2:
        return (float("nan"), float("nan"), float("nan"), n)
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    if sxx <= 0.0:
        return (float("nan"), float("nan"), float("nan"), n)
    b = sxy / sxx
    a = my - b * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    return (a, b, r2, n)


def analyze_carrier(bias, nodes, tris, node_cells, coords, is_boundary,
                    phi_s_key, phi_v_key, a_s_key, a_v_key, band, sig_frac, label):
    phi_s = {n: nodes[n][phi_s_key] for n in nodes}
    phi_v = {n: nodes[n][phi_v_key] for n in nodes}
    F_s = node_cell_gradient_magnitude(coords, phi_s, node_cells, tris)
    F_v = node_cell_gradient_magnitude(coords, phi_v, node_cells, tris)

    band_nodes = [n for n in nodes if band[0] - 1e-9 <= coords[n][0] <= band[1] + 1e-9]
    a_s = {n: nodes[n][a_s_key] for n in nodes}
    a_v = {n: nodes[n][a_v_key] for n in nodes}

    # Regression ln(alpha) vs 1/F over significant band nodes, per code.
    def fit(Fmap, amap):
        amax = max((amap[n] for n in band_nodes if math.isfinite(amap[n])), default=0.0)
        floor = sig_frac * amax
        xs, ys, used = [], [], []
        for n in band_nodes:
            F = Fmap[n]
            a = amap[n]
            if F > 0.0 and math.isfinite(a) and a > floor and a > 0.0:
                xs.append(1.0 / F)
                ys.append(math.log(a))
                used.append(n)
        return linear_fit(xs, ys), used, floor

    (a_s_fit, s_used, s_floor) = fit(F_s, a_s)
    (a_v_fit, v_used, v_floor) = fit(F_v, a_v)

    print(f"  --- {label} ---")
    print(f"    ln(alpha) = intercept + slope*(1/F)   [slope=-gamma*B, intercept=ln(gamma*A)]")
    print(f"    Sentaurus fit: slope={a_s_fit[1]:.4g}  intercept={a_s_fit[0]:.4g}  R2={a_s_fit[2]:.5f}  n={a_s_fit[3]}")
    print(f"    Vela      fit: slope={a_v_fit[1]:.4g}  intercept={a_v_fit[0]:.4g}  R2={a_v_fit[2]:.5f}  n={a_v_fit[3]}")

    # (A) SOLUTION contribution: F_vela/F_sen (same operator) interior vs boundary.
    def med(vals):
        v = sorted(x for x in vals if math.isfinite(x))
        return v[len(v) // 2] if v else float("nan")

    ratios_int, ratios_bnd = [], []
    for n in band_nodes:
        if F_s[n] > 0.0 and F_v[n] > 0.0:
            r = F_v[n] / F_s[n]
            (ratios_bnd if is_boundary(n) else ratios_int).append(r)
    print(f"    (A) driving-field ratio F_vela/F_sen (same operator, phi_n effect only):")
    print(f"        interior median={med(ratios_int):.4f}   boundary median={med(ratios_bnd):.4f}")

    # Per-node table over the band.
    print(f"    {'node':>4}{'x':>6}{'y':>6}{'bnd':>4}"
          f"{'F_sen':>11}{'F_vela':>11}{'F_v/F_s':>9}"
          f"{'a_sen':>11}{'a_vela':>11}{'a_v/a_s':>9}{'vfit_dex':>9}")
    for n in sorted(band_nodes, key=lambda m: (coords[m][0], coords[m][1])):
        F = F_v[n]
        # residual of this node's exported Vela alpha from the Vela fit line (in decades)
        if F > 0.0 and a_v[n] > 0.0 and math.isfinite(a_v_fit[0]):
            pred = a_v_fit[0] + a_v_fit[1] * (1.0 / F)
            vdex = (math.log(a_v[n]) - pred) / math.log(10.0)
        else:
            vdex = float("nan")
        fr = (F_v[n] / F_s[n]) if (F_s[n] > 0 and F_v[n] > 0) else float("nan")
        ar = (a_v[n] / a_s[n]) if (a_s[n] > 0 and a_v[n] > 0) else float("nan")
        print(f"    {n:>4}{coords[n][0]:>6.2f}{coords[n][1]:>6.2f}{int(is_boundary(n)):>4}"
              f"{F_s[n]:>11.4g}{F_v[n]:>11.4g}{fr:>9.3f}"
              f"{a_s[n]:>11.4g}{a_v[n]:>11.4g}{ar:>9.3g}{vdex:>9.3f}")
    print()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned", default=os.path.join(DEFAULT_DIR, "coarse_node_field_compare_aligned.csv"))
    parser.add_argument("--elements", default=os.path.join(DEFAULT_DIR, "sentaurus_multibias", "sentaurus_-5v", "elements.csv"))
    parser.add_argument("--biases", type=float, nargs="*", default=[-5.0, -10.0, -16.0, -20.0])
    parser.add_argument("--band", type=float, nargs=2, default=[0.5, 1.25], help="Avalanche band x range [um].")
    parser.add_argument("--sig-frac", type=float, default=0.02, help="alpha significance floor for the fit.")
    parser.add_argument("--carrier", choices=["electron", "hole", "both"], default="electron")
    args = parser.parse_args()

    tris = parse_elements(args.elements)
    data = parse_aligned(args.aligned)

    for bias in args.biases:
        nodes = data.get(bias)
        if not nodes:
            continue
        coords = {n: (nodes[n]["x"], nodes[n]["y"]) for n in nodes}
        node_cells = defaultdict(list)
        for c, (i, j, k) in enumerate(tris):
            for node in (i, j, k):
                node_cells[node].append(c)
        for n in nodes:
            node_cells.setdefault(n, [])

        xs = [coords[n][0] for n in nodes]
        ys = [coords[n][1] for n in nodes]
        xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

        def is_boundary(n):
            x, y = coords[n]
            return (abs(x - xmin) < COORD_TOL or abs(x - xmax) < COORD_TOL
                    or abs(y - ymin) < COORD_TOL or abs(y - ymax) < COORD_TOL)

        print("=" * 96)
        print(f"bias = {bias:g} V   band x in [{args.band[0]},{args.band[1]}] um")
        if args.carrier in ("electron", "both"):
            analyze_carrier(bias, nodes, tris, node_cells, coords, is_boundary,
                            "phin_s", "phin_v", "alpha_e_s", "alpha_e_v",
                            args.band, args.sig_frac, "ELECTRON  F=|grad phi_n|")
        if args.carrier in ("hole", "both"):
            analyze_carrier(bias, nodes, tris, node_cells, coords, is_boundary,
                            "phip_s", "phip_v", "alpha_p_s", "alpha_p_v",
                            args.band, args.sig_frac, "HOLE  F=|grad phi_p|")


if __name__ == "__main__":
    main()
