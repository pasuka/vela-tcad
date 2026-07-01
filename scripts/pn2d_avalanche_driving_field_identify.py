#!/usr/bin/env python3
"""Identify which field Sentaurus's avalanche alpha actually uses.

Both codes share IDENTICAL van Overstraeten electron parameters (single branch):
    alpha[1/m] = gamma*a * exp(-gamma*b / E),  a=7.03e7 /m, b=1.231e6 V/cm, gamma=1.

So the *effective* driving field each code evaluated alpha at can be inverted:
    E_eff = gamma*b / ln(gamma*a / alpha)                (valid for 0 < alpha < gamma*a).

We compare E_eff(Sentaurus) and E_eff(Vela) against two candidate driving fields,
both reconstructed on the SAME mesh with the SAME area-weighted cell-gradient
(box) operator:
    F_phin = |grad phi_n|   (quasi-Fermi-potential gradient  -- Vela's config)
    F_psi  = |grad psi|     (electrostatic field             -- Sentaurus default?)

If E_eff(Sentaurus) tracks F_psi (broad depletion field) rather than F_phin
(sharp junction spike), the BV gap is a driving-force *definition* difference:
Vela drives avalanche with the too-peaked quasi-Fermi gradient while Sentaurus
uses the broader electrostatic (current-parallel) field.

Interior nodes only for the operator-exact comparison. Std lib only.
Units: phi/psi [V], coords [um], field [V/cm], alpha [1/m].
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

GRAD_TO_V_PER_CM = 1.0e4
COORD_TOL = 1.0e-9

DEFAULT_DIR = os.path.join(
    "build-release", "reference_tcad", "pn2d_sentaurus2018_coarse7x3",
    "reports", "coarse_previous_full20_vector_current_20260630",
)

# van Overstraeten electron single-branch constants (match Sentaurus .par).
GAMMA = 1.0
A_E_PER_M = 7.03e7        # gamma*a prefactor [1/m]
B_E_V_PER_CM = 1.231e6    # gamma*b critical field [V/cm]

WANT = {
    "potential": ("psi_s", "psi_v"),
    "electron_qf": ("phin_s", "phin_v"),
    "electron_alpha_avalanche": ("alpha_s", "alpha_v"),
    "electric_field_mag": ("emag_s", "emag_v"),
}


def _fnum(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return float("nan")


def parse_elements(path):
    tris = []
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle):
            tris.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
    return tris


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
            e = data[bias][nid]
            e["x"] = _fnum(row[idx["x_um"]])
            e["y"] = _fnum(row[idx["y_um"]])
            sk, vk = WANT[q]
            e[sk] = _fnum(row[idx["sentaurus_value"]])
            e[vk] = _fnum(row[idx["vela_value_scaled_to_sentaurus_units"]])
    return data


def cell_gradient(coords, values, tri):
    i, j, k = tri
    x0, y0 = coords[i]
    x1, y1 = coords[j]
    x2, y2 = coords[k]
    det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
    if abs(det) <= 0.0 or not math.isfinite(det):
        return None
    dv1 = values[j] - values[i]
    dv2 = values[k] - values[i]
    gx = (dv1 * (y2 - y0) - (y1 - y0) * dv2) / det
    gy = ((x1 - x0) * dv2 - dv1 * (x2 - x0)) / det
    return gx, gy, 0.5 * abs(det)


def node_cell_grad_mag(coords, values, node_cells, tris):
    out = {}
    for node, cells in node_cells.items():
        gx = gy = tot = 0.0
        for c in cells:
            cg = cell_gradient(coords, values, tris[c])
            if cg is None:
                continue
            gx += cg[2] * cg[0]
            gy += cg[2] * cg[1]
            tot += cg[2]
        out[node] = math.hypot(gx / tot, gy / tot) * GRAD_TO_V_PER_CM if tot > 0 else 0.0
    return out


def eff_field_from_alpha(alpha_per_m):
    """Invert alpha = gamma*a*exp(-gamma*b/E) -> E [V/cm]."""
    if not math.isfinite(alpha_per_m) or alpha_per_m <= 0.0 or alpha_per_m >= A_E_PER_M:
        return float("nan")
    return B_E_V_PER_CM / math.log(A_E_PER_M / alpha_per_m)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned", default=os.path.join(DEFAULT_DIR, "coarse_node_field_compare_aligned.csv"))
    parser.add_argument("--elements", default=os.path.join(DEFAULT_DIR, "sentaurus_multibias", "sentaurus_-5v", "elements.csv"))
    parser.add_argument("--biases", type=float, nargs="*", default=[-5.0, -10.0, -16.0, -20.0])
    parser.add_argument("--band", type=float, nargs=2, default=[0.5, 1.5])
    args = parser.parse_args()

    tris = parse_elements(args.elements)
    data = parse_aligned(args.aligned)

    for bias in args.biases:
        nodes = data.get(bias)
        if not nodes:
            continue
        coords = {n: (nodes[n]["x"], nodes[n]["y"]) for n in nodes}
        node_cells = defaultdict(list)
        for c, tri in enumerate(tris):
            for nd in tri:
                node_cells[nd].append(c)
        for n in nodes:
            node_cells.setdefault(n, [])

        xs = [coords[n][0] for n in nodes]
        ys = [coords[n][1] for n in nodes]
        xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

        def interior(n):
            x, y = coords[n]
            return not (abs(x - xmin) < COORD_TOL or abs(x - xmax) < COORD_TOL
                        or abs(y - ymin) < COORD_TOL or abs(y - ymax) < COORD_TOL)

        F_phin_s = node_cell_grad_mag(coords, {n: nodes[n]["phin_s"] for n in nodes}, node_cells, tris)
        F_psi_s = node_cell_grad_mag(coords, {n: nodes[n]["psi_s"] for n in nodes}, node_cells, tris)

        print("=" * 100)
        print(f"bias = {bias:g} V   (interior column; E_eff back-solved from exported alpha via shared vO model)")
        print(f"  {'node':>4}{'x':>6}"
              f"{'alpha_sen':>12}{'E*_sen':>10}"
              f"{'|grad phin|':>13}{'|grad psi|':>12}{'|E|_exp':>11}"
              f"{'E*/phin':>9}{'E*/psi':>9}")
        band_int = [n for n in nodes if interior(n)
                    and args.band[0] - 1e-9 <= coords[n][0] <= args.band[1] + 1e-9]
        for n in sorted(band_int, key=lambda m: coords[m][0]):
            a_s = nodes[n]["alpha_s"]
            estar = eff_field_from_alpha(a_s)
            fph = F_phin_s[n]
            fps = F_psi_s[n]
            eexp = nodes[n].get("emag_s", float("nan"))
            r_ph = estar / fph if fph > 0 else float("nan")
            r_ps = estar / fps if fps > 0 else float("nan")
            print(f"  {n:>4}{coords[n][0]:>6.2f}"
                  f"{a_s:>12.4g}{estar:>10.4g}"
                  f"{fph:>13.4g}{fps:>12.4g}{eexp:>11.4g}"
                  f"{r_ph:>9.3f}{r_ps:>9.3f}")
        print()


if __name__ == "__main__":
    main()
