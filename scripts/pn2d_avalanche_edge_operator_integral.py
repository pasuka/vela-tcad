#!/usr/bin/env python3
"""Verify Vela's ACTUAL edge-assembled ionization integral vs exported nodal alpha.

Vela's solver assembles the avalanche source per EDGE (computeSgEdgeCurrentAvalanche
SourceRecords): alpha is evaluated from edgeAveragedCellScalarGradient(phin) =
the norm of the mean P1 cell-gradient of phin over the cells adjacent to the edge.
The exported per-NODE alpha instead averages the cell-gradient over ALL cells
touching a node, which smooths the junction peak and UNDER-states the integral.

This script reconstructs Vela's real EDGE operator from the coarse7x3 triangle
mesh + exported nodal phin (electron_qf), then integrates alpha along the junction
current path (interior-row horizontal edges 4-7-10-13-16, x=0.5..1.5um) and
compares against:
  * Vela exported NODAL alpha trapezoid  (the under-stated diagnostic)
  * Sentaurus exported NODAL alpha trapezoid (target)
  * the same EDGE operator applied to Sentaurus's OWN phin (local upper check)

van Overstraeten electron single branch (== Sentaurus .par): gamma=1,
  alpha[1/m] = 7.03e7 * exp(-1.231e6 / E[V/cm]).
Std lib only. phin/psi [V], coords [um], grad*1e4 -> V/cm, length*1e-6 -> m.
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

BASE = os.path.join(
    "build-release", "reference_tcad", "pn2d_sentaurus2018_coarse7x3",
    "reports", "coarse_previous_full20_vector_current_20260630",
)
DEFAULT_ALIGNED = os.path.join(BASE, "coarse_node_field_compare_aligned.csv")
DEFAULT_MESH_DIR = os.path.join(BASE, "sentaurus_multibias", "sentaurus_-5v")

# junction current path along interior row y=0.25 (x = 0.5,0.75,1.0,1.25,1.5)
PATH_NODES = [4, 7, 10, 13, 16]
PATH_EDGES = [(4, 7), (7, 10), (10, 13), (13, 16)]


def _f(t):
    try:
        return float(t)
    except (TypeError, ValueError):
        return float("nan")


def read_nodes(mesh_dir):
    coord = {}
    with open(os.path.join(mesh_dir, "nodes.csv"), newline="") as h:
        r = csv.DictReader(h)
        for row in r:
            coord[int(row["id"])] = (float(row["x_um"]), float(row["y_um"]))
    return coord


def read_triangles(mesh_dir):
    tris = []
    with open(os.path.join(mesh_dir, "elements.csv"), newline="") as h:
        r = csv.DictReader(h)
        for row in r:
            tris.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
    return tris


def read_aligned(path):
    """Return {bias: {node: {phin_v, phin_s, alpha_v, alpha_s}}}."""
    out = defaultdict(lambda: defaultdict(dict))
    with open(path, newline="") as h:
        r = csv.reader(h)
        header = [c.strip() for c in next(r)]
        idx = {n: i for i, n in enumerate(header)}
        for row in r:
            if len(row) < len(header):
                continue
            q = row[idx["quantity"]].strip()
            if q not in ("electron_qf", "electron_alpha_avalanche"):
                continue
            bias = _f(row[idx["bias_V"]])
            node = int(_f(row[idx["node_id"]]))
            e = out[bias][node]
            if q == "electron_qf":
                e["phin_s"] = _f(row[idx["sentaurus_value"]])
                e["phin_v"] = _f(row[idx["vela_value_scaled_to_sentaurus_units"]])
            else:
                e["alpha_s"] = _f(row[idx["sentaurus_value"]])
                e["alpha_v"] = _f(row[idx["vela_value_scaled_to_sentaurus_units"]])
    return out


def tri_gradient(coord, tri, values):
    (n1, n2, n3) = tri
    (x1, y1), (x2, y2), (x3, y3) = coord[n1], coord[n2], coord[n3]
    f1, f2, f3 = values[n1], values[n2], values[n3]
    det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if abs(det) < 1e-18:
        return (0.0, 0.0)
    dfdx = ((f2 - f1) * (y3 - y1) - (f3 - f1) * (y2 - y1)) / det
    dfdy = ((f3 - f1) * (x2 - x1) - (f2 - f1) * (x3 - x1)) / det
    return (dfdx, dfdy)


def edge_to_cells(tris):
    m = defaultdict(list)
    for ci, (a, b, c) in enumerate(tris):
        for u, v in ((a, b), (b, c), (c, a)):
            m[(min(u, v), max(u, v))].append(ci)
    return m


def local_alpha(field_v_per_cm):
    if field_v_per_cm <= 0.0:
        return 0.0
    return A_E_PER_M * math.exp(-B_E_V_PER_CM / field_v_per_cm)


def edge_field_vcm(coord, tris, e2c, grads, edge):
    key = (min(edge), max(edge))
    cells = e2c.get(key, [])
    if not cells:
        return 0.0
    gx = sum(grads[c][0] for c in cells) / len(cells)
    gy = sum(grads[c][1] for c in cells) / len(cells)
    return math.hypot(gx, gy) * GRAD_V_PER_UM_TO_V_PER_CM


def node_field_vcm(coord, tris, node_cells, grads, node):
    cells = node_cells.get(node, [])
    if not cells:
        return 0.0
    gx = sum(grads[c][0] for c in cells) / len(cells)
    gy = sum(grads[c][1] for c in cells) / len(cells)
    return math.hypot(gx, gy) * GRAD_V_PER_UM_TO_V_PER_CM


def trapz_nodes(coord, nodes, node_val):
    total = 0.0
    for i in range(len(nodes) - 1):
        dx = abs(coord[nodes[i + 1]][0] - coord[nodes[i]][0]) * UM_TO_M
        total += 0.5 * (node_val[nodes[i]] + node_val[nodes[i + 1]]) * dx
    return total


def edge_integral(coord, alpha_edge):
    total = 0.0
    for (a, b), al in alpha_edge.items():
        (x0, y0), (x1, y1) = coord[a], coord[b]
        length = math.hypot(x1 - x0, y1 - y0) * UM_TO_M
        total += al * length
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--aligned", default=DEFAULT_ALIGNED)
    ap.add_argument("--mesh-dir", default=DEFAULT_MESH_DIR)
    ap.add_argument("--biases", type=float, nargs="*", default=[-5.0, -10.0, -16.0, -20.0])
    args = ap.parse_args()

    coord = read_nodes(args.mesh_dir)
    tris = read_triangles(args.mesh_dir)
    e2c = edge_to_cells(tris)
    node_cells = defaultdict(list)
    for ci, (a, b, c) in enumerate(tris):
        for n in (a, b, c):
            node_cells[n].append(ci)
    data = read_aligned(args.aligned)

    print("# Vela ACTUAL edge-operator vs exported nodal alpha, junction path "
          "nodes 4-7-10-13-16 (x=0.5..1.5um)")
    print(f"# {'bias':>6}{'I_edge_vela':>13}{'I_node_vela':>13}{'I_exp_vela':>12}"
          f"{'I_exp_sen':>11}{'I_edge_senPhin':>16}{'edge/exp_sen':>13}{'exp_v/exp_s':>12}")

    for bias in args.biases:
        nd = data.get(bias)
        if not nd:
            continue
        phin_v = {n: nd[n].get("phin_v", float("nan")) for n in nd}
        phin_s = {n: nd[n].get("phin_s", float("nan")) for n in nd}
        alpha_v = {n: (nd[n].get("alpha_v") if math.isfinite(nd[n].get("alpha_v", float("nan"))) else 0.0) for n in nd}
        alpha_s = {n: (nd[n].get("alpha_s") if math.isfinite(nd[n].get("alpha_s", float("nan"))) else 0.0) for n in nd}

        grads_v = [tri_gradient(coord, t, phin_v) for t in tris]
        grads_s = [tri_gradient(coord, t, phin_s) for t in tris]

        # Vela edge-operator alpha along path edges
        alpha_edge_v = {}
        alpha_edge_sphin = {}
        for e in PATH_EDGES:
            fe_v = edge_field_vcm(coord, tris, e2c, grads_v, e)
            fe_s = edge_field_vcm(coord, tris, e2c, grads_s, e)
            alpha_edge_v[e] = local_alpha(fe_v)
            alpha_edge_sphin[e] = local_alpha(fe_s)

        # Vela node-operator alpha (box over all cells at node) reconstructed
        node_alpha_v = {n: local_alpha(node_field_vcm(coord, tris, node_cells, grads_v, n))
                        for n in PATH_NODES}

        I_edge_v = edge_integral(coord, alpha_edge_v)
        I_node_v = trapz_nodes(coord, PATH_NODES, node_alpha_v)
        I_exp_v = trapz_nodes(coord, PATH_NODES, alpha_v)
        I_exp_s = trapz_nodes(coord, PATH_NODES, alpha_s)
        I_edge_sphin = edge_integral(coord, alpha_edge_sphin)

        edge_ratio = I_edge_v / I_exp_s if I_exp_s > 0 else float("nan")
        exp_ratio = I_exp_v / I_exp_s if I_exp_s > 0 else float("nan")
        print(f"  {bias:>6.0f}{I_edge_v:>13.4g}{I_node_v:>13.4g}{I_exp_v:>12.4g}"
              f"{I_exp_s:>11.4g}{I_edge_sphin:>16.4g}{edge_ratio:>13.3f}{exp_ratio:>12.3f}")

    print()
    print("# Per-edge detail (Vela edge operator, field V/cm & alpha 1/m):")
    for bias in args.biases:
        nd = data.get(bias)
        if not nd:
            continue
        phin_v = {n: nd[n].get("phin_v", float("nan")) for n in nd}
        phin_s = {n: nd[n].get("phin_s", float("nan")) for n in nd}
        grads_v = [tri_gradient(coord, t, phin_v) for t in tris]
        grads_s = [tri_gradient(coord, t, phin_s) for t in tris]
        print(f"  bias={bias:.0f}V")
        for e in PATH_EDGES:
            fv = edge_field_vcm(coord, tris, e2c, grads_v, e)
            fs = edge_field_vcm(coord, tris, e2c, grads_s, e)
            ratio = fv / fs if fs > 0 else float("nan")
            print(f"    edge {e}: F_velaPhin={fv:.4g}  F_senPhin={fs:.4g} V/cm  "
                  f"Fv/Fs={ratio:.3f}  alpha_v={local_alpha(fv):.4g} alpha_s={local_alpha(fs):.4g}")


if __name__ == "__main__":
    main()
