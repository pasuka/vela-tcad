#!/usr/bin/env python3
"""Compare matched BVmethods device rows in native population coordinates.

Each row and its RHS are divided by that row's own diagonal. This removes
constant row scales, not state-dependent scaling derivatives or stabilization.
Contact potential rows/columns and the native circuit aliases are excluded.
The input state-match audit must be interpreted alongside these differences.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags, load_npz


def population_transform(nodes, blocks, population, g, n, voltage_ratio, density_scale):
    """dx/dz, x=(psi,fn,fp)/Vt; z=(psi/Vnative,n/Nref,p/Nref)."""
    if not (np.isfinite(voltage_ratio) and voltage_ratio > 0 and
            np.isfinite(density_scale) and density_scale > 0):
        raise ValueError("invalid coordinate scales")
    keys = list(zip(map(int, blocks), map(int, nodes)))
    if len(set(keys)) != len(keys) or any(b not in (0, 1, 2) or not 0 <= i < n for b, i in keys):
        raise ValueError("duplicate or invalid population coordinate")
    lookup = set(keys)
    rows, cols, values = [], [], []
    for j, (b, i) in enumerate(keys):
        if b == 0:
            for c in (0, 1, 2):
                if (c, i) in lookup:
                    rows.append(c * n + i); cols.append(j); values.append(voltage_ratio)
        else:
            denominator = population[b - 1, i] * g[b - 1, i]
            if not np.isfinite(denominator) or denominator <= 0:
                raise ValueError("nonpositive/nonfinite population derivative")
            rows.append(b * n + i); cols.append(j)
            values.append((-1 if b == 1 else 1) * density_scale / denominator)
    return coo_matrix((values, (rows, cols)), shape=(3 * n, len(keys))).tocsr()


def normalize_system(matrix, rhs):
    diagonal = matrix.diagonal()
    if (np.any(~np.isfinite(matrix.data)) or np.any(~np.isfinite(diagonal)) or
            np.any(diagonal == 0) or np.any(~np.isfinite(rhs))):
        raise ValueError("invalid diagonal or RHS")
    return diags(1 / diagonal) @ matrix, rhs / diagonal


def row_norms(a):
    return np.sqrt(np.asarray(a.multiply(a).sum(axis=1)).ravel())


def stats(a):
    return dict(zip(("min", "median", "p90", "p99", "max"), map(float, np.quantile(a, [0, .5, .9, .99, 1]))))


def compare(root, label):
    v = np.load(root / (label + "_native_vectors.npz"))
    m = np.load(root / (label + "_mapping.npz"))
    config = json.loads((root / (label + "_matrix/config.json")).read_text())
    mesh = json.loads(Path(config["mesh_file"]).read_text())
    n = len(mesh["nodes"])
    contacts = {i for c in mesh["contacts"] for i in c["node_ids"]}
    native_ids = np.flatnonzero(v["block"] != 1)
    blocks = np.where(v["block"][native_ids] == 0, 0, v["block"][native_ids] - 1)
    nodes = v["node"][native_ids]
    ids = blocks * n + nodes
    keep = np.array([not (b == 0 and i in contacts) for b, i in zip(blocks, nodes)])
    selected = np.flatnonzero(keep)
    t = population_transform(nodes, blocks, m["target"], m["g"], n,
                             float(m["native_potential_scale"]) / .025851999786435535, 14824351000.)
    with (root / (label + "_matrix/jacobian.csv")).open(newline="") as f:
        entries = list(csv.DictReader(f))
    j = coo_matrix(([float(r["value"]) for r in entries],
                    ([int(r["row"]) for r in entries], [int(r["col"]) for r in entries])), shape=(3*n, 3*n)).tocsr()
    x = np.genfromtxt(root / (label + "_matrix/rhs.csv"), delimiter=",", names=True)
    a = (j @ t)[ids, :][selected, :][:, selected].tocsr()
    c = load_npz(root / (label + "_native_matrix.npz"))[native_ids, :][:, native_ids][selected, :][:, selected].tocsr()
    ap, ba = normalize_system(a, x["rhs"][ids[keep]])
    cp, bc = normalize_system(c, v["rhs"][native_ids[keep]])
    difference = ap - cp
    row_error = row_norms(difference) / row_norms(cp)
    kb, kn = blocks[keep], nodes[keep]
    report = {"label": label, "common_free_unknowns": len(selected), "excluded_contact_psi": int((~keep).sum()),
              "scope": "fixed native drain voltage; free psi/n/p matched; contact constitutive and derived QF may differ; not a full circuit or nonlinear equivalence certificate", "blocks": {}, "focus": {}}
    for b, name in enumerate(("psi", "electron", "hole")):
        mask = kb == b
        report["blocks"][name] = {"row_relative_L2_difference": stats(row_error[mask]),
            "rhs_native_diag_normalized_L2": float(np.linalg.norm(bc[mask])),
            "rhs_vela_diag_normalized_L2": float(np.linalg.norm(ba[mask])),
            "rhs_difference_L2": float(np.linalg.norm(ba[mask]-bc[mask])),
            "subblock_relative_L2": {}}
        for col, cn in enumerate(("psi", "electron", "hole")):
            den = np.linalg.norm(row_norms(cp[mask, :][:, kb == col]))
            num = np.linalg.norm(row_norms(difference[mask, :][:, kb == col]))
            report["blocks"][name]["subblock_relative_L2"][cn] = float(num / den) if den else None
    for node in (129, 243, 785, 786, 789, 790, 802, 821, 850):
        for b, name in enumerate(("psi", "electron", "hole")):
            hit = np.flatnonzero((kb == b) & (kn == node))
            if len(hit):
                i = int(hit[0])
                report["focus"][f"{node}_{name}"] = {"row_relative_L2_difference": float(row_error[i]),
                    "native_rhs_over_diagonal": float(bc[i]), "vela_rhs_over_diagonal": float(ba[i])}
    with (root / (label + "_rows.csv")).open("w", newline="") as f:
        w = csv.writer(f); w.writerow(("block", "node", "native_diagonal", "vela_chained_diagonal", "native_rhs_over_diagonal", "vela_rhs_over_diagonal", "row_relative_L2_difference"))
        w.writerows(zip(kb, kn, c.diagonal(), a.diagonal(), bc, ba, row_error))
    # Store every element in the union of the two sparse patterns, including exact zeros on one side.
    pattern = ap.copy(); pattern.data[:] = 1
    pattern2 = cp.copy(); pattern2.data[:] = 1
    ii, jj = (pattern + pattern2).nonzero()
    av, cv = np.asarray(ap[ii, jj]).ravel(), np.asarray(cp[ii, jj]).ravel()
    with (root / (label + "_elements.csv")).open("w", newline="") as f:
        w = csv.writer(f); w.writerow(("row_block", "row_node", "col_block", "col_node", "native_normalized", "vela_normalized", "difference"))
        w.writerows(zip(kb[ii], kn[ii], kb[jj], kn[jj], cv, av, av-cv))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    args = p.parse_args()
    output = {"method": __doc__, "comparisons": [compare(args.root, s) for s in ("zero", "first", "second")]}
    (args.root / "element_comparison.json").write_text(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
