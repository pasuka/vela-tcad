#!/usr/bin/env python3
"""Boundary-neighbor weighted least-squares nodal electric-field scan (pn2d coarse7x3).

Motivation
----------
At bias -5V the coarse7x3 aligned comparison shows nodal *potential* differing
by <1% everywhere, yet the reconstructed nodal *electric field* differs from the
Sentaurus reference by up to ~8% at the boundary nodes. Interior nodes are
essentially exact (P1-linear reproduction) while the 1/d least-squares operator
aliases the strong x-gradient into the boundary rows.

This experiment separates interior and boundary nodes and re-weights the
least-squares stencil by node class: a neighbor that is an *interior* node keeps
weight coefficient 1, while a neighbor that is a *boundary* node has its 1/d
weight multiplied by a scan coefficient c (0.5, 1.0, 1.5, 2.0). Every node's
gradient is recomputed with this rule, and the field difference vs the Sentaurus
reference is reported split into interior / boundary / all.

The nodal potential fed to the operator is Vela's own full-precision potential
(``vela``, the default) so that c=1 reproduces Vela's shipped nodal field and the
scan shows how the boundary weighting moves the *actual* Vela field difference.

Reuses the validated parsing / topology / metric helpers from
``tryout_pn2d_efield_recovery.py``.

Units: potential [V], coordinates [um], field via grad(V/um) * 1e4 -> V/cm,
matching the aligned CSV field columns (validated by the c=1 reproduction check).
"""

from __future__ import annotations

import argparse
import csv
import math
import os

from tryout_pn2d_efield_recovery import (
    COORD_TOL,
    GRAD_TO_V_PER_CM,
    build_topology,
    parse_aligned,
    parse_elements,
    vector_error,
    weighted_rms_rel_error,
)

SCAN_COEFFS = (0.5, 1.0, 1.5, 2.0)


def recover_ls_boundary_weighted(coords, psi, neighbors, is_boundary, boundary_coeff):
    """1/d least-squares gradient with a per-neighbor boundary-class coefficient.

    A neighbor contributes with weight (1/d) * c, where c == 1 if the neighbor is
    an interior node and c == ``boundary_coeff`` if the neighbor is a boundary
    node. The center node's own class does not enter its stencil (a node is never
    its own neighbor); only the class of each *neighbor* re-weights the fit.
    """
    fields = {}
    for node, nb in neighbors.items():
        cx, cy = coords[node]
        sxx = sxy = syy = sxv = syv = 0.0
        for other in nb:
            ox, oy = coords[other]
            dx, dy = ox - cx, oy - cy
            dist = math.hypot(dx, dy)
            if dist <= 1.0e-30:
                continue
            weight = 1.0 / dist
            if is_boundary(other):
                weight *= boundary_coeff
            dv = psi[other] - psi[node]
            sxx += weight * dx * dx
            sxy += weight * dx * dy
            syy += weight * dy * dy
            sxv += weight * dx * dv
            syv += weight * dy * dv
        det = sxx * syy - sxy * sxy
        if abs(det) <= 1.0e-30 or not math.isfinite(det):
            fields[node] = None
            continue
        grad_x = (sxv * syy - syv * sxy) / det
        grad_y = (sxx * syv - sxy * sxv) / det
        fields[node] = (-grad_x * GRAD_TO_V_PER_CM, -grad_y * GRAD_TO_V_PER_CM)
    return fields


def make_is_boundary(coords, nodes):
    xs = [coords[n][0] for n in nodes]
    ys = [coords[n][1] for n in nodes]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    def is_boundary(node):
        x, y = coords[node]
        return (
            abs(x - xmin) < COORD_TOL
            or abs(x - xmax) < COORD_TOL
            or abs(y - ymin) < COORD_TOL
            or abs(y - ymax) < COORD_TOL
        )

    return is_boundary


def scope_error(field, ref, ref_mag, nodes, is_boundary, scope):
    errs, mags = [], []
    for n in nodes:
        if scope == "interior" and is_boundary(n):
            continue
        if scope == "boundary" and not is_boundary(n):
            continue
        err = vector_error(field[n], ref[n])
        if err is None or any(math.isnan(v) for v in ref[n]):
            continue
        errs.append(err)
        mags.append(ref_mag[n])
    return weighted_rms_rel_error(errs, mags)


def run_bias(bias, nodes, tris, potential_key, csv_rows):
    num_nodes = len(nodes)
    coords = {n: (nodes[n]["x"], nodes[n]["y"]) for n in nodes}
    psi = {n: nodes[n][potential_key] for n in nodes}
    _, neighbors = build_topology(tris, num_nodes)
    is_boundary = make_is_boundary(coords, nodes)

    ref = {n: (nodes[n]["ex_s"], nodes[n]["ey_s"]) for n in nodes}
    ref_mag = {n: math.hypot(*ref[n]) for n in nodes}
    max_ref = max((m for m in ref_mag.values() if math.isfinite(m)), default=0.0)

    # c == 1.0 reproduction check against the shipped Vela field (ex_v, ey_v).
    base = recover_ls_boundary_weighted(coords, psi, neighbors, is_boundary, 1.0)
    repro_max = 0.0
    for n in nodes:
        if base[n] is None or math.isnan(nodes[n].get("ex_v", float("nan"))):
            continue
        repro_max = max(
            repro_max,
            math.hypot(base[n][0] - nodes[n]["ex_v"], base[n][1] - nodes[n]["ey_v"]),
        )

    per_coeff = {}
    for coeff in SCAN_COEFFS:
        field = recover_ls_boundary_weighted(coords, psi, neighbors, is_boundary, coeff)
        per_coeff[coeff] = {
            "field": field,
            "all": scope_error(field, ref, ref_mag, nodes, is_boundary, "all"),
            "interior": scope_error(field, ref, ref_mag, nodes, is_boundary, "interior"),
            "boundary": scope_error(field, ref, ref_mag, nodes, is_boundary, "boundary"),
        }

    # Per-node CSV: boundary nodes only, one row per (node, coeff).
    for n in nodes:
        if not is_boundary(n):
            continue
        for coeff in SCAN_COEFFS:
            f = per_coeff[coeff]["field"][n]
            err = vector_error(f, ref[n])
            csv_rows.append(
                {
                    "bias_V": bias,
                    "node_id": n,
                    "x_um": coords[n][0],
                    "y_um": coords[n][1],
                    "boundary_coeff": coeff,
                    "ref_ex": ref[n][0],
                    "ref_ey": ref[n][1],
                    "ref_mag": ref_mag[n],
                    "field_ex": "" if f is None else f[0],
                    "field_ey": "" if f is None else f[1],
                    "abs_err": "" if err is None else err,
                    "rel_err": (
                        ""
                        if err is None or ref_mag[n] <= 0
                        else err / ref_mag[n]
                    ),
                }
            )

    return per_coeff, repro_max, ref_mag, max_ref, coords, is_boundary


def main():
    default_dir = os.path.join(
        "build-release",
        "reference_tcad",
        "pn2d_sentaurus2018_coarse7x3",
        "reports",
        "coarse_previous_full20_vector_current_20260630",
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--elements",
        default=os.path.join(
            default_dir, "sentaurus_multibias", "sentaurus_-5v", "elements.csv"
        ),
    )
    parser.add_argument(
        "--aligned",
        default=os.path.join(default_dir, "coarse_node_field_compare_aligned.csv"),
    )
    parser.add_argument(
        "--potential",
        choices=["sentaurus", "vela"],
        default="vela",
        help="Nodal potential fed to the operator (default: Vela full-precision).",
    )
    parser.add_argument(
        "--bias",
        type=float,
        default=-5.0,
        help="Bias to evaluate (default -5.0 V). Use 'all' via --all-biases.",
    )
    parser.add_argument("--all-biases", action="store_true")
    parser.add_argument(
        "--out",
        default=os.path.join("build", "diagnostics", "efield_boundary_weight_scan.csv"),
    )
    args = parser.parse_args()

    tris = parse_elements(args.elements)
    data = parse_aligned(args.aligned)
    potential_key = "psi_s" if args.potential == "sentaurus" else "psi_v"

    if args.all_biases:
        biases = sorted(data.keys())
    else:
        # match the requested bias against available keys with a tolerance
        biases = [
            b for b in sorted(data.keys()) if abs(b - args.bias) < 1.0e-6
        ] or [args.bias]

    csv_rows = []
    print(f"Potential input : {args.potential} ({potential_key})")
    print(f"Triangles       : {len(tris)}")
    print(f"Boundary coeffs : {', '.join(f'{c:g}' for c in SCAN_COEFFS)}")
    print("Rule: neighbor weight = (1/d) * c_boundary if neighbor is boundary else (1/d)\n")

    for bias in biases:
        nodes = data.get(bias)
        if not nodes:
            print(f"(no data for bias {bias:g} V)")
            continue
        per_coeff, repro_max, ref_mag, max_ref, coords, is_boundary = run_bias(
            bias, nodes, tris, potential_key, csv_rows
        )
        n_bnd = sum(1 for n in nodes if is_boundary(n))
        n_int = len(nodes) - n_bnd
        print("=" * 74)
        print(
            f"bias = {bias:g} V   nodes = {len(nodes)} "
            f"(interior {n_int}, boundary {n_bnd})   peak|E| = {max_ref:.4g} V/cm"
        )
        print(f"  (check) max |c=1 field - shipped Vela field| = {repro_max:.4g} V/cm")
        print("  field-normalized L2 error vs Sentaurus reference:")
        print(f"    {'c_boundary':>10}{'all':>12}{'interior':>12}{'boundary':>12}")
        for coeff in SCAN_COEFFS:
            s = per_coeff[coeff]
            print(
                f"    {coeff:>10.2f}{s['all']:>12.4%}"
                f"{s['interior']:>12.4%}{s['boundary']:>12.4%}"
            )
        # best boundary coefficient
        best = min(SCAN_COEFFS, key=lambda c: per_coeff[c]["boundary"])
        print(
            f"  -> boundary L2 error minimized at c_boundary = {best:g} "
            f"({per_coeff[best]['boundary']:.4%})"
        )

    if csv_rows:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"\nPer-node boundary table written to {args.out} ({len(csv_rows)} rows)")


if __name__ == "__main__":
    main()
