#!/usr/bin/env python3
"""Nodal E-vector angle difference vs triangle interior angles (pn2d coarse7x3).

For every node the electric-field *direction* is theta = atan2(E_y, E_x). This
script computes theta from Vela's shipped nodal field (ex_v, ey_v) and from the
Sentaurus reference (ex_s, ey_s), forms the wrapped angular difference
dtheta = wrap(theta_vela - theta_sentaurus) in degrees, and analyses how that
difference relates to the interior angles of the triangles incident to each node.

Interior-angle metrics per node (from elements.csv geometry):
* vertex_angle_min/max  -- the interior angle *subtended at this node* in each
                           incident triangle (min and max across those triangles).
* angle_sum             -- sum of the subtended vertex angles at the node
                           (~360 deg interior, <360 deg on the boundary; a direct
                           geometric boundary/corner detector).
* tri_min_angle         -- the smallest interior angle of *any* incident triangle
                           (a per-node mesh-skewness / quality indicator).

Because the field angle is meaningless where |E| ~ 0, nodes are filtered by a
significance gate |E| >= sig_frac * peak|E| (per bias) before correlation.

Units: coordinates [um], field via the aligned CSV columns, angles [deg].
Pure standard library (Pearson correlation computed manually).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

from tryout_pn2d_efield_recovery import parse_aligned, parse_elements


def wrap_deg(delta):
    """Wrap an angle difference to (-180, 180] degrees."""
    return (delta + 180.0) % 360.0 - 180.0


def vertex_angle(coords, apex, b, c):
    """Interior angle (deg) at vertex `apex` in triangle (apex, b, c)."""
    ax, ay = coords[apex]
    bx, by = coords[b]
    cx, cy = coords[c]
    ux, uy = bx - ax, by - ay
    vx, vy = cx - ax, cy - ay
    dot = ux * vx + uy * vy
    cross = ux * vy - uy * vx
    return math.degrees(math.atan2(abs(cross), dot))


def pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0.0 or syy <= 0.0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def build_node_angles(tris, num_nodes, coords):
    """Return per-node lists of subtended vertex angles and incident-tri min angles."""
    vertex_angles = defaultdict(list)
    tri_min_angle_at = defaultdict(list)
    for (i, j, k) in tris:
        angs = {
            i: vertex_angle(coords, i, j, k),
            j: vertex_angle(coords, j, k, i),
            k: vertex_angle(coords, k, i, j),
        }
        tmin = min(angs.values())
        for node, ang in angs.items():
            vertex_angles[node].append(ang)
            tri_min_angle_at[node].append(tmin)
    for n in range(num_nodes):
        vertex_angles.setdefault(n, [])
        tri_min_angle_at.setdefault(n, [])
    return vertex_angles, tri_min_angle_at


def run_bias(bias, nodes, tris, sig_frac, csv_rows):
    coords = {n: (nodes[n]["x"], nodes[n]["y"]) for n in nodes}
    num_nodes = len(nodes)
    vertex_angles, tri_min_angle_at = build_node_angles(tris, num_nodes, coords)

    xs = [coords[n][0] for n in nodes]
    ys = [coords[n][1] for n in nodes]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)

    def is_boundary(n):
        x, y = coords[n]
        return (
            abs(x - xmin) < 1e-9
            or abs(x - xmax) < 1e-9
            or abs(y - ymin) < 1e-9
            or abs(y - ymax) < 1e-9
        )

    mags = {}
    for n in nodes:
        es = math.hypot(nodes[n].get("ex_s", float("nan")), nodes[n].get("ey_s", float("nan")))
        mags[n] = es
    peak = max((m for m in mags.values() if math.isfinite(m)), default=0.0)
    gate = sig_frac * peak

    records = []
    for n in nodes:
        exs, eys = nodes[n].get("ex_s"), nodes[n].get("ey_s")
        exv, eyv = nodes[n].get("ex_v"), nodes[n].get("ey_v")
        if any(v is None or (isinstance(v, float) and math.isnan(v)) for v in (exs, eys, exv, eyv)):
            continue
        theta_s = math.degrees(math.atan2(eys, exs))
        theta_v = math.degrees(math.atan2(eyv, exv))
        dtheta = wrap_deg(theta_v - theta_s)
        va = vertex_angles[n]
        rec = {
            "bias_V": bias,
            "node_id": n,
            "x_um": coords[n][0],
            "y_um": coords[n][1],
            "boundary": int(is_boundary(n)),
            "emag_s": mags[n],
            "significant": int(mags[n] >= gate),
            "theta_sentaurus_deg": theta_s,
            "theta_vela_deg": theta_v,
            "dtheta_deg": dtheta,
            "abs_dtheta_deg": abs(dtheta),
            "vertex_angle_min": min(va) if va else float("nan"),
            "vertex_angle_max": max(va) if va else float("nan"),
            "angle_sum": sum(va) if va else float("nan"),
            "tri_min_angle": min(tri_min_angle_at[n]) if tri_min_angle_at[n] else float("nan"),
        }
        records.append(rec)
        csv_rows.append(rec)

    return records, is_boundary, peak


def summarize(records, label):
    sig = [r for r in records if r["significant"]]
    if not sig:
        print(f"    {label}: (no significant-field nodes)")
        return
    absd = [r["abs_dtheta_deg"] for r in sig]
    print(
        f"    {label:<10} n={len(sig):<3} "
        f"|dtheta| mean={sum(absd)/len(absd):7.3f}  max={max(absd):7.3f} deg"
    )
    for metric in ("vertex_angle_min", "vertex_angle_max", "angle_sum", "tri_min_angle"):
        xs = [r[metric] for r in sig]
        r = pearson(xs, absd)
        print(f"        corr(|dtheta|, {metric:<16}) = {r:+.3f}")


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
    parser.add_argument("--bias", type=float, default=-5.0)
    parser.add_argument("--all-biases", action="store_true")
    parser.add_argument("--sig-frac", type=float, default=0.01)
    parser.add_argument(
        "--out",
        default=os.path.join("build", "diagnostics", "efield_angle_vs_interior_angle.csv"),
    )
    args = parser.parse_args()

    tris = parse_elements(args.elements)
    data = parse_aligned(args.aligned)

    if args.all_biases:
        biases = sorted(data.keys())
    else:
        biases = [b for b in sorted(data.keys()) if abs(b - args.bias) < 1e-6] or [args.bias]

    csv_rows = []
    print(f"Triangles       : {len(tris)}")
    print(f"Significance gate: |E| >= {args.sig_frac:g} * peak|E| per bias")
    print("theta = atan2(E_y, E_x); dtheta = wrap(theta_vela - theta_sentaurus)\n")

    for bias in biases:
        nodes = data.get(bias)
        if not nodes:
            continue
        records, is_boundary, peak = run_bias(bias, nodes, tris, args.sig_frac, csv_rows)
        print("=" * 74)
        print(f"bias = {bias:g} V   peak|E| = {peak:.4g} V/cm")
        summarize(records, "all")
        summarize([r for r in records if r["boundary"]], "boundary")
        summarize([r for r in records if not r["boundary"]], "interior")

        # top angular-deviation nodes (significant only)
        top = sorted(
            (r for r in records if r["significant"]),
            key=lambda r: r["abs_dtheta_deg"],
            reverse=True,
        )[:6]
        if top:
            print("    top |dtheta| significant nodes:")
            print(
                f"      {'node':>4}{'x':>7}{'y':>7}{'bnd':>4}"
                f"{'|dtheta|':>10}{'vtx_min':>9}{'vtx_max':>9}{'ang_sum':>9}{'tri_min':>9}"
            )
            for r in top:
                print(
                    f"      {r['node_id']:>4}{r['x_um']:>7.3f}{r['y_um']:>7.3f}"
                    f"{r['boundary']:>4}{r['abs_dtheta_deg']:>10.3f}"
                    f"{r['vertex_angle_min']:>9.2f}{r['vertex_angle_max']:>9.2f}"
                    f"{r['angle_sum']:>9.2f}{r['tri_min_angle']:>9.2f}"
                )

    if csv_rows:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0].keys()))
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"\nPer-node table written to {args.out} ({len(csv_rows)} rows)")


if __name__ == "__main__":
    main()
