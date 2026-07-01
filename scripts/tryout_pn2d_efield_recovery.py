#!/usr/bin/env python3
"""Try-out reconstruction of nodal electric field on the pn2d coarse7x3 mesh.

This standalone script compares three node electric-field recovery operators
against the Sentaurus reference nodal ``ElectricField`` using only the two CSV
inputs that ship with the coarse7x3 comparison report:

* ``elements.csv``        -- triangle connectivity (node0, node1, node2)
* ``coarse_node_field_compare_aligned.csv`` -- per-node coordinates, the
  full-precision Sentaurus nodal potential, the Sentaurus reference field
  vector, and Vela's current least-squares field vector.

Operators evaluated (all fed the SAME nodal potential so the comparison
isolates the reconstruction operator, not any solution difference):

* ``LS_node``  -- Vela's current node-neighbor weighted (1/d) least-squares
                  gradient. Reproduced here as the baseline.
* ``Method1``  -- area-weighted average of the exact P1 constant per-cell
                  gradients (Sentaurus box-method style).
* ``Method2``  -- superconvergent patch recovery (SPR / Zienkiewicz-Zhu): a
                  linear field fitted to the incident cell gradients and
                  evaluated at the node, with area-average fallback.

Units: potential [V], coordinates [um], field [V/cm]. A potential gradient in
V/um converts to V/cm by a factor of 1e4 (1 um = 1e-4 cm).

Pure standard library; no numpy required.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict

# um-gradient (V/um) -> V/cm
GRAD_TO_V_PER_CM = 1.0e4
COORD_TOL = 1.0e-9


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------
def parse_elements(path):
    tris = []
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            tris.append((int(row["node0"]), int(row["node1"]), int(row["node2"])))
    return tris


def _norm_zero(value):
    # Collapse -0.0 (and tiny noise) to 0.0 so boundary detection is clean.
    return 0.0 if abs(value) < COORD_TOL else value


def parse_aligned(path):
    """Return {bias: {node_id: {x, y, psi_s, psi_v, ex_s, ey_s, ex_v, ey_v}}}."""
    data = defaultdict(lambda: defaultdict(dict))
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = [col.strip() for col in next(reader)]
        idx = {name: i for i, name in enumerate(header)}

        def cell(row, name):
            return row[idx[name]].strip()

        def fnum(text):
            try:
                return float(text)
            except ValueError:
                return float("nan")

        for row in reader:
            if not row or len(row) < len(header):
                continue
            bias = float(cell(row, "bias_V"))
            nid = int(cell(row, "node_id"))
            quantity = cell(row, "quantity")
            entry = data[bias][nid]
            entry["x"] = _norm_zero(float(cell(row, "x_um")))
            entry["y"] = _norm_zero(float(cell(row, "y_um")))
            s_val = fnum(cell(row, "sentaurus_value"))
            v_val = fnum(cell(row, "vela_value_scaled_to_sentaurus_units"))
            if quantity == "potential":
                entry["psi_s"] = s_val
                entry["psi_v"] = v_val
            elif quantity == "electric_field_x":
                entry["ex_s"] = s_val
                entry["ex_v"] = v_val
            elif quantity == "electric_field_y":
                entry["ey_s"] = s_val
                entry["ey_v"] = v_val
    return data


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------
def cell_gradient(coords, psi, tri):
    """Return (Ex, Ey, area) for a P1 triangle; field E = -grad(psi) in V/cm."""
    i, j, k = tri
    x0, y0 = coords[i]
    x1, y1 = coords[j]
    x2, y2 = coords[k]
    dx1, dy1 = x1 - x0, y1 - y0
    dx2, dy2 = x2 - x0, y2 - y0
    det = dx1 * dy2 - dy1 * dx2
    area = 0.5 * abs(det)
    if area <= 0.0 or not math.isfinite(det):
        return None
    dv1 = psi[j] - psi[i]
    dv2 = psi[k] - psi[i]
    grad_x = (dv1 * dy2 - dy1 * dv2) / det
    grad_y = (dx1 * dv2 - dv1 * dx2) / det
    ex = -grad_x * GRAD_TO_V_PER_CM
    ey = -grad_y * GRAD_TO_V_PER_CM
    return ex, ey, area


def solve3(matrix, rhs):
    """Solve a symmetric 3x3 system by Cramer's rule; None if near-singular."""
    a = matrix

    def det3(m):
        return (
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        )

    base = det3(a)
    if abs(base) <= 1.0e-30 or not math.isfinite(base):
        return None
    out = []
    for col in range(3):
        m = [list(row) for row in a]
        for r in range(3):
            m[r][col] = rhs[r]
        out.append(det3(m) / base)
    return out


# ---------------------------------------------------------------------------
# Recovery operators (all take a nodal potential array `psi`)
# ---------------------------------------------------------------------------
def recover_ls_node(coords, psi, neighbors):
    """Vela's current node-neighbor weighted (1/d) least-squares gradient."""
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
            dv = psi[other] - psi[node]
            sxx += weight * dx * dx
            sxy += weight * dx * dy
            syy += weight * dy * dy
            sxv += weight * dx * dv
            syv += weight * dy * dv
        det = sxx * syy - sxy * sxy
        if abs(det) <= 1.0e-30:
            fields[node] = None
            continue
        grad_x = (sxv * syy - syv * sxy) / det
        grad_y = (sxx * syv - sxy * sxv) / det
        fields[node] = (-grad_x * GRAD_TO_V_PER_CM, -grad_y * GRAD_TO_V_PER_CM)
    return fields


def recover_area_average(coords, psi, node_cells, tris):
    """Method 1: area-weighted average of constant per-cell fields."""
    cell_fields = {c: cell_gradient(coords, psi, tri) for c, tri in enumerate(tris)}
    fields = {}
    for node, cells in node_cells.items():
        ex = ey = total = 0.0
        for c in cells:
            cf = cell_fields[c]
            if cf is None:
                continue
            fex, fey, area = cf
            ex += area * fex
            ey += area * fey
            total += area
        fields[node] = (ex / total, ey / total) if total > 0.0 else None
    return fields, cell_fields


def fit_patch_polynomial(coords, tris, cell_fields, cells):
    """Fit a linear field E(x,y) = a0 + a1*x + a2*y to the incident cell fields.

    Returns (coeff_x, coeff_y) or None if fewer than 3 samples / singular.
    """
    rows = []
    for c in cells:
        cf = cell_fields[c]
        if cf is None:
            continue
        i, j, k = tris[c]
        cx = (coords[i][0] + coords[j][0] + coords[k][0]) / 3.0
        cy = (coords[i][1] + coords[j][1] + coords[k][1]) / 3.0
        rows.append((cx, cy, cf[0], cf[1]))
    if len(rows) < 3:
        return None
    ata = [[0.0] * 3 for _ in range(3)]
    atbx = [0.0, 0.0, 0.0]
    atby = [0.0, 0.0, 0.0]
    for cx, cy, fex, fey in rows:
        p = (1.0, cx, cy)
        for r in range(3):
            for col in range(3):
                ata[r][col] += p[r] * p[col]
            atbx[r] += p[r] * fex
            atby[r] += p[r] * fey
    coeff_x = solve3(ata, atbx)
    coeff_y = solve3(ata, atby)
    if coeff_x is None or coeff_y is None:
        return None
    return coeff_x, coeff_y


def _eval_poly(poly, x, y):
    coeff_x, coeff_y = poly
    return (
        coeff_x[0] + coeff_x[1] * x + coeff_x[2] * y,
        coeff_y[0] + coeff_y[1] * x + coeff_y[2] * y,
    )


def recover_direction_split(coords, psi, neighbors):
    """Method 4: axis-decoupled 1D differencing on a structured grid.

    E_x is fitted with a 1D weighted least squares using ONLY the horizontal
    (same-y) edge neighbors, and E_y using ONLY the vertical (same-x) edge
    neighbors. Because the two directions never share a stencil, x-curvature
    can no longer alias into E_y (and vice versa). For an interior node this
    reduces to a central difference along each axis; on a boundary it becomes a
    one-sided difference along the missing direction.
    """
    fields = {}
    for node, nb in neighbors.items():
        cx, cy = coords[node]
        sxx = sxv = 0.0
        syy = syv = 0.0
        for other in nb:
            ox, oy = coords[other]
            dx, dy = ox - cx, oy - cy
            dv = psi[other] - psi[node]
            if abs(dy) < COORD_TOL and abs(dx) > COORD_TOL:  # horizontal neighbor
                w = 1.0 / abs(dx)
                sxx += w * dx * dx
                sxv += w * dx * dv
            elif abs(dx) < COORD_TOL and abs(dy) > COORD_TOL:  # vertical neighbor
                w = 1.0 / abs(dy)
                syy += w * dy * dy
                syv += w * dy * dv
        grad_x = sxv / sxx if sxx > 0.0 else 0.0
        grad_y = syv / syy if syy > 0.0 else 0.0
        fields[node] = (-grad_x * GRAD_TO_V_PER_CM, -grad_y * GRAD_TO_V_PER_CM)
    return fields


def recover_spr(coords, psi, node_cells, tris, cell_fields, area_fields):
    """Method 2: superconvergent patch recovery of the cell fields."""
    fields = {}
    for node, cells in node_cells.items():
        poly = fit_patch_polynomial(coords, tris, cell_fields, cells)
        if poly is None:
            fields[node] = area_fields.get(node)
            continue
        nx, ny = coords[node]
        fields[node] = _eval_poly(poly, nx, ny)
    return fields


def recover_spr_boundary_aware(
    coords, node_cells, tris, cell_fields, area_fields, neighbors, is_boundary
):
    """Method 3: boundary-aware SPR (ZZ boundary treatment).

    Interior nodes use their own patch polynomial (already exact for a linear
    field). Boundary/corner nodes -- whose own patch has too few, one-sided
    cells and therefore extrapolates badly -- instead borrow the recovered
    polynomial from each adjacent interior node's patch and evaluate it at the
    boundary node location, averaging over the adjacent interior nodes.
    """
    node_poly = {
        node: fit_patch_polynomial(coords, tris, cell_fields, cells)
        for node, cells in node_cells.items()
    }

    fields = {}
    for node in node_cells:
        nx, ny = coords[node]
        if not is_boundary(node):
            fields[node] = (
                _eval_poly(node_poly[node], nx, ny)
                if node_poly[node] is not None
                else area_fields.get(node)
            )
            continue

        # Boundary node: borrow adjacent interior-node polynomials.
        donors = [
            node_poly[nb]
            for nb in neighbors[node]
            if not is_boundary(nb) and node_poly[nb] is not None
        ]
        if donors:
            ex = sum(_eval_poly(p, nx, ny)[0] for p in donors) / len(donors)
            ey = sum(_eval_poly(p, nx, ny)[1] for p in donors) / len(donors)
            fields[node] = (ex, ey)
        elif node_poly[node] is not None:
            fields[node] = _eval_poly(node_poly[node], nx, ny)
        else:
            fields[node] = area_fields.get(node)
    return fields


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------
def vector_error(field, ref):
    if field is None:
        return None
    ex, ey = field
    rx, ry = ref
    return math.hypot(ex - rx, ey - ry)


def weighted_rms_rel_error(errors, mags):
    """L2 field error normalized by the L2 of the reference magnitudes."""
    num = sum(e * e for e in errors)
    den = sum(m * m for m in mags)
    return math.sqrt(num / den) if den > 0.0 else float("nan")


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def build_topology(tris, num_nodes):
    node_cells = defaultdict(list)
    neighbors = defaultdict(set)
    for c, (i, j, k) in enumerate(tris):
        for node in (i, j, k):
            node_cells[node].append(c)
        for a in (i, j, k):
            for b in (i, j, k):
                if a != b:
                    neighbors[a].add(b)
    for node in range(num_nodes):
        node_cells.setdefault(node, [])
        neighbors.setdefault(node, set())
    return node_cells, neighbors


def run_bias(bias, nodes, tris, potential_key, sig_frac, csv_rows):
    num_nodes = len(nodes)
    coords = {n: (nodes[n]["x"], nodes[n]["y"]) for n in nodes}
    psi = {n: nodes[n][potential_key] for n in nodes}
    node_cells, neighbors = build_topology(tris, num_nodes)

    ls = recover_ls_node(coords, psi, neighbors)
    area, cell_fields = recover_area_average(coords, psi, node_cells, tris)

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

    spr = recover_spr(coords, psi, node_cells, tris, cell_fields, area)
    spr_bd = recover_spr_boundary_aware(
        coords, node_cells, tris, cell_fields, area, neighbors, is_boundary
    )
    dirsplit = recover_direction_split(coords, psi, neighbors)

    ref = {n: (nodes[n]["ex_s"], nodes[n]["ey_s"]) for n in nodes}
    ref_mag = {n: math.hypot(*ref[n]) for n in nodes}
    max_ref = max(ref_mag.values()) if ref_mag else 0.0
    sig_threshold = sig_frac * max_ref

    methods = {
        "LS_node": ls,
        "Method1_area": area,
        "Method2_spr": spr,
        "Method3_spr_bd": spr_bd,
        "Method4_dirsplit": dirsplit,
    }

    # Physical-truth surrogate: this is a 1D device (field is constant along y),
    # so the true field at a column equals its interior-row node value with
    # E_y == 0. Interior nodes reproduce Sentaurus exactly, so use them.
    col_interior = {}
    for n in nodes:
        if not is_boundary(n):
            col_interior[round(coords[n][0], 9)] = (nodes[n]["ex_s"], 0.0)
    truth = {n: col_interior.get(round(coords[n][0], 9)) for n in nodes}

    # --- global field-normalized L2 error, split interior/boundary ---
    stats = {}
    for name, field in methods.items():
        for scope in ("all", "interior", "boundary"):
            errs, mags = [], []
            for n in nodes:
                if scope == "interior" and is_boundary(n):
                    continue
                if scope == "boundary" and not is_boundary(n):
                    continue
                err = vector_error(field[n], ref[n])
                if err is None:
                    continue
                errs.append(err)
                mags.append(ref_mag[n])
            stats[(name, scope)] = weighted_rms_rel_error(errs, mags)

    # --- boundary error vs the physical 1D truth (E_y == 0) ---
    # Include Sentaurus's own reference field as a "method" to expose that it
    # also carries a boundary artifact relative to the 1D truth.
    truth_stats = {}
    for name, field in list(methods.items()) + [("Sentaurus_ref", ref)]:
        errs, mags = [], []
        for n in nodes:
            if not is_boundary(n) or truth[n] is None:
                continue
            f = field[n]
            if f is None:
                continue
            tx, ty = truth[n]
            errs.append(math.hypot(f[0] - tx, f[1] - ty))
            mags.append(math.hypot(tx, ty))
        truth_stats[name] = weighted_rms_rel_error(errs, mags)

    # --- validation: LS reconstructed here vs the CSV's stored Vela field ---
    val_max = 0.0
    for n in nodes:
        if ls[n] is None or math.isnan(nodes[n].get("ex_v", float("nan"))):
            continue
        d = math.hypot(ls[n][0] - nodes[n]["ex_v"], ls[n][1] - nodes[n]["ey_v"])
        val_max = max(val_max, d)

    # --- per-node table for the significant boundary nodes ---
    boundary_sig = sorted(
        (n for n in nodes if is_boundary(n) and ref_mag[n] >= sig_threshold),
        key=lambda n: ref_mag[n],
        reverse=True,
    )

    for n in nodes:
        row = {
            "bias_V": bias,
            "node_id": n,
            "x_um": coords[n][0],
            "y_um": coords[n][1],
            "boundary": int(is_boundary(n)),
            "ref_ex": ref[n][0],
            "ref_ey": ref[n][1],
            "ref_mag": ref_mag[n],
        }
        for name, field in methods.items():
            f = field[n]
            row[f"{name}_ex"] = "" if f is None else f[0]
            row[f"{name}_ey"] = "" if f is None else f[1]
            err = vector_error(f, ref[n])
            row[f"{name}_abs_err"] = "" if err is None else err
            row[f"{name}_rel_err"] = (
                "" if err is None or ref_mag[n] <= 0 else err / ref_mag[n]
            )
        csv_rows.append(row)

    return stats, truth_stats, val_max, boundary_sig, methods, ref, ref_mag, coords, is_boundary, truth


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
        default="sentaurus",
        help="Nodal potential fed to all operators (default: full-precision Sentaurus).",
    )
    parser.add_argument(
        "--sig-frac",
        type=float,
        default=0.01,
        help="Significant-field threshold as a fraction of the peak |E| per bias.",
    )
    parser.add_argument("--bias", type=float, default=None, help="Restrict to one bias.")
    parser.add_argument(
        "--out",
        default=os.path.join("build", "diagnostics", "efield_recovery_tryout.csv"),
    )
    args = parser.parse_args()

    tris = parse_elements(args.elements)
    data = parse_aligned(args.aligned)
    potential_key = "psi_s" if args.potential == "sentaurus" else "psi_v"

    biases = sorted(data.keys()) if args.bias is None else [args.bias]
    csv_rows = []

    print(f"Potential input : {args.potential} ({potential_key})")
    print(f"Triangles       : {len(tris)}")
    print(f"Significant band : |E| >= {args.sig_frac:.3g} * peak(|E|) per bias\n")

    for bias in biases:
        nodes = data.get(bias)
        if not nodes:
            continue
        (
            stats,
            truth_stats,
            val_max,
            boundary_sig,
            methods,
            ref,
            ref_mag,
            coords,
            is_boundary,
            truth,
        ) = run_bias(bias, nodes, tris, potential_key, args.sig_frac, csv_rows)

        print("=" * 78)
        print(f"bias = {bias:g} V   nodes = {len(nodes)}   peak|E| = {max(ref_mag.values()):.4g} V/cm")
        print(f"  (validation) max |LS_here - CSV Vela field| = {val_max:.4g} V/cm")
        print("  field-normalized L2 error vs Sentaurus reference:")
        print(f"    {'method':<18}{'all':>12}{'interior':>12}{'boundary':>12}")
        for name in (
            "LS_node",
            "Method1_area",
            "Method2_spr",
            "Method3_spr_bd",
            "Method4_dirsplit",
        ):
            print(
                f"    {name:<18}"
                f"{stats[(name, 'all')]:>12.4%}"
                f"{stats[(name, 'interior')]:>12.4%}"
                f"{stats[(name, 'boundary')]:>12.4%}"
            )

        print("\n  boundary L2 error vs PHYSICAL 1D truth (E_y == 0):")
        for name in (
            "Sentaurus_ref",
            "LS_node",
            "Method1_area",
            "Method2_spr",
            "Method3_spr_bd",
            "Method4_dirsplit",
        ):
            print(f"    {name:<18}{truth_stats[name]:>12.4%}")

        if boundary_sig:
            print("\n  significant boundary nodes -- transverse E_y [V/cm] (true E_y = 0):")
            print(
                f"    {'node':>4}{'(x,y)':>13}{'|E_ref|':>11}"
                f"{'Ey_sede':>11}{'Ey_LS':>11}{'Ey_dsplit':>11}"
            )
            for n in boundary_sig:
                x, y = coords[n]
                ey_ls = methods["LS_node"][n]
                ey_ds = methods["Method4_dirsplit"][n]
                print(
                    f"    {n:>4}{f'({x:g},{y:g})':>13}{ref_mag[n]:>11.4g}"
                    f"{ref[n][1]:>11.4g}"
                    f"{(ey_ls[1] if ey_ls else float('nan')):>11.4g}"
                    f"{(ey_ds[1] if ey_ds else float('nan')):>11.4g}"
                )

            print("\n  significant boundary nodes -- vector error vs Sentaurus / vs truth:")
            print(
                f"    {'node':>4}{'(x,y)':>13}"
                f"{'LS/sede':>10}{'dsplit/sede':>13}{'LS/truth':>11}{'dsplit/truth':>14}"
            )
            for n in boundary_sig:
                x, y = coords[n]

                def rel_to(field, target):
                    f = field[n]
                    t = target[n] if isinstance(target, dict) else target
                    if f is None or t is None:
                        return float("nan")
                    mag = math.hypot(*t)
                    return float("nan") if mag <= 0 else math.hypot(f[0] - t[0], f[1] - t[1]) / mag

                print(
                    f"    {n:>4}{f'({x:g},{y:g})':>13}"
                    f"{rel_to(methods['LS_node'], ref):>10.2%}"
                    f"{rel_to(methods['Method4_dirsplit'], ref):>13.2%}"
                    f"{rel_to(methods['LS_node'], truth):>11.2%}"
                    f"{rel_to(methods['Method4_dirsplit'], truth):>14.2%}"
                )
        print()

    if csv_rows:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        fieldnames = list(csv_rows[0].keys())
        with open(args.out, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)
        print(f"Per-node detail written to: {args.out}")


if __name__ == "__main__":
    main()
