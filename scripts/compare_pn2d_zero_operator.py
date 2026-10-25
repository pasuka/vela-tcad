#!/usr/bin/env python3
"""Compare an exported Vela augmented system with one native debug iteration.

Diagnostic only. Requires NumPy/SciPy. State populations must already match;
this script changes coordinates and eliminates constraints, never fits physics
parameters or alters acceptance gates. See the PN2D zero-operator report.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from compare_sentaurus_debug_matrix import iterations, super_iterations, super_closure
from constraint_elimination import eliminate_constraints
from native_debug_precision import variable_radii, propagate_state_rounding


def selected(stream, index):
    try:
        return next(itertools.islice(stream, index, None))
    except StopIteration as exc:
        raise ValueError(f"native iteration {index} not found") from exc


def canonical_contacts(matrix, rhs, contacts):
    """Use retained contact-potential constraints as invertible row operations."""
    matrix = matrix.tocsc()
    block = matrix[contacts, :][:, contacts]
    remainder = block - sparse.diags(block.diagonal())
    remainder.eliminate_zeros()
    if remainder.nnz or np.any(block.diagonal() == 0):
        raise ValueError("contact-potential constraint block must be diagonal")
    inverse = sparse.diags(1 / block.diagonal())
    embed = sparse.coo_matrix(
        (np.ones(len(contacts)), (contacts, np.arange(len(contacts)))),
        shape=(matrix.shape[0], len(contacts)),
    ).tocsc()
    operation = embed - matrix[:, contacts]
    return ((matrix + operation @ inverse @ matrix[contacts, :]).tocsr(),
            rhs + operation @ (inverse @ rhs[contacts]))


def equation_rhs_difference(vela_rhs, native_rhs, native_diagonal, row_ratios):
    """Use independent equation units, not each candidate Jacobian diagonal."""
    return (row_ratios * vela_rhs - native_rhs) / native_diagonal


def compare(vela_dir, debug, super_debug, index, populations, potential_v, density_cm3,
            matrix_output=None, row_scale_ratios=None):
    meta = json.loads((vela_dir / "identity.json").read_text())
    n = meta["nodes"]
    d = 3 * n + 1
    contacts = np.asarray(sorted(set(meta["contact_nodes"])), dtype=int)
    population = np.asarray(json.loads(populations.read_text()), dtype=float)
    if (population.shape != (n, 4) or not np.array_equal(population[:, 0], np.arange(n))
            or np.any(population[:, 2:] <= 0)):
        raise ValueError("populations must be ordered [node, psi_V, n_cm3, p_cm3]")
    with debug.open() as stream:
        labelled = selected(iterations(stream), index)
    native = selected(super_iterations(super_debug.read_text()), index)
    potential_error = max(abs(population[node, 1] - value * potential_v)
                          for (block, node), value in labelled["variables"].items()
                          if block == 0)
    population_error = max(abs(population[node, block] / (value * density_cm3) - 1)
                           for (block, node), value in labelled["variables"].items()
                           if block in (2, 3))
    if potential_error > 1e-12 or population_error > 1e-12:
        raise ValueError("exported population state does not match the selected native VAR/scales")
    actual = np.genfromtxt(vela_dir / "physical_state.csv", delimiter=",", names=True)
    assembled_psi_error = float(max(abs(actual["psi_V"] - population[:, 1])))
    assembled_density_error = float(max(
        max(abs(actual["n_cm3"] / population[:, 2] - 1)),
        max(abs(actual["p_cm3"] / population[:, 3] - 1))))
    if assembled_psi_error > 1e-12 or assembled_density_error > 1e-12:
        raise ValueError("Vela actually assembled a different physical state")
    free_e = sorted(node for block, node in labelled["variables"] if block == 2)
    free_h = sorted(node for block, node in labelled["variables"] if block == 3)
    if free_e != free_h or free_e != sorted(set(range(n)) - set(contacts)):
        raise ValueError("native free-carrier nodes disagree with Vela contacts")
    dimension = n + 1 + len(free_e) + len(free_h)
    if len(native["rhs"]) != dimension:
        raise ValueError("native internal dimension does not match device plus circuit")
    circuit = [value for (block, _), value in labelled["solution"].items() if block == 1]
    if not circuit or not np.allclose(circuit, circuit[0], rtol=1e-12, atol=1e-25):
        raise ValueError("native electrode aliases are not one independent circuit variable")
    native_direction = np.asarray(native["solution"])
    labelled_direction = np.asarray(
        [labelled["solution"][(0, i)] for i in range(n)] + [circuit[0]]
        + [labelled["solution"][(2, i)] for i in free_e]
        + [labelled["solution"][(3, i)] for i in free_h])
    order_error = np.max(abs(native_direction - labelled_direction)) / max(
        np.max(abs(native_direction)), np.finfo(float).tiny)
    if order_error > 1e-10:
        raise ValueError("native internal ordering is not verified by labelled increments")

    entries = np.genfromtxt(vela_dir / "augmented.csv", delimiter=",", names=True)
    rhs = np.genfromtxt(vela_dir / "rhs.csv", delimiter=",", names=True)["rhs"]
    matrix = sparse.coo_matrix(
        (entries["value"], (entries["row"].astype(int), entries["col"].astype(int))),
        shape=(d, d)).tocsc()
    # Full population-coordinate transform BEFORE eliminating contact carriers.
    vt_ratio = potential_v / meta["potential_scale_V"]
    nodes = np.arange(n)
    rows = np.r_[nodes, n + nodes, 2*n + nodes, n + nodes, 2*n + nodes, 3*n]
    cols = np.r_[nodes, nodes, nodes, n + nodes, 2*n + nodes, 3*n]
    values = np.r_[np.full(3*n, vt_ratio), -density_cm3 / population[:, 2],
                   density_cm3 / population[:, 3], potential_v]
    transform = sparse.coo_matrix((values, (rows, cols)), shape=(d, d)).tocsc()
    eliminated = np.r_[n + contacts, 2*n + contacts]
    reduced, reduced_rhs, kept, _ = eliminate_constraints(matrix @ transform, rhs, eliminated)
    lookup = {int(node): i for i, node in enumerate(kept)}
    order = [lookup[i] for i in list(range(n)) + [3*n]
             + [n+i for i in free_e] + [2*n+i for i in free_h]]
    vela, vb = canonical_contacts(reduced[order, :][:, order], reduced_rhs[order], contacts)
    entries = native["matrix"]
    nm = sparse.coo_matrix(
        (list(entries.values()), ([r for r, _ in entries], [c for _, c in entries])),
        shape=(dimension, dimension)).tocsc()
    nm, nb = canonical_contacts(nm, np.asarray(native["rhs"]), contacts)
    vd, nd = vela.diagonal(), nm.diagonal()
    if matrix_output is not None:
        matrix_output.mkdir(parents=True, exist_ok=True)
        sparse.save_npz(matrix_output / "vela.npz", vela)
        sparse.save_npz(matrix_output / "native.npz", nm)
        np.savez(matrix_output / "vectors.npz", vela_rhs=vb, native_rhs=nb)
    row_ratios = None
    if row_scale_ratios is not None:
        row_ratios = np.load(row_scale_ratios)
        if (row_ratios.shape != (dimension,) or not np.all(np.isfinite(row_ratios))
                or np.any(row_ratios == 0)):
            raise ValueError("invalid independently calibrated equation row scales")
    if np.any(vd == 0) or np.any(nd == 0):
        raise ValueError("diagonal normalization requires nonzero diagonals")
    vn, nn = sparse.diags(1/vd) @ vela, sparse.diags(1/nd) @ nm
    with debug.open() as stream:
        radii = variable_radii(stream, index)
    circuit_key = next(key for key in radii if key[0] == 1)
    ordered_keys = ([(0, i) for i in range(n)] + [circuit_key]
                    + [(2, i) for i in free_e] + [(3, i) for i in free_h])
    rounding_bound = propagate_state_rounding(nn, [radii[key] for key in ordered_keys])
    rhs_delta = abs(vb/vd - nb/nd)
    row_norm = lambda a: np.sqrt(np.asarray(a.multiply(a).sum(axis=1)).ravel())
    error = row_norm(vn - nn) / row_norm(nn)
    blocks = np.r_[np.zeros(n, int), 1, np.full(len(free_e), 2), np.full(len(free_h), 3)]
    denominator = abs(vela) @ abs(native_direction) + abs(vb)
    backward = np.divide(abs(vela @ native_direction - vb), denominator,
                         out=np.zeros(dimension), where=denominator > 0)
    result = {"iteration": index, "dimension": dimension,
              "row_scale_calibration_sha256": (hashlib.sha256(row_scale_ratios.read_bytes()).hexdigest()
                                                if row_scale_ratios is not None else None),
              "state_potential_mapping_error_V": float(potential_error),
              "state_population_mapping_relative_error": float(population_error),
              "assembled_potential_error_V": assembled_psi_error,
              "assembled_population_relative_error": assembled_density_error,
              "internal_order_relative_error": float(order_error),
              "native_linear_closure": super_closure(native), "blocks": {}}
    for block, name in enumerate(["psi", "circuit", "electron", "hole"]):
        sel = blocks == block
        native_norm = float(np.linalg.norm(nb[sel] / nd[sel]))
        difference = float(np.linalg.norm(vb[sel] / vd[sel] - nb[sel] / nd[sel]))
        result["blocks"][name] = {
            "VAR_rounding_first_order_radius_L2": float(np.linalg.norm(rounding_bound[sel])),
            "rhs_rows_outside_VAR_rounding_first_order_radius": int(np.count_nonzero(
                rhs_delta[sel] > rounding_bound[sel])),
            "matrix_relative_row_L2_max": float(max(error[sel])),
            "matrix_relative_row_L2_median": float(np.median(error[sel])),
            "rhs_native_diag_normalized_L2": native_norm,
            "rhs_difference_diag_normalized_L2": difference,
            # A ratio to roundoff-sized RHS has no acceptance meaning; retain
            # the actual denominator and use null only for an exactly zero RHS.
            "rhs_relative_L2": difference / native_norm if native_norm else None,
            "native_direction_vela_max_backward_error": float(max(backward[sel])),
        }
        if row_ratios is not None:
            # Own-diagonal RHS normalization mixes operator and RHS errors when
            # the two Jacobians differ. Use independently calibrated equation
            # scales for a physical RHS comparison in that case.
            common = equation_rhs_difference(vb[sel], nb[sel], nd[sel], row_ratios[sel])
            result["blocks"][name]["rhs_common_row_scale_difference_L2"] = float(np.linalg.norm(common))
            result["blocks"][name]["rhs_common_row_scale_relative_L2"] = (
                float(np.linalg.norm(common)) / native_norm if native_norm else None)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-dir", required=True, type=Path)
    parser.add_argument("--debug", required=True, type=Path)
    parser.add_argument("--super-debug", required=True, type=Path)
    parser.add_argument("--index", required=True, type=int)
    parser.add_argument("--populations", required=True, type=Path)
    parser.add_argument("--potential-scale-V", required=True, type=float)
    parser.add_argument("--density-scale-cm3", required=True, type=float)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--matrix-output-dir", type=Path)
    parser.add_argument("--row-scale-ratios", type=Path)
    args = parser.parse_args()
    if args.index < 0 or args.potential_scale_V <= 0 or args.density_scale_cm3 <= 0:
        parser.error("iteration and scales must be valid")
    result = compare(args.vela_dir, args.debug, args.super_debug, args.index,
                     args.populations, args.potential_scale_V, args.density_scale_cm3,
                     args.matrix_output_dir, args.row_scale_ratios)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
