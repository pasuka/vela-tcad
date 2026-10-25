#!/usr/bin/env python3
"""Compare independent Poisson then coupled zero-bias fields by mesh node ID."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy import sparse

from compare_sentaurus_debug_matrix import iterations, super_iterations


def field_errors(candidate, reference):
    if (candidate.shape != reference.shape or not np.all(np.isfinite(candidate))
            or not np.all(np.isfinite(reference))):
        raise ValueError("matching finite field arrays required")
    delta = candidate - reference
    return {"max_absolute": float(max(abs(delta))),
            "rms_absolute": float(np.sqrt(np.mean(delta**2))),
            "worst_node": int(np.argmax(abs(delta)))}


def density_errors(candidate, reference):
    if np.any(candidate <= 0) or np.any(reference <= 0):
        raise ValueError("positive carrier densities required")
    raw = field_errors(np.log10(candidate), np.log10(reference))
    result = {"max_log10_error": raw["max_absolute"],
              "rms_log10_error": raw["rms_absolute"], "worst_node": raw["worst_node"]}
    result["max_relative"] = float(max(abs(candidate / reference - 1)))
    return result


def compare(vela_dir, debug, super_debug, native_fields, vt):
    with debug.open() as stream:
        trace = list(iterations(stream))
    native_systems = list(super_iterations(super_debug.read_text()))
    first_dd = next(i for i, step in enumerate(trace) if len(step["names"]) > 1)
    if first_dd == 0 or any(len(s["names"]) != 1 for s in trace[:first_dd]):
        raise ValueError("native input must solve Poisson before the first coupled problem")
    field = lambda name: np.genfromtxt(native_fields / (name + "_region0.csv"),
                                       delimiter=",", names=True)["component0"]
    # All potentials are present before the first coupled Newton update.
    n = len(field("ElectrostaticPotential"))
    poisson_reference = np.array([trace[first_dd]["variables"][(0, i)] * vt for i in range(n)])
    final_reference = {"psi_V": field("ElectrostaticPotential"),
                       "n_cm3": field("eDensity"), "p_cm3": field("hDensity")}
    result = {"native_poisson_linear_solves": first_dd,
              "native_first_coupled_debug_index": first_dd,
              "native_P_to_final_DD_potential_V": field_errors(
                  final_reference["psi_V"], poisson_reference),
              "independent_solver": json.loads((vela_dir / "summary.json").read_text()),
              "fields": {}, "poisson_same_state_operators": {}}
    for phase in ["poisson", "coupled"]:
        actual = np.genfromtxt(vela_dir / (phase + ".csv"), delimiter=",", names=True)
        if not np.array_equal(actual["node"], np.arange(n)):
            raise ValueError("Vela fields must contain each mesh node in order")
        reference = poisson_reference if phase == "poisson" else final_reference["psi_V"]
        result["fields"][phase] = {"potential_V": field_errors(actual["psi_V"], reference)}
        if phase == "coupled":
            for carrier in ["n_cm3", "p_cm3"]:
                result["fields"][phase][carrier] = density_errors(actual[carrier], final_reference[carrier])
    for path in sorted(vela_dir.glob("P*_matrix.csv")):
        idx = int(path.name.split("_")[0][1:])
        if idx >= first_dd:
            raise ValueError("Poisson operator compared against a coupled-stage iteration")
        data = np.genfromtxt(path, delimiter=",", names=True)
        matrix = sparse.coo_matrix((data["value"],
            (data["row"].astype(int), data["col"].astype(int))), shape=(n, n)).tocsr()
        rhs = np.genfromtxt(vela_dir / f"P{idx}_rhs.csv", delimiter=",", names=True)["rhs"]
        item = native_systems[idx]
        if len(item["rhs"]) != n:
            raise ValueError("pure Poisson native system must have one unknown per node")
        entries = item["matrix"]
        native = sparse.coo_matrix((list(entries.values()),
            ([r for r, _ in entries], [c for _, c in entries])), shape=(n, n)).tocsr()
        normalized = sparse.diags(1/matrix.diagonal()) @ matrix
        native_normalized = sparse.diags(1/native.diagonal()) @ native
        row_norm = lambda m: np.sqrt(np.asarray(m.multiply(m).sum(axis=1)).ravel())
        nb = np.asarray(item["rhs"]) / native.diagonal()
        difference = rhs / matrix.diagonal() - nb
        result["poisson_same_state_operators"][str(idx)] = {
            "dimension": n,
            "matrix_max_relative_row_L2": float(max(row_norm(normalized-native_normalized)
                                                   / row_norm(native_normalized))),
            "rhs_difference_diag_normalized_L2": float(np.linalg.norm(difference)),
            "rhs_native_diag_normalized_L2": float(np.linalg.norm(nb)),
            "rhs_relative_L2": float(np.linalg.norm(difference)/np.linalg.norm(nb)),
        }
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vela-dir", type=Path, required=True)
    ap.add_argument("--native-debug", type=Path, required=True)
    ap.add_argument("--native-super-debug", type=Path, required=True)
    ap.add_argument("--native-fields", type=Path, required=True)
    ap.add_argument("--native-potential-scale", type=float, default=.025851995266484913)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    result = compare(a.vela_dir, a.native_debug, a.native_super_debug, a.native_fields,
                     a.native_potential_scale)
    a.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["fields"], indent=2))
