#!/usr/bin/env python3
"""Read labeled sdevice -d matrices without assuming Vela's variable ordering.

RHS1 and SOL are retained with their native signs. The closure check tests both
signs, and its denominator is |J| |SOL| + |RHS1|, not an arbitrary unit floor.
Native matrix entries are indexed by (variable block, global mesh node).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import re
import statistics

VARIABLE = re.compile(r"VAR:Variable (.*):$")
VALUE = re.compile(r"(VAR|RHS1|RHS2|SOL):#(\d+)\s*[:=]\s*(\S+)$")
ENTRY = re.compile(r"JAC\((\d+),(\d+)\)([dul]):#(?:\((\d+),(\d+)\)|(\d+))\s*:\s*(\S+)$")


def iterations(lines):
    current = None
    block = None
    for raw in lines:
        line = raw.strip()
        if line.startswith("Newton loop "):
            if current is not None:
                yield current
            current = {"loop": int(line.split()[-1]), "names": [], "variables": {},
                       "rhs1": {}, "rhs2": {}, "solution": {}, "matrix": {}}
            block = None
        elif current is not None:
            match = VARIABLE.fullmatch(line)
            if match:
                block = len(current["names"])
                current["names"].append(match[1])
                continue
            match = re.fullmatch(r"(?:RHS[12]:)?Block (\d+)", line)
            if match:
                block = int(match[1])
                continue
            match = VALUE.fullmatch(line)
            if match:
                if block is None:
                    raise ValueError("value without a variable block")
                key = (block, int(match[2]))
                value = float(match[3])
                target = {"VAR": "variables", "RHS1": "rhs1", "RHS2": "rhs2", "SOL": "solution"}[match[1]]
                if key in current[target] or not math.isfinite(value):
                    raise ValueError("duplicate or nonfinite vector value")
                current[target][key] = value
                continue
            match = ENTRY.fullmatch(line)
            if match:
                r, c = int(match[1]), int(match[2])
                i, j = (int(match[6]), int(match[6])) if match[3] == "d" else (int(match[4]), int(match[5]))
                key = ((r, i), (c, j))
                value = float(match[7])
                if key in current["matrix"] or not math.isfinite(value):
                    raise ValueError("duplicate or nonfinite matrix entry")
                current["matrix"][key] = value
            elif line.startswith(("VAR:#", "RHS1:#", "RHS2:#", "SOL:#", "JAC(")):
                raise ValueError("unrecognized native matrix/vector syntax")
    if current is not None:
        yield current


def closure(item):
    keys = set(item["variables"])
    if any(set(item[k]) != keys for k in ("rhs1", "rhs2", "solution")):
        raise ValueError("incomplete variable/RHS/solution mapping")
    products = {k: [] for k in keys}
    for (row, col), value in item["matrix"].items():
        if row not in keys or col not in keys:
            raise ValueError("matrix entry references an absent variable")
        products[row].append(value * item["solution"][col])
    by_block = {}
    for block, name in enumerate(item["names"]):
        rows = [k for k in keys if k[0] == block]
        den = math.fsum(math.fsum(abs(v) for v in products[k]) + abs(item["rhs1"][k]) for k in rows)
        errors = [math.fsum(abs(math.fsum(products[k]) + sign * item["rhs1"][k]) for k in rows) / den if den else 0.0
                  for sign in (-1, 1)]
        by_block[name] = {"rows": len(rows), "J_SOL_minus_RHS1_relative": errors[0],
                          "J_SOL_plus_RHS1_relative": errors[1],
                          "rhs1_l1": math.fsum(abs(item["rhs1"][k]) for k in rows),
                          "rhs2_l1": math.fsum(abs(item["rhs2"][k]) for k in rows),
                          "max_abs_solution": max(abs(item["solution"][k]) for k in rows)}
    return by_block


def row_summary(matrix, rhs, solution, key):
    row = {col: value for (r, col), value in matrix.items() if r == key and value != 0.0}
    diagonal = row.get(key, 0.0)
    return {"nonzero_columns": len(row), "diagonal": diagonal,
            "absolute_row_sum": math.fsum(abs(v) for v in row.values()),
            "rhs": rhs.get(key), "solution": solution.get(key) if solution else None,
            "columns_by_block": {str(b): sorted(k[1] for k in row if k[0] == b)
                                 for b in sorted({k[0] for k in row})}}


def super_iterations(text):
    """Decode the native symmetric-pattern, asymmetric-value MSR storage.

    ja uses one-based pointers/columns; matrix[n] stores the edge count, followed
    by upper values and lower values. Validate lengths rather than guessing.
    """
    for chunk in text.split("new matrix:")[1:]:
        arrays = {}
        for match in re.finditer(r"(ja|matrix|rhs|result):\s*size = (\d+)\s*([\s\S]*?)(?=\n(?:ja|matrix|rhs|result):|\Z)", chunk):
            values = [float(v) for v in match[3].split()]
            if len(values) != int(match[2]) or not all(math.isfinite(v) for v in values):
                raise ValueError("invalid super.debug array length/value")
            if match[1] in arrays:
                raise ValueError("duplicate super.debug array")
            arrays[match[1]] = values
        if set(arrays) != {"ja", "matrix", "rhs", "result"}:
            raise ValueError("incomplete super.debug iteration")
        n = len(arrays["rhs"])
        if any(v != int(v) for v in arrays["ja"]):
            raise ValueError("noninteger sparse index")
        ja = [int(v) for v in arrays["ja"]]
        edges = len(ja) - n - 1
        a = arrays["matrix"]
        if edges < 0 or len(a) != n + 1 + 2 * edges or len(arrays["result"]) != n or a[n] != edges:
            raise ValueError("unsupported super.debug storage layout")
        if ja[0] != n + 2 or ja[n] != len(ja) + 1 or any(ja[i] > ja[i+1] for i in range(n)):
            raise ValueError("invalid MSR row pointers")
        matrix = {(i, i): a[i] for i in range(n)}
        for i in range(n):
            for k in range(ja[i]-1, ja[i+1]-1):
                j = ja[k] - 1
                if not i < j < n:
                    raise ValueError("invalid MSR upper column")
                if (i,j) in matrix:
                    raise ValueError("duplicate MSR column")
                matrix[i,j] = a[k]
                matrix[j,i] = a[k+edges]
        yield {"matrix": matrix, "rhs": arrays["rhs"], "solution": arrays["result"]}


def super_closure(item):
    n = len(item["rhs"])
    products = [[] for _ in range(n)]
    for (i,j), v in item["matrix"].items():
        products[i].append(v * item["solution"][j])
    denominator = math.fsum(math.fsum(abs(v) for v in row) + abs(item["rhs"][i]) for i,row in enumerate(products))
    error = math.fsum(abs(math.fsum(row) - item["rhs"][i]) for i,row in enumerate(products))
    # The historical field is an aggregate ratio, not the worst component.
    # A well-solved large row can hide a poorly solved weak row in that ratio.
    row_ratios = []
    for i, row in enumerate(products):
        scale = math.fsum([abs(item["rhs"][i]), *(abs(v) for v in row)])
        defect = abs(math.fsum([*row, -item["rhs"][i]]))
        row_ratios.append(defect / scale if scale else 0.0)
    worst = max(range(n), key=row_ratios.__getitem__) if n else None
    aggregate = error / denominator if denominator else 0.0
    return {"dimension": n, "nonzero_entries": sum(v != 0 for v in item["matrix"].values()),
            "linear_componentwise_l1_backward_error": aggregate,
            "linear_global_l1_backward_error": aggregate,
            "linear_componentwise_max_backward_error": row_ratios[worst] if n else 0.0,
            "linear_componentwise_worst_row": worst}


def read_vela(root):
    with (root / "rhs.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    mapping = {int(r["row"]): (int(r["row"]) // (len(rows) // 3), int(r["node_id"])) for r in rows}
    with (root / "jacobian.csv").open(newline="", encoding="utf-8") as f:
        matrix = {(mapping[int(r["row"])], mapping[int(r["col"])]): float(r["value"]) for r in csv.DictReader(f)}
    return rows, mapping, matrix


def scalar_fields(directory, name):
    result = {}
    for path in sorted(directory.glob(name + "_region*.csv")):
        with path.open(newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                key, value = int(r["node_id"]), float(r["component0"])
                if key in result and not math.isclose(result[key], value, abs_tol=1e-12, rel_tol=1e-12):
                    raise ValueError("inconsistent values at shared region node")
                result[key] = value
    if not result:
        raise ValueError("missing native scalar fields: " + name)
    return result


def density_update_audit(old, delta, accepted, scale):
    """Test an observed rule; do not claim access to Sentaurus implementation."""
    nonpositive = [k for k in old if old[k] + delta[k] <= 0]
    reconstructed = {k: scale * (old[k] + delta[k] if old[k] + delta[k] > 0 else 0.01 * old[k]) for k in old}
    errors = [abs(accepted[k] / reconstructed[k] - 1) for k in old]
    return {"nodes": len(old), "raw_update_nonpositive_count": len(nonpositive),
            "accepted_nonpositive_count": sum(accepted[k] <= 0 for k in old),
            "negative_update_accepted_over_old": [accepted[k] / (scale * old[k]) for k in nonpositive],
            "observed_rule": "if old+SOL > 0 use old+SOL; otherwise use 0.01*old",
            "observed_rule_max_relative_error": max(errors)}


def physical_comparison(root, selected):
    """Compare accepted fields and dimensionless local shapes, not raw units.

    J_native[:,p] p is the exact chain rule to log(p). Comparing it with Vela's
    phip block uses the nondegenerate minority-hole approximation and different
    accepted states, so it is a diagnostic, not a Jacobian correctness gate.
    """
    result = {}
    mesh = json.loads((root / "mesh.json").read_text())
    for label, item in selected.items():
        export = root / ("native_export_" + label)
        with (export / "nodes.csv").open(newline="", encoding="utf-8") as f:
            points = {int(r["id"]): (float(r["x_um"]), float(r["y_um"])) for r in csv.DictReader(f)}
        if set(points) != {r["id"] for r in mesh["nodes"]}:
            raise ValueError("native and Vela node IDs differ")
        geometry_error = max(max(abs(points[r["id"]][0]-r["x"]), abs(points[r["id"]][1]-r["y"])) for r in mesh["nodes"])
        if geometry_error > 1e-12:
            raise ValueError("native and Vela coordinates differ")
        fields = {name: scalar_fields(export / "fields", name) for name in (
            "ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential", "eDensity", "hDensity")}
        # Field CSVs retain the TDR units. Densities are cm^-3; Vela CSV is m^-3.
        manifest = json.loads((export / "field_manifest.json").read_text())
        if any(f["unit"] != "cm^-3" for f in manifest["fields"] if f["name"] in ("eDensity", "hDensity")):
            raise ValueError("unexpected native density units")
        rows, mapping, matrix = read_vela(root / ("vela_" + label))
        state = next(r for r in rows if mapping[int(r["row"])] == (0,790))
        native_shape = {col[1]: value * item["variables"][col] for (row,col), value in item["matrix"].items()
                        if row == (3,790) and col[0] == 3 and value != 0}
        native_diagonal = native_shape[790]
        native_shape = {k: v/native_diagonal for k,v in native_shape.items()}
        vela_shape = {col[1]: value for (row,col), value in matrix.items() if row == (2,790) and col[0] == 2 and value != 0}
        vela_diagonal = vela_shape[790]
        vela_shape = {k: v/vela_diagonal for k,v in vela_shape.items()}
        keys = native_shape.keys() | vela_shape.keys()
        difference = math.fsum(abs(native_shape.get(k,0)-vela_shape.get(k,0)) for k in keys)
        denominator = math.fsum(abs(v) for v in native_shape.values())
        normalization = {}
        updates = {}
        for b,name in ((2,"eDensity"),(3,"hDensity")):
            ratios = [fields[name][node] / (value + item["solution"][b,node]) for (block,node),value in item["variables"].items()
                      if block == b and value + item["solution"][b,node] > 0]
            median = statistics.median(ratios)
            normalization[name] = {"median_cm3_per_debug_unit": median,
                                   "max_relative_deviation": max(abs(v/median-1) for v in ratios), "samples": len(ratios)}
            old = {node: value for (block,node),value in item["variables"].items() if block == b}
            delta = {node: item["solution"][b,node] for node in old}
            updates[name] = density_update_audit(old, delta, fields[name], median)
        result[label] = {"geometry_max_error_um": geometry_error, "geometry_nodes": len(points),
                         "native_density_normalization": normalization,
                         "native_density_update_audit": updates,
                         "node790_native_accepted": {name: values[790] * (1e6 if name.endswith("Density") else 1) for name,values in fields.items()},
                         "node790_vela_accepted": {k: float(state[k]) for k in ("psi_V","phin_V","phip_V","electrons_m3","holes_m3")},
                         "node790_native_pre_update_logp_shape": native_shape,
                         "node790_vela_accepted_phip_shape": vela_shape,
                         "shape_relative_l1_difference": difference/denominator,
                         "shape_is_correctness_gate": False,
                         "limitation": "native pre-update versus Vela accepted states; density versus quasi-Fermi variables; minority-hole logp approximation; full Fermi cross-block transformation not performed"}
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--debug", type=Path, required=True)
    ap.add_argument("--vela-root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--super-debug", type=Path)
    ap.add_argument("--physical-comparison", action="store_true", help="BVmethods 2026-10-09 fixture only: compare native_export_{zero,first,second} and iterations 7/11/13")
    args = ap.parse_args()
    native = []
    selected = {}
    with args.debug.open(encoding="utf-8") as f:
        for index, item in enumerate(iterations(f)):
            if index in (7,11,13):
                selected[{7:"zero",11:"first",13:"second"}[index]] = item
            blocks = closure(item)
            native.append({"index": index, "loop": item["loop"], "blocks": blocks,
                           "dimension": len(item["variables"]),
                           "stored_entries": len(item["matrix"]),
                           "nonzero_entries": sum(v != 0 for v in item["matrix"].values()),
                           "node790": {name: {"variable": item["variables"].get((b,790)),
                                               **row_summary(item["matrix"], item["rhs1"], item["solution"], (b,790))}
                                       for b, name in enumerate(item["names"]) if (b,790) in item["variables"]}})
    vela = {}
    for label in ("zero", "first", "second"):
        rows, mapping, matrix = read_vela(args.vela_root / ("vela_" + label))
        rhs = {mapping[int(r["row"])]: float(r["rhs"]) for r in rows}
        vela[label] = {"dimension": len(rows), "stored_entries": len(matrix),
                       "nonzero_entries": sum(v != 0 for v in matrix.values()),
                       "node790": {name: {**row_summary(matrix, rhs, None, (b,790)),
                                             "physical_state": next(r for r in rows if mapping[int(r["row"])] == (b,790))}
                                   for b,name in enumerate(("psi", "phin", "phip"))}}
    result = {"native_iterations": native, "vela": vela}
    if args.super_debug:
        result["native_linear_systems"] = [super_closure(item) for item in super_iterations(args.super_debug.read_text())]
        if len(result["native_linear_systems"]) != len(native):
            raise ValueError("debug and super.debug iteration counts differ")
    if args.physical_comparison:
        if set(selected) != {"zero", "first", "second"} or [x["loop"] for x in native] != [0,1,2,3,4,5,6,0,0,0,1,2,0,0]:
            raise ValueError("BVmethods diagnostic iteration mapping incomplete")
        result["accepted_state_comparison"] = physical_comparison(args.vela_root, selected)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(json.dumps({"native_iterations": len(native), "dimensions": sorted({r["dimension"] for r in native})}))


if __name__ == "__main__":
    main()
