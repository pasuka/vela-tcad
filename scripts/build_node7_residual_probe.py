#!/usr/bin/env python3
"""Build a newton_carrier_term_probe deck from the baseline sandbox -20V state.

Extracts the converged -20V state (Potential/ElectronQuasiFermi/HoleQuasiFermi)
from the baseline VTK into the node_id,component0 CSVs that readExternalState
wants, copies mesh/doping/materials, and writes a newton_carrier_term_probe
config reusing the baseline solver physics. Running the produced deck dumps the
COMPLETE per-node continuity residual + term decomposition (all edges incl
couple<=0, plus SRH and the avalanche source), giving the exact node7 balance.
Std lib only.
"""

from __future__ import annotations

import json
import os
import shutil

BASE = "build/qf_cap_exp/baseline"
OUT = "build/qf_cap_exp/residual_probe"
VTK = os.path.join(BASE, "vtk", "bv_0400_-20V.vtk")

FIELD_MAP = {
    "Potential": "ElectrostaticPotential",
    "ElectronQuasiFermi": "eQuasiFermiPotential",
    "HoleQuasiFermi": "hQuasiFermiPotential",
}


def parse_vtk_scalars(path, wanted):
    with open(path) as h:
        lines = h.read().splitlines()
    # find POINTS count
    npts = None
    for ln in lines:
        if ln.startswith("POINTS"):
            npts = int(ln.split()[1])
            break
    out = {}
    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if ln.startswith("SCALARS"):
            name = ln.split()[1]
            if name in wanted:
                # next line is LOOKUP_TABLE, then npts values
                vals = []
                j = i + 2
                while len(vals) < npts and j < len(lines):
                    for tok in lines[j].split():
                        vals.append(float(tok))
                    j += 1
                out[name] = vals[:npts]
        i += 1
    return out, npts


def main():
    os.makedirs(os.path.join(OUT, "state_fields"), exist_ok=True)
    scalars, npts = parse_vtk_scalars(VTK, set(FIELD_MAP))
    for vtk_name, csv_base in FIELD_MAP.items():
        vals = scalars[vtk_name]
        p = os.path.join(OUT, "state_fields", csv_base + "_region0.csv")
        with open(p, "w") as h:
            h.write("node_id,component0\n")
            for nid, v in enumerate(vals):
                h.write(f"{nid},{v!r}\n")
    # copy mesh/doping/materials
    for f in ("mesh.json", "doping.csv", "pn2d_sentaurus2018_iv_materials.json"):
        shutil.copy(os.path.join(BASE, f), os.path.join(OUT, f))

    # reuse baseline solver physics
    with open(os.path.join(BASE, "simulation_bv.json")) as h:
        base = json.load(h)
    solver = dict(base["solver"])

    cfg = {
        "_comment": "newton_carrier_term_probe: exact per-node continuity residual at -20V",
        "simulation_type": "newton_carrier_term_probe",
        "mesh_file": "mesh.json",
        "node_doping_file": "doping.csv",
        "materials_file": "pn2d_sentaurus2018_iv_materials.json",
        "scaling": {"mode": "unit_scaling"},
        "contacts": [
            {"name": "Anode", "bias": -20.0},
            {"name": "Cathode", "bias": 0.0},
        ],
        "solver": solver,
        "state_fields_dir": "state_fields",
        "output_csv": "carrier_terms.csv",
    }
    with open(os.path.join(OUT, "simulation_probe.json"), "w") as h:
        json.dump(cfg, h, indent=2)
    print(f"npts={npts}; wrote {OUT}/simulation_probe.json + state_fields + mesh/doping/materials")


if __name__ == "__main__":
    main()
