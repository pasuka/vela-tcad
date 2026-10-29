"""WP4a: gate-contact reference-level check (Sentaurus manual Eq. 104 vs Vela ``metal_gate``).

Sentaurus (contacts on insulators):  psi_gate = phi_F - Phi_MS          (Eq. 104)
Vela ``metal_gate``:                 psi_gate = bias  - work_function_eV (BoundaryCondition.cpp)

Both engines reference psi = 0 to the intrinsic semiconductor (ohmic contacts use
psi = V + Vt*asinh(N/2ni) in both), so the exact mapping is

    work_function_eV = Phi_MS            for ``Barrier = Phi_MS``
    work_function_eV = Phi_M - Phi_S,ref for ``WorkFunction = Phi_M`` with
        Phi_S,ref = chi + Eg(0)/2 + (kT/2) ln(Nc(T)/Nv(T))
    Native WP14 T-2022.03 controls establish the Eg(0) contact convention.
    The physical intrinsic work function uses Eg(T), and remains separate.

and the flat-band voltage of an oxide-on-uniform-semiconductor capacitor with an
ohmic back contact at 0 V is

    V_FB = Phi_MS + psi_bulk,   psi_bulk = Vt*asinh(N_net / (2 ni_eff))  (sign of N_net).

This module builds a 1-D-equivalent MOS capacitor deck for ``vela_example_runner``,
sweeps the gate with mandatory bias points around the prediction, extracts the
flat-band voltage from the surface/bulk potential difference (root of a
monotone function, secant-refined) and reports the deviation from Eq. 104.

Usage::

    python scripts/check_sdevice_gate_contact_mapping.py --runner build/vela_example_runner.exe \
        --out build/gate_contact_check [--barriers 0.0 0.55 -0.3]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Sequence

REPO = Path(__file__).resolve().parents[1]

KB = 1.380649e-23
Q = 1.602176634e-19
T = 300.0
VT = KB / Q * T

#: C++ built-ins are checked against this emitted no-BGN 300 K snapshot.
_SI = next(entry for entry in json.loads(
    (REPO / "share/vela/materials/default_materials.json").read_text(encoding="utf-8")
)["materials"] if entry["name"] == "Si")
SI_NI_CM3 = _SI["intrinsic_carrier_density_cm3"]
SI_CHI_EV = _SI["electron_affinity_eV"]
SI_EG_EV = _SI["bandgap_eV"]
SI_NC_CM3 = _SI["conduction_band_density_of_states_cm3"]
SI_NV_CM3 = _SI["valence_band_density_of_states_cm3"]

RUNNER_CANDIDATES = [REPO / "build" / "vela_example_runner.exe", REPO / "build" / "vela_example_runner",
                     REPO / "build-release" / "vela_example_runner.exe",
                     REPO / "build-release" / "vela_example_runner"]


def find_runner(explicit: str | None = None) -> Path | None:
    candidates = [Path(explicit)] if explicit else []
    env = os.environ.get("VELA_EXAMPLE_RUNNER")
    if env:
        candidates.append(Path(env))
    candidates.extend(RUNNER_CANDIDATES)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def intrinsic_work_function_eV(chi_eV: float = SI_CHI_EV, eg_eV: float = SI_EG_EV,
                               nc_cm3: float = SI_NC_CM3, nv_cm3: float = SI_NV_CM3) -> float:
    """Phi_S,i = chi + Eg/2 + (kT/2) ln(Nc/Nv): metal work function of an intrinsic semiconductor."""
    return chi_eV + 0.5 * eg_eV + 0.5 * VT * math.log(nc_cm3 / nv_cm3)


def phi_ms_from_workfunction(workfunction_eV: float, **material) -> float:
    return workfunction_eV - sentaurus_gate_reference_eV(**material)


def sentaurus_gate_reference_eV(chi_eV: float = SI_CHI_EV, eg_zero_eV: float = 1.16964,
                              nc_cm3: float = SI_NC_CM3, nv_cm3: float = SI_NV_CM3) -> float:
    """T-2022.03 gate WorkFunction reference: Eg(0) and DOS at 300 K.

    Native WP14 Barrier/WorkFunction controls establish this normalization;
    it does not replace the 300 K bandgap in the material or ni calculation.
    """
    return chi_eV + 0.5 * eg_zero_eV + 0.5 * VT * math.log(nc_cm3 / nv_cm3)


def bulk_potential_V(donors_cm3: float, acceptors_cm3: float, ni_cm3: float = SI_NI_CM3) -> float:
    """Ohmic-contact built-in potential Vt*asinh(N/2ni) (Eq. 101 / Vela ohmic BC)."""
    return VT * math.asinh((donors_cm3 - acceptors_cm3) / (2.0 * ni_cm3))


def predicted_flatband_V(phi_ms_V: float, donors_cm3: float, acceptors_cm3: float) -> float:
    return phi_ms_V + bulk_potential_V(donors_cm3, acceptors_cm3)


# ---------------------------------------------------------------------------
# Fixture: oxide strip on a uniformly doped silicon column (coordinates in um)
# ---------------------------------------------------------------------------

def build_mesh(si_thickness_um: float = 1.0, ox_thickness_um: float = 0.01, width_um: float = 0.05,
               si_layers: int = 24, ox_layers: int = 2) -> dict:
    """Two-column strip: x in [0, width], y from 0 (substrate) to si+ox (gate)."""
    ys: List[float] = []
    # Graded silicon layers (finer near the interface) then uniform oxide layers.
    for i in range(si_layers + 1):
        s = i / si_layers
        ys.append(si_thickness_um * (1.0 - (1.0 - s) ** 2))  # dense near y = si_thickness
    ys = sorted(set(round(y, 12) for y in ys))
    for j in range(1, ox_layers + 1):
        ys.append(round(si_thickness_um + ox_thickness_um * j / ox_layers, 12))
    nodes = []
    node_id = 0
    index = {}
    for row, y in enumerate(ys):
        for col, x in enumerate((0.0, width_um)):
            nodes.append({"id": node_id, "x": x, "y": y})
            index[(row, col)] = node_id
            node_id += 1
    triangles = []
    regions = {0: [], 1: []}
    tri_id = 0
    for row in range(len(ys) - 1):
        y_mid = 0.5 * (ys[row] + ys[row + 1])
        region_id = 0 if y_mid < si_thickness_um else 1
        a, b = index[(row, 0)], index[(row, 1)]
        c, d = index[(row + 1, 1)], index[(row + 1, 0)]
        for tri in ((a, b, c), (a, c, d)):
            triangles.append({"id": tri_id, "region_id": region_id, "node_ids": list(tri)})
            regions[region_id].append(tri_id)
            tri_id += 1
    top = len(ys) - 1
    interface_row = ys.index(round(si_thickness_um, 12))
    mesh = {
        "nodes": nodes,
        "triangles": triangles,
        "regions": [
            {"id": 0, "name": "R.Si", "material": "Si", "cell_ids": regions[0]},
            {"id": 1, "name": "R.Ox", "material": "SiO2", "cell_ids": regions[1]},
        ],
        "contacts": [
            {"id": 0, "name": "substrate", "region_id": 0, "node_ids": [index[(0, 0)], index[(0, 1)]]},
            {"id": 1, "name": "gate", "region_id": 1, "node_ids": [index[(top, 0)], index[(top, 1)]]},
        ],
    }
    probes = {
        "surface": [index[(interface_row, 0)], index[(interface_row, 1)]],
        "bulk": [index[(1, 0)], index[(1, 1)]],
        "gate": [index[(top, 0)], index[(top, 1)]],
        "substrate": [index[(0, 0)], index[(0, 1)]],
    }
    return {"mesh": mesh, "probes": probes}


def build_deck(mesh_file: str, phi_ms_V: float, bias_points: Sequence[float], acceptors_cm3: float,
               donors_cm3: float, csv_file: str, vtk_prefix: str) -> dict:
    return {
        "simulation_type": "dc_sweep",
        "mesh_file": mesh_file,
        "scaling": {"mode": "unit_scaling"},
        "doping": [{"region": "R.Si", "donors": donors_cm3, "acceptors": acceptors_cm3}],
        "contacts": [
            {"name": "substrate", "type": "ohmic", "bias": 0.0},
            {"name": "gate", "type": "metal_gate", "bias": 0.0, "work_function_eV": phi_ms_V},
        ],
        "solver": {"method": "newton", "max_iter": 200, "reltol": 1.0e-10, "abstol": 1.0e-14,
                   "recombination": [], "verbose": False},
        "sweep": {
            "mode": "iv", "contact": "gate", "current_contact": "substrate",
            "start": float(bias_points[0]), "stop": float(bias_points[-1]), "step": 0.05,
            "bias_points": [float(v) for v in bias_points],
            "csv_file": csv_file, "write_vtk": True, "vtk_prefix": vtk_prefix,
            "vtk_format": "vtu", "vtk_fields": ["ElectrostaticPotential"],
        },
        "runtime_log": {"enabled": False},
    }


_ARRAY_RE = re.compile(r'<DataArray[^>]*Name="ElectrostaticPotential"[^>]*>([\s\S]*?)</DataArray>')
_BIAS_RE = re.compile(r'<DataArray[^>]*Name="Bias_V"[^>]*>\s*([-+0-9.eE]+)')


def read_vtu_potential(path: Path) -> tuple[float, List[float]]:
    text = path.read_text(encoding="utf-8")
    bias = _BIAS_RE.search(text)
    array = _ARRAY_RE.search(text)
    if not bias or not array:
        raise RuntimeError(f"{path}: ElectrostaticPotential/Bias_V not found")
    return float(bias.group(1)), [float(v) for v in array.group(1).split()]


def run_points(runner: Path, workdir: Path, phi_ms_V: float, bias_points: Sequence[float],
               probes: dict, acceptors_cm3: float, donors_cm3: float, tag: str) -> Dict[float, float]:
    """Return {gate bias: psi_surface - psi_bulk} for the requested points."""
    deck = build_deck("mesh.json", phi_ms_V, bias_points, acceptors_cm3, donors_cm3,
                      f"{tag}.csv", f"{tag}_vtk/{tag}")
    deck_path = workdir / f"{tag}.json"
    deck_path.write_text(json.dumps(deck, indent=2), encoding="utf-8")
    env = {**os.environ, "Path": r"D:\msys64\ucrt64\bin;D:\msys64\usr\bin;" + os.environ.get("Path", "")}
    proc = subprocess.run([str(runner), "--config", deck_path.name], cwd=str(workdir), capture_output=True,
                          text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f"runner failed ({proc.returncode}) for {tag}:\n{proc.stdout[-3000:]}\n{proc.stderr[-3000:]}")
    result: Dict[float, float] = {}
    for vtu in sorted((workdir / f"{tag}_vtk").glob("*.vtu")):
        bias, psi = read_vtu_potential(vtu)
        surface = sum(psi[i] for i in probes["surface"]) / len(probes["surface"])
        bulk = sum(psi[i] for i in probes["bulk"]) / len(probes["bulk"])
        result[bias] = surface - bulk
    if len(result) != len(bias_points):
        raise RuntimeError(f"{tag}: expected {len(bias_points)} VTU files, found {len(result)}")
    return result


def flatband_from_samples(samples: Dict[float, float]) -> float:
    """Root of the monotone psi_surface - psi_bulk (V_g) by linear interpolation of the bracketing pair."""
    points = sorted(samples.items())
    for (v0, f0), (v1, f1) in zip(points, points[1:]):
        if f0 == 0.0:
            return v0
        if f0 * f1 < 0.0:
            return v0 - f0 * (v1 - v0) / (f1 - f0)
    raise RuntimeError(f"flat band not bracketed by samples {points}")


def check_barrier(runner: Path, workdir: Path, phi_ms_V: float, probes: dict, acceptors_cm3: float,
                  donors_cm3: float, half_window_V: float = 0.02, refine: int = 2) -> dict:
    predicted = predicted_flatband_V(phi_ms_V, donors_cm3, acceptors_cm3)
    tag = f"barrier_{phi_ms_V:+.3f}".replace("+", "p").replace("-", "m").replace(".", "_")
    window = [predicted - half_window_V, predicted - 0.5 * half_window_V, predicted,
              predicted + 0.5 * half_window_V, predicted + half_window_V]
    samples = run_points(runner, workdir, phi_ms_V, window, probes, acceptors_cm3, donors_cm3, tag)
    estimate = flatband_from_samples(samples)
    history = [estimate]
    for step in range(refine):
        # Tighten the bracket around the current estimate; a 1 mV window gives a
        # linear-interpolation error far below the 1 mV acceptance.
        width = half_window_V / (10.0 ** (step + 1))
        window = [estimate - width, estimate, estimate + width]
        samples = run_points(runner, workdir, phi_ms_V, window, probes, acceptors_cm3, donors_cm3,
                             f"{tag}_r{step}")
        estimate = flatband_from_samples(samples)
        history.append(estimate)
    return {
        "phi_ms_V": phi_ms_V,
        "predicted_flatband_V": predicted,
        "measured_flatband_V": estimate,
        "deviation_V": estimate - predicted,
        "bulk_potential_V": bulk_potential_V(donors_cm3, acceptors_cm3),
        "estimates": history,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runner", help="vela_example_runner executable")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--barriers", type=float, nargs="+", default=[0.0, 0.55, -0.3],
                        help="Phi_MS values (V) to check")
    parser.add_argument("--acceptors-cm3", type=float, default=1.0e17)
    parser.add_argument("--donors-cm3", type=float, default=0.0)
    parser.add_argument("--workfunction-eV", type=float, nargs="*", default=[4.05, 5.1],
                        help="metal work functions to convert through Phi_S,i (reported, not simulated)")
    args = parser.parse_args(argv)

    runner = find_runner(args.runner)
    if runner is None:
        print("vela_example_runner not found", file=sys.stderr)
        return 2
    args.out.mkdir(parents=True, exist_ok=True)
    fixture = build_mesh()
    (args.out / "mesh.json").write_text(json.dumps(fixture["mesh"]), encoding="utf-8")
    report = {
        "schema": "vela.sdevice_gate_contact_check.v1",
        "runner": str(runner),
        "temperature_K": T,
        "kT_eV": VT,
        "material": {"ni_cm3": SI_NI_CM3, "chi_eV": SI_CHI_EV, "eg_eV": SI_EG_EV, "nc_cm3": SI_NC_CM3,
                     "nv_cm3": SI_NV_CM3},
        "doping": {"acceptors_cm3": args.acceptors_cm3, "donors_cm3": args.donors_cm3},
        "intrinsic_work_function_eV": intrinsic_work_function_eV(),
        "sentaurus_gate_reference_eV": sentaurus_gate_reference_eV(),
        "mapping": {
            "Barrier": "work_function_eV = Barrier (Phi_MS)",
            "WorkFunction": "work_function_eV = WorkFunction - (chi + Eg(0)/2 + kT/2 ln(Nc(T)/Nv(T)))",
            "flatband": "V_FB = Phi_MS + Vt asinh(N_net/(2 ni))",
        },
        "workfunction_examples": [
            {"workfunction_eV": wf, "phi_ms_V": phi_ms_from_workfunction(wf)} for wf in args.workfunction_eV
        ],
        "barriers": [check_barrier(runner, args.out, phi, fixture["probes"], args.acceptors_cm3, args.donors_cm3)
                     for phi in args.barriers],
    }
    report["max_abs_deviation_V"] = max(abs(b["deviation_V"]) for b in report["barriers"])
    (args.out / "gate_contact_check.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    for entry in report["barriers"]:
        print(f"Phi_MS={entry['phi_ms_V']:+.3f} V  predicted V_FB={entry['predicted_flatband_V']:+.6f}  "
              f"measured={entry['measured_flatband_V']:+.6f}  deviation={entry['deviation_V'] * 1e3:+.4f} mV")
    print(f"max |deviation| = {report['max_abs_deviation_V'] * 1e3:.4f} mV")
    return 0 if report["max_abs_deviation_V"] <= 1.0e-3 else 1


if __name__ == "__main__":
    raise SystemExit(main())
