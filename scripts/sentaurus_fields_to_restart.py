#!/usr/bin/env python3
"""Merge a multi-region Sentaurus TDR CSV export into a Vela restart state."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


FIELDS = {
    "psi": "ElectrostaticPotential",
    "phin": "eQuasiFermiPotential",
    "phip": "hQuasiFermiPotential",
    "electrons_m3": "eDensity",
    "holes_m3": "hDensity",
    "electron_quantum_potential_V": "eQuantumPotential",
}
CONDUCTION_BAND_FIELD = "ConductionBandEnergy"
AFFINITY_FIELD = "ElectronAffinity"
ELECTROSTATIC_FIELD = "ElectrostaticPotential"
THERMAL_VOLTAGE_300K_V = 0.025851999786
SILICON_ELECTRON_DOS_MASS_RATIO = 1.0618016171622988


def read_field(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node = int(row["node_id"])
            value = float(row["component0"])
            if node in values and not math.isclose(values[node], value, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"inconsistent duplicate node {node} in {path}")
            values[node] = value
    return values


def merge_field(export_dir: Path, manifest: dict, name: str) -> tuple[dict[int, float], str]:
    merged: dict[int, float] = {}
    unit = ""
    for entry in manifest.get("fields", []):
        if entry.get("name") != name or entry.get("mapping_status") != "complete":
            continue
        unit = str(entry.get("unit", unit))
        values = read_field(export_dir / "fields" / str(entry["csv_file"]))
        for node, value in values.items():
            if node in merged and not math.isclose(merged[node], value, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"field {name} differs across regions at shared node {node}")
            merged[node] = value
    return merged, unit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--preserve-insulator-quantum-potential",
        action="store_true",
        help=("keep the global Sentaurus eQuantumPotential field on nodes "
              "without carrier support (required by include_insulators)"),
    )
    parser.add_argument(
        "--allow-missing-electron-quantum-potential",
        action="store_true",
        help=("treat an entirely absent eQuantumPotential field as zero; "
              "use only for classical DD Sentaurus states"),
    )
    args = parser.parse_args()

    manifest = json.loads((args.export_dir / "field_manifest.json").read_text())
    node_count = sum(1 for _ in (args.export_dir / "nodes.csv").open()) - 1
    columns: dict[str, dict[int, float]] = {}
    units: dict[str, str] = {}
    for output_name, sentaurus_name in FIELDS.items():
        columns[output_name], units[output_name] = merge_field(
            args.export_dir, manifest, sentaurus_name)
    conduction_band, conduction_band_unit = merge_field(
        args.export_dir, manifest, CONDUCTION_BAND_FIELD)
    affinity, affinity_unit = merge_field(
        args.export_dir, manifest, AFFINITY_FIELD)

    quantum = columns["electron_quantum_potential_V"]
    if not quantum and args.allow_missing_electron_quantum_potential:
        quantum.update({node: 0.0 for node in range(node_count)})
        units["electron_quantum_potential_V"] = "V"

    required_all_nodes = ("psi", "phin", "phip", "electron_quantum_potential_V")
    for name in required_all_nodes:
        missing = set(range(node_count)) - columns[name].keys()
        if missing:
            raise ValueError(f"{name} is missing {len(missing)} nodes")
    missing_conduction_band = set(range(node_count)) - conduction_band.keys()
    if missing_conduction_band:
        raise ValueError(
            "ConductionBandEnergy is missing "
            f"{len(missing_conduction_band)} nodes")
    if conduction_band_unit not in {"", "eV", "V"}:
        raise ValueError(
            "unsupported ConductionBandEnergy unit: " + conduction_band_unit)
    if affinity_unit not in {"", "eV", "V"}:
        raise ValueError("unsupported ElectronAffinity unit: " + affinity_unit)

    # Sentaurus band energies carry an arbitrary global energy origin. Recover
    # and remove it before constructing Phi/q; otherwise the constant appears
    # as a several-volt quantum correction in Vela's reaction term.
    band_origins = [
        conduction_band[node] + columns["psi"][node] + affinity[node]
        for node in range(node_count)
    ]
    band_origin = statistics.median(band_origins)
    if max(abs(value - band_origin) for value in band_origins) > 1.0e-8:
        raise ValueError(
            "ConductionBandEnergy + ElectrostaticPotential + "
            "ElectronAffinity is not a single global energy origin")

    transport_nodes: set[int] = set()
    with (args.export_dir / "elements.csv").open(newline="") as stream:
        for element in csv.DictReader(stream):
            if element["material"] not in {"Silicon", "Si", "PolySilicon"}:
                continue
            transport_nodes.update(
                int(element[f"node{local}"]) for local in range(3))


    # Densities exist only on carrier-supporting regions. Insulator values are
    # unused by the DD assembler and are written as zero.
    density_scale = {
        name: 1.0e6 if units[name].replace(" ", "") in {"cm^-3", "1/cm^3"} else 1.0
        for name in ("electrons_m3", "holes_m3")
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "node_id", *FIELDS, "electron_quantum_potential_like_V"])
        for node in range(node_count):
            # Sentaurus may carry arbitrary quasi-Fermi placeholders through
            # insulators.  Vela represents those non-transport rows as pinned
            # algebraic unknowns, so normalize nodes that have no density
            # support instead of feeding the placeholders to Newton.
            has_carrier_support = (
                node in columns["electrons_m3"] or
                node in columns["holes_m3"]
            )
            writer.writerow([
                node,
                columns["psi"][node],
                columns["phin"][node] if has_carrier_support else 0.0,
                columns["phip"][node] if has_carrier_support else 0.0,
                columns["electrons_m3"].get(node, 0.0) * density_scale["electrons_m3"],
                columns["holes_m3"].get(node, 0.0) * density_scale["holes_m3"],
                (columns["electron_quantum_potential_V"][node]
                 if (has_carrier_support or
                     args.preserve_insulator_quantum_potential) else 0.0),
                # Phi/q = Ec/q + Phi_m/q + Lambda. Vela uses the transport
                # trace at shared nodes.  Nc scales as m_DOS^(3/2), hence
                # Phi_m/q = -1.5*Vt*ln(m_DOS/m0); the stored primary subtracts
                # the corresponding +1.5*Vt*ln(m_DOS/m0) band-drive term.
                # Ec+Lambda has a conforming node trace in the TDR. Vela's
                # cell-side material drive supplies affinity, BGN, and the DOS
                # mass term; only the arbitrary energy origin is removed here.
                conduction_band[node] +
                columns["electron_quantum_potential_V"][node] -
                band_origin -
                1.5 * THERMAL_VOLTAGE_300K_V * math.log(
                    (SILICON_ELECTRON_DOS_MASS_RATIO
                     if node in transport_nodes else 0.42)),
            ])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
