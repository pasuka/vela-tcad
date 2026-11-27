#!/usr/bin/env python3
"""Prepare a strict frozen-state replay of the Sentaurus SLOT-LDMOS 1 V state."""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


SILICON_REGION = "Silicon_1"
OXIDE_REGION = "Oxide_1"


def _absolute(source_bundle: Path, value: str) -> str:
    return (source_bundle / value).resolve().as_posix()


def build_config(
    source: dict,
    source_bundle: Path,
    state: Path,
    output_directory: Path | None = None,
    bias_V: float = 1.0,
    mobility_model: str = "masetti_field_lombardi",
) -> dict:
    deck = copy.deepcopy(source)
    deck["_comment"] = (
        f"Diagnostic-only frozen replay of the Sentaurus {bias_V:g} V state. "
        "No nonlinear solve and no avalanche source are evaluated."
    )
    deck["mesh_file"] = _absolute(source_bundle, source["mesh_file"])
    deck["node_doping_file"] = _absolute(
        source_bundle, source["node_doping_file"])
    deck["materials_file"] = _absolute(
        source_bundle, source["materials_file"])
    output_root = (
        (output_directory / "outputs").resolve()
        if output_directory is not None
        else Path("outputs")
    )
    deck["output_csv"] = (output_root / "iv.csv").as_posix()

    solver = deck["solver"]
    solver["method"] = "frozen_state"
    solver["max_iter"] = 0
    solver.pop("handoff", None)
    solver.pop("impact_ionization", None)
    solver["recombination"] = ["srh", "auger"]
    solver["taun"] = 1.0e-5
    solver["taup"] = 3.0e-6
    solver["srh_doping_dependence"] = {
        "enabled": True,
        "concentration_basis": "total_impurity",
        "electron": {
            "tau_min_s": 0.0,
            "tau_max_s": 1.0e-5,
            # unit_scaling decks use the internal TCAD concentration unit cm^-3.
            "reference_doping_m3": 1.0e16,
            "gamma": 1.0,
        },
        "hole": {
            "tau_min_s": 0.0,
            "tau_max_s": 3.0e-6,
            "reference_doping_m3": 1.0e16,
            "gamma": 1.0,
        },
    }
    # unit_scaling decks use cm^6/s for Auger coefficient overrides.
    solver["auger_cn_m6_per_s"] = 2.90e-31
    solver["auger_cp_m6_per_s"] = 1.028e-31
    solver["mobility"] = {
        "model": mobility_model,
        "doping_concentration_basis": "total_impurity",
        "high_field_driving_force": "quasi_fermi_gradient",
        "high_field_gradient_discretization": "transport_cell_vector",
    }
    if mobility_model == "masetti_field_lombardi":
        solver["mobility"]["surface"] = {
            "surface_interface": [SILICON_REGION, OXIDE_REGION],
        }

    sweep = deck["sweep"]
    sweep["bias_points"] = [bias_V]
    sweep["start"] = bias_V
    sweep["stop"] = bias_V
    sweep["step"] = 1.0
    sweep["initial_state_file"] = state.resolve().as_posix()
    sweep["frozen_state_compute_current"] = True
    sweep["write_vtk"] = False
    sweep.pop("external_circuit", None)
    sweep.pop("boundary_control", None)
    sweep.pop("breakdown", None)
    sweep.pop("write_state_file", None)
    sweep.pop("write_state_every_point_prefix", None)
    sweep["diagnostics"] = {
        "terminal_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "substrate"],
            "csv_file": (output_root / "terminal_balance.csv").as_posix(),
        },
        "srh_balance": {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": (output_root / "recombination_balance.csv").as_posix(),
        },
        "contact_edge": {
            "enabled": True,
            "contacts": ["source", "drain", "substrate"],
            "csv_file": (output_root / "contact_edges.csv").as_posix(),
        },
        "continuity_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "substrate"],
            "csv_file": (output_root / "continuity_balance.csv").as_posix(),
        },
    }
    return deck


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-bundle", type=Path, required=True)
    parser.add_argument(
        "--source-config",
        default="simulation_04_avalanche_activation_1v.json",
    )
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=1.0)
    parser.add_argument(
        "--mobility-model",
        choices=("masetti_field_lombardi", "masetti_field", "masetti"),
        default="masetti_field_lombardi",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    source_bundle = args.source_bundle.resolve()
    source = json.loads(
        (source_bundle / args.source_config).read_text(encoding="utf-8"))
    audit_directory = args.output.parent / f"{args.output.name}_outputs"
    deck = build_config(
        source, source_bundle, args.state, audit_directory, args.bias,
        args.mobility_model)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    (audit_directory / "outputs").mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(deck, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
