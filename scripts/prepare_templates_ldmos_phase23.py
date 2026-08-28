#!/usr/bin/env python3
"""Prepare exact-mesh Templates/LDMOS phase-2 and G3 phase-3 decks.

The generator consumes the sealed stage-1 topology and the stage-1.5
Sentaurus restart.  It derives the PolySi(N) electrostatic offset from the
actual constant gate-boundary potential instead of embedding a handbook
value.  Generated run products stay below ignored ``reference_staging``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Sequence


CONTACTS = ("source", "drain", "gate", "substrate")

# Accepted drain voltages from the sealed T-2022.03-SP2 G3 quasistationary
# trajectory (InitialStep=0.01, Increment=1.35, MaxStep=0.2, Goal=0.1 V).
# Preserve this path as an exact-mesh WP1.5 diagnostic: a direct 0->0.1 V jump
# is not the Sentaurus oracle's numerical contract.  It remains separate from
# the production prebias until Vela can traverse and reclose the full path.
SENTAURUS_G3_DRAIN_PREBIAS_POINTS_V = [
    0.0,
    0.001,
    0.00214466666666667,
    0.00366859955555556,
    0.00569746220829630,
    0.00839855468664514,
    0.0119946091394869,
    0.0167821563010369,
    0.0231559774221138,
    0.0316416579413075,
    0.0429389272725274,
    0.0579793585088248,
    0.0777224312450045,
    0.0977224312450045,
    0.1,
]

# The first accepted Sentaurus point at which the G3 transition discrepancy is
# audited.  Keep this as a prefix of the sealed trajectory so the checkpoint is
# reproducible without relying on a later-state restart or interpolation.
SENTAURUS_G3_VD0P023_CHECKPOINT_POINTS_V = (
    SENTAURUS_G3_DRAIN_PREBIAS_POINTS_V[:9]
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def read_scalar_field(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(stream)
        }


def derive_polysi_flatband(mesh_path: Path, sentaurus_export: Path) -> dict[str, Any]:
    mesh = read_json(mesh_path)
    gate = next(contact for contact in mesh["contacts"] if contact["name"] == "gate")
    field_path = (
        sentaurus_export / "fields" /
        f"ElectrostaticPotential_region{int(gate['region_id'])}.csv"
    )
    values = read_scalar_field(field_path)
    boundary = [values[int(node)] for node in gate["node_ids"]]
    if not boundary:
        raise ValueError("gate contact has no values in the Sentaurus potential field")
    spread = max(boundary) - min(boundary)
    if spread > 1.0e-12:
        raise ValueError(
            f"Sentaurus gate potential is not constant (spread={spread:.17g} V)")
    gate_potential = math.fsum(boundary) / len(boundary)
    return {
        "sentaurus_gate_potential_V": gate_potential,
        "vela_flatband_voltage_V": -gate_potential,
        "node_count": len(boundary),
        "spread_V": spread,
        "mapping": "psi_gate = bias - flatband_voltage",
        "source": str(field_path.resolve()),
    }


def classical_solver(physics: dict[str, Any], *, high_field: bool) -> dict[str, Any]:
    mobility = physics["mobility"]
    electron = mobility["electron"]
    hole = mobility["hole"]
    recombination = physics["recombination"]
    srh = recombination["srh_doping_dependence"]
    bgn = physics["bandgap_narrowing"]
    mobility_config: dict[str, Any] = {
        "model": mobility["g3_model"] if high_field else "constant",
    }
    if high_field:
        if not mobility["high_field_saturation_enabled"]:
            raise ValueError("G3 requires high_field_saturation_enabled=true")
        if mobility["doping_dependence_enabled"]:
            raise ValueError("G3-no-IALMob must not enable DopingDependence")
        mobility_config.update({
            "high_field_driving_force": mobility["high_field_driving_force"],
            "high_field_gradient_discretization": mobility["high_field_gradient_discretization"],
            # Legacy unit_scaling keys consume TCAD internal cm-based values.
            # Do not convert these contract values to SI despite the historical
            # key suffixes; format-version-2 migration will rename them later.
            "electron_saturation_velocity_m_s": electron["saturation_velocity_cm_per_s"],
            "electron_high_field_beta": electron["high_field_exponent"],
            "hole_saturation_velocity_m_s": hole["saturation_velocity_cm_per_s"],
            "hole_high_field_beta": hole["high_field_exponent"],
        })
    return {
        "method": "newton",
        "max_iter": 100,
        "reltol": 1.0e-7,
        "abstol": 1.0e-12,
        "damping_factor": 1.0,
        "warm_start": True,
        "carrier_statistics": {"model": physics["carrier_statistics"]["model"]},
        "bandgap_narrowing": {
            "model": bgn["model"],
            # Legacy unit_scaling decks store concentrations in the TCAD
            # internal cm^-3 unit even though the transitional key retains
            # its historical ``_m3`` suffix.
            "reference_doping_m3": bgn["reference_doping_cm3"],
            "coefficient_eV": bgn["coefficient_eV"],
            "smoothing": bgn["smoothing"],
            "offset_eV": bgn["offset_eV"],
            "fermi_statistics_correction": bgn["fermi_statistics_correction"],
        },
        "mobility": mobility_config,
        "recombination": list(recombination["mechanisms"]),
        "taun": recombination["taun_s"],
        "taup": recombination["taup_s"],
        "srh_doping_dependence": {
            "enabled": srh["enabled"],
            "concentration_basis": srh["concentration_basis"],
            "density_coupling": srh["density_coupling"],
            "temperature_dependence": srh["temperature_dependence"],
            "temperature_K": srh["temperature_K"],
            "reference_temperature_K": srh["reference_temperature_K"],
            "electron_temperature_exponent": srh["electron_temperature_exponent"],
            "hole_temperature_exponent": srh["hole_temperature_exponent"],
            "electron": {
                "tau_min_s": srh["electron"]["tau_min_s"],
                "tau_max_s": srh["electron"]["tau_max_s"],
                "reference_doping_m3": srh["electron"]["reference_doping_cm3"],
                "gamma": srh["electron"]["gamma"],
            },
            "hole": {
                "tau_min_s": srh["hole"]["tau_min_s"],
                "tau_max_s": srh["hole"]["tau_max_s"],
                "reference_doping_m3": srh["hole"]["reference_doping_cm3"],
                "gamma": srh["hole"]["gamma"],
            },
        },
        "auger_cn_m6_per_s": recombination["auger_cn_m6_per_s"],
        "auger_cp_m6_per_s": recombination["auger_cp_m6_per_s"],
        "quasi_fermi_update_limit_V": 0.1,
        # The exact-mesh qualification matrix showed that block_filter can
        # force carrier-block decrease after those rows have reached their
        # numerical floor, reducing an otherwise full Newton step to tiny
        # damping factors.  The standard merit globalization is qualified for
        # this template after a fixed-QF Poisson equilibrium bootstrap.
        "line_search_mode": "merit",
        "residual_filter_gamma": 1.0e-4,
        "residual_filter_envelope_factor": 2.0,
        # The exact imported mesh has a measured scaled residual floor near
        # 1e-8--1e-7 on a same-bias equilibrium replay.  Keep this explicit so
        # a line-search stall below the qualified floor is classified, while
        # larger stalls remain hard failures.
        "stall_residual_floor": 1.0e-7,
        "continuity_row_scaling": {
            "enabled": True,
            "flux_fraction": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale": 1.0e-18,
            "min_weight": 1.0e-12,
            "max_weight": 1.0e12,
        },
        "impact_ionization": {"model": "none"},
    }


def contacts(flatband_V: float, gate_V: float, drain_V: float) -> list[dict[str, Any]]:
    return [
        {"name": "gate", "type": "metal_gate", "bias": gate_V,
         "flatband_voltage": flatband_V},
        {"name": "drain", "type": "ohmic", "bias": drain_V},
        {"name": "source", "type": "ohmic", "bias": 0.0},
        {"name": "substrate", "type": "ohmic", "bias": 0.0},
    ]


def sweep_diagnostics(stem: str, *, transport: bool = False) -> dict[str, Any]:
    result: dict[str, Any] = {
        "terminal_balance": {
            "enabled": True,
            "contacts": list(CONTACTS),
            "csv_file": f"{stem}_terminal_balance.csv",
        },
        "srh_balance": {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": list(CONTACTS),
            "resolution_margin_ratio": 10.0,
            "csv_file": f"{stem}_srh_balance.csv",
        },
    }
    if transport:
        result["transport"] = {"enabled": True}
    return result


def dc_deck(
    *, name: str, mesh: Path, doping: Path, materials: Path,
    physics: dict[str, Any], flatband_V: float, gate_V: float, drain_V: float,
    swept_contact: str, bias_points: list[float], initial_state: Path,
    output_dir: Path, high_field: bool, write_every_point: bool = False,
) -> dict[str, Any]:
    sweep: dict[str, Any] = {
        "mode": "iv",
        "contact": swept_contact,
        "current_contact": "drain",
        "start": bias_points[0],
        "stop": bias_points[-1],
        "step": max(abs(bias_points[-1] - bias_points[0]), 1.0),
        "bias_points": bias_points,
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str((output_dir / f"{name}_state.csv").resolve()),
        "write_vtk": False,
        "diagnostics": sweep_diagnostics(name, transport=high_field),
    }
    if write_every_point:
        sweep["write_state_every_point_prefix"] = str(
            (output_dir / f"{name}_point").resolve())
    return {
        "_comment": (
            "Templates/LDMOS phase-2/3 exact-mesh deck. Physics values are "
            "materialized from the versioned WP1.75 contracts."
        ),
        "simulation_type": "dc_sweep",
        "mesh_file": str(mesh.resolve()),
        "node_doping_file": str(doping.resolve()),
        "materials_file": str(materials.resolve()),
        "output_csv": str((output_dir / f"{name}.csv").resolve()),
        "scaling": {"mode": "unit_scaling"},
        "contacts": contacts(flatband_V, gate_V, drain_V),
        "solver": classical_solver(physics, high_field=high_field),
        "sweep": sweep,
    }


def exact_bias_points(curve: Path) -> list[float]:
    with curve.open(newline="", encoding="utf-8") as stream:
        points = [float(row["bias_V"]) for row in csv.DictReader(stream)]
    if not points or any(right <= left for left, right in zip(points, points[1:])):
        raise ValueError("Sentaurus CurrentPlot bias points must be strictly increasing")
    return points


def probe_deck(base: dict[str, Any], simulation_type: str, state: Path,
               output: Path) -> dict[str, Any]:
    result = {key: value for key, value in base.items() if key != "sweep"}
    result["simulation_type"] = simulation_type
    result["state_file"] = str(state.resolve())
    result["output_csv"] = str(output.resolve())
    result["_comment"] = "Fixed-state replay through the Vela production operator."
    return result


def prepare(stage1_dir: Path, oracle_dir: Path, contracts_dir: Path,
            output_dir: Path | None = None) -> dict[str, Any]:
    output = output_dir or stage1_dir / "phase23"
    output.mkdir(parents=True, exist_ok=True)
    exact = stage1_dir / "vela_exact_topology"
    qualification = stage1_dir / "qualification"
    physics = read_json(contracts_dir / "physics_contract.json")
    mapping = derive_polysi_flatband(
        exact / "mesh.json", qualification / "sentaurus_eq_0v_export")
    flatband = float(mapping["vela_flatband_voltage_V"])
    sentaurus_g0_eq = qualification / "sentaurus_eq_0v_state.csv"
    sentaurus_g4_eq = output / "G4_equilibrium_state.csv"
    classical_eq = sentaurus_g4_eq if sentaurus_g4_eq.is_file() else sentaurus_g0_eq
    materials = contracts_dir / physics["materials_file"]

    decks: dict[str, dict[str, Any]] = {}
    decks["g_contact_provisional_eq"] = dc_deck(
        name="g_contact_provisional_eq", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=0.0, gate_V=0.0, drain_V=0.0, swept_contact="drain",
        bias_points=[0.0], initial_state=classical_eq, output_dir=output,
        high_field=False)
    decks["g_contact_polysi_poisson_eq"] = dc_deck(
        name="g_contact_polysi_poisson_eq", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.0,
        swept_contact="drain", bias_points=[0.0], initial_state=classical_eq,
        output_dir=output, high_field=False)
    decks["g_contact_polysi_poisson_eq"]["solver"]["method"] = "poisson_only"
    decks["g_contact_polysi_eq"] = dc_deck(
        name="g_contact_polysi_eq", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.0,
        swept_contact="drain", bias_points=[0.0],
        initial_state=output / "g_contact_polysi_poisson_eq_state.csv",
        output_dir=output, high_field=False)
    # At zero applied bias, the Poisson bootstrap already fixes the physical
    # equilibrium carrier basins.  Keep any coupled cleanup at the serialization
    # noise scale so a depleted minority-QF direction cannot manufacture a
    # spurious contact majority-carrier drop.
    decks["g_contact_polysi_eq"]["solver"]["quasi_fermi_update_limit_V"] = 1.0e-12
    decks["g0_contact_polysi_replay_diagnostic"] = dc_deck(
        name="g0_contact_polysi_replay_diagnostic", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.0,
        swept_contact="drain", bias_points=[0.0], initial_state=sentaurus_g0_eq,
        output_dir=output, high_field=False)
    decks["g_contact_polysi_eq_repeat"] = dc_deck(
        name="g_contact_polysi_eq_repeat", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.0,
        swept_contact="drain", bias_points=[0.0],
        initial_state=output / "g_contact_polysi_eq_state.csv",
        output_dir=output, high_field=False)
    decks["g_contact_polysi_eq_repeat"]["solver"][
        "quasi_fermi_update_limit_V"
    ] = 1.0e-12
    decks["g3_drain_prebias"] = dc_deck(
        name="g3_drain_prebias", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.0,
        swept_contact="drain", bias_points=[0.0, 0.1],
        initial_state=output / "g_contact_polysi_eq_repeat_state.csv",
        output_dir=output, high_field=True)
    decks["g3_drain_prebias_sentaurus_path"] = dc_deck(
        name="g3_drain_prebias_sentaurus_path", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.0,
        swept_contact="drain",
        bias_points=list(SENTAURUS_G3_DRAIN_PREBIAS_POINTS_V),
        initial_state=output / "g_contact_polysi_eq_repeat_state.csv",
        output_dir=output, high_field=True)
    decks["g3_drain_prebias_vd0p023_checkpoint"] = dc_deck(
        name="g3_drain_prebias_vd0p023_checkpoint",
        mesh=exact / "mesh.json", doping=exact / "doping.csv",
        materials=materials, physics=physics, flatband_V=flatband,
        gate_V=0.0, drain_V=0.0, swept_contact="drain",
        bias_points=list(SENTAURUS_G3_VD0P023_CHECKPOINT_POINTS_V),
        initial_state=output / "g_contact_polysi_eq_repeat_state.csv",
        output_dir=output, high_field=True)
    decks["g3_drain_prebias_repeat"] = dc_deck(
        name="g3_drain_prebias_repeat", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.1,
        swept_contact="gate", bias_points=[0.0],
        initial_state=output / "g3_drain_prebias_state.csv",
        output_dir=output, high_field=True)
    # A second, independently qualified route starts from the exact Sentaurus
    # Vg=0/Vd=0.1 state.  It is the production seed when the deep-off drain
    # ramp exposes a continuation defect; the first point must still reclose
    # self-consistently under the frozen G3 contract.
    decks["g3_idvg_seed"] = dc_deck(
        name="g3_idvg_seed", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.1,
        swept_contact="gate", bias_points=[0.0],
        initial_state=qualification / "sentaurus_idvg_vg0_vd0p1_state.csv",
        output_dir=output, high_field=True)
    decks["g3_idvg_seed_repeat"] = dc_deck(
        name="g3_idvg_seed_repeat", mesh=exact / "mesh.json",
        doping=exact / "doping.csv", materials=materials, physics=physics,
        flatband_V=flatband, gate_V=0.0, drain_V=0.1,
        swept_contact="gate", bias_points=[0.0],
        initial_state=output / "g3_idvg_seed_state.csv",
        output_dir=output, high_field=True)
    idvg_points = exact_bias_points(oracle_dir / "IdVg_n2_des_drain_curve.csv")
    decks["g3_idvg"] = dc_deck(
        name="g3_idvg", mesh=exact / "mesh.json", doping=exact / "doping.csv",
        materials=materials, physics=physics, flatband_V=flatband,
        gate_V=0.0, drain_V=0.1, swept_contact="gate",
        bias_points=idvg_points,
        initial_state=output / "g3_drain_prebias_repeat_state.csv",
        output_dir=output, high_field=True, write_every_point=True)

    paths: dict[str, str] = {}
    for name, deck in decks.items():
        path = output / f"{name}.json"
        write_json(path, deck)
        paths[name] = str(path.resolve())

    for simulation_type in (
        "edge_mobility_probe", "sg_edge_flux_probe", "newton_carrier_term_probe"):
        name = f"sentaurus_eq_{simulation_type.removesuffix('_probe')}"
        path = output / f"{name}.json"
        write_json(path, probe_deck(
            decks["g_contact_polysi_eq"], simulation_type, classical_eq,
            output / f"{name}.csv"))
        paths[name] = str(path.resolve())

    manifest = {
        "schema": "vela.templates_ldmos.phase23_decks.v1",
        "physics_layer": "G3 classical bulk plus high-field; no IALMob/QP/thermal/Okuto",
        "gate_mapping": mapping,
        "classical_equilibrium_oracle": {
            "path": str(classical_eq.resolve()),
            "layer": "G4" if classical_eq == sentaurus_g4_eq else "G0_fallback_diagnostic",
            "qualified_for_spatial_scoring": classical_eq == sentaurus_g4_eq,
        },
        "idvg_bias_points_V": idvg_points,
        "decks": paths,
        "limitations": [
            "T-2022.03 parameter values still marked pending in physics_contract.json are not inferred here.",
            "The sealed G0 oracle curve is diagnostic only until a G3 Sentaurus ablation is captured.",
        ],
    }
    write_json(output / "deck_manifest.json", manifest)
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1-dir", type=Path, required=True)
    parser.add_argument("--oracle-dir", type=Path, required=True)
    parser.add_argument("--contracts-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    print(json.dumps(prepare(
        args.stage1_dir.resolve(), args.oracle_dir.resolve(),
        args.contracts_dir.resolve(),
        args.output_dir.resolve() if args.output_dir else None), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
