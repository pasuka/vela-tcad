#!/usr/bin/env python3
"""Prepare a fail-closed Vela SLOT-LDMOS mesh/contact policy bundle.

The input is the neutral CSV directory exported by ``sentaurus_import``.  The
script preserves the TDR topology and nodal doping, classifies gate and SLOT as
electrostatic-only ``metal_gate`` contacts, and selects one atomic avalanche
mesh profile.  The SG/GSS-Laux profile is rejected on any obtuse triangle.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


EXPECTED_CONTACTS = {"drain", "gate", "SLOT", "source", "substrate"}
ELECTROSTATIC_CONTACTS = {"gate", "SLOT"}
MATERIAL_MAP = {
    "Silicon": "Si",
    "Oxide": "SiO2",
    "Si": "Si",
    "SiO2": "SiO2",
    "PolySilicon": "PolySilicon",
    "Nitride": "Nitride",
}
INSULATORS = {"SiO2", "Nitride"}
PROFILE_BUNDLES = {
    "legacy_cell_reconstructed": {
        "mesh_geometry": {
            "node_volume_policy": "barycentric",
            "require_non_obtuse": False,
        },
        "impact_ionization": {
            "current_approximation": "cell_reconstructed",
            "source_mapping_mode": "triangle_gss_gradqf_truncated",
            "cell_reconstructed_midpoint_density": "gss_logistic",
        },
    },
    "element_edge_sg_gss_laux": {
        "mesh_geometry": {
            "node_volume_policy": "mixed_voronoi",
            "require_non_obtuse": True,
        },
        "impact_ionization": {
            "current_approximation": "element_edge_sg_gss_laux",
            "source_mapping_mode": "element_vertex_box_measure",
            "cell_reconstructed_midpoint_density": "bernoulli",
        },
    },
}

SENTAURUS_SERIES_RESISTANCE_OHM_UM = 1.0e12
SENTAURUS_BREAKDOWN_CURRENT_A_PER_UM = 1.0e-7
SENTAURUS_REFERENCE_BVDS_V = 38.520901203613384
SLOT_LDMOS_RESTART_RESIDUAL_SCALES = {
    "psi": 1.189625,
    "phin": 8476.074054,
    "phip": 3518.300497,
}
VELA_UNIT_RESISTOR_INNER_V = 0.008374398259206585


class PreparationError(ValueError):
    """Raised when the imported structure violates the SLOT/BVDS contract."""


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_ids(value: str) -> list[int]:
    return [int(item) for item in value.replace("|", ";").split(";") if item]


def triangle_angles(points: Iterable[tuple[float, float]]) -> list[float]:
    p = list(points)
    result: list[float] = []
    for index in range(3):
        origin = p[index]
        left = p[(index + 1) % 3]
        right = p[(index + 2) % 3]
        u = (left[0] - origin[0], left[1] - origin[1])
        v = (right[0] - origin[0], right[1] - origin[1])
        denom = math.hypot(*u) * math.hypot(*v)
        if denom == 0.0:
            raise PreparationError("mesh contains a zero-length triangle edge")
        cosine = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / denom))
        result.append(math.degrees(math.acos(cosine)))
    return result


def mesh_quality(nodes: list[dict[str, Any]], triangles: list[dict[str, Any]]) -> dict[str, Any]:
    coords = {node["id"]: (node["x"], node["y"]) for node in nodes}
    maximum_angle = 0.0
    minimum_area = math.inf
    maximum_aspect = 0.0
    obtuse: list[int] = []
    runtime_rejected: list[int] = []
    for cell in triangles:
        points = [coords[node_id] for node_id in cell["node_ids"]]
        twice_area = abs(
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[1][1] - points[0][1]) * (points[2][0] - points[0][0])
        )
        if twice_area == 0.0:
            raise PreparationError(f"mesh contains zero-area triangle {cell['id']}")
        minimum_area = min(minimum_area, 0.5 * twice_area)
        lengths = [
            math.dist(points[0], points[1]),
            math.dist(points[1], points[2]),
            math.dist(points[2], points[0]),
        ]
        maximum_aspect = max(maximum_aspect, max(lengths) / min(lengths))
        cell_maximum = max(triangle_angles(points))
        maximum_angle = max(maximum_angle, cell_maximum)
        if cell_maximum > 90.0:
            obtuse.append(cell["id"])
        if cell_maximum > 90.0 + 1.0e-10:
            runtime_rejected.append(cell["id"])
    return {
        "cell_count": len(triangles),
        "minimum_area_um2": minimum_area,
        "maximum_edge_aspect_ratio": maximum_aspect,
        "maximum_triangle_angle_deg": maximum_angle,
        "obtuse_cell_count": len(obtuse),
        "runtime_non_obtuse_rejected_cell_count": len(runtime_rejected),
        "non_obtuse": not runtime_rejected,
    }


def load_neutral_export(
    input_dir: Path, coordinate_scale_to_um: float
) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    node_rows = read_csv(input_dir / "nodes.csv")
    element_rows = read_csv(input_dir / "elements.csv")
    contact_rows = read_csv(input_dir / "contacts.csv")
    raw_coordinates = [
        (float(row["x_um"]), float(row["y_um"])) for row in node_rows
    ]
    nodes = [
        {
            "id": int(row["id"]),
            "x": raw_coordinates[index][0] * coordinate_scale_to_um,
            "y": raw_coordinates[index][1] * coordinate_scale_to_um,
        }
        for index, row in enumerate(node_rows)
    ]
    if [node["id"] for node in nodes] != list(range(len(nodes))):
        raise PreparationError("node ids must be zero-based and contiguous")

    region_cells: dict[str, list[int]] = defaultdict(list)
    region_material: dict[str, str] = {}
    triangles: list[dict[str, Any]] = []
    for expected_id, row in enumerate(element_rows):
        cell_id = int(row["id"])
        if cell_id != expected_id:
            raise PreparationError("element ids must be zero-based and contiguous")
        region = row["region"]
        raw_material = row["material"]
        if raw_material not in MATERIAL_MAP:
            raise PreparationError(f"unsupported material {raw_material!r}")
        material = MATERIAL_MAP[raw_material]
        previous = region_material.setdefault(region, material)
        if previous != material:
            raise PreparationError(f"region {region!r} has inconsistent materials")
        region_cells[region].append(cell_id)
        triangles.append({
            "id": cell_id,
            "region_id": None,
            "node_ids": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
        })

    region_ids = {name: index for index, name in enumerate(region_cells)}
    for triangle, row in zip(triangles, element_rows, strict=True):
        triangle["region_id"] = region_ids[row["region"]]
    regions = [
        {
            "id": region_ids[name],
            "name": name,
            "material": region_material[name],
            "cell_ids": cell_ids,
        }
        for name, cell_ids in region_cells.items()
    ]

    names = {row["name"] for row in contact_rows}
    if names != EXPECTED_CONTACTS:
        raise PreparationError(
            f"expected contacts {sorted(EXPECTED_CONTACTS)}, got {sorted(names)}"
        )
    incident_counts: dict[str, Counter[str]] = {}
    contacts: list[dict[str, Any]] = []
    for contact_id, row in enumerate(contact_rows):
        name = row["name"]
        owner = row["region"]
        if owner not in region_ids:
            raise PreparationError(f"contact {name!r} references unknown region {owner!r}")
        node_ids = sorted(set(parse_ids(row["node_ids"])))
        wanted = set(node_ids)
        counts: Counter[str] = Counter()
        for triangle in triangles:
            if wanted.intersection(triangle["node_ids"]):
                material = regions[triangle["region_id"]]["material"]
                counts[material] += 1
        if not counts:
            raise PreparationError(f"contact {name!r} touches no material cell")
        incident_counts[name] = counts
        contacts.append({
            "id": contact_id,
            "name": name,
            "region_id": region_ids[owner],
            "node_ids": node_ids,
        })

    slot_materials = set(incident_counts["SLOT"])
    if not slot_materials or not slot_materials.issubset(INSULATORS):
        raise PreparationError(
            "SLOT must be a collapsed conductor boundary touching only insulators; "
            f"got {sorted(slot_materials)}"
        )

    mesh = {
        "_comment": "Exact topology imported from the sealed SLOT-LDMOS Sentaurus TDR.",
        "nodes": nodes,
        "triangles": triangles,
        "regions": regions,
        "contacts": contacts,
    }
    contact_audit = {
        name: {
            "node_count": len(next(item for item in contacts if item["name"] == name)["node_ids"]),
            "owner_region": next(row["region"] for row in contact_rows if row["name"] == name),
            "incident_material_cell_counts": dict(sorted(counts.items())),
            "vela_type": "metal_gate" if name in ELECTROSTATIC_CONTACTS else "ohmic",
            "carrier_dirichlet": name not in ELECTROSTATIC_CONTACTS,
        }
        for name, counts in incident_counts.items()
    }
    coordinate_audit = {
        "input_column_labels": ["x_um", "y_um"],
        "input_physical_unit": "cm",
        "output_physical_unit": "um",
        "scale_to_um": coordinate_scale_to_um,
        "input_bounds": {
            "x_min": min(item[0] for item in raw_coordinates),
            "x_max": max(item[0] for item in raw_coordinates),
            "y_min": min(item[1] for item in raw_coordinates),
            "y_max": max(item[1] for item in raw_coordinates),
        },
        "output_bounds_um": {
            "x_min": min(node["x"] for node in nodes),
            "x_max": max(node["x"] for node in nodes),
            "y_min": min(node["y"] for node in nodes),
            "y_max": max(node["y"] for node in nodes),
        },
    }
    return mesh, contact_audit, coordinate_audit


def contact_specs(mesh: dict[str, Any]) -> list[dict[str, Any]]:
    specs: list[dict[str, Any]] = []
    for contact in mesh["contacts"]:
        name = contact["name"]
        spec: dict[str, Any] = {"name": name, "bias": 0.0}
        spec["type"] = "metal_gate" if name in ELECTROSTATIC_CONTACTS else "ohmic"
        if name in ELECTROSTATIC_CONTACTS:
            spec["flatband_voltage"] = 0.0
        specs.append(spec)
    return specs


def region_doping(input_dir: Path, mesh: dict[str, Any]) -> list[dict[str, Any]]:
    doping_rows = {
        int(row["node_id"]): (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
        for row in read_csv(input_dir / "doping.csv")
    }
    result: list[dict[str, Any]] = []
    for region in mesh["regions"]:
        if region["material"] != "Si":
            result.append({
                "region": region["name"],
                "donors": 0.0,
                "acceptors": 0.0,
            })
            continue
        node_ids: set[int] = set()
        for cell_id in region["cell_ids"]:
            node_ids.update(mesh["triangles"][cell_id]["node_ids"])
        samples = [doping_rows.get(node_id, (0.0, 0.0)) for node_id in node_ids]
        result.append({
            "region": region["name"],
            "donors": sum(item[0] for item in samples) / max(1, len(samples)),
            "acceptors": sum(item[1] for item in samples) / max(1, len(samples)),
        })
    return result


def write_transport_doping(
    input_path: Path, output_path: Path, mesh: dict[str, Any]
) -> dict[str, Any]:
    silicon_nodes: set[int] = set()
    for region in mesh["regions"]:
        if region["material"] != "Si":
            continue
        for cell_id in region["cell_ids"]:
            silicon_nodes.update(mesh["triangles"][cell_id]["node_ids"])

    rows = read_csv(input_path)
    if not rows:
        raise PreparationError("doping.csv must contain at least one node row")
    fieldnames = list(rows[0])
    required = {"node_id", "donors_cm3", "acceptors_cm3"}
    if not required.issubset(fieldnames):
        raise PreparationError(
            f"doping.csv is missing columns {sorted(required - set(fieldnames))}"
        )

    zeroed_nodes = 0
    zeroed_nonzero_nodes = 0
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            node_id = int(row["node_id"])
            if node_id not in silicon_nodes:
                donors = float(row["donors_cm3"])
                acceptors = float(row["acceptors_cm3"])
                zeroed_nodes += 1
                if donors != 0.0 or acceptors != 0.0:
                    zeroed_nonzero_nodes += 1
                row["donors_cm3"] = "0"
                row["acceptors_cm3"] = "0"
            writer.writerow(row)
    return {
        "policy": "retain doping only on nodes incident to Silicon cells",
        "silicon_transport_node_count": len(silicon_nodes),
        "nontransport_node_count": len(rows) - len(silicon_nodes),
        "zeroed_node_count": zeroed_nodes,
        "zeroed_previously_nonzero_node_count": zeroed_nonzero_nodes,
    }


def materials_document() -> dict[str, Any]:
    return {
        "_comment": (
            "Nitride and PolySilicon entries required to load the imported process mesh. "
            "PolySilicon is electrostatic-only: omitted ni/mobility/DOS fields make its "
            "carrier rows use the same nontransport gauge semantics as SiO2."
        ),
        "materials": [
            {
                "name": "PolySilicon",
                "eps_r": 11.7,
                "bandgap_eV": 1.12,
                "electron_affinity_eV": 4.05,
                "temperature_K": 300.0,
            },
            {
                "name": "Nitride",
                "eps_r": 7.5,
                "bandgap_eV": 5.0,
                "electron_affinity_eV": 1.9,
                "temperature_K": 300.0,
            },
        ],
    }


def avalanche_config(profile: str) -> dict[str, Any]:
    config: dict[str, Any] = {
        "model": "van_overstraeten",
        "coupling_mode": "self_consistent",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "current_magnitude_mode": "edge_scalar_abs",
        "driving_force_interpolation": "none",
        "edge_source_partition": "symmetric",
        "electron_driving_force_ref_density_m3": 0.0,
        "hole_driving_force_ref_density_m3": 0.0,
        "minimum_field_V_m": 0.0,
        "quasi_fermi_carrier_truncation": 0.0,
        "quasi_fermi_gradient_discretization": "cell_gradient",
        "source_geometry_scale": 1.0,
        "source_volume_factor": 0.0,
        "source_volume_policy": "genius_truncated",
    }
    config.update(PROFILE_BUNDLES[profile]["impact_ionization"])
    return config


def solver_config(profile: str, avalanche_coupling: str | None) -> dict[str, Any]:
    solver: dict[str, Any] = {
        "method": "gummel_newton",
        "max_iter": 80,
        "reltol": 1.0e-8,
        "damping_psi": 0.2,
        "bandgap_narrowing": "old_slotboom",
        "mobility": {
            "model": "masetti_field",
            "doping_concentration_basis": "total_impurity",
            "high_field_driving_force": "quasi_fermi_gradient",
        },
        "recombination": ["srh", "auger"],
    }
    if avalanche_coupling is not None:
        if avalanche_coupling not in {"postprocess_only", "self_consistent"}:
            raise PreparationError(
                f"unsupported avalanche coupling mode {avalanche_coupling!r}"
            )
        impact = avalanche_config(profile)
        impact["coupling_mode"] = avalanche_coupling
        solver["impact_ionization"] = impact
    return solver


def external_resistor_config(
    initial_inner_voltage_V: float, max_inner_voltage_step_V: float
) -> dict[str, Any]:
    return {
        "mode": "series_resistor",
        "resistance_ohm_um": SENTAURUS_SERIES_RESISTANCE_OHM_UM,
        "current_direction": 1.0,
        "initial_inner_voltage_V": initial_inner_voltage_V,
        "residual_tolerance_V": 1.0e-6,
        "voltage_tolerance_V": 1.0e-8,
        "max_inner_voltage_step_V": max_inner_voltage_step_V,
        "max_bracket_steps": 400,
        "max_iterations": 60,
    }


def stage_sweep(
    stage_id: str,
    bias_points: list[float],
    *,
    initial_state_file: str | None,
    external_resistor: bool,
    initial_inner_voltage_V: float = 0.0,
    max_inner_voltage_step_V: float = 0.01,
    predictor_max_step_factor: float = 2.0,
) -> dict[str, Any]:
    if not bias_points:
        raise PreparationError(f"stage {stage_id!r} must contain at least one bias point")
    sweep: dict[str, Any] = {
        "mode": "bv_reverse",
        "contact": "drain",
        "current_contact": "drain",
        "start": bias_points[0],
        "stop": bias_points[-1],
        "step": 1.0,
        "bias_points": bias_points,
        "initial_step": 1.0e-4,
        "min_step": 1.0e-10,
        "max_step": 1.0,
        "growth_factor": 1.2,
        "shrink_factor": 0.5,
        "max_retries": 20,
        "stop_on_failure": True,
        "write_vtk": False,
        "write_state_file": f"outputs/stages/{stage_id}/final_state.h5",
        "write_state_every_point_prefix": (
            f"outputs/stages/{stage_id}/states/state"
        ),
        "breakdown": {
            "max_electric_field_V_per_m": 1.0e12,
            "current_jump_ratio": 1.0e12,
            "non_convergence": True,
        },
    }
    if initial_state_file is not None:
        sweep["initial_state_file"] = initial_state_file
    if external_resistor:
        sweep["external_circuit"] = external_resistor_config(
            initial_inner_voltage_V, max_inner_voltage_step_V
        )
        sweep["boundary_control"] = {
            "evaluation_csv": (
                f"outputs/stages/{stage_id}/boundary_control_evaluations.csv"
            ),
            "checkpoint_directory": (
                f"outputs/stages/{stage_id}/boundary_control_checkpoints"
            ),
            "resume": True,
            "predictor_max_step_factor": predictor_max_step_factor,
            "preferred_max_evaluations": 3,
        }
    return sweep


def stage_document(
    common: dict[str, Any],
    profile: str,
    stage_id: str,
    purpose: str,
    bias_points: list[float],
    *,
    avalanche_coupling: str | None,
    initial_state_file: str | None,
    external_resistor: bool,
    initial_inner_voltage_V: float = 0.0,
    max_inner_voltage_step_V: float = 0.01,
    predictor_max_step_factor: float = 2.0,
) -> dict[str, Any]:
    solver = solver_config(profile, avalanche_coupling)
    if initial_state_file is not None:
        # Restarted stages already have a physically consistent state.  A fresh
        # Gummel pre-pass can stall at the tiny inner-voltage updates selected
        # by a 1e12 ohm*um load line, so hand the restart directly to Newton.
        solver.update({
            "abstol": 1.0e-9,
            "line_search": True,
            "warm_start": True,
            "residual_scales": SLOT_LDMOS_RESTART_RESIDUAL_SCALES,
            "handoff": {
                "fallback": "none",
                "require_gummel_convergence": False,
                "gummel_max_iter": 0,
                "newton_max_iter": 80,
            },
        })
    if avalanche_coupling == "self_consistent":
        solver["handoff"]["gummel_max_iter"] = 50
        # Sentaurus includes avalanche derivatives unless -AvalDerivatives is
        # explicitly selected.  Use the complete nine-column local derivative
        # of the legacy triangle GSS source for the production baseline.
        solver["impact_ionization"]["source_jacobian"] = "local_ad"
        solver.update({
            "carrier_regularization_scale": 1.0e-8,
            "quasi_fermi_update_limit_V": 0.1,
            "quasi_fermi_update_limit_minority_V": 0.05,
        })
    sweep = stage_sweep(
        stage_id,
        bias_points,
        initial_state_file=initial_state_file,
        external_resistor=external_resistor,
        initial_inner_voltage_V=initial_inner_voltage_V,
        max_inner_voltage_step_V=max_inner_voltage_step_V,
        predictor_max_step_factor=predictor_max_step_factor,
    )
    if avalanche_coupling is not None:
        sweep["diagnostics"] = {
            "release_bv_config_audit": {
                "enabled": True,
                "csv_file": f"outputs/stages/{stage_id}/avalanche_summary.csv",
                "summary_file": (
                    f"outputs/stages/{stage_id}/avalanche_summary.md"
                ),
            }
        }
    return {
        "_comment": purpose,
        "simulation_type": "dc_sweep",
        "output_csv": f"outputs/stages/{stage_id}/iv.csv",
        **copy.deepcopy(common),
        "solver": solver,
        "sweep": sweep,
    }


def write_json(path: Path, document: dict[str, Any]) -> None:
    path.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(
    input_dir: Path,
    output_dir: Path,
    profile: str,
    stop_voltage: float,
    coordinate_scale_to_um: float = 1.0e4,
) -> dict[str, Any]:
    if not math.isfinite(coordinate_scale_to_um) or coordinate_scale_to_um <= 0.0:
        raise PreparationError("coordinate scale to um must be positive and finite")
    mesh, contact_audit, coordinate_audit = load_neutral_export(
        input_dir, coordinate_scale_to_um
    )
    quality = mesh_quality(mesh["nodes"], mesh["triangles"])
    if profile == "element_edge_sg_gss_laux" and not quality["non_obtuse"]:
        raise PreparationError(
            "element_edge_sg_gss_laux is forbidden on this mesh: "
            f"{quality['obtuse_cell_count']} obtuse cells, maximum angle "
            f"{quality['maximum_triangle_angle_deg']:.9g} degrees"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "mesh.json", mesh)
    doping_audit = write_transport_doping(
        input_dir / "doping.csv", output_dir / "doping.csv", mesh
    )
    write_json(output_dir / "materials.json", materials_document())

    contacts = contact_specs(mesh)
    doping = region_doping(output_dir, mesh)
    common = {
        "mesh_file": "mesh.json",
        "materials_file": "materials.json",
        "node_doping_file": "doping.csv",
        "scaling": {"mode": "unit_scaling"},
        "doping": doping,
        "contacts": contacts,
        "mesh_geometry": PROFILE_BUNDLES[profile]["mesh_geometry"],
    }
    poisson = {
        "_comment": "SLOT electrostatic-boundary and obtuse-mesh fallback check.",
        "simulation_type": "poisson",
        "output_vtk": "outputs/slot_boundary_check.vtk",
        **common,
    }
    write_json(output_dir / "simulation_slot_boundary_check.json", poisson)

    if not math.isfinite(stop_voltage) or stop_voltage <= 0.0:
        raise PreparationError("stop voltage must be positive and finite")

    low_voltage_points = [
        value
        for value in [0.0, 1.0, 10.0, 20.0, 30.0, 35.0, 37.0, 38.0,
                      38.25, 38.4, 38.5, 39.0, 40.0, 50.0, stop_voltage]
        if value <= stop_voltage
    ]
    if low_voltage_points[-1] != stop_voltage:
        low_voltage_points.append(stop_voltage)
    low_voltage_points = sorted(set(low_voltage_points))

    equilibrium_state = "outputs/stages/00_equilibrium/final_state.h5"
    unit_resistor_state = "outputs/stages/01_unit_resistor_1v/final_state.h5"
    avalanche_activation_state = (
        "outputs/stages/04_avalanche_activation_1v/final_state.h5"
    )
    avalanche_state = "outputs/stages/05_avalanche_on_60v/final_state.h5"
    threshold_outer_voltage = (
        SENTAURUS_REFERENCE_BVDS_V
        + SENTAURUS_SERIES_RESISTANCE_OHM_UM
        * SENTAURUS_BREAKDOWN_CURRENT_A_PER_UM
    )
    final_outer_voltage = (
        SENTAURUS_REFERENCE_BVDS_V
        + 1.1
        * SENTAURUS_SERIES_RESISTANCE_OHM_UM
        * SENTAURUS_BREAKDOWN_CURRENT_A_PER_UM
    )
    final_outer_points = [
        stop_voltage,
        1.0e3,
        1.0e4,
        2.5e4,
        5.0e4,
        7.5e4,
        9.0e4,
        9.9e4,
        threshold_outer_voltage,
        final_outer_voltage,
    ]
    final_outer_points = sorted(
        set(value for value in final_outer_points if value >= stop_voltage)
    )

    stage_specs = [
        {
            "id": "00_equilibrium",
            "filename": "simulation_00_equilibrium.json",
            "purpose": (
                "Zero-bias full drift-diffusion equilibrium; avalanche and the "
                "external resistor are disabled to produce the common initial state."
            ),
            "depends_on": [],
            "bias_points": [0.0],
            "avalanche_coupling": None,
            "initial_state_file": None,
            "external_resistor": False,
            "initial_inner_voltage_V": 0.0,
            "max_inner_voltage_step_V": 0.01,
            "predictor_max_step_factor": 2.0,
        },
        {
            "id": "01_unit_resistor_1v",
            "filename": "simulation_01_unit_resistor_1v.json",
            "purpose": (
                "One-volt load-line smoke test for Vouter = Vinner + R*Id with "
                "R = 1e12 ohm*um; avalanche is disabled."
            ),
            "depends_on": ["00_equilibrium"],
            "bias_points": [0.0, 0.1, 0.2, 0.5, 1.0],
            "avalanche_coupling": None,
            "initial_state_file": equilibrium_state,
            "external_resistor": True,
            "initial_inner_voltage_V": 0.0,
            "max_inner_voltage_step_V": 0.01,
            "predictor_max_step_factor": 2.0,
        },
        {
            "id": "02_avalanche_off_60v",
            "filename": "simulation_02_avalanche_off_60v.json",
            "purpose": (
                "Leakage/load-line baseline through the low-voltage range with "
                "impact ionization disabled."
            ),
            "depends_on": ["00_equilibrium"],
            "bias_points": low_voltage_points,
            "avalanche_coupling": None,
            "initial_state_file": equilibrium_state,
            "external_resistor": True,
            "initial_inner_voltage_V": 0.0,
            "max_inner_voltage_step_V": 0.01,
            "predictor_max_step_factor": 2.0,
        },
        {
            "id": "03_iic_postprocess_60v",
            "filename": "simulation_03_iic_postprocess_60v.json",
            "purpose": (
                "Frozen-feedback ionization-integral diagnostic: avalanche is "
                "evaluated but not coupled into the carrier equations."
            ),
            "depends_on": ["00_equilibrium"],
            "bias_points": low_voltage_points,
            "avalanche_coupling": "postprocess_only",
            "initial_state_file": equilibrium_state,
            "external_resistor": True,
            "initial_inner_voltage_V": 0.0,
            "max_inner_voltage_step_V": 0.01,
            "predictor_max_step_factor": 2.0,
        },
        {
            "id": "04_avalanche_activation_1v",
            "filename": "simulation_04_avalanche_activation_1v.json",
            "purpose": (
                "Activate the self-consistent avalanche Jacobian at the verified "
                "1 V load-line state using a 1e-4 V maximum inner-voltage step."
            ),
            "depends_on": ["01_unit_resistor_1v"],
            "bias_points": [1.0],
            "avalanche_coupling": "self_consistent",
            "initial_state_file": unit_resistor_state,
            "external_resistor": True,
            "initial_inner_voltage_V": VELA_UNIT_RESISTOR_INNER_V,
            "max_inner_voltage_step_V": 1.0e-4,
            "predictor_max_step_factor": 2.0,
        },
        {
            "id": "05_avalanche_on_60v",
            "filename": "simulation_05_avalanche_on_60v.json",
            "purpose": (
                "Self-consistent avalanche sweep from the activated 1 V state "
                "through the intrinsic-knee region with the series resistor enabled."
            ),
            "depends_on": ["04_avalanche_activation_1v"],
            "bias_points": [value for value in low_voltage_points if value >= 10.0],
            "avalanche_coupling": "self_consistent",
            "initial_state_file": avalanche_activation_state,
            "external_resistor": True,
            "initial_inner_voltage_V": VELA_UNIT_RESISTOR_INNER_V,
            "max_inner_voltage_step_V": 0.01,
            "predictor_max_step_factor": 2.0,
        },
        {
            "id": "06_bvds_external_resistor_final",
            "filename": "simulation_06_bvds_external_resistor_final.json",
            "purpose": (
                "Final self-consistent load-line continuation in outer voltage, "
                "including the Sentaurus 1e-7 A/um threshold-equivalent point."
            ),
            "depends_on": ["05_avalanche_on_60v"],
            "bias_points": final_outer_points,
            "avalanche_coupling": "self_consistent",
            "initial_state_file": avalanche_state,
            "external_resistor": True,
            "initial_inner_voltage_V": SENTAURUS_REFERENCE_BVDS_V,
            "max_inner_voltage_step_V": 0.01,
            "predictor_max_step_factor": 2.0,
        },
    ]

    manifest_stages: list[dict[str, Any]] = []
    stage_documents: dict[str, dict[str, Any]] = {}
    for spec in stage_specs:
        stage_output_dir = output_dir / "outputs" / "stages" / spec["id"]
        (stage_output_dir / "states").mkdir(parents=True, exist_ok=True)
        if spec["external_resistor"]:
            (stage_output_dir / "boundary_control_checkpoints").mkdir(
                parents=True, exist_ok=True
            )
        document = stage_document(
            common,
            profile,
            spec["id"],
            spec["purpose"],
            spec["bias_points"],
            avalanche_coupling=spec["avalanche_coupling"],
            initial_state_file=spec["initial_state_file"],
            external_resistor=spec["external_resistor"],
            initial_inner_voltage_V=spec["initial_inner_voltage_V"],
            max_inner_voltage_step_V=spec["max_inner_voltage_step_V"],
            predictor_max_step_factor=spec["predictor_max_step_factor"],
        )
        stage_documents[spec["id"]] = document
        config_path = output_dir / spec["filename"]
        write_json(config_path, document)
        manifest_stages.append({
            "id": spec["id"],
            "config": spec["filename"],
            "purpose": spec["purpose"],
            "depends_on": spec["depends_on"],
            "output_csv": document["output_csv"],
            "final_state": document["sweep"]["write_state_file"],
            "sha256": file_sha256(config_path),
        })

    # Backward-compatible filename: now points at the external-resistor,
    # self-consistent low-voltage stage instead of the old direct-voltage deck.
    legacy_alias = copy.deepcopy(stage_documents["05_avalanche_on_60v"])
    legacy_alias["_comment"] = (
        "Compatibility alias for simulation_05_avalanche_on_60v.json. The "
        "authoritative execution order is slot_ldmos_bvds_stages_manifest.json."
    )
    write_json(output_dir / "simulation_bvds_legacy.json", legacy_alias)

    stages_manifest = {
        "schema": "vela.slot_ldmos_bvds_stages.v1",
        "selected_profile": profile,
        "execution_order": [stage["id"] for stage in manifest_stages],
        "external_circuit_contract": {
            "equation": "Vouter_V = Vinner_V + resistance_ohm_um * Id_A_per_um",
            "current_direction": 1.0,
            "resistance_ohm_um": SENTAURUS_SERIES_RESISTANCE_OHM_UM,
            "current_unit": "A/um",
            "resistance_unit": "ohm*um",
            "voltage_unit": "V",
        },
        "reference_stop_contract": {
            "sentaurus_reference_bvds_V": SENTAURUS_REFERENCE_BVDS_V,
            "breakdown_current_A_per_um": SENTAURUS_BREAKDOWN_CURRENT_A_PER_UM,
            "threshold_outer_voltage_V": threshold_outer_voltage,
            "final_outer_voltage_V": final_outer_voltage,
            "criterion": "first converged point with abs(Id) >= breakdown current",
        },
        "restart_solver_contract": {
            "strategy": "equilibrium state direct to Newton",
            "fixed_residual_scales": SLOT_LDMOS_RESTART_RESIDUAL_SCALES,
            "scale_source": (
                "raw block norms from the first converged equilibrium-restart "
                "load-line device solve"
            ),
        },
        "stages": manifest_stages,
    }
    write_json(output_dir / "slot_ldmos_bvds_stages_manifest.json", stages_manifest)

    report = {
        "schema": "vela.slot_ldmos_bvds_preparation.v1",
        "selected_profile": profile,
        "exact_tdr_topology_preserved": True,
        "coordinate_conversion": coordinate_audit,
        "contact_policy": contact_audit,
        "nontransport_region_policy": {
            "materials": ["PolySilicon", "SiO2", "Nitride"],
            "carrier_transport": False,
            "doping": doping_audit,
        },
        "mesh_quality": quality,
        "mesh_geometry": PROFILE_BUNDLES[profile]["mesh_geometry"],
        "impact_ionization_current_support": PROFILE_BUNDLES[profile]["impact_ionization"],
        "staged_bvds_manifest": "slot_ldmos_bvds_stages_manifest.json",
        "external_circuit_contract": stages_manifest["external_circuit_contract"],
        "forbidden_profiles": (
            ["element_edge_sg_gss_laux"] if not quality["non_obtuse"] else []
        ),
        "decision": (
            "Use positive barycentric control volumes and the box builder's "
            "negative-cotangent fallback on the exact TDR mesh. Do not claim the "
            "PN2D non-obtuse SG/GSS-Laux qualification for this device."
        ),
    }
    write_json(output_dir / "slot_ldmos_policy_report.json", report)
    return report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--profile",
        choices=sorted(PROFILE_BUNDLES),
        default="legacy_cell_reconstructed",
    )
    parser.add_argument("--stop-voltage", type=float, default=60.0)
    parser.add_argument(
        "--coordinate-scale-to-um",
        type=float,
        default=1.0e4,
        help=(
            "Scale applied to neutral-export coordinates. The sealed SProcess TDR "
            "stores centimetres, so its required cm-to-um value is 1e4."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        report = prepare(
            args.input_dir,
            args.output_dir,
            args.profile,
            args.stop_voltage,
            args.coordinate_scale_to_um,
        )
    except PreparationError as error:
        print(f"prepare_slot_ldmos_bvds: {error}", file=sys.stderr)
        return 2
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
