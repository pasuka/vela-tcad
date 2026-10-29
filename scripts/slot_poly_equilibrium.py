"""Guarded candidate restart for isolated, source-free SLOT PolySilicon.

This operation supplies an initial state only. The caller must still run the
original coupled convergence, update-history and all-active continuity gates.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slot_state_archive as states


_REQUIRED = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
_SPLIT = [f"{carrier}_qf_{part}_V" for carrier in ("electron", "hole")
          for part in ("increment", "reference")]
_QUANTUM = ["electron_quantum_potential_V", "electron_quantum_potential_like_V"]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError("Poly equilibrium restart: " + message)


def _finite(value, label: str) -> float:
    _require(not isinstance(value, bool), label + " must be finite numeric data")
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Poly equilibrium restart: invalid " + label) from exc
    _require(math.isfinite(number), label + " must be finite")
    return number


def _model(value, default="none"):
    return value.get("model", default) if isinstance(value, dict) else value


def _profile(document: dict) -> None:
    solver = document.get("solver", {})
    _require(isinstance(solver, dict), "solver must be an object")
    _require(solver.get("method", solver.get("type")) == "newton", "requires Newton SLOT profile")
    _require(_model(solver.get("carrier_statistics", "boltzmann"), "boltzmann") == "boltzmann"
             and solver.get("statistics", "boltzmann") == "boltzmann", "requires Boltzmann statistics")
    _require(_finite(solver.get("temperature_K", 300), "temperature_K") == 300, "requires 300 K")
    quantum = solver.get("electron_quantum_potential", False)
    _require((quantum.get("enabled", False) if isinstance(quantum, dict) else quantum) is False,
             "active quantum model is unsupported")
    scope = solver.get("recombination_materials")
    _require(isinstance(scope, list) and bool(scope) and all(isinstance(x, str) for x in scope)
             and "PolySilicon" not in scope, "requires explicit recombination scope excluding PolySilicon")
    _require(_model(solver.get("band_to_band", "none")) == "none", "band-to-band generation is unsupported")
    for owner in (document, solver):
        _require(not any(key in owner for key in ("generation", "optical_generation", "btbt",
                                                  "thermal", "electrothermal")),
                 "unexpected generation or thermal configuration")
    _require(document.get("simulation_type", "dc_sweep") in ("dc_sweep", "newton"),
             "unsupported simulation type")
    _require(not document.get("boundaries"), "explicit boundaries are unsupported")
    circuit = document.get("sweep", {}).get("external_circuit", {})
    _require(isinstance(circuit, dict) and circuit.get("enabled", False) is False,
             "external circuit bias is unsupported; requires direct contact bias")
    mobility = solver.get("mobility")
    _require(isinstance(mobility, dict) and mobility.get("model") == "ialmob", "requires IALMob")
    ialmob = mobility.get("ialmob", {})
    _require(isinstance(ialmob, dict), "IALMob configuration must be an object")
    overrides = ialmob.get("constant_material_mobility_cm2_per_Vs", {})
    _require(isinstance(overrides, dict), "constant mobility overrides must be an object")
    poly = overrides.get("PolySilicon", {})
    _require(isinstance(poly, dict) and all(k in poly for k in ("electron", "hole")),
             "requires explicit constant PolySilicon electron and hole mobility")
    for carrier in ("electron", "hole"):
        _require(_finite(poly[carrier], "PolySilicon " + carrier + " mobility") > 0,
                 "PolySilicon mobility must be positive")
    av = solver.get("impact_ionization", "none")
    if _model(av) == "none":
        return
    _require(isinstance(av, dict), "avalanche requires explicit canonical SG/GradQF configuration")
    expected = {"model": "van_overstraeten", "generation": "current_density",
                "driving_force": "quasi_fermi_gradient",
                "current_approximation": "element_edge_sg_gss_laux",
                "source_mapping_mode": "element_vertex_box_measure",
                "cell_reconstructed_midpoint_density": "bernoulli", "source_jacobian": "local_ad"}
    if av.get("current_approximation") == "ialmob_element_density_grad_qf":
        # AssemblerUtils::sgAvalancheSourceComponentIntegrals explicitly visits
        # Si cells only for this profile. _components additionally rejects any
        # Poly node shared with Si/another transport material. Thus the contact
        # electric-field fallback cannot introduce a Poly source here.
        expected.update(current_approximation="ialmob_element_density_grad_qf",
                        quasi_fermi_gradient_discretization="cell_gradient",
                        source_volume_policy="genius_truncated", source_volume_factor=0.,
                        source_geometry_scale=1., ialmob_boundary_field_projection=True,
                        contact_electric_field_fallback=True)
        _require(all(av.get(key) == value for key, value in expected.items()),
                 "observed Si-only avalanche requires its explicit canonical profile")
        _require(mobility.get("jacobian_field_derivatives") is True,
                 "observed Si-only avalanche requires mobility derivatives")
        _require(av.get("contact_electric_field_fallback_scope", "contact_node_cell") == "contact_node_cell"
                 and av.get("contact_electric_field_fallback_mode", "cell_gradient_magnitude") == "cell_gradient_magnitude"
                 and not av.get("debug_raw_van_overstraeten", False),
                 "unsupported observed avalanche contact/raw-field profile")
        for key in ("temperature_K", "reference_temperature_K"):
            _require(_finite(av.get(key, 300), key) == 300, "observed avalanche requires 300 K")
    else:
        _require(all(av.get(key) == value for key, value in expected.items())
                 and av.get("contact_electric_field_fallback") is False,
                 "avalanche requires canonical SG/current-density/GradQF with electric fallback false")
    _require(av.get("driving_force_interpolation", "none") == "none",
             "avalanche driving-force interpolation is unsupported")
    for key in ("electron_driving_force_ref_density_m3", "hole_driving_force_ref_density_m3",
                "quasi_fermi_carrier_truncation", "quasi_fermi_carrier_trucation"):
        _require(_finite(av.get(key, 0), key) == 0, "avalanche density interpolation/truncation is unsupported")


def _components(document: dict, directory: Path, bias: float):
    mesh_path = (directory / document["mesh_file"]).resolve()
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = [node["id"] for node in mesh["nodes"]]
    _require(all(type(node) is int for node in nodes) and sorted(nodes) == list(range(len(nodes))),
             "mesh node IDs must be unique and contiguous")
    node_set = set(nodes)
    regions = {region["id"]: region["material"] for region in mesh["regions"]}
    _require(len(regions) == len(mesh["regions"]), "duplicate mesh regions")
    # Unknown materials are conservatively treated as transport materials.
    # Only known oxide/nitride names with no transport override are exempt.
    insulators = {"SiO2", "Si3N4", "Nitride"}
    if document.get("materials_file"):
        material_path = directory / document["materials_file"]
        materials = json.loads(material_path.read_text(encoding="utf-8"))["materials"]
        for material in materials:
            if material["name"] == "PolySilicon":
                _require(_finite(material.get("temperature_K", 300), "PolySilicon material temperature_K") == 300,
                         "PolySilicon material requires 300 K")
            # MaterialDatabase accepts both the legacy names and the explicit
            # centimeter-unit names in vela.templates_ldmos.materials.v1.
            # Zero/nonzero transport classification is independent of units.
            transport_keys = ("ni", "mun", "mup", "intrinsic_carrier_density_cm3",
                              "electron_mobility_cm2_per_V_s", "hole_mobility_cm2_per_V_s")
            if material["name"] in insulators and any(
                    _finite(material.get(key, 0), "material " + key) != 0 for key in transport_keys):
                insulators.remove(material["name"])
    adjacency = {}
    other_nodes = set()
    _require(not mesh.get("quadrilaterals") and not mesh.get("quads"), "requires Tri3 mesh")
    for cell in mesh["triangles"]:
        ids = cell["node_ids"]
        _require(len(ids) == 3 and len(set(ids)) == 3 and set(ids) <= node_set,
                 "malformed Tri3 node coverage")
        material = regions[cell["region_id"]]
        if material == "PolySilicon":
            for node in ids:
                adjacency.setdefault(node, set()).update(ids)
        elif material not in insulators:
            other_nodes.update(ids)
    _require(bool(adjacency), "mesh has no PolySilicon component")
    _require(not (set(adjacency) & other_nodes), "PolySilicon shares nodes with another transport material")
    contacts = document.get("contacts", [])
    configured = {contact["name"]: contact for contact in contacts}
    _require(len(configured) == len(contacts), "duplicate configured contacts")
    driven = document.get("sweep", {}).get("contact")
    _require(driven is None or driven in configured, "unknown driven sweep contact")
    electrodes = {}
    for contact in mesh["contacts"]:
        _require(set(contact["node_ids"]) <= node_set, "contact has unknown nodes")
        touched = set(contact["node_ids"]) & set(adjacency)
        if not touched:
            continue
        name = contact["name"]
        _require(name in configured, "PolySilicon contact lacks configuration: " + name)
        cfg = configured[name]
        _require(cfg.get("type") == "ohmic", "PolySilicon contacts must be ohmic: " + name)
        voltage = bias if name == driven else _finite(cfg.get("bias", 0), "contact bias")
        for node in touched:
            electrodes.setdefault(node, []).append((name, voltage))
    remaining = set(adjacency)
    targets = {}
    components = []
    while remaining:
        pending = [min(remaining)]
        component = set()
        while pending:
            node = pending.pop()
            if node not in component:
                component.add(node)
                pending.extend(adjacency[node] - component)
        remaining -= component
        attached = [entry for node in component for entry in electrodes.get(node, [])]
        voltages = {entry[1] for entry in attached}
        _require(len(voltages) == 1, "floating or conflicting PolySilicon component contacts")
        voltage = voltages.pop()
        targets.update(dict.fromkeys(component, voltage))
        components.append({"nodes": len(component), "first_node_id": min(component),
                           "qf_V": voltage, "contacts": sorted({entry[0] for entry in attached})})
    return targets, components, node_set, mesh_path


def write_poly_equilibrium_restart(document: dict, input_path: Path, output_path: Path,
                                   config_directory: Path, bias: float) -> dict:
    """Validate the narrow profile and atomically write an unqualified candidate.

    HDF5 state quantities are SI regardless of deck input scaling. At fixed psi,
    temperature, doping and quantum potential, Boltzmann n changes by
    exp((old_fn-new_fn)/Vt), and p by exp((new_fp-old_fp)/Vt).
    """
    source, output = Path(input_path).resolve(), Path(output_path).resolve()
    _require(source != output and not (output.exists() and source.samefile(output)),
             "input and output must be different files")
    bias = _finite(bias, "sweep bias")
    _profile(document)
    try:
        targets, components, nodes, mesh_path = _components(document, Path(config_directory), bias)
    except (KeyError, TypeError, AttributeError) as exc:
        raise ValueError("Poly equilibrium restart: malformed mesh/material/contact configuration") from exc
    source_bytes = source.read_bytes()
    arrays, metadata = states.read(source, document, config_directory)
    rows = states.archive.fields_to_rows(arrays)
    split = states.archive.SPLIT <= arrays.keys()
    seen = set()
    maximum = 0.0
    vt = 1.380649e-23 * 300 / 1.602176634e-19
    for row in rows:
        _require(None not in row and all(value is not None for value in row.values()), "malformed HDF5 row")
        try:
            node = int(row["node_id"])
        except ValueError as exc:
            raise ValueError("Poly equilibrium restart: invalid HDF5 node ID") from exc
        _require(node in nodes and node not in seen, "unknown or duplicate HDF5 node ID")
        seen.add(node)
        values = {key: _finite(value, key) for key, value in row.items() if key != "node_id"}
        for key in ("electrons_m3", "holes_m3"):
            _require(values[key] >= 0 and (node not in targets or values[key] > 0),
                     "carrier densities must be nonnegative and positive in PolySilicon")
        deltas = {}
        for carrier, field in (("electron", "phin"), ("hole", "phip")):
            if split:
                reference = values[f"{carrier}_qf_reference_V"]
                increment = values[f"{carrier}_qf_increment_V"]
                physical = math.fsum((reference, increment))
                tolerance = 32 * sys.float_info.epsilon * max(1, abs(physical), abs(values[field]))
                _require(abs(physical - values[field]) <= tolerance, "inconsistent physical and split QF state")
                delta = math.fsum((reference, -targets.get(node, 0), increment))
            else:
                delta = values[field] - targets.get(node, 0)
            deltas[carrier] = delta
        if node not in targets:
            continue
        target = targets[node]
        maximum = max(maximum, *(abs(delta) for delta in deltas.values()))
        for carrier, field, density, sign in (("electron", "phin", "electrons_m3", 1),
                                              ("hole", "phip", "holes_m3", -1)):
            try:
                transformed = values[density] * math.exp(sign * deltas[carrier] / vt)
            except OverflowError as exc:
                raise ValueError("Poly equilibrium restart: density transformation overflow") from exc
            _require(math.isfinite(transformed) and transformed > 0, "density transformation overflow/underflow")
            row[density] = transformed
            row[field] = target
            if split:
                row[f"{carrier}_qf_reference_V"] = target
                row[f"{carrier}_qf_increment_V"] = 0.
                # The new reference is the complete target QF. Retained packed
                # increments must describe that same zero increment, including
                # an optional low component. Preserve psi and non-Poly hints.
                packed = f"packed_{carrier}_qf_increment"
                for key in (packed, packed + "_low"):
                    if key in row:
                        row[key] = 0.
    _require(seen == nodes, "HDF5 state does not cover every mesh node")
    provenance = {"seed": str(source), "input_sha256": hashlib.sha256(source_bytes).hexdigest(),
                  "output": str(output), "mesh_file": str(mesh_path), "bias_V": bias,
                  "temperature_K": 300, "poly_nodes": len(targets), "components": components,
                  "max_qf_change_V": maximum, "qualified": False,
                  "scope": "isolated source-free constant-mobility PolySilicon candidate; original gates required"}
    states.write_candidate(output, states.archive.rows_to_fields(rows), metadata,
                           document, [source], bias, "poly_equilibrium_candidate")
    return provenance
