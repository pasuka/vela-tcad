"""Check a native SLOT IALMob request against a prepared Vela element deck."""

from decimal import Decimal
from copy import deepcopy
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
from tempfile import TemporaryDirectory

from sentaurus_parameter_ir import active_blocks, parameter_number_lexemes, parse_parameter_ir

from ..mesh_import import _load_dfise, _load_tdr
from ..materials import evaluate_semiconductor, project_material_arithmetic
from .support import RunError


def _require(condition, message):
    if not condition:
        raise RunError("E407", message)


def _blocks(par, name, material, section):
    entries = [block for block in active_blocks(par) if block["section"] == section
               and block["scope"] == {"kind": "material", "name": material}
               and "variant" not in block]
    _require(len(entries) == 1, f"{name} requires one base {material}.{section} block")
    return {entry["raw_name"]: parameter_number_lexemes(entry) for entry in entries[0]["parameters"]}


def _path(reference, value):
    return (reference.parent / value).resolve(strict=True)


def _topology(mesh):
    cells = {tuple(sorted(row["node_ids"])): mesh["regions"][row["region_id"]]["material"]
             for row in mesh["triangles"]}
    _require(len(cells) == len(mesh["triangles"]), "SLOT grid contains duplicate triangles")
    contacts = {row["name"]: set(row["node_ids"]) for row in mesh["contacts"]}
    return cells, contacts


def check_slot_ialmob(ir, cmd, grid, reference):
    # The generic emitter continues to reject unsupported physics. This
    # entry proves a frozen, explicitly prepared SLOT contract without solving.
    fatal = [entry for entry in ir["errors"] if entry["severity"] == "fatal"]
    _require(not ir["preprocessor_findings"] and len(fatal) == 1 and
             fatal[0]["id"] == "E106" and fatal[0]["text"] == "Enormal(IALMob)",
             "SLOT IALMob check expects the native Enormal(IALMob) command")
    physics = ir["physics"]
    models = [*physics["global"], *(model for scope in physics["scoped"] for model in scope["models"])]
    mobility = [m for m in models if m["model"].lower() == "mobility"]
    _require(len(mobility) == 1, "SLOT IALMob check requires one mobility block")
    args = {row["name"].lower(): row.get("args", []) for row in mobility[0]["args"]}
    _require(set(args) == {"enormal", "highfieldsaturation"} and len(mobility[0]["args"]) == 2 and
             len(args["enormal"]) == 1 and args["enormal"][0]["name"].lower() == "ialmob" and
             not args["enormal"][0].get("args") and not args["highfieldsaturation"],
             "SLOT IALMob check requires unoriented IALMob and default high-field drive")
    _require(physics["temperature_K"] == 300., "SLOT IALMob check requires 300 K")
    par_file = _path(cmd, ir["files"]["parameter"]["resolved"])
    par = parse_parameter_ir(par_file)
    low = _blocks(par, "IALMob check", "Silicon", "IALMob")
    high = _blocks(par, "IALMob check", "Silicon", "HighFieldDependence")
    _require(low["EnormMinimum"] and [Decimal(x) for x in low["EnormMinimum"]] == [0],
             "SLOT requires the supported zero EnormMinimum")
    _require([Decimal(x) for x in high["alpha"]] == [0, 0] and
             [Decimal(x) for x in high["Vsat_Formula"]] == [1, 1],
             "SLOT IALMob check requires high-field alpha=0 and velocity formula 1")
    masses = [low["me_over_m0"][0], low["mh_over_m0"][0]]
    low_fields = set(low) - {"me_over_m0", "mh_over_m0", "EnormMinimum"}
    _require(all(len(low[key]) == 2 for key in low_fields), "IALMob requires electron/hole parameter pairs")

    deck = json.loads(reference.read_text(encoding="utf-8"), parse_float=Decimal)
    _require(deck.get("format_version") == 2 or deck.get("scaling", {}).get("mode") == "unit_scaling",
             "SLOT IALMob reference must use TCAD input units")
    solver = deck["solver"]
    mob = solver["mobility"]
    ial = mob["ialmob"]
    _require(solver.get("state_variables") == "potential_density" and
             solver.get("coupled_operator") == "legacy" and mob["model"] == "ialmob" and
             mob["high_field_driving_force"] == "quasi_fermi_gradient" and ial["high_field"],
             "SLOT reference must use the validated legacy element IALMob and high-field path")

    coefficients = {}
    mapped = {"ialmob": {}}
    for index, carrier in enumerate(("electron", "hole")):
        expected = {key: low[key][index] for key in sorted(low_fields)}
        expected.update(mass=masses[index], other_mass=masses[1 - index])
        families = ial[f"{carrier}_parameters_cm"]
        _require(set(families) == {"100", "110", "111"}, f"Invalid {carrier} orientation families")
        for family, payload in families.items():
            _require(payload == {key: Decimal(val) for key, val in expected.items()},
                     f"{carrier} IALMob coefficients for family {family} differ from the base PAR block")
        mapped["ialmob"][f"{carrier}_parameters_cm"] = {
            family: {key: float(literal) for key, literal in expected.items()} for family in families
        }
        for target, source in (("field_beta", "beta0"), ("saturation_velocity_m_s", "vsat0")):
            literal = high[source][index]
            key = f"{carrier}_{target}"
            _require(mob[key] == Decimal(literal), f"{key} differs from PAR in TCAD deck units")
            expected[key] = literal
            mapped[key] = float(literal)
        coefficients[carrier] = expected

    mesh_file = _path(reference, deck["mesh_file"])
    mesh = json.loads(mesh_file.read_text(encoding="utf-8"))
    native = _load_dfise(grid)
    grid_file = _path(cmd, ir["files"]["grid"]["resolved"])
    reference_coords = {node["id"]: (float(node["x"]), float(node["y"])) for node in mesh["nodes"]}
    _require(reference_coords.keys() == native["coords"].keys(), "SLOT DF-ISE node IDs differ")
    if grid_file.suffix.lower() == ".tdr":
        with TemporaryDirectory(prefix="vela_slot_ialmob_") as tmp:
            authoritative = _load_tdr(grid_file, None, Path(tmp))
    else:
        _require(grid_file == grid, "The command's DF-ISE grid must be the checked grid")
        authoritative = native
    _require(reference_coords == authoritative["coords"], "SLOT reference and authoritative grid coordinates differ")
    differences = [abs(a - b) for node, point in reference_coords.items()
                   for a, b in zip(point, native["coords"][node]) if a != b]
    reference_cells, contacts = _topology(mesh)
    _require((reference_cells, contacts) == _topology(authoritative) == _topology(native),
             "SLOT reference, native grid and DF-ISE topology differ")
    _require(set(contacts) == {entry["name"] for entry in ir["electrodes"]} and
             ial["effective_electrodes"] and set(ial["effective_electrodes"]) <= contacts.keys(),
             "SLOT mesh and CMD electrode identity differ")
    for field in ("crystal_x", "crystal_y"):
        _require(len(ial[field]) == 3 and any(ial[field]), f"Missing SLOT crystal axis {field}")

    silicon_cells = {row["id"]: row["node_ids"] for row in mesh["triangles"]
                     if mesh["regions"][row["region_id"]]["material"] == "Si"}
    paths = {"cmd": cmd, "par": par_file, "native_grid": grid_file, "dfise_grid": grid,
             "dfise_doping": native["companion"],
             "deck": reference, "mesh": mesh_file, "doping": _path(reference, deck["node_doping_file"]),
             "material": _path(reference, deck["materials_file"])}
    for kind, schema in (("geometry_file", "vela.ialmob.transport_geometry.v1"),
                         ("avalanche_geometry_file", "vela.ialmob.observed_avalanche_geometry.v1")):
        path = _path(reference, ial[kind])
        paths[kind] = path
        geometry = json.loads(path.read_text(encoding="utf-8"))
        _require(geometry["schema"] == schema and geometry["node_count"] == len(mesh["nodes"]) and
                 geometry["cell_count"] == len(mesh["triangles"]) and
                 {row["cell_id"]: row["node_ids"] for row in geometry["cells"]} == silicon_cells and
                 len(geometry["cells"]) == len(silicon_cells), f"SLOT {kind} topology differs")
    inputs = {key: {"path": str(path), "sha256": sha256(path.read_bytes()).hexdigest()}
              for key, path in paths.items()}
    return {"schema": "vela.slot_ialmob_input_check.v1", "inputs": inputs, "coefficients": coefficients,
            "solver_mobility_overrides": mapped,
            "mesh_nodes": len(mesh["nodes"]), "mesh_cells": len(mesh["triangles"]), "si_cells": len(silicon_cells),
            "effective_electrodes": ial["effective_electrodes"], "crystal_x": ial["crystal_x"],
            "crystal_y": ial["crystal_y"], "reference_uses_tcad_units": True,
            "dfise_coordinate_differences": {"components": len(differences), "max_um": max(differences, default=0.)},
            "simulation_started": False}


def bind_slot_stage_inputs(ir, cmd, profile):
    reference = _path(cmd, profile["reference_deck"])
    report = check_slot_ialmob(ir, cmd, _path(cmd, profile["dfise_grid"]), reference)
    _require(report["inputs"]["par"]["sha256"] == profile["parameter_sha256"],
             "SLOT parameter file differs from the explicitly bound PAR")
    native = {model["key"] for model in ir["physics"]["effective"]}
    _require(native == {"Mobility.Enormal", "Mobility.HighFieldSaturation",
                        "EffectiveIntrinsicDensity.OldSlotboom", "Recombination.SRH",
                        "Recombination.Auger", "Recombination.Avalanche"},
             "SLOT stage binding requires the validated physical model set")
    reference_text = reference.read_text(encoding="utf-8")
    deck = json.loads(reference_text)
    solver = deck["solver"]
    _require(solver.get("carrier_statistics", {}).get("model") == "boltzmann" and
             solver.get("bandgap_narrowing", {}).get("model") == "old_slotboom" and
             set(solver.get("recombination", [])) == {"srh", "auger"} and
             solver.get("impact_ionization", {}).get("model") == "van_overstraeten",
             "SLOT reference physics differs from the validated BVDS configuration")
    srh = next(m["parameters"] for m in ir["physics"]["effective"] if m["key"] == "Recombination.SRH")
    avalanche = next(m["parameters"] for m in ir["physics"]["effective"] if m["key"] == "Recombination.Avalanche")
    _require(srh.get("doping_dependence") and avalanche.get("driving_force") == "GradQuasiFermi" and
             avalanche.get("model", "VanOverstraeten") == "VanOverstraeten",
             "SLOT binds SRH(DopingDep) and VanOverstraeten GradQuasiFermi only")
    _require(not ir["physics"]["fixed_charges"] and
             all(not (set(e["attributes"]) - {"resist_ohm_um"}) for e in ir["electrodes"]),
             "SLOT requires the validated Ohmic contacts and no additional fixed charges")
    mapped = report["solver_mobility_overrides"]
    solver["mobility"].update({key: val for key, val in mapped.items() if key != "ialmob"})
    solver["mobility"]["ialmob"].update(mapped["ialmob"])
    source_coefficients = {}
    if ir["math"].get("arithmetic_precision") == "extended80":
        exact_solver = json.loads(reference_text, parse_float=Decimal)["solver"]

        def retain(mapping, compatibility, keys, prefix=""):
            for key in keys:
                if key in mapping:
                    _require(float(mapping[key]) == compatibility[key],
                             f"SLOT bound source coefficient differs from its compatibility value: {prefix}{key}")
                    source_coefficients[prefix + key] = str(mapping[key])

        retain(exact_solver, solver, ("taun", "taup", "auger_cn_m6_per_s", "auger_cp_m6_per_s"))
        for group, keys in (
                ("srh_doping_dependence", ("tau_min_s", "tau_max_s", "reference_doping_m3", "gamma")),
                ("auger_density_dependence", ("enhancement", "reference_density_m3"))):
            for carrier in ("electron", "hole"):
                retain(exact_solver.get(group, {}).get(carrier, {}),
                       solver.get(group, {}).get(carrier, {}), keys, f"{group}.{carrier}.")
        exact_ial = exact_solver["mobility"]["ialmob"]
        bound_ial = solver["mobility"]["ialmob"]
        constant_key = "constant_material_mobility_cm2_per_Vs"
        for name, values in exact_ial.get(constant_key, {}).items():
            retain(values, bound_ial[constant_key][name], ("electron", "hole"),
                   f"ialmob.{constant_key}.{name}.")
        retain(exact_solver["bandgap_narrowing"], solver["bandgap_narrowing"],
               ("reference_doping_m3", "coefficient_eV", "smoothing", "offset_eV"), "bandgap_narrowing.")
        impact = exact_solver["impact_ionization"]
        retain(impact, solver["impact_ionization"],
               ("electron_a_low_m_inv", "electron_a_high_m_inv", "electron_b_low_V_m",
                "electron_b_high_V_m", "hole_a_low_m_inv", "hole_a_high_m_inv", "hole_b_low_V_m",
                "hole_b_high_V_m", "switch_field_V_m", "minimum_field_V_m"), "impact_ionization.")
        high_field = {}
        decimal = {"high_field": high_field}
        for carrier, coefficients in report["coefficients"].items():
            prefix = f"{carrier}_"
            high_field.update({key: coefficients[key] for key in (
                prefix + "field_beta", prefix + "saturation_velocity_m_s")})
            family_key = prefix + "parameters_cm"
            decimal[family_key] = {
                family: {key: coefficients[key] for key in fields}
                for family, fields in mapped["ialmob"][family_key].items()
            }
        solver["mobility"]["ialmob"]["decimal_coefficients"] = decimal
    for key in ("geometry_file", "avalanche_geometry_file"):
        solver["mobility"]["ialmob"][key] = report["inputs"][key]["path"]
    mesh_geometry = deck["mesh_geometry"]
    if "external_averagebox_couples_file" in mesh_geometry:
        path = _path(reference, mesh_geometry["external_averagebox_couples_file"])
        mesh_geometry["external_averagebox_couples_file"] = str(path)
        report["inputs"]["transport_couples"] = {"path": str(path), "sha256": sha256(path.read_bytes()).hexdigest()}
    for key in ("mesh", "material", "doping"):
        report["inputs"][key]["copy_as"] = {"mesh": "mesh.json", "material": "materials.json",
                                           "doping": "doping.csv"}[key]
    if source := profile.get("poly_material_parameters"):
        path = _path(cmd, source["file"])
        digest = sha256(path.read_bytes()).hexdigest()
        _require(digest == source["sha256"], "SLOT PolySi parameters differ from the explicitly bound PAR")
        report["inputs"]["poly_material_par"] = {"path": str(path), "sha256": digest,
                                                  "copy_as": "slot_poly_material.par"}
    seeds = {}
    for index, (prefix, value) in enumerate(sorted(profile.get("load_states", {}).items())):
        path = _path(cmd, value)
        _require(path.suffix.lower() == ".h5", "SLOT Load requires an explicit Vela HDF5 state")
        filename = f"slot_input_{index:02d}.h5"
        report["inputs"][f"load:{prefix}"] = {"path": str(path), "sha256": sha256(path.read_bytes()).hexdigest(),
                                             "copy_as": filename}
        seeds[prefix] = filename
    qualified = deepcopy(ir)
    qualified["errors"] = [e for e in ir["errors"] if e["id"] != "E106"]
    qualified["verdict"] = {"accepted": True, "exit_code": 0, "fatal_ids": []}
    model = next(m for m in qualified["physics"]["effective"] if m["key"] == "Mobility.Enormal")
    model.update(status="supported", vela="solver.mobility.model=ialmob")
    model.pop("finding_id", None)
    return qualified, report, {"reference": deck, "load_states": seeds,
                               "source_coefficients": source_coefficients}


def check_slot_effective_doping(resolved, bound):
    def load(path):
        with path.open(newline="", encoding="utf-8") as handle:
            return {int(row["node_id"]): (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
                    for row in csv.DictReader(handle)}

    effective, reference = load(resolved), load(bound)
    _require(effective.keys() == reference.keys(), "SLOT effective doping has a different node set from the bound deck")
    max_relative = 0.0
    for node, expected in reference.items():
        for value, target in zip(effective[node], expected):
            # The qualified bound deck may originate in a rounded TDX text export.
            _require(math.isclose(value, target, rel_tol=8e-15),
                     f"SLOT effective doping disagrees with the bound deck at node {node}")
            if value != target:
                max_relative = max(max_relative, abs(value - target) / max(abs(value), abs(target)))
    return {"matched_bound_nodes": len(reference), "max_relative_difference": max_relative}


def _decimal_material_binding(materials, name, parameter, temperature_literal, profile):
    entry = next((item for item in materials["materials"] if item["name"] == name), None)
    binding = {"status": "reference_only", "material": name}
    poly = name == "PolySilicon"
    try:
        parameters = project_material_arithmetic(parameter, name, fold_old_slotboom=True)
        sections = {
            name: {key: [float(value) for value in operand["literals"]]
                   for key, operand in values.items() if operand["literals"]}
            for name, values in parameters["sections"].items()
        }
        temperature = float(temperature_literal)
        values = evaluate_semiconductor(sections, temperature, True, constants_profile=profile)
        if poly:
            band = values["derivation"]
            values["electron_affinity_eV"] += sections["Bandgap"]["Bgn2Chi"][0] * (band["Eg_0_eV"] - band["Eg_T_eV"])
            values["bandgap_eV"] = band["effective_gap_eV"]
            parameters["poly_bound_convention"] = True
        fields = {"ni": "intrinsic_carrier_density_cm3", "Nc_m3": "conduction_band_density_of_states_cm3",
                  "Nv_m3": "valence_band_density_of_states_cm3", "bandgap_eV": "bandgap_eV",
                  "electron_affinity_eV": "electron_affinity_eV"}
        limit = 4 if poly else 0
        if entry is None or entry.get("temperature_K") != temperature:
            binding["reason"] = "Material or temperature differs from the bound deck"
            return binding
        errors = {target: abs(entry[target] - values[source]) / math.ulp(entry[target])
                  for target, source in fields.items()}
        if any(error > limit for error in errors.values()):
            binding.update(reason="Native PAR projection differs from the bound material", errors_in_ulps=errors)
            return binding
        entry["decimal_semiconductor"] = {**parameters, "temperature_K": temperature_literal}
        binding.update(status="verified_binary64_projection", fields=list(fields))
        if poly:
            binding.update(errors_in_ulps=errors, reference_potential_source="bound_material")
    except (KeyError, ValueError) as error:
        binding["reason"] = str(error)
    return binding


def stage_slot_decimal_material(path, report, reference, temperature_literal):
    """Bind decimal material inputs after checking their existing operating values."""
    text = path.read_text(encoding="utf-8")
    materials = json.loads(text)
    source = json.loads(text, parse_float=Decimal, parse_int=Decimal)
    exact = {entry["name"]: entry for entry in source["materials"]}
    retained = {}
    for entry in materials["materials"]:
        name = entry["name"]
        fields = []
        for key in ("eps_r", "intrinsic_reference_potential_V"):
            if key in exact[name]:
                value = exact[name][key]
                _require(isinstance(value, Decimal) and value.is_finite() and float(value) == entry[key],
                         f"SLOT material {name}.{key} has incompatible precision inputs")
                entry["decimal_" + key] = str(value)
                fields.append(key)
        if fields:
            retained[name] = fields
    profile = reference.get("physical_constants_profile", "codata2018")
    inputs = report["inputs"]
    bindings = {"decimal_material": _decimal_material_binding(
        materials, "Si", Path(inputs["par"]["path"]), temperature_literal, profile)}
    if "poly_material_par" in inputs:
        bindings["decimal_poly_material"] = _decimal_material_binding(
            materials, "PolySilicon", path.with_name(inputs["poly_material_par"]["copy_as"]),
            temperature_literal, profile)
    if retained or any(entry["status"] == "verified_binary64_projection" for entry in bindings.values()):
        path.write_text(json.dumps(materials, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        digest = sha256(path.read_bytes()).hexdigest()
        for entry in bindings.values():
            if entry["status"] == "verified_binary64_projection":
                entry["staged_sha256"] = digest
    report.update(bindings, decimal_bound_material_fields=retained)
