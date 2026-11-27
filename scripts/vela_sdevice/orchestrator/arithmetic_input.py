"""Bind source numeric text for the C++ arithmetic input check."""

import csv
from hashlib import sha256
import json

from compare_sentaurus_tdr_tdx import iter_dat_datasets, parse_grd
from sentaurus_parameter_ir import active_blocks, parameter_number_lexemes, parse_parameter_ir

from ..mesh_import import ACCEPTOR_DATASETS, DONOR_DATASETS
from .support import RunError


def _source(path, numbers):
    return {"path": str(path), "sha256": sha256(path.read_bytes()).hexdigest(), "numbers": numbers}


class _Number(str):
    """Keep JSON numeric tokens distinct from string values."""


def _json_numbers(value, pointer=""):
    if isinstance(value, _Number):
        yield pointer, str(value)
    elif isinstance(value, dict):
        for key, child in value.items():
            escaped = key.replace("~", "~0").replace("/", "~1")
            yield from _json_numbers(child, f"{pointer}/{escaped}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            yield from _json_numbers(child, f"{pointer}/{index}")


def _read_json_tokens(path):
    return json.loads(path.read_text(encoding="utf-8"), parse_float=_Number, parse_int=_Number)


def _decimal_strings(value):
    if isinstance(value, dict):
        return {key: _decimal_strings(child) for key, child in value.items()}
    if isinstance(value, list):
        return [_decimal_strings(child) for child in value]
    return _Number(value)


def _csv_numbers(path, id_field, columns):
    with path.open(newline="", encoding="utf-8") as handle:
        return {f"{int(row[id_field])}/{column}": row[column]
                for row in csv.DictReader(handle) for column in columns}


def collect_input(cmd, cmd_input, par, grid, sidecar, profile, *, tdr_export=None, doping_fields=()):
    par_input = {}
    ir = parse_parameter_ir(par)
    for block in active_blocks(ir):
        for entry in block["parameters"]:
            if entry.get("shadowed_by") is not None:
                continue
            source = entry["source"]
            key = f"{source['file']}:{source['line']}:{block['section']}.{entry['raw_name']}"
            par_input.update((f"{key}/{index}", val)
                             for index, val in enumerate(parameter_number_lexemes(entry)))
    sources = {"cmd": _source(cmd, cmd_input["numbers"]), "par": _source(par, par_input)}
    if sidecar.is_file():
        sources["sidecar"] = _source(sidecar, dict(_json_numbers(_read_json_tokens(sidecar))))
    if grid.suffix.lower() == ".tdr":
        sources.update(collect_tdr_export(grid, tdr_export))
    elif grid.suffix.lower() == ".grd":
        sources.update(_collect_dfise(grid, doping_fields))
    else:
        raise RunError("E402", "Arithmetic input check requires DF-ISE or an imported HDF5 TDR")
    return {"schema": "vela.arithmetic_input.v1", "physical_constants_profile": profile,
            "temperature_K": cmd_input["temperature_K"], "sources": sources}


def collect_slot_input(work, binding, reference_deck, stage_decks):
    """Audit the files selected by SLOT binding, including the staged replacements."""
    mesh = work / "mesh.json"
    coordinates = {}
    for node in _read_json_tokens(mesh)["nodes"]:
        for component in ("x", "y"):
            value = node[component]
            if not isinstance(value, _Number):
                raise RunError("E402", f"SLOT bound coordinate {node['id']}/{component} is not numeric")
            coordinates[f"{int(node['id'])}/{component}"] = str(value)
    doping = work / "doping.csv"
    concentrations = _csv_numbers(doping, "node_id", ("donors_cm3", "acceptors_cm3"))
    if not coordinates or not concentrations:
        raise RunError("E402", "SLOT bound mesh or doping has no values to audit")

    materials_file = work / "materials.json"
    materials = {row["name"]: row["decimal_semiconductor"]
                 for row in json.loads(materials_file.read_text(encoding="utf-8"))["materials"]
                 if "decimal_semiconductor" in row}
    sources = {
        "slot_reference_deck": _source(reference_deck, dict(_json_numbers(_read_json_tokens(reference_deck)))),
        "slot_bound_coordinates": _source(mesh, coordinates),
        "slot_bound_doping": _source(doping, concentrations),
        "slot_staged_materials": _source(materials_file, dict(_json_numbers(_read_json_tokens(materials_file)))),
    }
    for key, source in (("deck", "slot_reference_deck"), ("mesh", "slot_bound_coordinates"),
                        ("doping", "slot_bound_doping")):
        if binding["inputs"][key]["sha256"] != sources[source]["sha256"]:
            raise RunError("E409", f"SLOT {key} differs from the bound input during arithmetic audit")
    for key in ("decimal_material", "decimal_poly_material"):
        if "staged_sha256" in binding.get(key, {}):
            if binding[key]["staged_sha256"] != sources["slot_staged_materials"]["sha256"]:
                raise RunError("E409", "SLOT staged material differs from the bound input during arithmetic audit")

    for stage_deck in stage_decks:
        config = json.loads(stage_deck.read_text(encoding="utf-8"))
        fields = dict(config.get("decimal_stage_controls", {}))
        ialmob = config["solver"].get("mobility", {}).get("ialmob", {}) if "solver" in config else {}
        if "decimal_coefficients" in ialmob:
            fields["ialmob"] = ialmob["decimal_coefficients"]
        if fields:
            sources[stage_deck.stem] = _source(stage_deck, dict(_json_numbers(_decimal_strings(fields))))
    return sources, materials


def collect_effective_doping(path):
    return _source(path, _csv_numbers(path, "node_id", ("donors_cm3", "acceptors_cm3")))


def _collect_dfise(grid, doping_fields):
    coordinates = parse_grd(grid, coordinate=str)["vertices"]
    if not coordinates:
        raise RunError("E402", "Arithmetic input check found no DF-ISE coordinates")
    position = _source(grid, {f"{node}/{component}": value
                              for node, row in enumerate(coordinates)
                              for component, value in enumerate(row)})
    dat = grid.with_suffix(".dat")
    names = set(DONOR_DATASETS) | set(ACCEPTOR_DATASETS) | {"DopingConcentration", "TotalConcentration"} | set(doping_fields)
    doping = {f"{name}/{block}/{index}": value
              for block, (name, record) in enumerate(iter_dat_datasets(dat, number=str)) if name in names
              for index, value in enumerate(record["values"])}
    if not doping:
        raise RunError("E402", "Arithmetic input check found no DF-ISE doping")
    return {"grd": position, "dat": _source(dat, doping)}


def collect_tdr_export(grid, export):
    """The native importer writes enough digits to recover its binary64 values."""
    if export is None:
        raise RunError("E402", "Arithmetic input check requires sentaurus_import output")
    original_sha = sha256(grid.read_bytes()).hexdigest()
    sources = {}
    for role, filename, id_field, columns in (
        ("tdr_coordinates", "nodes.csv", "id", ("x_um", "y_um")),
        ("tdr_doping", "doping.csv", "node_id", ("donors_cm3", "acceptors_cm3")),
    ):
        path = export / filename
        if not path.is_file():
            raise RunError("E402", f"TDR arithmetic input requires native {filename}")
        numbers = _csv_numbers(path, id_field, columns)
        if not numbers:
            raise RunError("E402", f"TDR arithmetic input has no values in {filename}")
        sources[role] = {**_source(path, numbers), "number_encoding": "binary64_roundtrip",
                         "original_tdr_sha256": original_sha}
    return sources
