"""WP4: pure execution-IR to v2 decks and C4 stage plans.

The emitter neither reads command text nor runs a solver. Input paths in the
runner decks are ASCII basenames staged by the WP5 orchestrator.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
import json
import math
from pathlib import Path
import re

from .. import schema_validator as sv
from ..resources import resource_root

REPO = resource_root()
SCHEMA_ID = "vela.stage_plan.v1"
IMPLICIT_FIELDS = ("ElectrostaticPotential", "eDensity", "hDensity", "DopingConcentration")
AVALANCHE_PROFILES = {
    "element_edge_sg_gss_laux": {
        "impact_ionization": {"current_approximation": "element_edge_sg_gss_laux",
            "source_mapping_mode": "element_vertex_box_measure", "cell_reconstructed_midpoint_density": "bernoulli"},
        "mesh_geometry": {"node_volume_policy": "mixed_voronoi", "require_non_obtuse": True}},
    "legacy_cell_reconstructed": {
        "impact_ionization": {"current_approximation": "cell_reconstructed",
            "source_mapping_mode": "triangle_gss_gradqf_truncated", "cell_reconstructed_midpoint_density": "gss_logistic"},
        "mesh_geometry": {"node_volume_policy": "barycentric", "require_non_obtuse": False}},
}
# src/io/VtuWriter.cpp owns the C++ catalogue (the header declares its API).
TDR_FIELDS = (
    "ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential",
    "eDensity", "hDensity", "DopingConcentration", "DonorConcentration",
    "AcceptorConcentration", "SpaceCharge", "ElectricField", "eCurrentDensity",
    "hCurrentDensity", "TotalCurrentDensity", "ConductionCurrentDensity",
    "eMobility", "hMobility", "eVelocity", "hVelocity", "srhRecombination",
    "AugerRecombination", "Band2BandGeneration", "AvalancheGeneration",
    "TotalRecombination", "eAlphaAvalanche", "hAlphaAvalanche", "eIonIntegral",
    "hIonIntegral", "MeanIonIntegral", "BandGap", "BandgapNarrowing",
    "EffectiveBandGap", "ElectronAffinity", "ConductionBandEnergy", "ValenceBandEnergy",
    "EffectiveIntrinsicDensity", "eEffectiveStateDensity", "hEffectiveStateDensity",
    "eEparallel", "hEparallel", "eEnormal", "hEnormal", "eGradQuasiFermi", "hGradQuasiFermi",
)
_SUFFIXES = (
    ("_m2_V_s", "_cm2_V_s"), ("_m6_per_s", "_cm6_per_s"),
    ("_m_per_s", "_cm_per_s"), ("_m_per_V", "_cm_per_V"),
    ("_V_per_m", "_V_per_cm"), ("_m_inv", "_cm_inv"),
    ("_m_s", "_cm_s"), ("_V_m", "_V_per_cm"), ("_m3", "_cm3"), ("_m2", "_cm2"),
)
_OLD_KEYS = {"ni": "ni_cm3", "mun": "mun_cm2_V_s", "mup": "mup_cm2_V_s",
             "donors": "donors_cm3", "acceptors": "acceptors_cm3",
             "contact_radius": "contact_radius_um", "depth_m": "depth_um",
             "surface_recombination_velocity": "surface_recombination_velocity_cm_per_s"}
_MESSAGES = {
    "E401": "Electrode '{name}' has no matching contact in the mesh (contacts: {list})",
    "E406": "Generated deck '{deck}' fails schema validation: {detail}",
    "E407": "Stage chain is inconsistent: {detail}",
    "E408": "Node-volume policy '{policy}' requires a non-obtuse mesh ({n} obtuse elements)",
    "E409": "Sidecar '{file}' is invalid: {detail}",
}


class DeckEmitError(ValueError):
    def __init__(self, message_id: str, **fields: object) -> None:
        self.message_id = message_id
        self.fields = fields
        self.message = _MESSAGES.get(message_id, "{detail}").format(
            **{"detail": message_id, **fields})
        super().__init__(f"{message_id}: {self.message}")


@dataclass(frozen=True)
class EmitOptions:
    exit_on_failure: bool = False
    threads: int | None = None
    plot_final_only: bool = False
    vtu_binary: bool = False
    keep_intermediate: bool = True
    allow_lossy: bool = False
    defer_current_times: bool = False


@dataclass
class StagePlanResult:
    stage_plan: dict
    decks: list[dict]
    decisions: list[dict]
    # Regeneration preserves save/load snapshots, controls and sidecar decisions
    # when an earlier stage stops before its goal. These are not part of C4.
    _context: dict = field(default_factory=dict, repr=False)


def v2_key(path: str) -> str:
    if path.startswith("solver.band_to_band."):
        return path
    parts = path.split(".")
    key = parts[-1]
    if key in _OLD_KEYS:
        key = _OLD_KEYS[key]
    else:
        for old, new in _SUFFIXES:
            if key.endswith(old):
                key = key[:-len(old)] + new
                break
    return ".".join(parts[:-1] + [key])


def _set(doc: dict, path: str, value: object) -> None:
    keys = path.split(".")
    for key in keys[:-1]:
        if key in doc and not isinstance(doc[key], dict):
            raise DeckEmitError("E407", detail=f"conflicting parameter path {path}")
        doc = doc.setdefault(key, {})
    doc[keys[-1]] = deepcopy(value)


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise DeckEmitError("E407", detail=f"{label} must be a finite number")
    return float(value)


def _basename(value: str, label: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9._-]+", value) or value in (".", ".."):
        raise DeckEmitError("E407", detail=f"{label} must be an ASCII basename")
    return value


def expand_times(spec: dict | None) -> list[float]:
    """Expand a C1 Time list/range with exact declared endpoints."""
    if spec is None:
        return []
    if "list" in spec:
        values = spec["list"]
    else:
        a, b = spec["range"]
        count = spec["intervals"]
        if isinstance(count, bool) or not isinstance(count, int) or count < 1 or a > b:
            raise DeckEmitError("E407", detail="invalid Time range or interval count")
        values = [a + (b - a) * i / count for i in range(count + 1)]
        values[0], values[-1] = a, b
    result = sorted(set(_finite(v, "Time") for v in values))
    if not result or result[0] < 0 or result[-1] > 1:
        raise DeckEmitError("E407", detail="Time points must lie in [0,1]")
    return result


def _bias(v0: float, v1: float, t: float) -> float:
    return v0 if t == 0 else v1 if t == 1 else v0 + t * (v1 - v0)


def _decision(decisions: list, keyword: str, value: object, key: str | None,
              mapped: object, note: str, message_id: str | None = None,
              stage: int | None = None) -> None:
    entry = {"keyword": keyword, "value": deepcopy(value), "vela_key": key,
             "vela_value": deepcopy(mapped), "note": note}
    if message_id:
        entry["message_id"] = message_id
    if stage is not None:
        entry["stage"] = stage
    decisions.append(entry)


def _reference_phi(report: dict) -> float:
    # All inputs come from WP3's evaluated report; no C++ built-in values.
    row = next((r for r in report["materials"] if r["name"] == "Si"), None)
    if row is None or not row.get("derivation"):
        raise DeckEmitError("E407", detail="WorkFunction requires the WP3 Si derivation")
    derivation = row["derivation"]
    chi = next((p["values"][0] for p in row["parameters"]
                if p["section"] == "Bandgap" and p["parameter"] == "Chi0"), None)
    if chi is None:
        raise DeckEmitError("E407", detail="WP3 report is missing Silicon Bandgap.Chi0")
    nc_nv = (derivation["electron_dos_mass_ratio"] / derivation["hole_dos_mass_ratio"]) ** 1.5
    return float(chi) + derivation["Eg_T_eV"] / 2 + derivation["kT_eV"] / 2 * math.log(nc_nv)


def _physics(ir: dict, report: dict, decisions: list) -> dict:
    deck = {"solver": {"mobility": {"model": "constant", "doping_concentration_basis": "total_impurity"}, "recombination": [],
                       "impact_ionization": {"model": "none"}}}
    solver = deck["solver"]
    # Ordinary SDevice Ohmic electrodes pin both carrier quasi-Fermi levels
    # to the applied voltage, including the minority carrier at high bias.
    solver["contact_boundary_minority_electron_relaxation"] = False
    suffixes = set()
    for model in ir["physics"]["effective"]:
        if model["status"] == "unsupported":
            raise DeckEmitError(model.get("finding_id", "E106"), detail=model["key"])
        mapping = model.get("vela")
        if model["status"] != "supported" or not mapping:
            _decision(decisions, model["key"], model["parameters"], None, None, "Metadata only")
            continue
        for operation in mapping.split(";"):
            operation = operation.strip()
            match = re.fullmatch(r"(solver\.[A-Za-z0-9_.]+)(\+=|=)([^=]+)", operation)
            if not match:
                raise DeckEmitError("E407", detail=f"unsupported Physics mapping {operation}")
            path, op, raw = match.groups()
            value = {"true": True, "false": False}.get(raw, raw)
            if path == "solver.mobility.model" and op == "+=":
                suffixes.add(raw)
            elif path == "solver.recombination" and op == "+=":
                if raw not in solver["recombination"]:
                    solver["recombination"].append(raw)
            elif path == "solver.bandgap_narrowing":
                _set(deck, path + ".model", value)
            elif path == "solver.srh_doping_dependence":
                _set(deck, path + ".enabled", value)
            else:
                _set(deck, path, value)
    solver["mobility"]["model"] += "".join(s for s in ("_field", "_lombardi") if s in suffixes)
    if "_field" in suffixes:
        solver["mobility"].setdefault("high_field_driving_force", "quasi_fermi_gradient")
    if solver["impact_ionization"]["model"] != "none":
        avalanche = solver["impact_ionization"]
        avalanche.setdefault("driving_force", "quasi_fermi_gradient")
        avalanche.update(generation="current_density", current_approximation="density_gradient")
        _decision(decisions, "Avalanche", "SDevice drift-diffusion defaults",
                  "solver.impact_ionization", avalanche,
                  "Current-density source with explicit driving force; generic SG edge support")
    values = {}
    for row in report["solver_parameters"]:
        if row.get("status") == "approximated" and not report.get("allow_lossy"):
            raise DeckEmitError("E407", detail="lossy materials report is not authorized")
        path = v2_key(row["target"])
        value = row["value_si"] if path.startswith("solver.band_to_band.") else row["value_par"]
        if row["target"] in ("solver.auger_cn_m6_per_s", "solver.auger_cp_m6_per_s"):
            # Older reports include the raw A/B/C rows as well as the combined row.
            if row["parameter"] != "A+B(T/300)+C(T/300)^2":
                continue
            value = row["value_si"] * 1e12
        _finite(value, path)
        if path in values and values[path] != value:
            raise DeckEmitError("E407", detail=f"material-dependent {path} cannot use a global solver value")
        values[path] = value
        _set(deck, path, value)
        if row["section"] == "Scharfetter" and row["parameter"] == "taumax":
            _set(deck, "solver.tau_n_s" if row["carrier"] == "electron" else "solver.tau_p_s", value)
    # A parameter declaration does not activate a Physics option. The C++
    # schema otherwise enables doping dependence whenever this object exists.
    if "srh_doping_dependence" in solver:
        solver["srh_doping_dependence"].setdefault("enabled", False)
    solver["temperature_K"] = ir["physics"]["temperature_K"] or 300.0
    if report.get("temperature_K") != solver["temperature_K"]:
        raise DeckEmitError("E407", detail="IR and materials temperatures differ")
    # OldSlotboom dEg0 is already in the material ni; retain only the doping law here.
    return solver


def _contacts(ir: dict, mesh: dict, report: dict, decisions: list) -> tuple[list, list]:
    names = [c["name"] for c in mesh["contacts"]]
    contacts, electrodes = [], []
    for electrode in ir["electrodes"]:
        name, attr = electrode["name"], electrode["attributes"]
        if name not in names:
            raise DeckEmitError("E401", name=name, list=", ".join(names))
        contact = {"name": name, "bias": _finite(electrode["voltage_V"], name), "type": "ohmic"}
        planned = {"name": name, "contact": name, "initial_voltage_V": contact["bias"]}
        schottky = attr.get("schottky", False)
        if schottky:
            contact.update(type="schottky", emission_model="thermionic_robin")
            if "barrier_V" in attr:
                contact["barrier_eV"] = attr["barrier_V"]
            for carrier, short, default in (("electron", "e", 2.573e6), ("hole", "h", 1.93e6)):
                key = short + "_rec_velocity_cm_s"
                velocity = attr.get(key, default)
                if _finite(velocity, key) <= 0:
                    raise DeckEmitError("E407", detail=f"{name}: {key} must be positive")
                contact[carrier + "_surface_recombination_velocity_cm_per_s"] = velocity
            if "e_rec_velocity_cm_s" not in attr or "h_rec_velocity_cm_s" not in attr:
                _decision(decisions, "Schottky", name, "contacts", contact,
                          "Manual Table 206 default recombination velocities", "W209")
        elif any(k in attr for k in ("barrier_V", "workfunction_eV", "material")):
            contact["type"] = "metal_gate"
        if not schottky and "barrier_V" in attr:
            contact["work_function_eV"] = attr["barrier_V"]
        if "workfunction_eV" in attr or "material" in attr:
            if "workfunction_eV" in attr:
                workfunction = attr["workfunction_eV"]
            else:
                row = next((r for r in report["electrodes"] if r["material"] == attr["material"]), None)
                if row is None:
                    raise DeckEmitError("E407", detail=f"missing electrode material {attr['material']}")
                workfunction = row["work_function_eV"]
            if "barrier_V" in attr:
                raise DeckEmitError("E407", detail=f"{name}: Barrier and WorkFunction are mutually exclusive")
            contact["work_function_eV"] = workfunction if schottky else (
                workfunction - _reference_phi(report))
        if schottky and not any(k in contact for k in ("barrier_eV", "work_function_eV")):
            raise DeckEmitError("E407", detail=f"Schottky contact {name} has no barrier or WorkFunction")
        if "resist_ohm_um" in attr:
            resistance = _finite(attr["resist_ohm_um"], "Resist")
            if resistance < 0:
                raise DeckEmitError("E407", detail="Resist must be non-negative")
            planned["series_resistance_ohm_um"] = resistance
        contacts.append(contact)
        electrodes.append(planned)
    return contacts, electrodes


def _charges(ir: dict, mesh: dict) -> dict:
    result = {}
    for charge in ir["physics"]["fixed_charges"]:
        scope, target = charge["scope"], charge["target"]
        if "conc_cm3" in charge:
            regions = [r["name"] for r in mesh["regions"] if
                       (scope == "Region" and r["name"] == target) or
                       (scope == "Material" and target.lower() in
                        (r["material"].lower(), r.get("source_material", "").lower()))]
            if not regions:
                raise DeckEmitError("E407", detail=f"fixed-charge scope {scope}={target} has no region")
            result.setdefault("regions", []).extend(
                {"name": name, "fixed_charge_cm3": charge["conc_cm3"]} for name in regions)
        else:
            pair = target.split("/")
            if len(pair) != 2:
                raise DeckEmitError("E407", detail="interface fixed charge requires an A/B pair")
            if scope == "RegionInterface":
                if pair[0] == pair[1] or any(name not in {r["name"] for r in mesh["regions"]} for name in pair):
                    raise DeckEmitError("E407", detail=f"fixed-charge interface {target} has no region pair")
                pairs = [pair]
            else:
                aliases = {"silicon": "Si", "oxide": "SiO2", "nitride": "Si3N4"}
                pair = [aliases.get(v.lower(), v) for v in pair]
                left = [r["name"] for r in mesh["regions"] if r["material"] == pair[0]]
                right = [r["name"] for r in mesh["regions"] if r["material"] == pair[1]]
                pairs = [[a, b] for a in left for b in right if a != b]
                if not pairs:
                    raise DeckEmitError("E407", detail=f"fixed-charge interface {target} has no region pair")
            result.setdefault("interfaces", []).extend(
                {"regions": regions, "fixed_charge_cm2": charge["conc_cm2"]} for regions in pairs)
    return result


def _merge(dst: dict, src: dict, prefix: str, overridden: list, decisions: list, index: int) -> None:
    for key, value in src.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            if not isinstance(dst.get(key), dict):
                dst[key] = {}
            _merge(dst[key], value, path, overridden, decisions, index)
        else:
            dst[key] = deepcopy(value)
            overridden.append(path)
            _decision(decisions, "sidecar", value, path, value, "Sidecar override", "W207", index)


def _validate_result(result: StagePlanResult) -> None:
    try:
        if len(result.stage_plan["stages"]) != len(result.decks):
            raise ValueError("stage/deck counts differ")
        sv.validate(sv.load_schema(REPO / "schemas/vela.stage_plan.v1.schema.json"), result.stage_plan)
        schema = sv.load_schema(REPO / "configs/schema/vela-simulation.schema.json")
        for stage, deck in zip(result.stage_plan["stages"], result.decks):
            sv.validate(schema, deck)
            def check(node: object, prefix: str = "") -> None:
                if isinstance(node, dict):
                    for key, value in node.items():
                        path = f"{prefix}.{key}" if prefix else key
                        if key == "scaling" or v2_key(path) != path:
                            raise ValueError(f"legacy v1 key {path}")
                        check(value, path)
                elif isinstance(node, list):
                    for value in node:
                        check(value, prefix)
                elif isinstance(node, float) and not math.isfinite(node):
                    raise ValueError(f"non-finite value at {prefix}")
            check(deck)
    except (sv.SchemaError, sv.ValidationError, ValueError) as error:
        raise DeckEmitError("E406", deck="stage_plan/decks", detail=str(error)) from error


def emit_stage_plan(ir: dict, *, mesh_qualification: dict, materials_report: dict,
                    work_dir: str, cmd_basename: str, cmd_sha256: str,
                    options: EmitOptions | None = None, sidecar: dict | None = None) -> StagePlanResult:
    """Generate a deterministic C4 plan and runner decks without mutating inputs."""
    try:
        return _emit(ir, mesh_qualification, materials_report, work_dir, cmd_basename,
                     cmd_sha256, options or EmitOptions(), sidecar, {})
    except DeckEmitError:
        raise
    except (KeyError, TypeError, IndexError, ValueError, OverflowError) as error:
        raise DeckEmitError("E407", detail=f"malformed emitter input: {error}") from error


def _emit(ir: dict, mesh: dict, report: dict, work_dir: str, cmd_name: str, sha: str,
          options: EmitOptions, sidecar: dict | None, rebases: dict,
          current_starts: dict | None = None, saved_voltages: dict | None = None) -> StagePlanResult:
    current_starts = current_starts or {}
    saved_voltages = saved_voltages or {}
    if not ir["verdict"]["accepted"]:
        raise DeckEmitError(ir["verdict"]["fatal_ids"][0], detail="rejected execution IR")
    try:
        sv.validate(sv.load_schema(REPO / "schemas/vela.sentaurus_execution_ir.v2.schema.json"), ir)
        sv.validate(sv.load_schema(REPO / "schemas/vela.mesh_qualification.v1.schema.json"), mesh)
    except (sv.SchemaError, sv.ValidationError) as error:
        raise DeckEmitError("E407", detail=str(error)) from error
    if sha != ir["source"]["sha256"]:
        raise DeckEmitError("E407", detail="cmd SHA256 differs from execution IR")
    if report.get("schema") != "vela.sdevice_materials_report.v1" or report.get("blocking"):
        raise DeckEmitError("E407", detail="invalid or blocking materials report")
    if report.get("lossy") and not options.allow_lossy:
        raise DeckEmitError("E407", detail="lossy materials require allow_lossy")
    _basename(work_dir, "work_dir")
    if mesh["dimension"] != 2 or not all(mesh["qualification"][k]
            for k in ("two_dimensional", "triangles_only", "coordinates_plausible")):
        raise DeckEmitError("E407", detail="a 2-D triangle mesh is required")
    if sidecar is not None:
        try:
            sv.validate(sv.load_schema(REPO / "schemas/vela.sdevice_sidecar.v1.schema.json"), sidecar)
        except (sv.SchemaError, sv.ValidationError) as error:
            raise DeckEmitError("E409", file="sidecar", detail=str(error)) from error
        protected = {"tau_n_s", "tau_p_s", "auger_cn_cm6_per_s", "auger_cp_cm6_per_s",
                     "srh_doping_dependence", "hole_quantum_potential"}
        if protected.intersection(sidecar.get("solver", {})):
            raise DeckEmitError("E409", file="sidecar", detail="physical parameters belong to the cmd/par")
    decisions = []
    solver = _physics(ir, report, decisions)
    contacts, electrodes = _contacts(ir, mesh, report, decisions)
    if mesh["angles"]["obtuse_element_count"]:
        _decision(decisions, "mesh.angles", mesh["angles"], "mesh_geometry.node_volume_policy",
                  "barycentric", "Imported obtuse elements retained with barycentric volumes")
    fields = [f["data_name"] + ("/Vector" if f["component"] == "vector" else "")
              for f in ir["plot"]["fields"]]
    if not ir["plot"]["explicit"]:
        added = [name for name in IMPLICIT_FIELDS if name not in fields]
        fields.extend(added)
        _decision(decisions, "PlotExplicit", False, "sweep.vtk_fields", added,
                  "Implicit plot fields appended", "W203")
    fields = list(dict.fromkeys(fields))
    math_cfg = ir["math"]
    for keyword in ("digits", "rel_err_control", "err_ref_electron_cm3", "err_ref_hole_cm3"):
        if math_cfg.get(keyword) is not None:
            _decision(decisions, keyword, math_cfg[keyword], None, None,
                      "Recorded without loosening Vela convergence tolerances", "W201")
    for ignored in math_cfg["ignored"]:
        _decision(decisions, ignored["keyword"], ignored["value"], None, None,
                  "Recorded and ignored", ignored["finding_id"])
    threads = options.threads if options.threads is not None else math_cfg["number_of_threads"]
    plan = {"schema": SCHEMA_ID, "cmd": {"path": cmd_name, "sha256": sha, "execution_ir": "execution_ir.json"},
            "work_dir": work_dir, "exit_on_failure": options.exit_on_failure or math_cfg["exit_on_failure"],
            "threads": threads, "inputs": {"mesh_json": "mesh.json", "doping_csv": "doping.csv",
                "materials_json": "materials.json", "mesh_qualification": "mesh_qualification.json",
                "sidecar": None if sidecar is None else Path(cmd_name).stem + ".vela.json"},
            "outputs": {"log": ir["files"]["output"]["resolved"],
                "final_vtu": ir["files"]["plot"]["vela"],
                "pvd": str(Path(ir["files"]["plot"]["vela"]).with_suffix(".pvd")),
                "plot_final_only": options.plot_final_only, "vtu_binary": options.vtu_binary,
                "keep_intermediate": options.keep_intermediate},
            "electrodes": electrodes, "stages": []}
    base = {"format_version": 2, "simulation_type": "dc_sweep", "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv", "materials_file": "materials.json", "doping": [],
            "mesh_geometry": {"node_volume_policy": "barycentric", "require_non_obtuse": False},
            "solver": solver, **_charges(ir, mesh)}
    voltages = {e["name"]: e["initial_voltage_V"] for e in electrodes}
    sources = deepcopy(ir["stages"])
    groups = []
    for source in sources:
        if source["kind"] == "initial" and groups and groups[-1][-1]["kind"] == "initial":
            groups[-1].append(source)
        else:
            groups.append([source])
    decks, snapshots = [], {}
    global_plot = 0
    for group in groups:
        source = group[-1]
        kind = source["kind"]
        if kind == "save":
            if not plan["stages"]:
                raise DeckEmitError("E407", detail="Save precedes a solved stage")
            statement = ir["solve"]["statements"][source["statement_index"]]
            prefix = statement["file_prefix"]
            state = _basename(prefix + ("" if "." in prefix else "_des.sav"), "Save")
            state = str(Path(state).with_suffix(".h5"))
            previous = plan["stages"][-1]
            if previous["output_state"] in [v[0] for v in snapshots.values()]:
                raise DeckEmitError("E407", detail="multiple Save aliases of one stage cannot be emitted")
            previous["output_state"] = state
            decks[-1]["sweep"]["write_state_file"] = state
            snapshots[prefix] = (state, deepcopy(saved_voltages.get(state, voltages)))
            _decision(decisions, "Save", prefix, "sweep.write_state_file", state, "HDF5 snapshot", "W205")
            continue
        index = len(decks)
        if index in rebases and kind != "load":
            if set(rebases[index]) - set(voltages):
                raise DeckEmitError("E407", detail="rebase contains an unknown electrode")
            voltages.update({k: _finite(v, k) for k, v in rebases[index].items()})
        input_state = plan["stages"][-1]["output_state"] if decks else None
        if kind == "load":
            prefix = ir["solve"]["statements"][source["statement_index"]]["file_prefix"]
            if prefix not in snapshots:
                raise DeckEmitError("E407", detail=f"Load has no earlier Save: {prefix}")
            input_state, voltages = deepcopy(snapshots[prefix])
            _decision(decisions, "Load", prefix, "contacts.bias", voltages, "Saved voltages restored", "W208", index)
        start = deepcopy(voltages)
        goal = source.get("goal")
        contact = goal["contact"] if goal else electrodes[0]["name"]
        if contact not in start:
            raise DeckEmitError("E407", detail=f"unknown goal contact {contact}")
        current = bool(goal and "current_A_per_um" in goal)
        if current and contact not in source.get("current_mode_contacts", []):
            raise DeckEmitError("E407", detail="Current goal requires Set(mode current)")
        v0 = start[contact]
        v1 = _finite(goal["voltage_V"], "Goal") if goal and not current else v0
        voltages[contact] = v1
        entry = {"index": index, "ir_stage_index": group[0]["index"],
            "kind": "load" if kind == "load" else "equilibrium" if kind == "initial" or (goal and v0 == v1 and not current) else "sweep",
            "deck": f"simulation_stage{index:02d}.json", "csv": f"stage{index:02d}.csv",
            "input_state": input_state, "output_state": f"stage{index:02d}_state.h5",
            "electrode_voltages_start_V": start, "electrode_voltages_goal_V": deepcopy(voltages),
            "current_file": source["outputs"]["current_file"],
            "plt_write": "all_accepted_steps", "vtu_write": {"per_accepted_step": not options.plot_final_only,
                "series": [], "after_stage": []}, "required_t": [], "sidecar_overrides": []}
        times = []
        for member in group:
            for plot in member["outputs"]["plot_files"]:
                if plot["time"] is None:
                    entry["vtu_write"]["after_stage"].append(plot["vela_pattern"])
                else:
                    t = expand_times(plot["time"])
                    entry["vtu_write"]["series"].append({"file_prefix": plot["file_prefix"], "t": t,
                        "global_index_start": global_plot, "pattern": plot["vela_pattern"]})
                    global_plot += len(t)
                    times.extend(t)
        current_times = source["outputs"]["current_plot_times"]
        if current_times is not None:
            entry["plt_write"] = expand_times(current_times)
            times.extend(entry["plt_write"])
        entry["required_t"] = sorted(set(times))
        deck = deepcopy(base)
        deck["output_csv"] = entry["csv"]
        deck["contacts"] = [dict(c, bias=start[c["name"]]) for c in contacts]
        method = "newton" if any(c["type"] == "schottky" for c in contacts) else "gummel_newton"
        if kind == "load":
            method = "frozen_state"
        elif source["equations"] == ["Poisson"]:
            method = "poisson_only"
        deck["solver"].update(method=method, max_iter=source["iterations"] or math_cfg["iterations"] or 100)
        if input_state and method in ("gummel_newton", "newton"):
            # Coupled continuation must preserve the accepted state. A fresh
            # Gummel pass can destroy a good high-injection Newton predictor.
            deck["solver"].update(method="newton", warm_start=True,
                                  quasi_fermi_reference="contact_majority")
            # A small normalized block residual alone cannot qualify leakage
            # currents. Use the existing source-aware row and closure gates.
            deck["solver"].update(
                quasi_fermi_update_limit_V=0.1,
                carrier_row_convergence={"mode": "enforce", "eps_row": 0.001,
                    "scale_floor": 1e-30, "min_source_scale_fraction": 0.0,
                    "min_source_scale": 1e-18, "min_newton_max_iter": deck["solver"]["max_iter"]},
                continuity_row_scaling={"enabled": True, "flux_fraction": 0.0,
                    "scale_floor": 1e-30, "min_source_scale": 1e-18,
                    "min_weight": 1e-12, "max_weight": 1e18},
                global_continuity_closure={"mode": "enforce", "tolerance": 0.01, "source_floor": 1e-14})
            _decision(decisions, "Coupled", source["equations"], "solver.method", "newton",
                      "Continue the accepted state with referenced quasi-Fermi coordinates", stage=index)
        _decision(decisions, "Iterations", source["iterations"] or math_cfg["iterations"],
                  "solver.max_iter", deck["solver"]["max_iter"], "Stage value takes precedence", "W201", index)
        sweep = {"mode": "iv", "contact": contact, "current_contact": contact,
                 "start": v0, "stop": v1, "step": 1.0,
                 "initial_step": 1.0, "min_step": 1.0, "max_step": 1.0,
                 "growth_factor": 2.0, "shrink_factor": 0.5,
                 "write_vtk": True, "vtk_prefix": f"stage{index:02d}",
                 "vtk_format": "vtu", "vtk_fields": fields, "vtk_binary": options.vtu_binary,
                 "csv_all_contacts": True, "write_state_file": entry["output_state"],
                 "predictor": "linear" if math_cfg["extrapolate"] else "none"}
        if input_state:
            sweep["initial_state_file"] = input_state
        elif method != "poisson_only" and any(s["equations"] == ["Poisson"] for s in group):
            sweep["initialization"] = {"mode": "poisson_block"}
        if entry["kind"] in ("equilibrium", "load") or current:
            sweep["bias_points"] = [v0]
        if kind == "sweep":
            control = source["t_control"]
            controls = {k: control[k] if control[k] is not None else default for k, default in
                        (("initial_step", .1), ("min_step", .001), ("max_step", 1.0), ("increment", 2.0), ("decrement", 2.0))}
            if not 0 < controls["min_step"] <= controls["initial_step"] <= controls["max_step"] <= 1 or controls["increment"] < 1 or controls["decrement"] <= 1:
                raise DeckEmitError("E407", detail="invalid Quasistationary step controls")
            delta = abs(v1 - v0)
            if current:
                if index not in current_starts and not options.defer_current_times and any(t != 1 for t in times):
                    raise DeckEmitError("E407", detail="current-goal Time points need the measured initial current (WP5)")
                target = goal["current_A_per_um"]
                direction = -1 if target < 0 else 1
                current_points = [abs(target)]
                if index in current_starts:
                    initial = _finite(current_starts[index], "initial current")
                    if direction * initial > abs(target):
                        raise DeckEmitError("E407", detail="runner current goals must increase in the selected direction")
                    current_points = sorted(set(direction * (initial + t * (target - initial))
                                                for t in set(times) | {1.0}))
                sweep["voltage_to_current"] = {"switch_voltage_V": v0,
                    "current_points_A_per_um": current_points, "current_direction": direction}
                _decision(decisions, "Current", goal["current_A_per_um"], "sweep.voltage_to_current",
                          sweep["voltage_to_current"], "Runner closes the final current goal; subsequent stages must rebase measured voltages", stage=index)
            elif delta:
                intervals = current_times.get("intervals") if current_times else None
                spacing = delta / intervals if intervals else delta * controls["max_step"]
                step = spacing if v1 > v0 else -spacing
                converted = {"contact": contact, "mode": "voltage", "start_V": v0, "stop_V": v1,
                    "step_V": step, "initial_step_V": delta * controls["initial_step"],
                    "min_step_V": delta * controls["min_step"], "max_step_V": delta * controls["max_step"],
                    "growth_factor": controls["increment"], "shrink_factor": controls["decrement"],
                    "t_to_bias": {"v0": v0, "v1": v1}}
                entry["sweep"] = converted
                sweep.update(step=step, initial_step=converted["initial_step_V"],
                    min_step=converted["min_step_V"], max_step=converted["max_step_V"],
                    growth_factor=converted["growth_factor"], shrink_factor=1.0 / converted["shrink_factor"],
                    vtk_bias_points=[_bias(v0, v1, t) for t in entry["required_t"]])
                _decision(decisions, "Decrement", controls["decrement"], "sweep.shrink_factor",
                          sweep["shrink_factor"], "C4 retains the divisor; runner uses its reciprocal", stage=index)
            else:
                _decision(decisions, "Goal", goal, "sweep.bias_points", [v0],
                          "Zero voltage span becomes a single equilibrium point", "W211", index)
            if any(e.get("series_resistance_ohm_um", 0) and e["name"] != contact for e in electrodes):
                raise DeckEmitError("E407", detail="runner supports series resistance on the swept contact only")
            resistance = next((e.get("series_resistance_ohm_um", 0) for e in electrodes if e["name"] == contact), 0)
            if resistance:
                if current:
                    raise DeckEmitError("E407", detail="current control and series resistance are mutually exclusive")
                sweep["external_circuit"] = {"mode": "series_resistor", "resistance_ohm_um": resistance}
        if current:
            entry["sweep"] = {"contact": contact, "mode": "current", "start_V": v0, "stop_V": v0,
                "stop_current_A_per_um": goal["current_A_per_um"], "step_V": 1.0,
                "initial_step_V": 1.0, "min_step_V": 1.0, "max_step_V": 1.0,
                "growth_factor": 1.0, "shrink_factor": 1.0}
        if math_cfg["compute_ionization_integrals"]:
            diag = {"enabled": True}
            integral = math_cfg["break_at_ion_integral"]
            if integral:
                diag.update(break_rank=integral["rank"], break_value=integral["value"])
                entry.setdefault("break", {})["ionization_integral"] = deepcopy(integral)
            sweep.setdefault("diagnostics", {})["path_ionization_integrals"] = diag
        criteria = math_cfg["break_criteria"]
        if criteria:
            if len(criteria) != 1:
                raise DeckEmitError("E407", detail="runner supports only one Current break contact")
            criterion = criteria[0]
            entry.setdefault("break", {})["current"] = [{k: v for k, v in criterion.items()
                if k in ("contact", "abs_val_A_per_um", "min_val_A_per_um", "max_val_A_per_um") and v is not None}]
            breakdown = {"current_contact": criterion["contact"]}
            for old, new in (("abs_val_A_per_um", "current_abs_A_per_um"),
                             ("min_val_A_per_um", "current_min_A_per_um"), ("max_val_A_per_um", "current_max_A_per_um")):
                if criterion.get(old) is not None:
                    breakdown[new] = criterion[old]
            sweep["breakdown"] = breakdown
        if not math_cfg["aval_derivatives"] and "impact_ionization" in deck["solver"]:
            deck["solver"]["impact_ionization"]["source_jacobian"] = "frozen"
        deck["sweep"] = sweep
        if sidecar is not None and (not sidecar.get("apply_to") or index in sidecar["apply_to"]):
            _merge(deck, {k: v for k, v in sidecar.items() if k not in ("schema", "description", "apply_to", "avalanche_current_support_profile")},
                   "", entry["sidecar_overrides"], decisions, index)
            if profile := sidecar.get("avalanche_current_support_profile"):
                avalanche = deck["solver"]["impact_ionization"]
                if avalanche["model"] == "none" or avalanche.get("driving_force") != "quasi_fermi_gradient":
                    raise DeckEmitError("E409", file="sidecar", detail="avalanche profile requires active GradQF avalanche")
                bundle = AVALANCHE_PROFILES[profile]
                for key, value in bundle["mesh_geometry"].items():
                    if key in sidecar.get("mesh_geometry", {}) and sidecar["mesh_geometry"][key] != value:
                        raise DeckEmitError("E409", file="sidecar", detail="partial avalanche profile geometry override")
                avalanche.update(bundle["impact_ionization"], quasi_fermi_gradient_discretization="cell_gradient",
                                 contact_electric_field_fallback=True,
                                 source_jacobian="local_ad" if math_cfg["aval_derivatives"] else "frozen")
                deck["mesh_geometry"].update(bundle["mesh_geometry"])
                entry["sidecar_overrides"].append("avalanche_current_support_profile")
                _decision(decisions, "avalanche_current_support_profile", profile,
                          "solver.impact_ionization + mesh_geometry", bundle,
                          "Explicit atomic discretization; coefficients and physical models unchanged", "W207", index)
        if any(c["type"] == "schottky" for c in contacts) and deck["solver"]["method"] not in ("newton", "frozen_state"):
            raise DeckEmitError("E407", detail="thermionic Schottky contacts require the Newton path")
        geometry = deck["mesh_geometry"]
        if geometry.get("node_volume_policy") == "mixed_voronoi":
            if mesh["angles"]["max_deg"] > 90 + 1e-8 or mesh["angles"]["obtuse_element_count"] or not mesh["qualification"]["non_obtuse"]:
                raise DeckEmitError("E408", policy="mixed_voronoi", n=mesh["angles"]["obtuse_element_count"])
            geometry["require_non_obtuse"] = True
        if sidecar and kind == "load" and deck["solver"]["method"] != "frozen_state":
            raise DeckEmitError("E409", file="sidecar", detail="Load must retain frozen_state")
        plan["stages"].append(entry)
        decks.append(deck)
    if not decks:
        raise DeckEmitError("E407", detail="no executable stages")
    if sidecar and any(i >= len(decks) for i in sidecar.get("apply_to", [])):
        raise DeckEmitError("E409", file="sidecar", detail="apply_to names an absent stage")
    result = StagePlanResult(plan, decks, decisions, deepcopy({"ir": ir, "mesh": mesh, "report": report,
        "work_dir": work_dir, "cmd_name": cmd_name, "sha": sha, "options": options,
        "sidecar": sidecar, "rebases": rebases, "current_starts": current_starts,
        "saved_voltages": saved_voltages}))
    try:
        _validate_result(result)
    except DeckEmitError as error:
        if sidecar:
            raise DeckEmitError("E409", file="sidecar", detail=error.message) from error
        raise
    return result


def rebase_stage(result: StagePlanResult, stage_index: int, reached_voltages_V: dict,
                 *, update_saved_state: bool = False) -> StagePlanResult:
    """Regenerate this and later stages from measured voltages, preserving snapshots."""
    if not isinstance(stage_index, int) or not 0 <= stage_index < len(result.decks):
        raise DeckEmitError("E407", detail="invalid rebase stage index")
    context = deepcopy(result._context)
    if not context:
        raise DeckEmitError("E407", detail="rebase requires the original emitter result")
    context["rebases"][stage_index] = deepcopy(reached_voltages_V)
    if update_saved_state and stage_index:
        previous = result.stage_plan["stages"][stage_index - 1]
        context["saved_voltages"][previous["output_state"]] = {
            **previous["electrode_voltages_start_V"], **reached_voltages_V}
    return _emit(**context)


def bind_current_start(result: StagePlanResult, stage_index: int,
                       initial_current_A_per_um: float, reached_voltages_V: dict) -> StagePlanResult:
    """Resolve Current-goal Time points from the preceding accepted CSV state."""
    if result.stage_plan["stages"][stage_index].get("sweep", {}).get("mode") != "current":
        raise DeckEmitError("E407", detail="current start binding requires a Current stage")
    context = deepcopy(result._context)
    context["current_starts"][stage_index] = _finite(initial_current_A_per_um, "initial current")
    context["rebases"][stage_index] = deepcopy(reached_voltages_V)
    return _emit(**context)


def write_stage_plan(out_dir: Path, result: StagePlanResult) -> list[Path]:
    """Validate before writing stage_plan, decks and mapping decisions."""
    _validate_result(result)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    documents = {"stage_plan.json": result.stage_plan, "decisions.json": result.decisions}
    documents.update({s["deck"]: d for s, d in zip(result.stage_plan["stages"], result.decks)})
    paths = []
    for name, doc in documents.items():
        path = out_dir / name
        path.write_text(json.dumps(doc, indent=2, ensure_ascii=True, allow_nan=False) + "\n", encoding="utf-8")
        paths.append(path)
    return paths
