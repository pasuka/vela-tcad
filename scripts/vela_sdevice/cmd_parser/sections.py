"""Section parsers: File, Electrode, Physics, Plot, CurrentPlot and Math."""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from .lexer import Token
from .syntax import (Block, ParseContext, SyntaxIssue, TokenStream, assignment_text, coerce_value,
                     is_number, parse_int, parse_number)

# --- File ---------------------------------------------------------------------------

_FILE_MAIN = {"grid": "grid", "plot": "plot", "parameter": "parameter", "parameters": "parameter",
              "current": "current", "output": "output"}
_FILE_EXTRA_OUTPUT = {"newtonplot", "acextract", "extraction", "nonlocalplot", "pmipath",
                      "savesubstructure", "currentplot", "plotexplicit", "spiceparam"}


def cmd_stem(path_name: str) -> str:
    for suffix in ("_des.cmd", ".cmd"):
        if path_name.lower().endswith(suffix):
            return path_name[: -len(suffix)]
    return path_name


def _complete(raw: str, extension: str) -> str:
    base = raw.replace("\\", "/").rsplit("/", 1)[-1]
    return raw if "." in base else raw + extension


def resolve_output_name(raw: Optional[str]) -> str:
    """Observed T-2022.03 behaviour: literal only when already ending in _des.log."""
    if raw is None:
        return "output_des.log"
    return raw if raw.endswith("_des.log") else raw + "_des.log"


def _read_assignments(stream: TokenStream, section: str) -> List[Tuple[Token, Optional[Token], List[Token]]]:
    """``key [= value] [(...)]`` triples; the trailing group tokens are returned raw."""
    items: List[Tuple[Token, Optional[Token], List[Token]]] = []
    while not stream.at_end():
        key = stream.next()
        if key.kind != "word":
            raise SyntaxIssue(f"unexpected '{key.value}' in {section} section", key.line, key.value)
        value: Optional[Token] = None
        group: List[Token] = []
        if stream.at_punct("="):
            stream.next()
            value = stream.take_value()
        if stream.at_punct("("):
            _, group, _ = stream.take_group()
        items.append((key, value, group))
    return items


def parse_file_block(block: Block, ctx: ParseContext, path_name: str) -> dict:
    section = "File"
    stem = cmd_stem(path_name)
    seen: Dict[str, Tuple[Token, Optional[Token]]] = {}
    ignored: List[dict] = []
    try:
        for key, value, _group in _read_assignments(TokenStream(block.body, ctx.source), section):
            name = key.value.lower()
            if name in _FILE_MAIN:
                canonical = _FILE_MAIN[name]
                if canonical in seen:
                    ctx.warn("W109", key.line, section, assignment_text(key, value), kw=key.value)
                if value is None:
                    raise SyntaxIssue(f"File keyword '{key.value}' needs a value", key.line, key.value)
                seen[canonical] = (key, value)
            elif name in _FILE_EXTRA_OUTPUT:
                ignored.append({"keyword": key.value, "value": value.value if value else None,
                                "line": key.line, "finding_id": "W107"})
                ctx.warn("W107", key.line, section, assignment_text(key, value), kw=key.value)
            else:
                ctx.error("E112", key.line, section, assignment_text(key, value), kw=key.value)
    except SyntaxIssue as issue:
        ctx.syntax_error(issue, section)

    files: dict = {"grid": None, "parameter": None, "current": None, "plot": None, "output": None,
                   "ignored": ignored}
    if "grid" in seen:
        key, value = seen["grid"]
        raw = value.value
        resolved = _complete(raw, ".tdr")
        lower = resolved.lower()
        kind = "tdr" if lower.endswith(".tdr") else "dfise" if lower.endswith(".grd") else "unknown"
        files["grid"] = {"raw": raw, "resolved": resolved, "kind": kind, "line": key.line}
        if kind == "unknown":
            ctx.error("E114", key.line, section, assignment_text(key, value), name=raw)
    if "parameter" in seen:
        key, value = seen["parameter"]
        files["parameter"] = {"raw": value.value, "resolved": _complete(value.value, ".par"),
                              "keyword": "Parameters" if key.value.lower() == "parameters" else "Parameter",
                              "line": key.line}
    else:
        ctx.error("E301", block.name.line, section, ctx.source.line_text(block.name.line))
    if "current" in seen:
        key, value = seen["current"]
        files["current"] = {"raw": value.value, "resolved": _complete(value.value, "_des.plt"), "line": key.line}
    else:
        files["current"] = {"raw": None, "resolved": f"{stem}_des.plt", "line": None}
    if "plot" in seen:
        key, value = seen["plot"]
        resolved = _complete(value.value, "_des.tdr")
        line: Optional[int] = key.line
        raw: Optional[str] = value.value
    else:
        resolved, line, raw = f"{stem}_des.tdr", None, None
    vela = resolved[: -len(".tdr")] + ".vtu" if resolved.lower().endswith(".tdr") else resolved + ".vtu"
    files["plot"] = {"raw": raw, "resolved": resolved, "vela": vela, "line": line}
    if "output" in seen:
        key, value = seen["output"]
        files["output"] = {"raw": value.value, "resolved": resolve_output_name(value.value), "line": key.line}
    else:
        files["output"] = {"raw": None, "resolved": resolve_output_name(None), "line": None}
    return files


def skeleton_files(path_name: str) -> dict:
    stem = cmd_stem(path_name)
    return {
        "grid": None,
        "parameter": None,
        "current": {"raw": None, "resolved": f"{stem}_des.plt", "line": None},
        "plot": {"raw": None, "resolved": f"{stem}_des.tdr", "vela": f"{stem}_des.vtu", "line": None},
        "output": {"raw": None, "resolved": "output_des.log", "line": None},
        "ignored": [],
    }


# --- Electrode -----------------------------------------------------------------------

_ELECTRODE_UNSUPPORTED = {"charge", "current", "distresist", "eqohmic", "srdoping", "fgcap", "extraction",
                          "cyclicnorm", "temperature", "surfaceresistance", "power", "schottkyresist",
                          "eschottkyresist", "hschottkyresist", "chargeneutral", "resistivity"}


def parse_electrode_block(block: Block, ctx: ParseContext) -> List[dict]:
    section = "Electrode"
    electrodes: List[dict] = []
    stream = TokenStream(block.body, ctx.source)
    try:
        while not stream.at_end():
            opener, inner, _closer = stream.take_group()
            if opener.value != "{":
                raise SyntaxIssue("electrode definitions must be enclosed in braces", opener.line, opener.value)
            electrodes.append(_parse_electrode(inner, opener.line, ctx))
    except SyntaxIssue as issue:
        ctx.syntax_error(issue, section)
    return electrodes


def _parse_electrode(tokens: List[Token], line: int, ctx: ParseContext) -> dict:
    section = "Electrode"
    stream = TokenStream(tokens, ctx.source)
    name: Optional[str] = None
    voltage: Optional[float] = None
    attributes: dict = {}
    while not stream.at_end():
        key = stream.next()
        if key.kind != "word":
            raise SyntaxIssue(f"unexpected '{key.value}' in electrode", key.line, key.value)
        lower = key.value.lower()
        value: Optional[Token] = None
        group: Optional[List[Token]] = None
        if stream.at_punct("="):
            stream.next()
            if lower == "name" and stream.peek() is not None and stream.peek().kind == "word" \
                    and stream.peek().value.lower() == "regexp":
                regexp = stream.next()
                _, inner, closer = stream.take_group()
                ctx.error("E117", key.line, section, ctx.source.slice(key.start, closer.end))
                name = inner[0].value if inner else regexp.value
                continue
            if stream.at_punct("("):
                _, group, _ = stream.take_group()
            else:
                value = stream.take_value()
        if stream.at_punct("(") and group is None:
            _, group, _ = stream.take_group()
        text = assignment_text(key, value)
        if lower == "name":
            name = value.value if value else None
        elif lower == "voltage":
            if group is not None or value is None:
                ctx.error("E108", key.line, section, text, kw=key.value)
            else:
                voltage = parse_number(value, "voltage")
        elif lower == "barrier":
            attributes["barrier_V"] = parse_number(value, "Barrier")
        elif lower == "workfunction":
            attributes["workfunction_eV"] = parse_number(value, "Workfunction")
        elif lower in ("resist", "resistor"):
            attributes["resist_ohm_um"] = parse_number(value, "Resist")
            if lower == "resistor":
                ctx.warn("W112", key.line, section, text)
        elif lower == "schottky":
            attributes["schottky"] = True
        elif lower == "erecvelocity":
            attributes["e_rec_velocity_cm_s"] = parse_number(value, "eRecVelocity")
        elif lower == "hrecvelocity":
            attributes["h_rec_velocity_cm_s"] = parse_number(value, "hRecVelocity")
        elif lower == "material":
            if group is not None or value is None:
                ctx.error("E108", key.line, section, text, kw=key.value)
            else:
                attributes["material"] = value.value
        elif lower == "poisson":
            if value is not None and value.value.lower() == "dirichlet":
                ctx.warn("W113", key.line, section, text, kw=key.value, value=value.value)
            else:
                ctx.error("E108", key.line, section, text, kw=key.value)
        elif lower == "areafactor":
            if value is not None and is_number(value.value) and float(value.value) == 1.0:
                ctx.warn("W113", key.line, section, text, kw=key.value, value=value.value)
            else:
                ctx.error("E108", key.line, section, text, kw=key.value)
        elif lower in _ELECTRODE_UNSUPPORTED:
            ctx.error("E108", key.line, section, text, kw=key.value)
        else:
            ctx.error("E111", key.line, section, text, kw=key.value)
    if name is None:
        raise SyntaxIssue("electrode without Name", line, ctx.source.line_text(line))
    return {"name": name, "voltage_V": voltage if voltage is not None else 0.0,
            "attributes": attributes, "line": line}


# --- Physics -------------------------------------------------------------------------

_INSULATORS = {"oxide", "sio2", "nitride", "si3n4", "hfo2", "al2o3", "gas", "ambient", "vacuum", "insulator"}
_DRIVING_FORCES = {"gradquasifermi": "GradQuasiFermi", "eparallel": "Eparallel", "electricfield": "ElectricField",
                   "carriertempdrive": "CarrierTempDrive"}
_AVALANCHE_MODELS = {"vanoverstraeten": "VanOverstraeten", "okuto": "Okuto", "lackner": "Lackner",
                     "unibo": "UniBo", "unibo2": "UniBo2", "hatakeyama": "Hatakeyama"}
_BAND2BAND_MODELS = {"e1": "E1", "e1_5": "E1_5", "e2": "E2", "schenk": "Schenk", "hurkx": "Hurkx",
                     "nonlocalpath": "NonLocalPath"}
_BGN_MODELS = {"oldslotboom": "OldSlotboom", "slotboom": "Slotboom", "bennettwilson": "BennettWilson",
               "delalamo": "delAlamo", "nobandgapnarrowing": "NoBandGapNarrowing", "tablebgn": "TableBGN",
               "jainroulston": "JainRoulston", "nofermi": "NoFermi"}
_UNSUPPORTED_TOP = {"thermodynamic", "hydrodynamic", "etemperature", "htemperature", "hquantumpotential",
                    "aniso", "piezo", "mechanics", "optics", "raytrace", "ldp", "montecarlo",
                    "heterointerface", "thermionic", "incompleteionization", "hqp", "eqp"}
_DEFAULT_VELA = {
    "Mobility.DopingDependence": "solver.mobility.model=masetti",
    "Mobility.ConstantMobility": "solver.mobility.model=constant",
    "Mobility.Enormal": "solver.mobility.model+=_lombardi",
    "Recombination.Auger": "solver.recombination+=auger",
    "EffectiveIntrinsicDensity.OldSlotboom": "solver.bandgap_narrowing=old_slotboom",
    "EffectiveIntrinsicDensity.NoBandGapNarrowing": "solver.bandgap_narrowing=none",
    "Fermi": "solver.carrier_statistics.model=fermi_dirac",
    "eQuantumPotential": "solver.electron_quantum_potential.enabled=true",
}
_DRIVING_FORCE_VELA = {"GradQuasiFermi": "quasi_fermi_gradient", "Eparallel": "electric_field_parallel",
                       "ElectricField": "electric_field"}


class PhysicsScope:
    def __init__(self, scope: Optional[str], target: Optional[str], line: int) -> None:
        self.scope = scope
        self.target = target
        self.line = line
        self.models: List[dict] = []          # verbatim model_entry trees
        self.effective: Dict[str, dict] = {}  # canonical key -> effective entry
        self.fixed_charges: List[dict] = []
        self.temperature_K: Optional[float] = None

    @property
    def is_interface(self) -> bool:
        return self.scope in ("MaterialInterface", "RegionInterface")

    @property
    def is_insulator(self) -> bool:
        return self.scope == "Material" and (self.target or "").lower() in _INSULATORS


def _parse_model_args(stream: TokenStream) -> List[dict]:
    args: List[dict] = []
    while not stream.at_end():
        token = stream.next()
        if token.kind == "punct":
            if token.value == "(":
                # bare nested group such as Traps((FixedCharge ...))
                stream.i -= 1
                _, inner, _ = stream.take_group()
                args.append({"name": "", "args": _parse_model_args(TokenStream(inner, stream.source))})
                continue
            raise SyntaxIssue(f"unexpected '{token.value}' in model arguments", token.line, token.value)
        entry: dict = {"name": token.value, "_line": token.line}
        if stream.at_punct("="):
            stream.next()
            if stream.at_punct("("):
                _, inner, _ = stream.take_group()
                entry["args"] = _parse_model_args(TokenStream(inner, stream.source))
            else:
                value = stream.take_value()
                entry["value"] = coerce_value(value)
        if stream.at_punct("("):
            _, inner, _ = stream.take_group()
            entry["args"] = _parse_model_args(TokenStream(inner, stream.source))
        args.append(entry)
    return args


def _arg_names(args: List[dict]) -> List[str]:
    return [a["name"] for a in args]


def _strip_lines(node: dict) -> dict:
    node.pop("_line", None)
    for arg in node.get("args", []):
        _strip_lines(arg)
    return node


def _line_of(node: dict, fallback: int) -> int:
    return int(node.get("_line", fallback))


def _add_effective(scope: PhysicsScope, key: str, parameters: dict, status: str = "supported",
                   vela: Optional[str] = None, finding_id: Optional[str] = None) -> None:
    entry = {"key": key, "parameters": parameters, "status": status, "vela": vela}
    if finding_id:
        entry["finding_id"] = finding_id
    scope.effective[key] = entry


def _canonicalise_model(model: dict, scope: PhysicsScope, ctx: ParseContext, line: int) -> None:
    section = "Physics"
    name = model["model"]
    lower = name.lower()
    args = model.get("args", [])
    if lower == "mobility":
        for arg in args:
            a = arg["name"].lower()
            sub = arg.get("args", [])
            if a in ("dopingdep", "dopingdependence"):
                _add_effective(scope, "Mobility.DopingDependence", {}, vela=_DEFAULT_VELA["Mobility.DopingDependence"])
            elif a == "constantmobility":
                _add_effective(scope, "Mobility.ConstantMobility", {}, vela=_DEFAULT_VELA["Mobility.ConstantMobility"])
            elif a in ("highfieldsaturation", "ehighfieldsaturation", "hhighfieldsaturation"):
                key = {"highfieldsaturation": "Mobility.HighFieldSaturation",
                       "ehighfieldsaturation": "Mobility.eHighFieldSaturation",
                       "hhighfieldsaturation": "Mobility.hHighFieldSaturation"}[a]
                params: dict = {}
                vela = "solver.mobility.model+=_field"
                for option in sub:
                    force = _DRIVING_FORCES.get(option["name"].lower())
                    if force is None:
                        ctx.error("E107", _line_of(option, line), section, option["name"], text=option["name"])
                        continue
                    params["driving_force"] = force
                    if force in _DRIVING_FORCE_VELA:
                        vela += f"; solver.mobility.high_field_driving_force={_DRIVING_FORCE_VELA[force]}"
                    else:
                        _add_effective(scope, key, params, "unsupported", None, "E106")
                        ctx.error("E106", _line_of(arg, line), section, arg["name"], key=key)
                        break
                else:
                    _add_effective(scope, key, params, vela=vela)
            elif a == "enormal":
                if sub and any(o["name"].lower() == "ialmob" for o in sub):
                    _add_effective(scope, "Mobility.Enormal", {"model": "IALMob"}, "unsupported", None, "E106")
                    ctx.error("E106", _line_of(arg, line), section, "Enormal(IALMob)", key="Mobility.Enormal.IALMob")
                else:
                    params = {"model": sub[0]["name"]} if sub else {}
                    _add_effective(scope, "Mobility.Enormal", params, vela=_DEFAULT_VELA["Mobility.Enormal"])
            elif a in ("ialmob", "carriercarrierscattering", "phumob", "philipsunifiedmobility",
                       "hydrodynamic", "thinlayer", "tomobility", "negativejunction", "monteCarlo".lower()):
                key = f"Mobility.{arg['name']}"
                _add_effective(scope, key, {}, "unsupported", None, "E106")
                ctx.error("E106", _line_of(arg, line), section, arg["name"], key=key)
            else:
                ctx.error("E107", _line_of(arg, line), section, arg["name"], text=arg["name"])
    elif lower == "recombination":
        for arg in args:
            a = arg["name"].lower()
            sub = arg.get("args", [])
            if a == "srh":
                params = {}
                vela = "solver.recombination+=srh"
                ok = True
                for option in sub:
                    o = option["name"].lower()
                    if o in ("dopingdep", "dopingdependence"):
                        params["doping_dependence"] = True
                        vela += "; solver.srh_doping_dependence=true"
                    elif o in ("tempdependence", "tempdep"):
                        params["temperature_dependence"] = True
                        vela += "; solver.srh_temperature_dependence=true"
                    elif o in ("electricfield", "tunneling", "hurkx", "nonlocalpath", "expdep"):
                        ok = False
                        _add_effective(scope, "Recombination.SRH", {"option": option["name"]}, "unsupported", None, "E106")
                        ctx.error("E106", _line_of(option, line), section, f"SRH({option['name']})", key="Recombination.SRH")
                        break
                    else:
                        ctx.error("E107", _line_of(option, line), section, option["name"], text=option["name"])
                        ok = False
                        break
                if ok:
                    _add_effective(scope, "Recombination.SRH", params, vela=vela)
            elif a == "auger":
                _add_effective(scope, "Recombination.Auger", {}, vela=_DEFAULT_VELA["Recombination.Auger"])
            elif a == "avalanche":
                params = {}
                unsupported = None
                for option in sub:
                    o = option["name"].lower()
                    if o in _AVALANCHE_MODELS:
                        params["model"] = _AVALANCHE_MODELS[o]
                        if params["model"] != "VanOverstraeten":
                            unsupported = params["model"]
                    elif o in _DRIVING_FORCES:
                        params["driving_force"] = _DRIVING_FORCES[o]
                    elif o in ("eparallel", "gradquasifermi"):
                        params["driving_force"] = _DRIVING_FORCES[o]
                    elif "value" in option and o in ("ebarrier", "hbarrier", "refdens", "a", "b"):
                        params[option["name"]] = option["value"]
                    else:
                        ctx.error("E107", _line_of(option, line), section, option["name"], text=option["name"])
                if unsupported:
                    _add_effective(scope, "Recombination.Avalanche", params, "unsupported", None, "E106")
                    ctx.error("E106", _line_of(arg, line), section, f"Avalanche({unsupported})", key="Recombination.Avalanche")
                else:
                    vela = "solver.impact_ionization.model=van_overstraeten"
                    force = params.get("driving_force")
                    if force in _DRIVING_FORCE_VELA:
                        vela += f"; solver.impact_ionization.driving_force={_DRIVING_FORCE_VELA[force]}"
                    _add_effective(scope, "Recombination.Avalanche", params, vela=vela)
            elif a == "band2band":
                params = {}
                for option in sub:
                    o = option["name"].lower()
                    if o == "model" and "value" in option:
                        params["model"] = _BAND2BAND_MODELS.get(str(option["value"]).lower(), str(option["value"]))
                    elif o in _BAND2BAND_MODELS:
                        params["model"] = _BAND2BAND_MODELS[o]
                    elif "value" in option:
                        params[option["name"]] = option["value"]
                    else:
                        ctx.error("E107", _line_of(option, line), section, option["name"], text=option["name"])
                if params.get("model") == "E2":
                    _add_effective(scope, "Recombination.Band2Band", params, vela="solver.band_to_band.model=e2")
                else:
                    _add_effective(scope, "Recombination.Band2Band", params, "unsupported", None, "E106")
                    ctx.error("E106", _line_of(arg, line), section, f"Band2Band({params.get('model', '')})", key="Recombination.Band2Band")
            elif a in ("radiative", "cdl", "trapassisted", "surfacesrh", "excitonic", "cdl1", "cdl2", "tat"):
                key = f"Recombination.{arg['name']}"
                _add_effective(scope, key, {}, "unsupported", None, "E106")
                ctx.error("E106", _line_of(arg, line), section, arg["name"], key=key)
            else:
                ctx.error("E107", _line_of(arg, line), section, arg["name"], text=arg["name"])
    elif lower == "effectiveintrinsicdensity":
        for arg in args:
            a = arg["name"].lower()
            canonical = _BGN_MODELS.get(a)
            if canonical is None:
                ctx.error("E107", _line_of(arg, line), section, arg["name"], text=arg["name"])
                continue
            key = f"EffectiveIntrinsicDensity.{canonical}"
            if key in _DEFAULT_VELA:
                _add_effective(scope, key, {}, vela=_DEFAULT_VELA[key])
            elif canonical == "NoFermi":
                continue
            else:
                _add_effective(scope, key, {}, "unsupported", None, "E106")
                ctx.error("E106", _line_of(arg, line), section, arg["name"], key=key)
    elif lower == "fermi":
        _add_effective(scope, "Fermi", {}, vela=_DEFAULT_VELA["Fermi"])
    elif lower == "-fermi":
        return
    elif lower == "equantumpotential":
        params = {}
        for option in args:
            if "value" in option:
                params[option["name"]] = option["value"]
        _add_effective(scope, "eQuantumPotential", params, vela=_DEFAULT_VELA["eQuantumPotential"])
    elif lower == "traps":
        _parse_traps(model, scope, ctx, line)
    elif lower in _UNSUPPORTED_TOP or lower.startswith("ialmob"):
        _add_effective(scope, name, {}, "unsupported", None, "E106")
        ctx.error("E106", line, section, name, key=name)
    elif lower == "areafactor":
        pass
    else:
        ctx.error("E107", line, section, name, text=name)


def _parse_traps(model: dict, scope: PhysicsScope, ctx: ParseContext, line: int) -> None:
    section = "Physics"
    # Accept Traps(Conc=x FixedCharge) and Traps((FixedCharge SpatialShape=Uniform Conc=x)).
    args = model.get("args", [])
    if len(args) == 1 and args[0]["name"] == "" and "args" in args[0]:
        args = args[0]["args"]
    names = {a["name"].lower() for a in args}
    conc = next((a.get("value") for a in args if a["name"].lower() == "conc"), None)
    shape = next((a.get("value") for a in args if a["name"].lower() == "spatialshape"), "Uniform")
    if "fixedcharge" in names and isinstance(conc, (int, float)) and str(shape).lower() == "uniform" \
            and names <= {"fixedcharge", "conc", "spatialshape"} and scope.scope is not None:
        entry = {"scope": scope.scope, "target": scope.target, "line": line}
        entry["conc_cm2" if scope.is_interface else "conc_cm3"] = float(conc)
        scope.fixed_charges.append(entry)
        return
    _add_effective(scope, "Traps", {}, "unsupported", None, "E106")
    ctx.error("E106", line, section, ctx.source.line_text(line), key="Traps")


def parse_physics_block(block: Block, ctx: ParseContext) -> PhysicsScope:
    section = "Physics"
    scope_kind: Optional[str] = None
    target: Optional[str] = None
    if block.scope is not None:
        try:
            stream = TokenStream(block.scope, ctx.source)
            key = stream.next()
            stream.expect_punct("=")
            value = stream.take_value()
            scope_kind = {"material": "Material", "region": "Region", "materialinterface": "MaterialInterface",
                          "regioninterface": "RegionInterface"}.get(key.value.lower())
            if scope_kind is None:
                raise SyntaxIssue(f"unknown Physics scope '{key.value}'", key.line, key.value)
            target = value.value
        except SyntaxIssue as issue:
            ctx.syntax_error(issue, section)
    scope = PhysicsScope(scope_kind, target, block.name.line)
    stream = TokenStream(block.body, ctx.source)
    try:
        while not stream.at_end():
            head = stream.next()
            if head.kind != "word":
                raise SyntaxIssue(f"unexpected '{head.value}' in Physics", head.line, head.value)
            if stream.at_punct("="):
                stream.next()
                value = stream.take_value()
                if head.value.lower() == "temperature":
                    scope.temperature_K = parse_number(value, "Temperature")
                elif head.value.lower() == "areafactor":
                    if not (is_number(value.value) and float(value.value) == 1.0):
                        ctx.error("E106", head.line, section, assignment_text(head, value), key="AreaFactor")
                else:
                    ctx.error("E107", head.line, section, head.value, text=head.value)
                continue
            model: dict = {"model": head.value, "args": [], "line": head.line}
            if stream.at_punct("("):
                _, inner, _ = stream.take_group()
                model["args"] = _parse_model_args(TokenStream(inner, ctx.source))
            scope.models.append(model)
            if scope.is_interface or scope.is_insulator:
                if head.value.lower() == "traps":
                    _parse_traps(model, scope, ctx, head.line)
                else:
                    ctx.error("E106", head.line, section, head.value, key=head.value)
            else:
                _canonicalise_model(model, scope, ctx, head.line)
            _strip_lines(model)
    except SyntaxIssue as issue:
        ctx.syntax_error(issue, section)
    return scope


def combine_physics(scopes: List[PhysicsScope], ctx: ParseContext) -> dict:
    section = "Physics"
    global_scope = PhysicsScope(None, None, 0)
    for scope in scopes:
        if scope.scope is None:
            global_scope.models.extend(scope.models)
            global_scope.effective.update(scope.effective)
            if scope.temperature_K is not None:
                global_scope.temperature_K = scope.temperature_K
    semiconductor_scopes = [s for s in scopes if s.scope is not None and not s.is_interface and not s.is_insulator]
    effective_sets: List[Tuple[PhysicsScope, Dict[str, dict]]] = []
    for scope in semiconductor_scopes:
        merged = dict(global_scope.effective)
        merged.update(scope.effective)
        effective_sets.append((scope, merged))
    effective = dict(global_scope.effective)
    if effective_sets:
        first_scope, first = effective_sets[0]
        for scope, other in effective_sets[1:]:
            if {k: (v["parameters"], v["status"]) for k, v in first.items()} != \
                    {k: (v["parameters"], v["status"]) for k, v in other.items()}:
                differing = sorted(set(first) ^ set(other) | {k for k in set(first) & set(other)
                                                              if first[k]["parameters"] != other[k]["parameters"]})
                ctx.error("E113", scope.line, section, ctx.source.line_text(scope.line),
                          detail=f"{first_scope.scope}={first_scope.target} vs {scope.scope}={scope.target}: "
                                 + ", ".join(differing))
        effective = first
    fixed_charges: List[dict] = []
    for scope in scopes:
        fixed_charges.extend(scope.fixed_charges)
    return {
        "temperature_K": global_scope.temperature_K if global_scope.temperature_K is not None else 300.0,
        "global": global_scope.models,
        "scoped": [{"scope": s.scope, "target": s.target, "models": s.models, "line": s.line}
                   for s in scopes if s.scope is not None],
        "effective": [effective[key] for key in sorted(effective)],
        "fixed_charges": fixed_charges,
    }


# --- Plot / CurrentPlot --------------------------------------------------------------

_PLOT_ALIASES = {
    "potential": "ElectrostaticPotential", "electrostaticpotential": "ElectrostaticPotential",
    "equasifermi": "eQuasiFermiPotential", "hquasifermi": "hQuasiFermiPotential",
    "equasifermipotential": "eQuasiFermiPotential", "hquasifermipotential": "hQuasiFermiPotential",
    "edensity": "eDensity", "hdensity": "hDensity",
    "doping": "DopingConcentration", "dopingconcentration": "DopingConcentration",
    "donorconcentration": "DonorConcentration", "acceptorconcentration": "AcceptorConcentration",
    "spacecharge": "SpaceCharge", "electricfield": "ElectricField",
    "ecurrent": "eCurrentDensity", "hcurrent": "hCurrentDensity", "totalcurrent": "TotalCurrentDensity",
    "ecurrentdensity": "eCurrentDensity", "hcurrentdensity": "hCurrentDensity",
    "totalcurrentdensity": "TotalCurrentDensity", "conductioncurrent": "ConductionCurrentDensity",
    "conductioncurrentdensity": "ConductionCurrentDensity", "displacementcurrent": "DisplacementCurrentDensity",
    "displacementcurrentdensity": "DisplacementCurrentDensity",
    "emobility": "eMobility", "hmobility": "hMobility", "evelocity": "eVelocity", "hvelocity": "hVelocity",
    "srh": "srhRecombination", "srhrecombination": "srhRecombination",
    "auger": "AugerRecombination", "augerrecombination": "AugerRecombination",
    "totalrecombination": "TotalRecombination",
    "band2band": "Band2BandGeneration", "band2bandgeneration": "Band2BandGeneration",
    "avalanchegeneration": "AvalancheGeneration", "eavalanchegeneration": "eAvalancheGeneration",
    "havalanchegeneration": "hAvalancheGeneration", "impactionization": "AvalancheGeneration",
    "eimpactionization": "eAvalancheGeneration", "himpactionization": "hAvalancheGeneration",
    "ealphaavalanche": "eAlphaAvalanche", "halphaavalanche": "hAlphaAvalanche",
    "eionintegral": "eIonIntegral", "hionintegral": "hIonIntegral", "meanionintegral": "MeanIonIntegral",
    "bandgap": "BandGap", "bandgapnarrowing": "BandgapNarrowing", "effectivebandgap": "EffectiveBandGap",
    "affinity": "ElectronAffinity", "electronaffinity": "ElectronAffinity",
    "conductionband": "ConductionBandEnergy", "valenceband": "ValenceBandEnergy",
    "conductionbandenergy": "ConductionBandEnergy", "valencebandenergy": "ValenceBandEnergy",
    "effectiveintrinsicdensity": "EffectiveIntrinsicDensity", "intrinsicdensity": "IntrinsicDensity",
    "eeffectivestatedensity": "eEffectiveStateDensity", "heffectivestatedensity": "hEffectiveStateDensity",
    "equantumpotential": "eQuantumPotential", "hquantumpotential": "hQuantumPotential",
    "eeparallel": "eEparallel", "heparallel": "hEparallel", "eenormal": "eEnormal", "henormal": "hEnormal",
    "egradquasifermi": "eGradQuasiFermi", "hgradquasifermi": "hGradQuasiFermi",
    "etemperature": "eTemperature", "htemperature": "hTemperature",
    "temperature": "LatticeTemperature", "latticetemperature": "LatticeTemperature",
    "elifetime": "eLifetime", "hlifetime": "hLifetime", "equasifermienergy": "eQuasiFermiEnergy",
    "hquasifermienergy": "hQuasiFermiEnergy", "etrappedcharge": "eTrappedCharge", "htrappedcharge": "hTrappedCharge",
}
_PLOT_OPTIONS = {"-loadable", "loadable", "datasetsfromgrid", "-plotloadable", "plotloadable", "plotexplicit",
                 "-plotexplicit", "explicit", "-explicit"}


def plot_data_name(name: str) -> Optional[str]:
    return _PLOT_ALIASES.get(name.lower())


def parse_plot_block(block: Block, ctx: ParseContext) -> dict:
    section = "Plot"
    fields: List[dict] = []
    explicit = False
    for token in block.body:
        if token.kind != "word":
            ctx.syntax_error(SyntaxIssue(f"unexpected '{token.value}' in Plot", token.line, token.value), section)
            continue
        lower = token.value.lower()
        if lower in _PLOT_OPTIONS:
            if "explicit" in lower and not lower.startswith("-"):
                explicit = True
            ctx.warn("W108", token.line, section, token.value, kw=token.value)
            continue
        name, component = token.value, "scalar"
        if "/" in token.value:
            name, suffix = token.value.split("/", 1)
            component = {"vector": "vector", "element": "element"}.get(suffix.lower(), "scalar")
        data_name = plot_data_name(name)
        if data_name is None:
            ctx.warn("W104", token.line, section, token.value, name=name)
        fields.append({"name": name, "data_name": data_name, "component": component, "line": token.line})
    return {"explicit": explicit, "fields": fields}


def parse_current_plot_block(block: Block, ctx: ParseContext) -> dict:
    section = "CurrentPlot"
    entries: List[dict] = []
    stream = TokenStream(block.body, ctx.source)
    try:
        while not stream.at_end():
            head = stream.next()
            if head.kind != "word":
                raise SyntaxIssue(f"unexpected '{head.value}' in CurrentPlot", head.line, head.value)
            last = head
            kind = None
            if stream.at_punct("("):
                opener, inner, closer = stream.take_group()
                last = closer
                inner_words = [t.value.lower() for t in inner if t.kind == "word"]
                if inner and inner[0].kind == "punct" and inner[0].value == "(":
                    kind = "W105"
                elif any(w in ("integrate", "maximum", "minimum", "average") for w in inner_words):
                    kind = "W106"
            if stream.at_punct("="):
                stream.next()
                last = stream.take_value()
            text = ctx.source.slice(head.start, last.end)
            if kind:
                ctx.warn(kind, head.line, section, text, text=text)
            entries.append({"text": text, "line": head.line, "finding_id": kind})
    except SyntaxIssue as issue:
        ctx.syntax_error(issue, section)
    return {"entries": entries}


# --- Math ----------------------------------------------------------------------------

_MATH_RECORDED = {"notdamped", "transient", "method", "submethod", "cnormprint", "wallclock",
                  "checkundefinedmodels", "parameterinheritance", "cylindrical", "rhsmin", "rhsfactor",
                  "rhsandupdateconvergence", "extendedprecision", "digits_dg", "errref_dg",
                  "numberofsolverthreads", "acmethod", "acsubmethod", "iterations_dg", "damped",
                  "stackednewtonweights", "newdiscretization", "elementvolumeavalanche", "avalpostprocessing",
                  "avalflatelementexclusion", "elementedgecurrent", "avalanche", "gradquasifermi",
                  "exitonfailure_dg", "spice_gmin", "traps", "currentweighting", "conductivemetal",
                  "-metalconductivity", "metalconductivity", "extrapolate_dg", "printvolume", "surfaceconductivity",
                  "elementcurrent", "eflux", "hflux", "quasifermitunneling"}


def parse_math_block(block: Block, ctx: ParseContext, math: dict) -> None:
    section = "Math"
    stream = TokenStream(block.body, ctx.source)
    try:
        while not stream.at_end():
            key = stream.next()
            if key.kind != "word":
                raise SyntaxIssue(f"unexpected '{key.value}' in Math", key.line, key.value)
            name = key.value
            lower = name.lower()
            negated = lower.startswith("-")
            base = lower[1:] if negated else lower
            group: Optional[List[Token]] = None
            value: Optional[Token] = None
            if stream.at_punct("("):
                _, group, _ = stream.take_group()
            if stream.at_punct("="):
                stream.next()
                value = stream.take_value()
            elif stream.at_punct("{"):
                _, group, _ = stream.take_group()
            text = assignment_text(key, value)
            if base == "iterations":
                math["iterations"] = parse_int(value, "Iterations")
            elif base == "digits":
                math["digits"] = coerce_value(value)
            elif base == "relerrcontrol":
                math["rel_err_control"] = not negated
            elif base == "errref":
                carrier = group[0].value.lower() if group else ""
                if carrier == "electron":
                    math["err_ref_electron_cm3"] = parse_number(value, "ErrRef")
                elif carrier == "hole":
                    math["err_ref_hole_cm3"] = parse_number(value, "ErrRef")
                else:
                    math["ignored"].append({"keyword": f"ErrRef({group[0].value if group else ''})",
                                            "value": value.value if value else None,
                                            "line": key.line, "finding_id": "W102"})
                    ctx.warn("W102", key.line, section, ctx.source.slice(key.start, (value or key).end),
                             kw=f"ErrRef({group[0].value if group else ''})")
            elif base == "extrapolate":
                math["extrapolate"] = not negated
                for option in group or []:
                    if option.kind == "word" and not is_number(option.value):
                        ctx.warn("W111", key.line, section, option.value, kw=option.value)
            elif base == "derivatives":
                math["derivatives"] = not negated
            elif base == "avalderivatives":
                math["aval_derivatives"] = not negated
            elif base == "computeionizationintegrals":
                math["compute_ionization_integrals"] = not negated
            elif base == "breakationintegral":
                numbers = [t for t in (group or []) if t.kind == "word"]
                if len(numbers) != 2:
                    raise SyntaxIssue("BreakAtIonIntegral expects (rank value)", key.line, text)
                math["break_at_ion_integral"] = {"rank": parse_int(numbers[0], "rank"),
                                                 "value": parse_number(numbers[1], "value")}
            elif base == "breakcriteria":
                _parse_break_criteria(group or [], ctx, math, key.line)
            elif base == "exitonfailure":
                math["exit_on_failure"] = not negated
            elif base == "numberofthreads":
                math["number_of_threads"] = parse_int(value, "NumberOfThreads")
            elif base in _MATH_RECORDED:
                math["ignored"].append({"keyword": name, "value": value.value if value else None,
                                        "line": key.line, "finding_id": "W102"})
                ctx.warn("W102", key.line, section, text, kw=name)
            else:
                math["ignored"].append({"keyword": name, "value": value.value if value else None,
                                        "line": key.line, "finding_id": "W101"})
                ctx.warn("W101", key.line, section, text, kw=name)
    except SyntaxIssue as issue:
        ctx.syntax_error(issue, section)


def _parse_break_criteria(tokens: List[Token], ctx: ParseContext, math: dict, line: int) -> None:
    section = "Math"
    stream = TokenStream(tokens, ctx.source)
    while not stream.at_end():
        head = stream.next()
        if head.kind != "word":
            raise SyntaxIssue(f"unexpected '{head.value}' in BreakCriteria", head.line, head.value)
        inner: List[Token] = []
        if stream.at_punct("("):
            _, inner, _ = stream.take_group()
        if head.value.lower() != "current":
            math["ignored"].append({"keyword": f"BreakCriteria.{head.value}", "value": None,
                                    "line": head.line, "finding_id": "W102"})
            ctx.warn("W102", head.line, section, head.value, kw=f"BreakCriteria.{head.value}")
            continue
        entry = {"kind": "Current", "contact": "", "abs_val_A_per_um": None, "min_val_A_per_um": None,
                 "max_val_A_per_um": None, "line": head.line}
        inner_stream = TokenStream(inner, ctx.source)
        while not inner_stream.at_end():
            key = inner_stream.next()
            inner_stream.expect_punct("=")
            value = inner_stream.take_value()
            lower = key.value.lower()
            if lower == "contact":
                entry["contact"] = value.value
            elif lower == "absval":
                entry["abs_val_A_per_um"] = parse_number(value, "AbsVal")
            elif lower == "minval":
                entry["min_val_A_per_um"] = parse_number(value, "MinVal")
            elif lower == "maxval":
                entry["max_val_A_per_um"] = parse_number(value, "MaxVal")
            else:
                raise SyntaxIssue(f"unknown BreakCriteria option '{key.value}'", key.line, key.value)
        math["break_criteria"].append(entry)


def default_math() -> dict:
    return {
        "iterations": 50, "digits": 5, "rel_err_control": True,
        "err_ref_electron_cm3": None, "err_ref_hole_cm3": None,
        "extrapolate": False, "derivatives": True, "aval_derivatives": True,
        "compute_ionization_integrals": False, "break_at_ion_integral": None,
        "break_criteria": [], "exit_on_failure": False, "number_of_threads": 1, "ignored": [],
    }
