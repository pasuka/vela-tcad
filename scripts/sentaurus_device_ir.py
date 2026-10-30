#!/usr/bin/env python3
"""Sentaurus SDE ``.cmd`` -> Vela device IR (``vela.sentaurus_device_ir.v1``).

The frontend is fail-closed: every ``sde:``/``sdegeo:``/``sdedr:`` command must
be classified by an explicit handler. Commands that are recognised but not
translated into IR semantics are recorded as ``metadata_only`` together with the
reason they cannot be honoured; anything else aborts with the source file, line,
column, and the original command text.

The IR is split into five sections so that downstream consumers can depend on a
stable contract:

``geometry``      regions with their material and shape
``materials``     material assignment per region
``contacts``      contact sets and their boundary pick points
``doping``        profile definitions and their placements
``mesh_control``  refinement sizes/windows and their placements
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Sequence

from sentaurus_sexp import (
    Atom,
    BodyHandle,
    PureFunction,
    Environment,
    Sexp,
    SexpError,
    eval_number,
    eval_scalar,
    expect_string,
    parse_file,
    parse_position,
)


DEVICE_IR_SCHEMA = "vela.sentaurus_device_ir.v1"

# Species names that map onto Vela's donor/acceptor split. Anything else is
# rejected rather than silently bucketed into one of the two.
DONOR_SPECIES = {
    "PhosphorusActiveConcentration",
    "ArsenicActiveConcentration",
    "AntimonyActiveConcentration",
    "PhosphorusConcentration",
    "ArsenicConcentration",
}
ACCEPTOR_SPECIES = {
    "BoronActiveConcentration",
    "IndiumActiveConcentration",
    "AluminumActiveConcentration",
    "BoronConcentration",
}

# SDE material names mapped to the Vela material database.
MATERIAL_ALIASES = {
    "Silicon": "Si",
    "PolySilicon": "PolySilicon",
    "SiO2": "SiO2",
    "Oxide": "SiO2",
    "Nitride": "Nitride",
}


class SdeParseError(SexpError):
    """Raised when an SDE file cannot be translated fail-closed."""


@dataclass
class MetadataOnlyRecord:
    command: str
    location: str
    text: str
    reason: str

    def to_json(self) -> dict[str, str]:
        return {
            "command": self.command,
            "location": self.location,
            "text": self.text,
            "reason": self.reason,
        }


@dataclass
class DeviceIrBuilder:
    """Accumulates IR state while walking the SDE forms in source order."""

    source: str
    frontend: str = "offline-static"
    env: Environment = field(default_factory=Environment)
    regions: list[dict[str, Any]] = field(default_factory=list)
    contact_sets: dict[str, dict[str, Any]] = field(default_factory=dict)
    contact_order: list[str] = field(default_factory=list)
    current_contact: str | None = None
    boolean_mode: str = "ABA"
    inserted_vertices: list[dict[str, Any]] = field(default_factory=list)
    windows: dict[str, dict[str, Any]] = field(default_factory=dict)
    profiles: dict[str, dict[str, Any]] = field(default_factory=dict)
    general_parameters: dict[str, float] = field(default_factory=dict)
    placements: list[dict[str, Any]] = field(default_factory=list)
    refinement_sizes: dict[str, dict[str, Any]] = field(default_factory=dict)
    refinement_placements: list[dict[str, Any]] = field(default_factory=list)
    refinement_functions: list[dict[str, Any]] = field(default_factory=list)
    mesh_build: dict[str, Any] | None = None
    mesh_operations: list[dict[str, Any]] = field(default_factory=list)
    metadata_only: list[MetadataOnlyRecord] = field(default_factory=list)

    def note_metadata_only(self, form: Sexp, reason: str) -> None:
        self.metadata_only.append(
            MetadataOnlyRecord(
                command=str(form.head()),
                location=form.location(),
                text=form.command_text(),
                reason=reason,
            )
        )

    def region_names(self) -> set[str]:
        return {region["name"] for region in self.regions}


def _require_args(form: Sexp, minimum: int, usage: str) -> Sequence[Any]:
    args = form.items[1:]
    if len(args) < minimum:
        raise SdeParseError(
            f"{form.location()}: '{form.head()}' expects at least {minimum} "
            f"arguments ({usage}), got {len(args)} in '{form.command_text()}'")
    return args


def _rectangle_bounds(lower: Sequence[float],
                      upper: Sequence[float]) -> dict[str, list[float]]:
    x0, x1 = sorted((float(lower[0]), float(upper[0])))
    y0, y1 = sorted((float(lower[1]), float(upper[1])))
    return {"lower_left": [x0, y0], "upper_right": [x1, y1]}


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------


def _handle_define(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 2, "(define name value)")
    if len(args) != 2:
        raise SdeParseError(f"{form.location()}: define requires one value/body")
    name_node, value = args
    if isinstance(name_node, Sexp):
        if not isinstance(value, Sexp) or not name_node.items or len(name_node.items) > 9:
            raise SdeParseError(f"{form.location()}: invalid pure function definition")
        symbols = name_node.items
        if any(not isinstance(item, Atom) or not item.is_symbol for item in symbols):
            raise SdeParseError(f"{form.location()}: pure function names must be symbols")
        name = str(symbols[0].value)
        parameters = tuple(str(item.value) for item in symbols[1:])
        if (len(parameters) != len(set(parameters)) or name in builder.env.functions or
                name in builder.env or
                name in {"position", "+", "-", "*", "/", "abs", "sqrt", "exp", "log", "min", "max"}):
            raise SdeParseError(f"{form.location()}: duplicate/reserved pure function binding")
        allowed = {"+", "-", "*", "/", "sqrt", "abs", "exp", "log", "min", "max",
                   "position", "position-2d", name, *builder.env.functions}

        def check(node, depth=0):
            if depth > 32:
                raise SdeParseError(f"{form.location()}: pure function expression exceeds nesting budget")
            if isinstance(node, Sexp):
                if node.head() not in allowed:
                    raise SdeParseError(f"{node.location()}: impure/unsupported function body")
                for item in node.items[1:]:
                    check(item, depth+1)

        check(value)
        builder.env.functions[name] = PureFunction(parameters, value, form.location())
        return
    if not isinstance(name_node, Atom) or not name_node.is_symbol:
        raise SdeParseError(f"{form.location()}: define requires a symbol name")
    if str(name_node.value) in builder.env.functions:
        raise SdeParseError(f"{form.location()}: cannot replace a pure function with a scalar/body")
    if isinstance(value, Sexp) and value.head() == "sdegeo:create-rectangle":
        _handle_create_rectangle(value, builder)
        body = builder.regions[-1]
        result = BodyHandle(f"body:{body['order']}", body['name'], value.location())
    elif (isinstance(value, Atom) and value.is_symbol and str(value.value) in builder.env and
          isinstance(builder.env.get(str(value.value)), BodyHandle)):
        result = builder.env.get(str(value.value))
    else:
        result = eval_scalar(value, builder.env)
    builder.env.define(str(name_node.value), result)


def _handle_ignore(form: Sexp, builder: DeviceIrBuilder) -> None:
    """Commands with no effect on the IR (session bookkeeping)."""
    return None


def _handle_evaluated_wrapper_control(form: Sexp,
                                      builder: DeviceIrBuilder) -> None:
    if builder.frontend != "sentaurus2022-evaluated":
        raise SdeParseError(
            f"{form.location()}: '{form.head()}' is outside the declarative "
            f"offline-static Scheme subset in '{form.command_text()}'; use the "
            "sentaurus2022-evaluated frontend")
    if form.head() == "set!":
        args = _require_args(form, 2, "(set! evaluate-log-file #t)")
        name = args[0]
        if not isinstance(name, Atom) or str(name.value) != "evaluate-log-file":
            raise SdeParseError(
                f"{form.location()}: unexpected set! survived evaluated log: "
                f"'{form.command_text()}'")
    builder.note_metadata_only(
        form, "authoritative Sentaurus wrapper control already executed by SDE")


def _handle_create_rectangle(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 4, "(sdegeo:create-rectangle lower upper material region)")
    lower = parse_position(args[0], builder.env)
    upper = parse_position(args[1], builder.env)
    material = expect_string(args[2], builder.env, "material name")
    name = expect_string(args[3], builder.env, "region name")
    if name in builder.region_names():
        raise SdeParseError(
            f"{form.location()}: region '{name}' is defined more than once in "
            f"'{form.command_text()}'")
    if material not in MATERIAL_ALIASES:
        raise SdeParseError(
            f"{form.location()}: unsupported SDE material '{material}' in "
            f"'{form.command_text()}'; supported materials: "
            f"{', '.join(sorted(MATERIAL_ALIASES))}")
    bounds = _rectangle_bounds(lower, upper)
    if bounds["lower_left"][0] == bounds["upper_right"][0] or \
            bounds["lower_left"][1] == bounds["upper_right"][1]:
        raise SdeParseError(
            f"{form.location()}: degenerate rectangle in '{form.command_text()}'")
    builder.regions.append({
        "name": name,
        "material": material,
        "vela_material": MATERIAL_ALIASES[material],
        "shape": "rectangle",
        "order": len(builder.regions),
        # Boolean selection is stateful in SDE.  Capture it at body creation;
        # consulting the final global mode later gives the wrong geometry.
        "boolean_mode": builder.boolean_mode,
        "location": form.location(),
        **bounds,
    })


def _handle_create_polygon(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: 'sdegeo:create-polygon' is not supported by the "
        f"Vela SDE frontend; only axis-aligned rectangles can be meshed "
        f"conformally today. Offending command: '{form.command_text()}'")


def _handle_define_contact_set(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:define-contact-set name ...)")
    name = expect_string(args[0], builder.env, "contact name")
    if name in builder.contact_sets:
        raise SdeParseError(
            f"{form.location()}: contact set '{name}' is defined more than once "
            f"in '{form.command_text()}'")
    display = []
    for node in args[1:]:
        if isinstance(node, Atom) and node.is_number:
            display.append(float(node.value))
        elif isinstance(node, Atom) and node.is_string:
            display.append(str(node.value))
        elif isinstance(node, Sexp) and node.head() == "color:rgb":
            display.append({
                "command": "color:rgb",
                "values": [eval_number(item, builder.env) for item in node.items[1:]],
            })
        else:
            display.append(eval_scalar(node, builder.env))
    builder.contact_sets[name] = {
        "name": name, "pick_points": [], "display": display,
        "location": form.location(),
    }
    builder.contact_order.append(name)


def _handle_set_current_contact_set(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:set-current-contact-set name)")
    name = expect_string(args[0], builder.env, "contact name")
    if name not in builder.contact_sets:
        # Sentaurus 2022 creates/activates the named set when an explicit name
        # has not previously been declared.
        builder.contact_sets[name] = {
            "name": name, "pick_points": [], "display": [],
            "location": form.location(),
        }
        builder.contact_order.append(name)
    builder.current_contact = name


def _pick_point_from_edge_form(form: Sexp,
                               node: Any,
                               builder: DeviceIrBuilder) -> tuple[float, float]:
    """Resolve ``(find-edge-id (position ...))`` to a 2-D pick point."""
    if not isinstance(node, Sexp) or node.head() != "find-edge-id":
        location = node.location() if isinstance(node, (Atom, Sexp)) else form.location()
        raise SdeParseError(
            f"{location}: 'sdegeo:define-2d-contact' requires a "
            f"(find-edge-id (position ...)) argument in '{form.command_text()}'")
    inner = node.items[1:]
    if len(inner) != 1:
        raise SdeParseError(
            f"{node.location()}: (find-edge-id ...) expects exactly one position "
            f"in '{form.command_text()}'")
    point = parse_position(inner[0], builder.env)
    return (point[0], point[1])


def _handle_define_2d_contact(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:define-2d-contact edges [name])")
    edge_nodes: list[Any]
    first = args[0]
    if isinstance(first, Sexp) and first.head() == "list":
        edge_nodes = list(first.items[1:])
    else:
        edge_nodes = [first]

    if len(args) >= 2:
        name = expect_string(args[1], builder.env, "contact name")
    elif builder.current_contact is not None:
        name = builder.current_contact
    else:
        raise SdeParseError(
            f"{form.location()}: 'sdegeo:define-2d-contact' has no contact name "
            f"and no current contact set is selected in '{form.command_text()}'")

    if name not in builder.contact_sets:
        raise SdeParseError(
            f"{form.location()}: contact '{name}' is used before "
            f"'sdegeo:define-contact-set' in '{form.command_text()}'")

    for edge_node in edge_nodes:
        x, y = _pick_point_from_edge_form(form, edge_node, builder)
        builder.contact_sets[name]["pick_points"].append({
            "x": x,
            "y": y,
            "location": form.location(),
        })


def _handle_set_contact(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:set-contact entities [name])")
    if len(args) >= 3:
        operation = expect_string(args[2], builder.env, "contact operation")
        raise SdeParseError(
            f"{form.location()}: unsupported 2-D body contact operation "
            f"'{operation}' in '{form.command_text()}'; body contact, remove, "
            "and 3-D imprint are outside the offline-static capability")
    if len(args) >= 2:
        name = expect_string(args[1], builder.env, "contact name")
        if name not in builder.contact_sets:
            builder.contact_sets[name] = {
                "name": name, "pick_points": [], "display": [],
                "location": form.location(),
            }
            builder.contact_order.append(name)
        builder.current_contact = name
    elif builder.current_contact is not None:
        name = builder.current_contact
    else:
        raise SdeParseError(
            f"{form.location()}: 'sdegeo:set-contact' has neither an explicit "
            f"name nor an active contact set in '{form.command_text()}'")

    entities = list(args[0].items[1:]) if isinstance(args[0], Sexp) and \
        args[0].head() == "list" else [args[0]]
    for entity in entities:
        if isinstance(entity, Sexp) and entity.head() == "find-body-id":
            raise SdeParseError(
                f"{entity.location()}: body contact is unsupported by the 2-D "
                f"offline-static frontend in '{form.command_text()}'")
        x, y = _pick_point_from_edge_form(form, entity, builder)
        builder.contact_sets[name]["pick_points"].append({
            "x": x, "y": y, "location": form.location(),
        })


def _handle_insert_vertex(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:insert-vertex position)")
    if len(args) != 1:
        raise SdeParseError(
            f"{form.location()}: 2-D 'sdegeo:insert-vertex' expects exactly one "
            f"position in '{form.command_text()}'")
    point = parse_position(args[0], builder.env)
    builder.inserted_vertices.append({
        "x": point[0], "y": point[1], "location": form.location(),
    })


def _handle_set_default_boolean(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:set-default-boolean mode)")
    mode = expect_string(args[0], builder.env, "boolean mode")
    if mode.upper() not in {"ABA", "BAB"}:
        raise SdeParseError(
            f"{form.location()}: boolean mode '{mode}' is not implemented by "
            f"the offline-static frontend; AB, ABiA, ABiB and XX fail closed")
    builder.boolean_mode = mode.upper()


def _handle_refeval_window(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 3, "(sdedr:define-refeval-window name shape geometry ...)")
    name = expect_string(args[0], builder.env, "window name")
    shape = expect_string(args[1], builder.env, "window shape")
    if name in builder.windows:
        raise SdeParseError(
            f"{form.location()}: window '{name}' is defined more than once in "
            f"'{form.command_text()}'")
    if shape in {"Rectangle", "Line"}:
        if len(args) != 4:
            raise SdeParseError(
                f"{form.location()}: {shape} Ref/Eval window requires exactly "
                "two positions")
        lower = parse_position(args[2], builder.env)
        upper = parse_position(args[3], builder.env)
        if lower[:2] == upper[:2]:
            raise SdeParseError(f"{form.location()}: degenerate {shape} window")
        builder.windows[name] = {
            "name": name, "shape": shape.lower(),
            "start": [lower[0], lower[1]], "end": [upper[0], upper[1]],
            "location": form.location(),
        }
        if shape == "Rectangle":
            builder.windows[name].update(_rectangle_bounds(lower, upper))
    elif shape == "Polygon":
        if len(args) not in {3} and not all(
                isinstance(node, Sexp) and node.head() == "position"
                for node in args[2:]):
            raise SdeParseError(
                f"{form.location()}: Polygon Ref/Eval window requires a position "
                "list or a sequence of positions")
        raw = args[2]
        nodes = list(raw.items[1:]) if isinstance(raw, Sexp) and raw.head() == "list" \
            else list(args[2:])
        points = [[p[0], p[1]] for p in (parse_position(node, builder.env) for node in nodes)]
        if len(points) < 4 or points[0] != points[-1]:
            raise SdeParseError(
                f"{form.location()}: Polygon window needs at least three vertices "
                "and must repeat its first position at the end")
        signed_area2 = sum(
            points[index][0] * points[index + 1][1] -
            points[index + 1][0] * points[index][1]
            for index in range(len(points) - 1))
        if signed_area2 <= 0.0:
            raise SdeParseError(
                f"{form.location()}: Polygon Ref/Eval vertices must be "
                "counterclockwise as required by Sentaurus 2022")
        builder.windows[name] = {
            "name": name, "shape": "polygon", "points": points[:-1],
            "location": form.location(),
        }
    else:
        raise SdeParseError(
            f"{form.location()}: unsupported refeval window shape '{shape}' in "
            f"'{form.command_text()}'; supported: Line, Rectangle, Polygon")


def _handle_constant_profile(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 3, "(sdedr:define-constant-profile name species value)")
    name = expect_string(args[0], builder.env, "profile name")
    species = expect_string(args[1], builder.env, "species name")
    value = eval_number(args[2], builder.env)
    if species in DONOR_SPECIES:
        kind = "donor"
    elif species in ACCEPTOR_SPECIES:
        kind = "acceptor"
    else:
        raise SdeParseError(
            f"{form.location()}: unsupported doping species '{species}' in "
            f"'{form.command_text()}'; supported species: "
            f"{', '.join(sorted(DONOR_SPECIES | ACCEPTOR_SPECIES))}")
    if value < 0.0:
        raise SdeParseError(
            f"{form.location()}: negative doping concentration in "
            f"'{form.command_text()}'")
    if name in builder.profiles:
        raise SdeParseError(
            f"{form.location()}: doping profile '{name}' is defined more than "
            f"once in '{form.command_text()}'")
    builder.profiles[name] = {
        "name": name,
        "type": "constant",
        "species": species,
        "carrier": kind,
        "value_cm3": value,
        "location": form.location(),
    }


def _profile_carrier(form: Sexp, species: str) -> str:
    if species in DONOR_SPECIES:
        return "donor"
    if species in ACCEPTOR_SPECIES:
        return "acceptor"
    raise SdeParseError(
        f"{form.location()}: unsupported doping species '{species}' in "
        f"'{form.command_text()}'; supported species: "
        f"{', '.join(sorted(DONOR_SPECIES | ACCEPTOR_SPECIES))}")


def _keyword_arguments(form: Sexp, nodes: Sequence[Any], builder: DeviceIrBuilder,
                       flag_keys: set[str] | None = None) -> dict[str, Any]:
    """Parse Sentaurus' alternating quoted-key/value argument convention."""
    flags = flag_keys or set()
    result: dict[str, Any] = {}
    index = 0
    while index < len(nodes):
        key = expect_string(nodes[index], builder.env, "parameter keyword")
        if key in flags:
            result[key] = True
            index += 1
            continue
        if index + 1 >= len(nodes):
            raise SdeParseError(
                f"{form.location()}: parameter '{key}' has no value in "
                f"'{form.command_text()}'")
        value_node = nodes[index + 1]
        if isinstance(value_node, Atom) and value_node.is_string:
            value: Any = str(value_node.value)
        else:
            value = eval_number(value_node, builder.env)
        result[key] = value
        index += 2
    return result


def _handle_general_profile(form: Sexp, builder: DeviceIrBuilder) -> None:
    import math
    from sentaurus_general import initialise, expression_tree

    args = _require_args(
        form, 6, "(sdedr:define-analytical-profile name species init function default General)")
    if len(args) != 6 or expect_string(args[5], builder.env, "function kind") != "General":
        raise SdeParseError(f"{form.location()}: only the six-argument General profile is supported")
    name = expect_string(args[0], builder.env, "profile name")
    species = expect_string(args[1], builder.env, "species")
    if name in builder.profiles or species not in DONOR_SPECIES | ACCEPTOR_SPECIES:
        raise SdeParseError(f"{form.location()}: duplicate profile or unsupported species")
    try:
        parameters = initialise(
            expect_string(args[2], builder.env, "initializer"), builder.general_parameters)
        expression = expect_string(args[3], builder.env, "function")
        tree = expression_tree(expression, set(parameters) | {"x", "y", "z"})
        default = eval_number(args[4], builder.env)
        if not math.isfinite(default) or default < 0:
            raise ValueError("invalid General default")
    except ValueError as error:
        raise SdeParseError(f"{form.location()}: {error}") from error
    builder.profiles[name] = {
        "name": name, "type": "general", "species": species,
        "carrier": "donor" if species in DONOR_SPECIES else "acceptor",
        "parameters": parameters, "expression": expression, "expression_tree": tree,
        "default_value_cm3": default, "location": form.location(),
    }


def _handle_gaussian_profile(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 6, "(sdedr:define-gaussian-profile name species key value ...)")
    name = expect_string(args[0], builder.env, "profile name")
    species = expect_string(args[1], builder.env, "species name")
    if name in builder.profiles:
        raise SdeParseError(f"{form.location()}: duplicate profile '{name}'")
    params = _keyword_arguments(form, args[2:], builder, {"Gauss", "Erf"})
    if "PeakPos" not in params:
        raise SdeParseError(f"{form.location()}: Gaussian profile requires PeakPos")
    if not ({"PeakVal", "Dose"} & params.keys()):
        raise SdeParseError(f"{form.location()}: Gaussian profile requires PeakVal or Dose")
    if "PeakVal" in params and not (
            "StdDev" in params or "Length" in params or
            ({"ValueAtDepth", "Depth"} <= params.keys())):
        raise SdeParseError(
            f"{form.location()}: Gaussian PeakVal requires StdDev, Length, or "
            "ValueAtDepth plus Depth")
    lateral = "Gauss" if params.pop("Gauss", False) else \
        "Erf" if params.pop("Erf", False) else "Gauss"
    builder.profiles[name] = {
        "name": name, "type": "gaussian", "species": species,
        "carrier": _profile_carrier(form, species), "parameters": params,
        "lateral_type": lateral.lower(), "location": form.location(),
    }


def _placement_target(form: Sexp,
                      builder: DeviceIrBuilder,
                      target: str,
                      target_kind: str) -> None:
    if target_kind == "region" and target not in builder.region_names():
        raise SdeParseError(
            f"{form.location()}: placement references unknown region '{target}' "
            f"in '{form.command_text()}'")
    if target_kind == "window" and target not in builder.windows:
        raise SdeParseError(
            f"{form.location()}: placement references unknown window '{target}' "
            f"in '{form.command_text()}'")
    if target_kind == "material" and target not in MATERIAL_ALIASES:
        raise SdeParseError(
            f"{form.location()}: placement references unsupported material "
            f"'{target}' in '{form.command_text()}'")


def _handle_constant_profile_placement(form: Sexp,
                                       builder: DeviceIrBuilder,
                                       target_kind: str) -> None:
    args = _require_args(
        form, 3, f"({form.head()} name profile {target_kind})")
    name = expect_string(args[0], builder.env, "placement name")
    profile = expect_string(args[1], builder.env, "profile name")
    target = expect_string(args[2], builder.env, f"{target_kind} name")
    if profile not in builder.profiles:
        raise SdeParseError(
            f"{form.location()}: placement references unknown doping profile "
            f"'{profile}' in '{form.command_text()}'")
    if any(item["name"] == name for item in builder.placements):
        raise SdeParseError(
            f"{form.location()}: doping placement '{name}' is defined more than "
            f"once in '{form.command_text()}'")
    _placement_target(form, builder, target, target_kind)
    # T-2022.03 pp. 544-546: [decay-length ["Gauss"]] [replace].
    # Constant placements express NoReplace by omitting ``replace``; unlike an
    # analytical placement, the literal "NoReplace" is not a valid argument.
    replacement = "NoReplace"
    decay_length = 0.0
    decay = "none"
    options = list(args[3:])
    if options and isinstance(options[0], Atom) and options[0].is_string:
        raise SdeParseError(
            f"{form.location()}: constant placement replacement requires a "
            "preceding numeric decay length (use 0 for no decay); omit all "
            "options for NoReplace")
    if options and not (isinstance(options[0], Atom) and options[0].is_string):
        decay_length = eval_number(options.pop(0), builder.env)
        if decay_length < 0.0:
            raise SdeParseError(f"{form.location()}: negative decay length")
        decay = "erf" if decay_length > 0.0 else "none"
        if options and isinstance(options[0], Atom) and options[0].is_string and \
                str(options[0].value) == "Gauss":
            options.pop(0)
            decay = "gaussian" if decay_length > 0.0 else "none"
    if options:
        replacement = expect_string(options.pop(0), builder.env, "replacement mode")
        if replacement not in {"Replace", "LocalReplace"}:
            raise SdeParseError(
                f"{form.location()}: constant placement replacement must be "
                "'Replace' or 'LocalReplace'; omit it for NoReplace")
    if options:
        raise SdeParseError(
            f"{form.location()}: too many constant placement options in "
            f"'{form.command_text()}'")
    builder.placements.append({
        "name": name,
        "profile": profile,
        "target_kind": target_kind,
        "target": target,
        # Later placements win where windows overlap, matching SDE's
        # last-definition-wins evaluation order.
        "priority": len(builder.placements),
        "replacement": replacement,
        "decay_length": decay_length,
        "decay": decay,
        "location": form.location(),
    })


def _handle_analytical_profile_placement(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 6, "(sdedr:define-analytical-profile-placement name profile "
        "baseline direction replacement evaluation [restriction ...])")
    name = expect_string(args[0], builder.env, "placement name")
    profile = expect_string(args[1], builder.env, "profile name")
    baseline = expect_string(args[2], builder.env, "baseline name")
    direction = expect_string(args[3], builder.env, "placement direction")
    replacement = expect_string(args[4], builder.env, "replacement mode")
    evaluation = expect_string(args[5], builder.env, "evaluation mode")
    if profile not in builder.profiles or builder.profiles[profile]["type"] not in {"gaussian", "general"}:
        raise SdeParseError(f"{form.location()}: unknown analytical profile '{profile}'")
    if baseline not in builder.windows or builder.windows[baseline]["shape"] != "line":
        raise SdeParseError(f"{form.location()}: '{baseline}' is not a Line window")
    direction = "Both" if direction == "Symm" else direction
    if direction not in {"Positive", "Negative", "Both"}:
        raise SdeParseError(f"{form.location()}: unsupported direction '{direction}'")
    if replacement not in {"Replace", "LocalReplace", "NoReplace"}:
        raise SdeParseError(f"{form.location()}: unsupported replacement '{replacement}'")
    if evaluation not in {"Eval", "NoEval"}:
        raise SdeParseError(f"{form.location()}: unsupported evaluation mode '{evaluation}'")
    restriction_kind = None
    restriction = None
    restriction_decay_length = 0.0
    restriction_decay = "none"
    if len(args) > 6:
        if len(args) not in {9, 10}:
            raise SdeParseError(
                f"{form.location()}: analytical restriction must be "
                "eval-window, decay-length, optional 'Gauss', and window type")
        restriction = expect_string(args[6], builder.env, "restriction name")
        restriction_decay_length = eval_number(args[7], builder.env)
        if restriction_decay_length < 0.0:
            raise SdeParseError(
                f"{form.location()}: negative analytical restriction decay length")
        kind_index = 8
        if len(args) == 10:
            decay_name = expect_string(args[8], builder.env, "restriction decay type")
            if decay_name != "Gauss":
                raise SdeParseError(
                    f"{form.location()}: analytical restriction decay type must be 'Gauss'")
            restriction_decay = (
                "gaussian" if restriction_decay_length > 0.0 else "none")
            kind_index = 9
        elif restriction_decay_length > 0.0:
            restriction_decay = "erf"
        raw_kind = expect_string(
            args[kind_index], builder.env, "restriction kind").lower()
        restriction_kind = {"region": "region", "material": "material",
                            "window": "window", "evalwin": "window"}.get(raw_kind)
        if restriction_kind is None:
            raise SdeParseError(f"{form.location()}: unsupported restriction '{raw_kind}'")
        _placement_target(form, builder, restriction, restriction_kind)
    if any(item["name"] == name for item in builder.placements):
        raise SdeParseError(f"{form.location()}: duplicate placement '{name}'")
    builder.placements.append({
        "name": name, "profile": profile, "target_kind": "window",
        "target": baseline, "priority": len(builder.placements),
        "replacement": replacement, "direction": direction,
        "evaluation": evaluation, "restriction_kind": restriction_kind,
        "restriction": restriction,
        "restriction_decay_length": restriction_decay_length,
        "restriction_decay": restriction_decay,
        "location": form.location(),
    })


def _handle_refinement_size(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 3, "(sdedr:define-refinement-size name xmax ymax [xmin ymin])")
    name = expect_string(args[0], builder.env, "refinement name")
    values = [eval_number(arg, builder.env) for arg in args[1:]]
    if len(values) not in (2, 4, 6):
        raise SdeParseError(
            f"{form.location()}: 'sdedr:define-refinement-size' expects 2, 4, or "
            f"6 numeric extents, got {len(values)} in '{form.command_text()}'")
    if any(value <= 0.0 for value in values):
        raise SdeParseError(
            f"{form.location()}: refinement extents must be positive in "
            f"'{form.command_text()}'")
    if name in builder.refinement_sizes:
        raise SdeParseError(
            f"{form.location()}: refinement size '{name}' is defined more than "
            f"once in '{form.command_text()}'")
    entry: dict[str, Any] = {
        "name": name,
        "max_x": values[0],
        "max_y": values[1],
        "location": form.location(),
    }
    if len(values) >= 4:
        entry["min_x"] = values[2]
        entry["min_y"] = values[3]
    builder.refinement_sizes[name] = entry


def _handle_refinement_placement(form: Sexp,
                                 builder: DeviceIrBuilder,
                                 target_kind: str) -> None:
    args = _require_args(form, 3, f"({form.head()} name refinement {target_kind})")
    name = expect_string(args[0], builder.env, "placement name")
    refinement = expect_string(args[1], builder.env, "refinement name")
    target = expect_string(args[2], builder.env, f"{target_kind} name")
    if refinement not in builder.refinement_sizes:
        raise SdeParseError(
            f"{form.location()}: placement references unknown refinement "
            f"'{refinement}' in '{form.command_text()}'")
    if any(item["name"] == name for item in builder.refinement_placements):
        raise SdeParseError(
            f"{form.location()}: refinement placement '{name}' is defined more "
            f"than once in '{form.command_text()}'")
    _placement_target(form, builder, target, target_kind)
    builder.refinement_placements.append({
        "name": name,
        "refinement": refinement,
        "target_kind": target_kind,
        "target": target,
        "location": form.location(),
    })


def _handle_refinement_function(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 3, "(sdedr:define-refinement-function refinement function ...)")
    refinement = expect_string(args[0], builder.env, "refinement name")
    function = expect_string(args[1], builder.env, "refinement function")
    value_nodes = list(args[2:])
    variable = None
    if function not in {"MaxTransDiff", "MaxGradient", "MaxLenInt"}:
        # Profile-driven functions use the dataset name before the actual
        # refinement function, e.g. "DopingConcentration" "MaxTransDiff" 1.
        variable = function
        if not value_nodes:
            raise SdeParseError(
                f"{form.location()}: missing refinement function after variable "
                f"{variable!r}")
        function = expect_string(
            value_nodes.pop(0), builder.env, "refinement function")
    if refinement not in builder.refinement_sizes:
        raise SdeParseError(f"{form.location()}: unknown refinement '{refinement}'")
    if function not in {"MaxTransDiff", "MaxGradient", "MaxLenInt"}:
        raise SdeParseError(
            f"{form.location()}: refinement function '{function}' is unsupported; "
            "MaxInterval, multibox and offsetting fail closed")
    values: list[Any] = []
    for node in value_nodes:
        if isinstance(node, Atom) and node.is_string:
            values.append(str(node.value))
        else:
            values.append(eval_number(node, builder.env))
    if function in {"MaxTransDiff", "MaxGradient"} and not values:
        raise SdeParseError(f"{form.location()}: {function} requires a threshold")
    entry: dict[str, Any] = {
        "refinement": refinement, "function": function, "arguments": values,
        "location": form.location(),
    }
    if variable is not None:
        entry["variable"] = variable
    if function == "MaxLenInt":
        if len(values) < 3 or not all(isinstance(value, str) for value in values[:2]) \
                or not isinstance(values[2], (int, float)):
            raise SdeParseError(
                f"{form.location()}: MaxLenInt requires two interface names and "
                "a positive target length")
        target_length = float(values[2])
        if target_length <= 0.0:
            raise SdeParseError(f"{form.location()}: MaxLenInt length must be positive")
        index = 3
        factor = None
        if index < len(values) and isinstance(values[index], (int, float)):
            factor = float(values[index])
            if factor < 1.0:
                raise SdeParseError(
                    f"{form.location()}: MaxLenInt factor must be at least 1")
            index += 1
        flags = values[index:]
        if any(flag not in {"DoubleSide", "UseRegionNames"} for flag in flags) or \
                len(flags) != len(set(flags)):
            raise SdeParseError(
                f"{form.location()}: invalid or duplicate MaxLenInt option in "
                f"'{form.command_text()}'")
        entry.update({
            "interface": [values[0], values[1]],
            "value": target_length,
            "factor": factor,
            "double_side": "DoubleSide" in flags,
            "use_region_names": "UseRegionNames" in flags,
        })
    builder.refinement_functions.append(entry)


def _handle_multibox(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: '{form.head()}' multibox refinement is not "
        f"implemented by the Vela mesh generator. Offending command: "
        f"'{form.command_text()}'")


def _handle_fillet(form: Sexp, builder: DeviceIrBuilder) -> None:
    import math
    args = _require_args(form, 2, "(sdegeo:fillet-2d (find-vertex-id position) radius)")
    if len(args) != 2 or not isinstance(args[0], Sexp) or args[0].head() != "find-vertex-id" or len(args[0].items) != 2:
        raise SdeParseError(f"{form.location()}: fillet requires one explicit vertex pick")
    x,y,z = parse_position(args[0].items[1], builder.env)
    radius = eval_number(args[1], builder.env)
    candidates=[]
    for body in builder.regions:
        x0,y0=body['lower_left'];x1,y1=body['upper_right']
        if any(math.hypot(x-a,y-b)<1e-10 for a,b in [(x0,y0),(x1,y0),(x1,y1),(x0,y1)]):
            candidates.append(body)
    if len(candidates)!=1 or z!=0 or not math.isfinite(radius) or radius<=0:
        raise SdeParseError(f"{form.location()}: fillet requires an unambiguous convex rectangle corner and positive radius")
    body=candidates[0];x0,y0=body['lower_left'];x1,y1=body['upper_right']
    if radius>min(x1-x0,y1-y0)+1e-12 or body.get('fillets'):
        raise SdeParseError(f"{form.location()}: overlapping or oversized fillet is unsupported")
    body['fillets']=[{'corner':[x,y], 'radius_um':radius,'location':form.location()}]


def _handle_shortcut(form: Sexp, builder: DeviceIrBuilder) -> None:
    import math
    args=form.items[1:]
    tag=f"Shortcut.{len(builder.refinement_sizes)}"
    location=form.location()
    head=form.head()
    extent=max(max(abs(v) for v in b['lower_left']+b['upper_right']) for b in builder.regions)*4+1
    size={'name':tag,'max_x':extent,'max_y':extent,'location':location}
    if head=='sdedr:refine-box':
        if len(args) not in {4,6}:raise SdeParseError(f"{location}: refine-box expects scalar or three directional spacings")
        tag='Shortcut.'+expect_string(args[0],builder.env,'refinement name')
        start=parse_position(args[1],builder.env);end=parse_position(args[2],builder.env)
        if start[2]!=0 or end[2]!=0:raise SdeParseError(f'{location}: refine-box is 2-D only')
        h=[eval_number(a,builder.env) for a in args[3:]]
        if any(not math.isfinite(v) or v<=0 for v in h):raise SdeParseError(f'{location}: invalid directional spacing')
        if len(h)==1:h=h*3
        builder.windows[tag]={'name':tag,'shape':'rectangle',**_rectangle_bounds(start,end),'location':location}
        size.update(name=tag,max_x=h[0],max_y=h[1])
        placement={'name':tag,'refinement':tag,'target_kind':'window','target':tag,'location':location}
    elif head=='sdedr:refine-doping':
        # P4 supports the source scalar form; other optional syntax stays explicit.
        if len(args)!=1:raise SdeParseError(f"{location}: refine-doping supports scalar spacing only")
        h=eval_number(args[0],builder.env)
        size.update(min_x=h,min_y=h)
        placement={'name':tag,'refinement':tag,'target_kind':'material','target':'Silicon','location':location}
        builder.refinement_functions.append({'refinement':tag,'function':'MaxTransDiff','variable':'DopingConcentration','arguments':[1.0],'location':location})
    else:
        if len(args) not in {5,7}:raise SdeParseError(f"{location}: unsupported refine-interface arguments")
        kind=expect_string(args[0],builder.env,'scope')
        first=expect_string(args[1],builder.env,'first interface');second=expect_string(args[2],builder.env,'second interface')
        if kind not in {'region','material'} or expect_string(args[3],builder.env,'hlocal')!='hlocal':
            raise SdeParseError(f"{location}: refine-interface requires region/material and hlocal")
        value=eval_number(args[4],builder.env);factor=2.0
        if len(args)==7:
            if expect_string(args[5],builder.env,'factor')!='factor':raise SdeParseError(f"{location}: expected factor")
            factor=eval_number(args[6],builder.env)
        if not math.isfinite(value) or value<=0 or not math.isfinite(factor) or factor<1:
            raise SdeParseError(f"{location}: invalid interface spacing/factor")
        placement={'name':tag,'refinement':tag,'target_kind':kind,'target':first,'location':location}
        builder.refinement_functions.append({'refinement':tag,'function':'MaxLenInt','arguments':[first,second,value,factor],
            'interface':[first,second],'value':value,'factor':factor,'double_side':False,'use_region_names':kind=='region','location':location})
    if tag in builder.refinement_sizes or any(not math.isfinite(float(size[k])) or float(size[k])<=0 for k in ('max_x','max_y')):
        raise SdeParseError(f"{location}: duplicate or invalid shortcut sizing")
    if 'min_x' in size and (not math.isfinite(size['min_x']) or size['min_x']<=0):raise SdeParseError(f"{location}: invalid minimum spacing")
    builder.refinement_sizes[tag]=size
    builder.refinement_placements.append(placement)


def _handle_mesh_mirror(form: Sexp, builder: DeviceIrBuilder) -> None:
    import re
    args=_require_args(form,1,'(system:command known-tdx-mirror)')
    text=expect_string(args[0],builder.env,'mirror command')
    match=re.fullmatch(r'tdx\s+-mtt\s+-y\s+-M\s+0\s+-S\s+0\s+-ren\s+drain=source\s+([A-Za-z0-9_@.-]+)\s+([A-Za-z0-9_@.-]+)',text)
    if len(args)!=1 or not match or builder.mesh_build is None or match[1]!=builder.mesh_build['prefix']:
        raise SdeParseError(f"{form.location()}: unsupported system command; only matching local y=0 mesh reflection is allowed")
    builder.mesh_operations.append({'kind':'mirror','axis':'y','coordinate_um':0.0,'rename':{'drain':'source'},
                                    'input':match[1],'output':match[2],'location':form.location()})


def _handle_build_mesh(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = form.items[1:]
    values = [expect_string(arg, builder.env, "build-mesh argument")
              for arg in args
              if isinstance(arg, Atom) and arg.is_string]
    prefix = values[-1] if values else ""
    engine = values[0] if len(values) >= 2 else "snmesh"
    if engine and engine not in {"snmesh", "mesh", "noise"}:
        raise SdeParseError(
            f"{form.location()}: unsupported mesh engine '{engine}' in "
            f"'{form.command_text()}'")
    builder.mesh_build = {"engine": engine or "snmesh", "prefix": prefix}


def _handle_process_up_direction(form: Sexp, builder: DeviceIrBuilder) -> None:
    builder.note_metadata_only(
        form,
        "process up-direction affects SDE's own coordinate bookkeeping only",
    )


Handler = Callable[[Sexp, DeviceIrBuilder], None]

HANDLERS: dict[str, Handler] = {
    "define": _handle_define,
    "set!": _handle_evaluated_wrapper_control,
    "load": _handle_evaluated_wrapper_control,
    "sde:clear": _handle_ignore,
    "sde:build-mesh": _handle_build_mesh,
    "system:command": _handle_mesh_mirror,
    "sde:set-process-up-direction": _handle_process_up_direction,
    "sdegeo:create-rectangle": _handle_create_rectangle,
    "sdegeo:create-polygon": _handle_create_polygon,
    "sdegeo:define-contact-set": _handle_define_contact_set,
    "sdegeo:set-current-contact-set": _handle_set_current_contact_set,
    "sdegeo:define-2d-contact": _handle_define_2d_contact,
    "sdegeo:set-contact": _handle_set_contact,
    "sdegeo:insert-vertex": _handle_insert_vertex,
    "sdegeo:set-default-boolean": _handle_set_default_boolean,
    "sdegeo:fillet-2d": _handle_fillet,
    "sdedr:define-refeval-window": _handle_refeval_window,
    # SDE accepts both spellings for the evaluation window.
    "sdedr:define-refinement-window": _handle_refeval_window,
    "sdedr:define-constant-profile": _handle_constant_profile,
    "sdedr:define-gaussian-profile": _handle_gaussian_profile,
    "sdedr:define-analytical-profile": _handle_general_profile,
    "sdedr:define-constant-profile-placement":
        lambda form, builder: _handle_constant_profile_placement(form, builder, "window"),
    "sdedr:define-constant-profile-region":
        lambda form, builder: _handle_constant_profile_placement(form, builder, "region"),
    "sdedr:define-constant-profile-material":
        lambda form, builder: _handle_constant_profile_placement(form, builder, "material"),
    "sdedr:define-analytical-profile-placement": _handle_analytical_profile_placement,
    "sdedr:refine-box": _handle_shortcut,
    "sdedr:refine-doping": _handle_shortcut,
    "sdedr:refine-interface": _handle_shortcut,
    "sdedr:define-refinement-size": _handle_refinement_size,
    "sdedr:define-refinement-placement":
        lambda form, builder: _handle_refinement_placement(form, builder, "window"),
    "sdedr:define-refinement-region":
        lambda form, builder: _handle_refinement_placement(form, builder, "region"),
    "sdedr:define-refinement-material":
        lambda form, builder: _handle_refinement_placement(form, builder, "material"),
    "sdedr:define-refinement-function": _handle_refinement_function,
    "sdedr:define-multibox-size": _handle_multibox,
    "sdedr:define-multibox-placement": _handle_multibox,
}


def parse_sde_device_ir(path: Path,
                        frontend: str = "offline-static",
                        parameters: dict[str,str] | None = None) -> dict[str, Any]:
    """Parse an SDE ``.cmd`` file into ``vela.sentaurus_device_ir.v1``."""
    source = Path(path)
    import re
    from sentaurus_sexp import parse_forms
    text=source.read_text(encoding='utf-8-sig')
    used=set(re.findall(r'@([A-Za-z_][A-Za-z_0-9]*)@',text))
    parameters=parameters or {}
    if used-set(parameters):raise SdeParseError(f"{source}: unresolved preprocessing parameters: {sorted(used-set(parameters))}")
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',str(v)) for v in (parameters[name] for name in used)):
        raise SdeParseError(f"{source}: unsafe preprocessing parameter value")
    for name in used:text=text.replace('@'+name+'@',str(parameters[name]))
    forms = parse_forms(text,str(source))
    if frontend not in {"offline-static", "sentaurus2022-evaluated"}:
        raise ValueError(f"unknown SDE frontend '{frontend}'")
    builder = DeviceIrBuilder(source=str(source), frontend=frontend)

    for form in forms:
        if builder.mesh_operations:
            raise SdeParseError(f"{form.location()}: commands after final mesh reflection are unsupported")
        head = form.head()
        if head is None:
            raise SdeParseError(
                f"{form.location()}: expected a command symbol at the head of "
                f"'{form.command_text()}'")
        handler = HANDLERS.get(head)
        if handler is None:
            hint = ""
            if head in {"set!", "if", "cond", "do", "for-each", "lambda", "load"}:
                hint = ("; offline-static supports only a declarative Scheme subset; "
                        "use the sentaurus2022-evaluated frontend")
            raise SdeParseError(
                f"{form.location()}: unsupported SDE command '{head}' in "
                f"'{form.command_text()}'{hint}")
        handler(form, builder)

    # A declared-but-unassigned contact set is legal SDE, so a missing
    # 'sdegeo:define-2d-contact' pick point is not a parse error. It only
    # becomes fatal when something tries to materialise the boundary nodes,
    # which is where scripts/sentaurus_mesh_builder.py raises.
    return build_device_ir(builder)


def build_device_ir(builder: DeviceIrBuilder) -> dict[str, Any]:
    regions = [
        {
            "name": region["name"],
            "shape": region["shape"],
            "lower_left": region["lower_left"],
            "upper_right": region["upper_right"],
            "order": region["order"],
            "boolean_mode": region["boolean_mode"],
            "location": region["location"],
            **({"fillets": region["fillets"]} if region.get("fillets") else {}),
        }
        for region in builder.regions
    ]
    materials = [
        {
            "region": region["name"],
            "sde_material": region["material"],
            "vela_material": region["vela_material"],
        }
        for region in builder.regions
    ]
    contacts = [
        {
            "name": name,
            "display": builder.contact_sets[name]["display"],
            "pick_points": [
                {"x": point["x"], "y": point["y"], "location": point["location"]}
                for point in builder.contact_sets[name]["pick_points"]
            ],
            "location": builder.contact_sets[name]["location"],
        }
        for name in builder.contact_order
    ]
    doping = {
        "profiles": [builder.profiles[name] for name in sorted(builder.profiles)],
        "placements": list(builder.placements),
        "windows": [builder.windows[name] for name in sorted(builder.windows)],
    }
    mesh_control = {
        "refinement_sizes": [
            builder.refinement_sizes[name] for name in sorted(builder.refinement_sizes)
        ],
        "refinement_placements": list(builder.refinement_placements),
        "refinement_functions": list(builder.refinement_functions),
        "windows": doping["windows"],
        "build": builder.mesh_build,
        **({"operations": builder.mesh_operations} if builder.mesh_operations else {}),
    }
    return {
        "schema": DEVICE_IR_SCHEMA,
        "source": builder.source,
        "defines": {
            **{key: (value.to_json() if isinstance(value, BodyHandle) else value)
               for key, value in sorted(builder.env.values.items())},
            **{key: value.to_json() for key, value in sorted(builder.env.functions.items())}
        },
        "geometry": {
            "regions": regions,
            "inserted_vertices": [
                {"x": item["x"], "y": item["y"], "location": item["location"]}
                for item in builder.inserted_vertices
            ],
        },
        "materials": materials,
        "contacts": contacts,
        "doping": doping,
        "mesh_control": mesh_control,
        "metadata_only": [record.to_json() for record in builder.metadata_only],
    }


def legacy_summary(device_ir: dict[str, Any]) -> dict[str, Any]:
    """Project the device IR onto the historical ``sde`` summary shape.

    Existing reference-import flows consume this dictionary; keeping it as a
    projection avoids a second parser while the newer IR consumers migrate.
    """
    material_by_region = {
        item["region"]: item["sde_material"] for item in device_ir["materials"]
    }
    rectangles = [
        {
            "region": region["name"],
            "material": material_by_region[region["name"]],
            "lower_left": list(region["lower_left"]),
            "upper_right": list(region["upper_right"]),
        }
        for region in device_ir["geometry"]["regions"]
    ]
    profiles = {
        profile["name"]: {
            "species": profile["species"],
            "value": profile["value_cm3"],
        }
        for profile in device_ir["doping"]["profiles"]
    }
    placements = {
        placement["name"]: {
            "profile": placement["profile"],
            placement["target_kind"]: placement["target"],
            # Historical consumers only ever read "region"; keep the key present
            # for window placements too so they are not silently invisible.
            "region": placement["target"],
        }
        for placement in device_ir["doping"]["placements"]
    }
    windows = [
        {
            "name": window["name"],
            "shape": "Rectangle",
            "lower_left": list(window["lower_left"]),
            "upper_right": list(window["upper_right"]),
        }
        for window in device_ir["doping"]["windows"]
    ]
    return {
        "defines": dict(device_ir["defines"]),
        "geometry": {"rectangles": rectangles},
        "doping_profiles": profiles,
        "doping_placements": placements,
        "contacts": [
            {"name": contact["name"], "pick_points": contact["pick_points"]}
            for contact in device_ir["contacts"]
        ],
        "refinement_windows": windows,
        "device_ir": device_ir,
    }


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
