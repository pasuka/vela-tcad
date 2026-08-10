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
    "SiO2": "SiO2",
    "Oxide": "SiO2",
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
    env: Environment = field(default_factory=Environment)
    regions: list[dict[str, Any]] = field(default_factory=list)
    contact_sets: dict[str, dict[str, Any]] = field(default_factory=dict)
    contact_order: list[str] = field(default_factory=list)
    current_contact: str | None = None
    windows: dict[str, dict[str, Any]] = field(default_factory=dict)
    profiles: dict[str, dict[str, Any]] = field(default_factory=dict)
    placements: list[dict[str, Any]] = field(default_factory=list)
    refinement_sizes: dict[str, dict[str, Any]] = field(default_factory=dict)
    refinement_placements: list[dict[str, Any]] = field(default_factory=list)
    mesh_build: dict[str, Any] | None = None
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
    name_node = args[0]
    if not isinstance(name_node, Atom) or not name_node.is_symbol:
        raise SdeParseError(
            f"{form.location()}: (define ...) requires a symbol name in "
            f"'{form.command_text()}'")
    builder.env.define(str(name_node.value), eval_scalar(args[1], builder.env))


def _handle_ignore(form: Sexp, builder: DeviceIrBuilder) -> None:
    """Commands with no effect on the IR (session bookkeeping)."""
    return None


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
    builder.contact_sets[name] = {"name": name, "pick_points": []}
    builder.contact_order.append(name)


def _handle_set_current_contact_set(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:set-current-contact-set name)")
    name = expect_string(args[0], builder.env, "contact name")
    if name not in builder.contact_sets:
        raise SdeParseError(
            f"{form.location()}: contact set '{name}' is selected before it is "
            f"defined in '{form.command_text()}'")
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
    raise SdeParseError(
        f"{form.location()}: 'sdegeo:set-contact' selects boundary entities by "
        f"interactive body/edge id, which the Vela SDE frontend cannot resolve "
        f"deterministically. Use 'sdegeo:define-2d-contact' with an explicit "
        f"(find-edge-id (position ...)) pick point. Offending command: "
        f"'{form.command_text()}'")


def _handle_set_default_boolean(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(form, 1, "(sdegeo:set-default-boolean mode)")
    mode = expect_string(args[0], builder.env, "boolean mode")
    # "ABA" keeps previously created bodies intact where they overlap, which is
    # what the ordered-rectangle geometry model already assumes.
    if mode.upper() not in {"ABA", "BAB"}:
        raise SdeParseError(
            f"{form.location()}: unsupported boolean mode '{mode}' in "
            f"'{form.command_text()}'")
    builder.note_metadata_only(
        form,
        f"boolean mode '{mode}' recorded only; the rectangle geometry model "
        f"resolves overlaps by creation order",
    )


def _handle_refeval_window(form: Sexp, builder: DeviceIrBuilder) -> None:
    args = _require_args(
        form, 4, "(sdedr:define-refeval-window name shape lower upper)")
    name = expect_string(args[0], builder.env, "window name")
    shape = expect_string(args[1], builder.env, "window shape")
    if shape != "Rectangle":
        raise SdeParseError(
            f"{form.location()}: unsupported refeval window shape '{shape}' in "
            f"'{form.command_text()}'; only \"Rectangle\" is supported")
    lower = parse_position(args[2], builder.env)
    upper = parse_position(args[3], builder.env)
    if name in builder.windows:
        raise SdeParseError(
            f"{form.location()}: window '{name}' is defined more than once in "
            f"'{form.command_text()}'")
    builder.windows[name] = {
        "name": name,
        "shape": "rectangle",
        **_rectangle_bounds(lower, upper),
    }


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
    }


def _handle_gaussian_profile(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: 'sdedr:define-gaussian-profile' is not yet "
        f"translated by the Vela SDE frontend; accepting it would silently drop "
        f"the analytic decay parameters. Offending command: "
        f"'{form.command_text()}'")


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
    builder.placements.append({
        "name": name,
        "profile": profile,
        "target_kind": target_kind,
        "target": target,
        # Later placements win where windows overlap, matching SDE's
        # last-definition-wins evaluation order.
        "priority": len(builder.placements),
    })


def _handle_analytical_profile_placement(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: 'sdedr:define-analytical-profile-placement' is not "
        f"yet translated by the Vela SDE frontend. Offending command: "
        f"'{form.command_text()}'")


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
    })


def _handle_refinement_function(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: 'sdedr:define-refinement-function' drives adaptive "
        f"refinement that the Vela mesh generator does not implement. Offending "
        f"command: '{form.command_text()}'")


def _handle_multibox(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: '{form.head()}' multibox refinement is not "
        f"implemented by the Vela mesh generator. Offending command: "
        f"'{form.command_text()}'")


def _handle_fillet(form: Sexp, builder: DeviceIrBuilder) -> None:
    raise SdeParseError(
        f"{form.location()}: 'sdegeo:fillet-2d' rounds corners and changes the "
        f"device outline, which the rectangle geometry model cannot represent. "
        f"Offending command: '{form.command_text()}'")


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
    "sde:clear": _handle_ignore,
    "sde:build-mesh": _handle_build_mesh,
    "sde:set-process-up-direction": _handle_process_up_direction,
    "sdegeo:create-rectangle": _handle_create_rectangle,
    "sdegeo:create-polygon": _handle_create_polygon,
    "sdegeo:define-contact-set": _handle_define_contact_set,
    "sdegeo:set-current-contact-set": _handle_set_current_contact_set,
    "sdegeo:define-2d-contact": _handle_define_2d_contact,
    "sdegeo:set-contact": _handle_set_contact,
    "sdegeo:set-default-boolean": _handle_set_default_boolean,
    "sdegeo:fillet-2d": _handle_fillet,
    "sdedr:define-refeval-window": _handle_refeval_window,
    # SDE accepts both spellings for the evaluation window.
    "sdedr:define-refinement-window": _handle_refeval_window,
    "sdedr:define-constant-profile": _handle_constant_profile,
    "sdedr:define-gaussian-profile": _handle_gaussian_profile,
    "sdedr:define-constant-profile-placement":
        lambda form, builder: _handle_constant_profile_placement(form, builder, "window"),
    "sdedr:define-constant-profile-region":
        lambda form, builder: _handle_constant_profile_placement(form, builder, "region"),
    "sdedr:define-analytical-profile-placement": _handle_analytical_profile_placement,
    "sdedr:define-refinement-size": _handle_refinement_size,
    "sdedr:define-refinement-placement":
        lambda form, builder: _handle_refinement_placement(form, builder, "window"),
    "sdedr:define-refinement-region":
        lambda form, builder: _handle_refinement_placement(form, builder, "region"),
    "sdedr:define-refinement-function": _handle_refinement_function,
    "sdedr:define-multibox-size": _handle_multibox,
    "sdedr:define-multibox-placement": _handle_multibox,
}


def parse_sde_device_ir(path: Path) -> dict[str, Any]:
    """Parse an SDE ``.cmd`` file into ``vela.sentaurus_device_ir.v1``."""
    source = Path(path)
    forms = parse_file(source)
    builder = DeviceIrBuilder(source=str(source))

    for form in forms:
        head = form.head()
        if head is None:
            raise SdeParseError(
                f"{form.location()}: expected a command symbol at the head of "
                f"'{form.command_text()}'")
        handler = HANDLERS.get(head)
        if handler is None:
            raise SdeParseError(
                f"{form.location()}: unsupported SDE command '{head}' in "
                f"'{form.command_text()}'")
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
            "pick_points": [
                {"x": point["x"], "y": point["y"]}
                for point in builder.contact_sets[name]["pick_points"]
            ],
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
        "windows": doping["windows"],
        "build": builder.mesh_build,
    }
    return {
        "schema": DEVICE_IR_SCHEMA,
        "source": builder.source,
        "defines": {
            key: builder.env.get(key) for key in sorted(builder.env.values)
        },
        "geometry": {"regions": regions},
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
