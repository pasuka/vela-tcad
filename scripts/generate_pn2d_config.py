#!/usr/bin/env python3
"""Render validated, reproducible PN2D and BVmethods configurations."""

from __future__ import annotations

import argparse
import copy
import json
import math
import re
from pathlib import Path, PureWindowsPath
from typing import Any


REPO = Path(__file__).resolve().parents[1]
TEMPLATE_DIR = REPO / "configs" / "templates"
TEMPLATES = {
    "bvmethods_nmos_external_resistor": (
        TEMPLATE_DIR / "bvmethods_nmos_external_resistor.template.json"
    ),
    "bvmethods_nmos_voltage_to_current": (
        TEMPLATE_DIR / "bvmethods_nmos_voltage_to_current.template.json"
    ),
    "pn2d_iv": TEMPLATE_DIR / "pn2d_iv.template.json",
    "pn2d_bv": TEMPLATE_DIR / "pn2d_bv.template.json",
}
PLACEHOLDER = re.compile(r"^\$\{([a-zA-Z_][a-zA-Z0-9_]*)\}$")
BV_AVALANCHE_CURRENT_SUPPORT_PROFILES = {
    "element_edge_sg_gss_laux": {
        "current_approximation": "element_edge_sg_gss_laux",
        "source_mapping_mode": "element_vertex_box_measure",
        "cell_reconstructed_midpoint_density": "bernoulli",
    },
    "legacy_cell_reconstructed": {
        "current_approximation": "cell_reconstructed",
        "source_mapping_mode": "triangle_gss_gradqf_truncated",
        "cell_reconstructed_midpoint_density": "gss_logistic",
    },
}
BV_CONFIGURATION_PROFILES = {
    "element_edge_sg_gss_laux": {
        "impact_ionization": BV_AVALANCHE_CURRENT_SUPPORT_PROFILES[
            "element_edge_sg_gss_laux"
        ],
        "mesh_geometry": {
            "node_volume_policy": "mixed_voronoi",
            "require_non_obtuse": True,
        },
    },
    "legacy_cell_reconstructed": {
        "impact_ionization": BV_AVALANCHE_CURRENT_SUPPORT_PROFILES[
            "legacy_cell_reconstructed"
        ],
        "mesh_geometry": {
            "node_volume_policy": "barycentric",
            "require_non_obtuse": False,
        },
    },
}


class TemplateError(ValueError):
    """Raised when a PN2D template or an override is invalid."""


def _parse_override(raw: str) -> tuple[str, Any]:
    if "=" not in raw:
        raise TemplateError(f"override must use NAME=VALUE: {raw!r}")
    name, encoded = raw.split("=", 1)
    name = name.strip()
    if not name:
        raise TemplateError("override name must not be empty")
    try:
        value = json.loads(encoded)
    except json.JSONDecodeError:
        value = encoded
    return name, value


def _is_absolute_path(value: str) -> bool:
    return Path(value).is_absolute() or PureWindowsPath(value).is_absolute()


def _check_parameter(name: str, value: Any, definition: dict[str, Any],
                     allow_absolute_paths: bool) -> None:
    expected = definition.get("type")
    valid = {
        "string": isinstance(value, str),
        "number": isinstance(value, (int, float)) and not isinstance(value, bool),
        "boolean": isinstance(value, bool),
    }.get(expected, False)
    if not valid:
        raise TemplateError(
            f"parameter {name!r} must have type {expected}, got {type(value).__name__}"
        )
    if expected == "number" and not math.isfinite(float(value)):
        raise TemplateError(f"parameter {name!r} must be finite")
    if "choices" in definition and value not in definition["choices"]:
        raise TemplateError(
            f"parameter {name!r} must be one of {definition['choices']}, got {value!r}"
        )
    if definition.get("path") and not allow_absolute_paths and _is_absolute_path(value):
        raise TemplateError(
            f"parameter {name!r} must be relative; use --allow-absolute-paths "
            "only for legacy/external workflows"
        )


def _substitute(value: Any, parameters: dict[str, Any]) -> Any:
    if isinstance(value, dict):
        return {key: _substitute(item, parameters) for key, item in value.items()}
    if isinstance(value, list):
        return [_substitute(item, parameters) for item in value]
    if not isinstance(value, str):
        return value
    match = PLACEHOLDER.fullmatch(value)
    if match:
        name = match.group(1)
        if name not in parameters:
            raise TemplateError(f"unresolved template parameter: {name}")
        return copy.deepcopy(parameters[name])
    return value


def validate_pn2d_config(config: dict[str, Any], template_name: str) -> None:
    required = (
        "simulation_type", "mesh_file", "node_doping_file", "materials_file",
        "output_csv", "contacts", "solver", "sweep",
    )
    missing = [name for name in required if name not in config]
    if missing:
        raise TemplateError(f"rendered config is missing fields: {missing}")
    if config["simulation_type"] != "dc_sweep":
        raise TemplateError("PN2D templates require simulation_type=dc_sweep")

    sweep = config["sweep"]
    for name in ("start", "stop", "step", "initial_step", "min_step", "max_step"):
        value = sweep.get(name)
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise TemplateError(f"sweep.{name} must be numeric")
        if not math.isfinite(float(value)):
            raise TemplateError(f"sweep.{name} must be finite")
    if sweep["step"] == 0:
        raise TemplateError("sweep.step must be non-zero")
    if not 0 < sweep["min_step"] <= sweep["initial_step"] <= sweep["max_step"]:
        raise TemplateError(
            "sweep step bounds must satisfy min_step <= initial_step <= max_step"
        )
    if sweep["stop"] > sweep["start"] and sweep["step"] < 0:
        raise TemplateError("forward sweep requires a positive sweep.step")
    if sweep["stop"] < sweep["start"] and sweep["step"] > 0:
        raise TemplateError("reverse sweep requires a negative sweep.step")

    solver = config["solver"]
    mobility = solver["mobility"]
    impact = solver["impact_ionization"]
    if template_name == "pn2d_iv":
        if sweep.get("mode") != "iv" or impact.get("model") != "none":
            raise TemplateError("pn2d_iv must use IV mode with impact ionization off")
        if mobility.get("model") != "masetti":
            raise TemplateError("pn2d_iv must use the low-field Masetti model")
    elif template_name == "pn2d_bv":
        if sweep.get("mode") != "bv_reverse":
            raise TemplateError("pn2d_bv must use bv_reverse mode")
        if impact.get("model") != "van_overstraeten":
            raise TemplateError("pn2d_bv must enable van_overstraeten")
        if mobility.get("model") != "masetti_field":
            raise TemplateError("pn2d_bv must use masetti_field")
        geometry = config.get("mesh_geometry", {})
        if not isinstance(geometry, dict):
            raise TemplateError("pn2d_bv mesh_geometry must be an object")
        if not isinstance(geometry.get("node_volume_policy"), str):
            raise TemplateError(
                "pn2d_bv mesh_geometry.node_volume_policy must be a string"
            )
        if not isinstance(geometry.get("require_non_obtuse"), bool):
            raise TemplateError(
                "pn2d_bv mesh_geometry.require_non_obtuse must be boolean"
            )
        active_profile = {
            "impact_ionization": {
                name: impact.get(name)
                for name in (
                    "current_approximation",
                    "source_mapping_mode",
                    "cell_reconstructed_midpoint_density",
                )
            },
            "mesh_geometry": {
                "node_volume_policy": geometry.get(
                    "node_volume_policy", "barycentric"
                ),
                "require_non_obtuse": geometry.get(
                    "require_non_obtuse", False
                ),
            },
        }
        if active_profile not in BV_CONFIGURATION_PROFILES.values():
            raise TemplateError(
                "pn2d_bv solver and mesh geometry must match one complete atomic "
                "profile; half-migrated or omitted profile fields are forbidden"
            )
    elif template_name in {
        "bvmethods_nmos_external_resistor",
        "bvmethods_nmos_voltage_to_current",
    }:
        if sweep.get("mode") != "bv_reverse":
            raise TemplateError(f"{template_name} must use bv_reverse mode")
        if impact.get("model") != "van_overstraeten":
            raise TemplateError(f"{template_name} must enable van_overstraeten")
        if impact.get("coupling_mode") != "self_consistent":
            raise TemplateError(f"{template_name} must use self_consistent avalanche")
        if mobility.get("model") != "masetti_field":
            raise TemplateError(f"{template_name} must use masetti_field")
        if template_name == "bvmethods_nmos_external_resistor":
            circuit = sweep.get("external_circuit")
            if not isinstance(circuit, dict) or circuit.get("mode") != "series_resistor":
                raise TemplateError(
                    "bvmethods_nmos_external_resistor requires a series_resistor circuit"
                )
            if "voltage_to_current" in sweep:
                raise TemplateError(
                    "external-resistor and voltage-to-current controls are mutually exclusive"
                )
            if float(circuit.get("resistance_ohm_um", 0.0)) <= 0.0:
                raise TemplateError("external resistance must be positive")
        else:
            control = sweep.get("voltage_to_current")
            if not isinstance(control, dict):
                raise TemplateError(
                    "bvmethods_nmos_voltage_to_current requires voltage_to_current"
                )
            if "external_circuit" in sweep:
                raise TemplateError(
                    "external-resistor and voltage-to-current controls are mutually exclusive"
                )
            if control.get("switch_voltage_V") != sweep.get("stop"):
                raise TemplateError(
                    "voltage-to-current switch voltage must equal the voltage sweep stop"
                )
            points = control.get("current_points_A_per_um")
            if not isinstance(points, list) or not points:
                raise TemplateError("voltage-to-current target list must not be empty")
            if any(not isinstance(value, (int, float)) or value <= 0 for value in points):
                raise TemplateError("voltage-to-current targets must be positive")
    else:
        raise TemplateError(f"unsupported PN2D template: {template_name}")


def _apply_compensated_pn2d(config: dict[str, Any], template_name: str,
                            mobility_basis_overridden: bool) -> None:
    solver = config["solver"]
    solver.pop("handoff")
    solver.update(
        method="newton", state_variables="potential_density",
        primary_variable_coordinates="unit", compensated_population_state=True,
        quasi_fermi_reference="none", quasi_fermi_update_limit_V=0.0,
        population_stopping_rule="native_approx",
        population_update_convergence={"digits": 5, "potential_reference_V": 0.0258},
        line_search=False, tau_n_s=1e-5, tau_p_s=3e-6,
    )
    solver["continuity_row_scaling"]["enabled"] = False
    solver["carrier_row_convergence"]["mode"] = "report"
    solver["global_continuity_closure"]["mode"] = "report"
    if not mobility_basis_overridden:
        solver["mobility"]["doping_concentration_basis"] = "total_impurity"

    reverse = template_name == "pn2d_bv"
    if reverse:
        solver["mobility"]["high_field_gradient_discretization"] = "transport_cell_vector"
        solver["impact_ionization"].update(
            contact_electric_field_fallback=True, source_jacobian="local_ad",
        )
    sweep = config["sweep"]
    span = sweep["stop"] - sweep["start"]
    if span * sweep["step"] < 0:
        direction = "positive" if span > 0 else "negative"
        raise TemplateError(f"{template_name} sweep requires a {direction} sweep.step")
    sweep.update(
        step=span or sweep["step"], step_growth_mode="newton_iterations",
        stop_below_minimum_step=True, preserve_step_after_output_clip=True,
        continuation={"predictor": {"mode": "linear"}},
        initialization={"mode": "poisson_block", "poisson_max_iter": 100},
        stop_on_failure=True, csv_all_contacts=True,
    )
    if span:
        for key in ("initial_step", "min_step", "max_step"):
            sweep[key] *= abs(span)
        sweep["max_retries"] = math.ceil(
            (math.log(sweep["max_step"]) - math.log(sweep["min_step"]))
            / -math.log(sweep["shrink_factor"])
        )
        intervals = 400 if reverse else 200
        sweep["vtk_bias_points"] = [
            sweep["start"] + span * (i / intervals) for i in range(intervals + 1)
        ]
        sweep["vtk_format"] = "vtu"


def render_named_template(
    template_name: str,
    overrides: dict[str, Any] | None = None,
    *,
    allow_absolute_paths: bool = False,
) -> tuple[dict[str, Any], dict[str, Any]]:
    if template_name not in TEMPLATES:
        raise TemplateError(
            f"unknown template {template_name!r}; choose from {sorted(TEMPLATES)}"
        )
    document = json.loads(TEMPLATES[template_name].read_text(encoding="utf-8"))
    if document.get("template_schema") != "vela.simulation-template.v1":
        raise TemplateError("unsupported template schema")
    definitions = document.get("parameters", {})
    supplied = dict(overrides or {})
    unknown = sorted(set(supplied) - set(definitions))
    if unknown:
        raise TemplateError(f"unknown template parameters: {unknown}")

    parameters: dict[str, Any] = {}
    for name, definition in definitions.items():
        if name in supplied:
            value = supplied[name]
        elif "default" in definition:
            value = copy.deepcopy(definition["default"])
        else:
            raise TemplateError(f"required template parameter is missing: {name}")
        _check_parameter(name, value, definition, allow_absolute_paths)
        parameters[name] = value

    config = _substitute(document["config"], parameters)
    if template_name == "pn2d_bv":
        profile = parameters["avalanche_current_support_profile"]
        bundle = BV_CONFIGURATION_PROFILES[profile]
        impact = config["solver"]["impact_ionization"]
        impact.update(copy.deepcopy(bundle["impact_ionization"]))
        config["mesh_geometry"] = copy.deepcopy(bundle["mesh_geometry"])
    if (parameters.get("solver_profile") == "compensated"
            and parameters.get("avalanche_current_support_profile") == "legacy_cell_reconstructed"):
        if "solver_profile" in supplied:
            raise TemplateError("compensated BV requires the element_edge_sg_gss_laux profile")
        parameters["solver_profile"] = "legacy"
    if parameters.get("solver_profile") == "compensated":
        _apply_compensated_pn2d(config, template_name, "mobility_doping_basis" in supplied)
        parameters["mobility_doping_basis"] = config["solver"]["mobility"]["doping_concentration_basis"]
    elif parameters.get("solver_profile") == "legacy" and template_name == "pn2d_iv" and "stop_voltage" not in supplied:
        config["sweep"]["stop"] = parameters["stop_voltage"] = 20.0
    validate_pn2d_config(config, template_name)
    manifest = {
        "generator": "scripts/generate_pn2d_config.py",
        "overrides": supplied,
        "parameters": parameters,
        "template": template_name,
        "template_schema": document["template_schema"],
        "template_version": document["version"],
    }
    if template_name == "pn2d_bv":
        manifest["resolved_profile"] = copy.deepcopy(
            BV_CONFIGURATION_PROFILES[parameters["avalanche_current_support_profile"]]
        )
    return config, manifest


def write_rendered_config(
    template_name: str,
    output: Path,
    overrides: dict[str, Any] | None = None,
    *,
    manifest_path: Path | None = None,
    allow_absolute_paths: bool = False,
) -> tuple[Path, Path]:
    config, manifest = render_named_template(
        template_name, overrides, allow_absolute_paths=allow_absolute_paths
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if manifest_path is None:
        manifest_path = output.with_suffix(".manifest.json")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return output, manifest_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", choices=sorted(TEMPLATES), required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument(
        "--set", action="append", default=[], metavar="NAME=VALUE",
        help="Override a declared template parameter; VALUE accepts JSON syntax.",
    )
    parser.add_argument(
        "--allow-absolute-paths", action="store_true",
        help="Permit absolute artifact paths for legacy/external workflows.",
    )
    args = parser.parse_args()
    try:
        overrides = dict(_parse_override(raw) for raw in args.set)
        write_rendered_config(
            args.template,
            args.output,
            overrides,
            manifest_path=args.manifest,
            allow_absolute_paths=args.allow_absolute_paths,
        )
    except (OSError, TemplateError, json.JSONDecodeError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
