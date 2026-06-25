#!/usr/bin/env python3
"""Audit and inject explicit PN2D BV impact-ionization coefficients.

The tool distinguishes two cases that are easy to conflate:

* Sentaurus command files can select ``Avalanche(VanOverstraeten)``.
* ``models.par`` may or may not contain local numeric coefficients for that model.

When local coefficients are present in a recognizable VanOverstraeten block, they are
converted to Vela SI JSON keys.  When the current Sentaurus 2018 silicon
``models.par`` has no such local block, ``--fallback vela-defaults`` writes Vela's
Sentaurus-2018-compatible built-in coefficients explicitly and records that they did
not come from ``models.par``.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_MODELS_PAR = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "models.par"
DEFAULT_CONFIG = REPO / "build" / "pn2d_bv_gradqf_full20" / "simulation_bv_genius_cell_gradient_full20.json"
DEFAULT_REPORT = (
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reports" /
    "impact_parameter_audit" / "impact_parameter_audit_summary.json"
)
DEFAULT_WRITE_CONFIG = DEFAULT_CONFIG.with_name("simulation_bv_genius_cell_gradient_full20_impact_explicit.json")

SENTAURUS2018_VELA_DEFAULTS: dict[str, float | str] = {
    "model": "van_overstraeten",
    "electron_A_m_inv": 7.03e7,
    "electron_B_V_m": 1.231e8,
    "hole_A_m_inv": 1.582e8,
    "hole_B_V_m": 2.036e8,
    "electron_a_low_m_inv": 7.03e7,
    "electron_a_high_m_inv": 7.03e7,
    "electron_b_low_V_m": 1.231e8,
    "electron_b_high_V_m": 1.231e8,
    "hole_a_low_m_inv": 1.582e8,
    "hole_a_high_m_inv": 6.71e7,
    "hole_b_low_V_m": 2.036e8,
    "hole_b_high_V_m": 1.693e8,
    "switch_field_V_m": 4.0e7,
    "phonon_energy_eV": 0.063,
    "reference_temperature_K": 300.0,
    "temperature_K": 300.0,
    "minimum_field_V_m": 0.0,
}

INVERSE_LENGTH_KEYS = {
    "electron_A_m_inv",
    "hole_A_m_inv",
    "electron_a_low_m_inv",
    "electron_a_high_m_inv",
    "hole_a_low_m_inv",
    "hole_a_high_m_inv",
}
ELECTRIC_FIELD_KEYS = {
    "electron_B_V_m",
    "hole_B_V_m",
    "electron_b_low_V_m",
    "electron_b_high_V_m",
    "hole_b_low_V_m",
    "hole_b_high_V_m",
    "switch_field_V_m",
    "minimum_field_V_m",
}

PARAM_ALIASES: dict[str, str] = {
    "electron_A": "electron_A_m_inv",
    "electron_a": "electron_A_m_inv",
    "e_A": "electron_A_m_inv",
    "n_A": "electron_A_m_inv",
    "electron_B": "electron_B_V_m",
    "electron_b": "electron_B_V_m",
    "e_B": "electron_B_V_m",
    "n_B": "electron_B_V_m",
    "hole_A": "hole_A_m_inv",
    "hole_a": "hole_A_m_inv",
    "h_A": "hole_A_m_inv",
    "p_A": "hole_A_m_inv",
    "hole_B": "hole_B_V_m",
    "hole_b": "hole_B_V_m",
    "h_B": "hole_B_V_m",
    "p_B": "hole_B_V_m",
    "electron_a_low": "electron_a_low_m_inv",
    "electron_a_high": "electron_a_high_m_inv",
    "electron_b_low": "electron_b_low_V_m",
    "electron_b_high": "electron_b_high_V_m",
    "hole_a_low": "hole_a_low_m_inv",
    "hole_a_high": "hole_a_high_m_inv",
    "hole_b_low": "hole_b_low_V_m",
    "hole_b_high": "hole_b_high_V_m",
    "e_Alow": "electron_a_low_m_inv",
    "e_Ahigh": "electron_a_high_m_inv",
    "e_Blow": "electron_b_low_V_m",
    "e_Bhigh": "electron_b_high_V_m",
    "h_Alow": "hole_a_low_m_inv",
    "h_Ahigh": "hole_a_high_m_inv",
    "h_Blow": "hole_b_low_V_m",
    "h_Bhigh": "hole_b_high_V_m",
    "n_Alow": "electron_a_low_m_inv",
    "n_Ahigh": "electron_a_high_m_inv",
    "n_Blow": "electron_b_low_V_m",
    "n_Bhigh": "electron_b_high_V_m",
    "p_Alow": "hole_a_low_m_inv",
    "p_Ahigh": "hole_a_high_m_inv",
    "p_Blow": "hole_b_low_V_m",
    "p_Bhigh": "hole_b_high_V_m",
    "switch_field": "switch_field_V_m",
    "E0": "switch_field_V_m",
    "phonon_energy": "phonon_energy_eV",
    "reference_temperature": "reference_temperature_K",
    "temperature": "temperature_K",
}
for key in list(INVERSE_LENGTH_KEYS | ELECTRIC_FIELD_KEYS):
    PARAM_ALIASES[key] = key
PARAM_ALIASES.update({
    key.replace("_m_inv", "_cm_inv"): key for key in INVERSE_LENGTH_KEYS
})
PARAM_ALIASES.update({
    key.replace("_V_m", "_V_cm"): key for key in ELECTRIC_FIELD_KEYS
})

PARAM_RE = re.compile(
    r"^\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^#\n]+?)(?:\s*#\s*\[([^\]]+)\])?\s*$"
)
NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")
SECTION_START_RE = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)(?::)?(?:\s*\{)?\s*(?:#.*)?$")


def parse_number(value: str) -> float | None:
    match = NUMBER_RE.search(value.replace("D", "e").replace("d", "e"))
    if not match:
        return None
    parsed = float(match.group(0))
    return parsed if math.isfinite(parsed) else None


def canonical_unit(raw: str | None, key: str, source_name: str) -> str:
    if raw:
        unit = raw.strip().replace(" ", "")
        if unit in {"1/cm", "cm^-1", "cm-1"}:
            return "1/cm"
        if unit in {"1/m", "m^-1", "m-1"}:
            return "1/m"
        if unit in {"V/cm", "Vcm^-1", "Vcm-1"}:
            return "V/cm"
        if unit in {"V/m", "Vm^-1", "Vm-1"}:
            return "V/m"
        if unit in {"eV", "K", "1"}:
            return unit
    if source_name.endswith("_cm_inv"):
        return "1/cm"
    if source_name.endswith("_V_cm"):
        return "V/cm"
    if key in INVERSE_LENGTH_KEYS:
        return "1/m"
    if key in ELECTRIC_FIELD_KEYS:
        return "V/m"
    return "1"


def convert_to_vela_si(value: float, key: str, unit: str) -> float:
    if key in INVERSE_LENGTH_KEYS:
        return value * 100.0 if unit == "1/cm" else value
    if key in ELECTRIC_FIELD_KEYS:
        return value * 100.0 if unit == "V/cm" else value
    return value


def parse_top_level_sections(text: str) -> dict[str, str]:
    sections: dict[str, str] = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if not stripped or stripped.startswith("*"):
            i += 1
            continue
        match = SECTION_START_RE.match(line)
        if not match:
            i += 1
            continue
        name = match.group(1)
        open_on_line = "{" in line
        if not open_on_line:
            j = i + 1
            while j < len(lines) and (not lines[j].strip() or lines[j].strip().startswith("*")):
                j += 1
            if j >= len(lines) or not lines[j].lstrip().startswith("{"):
                i += 1
                continue
            i = j
            line = lines[i]
        depth = line.count("{") - line.count("}")
        body: list[str] = []
        i += 1
        while i < len(lines) and depth > 0:
            current = lines[i]
            depth += current.count("{") - current.count("}")
            if depth >= 0 and not current.lstrip().startswith("}"):
                body.append(current)
            i += 1
        sections[name] = "\n".join(body)
    return sections


def parse_section_parameters(body: str) -> dict[str, tuple[float, str, str]]:
    params: dict[str, tuple[float, str, str]] = {}
    for line in body.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("*"):
            continue
        match = PARAM_RE.match(line)
        if not match:
            continue
        raw_name, raw_value, raw_unit = match.groups()
        canonical = PARAM_ALIASES.get(raw_name)
        if canonical is None:
            continue
        parsed = parse_number(raw_value)
        if parsed is None:
            continue
        unit = canonical_unit(raw_unit, canonical, raw_name)
        params[canonical] = (convert_to_vela_si(parsed, canonical, unit), unit, raw_name)
    return params


def extract_impact_parameters(models_par: Path) -> tuple[dict[str, float], str | None, list[str], list[str]]:
    text = models_par.read_text(encoding="utf-8", errors="replace")
    sections = parse_top_level_sections(text)
    warnings: list[str] = []
    observed = sorted(sections.keys())
    candidates: list[tuple[str, dict[str, tuple[float, str, str]]]] = []
    for name, body in sections.items():
        lower = name.lower()
        if "vanoverstraeten" not in lower and "impact" not in lower and lower != "avalanche":
            continue
        if lower in {"avalanchefactors", "ionization"}:
            continue
        params = parse_section_parameters(body)
        if params:
            candidates.append((name, params))
    if not candidates:
        warnings.append("no_local_impact_ionization_coefficients_in_models_par")
        return {}, None, warnings, observed
    name, params = max(candidates, key=lambda item: len(item[1]))
    extracted = {key: value for key, (value, _unit, _raw_name) in params.items()}
    return extracted, f"models_par:{name}", warnings, observed


def load_config(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def uses_unit_scaling(config: dict[str, Any]) -> bool:
    scaling = config.get("scaling")
    return isinstance(scaling, dict) and scaling.get("mode") == "unit_scaling"


def params_for_config_units(params: dict[str, float | str],
                            config: dict[str, Any]) -> tuple[dict[str, float | str], str]:
    if not uses_unit_scaling(config):
        return dict(params), "legacy_si"
    converted: dict[str, float | str] = {}
    for key, value in params.items():
        if isinstance(value, (int, float)) and key in (INVERSE_LENGTH_KEYS | ELECTRIC_FIELD_KEYS):
            converted[key] = float(value) / 100.0
        else:
            converted[key] = value
    return converted, "unit_scaling_input_cm"


def update_config(config: dict[str, Any], params: dict[str, float | str], source: str) -> dict[str, Any]:
    solver = config.setdefault("solver", {})
    impact = solver.setdefault("impact_ionization", {})
    if isinstance(impact, str):
        impact = {"model": impact}
        solver["impact_ionization"] = impact
    if not isinstance(impact, dict):
        raise TypeError("solver.impact_ionization must be a string or object")
    for key, value in params.items():
        if key == "model" and "model" in impact:
            continue
        impact[key] = value
    impact["_parameter_source"] = source
    return config


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_report(args: argparse.Namespace,
                 extracted: dict[str, float],
                 source: str,
                 warnings: list[str],
                 observed_sections: list[str],
                 written_config: Path | None,
                 written_parameter_units: str | None) -> dict[str, Any]:
    return {
        "models_par": str(args.models_par.resolve()),
        "config": str(args.config.resolve()) if args.config else None,
        "written_config": str(written_config.resolve()) if written_config else None,
        "impact_parameter_source": source,
        "written_parameter_units": written_parameter_units,
        "extracted_from_models_par": bool(extracted),
        "extracted_keys": sorted(extracted.keys()),
        "warnings": warnings,
        "observed_relevant_sections": [
            name for name in observed_sections
            if name in {"AvalancheFactors", "Ionization"} or "avalanche" in name.lower()
            or "impact" in name.lower() or "vanoverstraeten" in name.lower()
        ],
    }

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--models-par", type=Path, default=DEFAULT_MODELS_PAR)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--write-config", type=Path, default=DEFAULT_WRITE_CONFIG)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument(
        "--fallback",
        choices=["none", "vela-defaults"],
        default="none",
        help="Parameter source to use when models.par has no local impact-ionization coefficients.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    extracted, source, warnings, observed_sections = extract_impact_parameters(args.models_par)
    params: dict[str, float | str]
    if extracted:
        params = dict(extracted)
        source = source or "models_par"
    elif args.fallback == "vela-defaults":
        params = dict(SENTAURUS2018_VELA_DEFAULTS)
        source = "vela_default_sentaurus2018_compatible"
    else:
        source = "none"
        params = {}
        warnings.append("no_parameters_written_without_fallback")

    written_config: Path | None = None
    written_parameter_units: str | None = None
    if args.write_config:
        if not params:
            print("No impact-ionization parameters available to write.", file=sys.stderr)
            return 2
        config = load_config(args.config)
        writable_params, written_parameter_units = params_for_config_units(params, config)
        update_config(config, writable_params, source)
        write_json(args.write_config, config)
        written_config = args.write_config

    report = build_report(
        args,
        extracted,
        source,
        warnings,
        observed_sections,
        written_config,
        written_parameter_units,
    )
    if args.report:
        write_json(args.report, report)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0

if __name__ == "__main__":
    raise SystemExit(main())