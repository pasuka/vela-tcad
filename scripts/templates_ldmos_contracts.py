#!/usr/bin/env python3
"""Validate and render Templates/LDMOS phase 0/1 governance contracts."""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
SCHEMA_DIR = REPO / "schemas"
BENCHMARK = "sentaurus_t2022_03_sp2_templates_ldmos"
SCHEMA_FILES = {
    "vela.templates_ldmos.validation_summary.v1":
        "vela.templates_ldmos.validation_summary.v1.schema.json",
    "vela.templates_ldmos.known_difference_ledger.v1":
        "vela.templates_ldmos.known_difference_ledger.v1.schema.json",
    "vela.templates_ldmos.threshold_freeze.v1":
        "vela.templates_ldmos.threshold_freeze.v1.schema.json",
    "vela.templates_ldmos.budget_freeze.v1":
        "vela.templates_ldmos.budget_freeze.v1.schema.json",
    "vela.templates_ldmos.phase01_manifest.v1":
        "vela.templates_ldmos.phase01_manifest.v1.schema.json",
    "vela.templates_ldmos.materials.v1":
        "vela.templates_ldmos.materials.v1.schema.json",
    "vela.templates_ldmos.physics_contract.v1":
        "vela.templates_ldmos.physics_contract.v1.schema.json",
    "vela.templates_ldmos.discretization_contract.v1":
        "vela.templates_ldmos.discretization_contract.v1.schema.json",
    "vela.templates_ldmos.discretization_contract.v2":
        "vela.templates_ldmos.discretization_contract.v2.schema.json",
    "vela.templates_ldmos.restart_qualification.v1":
        "vela.templates_ldmos.restart_qualification.v1.schema.json",
    "vela.templates_ldmos_g3_interface_charge_volume_audit.v1":
        "vela.templates_ldmos_g3_interface_charge_volume_audit.v1.schema.json",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return payload


def _resolve_ref(root: dict[str, Any], ref: str) -> dict[str, Any]:
    if not ref.startswith("#/"):
        raise ValueError(f"unsupported external schema reference: {ref}")
    value: Any = root
    for token in ref[2:].split("/"):
        value = value[token.replace("~1", "/").replace("~0", "~")]
    if not isinstance(value, dict):
        raise ValueError(f"schema reference is not an object: {ref}")
    return value


def _type_matches(value: Any, expected: str) -> bool:
    if expected == "null":
        return value is None
    if expected == "boolean":
        return isinstance(value, bool)
    if expected == "integer":
        return isinstance(value, int) and not isinstance(value, bool)
    if expected == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected == "string":
        return isinstance(value, str)
    if expected == "array":
        return isinstance(value, list)
    if expected == "object":
        return isinstance(value, dict)
    raise ValueError(f"unsupported schema type: {expected}")


def validate_instance(value: Any,
                      schema: dict[str, Any],
                      *,
                      root: dict[str, Any] | None = None,
                      path: str = "$") -> None:
    root = root or schema
    if "$ref" in schema:
        validate_instance(value, _resolve_ref(root, schema["$ref"]), root=root, path=path)
        return
    if "const" in schema and value != schema["const"]:
        raise ValueError(f"{path}: expected constant {schema['const']!r}")
    if "enum" in schema and value not in schema["enum"]:
        raise ValueError(f"{path}: value {value!r} is outside enum")
    expected = schema.get("type")
    if expected is not None:
        choices = expected if isinstance(expected, list) else [expected]
        if not any(_type_matches(value, choice) for choice in choices):
            raise ValueError(f"{path}: expected type {choices}, got {type(value).__name__}")
    if isinstance(value, dict):
        properties = schema.get("properties", {})
        missing = [name for name in schema.get("required", []) if name not in value]
        if missing:
            raise ValueError(f"{path}: missing required keys {missing}")
        for name, child in value.items():
            if name in properties:
                validate_instance(child, properties[name], root=root, path=f"{path}.{name}")
            elif schema.get("additionalProperties") is False:
                raise ValueError(f"{path}: unexpected key {name!r}")
            elif isinstance(schema.get("additionalProperties"), dict):
                validate_instance(
                    child, schema["additionalProperties"], root=root, path=f"{path}.{name}")
    if isinstance(value, list):
        if len(value) < int(schema.get("minItems", 0)):
            raise ValueError(f"{path}: too few items")
        if "maxItems" in schema and len(value) > int(schema["maxItems"]):
            raise ValueError(f"{path}: too many items")
        if schema.get("uniqueItems") and len({json.dumps(item, sort_keys=True)
                                               for item in value}) != len(value):
            raise ValueError(f"{path}: array items must be unique")
        item_schema = schema.get("items")
        if isinstance(item_schema, dict):
            for index, item in enumerate(value):
                validate_instance(item, item_schema, root=root, path=f"{path}[{index}]")
    if isinstance(value, str):
        if len(value) < int(schema.get("minLength", 0)):
            raise ValueError(f"{path}: string is too short")
        if "pattern" in schema and re.fullmatch(schema["pattern"], value) is None:
            raise ValueError(f"{path}: value does not match {schema['pattern']!r}")
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            raise ValueError(f"{path}: value is below minimum {schema['minimum']}")
        if "maximum" in schema and value > schema["maximum"]:
            raise ValueError(f"{path}: value is above maximum {schema['maximum']}")
        if "exclusiveMinimum" in schema and value <= schema["exclusiveMinimum"]:
            raise ValueError(
                f"{path}: value must exceed {schema['exclusiveMinimum']}")
        if "exclusiveMaximum" in schema and value >= schema["exclusiveMaximum"]:
            raise ValueError(
                f"{path}: value must be below {schema['exclusiveMaximum']}")


def validate_document(document: dict[str, Any], schema_dir: Path = SCHEMA_DIR) -> None:
    schema_name = document.get("schema")
    if schema_name not in SCHEMA_FILES:
        raise ValueError(f"unsupported Templates/LDMOS schema: {schema_name!r}")
    schema = read_json(schema_dir / SCHEMA_FILES[str(schema_name)])
    validate_instance(document, schema)


def canonical_round_trip(document: dict[str, Any]) -> dict[str, Any]:
    """Validate a contract and prove canonical JSON serialization is lossless."""
    validate_document(document)
    restored = json.loads(json.dumps(
        document, sort_keys=True, separators=(",", ":"), allow_nan=False))
    if restored != document:
        raise ValueError("contract changed during canonical JSON round-trip")
    validate_document(restored)
    return restored


def migrate_materials_v0(document: dict[str, Any],
                         source_unit_system: str) -> dict[str, Any]:
    """Migrate an explicitly identified legacy material file to v1 units."""
    if source_unit_system not in {"legacy_si", "tcad_internal"}:
        raise ValueError(
            "source_unit_system must be 'legacy_si' or 'tcad_internal'")
    if set(document) != {"materials"} or not isinstance(document["materials"], list):
        raise ValueError(
            "legacy material migration accepts only an object containing a materials array")
    allowed = {
        "name", "eps_r", "ni", "mun", "mup", "bandgap_eV",
        "electron_affinity_eV", "Nc_m3", "Nv_m3", "temperature_K",
        "electron_quantum_gamma", "electron_quantum_dos_mass_ratio",
        "electron_quantum_coefficient_mass_ratio",
        "thermal_conductivity_W_per_m_K", "specific_heat_J_per_kg_K",
        "mass_density_kg_per_m3",
    }
    key_map = {
        "ni": "intrinsic_carrier_density_cm3",
        "mun": "electron_mobility_cm2_per_V_s",
        "mup": "hole_mobility_cm2_per_V_s",
        "Nc_m3": "conduction_band_density_of_states_cm3",
        "Nv_m3": "valence_band_density_of_states_cm3",
    }
    migrated_materials = []
    for index, source in enumerate(document["materials"]):
        if not isinstance(source, dict):
            raise ValueError(f"legacy material {index} must be an object")
        unexpected = sorted(set(source) - allowed)
        if unexpected:
            raise ValueError(f"legacy material {index} has unknown keys {unexpected}")
        name = source.get("name")
        if not isinstance(name, str) or not name:
            raise ValueError(f"legacy material {index} requires a non-empty name")
        target: dict[str, Any] = {
            "name": name,
            "role": "dielectric" if name in {"SiO2", "Oxide"}
                    else "transport_semiconductor",
            "provenance": f"migrated_from_{source_unit_system}_v0",
        }
        for key, value in source.items():
            if key == "name":
                continue
            target_key = key_map.get(key, key)
            converted = value
            if source_unit_system == "legacy_si":
                if key in {"ni", "Nc_m3", "Nv_m3"}:
                    converted = value / 1.0e6
                elif key in {"mun", "mup"}:
                    converted = value / 1.0e-4
            target[target_key] = converted
        migrated_materials.append(target)
    migrated = {
        "schema": "vela.templates_ldmos.materials.v1",
        "benchmark": BENCHMARK,
        "revision": 1,
        "unit_system": {
            "concentration": "cm^-3",
            "mobility": "cm^2/(V*s)",
            "energy": "eV",
            "temperature": "K",
            "thermal_conductivity": "W/(m*K)",
            "specific_heat": "J/(kg*K)",
            "mass_density": "kg/m^3",
        },
        "materials": migrated_materials,
    }
    return canonical_round_trip(migrated)


def migrate_discretization_draft(document: dict[str, Any]) -> dict[str, Any]:
    expected = {
        "schema", "profile_name", "applicable_mesh", "current_support",
        "control_volume", "field_recovery", "volume_source_mapping",
        "contact_edge_integration", "obtuse_policy", "require_non_obtuse",
        "non_delaunay_policy", "avalanche_profile", "forbidden_inference",
        "physics_use_authorized", "status",
    }
    unexpected = sorted(set(document) - expected)
    if unexpected:
        raise ValueError(f"discretization draft has unknown keys {unexpected}")
    if document.get("schema") != \
            "vela.templates_ldmos.discretization_contract.v1-draft-unvalidated":
        raise ValueError("input is not the recognized stage-1 discretization draft")
    migrated = {
        "schema": "vela.templates_ldmos.discretization_contract.v2",
        "benchmark": BENCHMARK,
        "revision": 2,
        "profile_name": "templates_ldmos_exact_topology_phase_a_classical_v2",
        "applicable_mesh": document["applicable_mesh"],
        "scope": ["stage_1_5", "L2", "L3"],
        "current_support":
            "scharfetter_gummel_edge_flux_external_averagebox_couple",
        "control_volume": "mesh_barycentric",
        "poisson_charge_volume": "transport_material_local_barycentric",
        "continuity_source_volume": "mesh_barycentric",
        "poisson_edge_support":
            "mesh_default_edge_average_epsilon_cotangent_fallback",
        "field_recovery": document["field_recovery"],
        "volume_source_mapping": "cell_reconstructed_unchanged",
        "contact_edge_integration": document["contact_edge_integration"],
        "obtuse_policy": document["obtuse_policy"],
        "non_delaunay_policy": "qualified_for_phase_a_classical_only",
        "avalanche_profile": "not_authorized_in_phase_a",
        "atomic_constraints": [
            "carrier transport uses the exact 16237-edge archived AverageBox profile only on the qualified mesh",
            "material-local volume changes only Poisson mobile-carrier and ionized-dopant charge residual and Jacobian terms",
            "continuity source volumes, fixed/interface charge, BTBT, SRH/Auger, impact ionization and stored charge remain on their existing mappings",
            "avalanche and thermal source mappings require a separately versioned phase-B profile",
        ],
        "forbidden_inference": document["forbidden_inference"],
        "physics_use_authorized": True,
        "status": "qualified_for_phase_a_classical_stage3",
    }
    return canonical_round_trip(migrated)


def draft_governance_contracts(vela_commit: str,
                               oracle_manifest_sha256: str | None = None
                               ) -> dict[str, dict[str, Any]]:
    approval = {
        "status": "draft",
        "benchmark_owner": None,
        "independent_reviewer": None,
        "approved_at": None,
    }
    ledger = {
        "schema": "vela.templates_ldmos.known_difference_ledger.v1",
        "benchmark": BENCHMARK,
        "entries": [],
    }
    thresholds = {
        "schema": "vela.templates_ldmos.threshold_freeze.v1",
        "benchmark": BENCHMARK,
        "revision": 1,
        "scope": "phase_a_decisions",
        "thresholds": [
            {
                "id": "hqp_relative_budget_fraction",
                "metric": "hQP delta divided by applicable final tolerance",
                "operator": "<",
                "value": 0.2,
                "unit": "fraction",
                "gate_type": "decision",
                "applies_to": ["G0-G1", "D2-D3"],
            },
            {
                "id": "hqp_vth_shift",
                "metric": "absolute fixed-current threshold-voltage shift",
                "operator": "<",
                "value": 10.0,
                "unit": "mV",
                "gate_type": "decision",
                "applies_to": ["G0-G1"],
            },
            {
                "id": "hqp_strong_inversion_current_delta",
                "metric": "absolute strong-inversion current change",
                "operator": "<",
                "value": 2.0,
                "unit": "percent",
                "gate_type": "decision",
                "applies_to": ["G0-G1"],
            },
            {
                "id": "hrec_relative_budget_fraction",
                "metric": "hRecVelocity delta divided by applicable final tolerance",
                "operator": "<",
                "value": 0.2,
                "unit": "fraction",
                "gate_type": "decision",
                "applies_to": ["D1-D2"],
            },
        ],
        "evidence": [],
        "approval": approval.copy(),
    }
    budget = {
        "schema": "vela.templates_ldmos.budget_freeze.v1",
        "benchmark": BENCHMARK,
        "revision": 1,
        "basis": {
            "vela_commit": vela_commit,
            "oracle_manifest_sha256": oracle_manifest_sha256,
            "probe": {"status": "not_run"},
        },
        "scenarios": [
            {"name": name, "wall_clock_hours": None, "peak_memory_gib": None,
             "storage_gib": None, "max_concurrency": None}
            for name in ("best", "base", "worst")
        ],
        "scope_budgets": [
            {
                "scope": scope,
                "status": "pending_evidence",
                "best_wall_clock_hours": None,
                "base_wall_clock_hours": None,
                "worst_wall_clock_hours": None,
                "estimated_bias_points": None,
            }
            for scope in ("phase_a", "phase_a_plus", "phase_b")
        ],
        "matrix_policy": "Pending exact-mesh phase-1 cost probe.",
        "approval": approval.copy(),
    }
    return {
        "known_difference_ledger.json": ledger,
        "threshold_freeze.json": thresholds,
        "budget_freeze.json": budget,
    }


def render_summary(document: dict[str, Any]) -> str:
    validate_document(document)
    lines = [
        "# Templates/LDMOS validation summary",
        "",
        f"- Status: `{document['status']}`",
        f"- Highest validated level: `{document['highest_level']}`",
        f"- Generated: `{document['generated_at']}`",
        "",
        "## Gates",
        "",
        "| Gate | Status | Summary |",
        "| --- | --- | --- |",
    ]
    for gate in document["gates"]:
        summary = str(gate["summary"]).replace("|", "\\|")
        lines.append(f"| `{gate['id']}` | `{gate['status']}` | {summary} |")
    lines.extend(["", "## Limitations", ""])
    limitations = document.get("limitations", [])
    lines.extend(f"- {item}" for item in limitations)
    if not limitations:
        lines.append("- None recorded.")
    lines.extend(["", "## Evidence", ""])
    for item in document["evidence"]:
        lines.append(f"- `{item['path']}` — `{item['sha256']}`")
    if not document["evidence"]:
        lines.append("- No evidence sealed yet.")
    return "\n".join(lines) + "\n"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    validate_document(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("documents", type=Path, nargs="+")
    render_parser = subparsers.add_parser("render")
    render_parser.add_argument("summary", type=Path)
    render_parser.add_argument("output", type=Path)
    migrate_parser = subparsers.add_parser("migrate")
    migrate_parser.add_argument("kind", choices=("materials_v0", "discretization_draft"))
    migrate_parser.add_argument("source", type=Path)
    migrate_parser.add_argument("output", type=Path)
    migrate_parser.add_argument(
        "--source-unit-system", choices=("legacy_si", "tcad_internal"))
    args = parser.parse_args(argv)
    if args.command == "validate":
        for path in args.documents:
            validate_document(read_json(path))
            print(f"validated {path}")
        return 0
    if args.command == "render":
        document = read_json(args.summary)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(render_summary(document), encoding="utf-8")
        print(args.output)
        return 0
    source = read_json(args.source)
    if args.kind == "materials_v0":
        if args.source_unit_system is None:
            parser.error("materials_v0 migration requires --source-unit-system")
        document = migrate_materials_v0(source, args.source_unit_system)
    else:
        document = migrate_discretization_draft(source)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
