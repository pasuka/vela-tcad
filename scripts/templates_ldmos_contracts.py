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


def validate_document(document: dict[str, Any], schema_dir: Path = SCHEMA_DIR) -> None:
    schema_name = document.get("schema")
    if schema_name not in SCHEMA_FILES:
        raise ValueError(f"unsupported Templates/LDMOS schema: {schema_name!r}")
    schema = read_json(schema_dir / SCHEMA_FILES[str(schema_name)])
    validate_instance(document, schema)


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
    args = parser.parse_args(argv)
    if args.command == "validate":
        for path in args.documents:
            validate_document(read_json(path))
            print(f"validated {path}")
        return 0
    document = read_json(args.summary)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render_summary(document), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
