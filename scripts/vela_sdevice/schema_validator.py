"""Minimal JSON Schema (draft-07 subset) validator with no third-party dependency.

The ``vela device`` orchestrator ships with a bundled CPython and no site
packages, so it cannot rely on ``jsonschema`` to check sidecars, stage plans or
artifact manifests.  This module implements the subset of draft-07 that the
Vela contract schemas use:

``type`` (incl. lists), ``enum``, ``const``, ``required``, ``properties``,
``additionalProperties`` (bool or schema), ``items`` (single schema),
``minItems``/``maxItems``, ``minimum``/``maximum``/``exclusiveMinimum``/
``exclusiveMaximum`` (draft-07 numeric form), ``minLength``/``maxLength``,
``pattern``, ``allOf``/``anyOf``/``oneOf``/``not`` and local ``$ref``
(``#/definitions/...`` or any ``#/...`` JSON pointer).  ``default``,
``description``, ``title``, ``$schema``, ``$id`` and ``format`` are ignored.

Unsupported keywords raise :class:`UnsupportedSchemaKeyword` so that a schema
silently accepted here can never be stricter under a full validator without
the test-suite noticing.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Iterator, List, Sequence

__all__ = [
    "SchemaError",
    "UnsupportedSchemaKeyword",
    "ValidationError",
    "load_schema",
    "iter_errors",
    "validate",
    "is_valid",
]

_IGNORED_KEYWORDS = {
    "$schema", "$id", "$comment", "title", "description", "default", "examples", "format",
    "definitions", "readOnly", "writeOnly",
}
_HANDLED_KEYWORDS = {
    "type", "enum", "const", "required", "properties", "additionalProperties", "items",
    "minItems", "maxItems", "minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum",
    "minLength", "maxLength", "pattern", "allOf", "anyOf", "oneOf", "not", "$ref",
}


class SchemaError(ValueError):
    """The schema itself is malformed."""


class UnsupportedSchemaKeyword(SchemaError):
    """The schema uses a keyword this validator does not implement."""


class ValidationError(ValueError):
    def __init__(self, path: str, message: str) -> None:
        super().__init__(f"{path or '$'}: {message}")
        self.path = path or "$"
        self.message = message


def load_schema(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _type_matches(value: Any, type_name: str) -> bool:
    if type_name == "object":
        return isinstance(value, dict)
    if type_name == "array":
        return isinstance(value, list)
    if type_name == "string":
        return isinstance(value, str)
    if type_name == "integer":
        return (isinstance(value, int) and not isinstance(value, bool)) or (
            isinstance(value, float) and value.is_integer())
    if type_name == "number":
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if type_name == "boolean":
        return isinstance(value, bool)
    if type_name == "null":
        return value is None
    raise SchemaError(f"unknown type name {type_name!r}")


def _resolve_ref(root: dict, ref: str) -> Any:
    if not ref.startswith("#"):
        raise UnsupportedSchemaKeyword(f"only local $ref is supported, got {ref!r}")
    node: Any = root
    for raw in ref[1:].split("/"):
        if raw == "":
            continue
        token = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(node, list):
            node = node[int(token)]
        elif isinstance(node, dict) and token in node:
            node = node[token]
        else:
            raise SchemaError(f"unresolvable $ref {ref!r}")
    return node


def _equal(a: Any, b: Any) -> bool:
    # JSON equality: 1 == 1.0 but True != 1.
    if isinstance(a, bool) or isinstance(b, bool):
        return isinstance(a, bool) and isinstance(b, bool) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return a == b
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_equal(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_equal(a[k], b[k]) for k in a)
    return type(a) is type(b) and a == b


def _join(path: str, token: Any) -> str:
    return f"{path}/{token}" if path else f"/{token}"


def _iter(root: dict, schema: Any, value: Any, path: str) -> Iterator[ValidationError]:
    if schema is True:
        return
    if schema is False:
        yield ValidationError(path, "schema forbids any value")
        return
    if not isinstance(schema, dict):
        raise SchemaError(f"schema node at {path or '$'} must be an object or boolean")

    for key in schema:
        if key not in _HANDLED_KEYWORDS and key not in _IGNORED_KEYWORDS:
            raise UnsupportedSchemaKeyword(f"keyword {key!r} at {path or '$'} is not supported")

    if "$ref" in schema:
        yield from _iter(root, _resolve_ref(root, schema["$ref"]), value, path)

    if "type" in schema:
        types = schema["type"]
        if isinstance(types, str):
            types = [types]
        if not any(_type_matches(value, t) for t in types):
            yield ValidationError(path, f"expected type {'/'.join(types)}, got {type(value).__name__}")
            return

    if "enum" in schema and not any(_equal(value, option) for option in schema["enum"]):
        yield ValidationError(path, f"value {json.dumps(value)} is not one of {json.dumps(schema['enum'])}")
    if "const" in schema and not _equal(value, schema["const"]):
        yield ValidationError(path, f"value {json.dumps(value)} is not the constant {json.dumps(schema['const'])}")

    if isinstance(value, dict):
        for key in schema.get("required", []):
            if key not in value:
                yield ValidationError(path, f"missing required property {key!r}")
        properties = schema.get("properties", {})
        for key, sub in properties.items():
            if key in value:
                yield from _iter(root, sub, value[key], _join(path, key))
        additional = schema.get("additionalProperties", True)
        for key in value:
            if key in properties:
                continue
            if additional is False:
                yield ValidationError(path, f"additional property {key!r} is not allowed")
            elif additional is not True:
                yield from _iter(root, additional, value[key], _join(path, key))

    if isinstance(value, list):
        if "minItems" in schema and len(value) < schema["minItems"]:
            yield ValidationError(path, f"expected at least {schema['minItems']} items, got {len(value)}")
        if "maxItems" in schema and len(value) > schema["maxItems"]:
            yield ValidationError(path, f"expected at most {schema['maxItems']} items, got {len(value)}")
        items = schema.get("items")
        if items is not None:
            if isinstance(items, list):
                raise UnsupportedSchemaKeyword("tuple-form 'items' is not supported")
            for index, element in enumerate(value):
                yield from _iter(root, items, element, _join(path, index))

    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if "minimum" in schema and value < schema["minimum"]:
            yield ValidationError(path, f"{value} is less than the minimum {schema['minimum']}")
        if "maximum" in schema and value > schema["maximum"]:
            yield ValidationError(path, f"{value} is greater than the maximum {schema['maximum']}")
        if "exclusiveMinimum" in schema and value <= schema["exclusiveMinimum"]:
            yield ValidationError(path, f"{value} is not greater than {schema['exclusiveMinimum']}")
        if "exclusiveMaximum" in schema and value >= schema["exclusiveMaximum"]:
            yield ValidationError(path, f"{value} is not less than {schema['exclusiveMaximum']}")

    if isinstance(value, str):
        if "minLength" in schema and len(value) < schema["minLength"]:
            yield ValidationError(path, f"string shorter than {schema['minLength']}")
        if "maxLength" in schema and len(value) > schema["maxLength"]:
            yield ValidationError(path, f"string longer than {schema['maxLength']}")
        if "pattern" in schema and re.search(schema["pattern"], value) is None:
            yield ValidationError(path, f"string {value!r} does not match pattern {schema['pattern']!r}")

    if "allOf" in schema:
        for index, sub in enumerate(schema["allOf"]):
            yield from _iter(root, sub, value, path)
    if "anyOf" in schema:
        if not any(_valid(root, sub, value, path) for sub in schema["anyOf"]):
            yield ValidationError(path, "value does not match any of the anyOf alternatives")
    if "oneOf" in schema:
        matches = sum(1 for sub in schema["oneOf"] if _valid(root, sub, value, path))
        if matches != 1:
            yield ValidationError(path, f"value matches {matches} oneOf alternatives, expected exactly 1")
    if "not" in schema and _valid(root, schema["not"], value, path):
        yield ValidationError(path, "value matches the forbidden 'not' schema")


def _valid(root: dict, schema: Any, value: Any, path: str) -> bool:
    for _ in _iter(root, schema, value, path):
        return False
    return True


def iter_errors(schema: dict, document: Any) -> Iterator[ValidationError]:
    """Yield every validation error of ``document`` against ``schema``."""
    yield from _iter(schema, schema, document, "")


def validate(schema: dict, document: Any) -> None:
    """Raise :class:`ValidationError` with the first problem found, if any."""
    for error in iter_errors(schema, document):
        raise error


def is_valid(schema: dict, document: Any) -> bool:
    return not any(True for _ in iter_errors(schema, document))


def collect_errors(schema: dict, document: Any) -> List[str]:
    return [str(error) for error in iter_errors(schema, document)]


def _main(argv: Sequence[str]) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Validate JSON documents against a Vela contract schema.")
    parser.add_argument("schema", type=Path)
    parser.add_argument("documents", type=Path, nargs="+")
    args = parser.parse_args(argv)
    schema = load_schema(args.schema)
    status = 0
    for document_path in args.documents:
        document = json.loads(document_path.read_text(encoding="utf-8"))
        errors = collect_errors(schema, document)
        if errors:
            status = 1
            print(f"{document_path}: INVALID")
            for line in errors:
                print(f"  {line}")
        else:
            print(f"{document_path}: ok")
    return status


if __name__ == "__main__":  # pragma: no cover
    import sys

    raise SystemExit(_main(sys.argv[1:]))
