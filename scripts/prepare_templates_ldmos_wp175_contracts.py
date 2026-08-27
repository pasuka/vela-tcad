#!/usr/bin/env python3
"""Prepare validated Templates/LDMOS WP1.75 contracts for one stage-1 run."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

from templates_ldmos_contracts import (
    canonical_round_trip,
    migrate_discretization_draft,
    read_json,
)


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = REPO / "reference_tcad" / "templates_ldmos_sentaurus2022" / "contracts"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(stage1_dir: Path,
            source_dir: Path = DEFAULT_SOURCE) -> dict[str, Any]:
    source_materials = source_dir / "materials.json"
    source_physics = source_dir / "physics_contract.json"
    draft_discretization = stage1_dir / "contracts" / "discretization_contract.json"
    for path in (source_materials, source_physics, draft_discretization):
        if not path.is_file():
            raise ValueError(f"WP1.75 input is missing: {path}")

    materials = canonical_round_trip(read_json(source_materials))
    physics = canonical_round_trip(read_json(source_physics))
    discretization = migrate_discretization_draft(read_json(draft_discretization))
    if physics["materials_file"] != "materials.json":
        raise ValueError("physics_contract.materials_file must be materials.json")
    if physics["discretization_profile"] != discretization["profile_name"]:
        raise ValueError(
            "physics and discretization contracts name different profiles")

    output = stage1_dir / "contracts" / "wp175"
    output.mkdir(parents=True, exist_ok=True)
    documents = {
        "materials.json": materials,
        "physics_contract.json": physics,
        "discretization_contract.json": discretization,
    }
    for name, document in documents.items():
        (output / name).write_text(
            json.dumps(document, indent=2) + "\n", encoding="utf-8")

    report = {
        "status": "pass",
        "output_directory": str(output),
        "contracts": [
            {
                "path": name,
                "schema": document["schema"],
                "sha256": sha256_file(output / name),
            }
            for name, document in sorted(documents.items())
        ],
        "checks": [
            "all contracts pass their versioned schema",
            "canonical JSON round-trip is lossless",
            "physics materials_file resolves to materials.json",
            "physics and discretization profile names match",
            "phase-A discretization explicitly forbids avalanche use",
        ],
    }
    (output / "qualification_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Templates/LDMOS WP1.75 contract qualification",
        "",
        "- Status: `pass`",
        "- Scope: versioned material, solver-physics and phase-A discretization contracts",
        "",
        "| Contract | Schema | SHA-256 |",
        "| --- | --- | --- |",
    ]
    for item in report["contracts"]:
        lines.append(
            f"| `{item['path']}` | `{item['schema']}` | `{item['sha256']}` |")
    lines.extend(["", "## Checks", ""])
    lines.extend(f"- {item}" for item in report["checks"])
    (output / "qualification_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1-dir", type=Path, required=True)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    args = parser.parse_args(argv)
    report = prepare(args.stage1_dir, args.source_dir)
    print(json.dumps({"status": report["status"], "output": report["output_directory"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
