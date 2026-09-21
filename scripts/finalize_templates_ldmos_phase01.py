#!/usr/bin/env python3
"""Assemble the Templates/LDMOS phase-A prerequisite validation summary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from run_templates_ldmos_sentaurus_vm import sha256_file, utc_now, write_json
from templates_ldmos_contracts import BENCHMARK, render_summary, validate_document


def read(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def evidence(path: Path, root: Path, description: str) -> dict[str, str]:
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": sha256_file(path),
        "description": description,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--stage1-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    stage1 = args.stage1_dir.resolve()
    oracle_path = run_dir / "reports" / "oracle_audit.json"
    structure_path = stage1 / "reports" / "structure_audit.json"
    repeat_path = stage1 / "reports" / "repeatability.json"
    budget_path = stage1 / "contracts" / "budget_freeze.json"
    cost_path = stage1 / "cost_probe" / "cost_probe.json"
    required = (oracle_path, structure_path, repeat_path, budget_path, cost_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing phase-0/1 evidence: {missing}")
    oracle = read(oracle_path)
    structure = read(structure_path)
    repeat = read(repeat_path)
    budget = read(budget_path)
    cost = read(cost_path)
    restart_path = stage1 / "qualification" / "qualification_summary.json"
    wp175_path = stage1 / "contracts" / "wp175" / "qualification_report.json"
    restart = read(restart_path) if restart_path.is_file() else None
    wp175 = read(wp175_path) if wp175_path.is_file() else None
    # Keep the run-level governance location canonical even when the long VM
    # runner completed with an older in-memory WP0 schema implementation.
    # The stage-1 copy remains evidence and is hashed below.
    write_json(run_dir / "manifest" / "budget_freeze.json", budget, validate=True)
    oracle_pass = all(oracle["gates"].values())
    structure_pass = all(structure["gates"].values()) and repeat["identical"]
    probe_pass = all(item["return_code"] == 0 for item in cost["runs"].values())
    budget_approved = budget["approval"]["status"] == "approved"
    gates = [
        {
            "id": "official_oracle_and_representative_states",
            "status": "pass" if oracle_pass else "fail",
            "summary": "All stage-0 exit, curve, log, final-structure, and representative-state gates.",
        },
        {
            "id": "exact_topology_structure",
            "status": "pass" if structure_pass else "fail",
            "summary": "Coordinates, triangles, contacts, doping, mesh anomalies, and deterministic repeat.",
        },
        {
            "id": "exact_mesh_cost_probe",
            "status": "pass" if probe_pass else "fail",
            "summary": "Poisson-only structural lower bound; not physics acceptance.",
        },
        {
            "id": "budget_double_approval",
            "status": "pass" if budget_approved else "unresolved",
            "summary": f"budget_freeze approval.status={budget['approval']['status']}",
        },
    ]
    if restart is not None:
        gates.append({
            "id": "exact_mesh_restart_qualification",
            "status": restart["status"],
            "summary": (
                "17-digit round-trip, frozen replay, equilibrium/IdVg/IdVd "
                "same-bias reclose, resolved current and KCL gates."
            ),
        })
    if wp175 is not None:
        gates.append({
            "id": "wp175_versioned_contracts",
            "status": wp175["status"],
            "summary": (
                "Strict material, solver-physics and phase-A discretization "
                "schemas, migration, unit and round-trip checks."
            ),
        })
    restart_pass = restart is None or restart["status"] == "pass"
    wp175_pass = wp175 is None or wp175["status"] == "pass"
    technical_pass = (
        oracle_pass and structure_pass and probe_pass and restart_pass and wp175_pass)
    if not technical_pass:
        status, highest = "fail", ("L0" if oracle_pass else "none")
    elif not budget_approved:
        status, highest = "unresolved", "L0"
    else:
        status, highest = "pass", "L1"
    limitations = [
        "No classic DD curve acceptance or new phase-A/phase-B physics acceptance was entered.",
        "The cost probe is a Poisson-only lower bound and must be re-frozen after the first qualified classic-DD run.",
    ]
    if restart is None:
        limitations.append("Stage 1.5 restart qualification has not been attached.")
    if wp175 is None:
        limitations.append("WP1.75 versioned contracts have not been attached.")
    if not budget_approved:
        limitations.append(
            "L1 remains governance-unresolved until budget_freeze has both required approvals."
        )
    summary = {
        "schema": "vela.templates_ldmos.validation_summary.v1",
        "benchmark": BENCHMARK,
        "generated_at": utc_now(),
        "status": status,
        "highest_level": highest,
        "oracle": gates[0],
        "structure": gates[1],
        "budget": gates[3],
        "gates": gates,
        "evidence": [
            evidence(oracle_path, run_dir, "stage-0 oracle audit"),
            evidence(structure_path, run_dir, "stage-1 structural audit"),
            evidence(repeat_path, run_dir, "deterministic repeat audit"),
            evidence(cost_path, run_dir, "exact-mesh structural cost probe"),
            evidence(budget_path, run_dir, "draft or approved execution budget"),
        ] + ([evidence(
            restart_path, run_dir, "exact-mesh restart and reclose qualification")]
             if restart is not None else []) + ([evidence(
            wp175_path, run_dir, "WP1.75 contract qualification")]
             if wp175 is not None else []),
        "limitations": limitations,
    }
    validate_document(summary)
    output = run_dir / "reports" / "validation_summary.json"
    write_json(output, summary, validate=True)
    (run_dir / "reports" / "validation_summary.md").write_text(
        render_summary(summary), encoding="utf-8"
    )
    print(json.dumps({"status": status, "highest_level": highest, "gates": gates}, indent=2))
    return 0 if technical_pass else 2


if __name__ == "__main__":
    raise SystemExit(main())
