#!/usr/bin/env python3
"""Run a physics-disabled exact-mesh Poisson cost probe and draft budgets.

The probe is a structural lower bound.  Its placeholder PolySilicon and
Nitride permittivities are not a materials contract and its results are never
used as electrical validation data.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
import hashlib
import json
import os
import shutil
import subprocess
import time
from ctypes import wintypes
from pathlib import Path
from typing import Any

from templates_ldmos_contracts import draft_governance_contracts, validate_document


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PROCESS_MEMORY_COUNTERS_EX(ctypes.Structure):
    _fields_ = [
        ("cb", wintypes.DWORD),
        ("PageFaultCount", wintypes.DWORD),
        ("PeakWorkingSetSize", ctypes.c_size_t),
        ("WorkingSetSize", ctypes.c_size_t),
        ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPagedPoolUsage", ctypes.c_size_t),
        ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
        ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
        ("PagefileUsage", ctypes.c_size_t),
        ("PeakPagefileUsage", ctypes.c_size_t),
        ("PrivateUsage", ctypes.c_size_t),
    ]


def windows_peak_rss(process: subprocess.Popen[str]) -> int | None:
    if os.name != "nt":
        return None
    counters = PROCESS_MEMORY_COUNTERS_EX()
    counters.cb = ctypes.sizeof(counters)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    if not psapi.GetProcessMemoryInfo(
        wintypes.HANDLE(int(process._handle)),  # type: ignore[attr-defined]
        ctypes.byref(counters),
        counters.cb,
    ):
        return None
    return int(counters.PeakWorkingSetSize)


def measured_run(command: list[str], cwd: Path) -> dict[str, Any]:
    started = time.perf_counter()
    process = subprocess.Popen(
        command,
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    peak = 0
    while process.poll() is None:
        sample = windows_peak_rss(process)
        if sample is not None:
            peak = max(peak, sample)
        time.sleep(0.02)
    output, _ = process.communicate()
    elapsed = time.perf_counter() - started
    sample = windows_peak_rss(process)
    if sample is not None:
        peak = max(peak, sample)
    return {
        "command": command,
        "return_code": process.returncode,
        "wall_clock_seconds": elapsed,
        "peak_working_set_bytes": peak or None,
        "stdout": output,
    }


def mesh_pattern(mesh: dict[str, Any]) -> dict[str, int]:
    nodes = len(mesh["nodes"])
    edges: set[tuple[int, int]] = set()
    for triangle in mesh["triangles"]:
        a, b, c = map(int, triangle["node_ids"])
        for left, right in ((a, b), (b, c), (c, a)):
            edges.add((left, right) if left < right else (right, left))
    contact_nodes = {int(node) for contact in mesh["contacts"] for node in contact["node_ids"]}
    poisson_upper = nodes + 2 * len(edges)
    poisson_dirichlet_row = poisson_upper - sum(
        sum(1 for edge in edges if node in edge) for node in contact_nodes
    )
    return {
        "nodes": nodes,
        "unique_edges": len(edges),
        "poisson_structural_nnz_upper_bound": poisson_upper,
        "poisson_dirichlet_row_nnz_estimate": poisson_dirichlet_row,
        "coupled_dd_unknowns_estimate": 3 * nodes,
    }


def count_iterations(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"rows": 0, "accepted_rows": 0, "note": "diagnostic file not emitted"}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    accepted = sum(
        str(row.get("accepted", "")).lower() in {"1", "true", "yes"}
        or "accepted" in str(row.get("row_type", row.get("event", ""))).lower()
        for row in rows
    )
    return {"rows": len(rows), "accepted_rows": accepted, "columns": list(rows[0]) if rows else []}


def prepare_deck(source: dict[str, Any], output_dir: Path, name: str, bias: float) -> Path:
    deck = json.loads(json.dumps(source))
    deck["_comment"] = (
        "Templates/LDMOS stage-1 structural cost probe only. All new physics is disabled; "
        "placeholder dielectric/poly properties are not a validation physics contract."
    )
    deck["mesh_file"] = "mesh.json"
    deck["node_doping_file"] = "doping.csv"
    deck["materials_file"] = "materials_probe_only.json"
    deck["simulation_type"] = "poisson"
    deck["output_vtk"] = f"{name}.vtk"
    deck.pop("output_csv", None)
    deck.pop("solver", None)
    deck.pop("sweep", None)
    deck.pop("node_doping_file", None)
    # This is deliberately a matrix/mesh lower bound, not a device solution.
    # Region-average values inherited from the converter would make the timing
    # depend on an unqualified materials/doping contract.
    deck["doping"] = [
        {"region": item["region"], "donors": 0.0, "acceptors": 0.0}
        for item in deck.get("doping", [])
    ]
    for contact in deck.get("contacts", []):
        contact["bias"] = bias if contact.get("name") == "drain" else 0.0
    path = output_dir / f"{name}.json"
    write_json(path, deck)
    return path


def measured_poisson_sequence(runner: Path,
                              decks: list[Path],
                              cwd: Path) -> dict[str, Any]:
    point_runs = [
        measured_run([str(runner), "--config", deck.name], cwd)
        for deck in decks
    ]
    return {
        "commands": [item["command"] for item in point_runs],
        "return_code": 0 if all(item["return_code"] == 0 for item in point_runs) else 1,
        "wall_clock_seconds": sum(item["wall_clock_seconds"] for item in point_runs),
        "peak_working_set_bytes": max(
            (item.get("peak_working_set_bytes") or 0 for item in point_runs),
            default=0,
        ) or None,
        "stdout": "\n".join(item["stdout"] for item in point_runs),
        "points": point_runs,
        "newton_iterations": {
            "rows": 0,
            "accepted_rows": 0,
            "note": "linear Poisson structural probe; Newton is not invoked",
        },
    }


def budget_from_probe(vela_commit: str,
                      oracle_manifest_sha256: str | None,
                      probe: dict[str, Any],
                      staging_bytes: int) -> dict[str, Any]:
    result = draft_governance_contracts(vela_commit, oracle_manifest_sha256)["budget_freeze.json"]
    single = probe["runs"]["single_bias"]
    short = probe["runs"]["five_point"]
    per_point = max(short["wall_clock_seconds"] / 5.0, single["wall_clock_seconds"] * 0.2)
    peak_gib = max(
        (item.get("peak_working_set_bytes") or 0) / (1024 ** 3)
        for item in probe["runs"].values()
    )
    point_estimates = {"phase_a": 360, "phase_a_plus": 180, "phase_b": 600}
    # Poisson-only is a lower bound.  Multipliers intentionally cover coupled
    # carrier/DG/thermal unknowns, continuation retries, and diagnostics.
    scope_multipliers = {"phase_a": 6.0, "phase_a_plus": 12.0, "phase_b": 30.0}
    scopes = []
    for scope in ("phase_a", "phase_a_plus", "phase_b"):
        base_hours = per_point * point_estimates[scope] * scope_multipliers[scope] / 3600.0
        scopes.append({
            "scope": scope,
            "status": "estimated",
            "best_wall_clock_hours": base_hours * 0.5,
            "base_wall_clock_hours": base_hours,
            "worst_wall_clock_hours": base_hours * 2.5,
            "estimated_bias_points": point_estimates[scope],
        })
    total_base = sum(item["base_wall_clock_hours"] for item in scopes)
    result["basis"]["probe"] = {
        "status": "structural_lower_bound_complete",
        "single_bias_wall_seconds": single["wall_clock_seconds"],
        "five_point_wall_seconds": short["wall_clock_seconds"],
        "derived_per_point_seconds": per_point,
        "peak_memory_gib": peak_gib,
        "matrix_nnz_upper_bound": probe["matrix_pattern"]["poisson_structural_nnz_upper_bound"],
        "physics_contract": "placeholder_poisson_only_not_for_acceptance",
    }
    result["scope_budgets"] = scopes
    storage_gib = staging_bytes / (1024 ** 3)
    result["scenarios"] = [
        {"name": "best", "wall_clock_hours": total_base * 0.5, "peak_memory_gib": peak_gib * 3.0,
         "storage_gib": max(storage_gib * 3.0, 5.0), "max_concurrency": 1},
        {"name": "base", "wall_clock_hours": total_base, "peak_memory_gib": peak_gib * 6.0,
         "storage_gib": max(storage_gib * 8.0, 15.0), "max_concurrency": 1},
        {"name": "worst", "wall_clock_hours": total_base * 2.5, "peak_memory_gib": peak_gib * 12.0,
         "storage_gib": max(storage_gib * 20.0, 40.0), "max_concurrency": 1},
    ]
    result["matrix_policy"] = (
        "Draft only: use exact common bias points; keep equilibrium, threshold/knee, and endpoint points; "
        "trim redundant smooth-curve interior points by information gain before execution, never afterward. "
        "Poisson timings are lower bounds; re-freeze after stage 1.5 and first classic-DD exact-mesh run."
    )
    validate_document(result)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage1-dir", type=Path, required=True)
    parser.add_argument("--vela-runner", type=Path, required=True)
    parser.add_argument("--vela-commit", required=True)
    parser.add_argument("--oracle-manifest", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    stage1 = args.stage1_dir.resolve()
    converted = stage1 / "vela_exact_topology"
    probe_dir = stage1 / "cost_probe"
    if probe_dir.exists():
        raise FileExistsError(f"refusing to overwrite cost probe: {probe_dir}")
    probe_dir.mkdir(parents=True)
    for name in ("mesh.json", "doping.csv"):
        shutil.copyfile(converted / name, probe_dir / name)
    write_json(probe_dir / "materials_probe_only.json", {
        "materials": [
            {"name": "PolySilicon", "eps_r": 11.7, "ni": 0.0, "mun": 0.0, "mup": 0.0},
            {"name": "Nitride", "eps_r": 7.5, "ni": 0.0, "mun": 0.0, "mup": 0.0},
        ],
        "warning": "Probe-only placeholders; not a Templates/LDMOS materials contract.",
    })
    source = json.loads((converted / "simulation_iv.json").read_text(encoding="utf-8"))
    single_decks = [prepare_deck(source, probe_dir, "single_bias", 0.0)]
    five_decks = [
        prepare_deck(source, probe_dir, f"five_point_{index}", bias)
        for index, bias in enumerate((0.0, 0.1, 0.2, 0.3, 0.4))
    ]
    runs = {
        "single_bias": measured_poisson_sequence(args.vela_runner, single_decks, probe_dir),
        "five_point": measured_poisson_sequence(args.vela_runner, five_decks, probe_dir),
    }
    mesh = json.loads((probe_dir / "mesh.json").read_text(encoding="utf-8"))
    probe = {
        "schema": "vela.templates_ldmos.cost_probe.v1",
        "classification": "structural_poisson_lower_bound_not_physics_acceptance",
        "runner": {"path": str(args.vela_runner), "sha256": sha256(args.vela_runner)},
        "matrix_pattern": mesh_pattern(mesh),
        "runs": runs,
    }
    write_json(probe_dir / "cost_probe.json", probe)
    oracle_hash = sha256(args.oracle_manifest) if args.oracle_manifest and args.oracle_manifest.is_file() else None
    staging_bytes = sum(path.stat().st_size for path in stage1.rglob("*") if path.is_file())
    budget = budget_from_probe(args.vela_commit, oracle_hash, probe, staging_bytes)
    write_json(stage1 / "contracts" / "budget_freeze.json", budget)
    print(json.dumps({
        "return_codes": {name: item["return_code"] for name, item in runs.items()},
        "matrix_pattern": probe["matrix_pattern"],
        "budget_approval": budget["approval"]["status"],
    }, indent=2))
    return 0 if all(item["return_code"] == 0 for item in runs.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
