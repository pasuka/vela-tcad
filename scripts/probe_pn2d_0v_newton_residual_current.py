#!/usr/bin/env python3
"""Probe whether relaxed Newton stopping leaves Sentaurus-like 0V residual current."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any

from diagnose_pn2d_0v_current_balance import (
    DEFAULT_CONTACTS,
    edge_report,
    resolve_deck_path,
    sentaurus_current_parity,
    sentaurus_current_reference,
    terminal_report,
    vtk_report,
    write_json,
)


REPORT_NAME = "pn2d_0v_newton_residual_probe.json"
MARKDOWN_NAME = "pn2d_0v_newton_residual_probe.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def run_command(command: list[str], cwd: Path) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{result.returncode}: {' '.join(command)}\nstdout={result.stdout}\nstderr={result.stderr}"
        )


def derive_candidate_deck(
    reference_root: Path,
    output_dir: Path,
    name: str,
    reltol: float,
    abstol: float,
    max_iter: int,
) -> tuple[Path, Path, Path]:
    base_path = reference_root / "vela" / "simulation_0v.json"
    deck = read_json(base_path)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        resolve_deck_path(deck, reference_root, key)

    candidate_dir = output_dir / name
    candidate_dir.mkdir(parents=True, exist_ok=True)
    terminal_csv = candidate_dir / "terminal_balance.csv"
    edge_csv = candidate_dir / "contact_edges.csv"

    solver = deck.setdefault("solver", {})
    solver["method"] = "gummel_newton"
    solver["reltol"] = reltol
    solver["abstol"] = abstol
    solver["max_iter"] = max_iter
    handoff = solver.setdefault("handoff", {})
    handoff["fallback"] = "none"
    handoff["require_gummel_convergence"] = False
    handoff["gummel_max_iter"] = 0
    handoff["newton_max_iter"] = max_iter

    deck["simulation_type"] = "dc_sweep"
    deck["output_csv"] = str((candidate_dir / "pn2d_0v_probe.csv").resolve())
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = sweep.get("current_contact", "Anode")
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((candidate_dir / "pn2d_0v_probe").resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": DEFAULT_CONTACTS,
        "csv_file": str(terminal_csv.resolve()),
    }
    diagnostics["contact_edge"] = {
        "enabled": True,
        "contacts": DEFAULT_CONTACTS,
        "csv_file": str(edge_csv.resolve()),
    }

    deck_path = candidate_dir / "simulation_0v_probe.json"
    write_json(deck_path, deck)
    return deck_path, terminal_csv, edge_csv


def run_candidate(
    reference_root: Path,
    output_dir: Path,
    runner: str,
    name: str,
    reltol: float,
    abstol: float,
    max_iter: int,
) -> dict[str, Any]:
    deck_path, terminal_csv, edge_csv = derive_candidate_deck(
        reference_root, output_dir, name, reltol, abstol, max_iter
    )
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    try:
        run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    except RuntimeError as exc:
        return {
            "name": name,
            "status": "runner_failed",
            "deck": str(deck_path),
            "solver": {
                "reltol": reltol,
                "abstol": abstol,
                "max_iter": max_iter,
            },
            "error": str(exc),
            "sentaurus_total_abs_ratio_median_log10_error": 300.0,
        }
    vtk_matches = sorted(deck_path.parent.glob("pn2d_0v_probe*.vtk"))
    if not vtk_matches:
        return {
            "name": name,
            "status": "runner_failed",
            "deck": str(deck_path),
            "solver": {
                "reltol": reltol,
                "abstol": abstol,
                "max_iter": max_iter,
            },
            "error": f"candidate {name} did not produce a VTK probe",
            "sentaurus_total_abs_ratio_median_log10_error": 300.0,
        }

    terminal = terminal_report(terminal_csv, relative_gate=5.0e-2, abs_gate=1.0e-24)
    sentaurus_reference = sentaurus_current_reference(reference_root)
    parity = sentaurus_current_parity(terminal, sentaurus_reference)
    edges = edge_report(edge_csv, terminal, top_n=5)
    vtk = vtk_report(vtk_matches[0])
    score = sentaurus_total_abs_ratio_median_log10_error(parity)
    return {
        "name": name,
        "status": "pass",
        "deck": str(deck_path),
        "solver": {
            "reltol": reltol,
            "abstol": abstol,
            "max_iter": max_iter,
        },
        "terminal_balance": terminal,
        "sentaurus_current_parity": parity,
        "contact_edges": edges,
        "vtk_fields": vtk,
        "sentaurus_total_abs_ratio_median_log10_error": score,
    }


def sentaurus_total_abs_ratio_median_log10_error(parity: dict[str, Any]) -> float:
    values: list[float] = []
    for row in parity.get("by_contact", {}).values():
        ratio = row.get("abs_ratio_sentaurus_to_vela")
        if ratio is None:
            continue
        if math.isinf(float(ratio)):
            values.append(300.0)
        elif float(ratio) > 0.0:
            values.append(abs(math.log10(float(ratio))))
    if not values:
        return 300.0
    values.sort()
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return 0.5 * (values[mid - 1] + values[mid])


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Newton Residual Current Probe",
        "",
        f"Status: {report.get('status')}",
        "",
        "| Candidate | reltol | abstol | max_iter | log10 ratio error | parity |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    for item in report.get("candidates", []):
        solver = item.get("solver", {})
        parity = item.get("sentaurus_current_parity", {})
        lines.append(
            f"| {item.get('name')} | {solver.get('reltol')} | {solver.get('abstol')} | "
            f"{solver.get('max_iter')} | {item.get('sentaurus_total_abs_ratio_median_log10_error')} | "
            f"{parity.get('status')} |"
        )
    best = report.get("best_candidate")
    if best:
        lines.extend([
            "",
            "## Best Candidate",
            "",
            f"- name: {best.get('name')}",
            f"- solver: {best.get('solver')}",
            f"- log10 ratio error: {best.get('sentaurus_total_abs_ratio_median_log10_error')}",
        ])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--reltol-candidate", action="append", type=float, default=[])
    parser.add_argument("--abstol-candidate", action="append", type=float, default=[])
    parser.add_argument("--max-iter-candidate", action="append", type=int, default=[])
    args = parser.parse_args()

    reltols = args.reltol_candidate or [1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2]
    abstols = args.abstol_candidate or [1.0e-24]
    max_iters = args.max_iter_candidate or [80, 10, 3, 1]
    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidates: list[dict[str, Any]] = []
    for index, (reltol, abstol, max_iter) in enumerate(itertools.product(reltols, abstols, max_iters)):
        name = f"candidate_{index:02d}_rtol_{reltol:.0e}_atol_{abstol:.0e}_iter_{max_iter}"
        candidates.append(
            run_candidate(args.reference_root, args.output_dir, args.runner, name, reltol, abstol, max_iter)
        )

    successful = [item for item in candidates if item.get("status") == "pass"]
    balanced = [
        item for item in successful
        if item.get("terminal_balance", {}).get("status") == "pass"
    ]
    best_pool = balanced or successful or candidates
    best = min(best_pool, key=lambda item: item["sentaurus_total_abs_ratio_median_log10_error"])
    report = {
        "status": "pass",
        "reference_root": str(args.reference_root),
        "sentaurus_current_reference": sentaurus_current_reference(args.reference_root),
        "best_selection_policy": "prefer_successful_terminal_balanced_then_lowest_log10_ratio_error",
        "best_candidate": best,
        "candidates": candidates,
    }
    write_json(args.output_dir / REPORT_NAME, report)
    (args.output_dir / MARKDOWN_NAME).write_text(markdown_report(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
