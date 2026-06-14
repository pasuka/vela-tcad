#!/usr/bin/env python3
"""Probe candidate drivers for the pn2d 0 V quasi-Fermi split."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any


QF_EQUILIBRIUM_TOL_V = 1.0e-6


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def resolve_deck_path(deck: dict[str, Any], reference_root: Path, key: str) -> None:
    if key not in deck:
        return
    candidate = Path(str(deck[key]))
    if not candidate.is_absolute():
        deck[key] = str((reference_root / "vela" / candidate).resolve())


def scalar_stats(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"points": 0, "min": None, "max": None, "span": None, "mean": None}
    return {
        "points": len(values),
        "min": min(values),
        "max": max(values),
        "span": max(values) - min(values),
        "mean": sum(values) / len(values),
    }


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    lines = path.read_text().splitlines()
    npoints = 0
    for line in lines:
        parts = line.split()
        if len(parts) >= 2 and parts[0] == "POINTS":
            npoints = int(parts[1])
            break

    fields: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < npoints:
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                head = line.split()[0]
                if head in {"SCALARS", "CELL_DATA", "POINT_DATA", "VECTORS"}:
                    break
                values.extend(float(item) for item in line.split())
                i += 1
            fields[name] = values[:npoints]
            continue
        i += 1
    return fields


def vtk_report(path: Path | None) -> dict[str, Any]:
    required = ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi"]
    if path is None or not path.exists():
        return {
            "path": None if path is None else str(path),
            "missing_fields": required,
            "fields": {name: scalar_stats([]) for name in required},
            "qf_max_span_V": None,
            "near_equilibrium_qf": False,
        }

    fields = parse_vtk_scalars(path)
    missing = [name for name in required if name not in fields]
    stats = {name: scalar_stats(fields.get(name, [])) for name in required}
    qf_max_span = max(
        (stats[name]["span"] or 0.0 for name in ["ElectronQuasiFermi", "HoleQuasiFermi"]),
        default=0.0,
    )
    return {
        "path": str(path),
        "missing_fields": missing,
        "fields": stats,
        "qf_max_span_V": qf_max_span,
        "near_equilibrium_qf": qf_max_span <= QF_EQUILIBRIUM_TOL_V,
        "qf_equilibrium_tolerance_V": QF_EQUILIBRIUM_TOL_V,
    }


def first_csv_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows[0] if rows else {}


def terminal_balance_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def runner_command(raw: str) -> list[str]:
    command = shlex.split(raw, posix=not sys.platform.startswith("win"))
    if not command:
        raise ValueError("--runner must not be empty")
    runner = Path(command[0])
    if not runner.is_absolute():
        command[0] = str((Path.cwd() / runner).resolve())
    return command


def variant_deck(
    base: dict[str, Any],
    reference_root: Path,
    work_dir: Path,
    name: str,
    solver_patch: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Path]]:
    deck = json.loads(json.dumps(base))
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        resolve_deck_path(deck, reference_root, key)

    sweep_csv = work_dir / f"{name}.csv"
    terminal_csv = work_dir / f"{name}_terminal_balance.csv"
    edge_csv = work_dir / f"{name}_contact_edges.csv"
    vtk_prefix = work_dir / name

    deck["simulation_type"] = "dc_sweep"
    deck["output_csv"] = str(sweep_csv.resolve())
    solver = deck.setdefault("solver", {})
    solver.update(solver_patch)
    solver["diagnostics"] = True
    handoff = solver.setdefault("handoff", {})
    handoff["fallback"] = "none"
    handoff["require_gummel_convergence"] = False
    handoff["gummel_max_iter"] = 0
    handoff["newton_max_iter"] = solver.get("max_iter", 80)

    sweep = deck.setdefault("sweep", {})
    sweep.update({
        "mode": "iv",
        "contact": "Anode",
        "current_contact": "Anode",
        "start": 0.0,
        "stop": 0.0,
        "step": 1.0,
        "write_vtk": True,
        "vtk_prefix": str(vtk_prefix.resolve()),
        "diagnostics": {
            "terminal_balance": {
                "enabled": True,
                "contacts": ["Anode", "Cathode"],
                "csv_file": str(terminal_csv.resolve()),
            },
            "contact_edge": {
                "enabled": True,
                "contacts": ["Anode", "Cathode"],
                "csv_file": str(edge_csv.resolve()),
            },
        },
    })

    return deck, {
        "sweep_csv": sweep_csv,
        "terminal_csv": terminal_csv,
        "edge_csv": edge_csv,
        "vtk": work_dir / f"{name}_0000_0V.vtk",
    }


def run_variant(command: list[str], deck_path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [*command, "--config", str(deck_path.resolve())],
        cwd=deck_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def to_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def summarize_variant(
    name: str,
    patch: dict[str, Any],
    paths: dict[str, Path],
    result: subprocess.CompletedProcess[str],
) -> dict[str, Any]:
    sweep = first_csv_row(paths["sweep_csv"])
    terminal = terminal_balance_rows(paths["terminal_csv"])
    vtk = vtk_report(paths["vtk"])
    currents = {
        row.get("contact", ""): to_float(row.get("current_total_A_per_um"))
        for row in terminal
        if row.get("contact")
    }
    current_sum = sum(value for value in currents.values() if value is not None)
    max_abs_current = max((abs(value) for value in currents.values() if value is not None), default=0.0)
    pair_balance_relative = (
        abs(current_sum) / max_abs_current if max_abs_current > 0.0 else None
    )
    return {
        "variant": name,
        "solver_patch": patch,
        "returncode": result.returncode,
        "converged": sweep.get("converged"),
        "newton_iterations": sweep.get("newton_iterations"),
        "handoff_stage": sweep.get("handoff_stage"),
        "failure_reason": sweep.get("failure_reason"),
        "newton_failure_class": sweep.get("newton_failure_class"),
        "current_total_A_per_um": to_float(sweep.get("current_total_A_per_um")),
        "terminal_current_sum_A_per_um": current_sum,
        "terminal_pair_balance_relative": pair_balance_relative,
        "qf_max_span_V": vtk["qf_max_span_V"],
        "near_equilibrium_qf": vtk["near_equilibrium_qf"],
        "vtk": vtk,
        "paths": {key: str(path) for key, path in paths.items()},
        "stdout_tail": result.stdout[-2000:],
        "stderr_tail": result.stderr[-2000:],
    }


def write_csv_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = [
        "variant",
        "returncode",
        "converged",
        "newton_iterations",
        "handoff_stage",
        "failure_reason",
        "newton_failure_class",
        "current_total_A_per_um",
        "terminal_current_sum_A_per_um",
        "terminal_pair_balance_relative",
        "qf_max_span_V",
        "near_equilibrium_qf",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fields})


def write_markdown_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    lines = [
        "# PN2D 0V QF Driver Probe",
        "",
        "| Variant | Converged | Newton iters | Handoff | QF max span (V) | I total (A/um) | Pair balance |",
        "| --- | ---: | ---: | --- | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            "| {variant} | {converged} | {newton_iterations} | {handoff_stage} | "
            "{qf_max_span_V} | {current_total_A_per_um} | {terminal_pair_balance_relative} |".format(**row)
        )
    lines.append("")
    lines.append(
        f"A variant with `qf_max_span_V <= {QF_EQUILIBRIUM_TOL_V:g}` is treated "
        "as near-equilibrium for this probe."
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    args.reference_root = args.reference_root.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    base_path = args.reference_root / "vela" / "simulation_0v.json"
    base = read_json(base_path)
    command = runner_command(args.runner)

    variants: dict[str, dict[str, Any]] = {
        "baseline": {},
        "no_recombination": {"recombination": []},
        "no_bgn": {"bandgap_narrowing": "none"},
        "l2_residual": {"residual_norm": "l2"},
        "tight_block_scales": {
            "residual_scales": {"psi": 1.0, "phin": 1.0e-24, "phip": 1.0e-24},
        },
    }

    summary: list[dict[str, Any]] = []
    for name, patch in variants.items():
        work_dir = args.output_dir / name
        work_dir.mkdir(parents=True, exist_ok=True)
        deck, paths = variant_deck(base, args.reference_root, work_dir, name, patch)
        deck_path = work_dir / f"{name}.json"
        write_json(deck_path, deck)
        result = run_variant(command, deck_path)
        summary.append(summarize_variant(name, patch, paths, result))

    report = {
        "reference_root": str(args.reference_root),
        "runner": command,
        "status": "complete" if all(row["returncode"] == 0 for row in summary)
        else "complete_with_variant_failures",
        "all_variants_returned_zero": all(row["returncode"] == 0 for row in summary),
        "variants": summary,
    }
    write_json(args.output_dir / "pn2d_0v_qf_driver_summary.json", report)
    write_csv_summary(args.output_dir / "pn2d_0v_qf_driver_summary.csv", summary)
    write_markdown_summary(args.output_dir / "pn2d_0v_qf_driver_summary.md", summary)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
