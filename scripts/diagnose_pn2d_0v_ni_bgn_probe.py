#!/usr/bin/env python3
"""Run PN2D 0V ni/BGN probe candidates and rank state parity against Sentaurus."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shlex
import statistics
import subprocess
import sys
from pathlib import Path
from typing import Any


REPORT_NAME = "pn2d_0v_ni_bgn_probe.json"
MARKDOWN_NAME = "pn2d_0v_ni_bgn_probe.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def field_files(root: Path, field: str) -> list[Path]:
    exact = root / f"{field}_region0.csv"
    if exact.is_file():
        return [exact]
    return sorted(root.glob(f"{field}_region*.csv"))


def read_sentaurus_scalar(fields_root: Path, field: str) -> dict[int, float]:
    values: dict[int, float] = {}
    for path in field_files(fields_root, field):
        for row in read_csv_rows(path):
            values[int(float(row["node_id"]))] = float(row["component0"])
    return values


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


def as_node_map(values: list[float], multiplier: float = 1.0) -> dict[int, float]:
    return {idx: value * multiplier for idx, value in enumerate(values)}


def common(reference: dict[int, float], candidate: dict[int, float]) -> list[int]:
    return sorted(set(reference) & set(candidate))


def value_stats(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    nodes = common(reference, candidate)
    diffs = [candidate[node] - reference[node] for node in nodes]
    abs_diffs = [abs(value) for value in diffs]
    rel_diffs = [
        abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
        for node in nodes
    ]
    return {
        "points": len(nodes),
        "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
        "median_abs_diff": statistics.median(abs_diffs) if abs_diffs else None,
        "rms_diff": math.sqrt(sum(value * value for value in diffs) / len(diffs)) if diffs else None,
        "max_abs_diff": max(abs_diffs) if abs_diffs else None,
        "max_rel_diff": max(rel_diffs) if rel_diffs else None,
    }


def log10_error_stats(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    values: list[float] = []
    for node in common(reference, candidate):
        ref = reference[node]
        cand = candidate[node]
        if ref > 0.0 and cand > 0.0:
            values.append(abs(math.log10(cand / ref)))
    ordered = sorted(values)
    p95 = ordered[int(round((len(ordered) - 1) * 0.95))] if ordered else None
    return {
        "points": len(values),
        "mean_log10_error": sum(values) / len(values) if values else None,
        "median_log10_error": statistics.median(values) if values else None,
        "p95_log10_error": p95,
        "max_log10_error": max(values) if values else None,
    }


def inferred_ni_from_np(electrons_cm3: dict[int, float],
                        holes_cm3: dict[int, float]) -> dict[int, float]:
    out: dict[int, float] = {}
    for node in common(electrons_cm3, holes_cm3):
        out[node] = math.sqrt(max(electrons_cm3[node] * holes_cm3[node], 0.0))
    return out


def median_value(values: dict[int, float]) -> float | None:
    raw = list(values.values())
    return statistics.median(raw) if raw else None


def resolve_deck_path(deck: dict[str, Any], reference_root: Path, key: str) -> None:
    if key not in deck:
        return
    candidate = Path(str(deck[key]))
    if not candidate.is_absolute():
        deck[key] = str((reference_root / "vela" / candidate).resolve())


def run_command(command: list[str], cwd: Path) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{result.returncode}: {' '.join(command)}\nstdout={result.stdout}\nstderr={result.stderr}"
        )


def parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return default
    out: list[float] = []
    for raw in values:
        out.extend(float(item) for item in raw.split(",") if item.strip())
    return out


def parse_string_list(values: list[str] | None, default: list[str]) -> list[str]:
    if not values:
        return default
    out: list[str] = []
    for raw in values:
        out.extend(item.strip() for item in raw.split(",") if item.strip())
    return list(dict.fromkeys(out))


def candidate_name(ni_cm3: float, bgn: str) -> str:
    return f"ni_{ni_cm3:.8g}_bgn_{bgn}".replace("+", "").replace(".", "p")


def run_candidate(reference_root: Path,
                  output_dir: Path,
                  runner: str,
                  ni_cm3: float,
                  bgn: str) -> Path:
    base_path = reference_root / "vela" / "simulation_0v.json"
    if not base_path.is_file():
        raise RuntimeError(f"missing base deck: {base_path}")
    deck = json.loads(json.dumps(read_json(base_path)))
    for key in ("mesh_file", "node_doping_file"):
        resolve_deck_path(deck, reference_root, key)
    name = candidate_name(ni_cm3, bgn)
    candidate_dir = output_dir / name
    candidate_dir.mkdir(parents=True, exist_ok=True)
    materials_path = candidate_dir / "materials.json"
    write_json(materials_path, {"materials": [{"name": "Si", "ni": ni_cm3}]})
    deck["materials_file"] = str(materials_path.resolve())
    solver = deck.setdefault("solver", {})
    if isinstance(solver, dict):
        solver["bandgap_narrowing"] = bgn
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
    deck_path = candidate_dir / "simulation_0v_probe.json"
    write_json(deck_path, deck)
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    matches = sorted(candidate_dir.glob("pn2d_0v_probe*.vtk"))
    if not matches:
        raise RuntimeError(f"candidate {name} did not produce VTK")
    return matches[0]


def evaluate_candidate(sentaurus: dict[str, dict[int, float]],
                       vtk_path: Path,
                       ni_cm3: float,
                       bgn: str) -> dict[str, Any]:
    vtk = parse_vtk_scalars(vtk_path)
    missing = [name for name in ["Potential", "Electrons", "Holes"] if name not in vtk]
    if missing:
        raise RuntimeError(f"candidate VTK missing fields: {', '.join(missing)}")
    psi = as_node_map(vtk["Potential"])
    n_cm3 = as_node_map(vtk["Electrons"], 1.0e-6)
    p_cm3 = as_node_map(vtk["Holes"], 1.0e-6)
    density_log = {
        "electron": log10_error_stats(sentaurus["eDensity"], n_cm3),
        "hole": log10_error_stats(sentaurus["hDensity"], p_cm3),
    }
    max_log = max(
        value for value in [
            density_log["electron"]["max_log10_error"],
            density_log["hole"]["max_log10_error"],
        ]
        if value is not None
    )
    s_ni = inferred_ni_from_np(sentaurus["eDensity"], sentaurus["hDensity"])
    v_ni = inferred_ni_from_np(n_cm3, p_cm3)
    s_ni_median = median_value(s_ni)
    v_ni_median = median_value(v_ni)
    ni_median_rel_error = (
        abs(v_ni_median - s_ni_median) / max(abs(s_ni_median), 1.0e-300)
        if s_ni_median is not None and v_ni_median is not None else None
    )
    potential = value_stats(sentaurus["ElectrostaticPotential"], psi)
    score = max_log
    if potential["rms_diff"] is not None:
        score += potential["rms_diff"]
    if ni_median_rel_error is not None:
        score += ni_median_rel_error
    return {
        "ni_cm3": ni_cm3,
        "bandgap_narrowing": bgn,
        "vtk": str(vtk_path),
        "score": score,
        "metrics": {
            "potential": potential,
            "density": {
                "max_log10_error": max_log,
                "electron": density_log["electron"],
                "hole": density_log["hole"],
            },
            "ni_eff_median_cm3": {
                "sentaurus": s_ni_median,
                "vela": v_ni_median,
                "relative_error": ni_median_rel_error,
            },
        },
    }


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V ni/BGN Probe",
        "",
        f"Status: {report.get('status')}",
    ]
    best = report.get("best_candidate", {})
    if best:
        lines.append(
            "Best: "
            f"ni={best.get('ni_cm3')} cm^-3, "
            f"BGN={best.get('bandgap_narrowing')}, "
            f"score={best.get('score')}"
        )
    lines.extend([
        "",
        "| ni cm^-3 | BGN | score | median density log10 err | p95 density log10 err | max density log10 err | potential RMS V | ni_eff rel err |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ])
    for candidate in report.get("candidates", []):
        metrics = candidate.get("metrics", {})
        density = metrics.get("density", {})
        electron = density.get("electron", {})
        hole = density.get("hole", {})
        median_log = max(
            value for value in [electron.get("median_log10_error"), hole.get("median_log10_error")]
            if value is not None
        )
        p95_log = max(
            value for value in [electron.get("p95_log10_error"), hole.get("p95_log10_error")]
            if value is not None
        )
        lines.append(
            f"| {candidate.get('ni_cm3')} | {candidate.get('bandgap_narrowing')} | "
            f"{candidate.get('score')} | "
            f"{median_log} | "
            f"{p95_log} | "
            f"{density.get('max_log10_error')} | "
            f"{metrics.get('potential', {}).get('rms_diff')} | "
            f"{metrics.get('ni_eff_median_cm3', {}).get('relative_error')} |"
        )
    if report.get("errors"):
        lines.extend(["", "Errors:"])
        lines.extend(f"- {item}" for item in report["errors"])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ni-cm3", action="append")
    parser.add_argument("--bgn", action="append")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME

    fields_root = args.reference_root / "sim_fields" / "0v" / "fields"
    sentaurus = {
        name: read_sentaurus_scalar(fields_root, name)
        for name in ["ElectrostaticPotential", "eDensity", "hDensity"]
    }
    missing = [name for name, values in sentaurus.items() if not values]
    if missing:
        report = {"status": "fail", "errors": [f"missing sentaurus field {name}" for name in missing]}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    ni_values = parse_float_list(args.ni_cm3, [1.0e10, 1.45e10, 1.6556207295e10])
    bgn_values = parse_string_list(args.bgn, ["none", "slotboom"])
    candidates: list[dict[str, Any]] = []
    errors: list[str] = []
    for ni_cm3 in ni_values:
        for bgn in bgn_values:
            try:
                vtk_path = run_candidate(args.reference_root, args.output_dir, args.runner, ni_cm3, bgn)
                candidates.append(evaluate_candidate(sentaurus, vtk_path, ni_cm3, bgn))
            except Exception as exc:
                errors.append(f"ni={ni_cm3:g}, bgn={bgn}: {exc}")

    candidates = sorted(candidates, key=lambda item: item["score"])
    report = {
        "status": "pass" if candidates else "fail",
        "reference_root": str(args.reference_root),
        "ni_candidates_cm3": ni_values,
        "bgn_candidates": bgn_values,
        "best_candidate": candidates[0] if candidates else None,
        "candidates": candidates,
        "errors": errors,
    }
    write_json(report_path, report)
    md_path.write_text(markdown_report(report))
    return 0 if candidates else 1


if __name__ == "__main__":
    raise SystemExit(main())
