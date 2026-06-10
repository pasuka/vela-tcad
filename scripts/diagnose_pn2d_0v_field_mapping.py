#!/usr/bin/env python3
"""Diagnose whether PN2D 0V field mismatches come from node ordering or coordinates."""

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


REPORT_NAME = "pn2d_0v_field_mapping.json"
MARKDOWN_NAME = "pn2d_0v_field_mapping.md"
DEFAULT_FIELDS = [
    "ElectrostaticPotential:Potential",
    "eQuasiFermiPotential:ElectronQuasiFermi",
    "hQuasiFermiPotential:HoleQuasiFermi",
    "eDensity:Electrons",
    "hDensity:Holes",
]
UNIT_CONVENTIONS = {
    "raw": 1.0,
    "vtk_m3_to_cm3": 1.0e-6,
    "vtk_cm3_to_m3": 1.0e6,
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_nodes(path: Path) -> list[dict[str, float]]:
    nodes: list[dict[str, float]] = []
    for row in read_csv_rows(path):
        nodes.append({
            "id": int(float(row["id"])),
            "x_um": float(row.get("x_um", row.get("x", 0.0))),
            "y_um": float(row.get("y_um", row.get("y", 0.0))),
        })
    return sorted(nodes, key=lambda row: row["id"])


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


def parse_vtk(path: Path) -> tuple[list[tuple[float, float, float]], dict[str, list[float]]]:
    lines = path.read_text().splitlines()
    points: list[tuple[float, float, float]] = []
    fields: dict[str, list[float]] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 3 and parts[0] == "POINTS":
            npoints = int(parts[1])
            i += 1
            raw: list[float] = []
            while i < len(lines) and len(raw) < 3 * npoints:
                raw.extend(float(item) for item in lines[i].split())
                i += 1
            points = [
                (raw[3 * idx], raw[3 * idx + 1], raw[3 * idx + 2])
                for idx in range(npoints)
            ]
            continue
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < len(points):
                line = lines[i].strip()
                if not line:
                    i += 1
                    continue
                head = line.split()[0]
                if head in {"SCALARS", "CELL_DATA", "POINT_DATA", "VECTORS"}:
                    break
                values.extend(float(item) for item in line.split())
                i += 1
            fields[name] = values[:len(points)]
            continue
        i += 1
    return points, fields


def stats(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    common = sorted(set(reference) & set(candidate))
    diffs = [candidate[node] - reference[node] for node in common]
    abs_diffs = [abs(value) for value in diffs]
    rel_diffs = [
        abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
        for node in common
    ]
    worst = None
    if common:
        worst_node = max(common, key=lambda node: abs(candidate[node] - reference[node]))
        worst = {
            "node_id": worst_node,
            "reference": reference[worst_node],
            "candidate": candidate[worst_node],
            "abs_diff": abs(candidate[worst_node] - reference[worst_node]),
        }
    return {
        "matched_points": len(common),
        "max_abs_diff": max(abs_diffs) if abs_diffs else None,
        "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
        "rms_diff": math.sqrt(sum(value * value for value in diffs) / len(diffs)) if diffs else None,
        "max_rel_diff": max(rel_diffs) if rel_diffs else None,
        "worst": worst,
    }


def median(values: list[float]) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return 0.5 * (ordered[mid - 1] + ordered[mid])


def scaled_candidate(candidate: dict[int, float], multiplier: float) -> dict[int, float]:
    return {node: value * multiplier for node, value in candidate.items()}


def median_ratio(reference: dict[int, float], candidate: dict[int, float]) -> float | None:
    ratios: list[float] = []
    for node in sorted(set(reference) & set(candidate)):
        ref = reference[node]
        if abs(ref) <= 1.0e-300:
            continue
        ratios.append(candidate[node] / ref)
    return median(ratios)


def unit_convention_report(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    reports: dict[str, Any] = {}
    for name, multiplier in UNIT_CONVENTIONS.items():
        transformed = scaled_candidate(candidate, multiplier)
        entry = stats(reference, transformed)
        entry["candidate_multiplier"] = multiplier
        entry["median_candidate_over_reference"] = median_ratio(reference, transformed)
        reports[name] = entry
    best = min(
        reports,
        key=lambda name: (
            reports[name]["max_rel_diff"] if reports[name]["max_rel_diff"] is not None else math.inf,
            reports[name]["max_abs_diff"] if reports[name]["max_abs_diff"] is not None else math.inf,
        ),
    )
    raw_ratio = reports["raw"].get("median_candidate_over_reference")
    classification = "raw"
    if best == "vtk_m3_to_cm3" and raw_ratio is not None and 1.0e5 <= abs(raw_ratio) <= 1.0e7:
        classification = "vtk_m3_reference_cm3"
    elif best == "vtk_cm3_to_m3" and raw_ratio is not None and 1.0e-7 <= abs(raw_ratio) <= 1.0e-5:
        classification = "vtk_cm3_reference_m3"
    elif best != "raw":
        classification = "scaled"
    return {
        "best_unit_convention": best,
        "unit_classification": classification,
        "unit_conventions": reports,
    }


def should_check_density_units(sentaurus_name: str, vtk_name: str) -> bool:
    tokens = f"{sentaurus_name} {vtk_name}".lower()
    return any(
        marker in tokens
        for marker in ("density", "concentration", "electrons", "holes")
    )


def raw_unit_report(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    entry = stats(reference, candidate)
    entry["candidate_multiplier"] = 1.0
    entry["median_candidate_over_reference"] = median_ratio(reference, candidate)
    return {
        "best_unit_convention": "raw",
        "unit_classification": "raw",
        "unit_conventions": {"raw": entry},
    }


def coordinate_candidates(nodes: list[dict[str, float]],
                          points: list[tuple[float, float, float]]) -> dict[str, Any]:
    candidates = {
        "um": 1.0,
        "m_to_um": 1.0e6,
    }
    out: dict[str, Any] = {}
    for name, scale in candidates.items():
        mapping, distances = nearest_mapping(nodes, points, scale)
        out[name] = {
            "scale_to_um": scale,
            "mapping": mapping,
            "matched_points": len(mapping),
            "max_distance_um": max(distances) if distances else None,
            "mean_distance_um": sum(distances) / len(distances) if distances else None,
        }
    return out


def nearest_mapping(nodes: list[dict[str, float]],
                    points: list[tuple[float, float, float]],
                    scale_to_um: float) -> tuple[dict[int, int], list[float]]:
    available = set(range(len(points)))
    mapping: dict[int, int] = {}
    distances: list[float] = []
    for node in nodes:
        if not available:
            break
        best_index = min(
            available,
            key=lambda idx: (
                (points[idx][0] * scale_to_um - node["x_um"]) ** 2
                + (points[idx][1] * scale_to_um - node["y_um"]) ** 2
            ),
        )
        distance = math.sqrt(
            (points[best_index][0] * scale_to_um - node["x_um"]) ** 2
            + (points[best_index][1] * scale_to_um - node["y_um"]) ** 2
        )
        mapping[int(node["id"])] = best_index
        distances.append(distance)
        available.remove(best_index)
    return mapping, distances


def parse_field_specs(values: list[str] | None) -> list[tuple[str, str]]:
    raw_specs = values or DEFAULT_FIELDS
    specs: list[tuple[str, str]] = []
    for raw in raw_specs:
        for item in raw.split(","):
            item = item.strip()
            if not item:
                continue
            if ":" in item:
                sentaurus, vtk = item.split(":", 1)
            else:
                sentaurus = vtk = item
            specs.append((sentaurus, vtk))
    return specs


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


def run_probe(reference_root: Path, output_dir: Path, runner: str) -> Path:
    base = read_json(reference_root / "vela" / "simulation_0v.json")
    deck = json.loads(json.dumps(base))
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        resolve_deck_path(deck, reference_root, key)
    probe_dir = output_dir / "probe"
    probe_dir.mkdir(parents=True, exist_ok=True)
    deck["simulation_type"] = "dc_sweep"
    deck["output_csv"] = str((probe_dir / "pn2d_0v_field_mapping.csv").resolve())
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = sweep.get("current_contact", "Anode")
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((probe_dir / "pn2d_0v_field_mapping").resolve())
    deck_path = probe_dir / "simulation_0v_field_mapping_probe.json"
    write_json(deck_path, deck)
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    vtk_matches = sorted(probe_dir.glob("pn2d_0v_field_mapping*.vtk"))
    if not vtk_matches:
        raise RuntimeError("runner did not produce a 0V VTK probe")
    return vtk_matches[0]


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Field Mapping Diagnosis",
        "",
        f"Status: {report.get('status')}",
        f"Best coordinate scale: {report.get('coordinate_alignment', {}).get('best_scale')}",
        "",
        "| Field | Best pairing | Best unit | Unit class | Raw median ratio | Direct max abs | Best-unit max abs |",
        "| --- | --- | --- | --- | ---: | ---: | ---: |",
    ]
    for name, field in report.get("fields", {}).items():
        unit_name = field.get("best_unit_convention")
        unit_stats = field.get("unit_conventions", {}).get(unit_name, {})
        raw_stats = field.get("unit_conventions", {}).get("raw", {})
        lines.append(
            f"| {name} | {field.get('best_pairing')} | "
            f"{unit_name} | "
            f"{field.get('unit_classification')} | "
            f"{raw_stats.get('median_candidate_over_reference')} | "
            f"{field.get('direct_node_id', {}).get('max_abs_diff')} | "
            f"{unit_stats.get('max_abs_diff')} |"
        )
    if report.get("missing"):
        lines.extend(["", "Missing:"])
        for item in report["missing"]:
            lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--existing-vtk", type=Path)
    parser.add_argument("--runner")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--fields", action="append")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME

    try:
        vtk_path = args.existing_vtk
        if vtk_path is None:
            if not args.runner:
                raise RuntimeError("--runner is required unless --existing-vtk is supplied")
            vtk_path = run_probe(args.reference_root, args.output_dir, args.runner)
        points, vtk_fields = parse_vtk(vtk_path)
        nodes = read_nodes(args.reference_root / "nodes.csv")
    except Exception as exc:
        report = {"status": "error", "error": str(exc)}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    coordinate_reports = coordinate_candidates(nodes, points)
    best_scale = min(
        coordinate_reports,
        key=lambda name: (
            coordinate_reports[name]["max_distance_um"]
            if coordinate_reports[name]["max_distance_um"] is not None else math.inf,
            coordinate_reports[name]["mean_distance_um"]
            if coordinate_reports[name]["mean_distance_um"] is not None else math.inf,
        ),
    )
    nearest = coordinate_reports[best_scale]["mapping"]

    fields_root = args.reference_root / "sim_fields" / "0v" / "fields"
    field_reports: dict[str, Any] = {}
    missing: list[str] = []
    for sentaurus_name, vtk_name in parse_field_specs(args.fields):
        reference = read_sentaurus_scalar(fields_root, sentaurus_name)
        if not reference:
            missing.append(f"sentaurus:{sentaurus_name}")
            continue
        vtk_values = vtk_fields.get(vtk_name)
        if vtk_values is None:
            missing.append(f"vtk:{vtk_name}")
            continue
        direct = {
            node: vtk_values[node]
            for node in reference
            if 0 <= node < len(vtk_values)
        }
        nearest_values = {
            node: vtk_values[vtk_index]
            for node, vtk_index in nearest.items()
            if node in reference and 0 <= vtk_index < len(vtk_values)
        }
        direct_stats = stats(reference, direct)
        nearest_stats = stats(reference, nearest_values)
        direct_score = direct_stats["max_abs_diff"] if direct_stats["max_abs_diff"] is not None else math.inf
        nearest_score = nearest_stats["max_abs_diff"] if nearest_stats["max_abs_diff"] is not None else math.inf
        field_reports[sentaurus_name] = {
            "vtk_field": vtk_name,
            "best_pairing": "nearest_coordinate" if nearest_score < direct_score else "direct_node_id",
            "direct_node_id": direct_stats,
            "nearest_coordinate": nearest_stats,
            **(
                unit_convention_report(
                    reference,
                    nearest_values if nearest_score < direct_score else direct,
                )
                if should_check_density_units(sentaurus_name, vtk_name)
                else raw_unit_report(
                    reference,
                    nearest_values if nearest_score < direct_score else direct,
                )
            ),
        }

    report = {
        "status": "fail" if missing else "pass",
        "reference_root": str(args.reference_root),
        "vtk": str(vtk_path),
        "missing": missing,
        "coordinate_alignment": {
            "best_scale": best_scale,
            "candidates": {
                name: {
                    key: value for key, value in entry.items()
                    if key != "mapping"
                }
                for name, entry in coordinate_reports.items()
            },
        },
        "fields": field_reports,
    }
    write_json(report_path, report)
    md_path.write_text(markdown_report(report))
    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())
