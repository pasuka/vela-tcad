#!/usr/bin/env python3
"""Diagnose PN2D 0V potential, quasi-Fermi, and carrier-density conventions."""

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


K_B_OVER_Q = 8.617333262145e-5
REPORT_NAME = "pn2d_0v_field_conventions.json"
MARKDOWN_NAME = "pn2d_0v_field_conventions.md"


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


def read_doping(path: Path) -> dict[int, float]:
    if not path.is_file():
        return {}
    doping: dict[int, float] = {}
    for row in read_csv_rows(path):
        node_id = int(float(row["node_id"]))
        donors = float(row.get("donors_cm3", 0.0) or 0.0)
        acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0)
        doping[node_id] = donors - acceptors
    return doping


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


def as_node_map(values: list[float]) -> dict[int, float]:
    return {idx: value for idx, value in enumerate(values)}


def common_nodes(reference: dict[int, float], candidate: dict[int, float]) -> list[int]:
    return sorted(set(reference) & set(candidate))


def value_stats(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    nodes = common_nodes(reference, candidate)
    diffs = [candidate[node] - reference[node] for node in nodes]
    abs_diffs = [abs(value) for value in diffs]
    rel_diffs = [
        abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
        for node in nodes
    ]
    rms = math.sqrt(sum(value * value for value in diffs) / len(diffs)) if diffs else None
    return {
        "points": len(nodes),
        "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
        "rms_diff": rms,
        "max_abs_diff": max(abs_diffs) if abs_diffs else None,
        "max_rel_diff": max(rel_diffs) if rel_diffs else None,
    }


def mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def affine_fit(reference: dict[int, float], candidate: dict[int, float]) -> tuple[float, float]:
    nodes = common_nodes(reference, candidate)
    xs = [candidate[node] for node in nodes]
    ys = [reference[node] for node in nodes]
    xbar = mean(xs)
    ybar = mean(ys)
    denom = sum((x - xbar) ** 2 for x in xs)
    if denom <= 0.0:
        return 0.0, ybar
    slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denom
    offset = ybar - slope * xbar
    return slope, offset


def potential_convention(reference: dict[int, float], vela: dict[int, float]) -> dict[str, Any]:
    nodes = common_nodes(reference, vela)
    same_offset = mean([reference[node] - vela[node] for node in nodes])
    opposite_offset = mean([reference[node] + vela[node] for node in nodes])
    affine_slope, affine_offset = affine_fit(reference, vela)
    candidates: dict[str, dict[str, Any]] = {
        "same_sign": {"slope": 1.0, "offset": 0.0},
        "opposite_sign": {"slope": -1.0, "offset": 0.0},
        "same_sign_with_offset": {"slope": 1.0, "offset": same_offset},
        "opposite_sign_with_offset": {"slope": -1.0, "offset": opposite_offset},
        "affine": {"slope": affine_slope, "offset": affine_offset},
    }
    priority = {
        "same_sign": 0,
        "opposite_sign": 1,
        "same_sign_with_offset": 2,
        "opposite_sign_with_offset": 3,
        "affine": 4,
    }
    for name, entry in candidates.items():
        transformed = {
            node: entry["slope"] * vela[node] + entry["offset"]
            for node in nodes
        }
        entry["stats"] = value_stats(reference, transformed)
        entry["max_abs_diff"] = entry["stats"]["max_abs_diff"]
        entry["rms_diff"] = entry["stats"]["rms_diff"]
        entry["max_rel_diff"] = entry["stats"]["max_rel_diff"]
    best = min(
        candidates,
        key=lambda name: (
            candidates[name]["stats"]["max_abs_diff"] if candidates[name]["stats"]["max_abs_diff"] is not None else math.inf,
            priority[name],
        ),
    )
    return {"best": best, "candidates": candidates}


def slotboom_delta_eg_eV(net_doping_cm3: float,
                         reference_doping_cm3: float = 1.0e17,
                         coefficient_eV: float = 9.0e-3,
                         smoothing: float = 0.5) -> float:
    effective = abs(net_doping_cm3)
    if effective <= 0.0:
        return 0.0
    x = math.log(effective / reference_doping_cm3)
    return max(coefficient_eV * (x + math.sqrt(x * x + smoothing)), 0.0)


def ni_eff_cm3(ni_cm3: float, vt: float, net_doping_cm3: float, bgn: str) -> float:
    if bgn == "none":
        return ni_cm3
    return ni_cm3 * math.exp(slotboom_delta_eg_eV(net_doping_cm3) / (2.0 * vt))


def solver_bgn(config: dict[str, Any]) -> str:
    solver = config.get("solver", {})
    if not isinstance(solver, dict):
        return "none"
    value = solver.get("bandgap_narrowing", "none")
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return str(value.get("model", "none"))
    return "none"


def density_formula_values(psi: dict[int, float],
                           qf: dict[int, float],
                           doping_cm3: dict[int, float],
                           ni_cm3: float,
                           bgn: str,
                           vt: float,
                           formula: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for node, psi_value in psi.items():
        if node not in qf:
            continue
        qf_value = qf[node]
        ni_node = ni_eff_cm3(ni_cm3, vt, doping_cm3.get(node, 0.0), bgn)
        if formula == "psi_minus_qf":
            exponent = (psi_value - qf_value) / vt
        elif formula == "psi_plus_qf":
            exponent = (psi_value + qf_value) / vt
        elif formula == "qf_minus_psi":
            exponent = (qf_value - psi_value) / vt
        elif formula == "minus_qf_minus_psi":
            exponent = (-qf_value - psi_value) / vt
        else:
            raise ValueError(f"unknown formula {formula}")
        out[node] = ni_node * math.exp(max(min(exponent, 500.0), -500.0))
    return out


def carrier_formula_report(reference: dict[int, float],
                           psi: dict[int, float],
                           qf: dict[int, float],
                           doping_cm3: dict[int, float],
                           ni_values: list[float],
                           bgn_values: list[str],
                           vt: float,
                           formulas: list[str]) -> dict[str, Any]:
    candidates: dict[str, Any] = {}
    for ni in ni_values:
        for bgn in bgn_values:
            for formula in formulas:
                values = density_formula_values(psi, qf, doping_cm3, ni, bgn, vt, formula)
                key = f"ni_{ni:g}_bgn_{bgn}_{formula}"
                candidates[key] = {
                    "ni_cm3": ni,
                    "bandgap_narrowing": bgn,
                    "formula": formula,
                    "stats": value_stats(reference, values),
                }
    best_key = min(
        candidates,
        key=lambda key: candidates[key]["stats"]["max_rel_diff"]
        if candidates[key]["stats"]["max_rel_diff"] is not None else math.inf,
    )
    return {"best_key": best_key, "best": candidates[best_key], "candidates": candidates}


def parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return default
    out: list[float] = []
    for raw in values:
        out.extend(float(item) for item in raw.split(",") if item.strip())
    return out


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
    deck["output_csv"] = str((probe_dir / "pn2d_0v_field_conventions.csv").resolve())
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = sweep.get("current_contact", "Anode")
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((probe_dir / "pn2d_0v_field_conventions").resolve())
    deck_path = probe_dir / "simulation_0v_field_conventions_probe.json"
    write_json(deck_path, deck)
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    vtk_matches = sorted(probe_dir.glob("pn2d_0v_field_conventions*.vtk"))
    if not vtk_matches:
        raise RuntimeError("runner did not produce a 0V VTK probe")
    return vtk_matches[0]


def markdown_report(report: dict[str, Any]) -> str:
    potential = report.get("potential_convention", {})
    electron = report.get("carrier_formula", {}).get("electron", {}).get("best", {})
    hole = report.get("carrier_formula", {}).get("hole", {}).get("best", {})
    lines = [
        "# PN2D 0V Field Convention Diagnosis",
        "",
        f"Status: {report.get('status')}",
        f"Potential convention: {potential.get('best')}",
        f"Electron formula: {electron.get('formula')}",
        f"Hole formula: {hole.get('formula')}",
        "",
    ]
    if report.get("missing_fields"):
        lines.append("Missing fields:")
        for item in report["missing_fields"]:
            lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner")
    parser.add_argument("--existing-vtk", type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ni-cm3", action="append")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME
    fields_root = args.reference_root / "sim_fields" / "0v" / "fields"
    config_path = args.reference_root / "vela" / "simulation_0v.json"
    config = read_json(config_path) if config_path.is_file() else {}

    missing: list[str] = []
    sentaurus = {}
    for name in [
        "ElectrostaticPotential",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eDensity",
        "hDensity",
    ]:
        values = read_sentaurus_scalar(fields_root, name)
        if not values:
            missing.append(f"sentaurus:{name}")
        sentaurus[name] = values

    try:
        vtk_path = args.existing_vtk
        if vtk_path is None:
            if not args.runner:
                raise RuntimeError("--runner is required unless --existing-vtk is supplied")
            vtk_path = run_probe(args.reference_root, args.output_dir, args.runner)
        vtk = parse_vtk_scalars(vtk_path)
    except Exception as exc:
        report = {"status": "error", "missing_fields": missing, "error": str(exc)}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    required_vtk = ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi"]
    for name in required_vtk:
        if name not in vtk:
            missing.append(f"vtk:{name}")
    if missing:
        report = {"status": "fail", "missing_fields": sorted(dict.fromkeys(missing))}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    vt = K_B_OVER_Q * args.temperature_k
    ni_values = parse_float_list(args.ni_cm3, [1.0e10, 1.45e10])
    primary_bgn = solver_bgn(config)
    bgn_values = list(dict.fromkeys([primary_bgn, "none", "slotboom"]))
    doping_cm3 = read_doping(args.reference_root / "doping.csv")

    vela_potential = as_node_map(vtk["Potential"])
    potential = potential_convention(sentaurus["ElectrostaticPotential"], vela_potential)
    best_potential = potential["candidates"][potential["best"]]
    transformed_psi = {
        node: best_potential["slope"] * value + best_potential["offset"]
        for node, value in vela_potential.items()
    }

    qfn = as_node_map(vtk["ElectronQuasiFermi"])
    qfp = as_node_map(vtk["HoleQuasiFermi"])
    qf_conventions = {
        "electron": potential_convention(sentaurus["eQuasiFermiPotential"], qfn),
        "hole": potential_convention(sentaurus["hQuasiFermiPotential"], qfp),
    }
    electron_formula = carrier_formula_report(
        sentaurus["eDensity"],
        transformed_psi,
        qfn,
        doping_cm3,
        ni_values,
        bgn_values,
        vt,
        ["psi_minus_qf", "psi_plus_qf"],
    )
    hole_formula = carrier_formula_report(
        sentaurus["hDensity"],
        transformed_psi,
        qfp,
        doping_cm3,
        ni_values,
        bgn_values,
        vt,
        ["qf_minus_psi", "minus_qf_minus_psi"],
    )

    np_nodes = sorted(set(sentaurus["eDensity"]) & set(sentaurus["hDensity"]))
    sentaurus_ni_from_np = {
        node: math.sqrt(max(sentaurus["eDensity"][node] * sentaurus["hDensity"][node], 0.0))
        for node in np_nodes
    }
    ni_stats = {
        "points": len(sentaurus_ni_from_np),
        "min_cm3": min(sentaurus_ni_from_np.values()) if sentaurus_ni_from_np else None,
        "max_cm3": max(sentaurus_ni_from_np.values()) if sentaurus_ni_from_np else None,
        "mean_cm3": mean(list(sentaurus_ni_from_np.values())) if sentaurus_ni_from_np else None,
    }

    report = {
        "status": "pass",
        "reference_root": str(args.reference_root),
        "vtk": str(vtk_path),
        "missing_fields": [],
        "thermal_voltage_V": vt,
        "potential_convention": potential,
        "qf_convention": qf_conventions,
        "carrier_formula": {
            "electron": electron_formula,
            "hole": hole_formula,
        },
        "sentaurus_ni_from_np": ni_stats,
        "diagnostic_priorities": [
            "potential sign/offset",
            "QF sign/offset",
            "carrier formula",
            "ni/BGN",
        ],
    }
    write_json(report_path, report)
    md_path.write_text(markdown_report(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
