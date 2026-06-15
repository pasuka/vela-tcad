#!/usr/bin/env python3
"""Decompose PN2D 0V density mismatch into units, formulas, and state differences."""

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


K_B_OVER_Q = 8.617333262145e-5
REPORT_NAME = "pn2d_0v_density_decomposition.json"
MARKDOWN_NAME = "pn2d_0v_density_decomposition.md"


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
    values: dict[int, float] = {}
    for row in read_csv_rows(path):
        node = int(float(row["node_id"]))
        donors = float(row.get("donors_cm3", 0.0) or 0.0)
        acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0)
        values[node] = donors - acceptors
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


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    return ordered[int(round((len(ordered) - 1) * p))]


def value_stats(reference: dict[int, float], candidate: dict[int, float]) -> dict[str, Any]:
    nodes = common(reference, candidate)
    diffs = [candidate[node] - reference[node] for node in nodes]
    abs_diffs = [abs(value) for value in diffs]
    rel_diffs = [
        abs(candidate[node] - reference[node]) / max(abs(reference[node]), 1.0e-300)
        for node in nodes
    ]
    worst = None
    if nodes:
        worst_node = max(nodes, key=lambda node: abs(candidate[node] - reference[node]))
        worst = {
            "node_id": worst_node,
            "reference": reference[worst_node],
            "candidate": candidate[worst_node],
            "abs_diff": abs(candidate[worst_node] - reference[worst_node]),
            "rel_diff": abs(candidate[worst_node] - reference[worst_node]) / max(abs(reference[worst_node]), 1.0e-300),
        }
    return {
        "points": len(nodes),
        "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
        "median_abs_diff": statistics.median(abs_diffs) if abs_diffs else None,
        "p95_abs_diff": quantile(abs_diffs, 0.95),
        "max_abs_diff": max(abs_diffs) if abs_diffs else None,
        "max_rel_diff": max(rel_diffs) if rel_diffs else None,
        "worst": worst,
    }


def distribution(values: dict[int, float]) -> dict[str, Any]:
    raw = list(values.values())
    return {
        "points": len(raw),
        "min": min(raw) if raw else None,
        "median": statistics.median(raw) if raw else None,
        "mean": sum(raw) / len(raw) if raw else None,
        "max": max(raw) if raw else None,
    }


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


def density_from_qf(psi: dict[int, float],
                    qf: dict[int, float],
                    doping_cm3: dict[int, float],
                    ni_cm3: float,
                    bgn: str,
                    vt: float,
                    carrier: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for node, psi_value in psi.items():
        if node not in qf:
            continue
        ni_node = ni_eff_cm3(ni_cm3, vt, doping_cm3.get(node, 0.0), bgn)
        if carrier == "electron":
            exponent = (psi_value - qf[node]) / vt
        else:
            exponent = (qf[node] - psi_value) / vt
        out[node] = ni_node * math.exp(max(min(exponent, 500.0), -500.0))
    return out


def inferred_ni_from_density(psi: dict[int, float],
                             qf: dict[int, float],
                             density_cm3: dict[int, float],
                             vt: float,
                             carrier: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for node, value in density_cm3.items():
        if node not in psi or node not in qf or value < 0.0:
            continue
        if carrier == "electron":
            exponent = -(psi[node] - qf[node]) / vt
        else:
            exponent = -(qf[node] - psi[node]) / vt
        out[node] = value * math.exp(max(min(exponent, 500.0), -500.0))
    return out


def inferred_ni_from_np(electrons_cm3: dict[int, float],
                        holes_cm3: dict[int, float]) -> dict[int, float]:
    out: dict[int, float] = {}
    for node in common(electrons_cm3, holes_cm3):
        out[node] = math.sqrt(max(electrons_cm3[node] * holes_cm3[node], 0.0))
    return out


def solver_bgn(config: dict[str, Any]) -> str:
    solver = config.get("solver", {})
    if not isinstance(solver, dict):
        return "none"
    value = solver.get("bandgap_narrowing", "none")
    if isinstance(value, dict):
        return str(value.get("model", "none"))
    return str(value)


def parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return default
    out: list[float] = []
    for raw in values:
        out.extend(float(item) for item in raw.split(",") if item.strip())
    return out


def best_formula(reference: dict[int, float],
                 psi: dict[int, float],
                 qf: dict[int, float],
                 doping_cm3: dict[int, float],
                 ni_values: list[float],
                 bgn_values: list[str],
                 vt: float,
                 carrier: str) -> dict[str, Any]:
    candidates: dict[str, Any] = {}
    for ni in ni_values:
        for bgn in bgn_values:
            values = density_from_qf(psi, qf, doping_cm3, ni, bgn, vt, carrier)
            key = f"ni_{ni:g}_bgn_{bgn}"
            candidates[key] = {
                "ni_cm3": ni,
                "bandgap_narrowing": bgn,
                "stats": value_stats(reference, values),
            }
    best_key = min(
        candidates,
        key=lambda key: (
            candidates[key]["stats"]["max_rel_diff"] if candidates[key]["stats"]["max_rel_diff"] is not None else math.inf,
            candidates[key]["stats"]["max_abs_diff"] if candidates[key]["stats"]["max_abs_diff"] is not None else math.inf,
        ),
    )
    return {"best_key": best_key, "best": candidates[best_key], "candidates": candidates}


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
    deck["output_csv"] = str((probe_dir / "pn2d_0v_density_decomposition.csv").resolve())
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = sweep.get("current_contact", "Anode")
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((probe_dir / "pn2d_0v_density_decomposition").resolve())
    deck_path = probe_dir / "simulation_0v_density_decomposition_probe.json"
    write_json(deck_path, deck)
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    vtk_matches = sorted(probe_dir.glob("pn2d_0v_density_decomposition*.vtk"))
    if not vtk_matches:
        raise RuntimeError("runner did not produce a 0V VTK probe")
    return vtk_matches[0]


def classify(report: dict[str, Any]) -> str:
    vela_self = report["vela_self_consistency"]
    sentaurus_vtk = report["sentaurus_vs_vela_density"]
    e_self = vela_self["electron"]["max_rel_diff"]
    h_self = vela_self["hole"]["max_rel_diff"]
    e_vs = sentaurus_vtk["electron_vtk_cm3"]["max_rel_diff"]
    h_vs = sentaurus_vtk["hole_vtk_cm3"]["max_rel_diff"]
    e_self = e_self if e_self is not None else math.inf
    h_self = h_self if h_self is not None else math.inf
    e_vs = e_vs if e_vs is not None else math.inf
    h_vs = h_vs if h_vs is not None else math.inf
    # ASCII VTK truncation leaves a small reconstruction mismatch even when the
    # solver density export and psi/QF formulas are internally consistent.
    if max(e_self, h_self) > 1.0e-4:
        return "vela_density_export_or_formula"
    if max(e_vs, h_vs) > 1.0e-2:
        return "vela_state_differs_from_sentaurus"
    return "density_parity"


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Density Decomposition",
        "",
        f"Status: {report.get('status')}",
        f"Classification: {report.get('classification')}",
        f"Thermal voltage: {report.get('thermal_voltage_V')}",
        "",
        "| Check | Electron max rel | Hole max rel |",
        "| --- | ---: | ---: |",
    ]
    checks = [
        ("Vela self consistency", "vela_self_consistency", "electron", "hole"),
        ("Sentaurus with primary ni", "sentaurus_self_consistency", "electron", "hole"),
        ("Sentaurus vs Vela VTK density", "sentaurus_vs_vela_density", "electron_vtk_cm3", "hole_vtk_cm3"),
        ("Sentaurus vs Vela formula", "sentaurus_vs_vela_formula", "electron_best", "hole_best"),
    ]
    for label, section, electron_key, hole_key in checks:
        item = report.get(section, {})
        electron = item.get(electron_key, {})
        hole = item.get(hole_key, {})
        if "stats" in electron:
            electron = electron["stats"]
        if "best" in electron:
            electron = electron["best"].get("stats", {})
        if "stats" in hole:
            hole = hole["stats"]
        if "best" in hole:
            hole = hole["best"].get("stats", {})
        lines.append(f"| {label} | {electron.get('max_rel_diff')} | {hole.get('max_rel_diff')} |")
    ni = report.get("inferred_ni_eff_cm3", {})
    lines.extend([
        "",
        f"Sentaurus sqrt(n*p) median ni_eff cm^-3: {ni.get('sentaurus_np', {}).get('median')}",
        f"Vela VTK sqrt(n*p) median ni_eff cm^-3: {ni.get('vela_vtk_np', {}).get('median')}",
    ])
    if report.get("missing_fields"):
        lines.extend(["", "Missing fields:"])
        lines.extend(f"- {item}" for item in report["missing_fields"])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--existing-vtk", type=Path)
    parser.add_argument("--runner")
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

    try:
        vtk_path = args.existing_vtk
        if vtk_path is None:
            if not args.runner:
                raise RuntimeError("--runner is required unless --existing-vtk is supplied")
            vtk_path = run_probe(args.reference_root, args.output_dir, args.runner)
        vtk = parse_vtk_scalars(vtk_path)
    except Exception as exc:
        report = {"status": "error", "error": str(exc)}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    required_sentaurus = [
        "ElectrostaticPotential",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eDensity",
        "hDensity",
    ]
    sentaurus = {name: read_sentaurus_scalar(fields_root, name) for name in required_sentaurus}
    missing = [f"sentaurus:{name}" for name, values in sentaurus.items() if not values]
    for name in ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"]:
        if name not in vtk:
            missing.append(f"vtk:{name}")
    if missing:
        report = {"status": "fail", "missing_fields": missing}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    vt = K_B_OVER_Q * args.temperature_k
    ni_values = parse_float_list(args.ni_cm3, [1.0e10, 1.45e10])
    primary_ni = ni_values[0]
    primary_bgn = solver_bgn(config)
    bgn_values = list(dict.fromkeys([primary_bgn, "none", "slotboom"]))
    doping_cm3 = read_doping(args.reference_root / "doping.csv")

    s_psi = sentaurus["ElectrostaticPotential"]
    s_qfn = sentaurus["eQuasiFermiPotential"]
    s_qfp = sentaurus["hQuasiFermiPotential"]
    s_n = sentaurus["eDensity"]
    s_p = sentaurus["hDensity"]
    v_psi = as_node_map(vtk["Potential"])
    v_qfn = as_node_map(vtk["ElectronQuasiFermi"])
    v_qfp = as_node_map(vtk["HoleQuasiFermi"])
    v_n_cm3 = as_node_map(vtk["Electrons"], 1.0e-6)
    v_p_cm3 = as_node_map(vtk["Holes"], 1.0e-6)

    v_n_formula = density_from_qf(v_psi, v_qfn, doping_cm3, primary_ni, primary_bgn, vt, "electron")
    v_p_formula = density_from_qf(v_psi, v_qfp, doping_cm3, primary_ni, primary_bgn, vt, "hole")
    s_n_formula = density_from_qf(s_psi, s_qfn, {}, primary_ni, "none", vt, "electron")
    s_p_formula = density_from_qf(s_psi, s_qfp, {}, primary_ni, "none", vt, "hole")

    formula_matrix = {
        "electron_best": best_formula(s_n, v_psi, v_qfn, doping_cm3, ni_values, bgn_values, vt, "electron"),
        "hole_best": best_formula(s_p, v_psi, v_qfp, doping_cm3, ni_values, bgn_values, vt, "hole"),
    }

    report: dict[str, Any] = {
        "status": "pass",
        "reference_root": str(args.reference_root),
        "vtk": str(vtk_path),
        "thermal_voltage_V": vt,
        "ni_candidates_cm3": ni_values,
        "primary_bandgap_narrowing": primary_bgn,
        "missing_fields": [],
        "sentaurus_vs_vela_density": {
            "electron_vtk_cm3": value_stats(s_n, v_n_cm3),
            "hole_vtk_cm3": value_stats(s_p, v_p_cm3),
        },
        "vela_self_consistency": {
            "electron": value_stats(v_n_cm3, v_n_formula),
            "hole": value_stats(v_p_cm3, v_p_formula),
        },
        "sentaurus_self_consistency": {
            "electron": value_stats(s_n, s_n_formula),
            "hole": value_stats(s_p, s_p_formula),
        },
        "sentaurus_vs_vela_formula": formula_matrix,
        "inferred_ni_eff_cm3": {
            "sentaurus_np": distribution(inferred_ni_from_np(s_n, s_p)),
            "sentaurus_e_density_psi_qf": distribution(inferred_ni_from_density(s_psi, s_qfn, s_n, vt, "electron")),
            "sentaurus_h_density_psi_qf": distribution(inferred_ni_from_density(s_psi, s_qfp, s_p, vt, "hole")),
            "vela_vtk_np": distribution(inferred_ni_from_np(v_n_cm3, v_p_cm3)),
            "vela_e_density_psi_qf": distribution(inferred_ni_from_density(v_psi, v_qfn, v_n_cm3, vt, "electron")),
            "vela_h_density_psi_qf": distribution(inferred_ni_from_density(v_psi, v_qfp, v_p_cm3, vt, "hole")),
        },
    }
    report["classification"] = classify(report)
    write_json(report_path, report)
    md_path.write_text(markdown_report(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
