#!/usr/bin/env python3
"""Compare Vela PN2D zero-bias state against Sentaurus 0V TDR exports."""

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
REPORT_NAME = "pn2d_0v_state_comparison.json"
MARKDOWN_NAME = "pn2d_0v_state_comparison.md"


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
    return nodes


def read_doping(path: Path) -> dict[int, float]:
    doping: dict[int, float] = {}
    for row in read_csv_rows(path):
        node_id = int(float(row["node_id"]))
        donors = float(row.get("donors_cm3", 0.0) or 0.0)
        acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0)
        doping[node_id] = donors - acceptors
    return doping


def read_contacts(path: Path) -> dict[str, list[int]]:
    if not path.is_file():
        return {}
    contacts: dict[str, list[int]] = {}
    for row in read_csv_rows(path):
        raw_ids = row.get("node_ids", "")
        contacts[row["name"]] = [int(item) for item in raw_ids.split(";") if item.strip()]
    return contacts


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


def as_node_map(values: list[float]) -> dict[int, float]:
    return {idx: value for idx, value in enumerate(values)}


def quantile(values: list[float], p: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    idx = int(round((len(ordered) - 1) * p))
    return ordered[idx]


def selected_node_groups(nodes: list[dict[str, float]],
                         contacts: dict[str, list[int]],
                         doping_cm3: dict[int, float]) -> dict[str, set[int]]:
    all_ids = {node["id"] for node in nodes}
    groups: dict[str, set[int]] = {"all": all_ids}
    if nodes:
        ys = [node["y_um"] for node in nodes]
        y_mid = 0.5 * (min(ys) + max(ys))
        y_span = max(ys) - min(ys)
        halfwidth = max(y_span * 0.02, 0.01)
        groups["centerline"] = {node["id"] for node in nodes if abs(node["y_um"] - y_mid) <= halfwidth}
        xs = [node["x_um"] for node in nodes]
        x_mid = 0.5 * (min(xs) + max(xs))
        x_halfwidth = max((max(xs) - min(xs)) * 0.02, 0.02)
        sign_change_nodes = {
            node_id for node_id, net in doping_cm3.items()
            if abs(net) <= 1.0e-30
        }
        groups["junction_near"] = sign_change_nodes or {
            node["id"] for node in nodes if abs(node["x_um"] - x_mid) <= x_halfwidth
        }
    for name, ids in contacts.items():
        groups[f"contact:{name}"] = set(ids)
    return groups


def compare_node_values(reference: dict[int, float],
                        candidate: dict[int, float],
                        groups: dict[str, set[int]]) -> dict[str, Any]:
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
            "rel_diff": abs(candidate[worst_node] - reference[worst_node]) / max(abs(reference[worst_node]), 1.0e-300),
        }

    group_stats: dict[str, Any] = {}
    for name, ids in groups.items():
        selected = [node for node in common if node in ids]
        selected_abs = [abs(candidate[node] - reference[node]) for node in selected]
        group_stats[name] = {
            "points": len(selected),
            "mean_abs_diff": sum(selected_abs) / len(selected_abs) if selected_abs else None,
            "max_abs_diff": max(selected_abs) if selected_abs else None,
        }

    return {
        "points_compared": len(common),
        "reference_points": len(reference),
        "candidate_points": len(candidate),
        "missing_candidate_nodes": len(set(reference) - set(candidate)),
        "missing_reference_nodes": len(set(candidate) - set(reference)),
        "mean_abs_diff": sum(abs_diffs) / len(abs_diffs) if abs_diffs else None,
        "median_abs_diff": statistics.median(abs_diffs) if abs_diffs else None,
        "p95_abs_diff": quantile(abs_diffs, 0.95),
        "max_abs_diff": max(abs_diffs) if abs_diffs else None,
        "max_rel_diff": max(rel_diffs) if rel_diffs else None,
        "worst": worst,
        "groups": group_stats,
    }


def slotboom_delta_eg_eV(net_doping_cm3: float,
                         n_cm3: float,
                         p_cm3: float,
                         reference_doping_cm3: float = 1.0e17,
                         coefficient_eV: float = 9.0e-3,
                         smoothing: float = 0.5) -> float:
    effective = max(abs(net_doping_cm3), n_cm3, p_cm3)
    if effective <= 0.0:
        return 0.0
    x = math.log(effective / reference_doping_cm3)
    return max(coefficient_eV * (x + math.sqrt(x * x + smoothing)), 0.0)


def ni_eff_cm3(ni_cm3: float,
               vt: float,
               net_doping_cm3: float,
               bgn: str) -> float:
    if bgn == "none":
        return ni_cm3
    delta = slotboom_delta_eg_eV(net_doping_cm3, 0.0, 0.0)
    return ni_cm3 * math.exp(delta / (2.0 * vt))


def density_from_qf(psi: dict[int, float],
                    qf: dict[int, float],
                    doping_cm3: dict[int, float],
                    ni_cm3: float,
                    bgn: str,
                    vt: float,
                    carrier: str) -> dict[int, float]:
    out: dict[int, float] = {}
    for node_id, psi_value in psi.items():
        qf_value = qf.get(node_id)
        if qf_value is None:
            continue
        ni_node = ni_eff_cm3(ni_cm3, vt, doping_cm3.get(node_id, 0.0), bgn)
        if carrier == "electron":
            exponent = max(min((psi_value - qf_value) / vt, 500.0), -500.0)
        else:
            exponent = max(min((qf_value - psi_value) / vt, 500.0), -500.0)
        out[node_id] = ni_node * math.exp(exponent)
    return out


def srh_rate_cm3_s(n_cm3: dict[int, float],
                   p_cm3: dict[int, float],
                   doping_cm3: dict[int, float],
                   ni_cm3: float,
                   bgn: str,
                   vt: float,
                   taun: float,
                   taup: float) -> dict[int, float]:
    out: dict[int, float] = {}
    for node_id, n_value in n_cm3.items():
        p_value = p_cm3.get(node_id)
        if p_value is None:
            continue
        ni_node = ni_eff_cm3(ni_cm3, vt, doping_cm3.get(node_id, 0.0), bgn)
        denominator = taup * (n_value + ni_node) + taun * (p_value + ni_node)
        excess = n_value * p_value - ni_node * ni_node
        if abs(excess) <= max(ni_node * ni_node, 1.0) * 1.0e-12:
            excess = 0.0
        out[node_id] = 0.0 if abs(denominator) < 1.0e-300 else excess / denominator
    return out


def solver_value(config: dict[str, Any], key: str, default: Any) -> Any:
    solver = config.get("solver", {})
    if isinstance(solver, dict):
        return solver.get(key, default)
    return default


def solver_bgn(config: dict[str, Any]) -> str:
    value = solver_value(config, "bandgap_narrowing", "none")
    if isinstance(value, str):
        return value
    if isinstance(value, dict):
        return str(value.get("model", "none"))
    return "none"


def derive_probe_deck(base: dict[str, Any],
                      output_dir: Path,
                      reference_root: Path,
                      contact: str) -> Path:
    deck = json.loads(json.dumps(base))
    deck["simulation_type"] = "dc_sweep"
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in deck:
            candidate = Path(str(deck[key]))
            if not candidate.is_absolute():
                deck[key] = str((reference_root / "vela" / candidate).resolve())
    deck["output_csv"] = f"pn2d_0v_{contact.lower()}_probe.csv"
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = contact
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = f"pn2d_0v_{contact.lower()}_probe"
    path = output_dir / "probes" / f"simulation_0v_{contact.lower()}_probe.json"
    write_json(path, deck)
    return path


def run_command(command: list[str], cwd: Path) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{result.returncode}: {' '.join(command)}\nstdout={result.stdout}\nstderr={result.stderr}"
        )


def run_probe(reference_root: Path,
              output_dir: Path,
              runner: str,
              contacts: list[str]) -> tuple[Path, Path]:
    base_path = reference_root / "vela" / "simulation_0v.json"
    base = read_json(base_path)
    probe_dir = output_dir / "probes"
    probe_dir.mkdir(parents=True, exist_ok=True)
    combined_terminal = output_dir / "terminal_currents.csv"
    rows: list[dict[str, str]] = []
    first_vtk: Path | None = None
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    for contact in contacts:
        deck_path = derive_probe_deck(base, output_dir, reference_root, contact)
        run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
        csv_path = deck_path.parent / f"pn2d_0v_{contact.lower()}_probe.csv"
        rows.extend(read_csv_rows(csv_path))
        vtk_matches = sorted(deck_path.parent.glob(f"pn2d_0v_{contact.lower()}_probe*.vtk"))
        if vtk_matches and first_vtk is None:
            first_vtk = vtk_matches[0]
    if not first_vtk:
        raise RuntimeError("runner did not produce a 0V VTK probe")
    if rows:
        with combined_terminal.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
    return first_vtk, combined_terminal


def terminal_report(
    path: Path,
    current_abs_gate_A_per_um: float,
    current_pair_relative_gate: float,
    current_pair_abs_gate_A_per_um: float,
) -> dict[str, Any]:
    rows = read_csv_rows(path)
    terminal_rows = [row for row in rows if row.get("current_contact")]
    by_contact: dict[str, dict[str, Any]] = {}
    total = 0.0
    max_abs = 0.0
    for row in terminal_rows:
        contact = row.get("current_contact", "")
        current_um = float(row.get("current_total_A_per_um", row.get("current_total", 0.0)) or 0.0)
        current_m = float(row.get("current_total", 0.0) or 0.0)
        total += current_um
        max_abs = max(max_abs, abs(current_um))
        by_contact[contact] = {
            "current_total_A_per_m": current_m,
            "current_total_A_per_um": current_um,
            "converged": row.get("converged"),
            "solver_method": row.get("solver_method"),
            "handoff_stage": row.get("handoff_stage"),
            "newton_iterations": int(float(row.get("newton_iterations", 0) or 0)),
        }
    pair_balance_relative = 0.0 if max_abs == 0.0 else abs(total) / max_abs
    failure_reasons: list[str] = []
    if not by_contact:
        failure_reasons.append("missing_terminal_currents")
    if max_abs > current_abs_gate_A_per_um:
        failure_reasons.append("terminal_current_above_abs_gate")
    if len(by_contact) != 2:
        failure_reasons.append("expected_two_terminal_currents")
    if (pair_balance_relative > current_pair_relative_gate and
            abs(total) > current_pair_abs_gate_A_per_um):
        failure_reasons.append("terminal_currents_not_equal_and_opposite")
    status = "pass"
    if failure_reasons:
        status = "fail"
    return {
        "status": status,
        "contacts": by_contact,
        "max_abs_A_per_um": max_abs,
        "sum_A_per_um": total,
        "gate_A_per_um": current_abs_gate_A_per_um,
        "pair_balance_relative": pair_balance_relative,
        "pair_relative_gate": current_pair_relative_gate,
        "pair_abs_gate_A_per_um": current_pair_abs_gate_A_per_um,
        "failure_reasons": failure_reasons,
    }


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V State Comparison",
        "",
        f"Status: {report['status']}",
        "",
        "| Field | Points | Max abs diff | Max rel diff |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name, stats in report.get("field_stats", {}).items():
        lines.append(
            f"| {name} | {stats.get('points_compared', 0)} | "
            f"{stats.get('max_abs_diff')} | {stats.get('max_rel_diff')} |"
        )
    terminal = report.get("terminal_currents", {})
    lines.extend([
        "",
        f"Terminal current status: {terminal.get('status')}",
        f"Terminal current sum A/um: {terminal.get('sum_A_per_um')}",
        f"Terminal current pair balance relative: {terminal.get('pair_balance_relative')}",
        "",
    ])
    if report.get("current_balance_report"):
        lines.extend([
            "Single-solution current balance report:",
            f"- {report['current_balance_report']}",
            "",
        ])
    if report.get("missing_fields"):
        lines.append("Missing fields:")
        for item in report["missing_fields"]:
            lines.append(f"- {item}")
    return "\n".join(lines) + "\n"


def parse_float_list(values: list[str] | None, default: list[float]) -> list[float]:
    if not values:
        return default
    out: list[float] = []
    for raw in values:
        out.extend(float(item) for item in raw.split(",") if item.strip())
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--existing-vtk", type=Path)
    parser.add_argument("--existing-terminal-csv", type=Path)
    parser.add_argument("--ni-cm3", action="append")
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--current-abs-gate-A-per-um", type=float, default=1.0e-12)
    parser.add_argument("--current-pair-relative-gate", type=float, default=5.0e-2)
    parser.add_argument("--current-pair-abs-gate-A-per-um", type=float, default=1.0e-24)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME
    config_path = args.reference_root / "vela" / "simulation_0v.json"
    config = read_json(config_path) if config_path.is_file() else {}

    missing: list[str] = []
    fields_root = args.reference_root / "sim_fields" / "0v" / "fields"
    sentaurus_names = [
        "ElectrostaticPotential",
        "eDensity",
        "hDensity",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "srhRecombination",
    ]
    sentaurus = {}
    for name in sentaurus_names:
        values = read_sentaurus_scalar(fields_root, name)
        if not values:
            missing.append(f"sentaurus:{name}")
        sentaurus[name] = values

    try:
        if args.existing_vtk and args.existing_terminal_csv:
            vtk_path = args.existing_vtk
            terminal_path = args.existing_terminal_csv
        else:
            if not args.runner:
                raise RuntimeError("--runner is required unless existing VTK/terminal CSV are supplied")
            vtk_path, terminal_path = run_probe(
                args.reference_root,
                args.output_dir,
                args.runner,
                ["Anode", "Cathode"],
            )
        vtk = parse_vtk_scalars(vtk_path)
    except Exception as exc:
        report = {"status": "fail", "missing_fields": missing, "error": str(exc)}
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1

    required_vtk = ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi"]
    for name in required_vtk:
        if name not in vtk:
            missing.append(f"vtk:{name}")

    nodes_path = args.reference_root / "nodes.csv"
    doping_path = args.reference_root / "doping.csv"
    contacts_path = args.reference_root / "contacts.csv"
    nodes = read_nodes(nodes_path) if nodes_path.is_file() else []
    doping_cm3 = read_doping(doping_path) if doping_path.is_file() else {}
    contacts = read_contacts(contacts_path)
    groups = selected_node_groups(nodes, contacts, doping_cm3)

    vt = K_B_OVER_Q * args.temperature_k
    ni_values = parse_float_list(args.ni_cm3, [1.0e10, 1.45e10])
    primary_ni = ni_values[0]
    primary_bgn = solver_bgn(config)
    bgn_values = list(dict.fromkeys([primary_bgn, "none", "slotboom"]))

    field_stats: dict[str, Any] = {}
    if not missing:
        psi = as_node_map(vtk["Potential"])
        phin = as_node_map(vtk["ElectronQuasiFermi"])
        phip = as_node_map(vtk["HoleQuasiFermi"])
        field_stats["ElectrostaticPotential"] = compare_node_values(
            sentaurus["ElectrostaticPotential"], psi, groups)
        field_stats["eQuasiFermiPotential"] = compare_node_values(
            sentaurus["eQuasiFermiPotential"], phin, groups)
        field_stats["hQuasiFermiPotential"] = compare_node_values(
            sentaurus["hQuasiFermiPotential"], phip, groups)
        e_formula = density_from_qf(psi, phin, doping_cm3, primary_ni, primary_bgn, vt, "electron")
        h_formula = density_from_qf(psi, phip, doping_cm3, primary_ni, primary_bgn, vt, "hole")
        field_stats["eDensity_formula"] = compare_node_values(sentaurus["eDensity"], e_formula, groups)
        field_stats["hDensity_formula"] = compare_node_values(sentaurus["hDensity"], h_formula, groups)
        if "Electrons" in vtk:
            field_stats["eDensity_vtk_raw_cm3"] = compare_node_values(
                sentaurus["eDensity"],
                {idx: value / 1.0e6 for idx, value in enumerate(vtk["Electrons"])},
                groups,
            )
        if "Holes" in vtk:
            field_stats["hDensity_vtk_raw_cm3"] = compare_node_values(
                sentaurus["hDensity"],
                {idx: value / 1.0e6 for idx, value in enumerate(vtk["Holes"])},
                groups,
            )
        taun = float(solver_value(config, "taun", 1.0e-6))
        taup = float(solver_value(config, "taup", 1.0e-6))
        srh = srh_rate_cm3_s(e_formula, h_formula, doping_cm3, primary_ni, primary_bgn, vt, taun, taup)
        field_stats["srhRecombination"] = compare_node_values(sentaurus["srhRecombination"], srh, groups)

    formula_matrix: dict[str, Any] = {}
    if not missing and "Potential" in vtk and "ElectronQuasiFermi" in vtk and "HoleQuasiFermi" in vtk:
        psi = as_node_map(vtk["Potential"])
        phin = as_node_map(vtk["ElectronQuasiFermi"])
        phip = as_node_map(vtk["HoleQuasiFermi"])
        for ni in ni_values:
            for bgn in bgn_values:
                e_values = density_from_qf(psi, phin, doping_cm3, ni, bgn, vt, "electron")
                h_values = density_from_qf(psi, phip, doping_cm3, ni, bgn, vt, "hole")
                key = f"ni_{ni:g}_bgn_{bgn}"
                formula_matrix[key] = {
                    "eDensity": compare_node_values(sentaurus["eDensity"], e_values, groups),
                    "hDensity": compare_node_values(sentaurus["hDensity"], h_values, groups),
                }

    terminal = terminal_report(
        terminal_path,
        args.current_abs_gate_A_per_um,
        args.current_pair_relative_gate,
        args.current_pair_abs_gate_A_per_um,
    )
    current_balance_report = (
        args.output_dir.parent
        / "0v_current_balance"
        / "pn2d_0v_current_balance.json"
    )
    status = "fail" if missing or terminal["status"] != "pass" else "pass"
    report = {
        "status": status,
        "reference_root": str(args.reference_root),
        "vtk": str(vtk_path),
        "current_balance_report": str(current_balance_report) if current_balance_report.is_file() else None,
        "missing_fields": sorted(dict.fromkeys(missing)),
        "thermal_voltage_V": vt,
        "diagnostic_matrix": {
            "ni_cm3": ni_values,
            "bandgap_narrowing": bgn_values,
            "primary": {"ni_cm3": primary_ni, "bandgap_narrowing": primary_bgn},
            "formula_variants": formula_matrix,
            "priorities": [
                "ni",
                "OldSlotboom/BGN",
                "Ohmic contact",
                "QF definitions",
                "carrier formulas",
                "current units",
            ],
        },
        "field_stats": field_stats,
        "terminal_currents": terminal,
    }
    write_json(report_path, report)
    md_path.write_text(markdown_report(report))
    return 0 if status == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
