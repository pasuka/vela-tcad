#!/usr/bin/env python3
"""Diagnose PN2D 0V contact/electrostatic state parity against Sentaurus.

Runs a strict coupled-Newton 0V probe (or consumes an existing VTK), then
compares contact-node and contact-adjacent ("first ring") electrostatic state
between Vela and Sentaurus to classify the dominant parity root cause.
"""

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


REPORT_NAME = "pn2d_0v_contact_state.json"
MARKDOWN_NAME = "pn2d_0v_contact_state.md"

K_B_OVER_Q = 8.617333262145e-5  # V/K

SENTAURUS_FIELDS = [
    "ElectrostaticPotential",
    "eQuasiFermiPotential",
    "hQuasiFermiPotential",
    "eDensity",
    "hDensity",
]
VTK_FIELDS = [
    "Potential",
    "ElectronQuasiFermi",
    "HoleQuasiFermi",
    "Electrons",
    "Holes",
]

ROOT_CAUSE_ORDER = [
    "contact_node_selection_mismatch",
    "contact_psi_pinning_mismatch",
    "qf_contact_boundary_mismatch",
    "built_in_potential_mismatch",
    "bulk_poisson_state_mismatch",
]


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


def parse_vtk_scalars(path: Path) -> tuple[int, dict[str, list[float]]]:
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
    return npoints, fields


def as_node_map(values: list[float], multiplier: float = 1.0) -> dict[int, float]:
    return {idx: value * multiplier for idx, value in enumerate(values)}


def read_doping(path: Path) -> dict[int, dict[str, float]]:
    doping: dict[int, dict[str, float]] = {}
    for row in read_csv_rows(path):
        node_id = int(float(row["node_id"]))
        donors = float(row.get("donors_cm3", 0.0) or 0.0)
        acceptors = float(row.get("acceptors_cm3", 0.0) or 0.0)
        doping[node_id] = {
            "donors_cm3": donors,
            "acceptors_cm3": acceptors,
            "net_cm3": donors - acceptors,
        }
    return doping


def read_contacts_csv(path: Path) -> dict[str, list[int]]:
    if not path.is_file():
        return {}
    contacts: dict[str, list[int]] = {}
    for row in read_csv_rows(path):
        raw_ids = row.get("node_ids", "")
        contacts[row["name"]] = [int(item) for item in raw_ids.split(";") if item.strip()]
    return contacts


def read_mesh(path: Path) -> dict[str, Any]:
    mesh = read_json(path)
    contacts: dict[str, list[int]] = {}
    for contact in mesh.get("contacts", []):
        contacts[str(contact.get("name"))] = [int(node) for node in contact.get("node_ids", [])]
    neighbors: dict[int, set[int]] = {}
    for triangle in mesh.get("triangles", []):
        node_ids = [int(node) for node in triangle.get("node_ids", [])]
        for a in node_ids:
            for b in node_ids:
                if a != b:
                    neighbors.setdefault(a, set()).add(b)
    return {
        "node_count": len(mesh.get("nodes", [])),
        "contacts": contacts,
        "neighbors": neighbors,
    }


def first_ring_nodes(contact_nodes: list[int],
                     all_contact_nodes: set[int],
                     neighbors: dict[int, set[int]]) -> list[int]:
    ring: set[int] = set()
    for node in contact_nodes:
        ring.update(neighbors.get(node, set()))
    return sorted(ring - all_contact_nodes)


def value_summary(node_map: dict[int, float], nodes: list[int]) -> dict[str, Any]:
    values = [node_map[node] for node in nodes if node in node_map]
    if not values:
        return {"points": 0, "min": None, "median": None, "max": None, "mean": None}
    return {
        "points": len(values),
        "min": min(values),
        "median": statistics.median(values),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def delta_summary(reference: dict[int, float],
                  candidate: dict[int, float],
                  nodes: list[int]) -> dict[str, Any]:
    diffs = [
        candidate[node] - reference[node]
        for node in nodes
        if node in candidate and node in reference
    ]
    if not diffs:
        return {"points": 0, "min": None, "median": None, "max": None, "rms": None}
    return {
        "points": len(diffs),
        "min": min(diffs),
        "median": statistics.median(diffs),
        "max": max(diffs),
        "rms": math.sqrt(sum(value * value for value in diffs) / len(diffs)),
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
    delta = slotboom_delta_eg_eV(net_doping_cm3)
    return ni_cm3 * math.exp(delta / (2.0 * vt))


def ohmic_state_summary(n_cm3: dict[int, float],
                        p_cm3: dict[int, float],
                        doping: dict[int, dict[str, float]],
                        nodes: list[int],
                        ni_cm3_value: float,
                        bgn: str,
                        vt: float) -> dict[str, Any]:
    neutrality: list[float] = []
    mass_action: list[float] = []
    for node in nodes:
        if node not in n_cm3 or node not in p_cm3 or node not in doping:
            continue
        net = doping[node]["net_cm3"]
        ni_eff = ni_eff_cm3(ni_cm3_value, vt, net, bgn)
        scale = max(abs(net), ni_eff)
        neutrality.append(abs(n_cm3[node] - p_cm3[node] - net) / scale)
        mass_action.append(abs(n_cm3[node] * p_cm3[node] - ni_eff * ni_eff) / (ni_eff * ni_eff))
    return {
        "points": len(neutrality),
        "neutrality_median_rel": statistics.median(neutrality) if neutrality else None,
        "neutrality_max_rel": max(neutrality) if neutrality else None,
        "mass_action_median_rel": statistics.median(mass_action) if mass_action else None,
        "mass_action_max_rel": max(mass_action) if mass_action else None,
    }


def plateau_nodes(doping: dict[int, dict[str, float]],
                  excluded: set[int],
                  side: str,
                  fraction: float = 0.25) -> list[int]:
    if side == "p":
        candidates = {node: -info["net_cm3"] for node, info in doping.items() if info["net_cm3"] < 0.0}
    else:
        candidates = {node: info["net_cm3"] for node, info in doping.items() if info["net_cm3"] > 0.0}
    if not candidates:
        return []
    peak = max(candidates.values())
    return sorted(
        node for node, magnitude in candidates.items()
        if magnitude >= fraction * peak and node not in excluded
    )


def potential_parity(sentaurus_psi: dict[int, float],
                     vela_psi: dict[int, float],
                     doping: dict[int, dict[str, float]]) -> dict[str, Any]:
    nodes = sorted(set(sentaurus_psi) & set(vela_psi))
    diffs = {node: vela_psi[node] - sentaurus_psi[node] for node in nodes}
    if not diffs:
        return {
            "points": 0,
            "rms_V": None,
            "global_shift_V": None,
            "rms_after_global_shift_V": None,
            "side_shifts_V": {"p_side": None, "n_side": None},
            "rms_after_side_shift_V": None,
        }
    values = list(diffs.values())
    rms = math.sqrt(sum(value * value for value in values) / len(values))
    global_shift = sum(values) / len(values)
    centered = [value - global_shift for value in values]
    rms_global = math.sqrt(sum(value * value for value in centered) / len(centered))

    side_groups: dict[str, list[float]] = {"p_side": [], "n_side": []}
    side_of: dict[int, str] = {}
    for node in nodes:
        net = doping.get(node, {}).get("net_cm3", 0.0)
        side = "p_side" if net < 0.0 else "n_side"
        side_of[node] = side
        side_groups[side].append(diffs[node])
    side_shifts = {
        side: (sum(group) / len(group) if group else 0.0)
        for side, group in side_groups.items()
    }
    side_residuals = [diffs[node] - side_shifts[side_of[node]] for node in nodes]
    rms_side = math.sqrt(sum(value * value for value in side_residuals) / len(side_residuals))
    return {
        "points": len(nodes),
        "rms_V": rms,
        "global_shift_V": global_shift,
        "rms_after_global_shift_V": rms_global,
        "side_shifts_V": {
            "p_side": side_shifts["p_side"] if side_groups["p_side"] else None,
            "n_side": side_shifts["n_side"] if side_groups["n_side"] else None,
        },
        "rms_after_side_shift_V": rms_side,
    }


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


def run_strict_newton_probe(reference_root: Path,
                            output_dir: Path,
                            runner: str,
                            ni_cm3_value: float,
                            bgn: str) -> tuple[Path, Path, Path]:
    base_path = reference_root / "vela" / "simulation_0v.json"
    if not base_path.is_file():
        raise RuntimeError(f"missing base deck: {base_path}")
    deck = json.loads(json.dumps(read_json(base_path)))
    for key in ("mesh_file", "node_doping_file"):
        resolve_deck_path(deck, reference_root, key)
    probe_dir = output_dir / "probe"
    probe_dir.mkdir(parents=True, exist_ok=True)
    materials_path = probe_dir / "materials.json"
    write_json(materials_path, {"materials": [{"name": "Si", "ni": ni_cm3_value}]})
    deck["materials_file"] = str(materials_path.resolve())
    solver = deck.setdefault("solver", {})
    if isinstance(solver, dict):
        solver["method"] = "newton"
        solver["bandgap_narrowing"] = bgn
        handoff = solver.setdefault("handoff", {})
        if isinstance(handoff, dict):
            handoff["fallback"] = "none"
    deck["simulation_type"] = "dc_sweep"
    deck["output_csv"] = str((probe_dir / "pn2d_0v_contact_state_probe.csv").resolve())
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = sweep.get("current_contact", "Anode")
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((probe_dir / "pn2d_0v_contact_state_probe").resolve())
    deck_path = probe_dir / "simulation_0v_contact_state_probe.json"
    write_json(deck_path, deck)
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    matches = sorted(probe_dir.glob("pn2d_0v_contact_state_probe*.vtk"))
    if not matches:
        raise RuntimeError("strict Newton probe did not produce VTK output")
    return matches[0], deck_path, materials_path


def contact_report(name: str,
                   contact_nodes: list[int],
                   ring_nodes: list[int],
                   sentaurus: dict[str, dict[int, float]],
                   vela: dict[str, dict[int, float]],
                   doping: dict[int, dict[str, float]],
                   ni_cm3_value: float,
                   bgn: str,
                   vt: float) -> dict[str, Any]:
    acceptors = [doping[node]["acceptors_cm3"] for node in contact_nodes if node in doping]
    donors = [doping[node]["donors_cm3"] for node in contact_nodes if node in doping]
    nets = [doping[node]["net_cm3"] for node in contact_nodes if node in doping]

    def side_state(psi: dict[int, float],
                   phin: dict[int, float],
                   phip: dict[int, float],
                   n_cm3: dict[int, float],
                   p_cm3: dict[int, float]) -> dict[str, Any]:
        contact_psi = value_summary(psi, contact_nodes)
        ring_psi = value_summary(psi, ring_nodes)
        delta = None
        if contact_psi["median"] is not None and ring_psi["median"] is not None:
            delta = ring_psi["median"] - contact_psi["median"]
        return {
            "contact_psi_V": contact_psi,
            "first_ring_psi_V": ring_psi,
            "contact_to_first_ring_delta_V": delta,
            "contact_phin_V": value_summary(phin, contact_nodes),
            "contact_phip_V": value_summary(phip, contact_nodes),
            "ohmic_state": ohmic_state_summary(
                n_cm3, p_cm3, doping, contact_nodes, ni_cm3_value, bgn, vt
            ),
        }

    return {
        "name": name,
        "contact_node_count": len(contact_nodes),
        "first_ring_node_count": len(ring_nodes),
        "doping": {
            "acceptors_median_cm3": statistics.median(acceptors) if acceptors else None,
            "donors_median_cm3": statistics.median(donors) if donors else None,
            "net_median_cm3": statistics.median(nets) if nets else None,
        },
        "sentaurus": side_state(
            sentaurus["psi"], sentaurus["phin"], sentaurus["phip"],
            sentaurus["n_cm3"], sentaurus["p_cm3"],
        ),
        "vela": side_state(
            vela["psi"], vela["phin"], vela["phip"],
            vela["n_cm3"], vela["p_cm3"],
        ),
        "delta_vela_minus_sentaurus": {
            "contact_psi_V": delta_summary(sentaurus["psi"], vela["psi"], contact_nodes),
            "first_ring_psi_V": delta_summary(sentaurus["psi"], vela["psi"], ring_nodes),
            "contact_phin_V": delta_summary(sentaurus["phin"], vela["phin"], contact_nodes),
            "contact_phip_V": delta_summary(sentaurus["phip"], vela["phip"], contact_nodes),
        },
    }


def classify(report: dict[str, Any],
             contact_delta_gate_v: float,
             qf_gate_v: float,
             built_in_gate_v: float,
             bulk_residual_floor_v: float) -> tuple[str, list[str], dict[str, bool]]:
    flags = {name: False for name in ROOT_CAUSE_ORDER}
    reasons: list[str] = []

    mesh = report["mesh"]
    for name, info in report["contacts"].items():
        if info["contact_node_count"] == 0:
            flags["contact_node_selection_mismatch"] = True
            reasons.append(f"{name}: empty contact node set in mesh.json")
    for name, status in mesh.get("contacts_csv_match", {}).items():
        if status is False:
            flags["contact_node_selection_mismatch"] = True
            reasons.append(f"{name}: mesh.json contact nodes differ from contacts.csv")
    if mesh.get("vtk_point_count") not in (None, mesh.get("node_count")):
        flags["contact_node_selection_mismatch"] = True
        reasons.append(
            f"VTK point count {mesh.get('vtk_point_count')} != mesh node count {mesh.get('node_count')}"
        )

    for name, info in report["contacts"].items():
        delta = info["delta_vela_minus_sentaurus"]["contact_psi_V"].get("median")
        if delta is not None and abs(delta) > contact_delta_gate_v:
            flags["contact_psi_pinning_mismatch"] = True
            reasons.append(
                f"{name}: contact psi median delta {delta:+.6g} V exceeds gate {contact_delta_gate_v:g} V"
            )
        for field in ("contact_phin_V", "contact_phip_V"):
            vela_median = info["vela"][field].get("median")
            qf_delta = info["delta_vela_minus_sentaurus"][field].get("median")
            if vela_median is not None and abs(vela_median) > qf_gate_v:
                flags["qf_contact_boundary_mismatch"] = True
                reasons.append(
                    f"{name}: Vela {field} median {vela_median:+.6g} V deviates from 0 V bias"
                )
            if qf_delta is not None and abs(qf_delta) > qf_gate_v:
                flags["qf_contact_boundary_mismatch"] = True
                reasons.append(
                    f"{name}: {field} Vela-Sentaurus median delta {qf_delta:+.6g} V exceeds gate"
                )

    built_in = report["built_in_potential"]
    estimate = built_in.get("estimate_V")
    s_delta = built_in.get("sentaurus_plateau_delta_V")
    v_delta = built_in.get("vela_plateau_delta_V")
    if s_delta is not None and v_delta is not None:
        if abs(v_delta - s_delta) > built_in_gate_v:
            flags["built_in_potential_mismatch"] = True
            reasons.append(
                f"plateau delta mismatch: Vela {v_delta:+.6g} V vs Sentaurus {s_delta:+.6g} V"
            )
        elif estimate is not None and abs(v_delta - estimate) > built_in_gate_v \
                and abs(s_delta - estimate) <= built_in_gate_v:
            flags["built_in_potential_mismatch"] = True
            reasons.append(
                f"Vela plateau delta {v_delta:+.6g} V deviates from built-in estimate {estimate:+.6g} V"
            )

    parity = report["potential_parity"]
    rms = parity.get("rms_V")
    rms_side = parity.get("rms_after_side_shift_V")
    if rms is not None and rms > contact_delta_gate_v and rms_side is not None:
        explained_gate = max(0.25 * rms, bulk_residual_floor_v)
        if rms_side > explained_gate:
            flags["bulk_poisson_state_mismatch"] = True
            reasons.append(
                f"potential RMS {rms:.6g} V not explained by contact-side offsets "
                f"(residual {rms_side:.6g} V > gate {explained_gate:.6g} V)"
            )

    for name in ROOT_CAUSE_ORDER:
        if flags[name]:
            return name, reasons, flags
    return "consistent", reasons, flags


def fmt(value: Any, precision: str = ".6g") -> str:
    if value is None:
        return "n/a"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return format(value, precision)
    return str(value)


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Contact/Electrostatic State Parity",
        "",
        f"Status: {report.get('status')}",
        f"Classification: {report.get('classification')}",
        "",
        f"- ni = {fmt(report.get('ni_cm3'))} cm^-3, BGN = {report.get('bandgap_narrowing')}, "
        f"Vt = {fmt(report.get('thermal_voltage_V'))} V",
        f"- VTK: {report.get('vtk')}",
        "",
        "## Contact-local psi parity (Vela - Sentaurus)",
        "",
        "| Contact | contact psi delta median V | contact psi delta rms V | "
        "first-ring psi delta median V | first-ring psi delta rms V |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for name, info in report.get("contacts", {}).items():
        delta = info["delta_vela_minus_sentaurus"]
        lines.append(
            f"| {name} | {fmt(delta['contact_psi_V'].get('median'))} | "
            f"{fmt(delta['contact_psi_V'].get('rms'))} | "
            f"{fmt(delta['first_ring_psi_V'].get('median'))} | "
            f"{fmt(delta['first_ring_psi_V'].get('rms'))} |"
        )
    lines.extend([
        "",
        "## Contact psi distributions",
        "",
        "| Contact | side | contact psi min/median/max V | first-ring psi min/median/max V | "
        "contact-to-first-ring delta V |",
        "| --- | --- | --- | --- | ---: |",
    ])
    for name, info in report.get("contacts", {}).items():
        for side in ("sentaurus", "vela"):
            state = info[side]
            contact_psi = state["contact_psi_V"]
            ring_psi = state["first_ring_psi_V"]
            lines.append(
                f"| {name} | {side} | "
                f"{fmt(contact_psi.get('min'))} / {fmt(contact_psi.get('median'))} / {fmt(contact_psi.get('max'))} | "
                f"{fmt(ring_psi.get('min'))} / {fmt(ring_psi.get('median'))} / {fmt(ring_psi.get('max'))} | "
                f"{fmt(state.get('contact_to_first_ring_delta_V'))} |"
            )
    built_in = report.get("built_in_potential", {})
    parity = report.get("potential_parity", {})
    lines.extend([
        "",
        "## Built-in potential",
        "",
        f"- Estimate Vt*ln(Na*Nd/ni_eff^2): {fmt(built_in.get('estimate_V'))} V "
        f"(Na={fmt(built_in.get('na_cm3'))} cm^-3, Nd={fmt(built_in.get('nd_cm3'))} cm^-3)",
        f"- Sentaurus quasi-neutral plateau delta (n-side minus p-side): "
        f"{fmt(built_in.get('sentaurus_plateau_delta_V'))} V",
        f"- Vela quasi-neutral plateau delta (n-side minus p-side): "
        f"{fmt(built_in.get('vela_plateau_delta_V'))} V",
        "",
        "## Potential RMS decomposition",
        "",
        f"- Global potential RMS (Vela - Sentaurus): {fmt(parity.get('rms_V'))} V",
        f"- After removing one global shift: {fmt(parity.get('rms_after_global_shift_V'))} V "
        f"(shift {fmt(parity.get('global_shift_V'))} V)",
        f"- After removing per-side (p/n) shifts: {fmt(parity.get('rms_after_side_shift_V'))} V "
        f"(p-side {fmt(parity.get('side_shifts_V', {}).get('p_side'))} V, "
        f"n-side {fmt(parity.get('side_shifts_V', {}).get('n_side'))} V)",
        f"- Contact-boundary offsets explain the potential RMS: "
        f"{fmt(report.get('contact_offset_explains_potential_rms'))}",
        "",
        "## Root cause flags",
        "",
    ])
    for name, value in report.get("root_cause_flags", {}).items():
        lines.append(f"- {name}: {fmt(value)}")
    if report.get("reasons"):
        lines.extend(["", "Reasons:"])
        lines.extend(f"- {item}" for item in report["reasons"])
    if report.get("errors"):
        lines.extend(["", "Errors:"])
        lines.extend(f"- {item}" for item in report["errors"])
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner", default=None)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--existing-vtk", type=Path, default=None)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--ni-cm3", type=float, default=1.6556207295e10)
    parser.add_argument("--bgn", default="none")
    parser.add_argument("--contact-delta-gate-v", type=float, default=0.02)
    parser.add_argument("--qf-gate-v", type=float, default=1.0e-3)
    parser.add_argument("--built-in-gate-v", type=float, default=0.02)
    parser.add_argument("--bulk-residual-floor-v", type=float, default=5.0e-3)
    parser.add_argument("--plateau-doping-fraction", type=float, default=0.25)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME
    vt = K_B_OVER_Q * args.temperature_k

    errors: list[str] = []
    report: dict[str, Any] = {
        "status": "fail",
        "reference_root": str(args.reference_root),
        "thermal_voltage_V": vt,
        "ni_cm3": args.ni_cm3,
        "bandgap_narrowing": args.bgn,
        "errors": errors,
    }

    try:
        if args.existing_vtk is not None:
            vtk_path = args.existing_vtk
            report["probe_deck"] = None
            report["materials_file"] = None
        else:
            if not args.runner:
                raise RuntimeError("either --runner or --existing-vtk is required")
            vtk_path, deck_path, materials_path = run_strict_newton_probe(
                args.reference_root, args.output_dir, args.runner, args.ni_cm3, args.bgn
            )
            report["probe_deck"] = str(deck_path)
            report["materials_file"] = str(materials_path)
        report["vtk"] = str(vtk_path)

        fields_root = args.reference_root / "sim_fields" / "0v" / "fields"
        sentaurus_raw = {
            name: read_sentaurus_scalar(fields_root, name) for name in SENTAURUS_FIELDS
        }
        missing = [name for name, values in sentaurus_raw.items() if not values]
        if missing:
            raise RuntimeError(f"missing sentaurus fields: {', '.join(missing)}")

        vtk_point_count, vtk = parse_vtk_scalars(vtk_path)
        missing_vtk = [name for name in VTK_FIELDS if name not in vtk]
        if missing_vtk:
            raise RuntimeError(f"missing VTK fields: {', '.join(missing_vtk)}")

        sentaurus = {
            "psi": sentaurus_raw["ElectrostaticPotential"],
            "phin": sentaurus_raw["eQuasiFermiPotential"],
            "phip": sentaurus_raw["hQuasiFermiPotential"],
            "n_cm3": sentaurus_raw["eDensity"],
            "p_cm3": sentaurus_raw["hDensity"],
        }
        vela = {
            "psi": as_node_map(vtk["Potential"]),
            "phin": as_node_map(vtk["ElectronQuasiFermi"]),
            "phip": as_node_map(vtk["HoleQuasiFermi"]),
            # Vela VTK densities are m^-3; convert to cm^-3 before comparison.
            "n_cm3": as_node_map(vtk["Electrons"], 1.0e-6),
            "p_cm3": as_node_map(vtk["Holes"], 1.0e-6),
        }

        mesh = read_mesh(args.reference_root / "vela" / "mesh.json")
        doping = read_doping(args.reference_root / "doping.csv")
        contacts_csv = read_contacts_csv(args.reference_root / "contacts.csv")

        contacts_csv_match: dict[str, bool] = {}
        for name, node_ids in mesh["contacts"].items():
            if name in contacts_csv:
                contacts_csv_match[name] = set(node_ids) == set(contacts_csv[name])

        all_contact_nodes: set[int] = set()
        for node_ids in mesh["contacts"].values():
            all_contact_nodes.update(node_ids)

        rings = {
            name: first_ring_nodes(node_ids, all_contact_nodes, mesh["neighbors"])
            for name, node_ids in mesh["contacts"].items()
        }

        report["mesh"] = {
            "node_count": mesh["node_count"],
            "vtk_point_count": vtk_point_count,
            "contact_node_counts": {
                name: len(node_ids) for name, node_ids in mesh["contacts"].items()
            },
            "first_ring_node_counts": {name: len(nodes) for name, nodes in rings.items()},
            "contacts_csv_match": contacts_csv_match,
        }

        report["contacts"] = {
            name: contact_report(
                name, node_ids, rings[name], sentaurus, vela,
                doping, args.ni_cm3, args.bgn, vt,
            )
            for name, node_ids in mesh["contacts"].items()
        }

        # Built-in potential estimate from contact doping and ni_eff.
        na = None
        nd = None
        for info in report["contacts"].values():
            net = info["doping"].get("net_median_cm3")
            if net is None:
                continue
            if net < 0.0:
                na = info["doping"].get("acceptors_median_cm3")
            else:
                nd = info["doping"].get("donors_median_cm3")
        estimate = None
        ni_eff_p = ni_eff_cm3(args.ni_cm3, vt, -(na or 0.0), args.bgn) if na else None
        ni_eff_n = ni_eff_cm3(args.ni_cm3, vt, nd or 0.0, args.bgn) if nd else None
        if na and nd and na > 0.0 and nd > 0.0:
            estimate = vt * math.log((na * nd) / (ni_eff_p * ni_eff_n))

        excluded = set(all_contact_nodes)
        for nodes in rings.values():
            excluded.update(nodes)
        p_nodes = plateau_nodes(doping, excluded, "p", args.plateau_doping_fraction)
        n_nodes = plateau_nodes(doping, excluded, "n", args.plateau_doping_fraction)

        def plateau_delta(psi: dict[int, float]) -> tuple[float | None, float | None, float | None]:
            p_summary = value_summary(psi, p_nodes)
            n_summary = value_summary(psi, n_nodes)
            delta = None
            if p_summary["median"] is not None and n_summary["median"] is not None:
                delta = n_summary["median"] - p_summary["median"]
            return p_summary["median"], n_summary["median"], delta

        s_p, s_n, s_delta = plateau_delta(sentaurus["psi"])
        v_p, v_n, v_delta = plateau_delta(vela["psi"])
        report["built_in_potential"] = {
            "estimate_V": estimate,
            "na_cm3": na,
            "nd_cm3": nd,
            "ni_eff_p_side_cm3": ni_eff_p,
            "ni_eff_n_side_cm3": ni_eff_n,
            "plateau_node_counts": {"p_side": len(p_nodes), "n_side": len(n_nodes)},
            "sentaurus_plateaus_V": {"p_side": s_p, "n_side": s_n},
            "vela_plateaus_V": {"p_side": v_p, "n_side": v_n},
            "sentaurus_plateau_delta_V": s_delta,
            "vela_plateau_delta_V": v_delta,
        }

        parity = potential_parity(sentaurus["psi"], vela["psi"], doping)
        report["potential_parity"] = parity
        rms = parity.get("rms_V")
        rms_side = parity.get("rms_after_side_shift_V")
        explains = None
        if rms is not None and rms > args.contact_delta_gate_v and rms_side is not None:
            explains = rms_side <= max(0.25 * rms, args.bulk_residual_floor_v)
        report["contact_offset_explains_potential_rms"] = explains

        classification, reasons, flags = classify(
            report,
            args.contact_delta_gate_v,
            args.qf_gate_v,
            args.built_in_gate_v,
            args.bulk_residual_floor_v,
        )
        report["classification"] = classification
        report["reasons"] = reasons
        report["root_cause_flags"] = flags
        report["status"] = "pass"
    except Exception as exc:
        errors.append(str(exc))
        report.setdefault("classification", "error")
        report.setdefault("root_cause_flags", {name: False for name in ROOT_CAUSE_ORDER})
        report.setdefault("reasons", [])

    write_json(report_path, report)
    md_path.write_text(markdown_report(report))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
