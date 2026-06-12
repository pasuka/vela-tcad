#!/usr/bin/env python3
"""Diagnose PN2D 0V terminal-current balance on one Vela Newton solution."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any


REPORT_NAME = "pn2d_0v_current_balance.json"
MARKDOWN_NAME = "pn2d_0v_current_balance.md"
DEFAULT_CONTACTS = ["Anode", "Cathode"]
SENTAURUS_CURRENT_FIELDS = {
    "eCurrent": "electron_current",
    "hCurrent": "hole_current",
    "TotalCurrent": "total_current",
}
DEFAULT_THERMAL_VOLTAGE_V = 0.025851999786435535
EXP_RESOLUTION_DELTA = 2.0e-16


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def float_value(row: dict[str, str], keys: list[str], default: float = 0.0) -> float:
    for key in keys:
        raw = row.get(key)
        if raw not in (None, ""):
            return float(raw)
    return default


def current_a_per_um(row: dict[str, str], base: str, default: float = 0.0) -> float:
    unit_key = f"{base}_A_per_um"
    if row.get(unit_key) not in (None, ""):
        return float(row[unit_key])
    if row.get(base) not in (None, ""):
        return float(row[base]) / 1.0e6
    return default


def optional_float(row: dict[str, str], key: str) -> float | None:
    raw = row.get(key)
    if raw in (None, ""):
        return None
    return float(raw)


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


def resolve_deck_path(deck: dict[str, Any], reference_root: Path, key: str) -> None:
    if key not in deck:
        return
    candidate = Path(str(deck[key]))
    if not candidate.is_absolute():
        deck[key] = str((reference_root / "vela" / candidate).resolve())


def derive_probe_deck(reference_root: Path, output_dir: Path) -> tuple[Path, Path, Path]:
    base_path = reference_root / "vela" / "simulation_0v.json"
    deck = read_json(base_path)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        resolve_deck_path(deck, reference_root, key)
    probe_dir = output_dir / "probe"
    probe_dir.mkdir(parents=True, exist_ok=True)
    terminal_csv = probe_dir / "terminal_balance.csv"
    edge_csv = probe_dir / "contact_edges.csv"
    deck["simulation_type"] = "dc_sweep"
    deck["output_csv"] = str((probe_dir / "pn2d_0v_current_balance.csv").resolve())
    sweep = deck.setdefault("sweep", {})
    sweep["mode"] = "iv"
    sweep["contact"] = sweep.get("contact", "Anode")
    sweep["current_contact"] = sweep.get("current_contact", "Anode")
    sweep["start"] = 0.0
    sweep["stop"] = 0.0
    sweep["step"] = 1.0
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((probe_dir / "pn2d_0v_current_balance").resolve())
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
    path = probe_dir / "simulation_0v_current_balance_probe.json"
    write_json(path, deck)
    return path, terminal_csv, edge_csv


def run_command(command: list[str], cwd: Path) -> None:
    result = subprocess.run(command, cwd=cwd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "command failed with exit code "
            f"{result.returncode}: {' '.join(command)}\nstdout={result.stdout}\nstderr={result.stderr}"
        )


def run_probe(reference_root: Path, output_dir: Path, runner: str) -> tuple[Path, Path, Path]:
    deck_path, terminal_csv, edge_csv = derive_probe_deck(reference_root, output_dir)
    runner_command = shlex.split(runner, posix=not sys.platform.startswith("win"))
    if runner_command:
        runner_path = Path(runner_command[0])
        if not runner_path.is_absolute():
            runner_command[0] = str((Path.cwd() / runner_path).resolve())
    run_command([*runner_command, "--config", str(deck_path.resolve())], cwd=deck_path.parent)
    vtk_matches = sorted(deck_path.parent.glob("pn2d_0v_current_balance*.vtk"))
    if not vtk_matches:
        raise RuntimeError("runner did not produce a 0V VTK probe")
    return vtk_matches[0], terminal_csv, edge_csv


def balance_metric(values: dict[str, float], relative_gate: float, abs_gate: float) -> dict[str, Any]:
    total = sum(values.values())
    max_abs = max((abs(value) for value in values.values()), default=0.0)
    relative = 0.0 if max_abs == 0.0 else abs(total) / max_abs
    absolute_floor_pass = abs(total) <= abs_gate
    relative_pass = relative <= relative_gate
    passed = relative_pass or absolute_floor_pass
    return {
        "by_contact_A_per_um": values,
        "sum_A_per_um": total,
        "abs_pair_sum_A_per_um": abs(total),
        "max_abs_A_per_um": max_abs,
        "pair_balance_relative": relative,
        "pair_balance_relative_gate": relative_gate,
        "pair_balance_abs_gate_A_per_um": abs_gate,
        "relative_gate_pass": relative_pass,
        "absolute_floor_gate_pass": absolute_floor_pass,
        "status": "pass" if passed else "fail",
    }


def terminal_report(path: Path, relative_gate: float, abs_gate: float) -> dict[str, Any]:
    rows = read_csv_rows(path)
    latest_by_contact: dict[str, dict[str, str]] = {}
    for row in rows:
        contact = row.get("contact") or row.get("current_contact") or ""
        if contact:
            latest_by_contact[contact] = row

    electron: dict[str, float] = {}
    hole: dict[str, float] = {}
    minus: dict[str, float] = {}
    plus: dict[str, float] = {}
    raw_m: dict[str, Any] = {}
    for contact, row in latest_by_contact.items():
        electron[contact] = current_a_per_um(row, "current_electron")
        hole[contact] = current_a_per_um(row, "current_hole")
        minus[contact] = current_a_per_um(
            row,
            "electron_minus_hole",
            current_a_per_um(row, "current_total"),
        )
        if row.get("electron_plus_hole_A_per_um") not in (None, "") or row.get("electron_plus_hole") not in (None, ""):
            plus[contact] = current_a_per_um(row, "electron_plus_hole")
        else:
            plus[contact] = electron[contact] + hole[contact]
        raw_m[contact] = {
            "current_electron_A_per_m": float_value(row, ["current_electron"], electron[contact] * 1.0e6),
            "current_hole_A_per_m": float_value(row, ["current_hole"], hole[contact] * 1.0e6),
            "electron_minus_hole_A_per_m": float_value(row, ["electron_minus_hole", "current_total"], minus[contact] * 1.0e6),
            "electron_plus_hole_A_per_m": float_value(row, ["electron_plus_hole"], plus[contact] * 1.0e6),
            "converged": row.get("converged"),
            "solver_method": row.get("solver_method"),
            "gummel_iterations": row.get("gummel_iterations"),
            "newton_iterations": row.get("newton_iterations"),
            "handoff_stage": row.get("handoff_stage"),
        }

    missing_contacts = [contact for contact in DEFAULT_CONTACTS if contact not in latest_by_contact]
    conventions = {
        "electron_minus_hole": balance_metric(minus, relative_gate, abs_gate),
        "electron_plus_hole": balance_metric(plus, relative_gate, abs_gate),
        "electron": balance_metric(electron, relative_gate, abs_gate),
        "hole": balance_metric(hole, relative_gate, abs_gate),
    }
    return {
        "path": str(path),
        "contacts": raw_m,
        "missing_contacts": missing_contacts,
        "conventions": conventions,
        "status": "pass" if not missing_contacts and conventions["electron_minus_hole"]["status"] == "pass" else "fail",
    }


def conservation_summary(terminal: dict[str, Any]) -> dict[str, Any]:
    conventions = terminal.get("conventions", {})
    minus = conventions.get("electron_minus_hole", {})
    plus = conventions.get("electron_plus_hole", {})
    electron = conventions.get("electron", {})
    hole = conventions.get("hole", {})
    return {
        "electron_minus_hole_sum_A_per_um": minus.get("sum_A_per_um"),
        "electron_plus_hole_sum_A_per_um": plus.get("sum_A_per_um"),
        "electron_sum_A_per_um": electron.get("sum_A_per_um"),
        "hole_sum_A_per_um": hole.get("sum_A_per_um"),
        "max_abs_terminal_current_A_per_um": minus.get("max_abs_A_per_um"),
        "abs_pair_sum_A_per_um": minus.get("abs_pair_sum_A_per_um"),
        "pair_balance_relative": minus.get("pair_balance_relative"),
        "pair_balance_relative_gate": minus.get("pair_balance_relative_gate"),
        "absolute_floor_gate_A_per_um": minus.get("pair_balance_abs_gate_A_per_um"),
        "relative_gate_pass": minus.get("relative_gate_pass"),
        "absolute_floor_gate_pass": minus.get("absolute_floor_gate_pass"),
        "status": minus.get("status"),
    }


def classify_zero_component_reason(
    row: dict[str, str],
    component: str,
    thermal_voltage: float = DEFAULT_THERMAL_VOLTAGE_V,
) -> dict[str, Any]:
    if component == "electron":
        branch = row.get("electron_branch")
        mobility = optional_float(row, "mun")
        carrier0 = optional_float(row, "n0")
        carrier1 = optional_float(row, "n1")
        qf0 = optional_float(row, "phin0")
        qf1 = optional_float(row, "phin1")
        flux = optional_float(row, "electron_continuity_flux")
    else:
        branch = row.get("hole_branch")
        mobility = optional_float(row, "mup")
        carrier0 = optional_float(row, "p0")
        carrier1 = optional_float(row, "p1")
        qf0 = optional_float(row, "phip0")
        qf1 = optional_float(row, "phip1")
        flux = optional_float(row, "hole_continuity_flux")

    qf_delta = None if qf0 is None or qf1 is None else qf0 - qf1
    qf_delta_over_vt = None if qf_delta is None else abs(qf_delta) / thermal_voltage
    mobility_nonzero = mobility is not None and abs(mobility) > 0.0
    carrier_nonzero = (
        carrier0 is not None and carrier1 is not None and
        (abs(carrier0) > 0.0 or abs(carrier1) > 0.0)
    )
    flux_zero = flux is None or abs(flux) == 0.0
    reason = "unknown"
    if not mobility_nonzero:
        reason = "zero_mobility"
    elif not carrier_nonzero:
        reason = "zero_carrier_density"
    elif branch == "quasi_fermi" and qf_delta_over_vt is not None and qf_delta_over_vt <= EXP_RESOLUTION_DELTA:
        reason = "qf_delta_below_exp_resolution"
    elif flux_zero:
        reason = "sg_flux_zero_after_evaluation"
    else:
        reason = "component_current_zero_after_scaling"
    return {
        "reason": reason,
        "branch": branch,
        "mobility": mobility,
        "mobility_nonzero": mobility_nonzero,
        "carrier0": carrier0,
        "carrier1": carrier1,
        "carrier_nonzero": carrier_nonzero,
        "qf_delta_V": qf_delta,
        "abs_qf_delta_over_Vt": qf_delta_over_vt,
        "exp_resolution_delta": EXP_RESOLUTION_DELTA,
        "continuity_flux": flux,
        "continuity_flux_zero": flux_zero,
    }


def zero_component_diagnostics(rows: list[dict[str, str]]) -> dict[str, Any]:
    by_contact: dict[str, dict[str, Any]] = {}
    examples: list[dict[str, Any]] = []
    for row in rows:
        contact = row.get("current_contact") or row.get("contact") or ""
        if not contact:
            continue
        contact_bucket = by_contact.setdefault(contact, {})
        for component in ("electron", "hole"):
            current = current_a_per_um(row, f"current_{component}")
            if current != 0.0:
                continue
            detail = classify_zero_component_reason(row, component)
            component_bucket = contact_bucket.setdefault(component, {
                "count": 0,
                "reason_counts": {},
                "dominant_reason": None,
                "max_abs_qf_delta_over_Vt": 0.0,
            })
            component_bucket["count"] += 1
            reason_counts = component_bucket["reason_counts"]
            reason_counts[detail["reason"]] = reason_counts.get(detail["reason"], 0) + 1
            if detail["abs_qf_delta_over_Vt"] is not None:
                component_bucket["max_abs_qf_delta_over_Vt"] = max(
                    component_bucket["max_abs_qf_delta_over_Vt"],
                    detail["abs_qf_delta_over_Vt"],
                )
            if len(examples) < 10:
                examples.append({
                    "contact": contact,
                    "component": component,
                    "edge_id": row.get("edge_id"),
                    **detail,
                })
    for contact_bucket in by_contact.values():
        for component_bucket in contact_bucket.values():
            reason_counts = component_bucket["reason_counts"]
            if reason_counts:
                component_bucket["dominant_reason"] = max(
                    reason_counts.items(),
                    key=lambda item: (item[1], item[0]),
                )[0]
    return {
        "thermal_voltage_V": DEFAULT_THERMAL_VOLTAGE_V,
        "exp_resolution_delta": EXP_RESOLUTION_DELTA,
        "by_contact": by_contact,
        "examples": examples,
    }


def edge_report(path: Path, terminal: dict[str, Any], top_n: int) -> dict[str, Any]:
    rows = read_csv_rows(path)
    sums: dict[str, float] = {}
    component_sums: dict[str, dict[str, float]] = {}
    counts: dict[str, int] = {}
    top_rows: list[dict[str, Any]] = []
    for row in rows:
        contact = row.get("current_contact") or row.get("contact") or ""
        if not contact:
            continue
        current = current_a_per_um(row, "current_total")
        electron_current = current_a_per_um(row, "current_electron")
        hole_current = current_a_per_um(row, "current_hole")
        contact_components = component_sums.setdefault(contact, {
            "electron": 0.0,
            "hole": 0.0,
            "electron_minus_hole": 0.0,
            "electron_plus_hole": 0.0,
        })
        contact_components["electron"] += electron_current
        contact_components["hole"] += hole_current
        contact_components["electron_minus_hole"] += current
        contact_components["electron_plus_hole"] += electron_current + hole_current
        sums[contact] = sums.get(contact, 0.0) + current
        counts[contact] = counts.get(contact, 0) + 1
        top_rows.append({
            "contact": contact,
            "edge_id": row.get("edge_id"),
            "node0": row.get("node0"),
            "node1": row.get("node1"),
            "current_total_A_per_um": current,
            "current_electron_A_per_um": electron_current,
            "current_hole_A_per_um": hole_current,
            "abs_current_total_A_per_um": abs(current),
            "psi0": row.get("psi0"),
            "psi1": row.get("psi1"),
            "phin0": row.get("phin0"),
            "phin1": row.get("phin1"),
            "phip0": row.get("phip0"),
            "phip1": row.get("phip1"),
            "electron_continuity_flux": row.get("electron_continuity_flux"),
            "hole_continuity_flux": row.get("hole_continuity_flux"),
        })
    top_rows.sort(key=lambda item: item["abs_current_total_A_per_um"], reverse=True)

    terminal_minus = terminal["conventions"]["electron_minus_hole"]["by_contact_A_per_um"]
    aggregation: dict[str, Any] = {}
    aggregation_fail = False
    for contact, terminal_value in terminal_minus.items():
        edge_sum = sums.get(contact, 0.0)
        abs_error = abs(edge_sum - terminal_value)
        scale = max(abs(edge_sum), abs(terminal_value), 1.0e-300)
        rel_error = abs_error / scale
        fail = rel_error > 1.0e-9 and abs_error > 1.0e-24
        aggregation_fail = aggregation_fail or fail
        aggregation[contact] = {
            "edge_sum_A_per_um": edge_sum,
            "terminal_electron_minus_hole_A_per_um": terminal_value,
            "abs_error_A_per_um": abs_error,
            "rel_error": rel_error,
            "status": "fail" if fail else "pass",
        }

    missing_contacts = [contact for contact in DEFAULT_CONTACTS if counts.get(contact, 0) == 0]
    return {
        "path": str(path),
        "flux_link_count_by_contact": counts,
        "edge_count_by_contact": counts,
        "flux_link_sum_A_per_um_by_contact": sums,
        "edge_sum_A_per_um_by_contact": sums,
        "component_sum_A_per_um_by_contact": component_sums,
        "zero_component_diagnostics": zero_component_diagnostics(rows),
        "missing_contacts": missing_contacts,
        "aggregation": aggregation,
        "aggregation_status": "fail" if aggregation_fail else "pass",
        "top_edges": top_rows[:top_n],
    }


def mesh_reference_report(reference_root: Path) -> dict[str, Any]:
    inventory = reference_root / "tdr_inventory" / "mesh.json"
    contact_boundary_counts: dict[str, int] = {}
    contact_node_counts: dict[str, int] = {}
    bulk_triangle_count: int | None = None
    tdr_total_element_count: int | None = None
    source = None
    if inventory.is_file():
        data = read_json(inventory)
        geometry = data.get("geometry", {})
        source = str(inventory)
        for region in geometry.get("regions", []):
            name = str(region.get("name", ""))
            region_type = int(region.get("type", 99))
            if region_type == 0:
                triangles = region.get("triangles", [])
                bulk_triangle_count = (bulk_triangle_count or 0) + len(triangles)
            elif region_type == 1:
                edges = region.get("edges", [])
                contact_boundary_counts[name] = len(edges)
                contact_node_counts[name] = len({int(node) for edge in edges for node in edge})
        if bulk_triangle_count is not None:
            tdr_total_element_count = bulk_triangle_count + sum(contact_boundary_counts.values())

    mesh_path = reference_root / "vela" / "mesh.json"
    vela_triangle_count: int | None = None
    vela_contact_node_counts: dict[str, int] = {}
    if mesh_path.is_file():
        mesh = read_json(mesh_path)
        vela_triangle_count = len(mesh.get("triangles", []))
        for contact in mesh.get("contacts", []):
            vela_contact_node_counts[str(contact.get("name", ""))] = len(contact.get("node_ids", []))

    return {
        "source": source,
        "sentaurus_total_elements_including_contacts": tdr_total_element_count,
        "sentaurus_bulk_triangle_elements": bulk_triangle_count,
        "sentaurus_contact_boundary_elements": contact_boundary_counts,
        "sentaurus_contact_boundary_node_counts": contact_node_counts,
        "vela_triangle_cells": vela_triangle_count,
        "vela_contact_node_counts": vela_contact_node_counts,
        "note": (
            "Sentaurus contact boundary elements are boundary segments from the grd/TDR. "
            "Vela contact-current flux_link_count counts finite-volume links from contact "
            "Dirichlet nodes into the bulk, so it may be node_count rather than segment_count."
        ),
    }


def vtk_report(path: Path) -> dict[str, Any]:
    fields = parse_vtk_scalars(path)
    required = ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi"]
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
        "near_equilibrium_qf": qf_max_span <= 1.0e-8,
    }


def sentaurus_log_reference(reference_root: Path) -> dict[str, Any]:
    repo = Path(__file__).resolve().parents[1]
    candidates = [
        reference_root / "source" / "pn2d_0v.log_des.log",
        repo / "reference_tcad" / "pn2d_sentaurus2018" / "source" / "pn2d_0v.log_des.log",
    ]
    for path in candidates:
        if path.is_file():
            lines = path.read_text(errors="replace").splitlines()
            snippets: dict[str, list[str]] = {}
            for start, stop in [(516, 518), (547, 549)]:
                snippets[f"{start}-{stop}"] = lines[start - 1:stop]
            return {
                "path": str(path),
                "line_ranges": ["516-518", "547-549"],
                "snippets": snippets,
                "expected_behavior": "Anode and Cathode conduction currents are equal magnitude and opposite sign at 0 V.",
            }
    return {
        "path": None,
        "line_ranges": ["516-518", "547-549"],
        "snippets": {},
        "expected_behavior": "Anode and Cathode conduction currents are equal magnitude and opposite sign at 0 V.",
    }


def parse_sentaurus_current_log(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    lines = path.read_text(errors="replace").splitlines()
    blocks: list[dict[str, Any]] = []
    for idx, line in enumerate(lines):
        if "electron current" not in line or "hole current" not in line or "conduction current" not in line:
            continue
        contacts: dict[str, dict[str, float]] = {}
        for contact_line in lines[idx + 1:idx + 3]:
            parts = contact_line.split()
            if len(parts) < 5:
                continue
            contacts[parts[0]] = {
                "voltage": float(parts[1]),
                "electron_current": float(parts[2]),
                "hole_current": float(parts[3]),
                "total_current": float(parts[4]),
            }
        if contacts:
            blocks.append({
                "source": "log",
                "line_range": f"{idx + 1}-{idx + 1 + len(contacts)}",
                "contacts": contacts,
            })
    if not blocks:
        return None
    return {
        "path": str(path),
        "initial_poisson": blocks[0],
        "final_coupled": blocks[-1],
        "blocks": blocks,
    }


def parse_quoted_tokens(block: str) -> list[str]:
    return re.findall(r'"([^"]+)"', block)


def sentaurus_plt_records(path: Path) -> list[dict[str, float]]:
    text = path.read_text(errors="replace")
    datasets_match = re.search(r"datasets\s*=\s*\[(.*?)\]", text, flags=re.S)
    data_match = re.search(r"Data\s*\{(.*?)\}", text, flags=re.S)
    if not datasets_match or not data_match:
        return []
    datasets = parse_quoted_tokens(datasets_match.group(1))
    if not datasets:
        return []
    values = [float(item) for item in re.findall(r"[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][+-]?\d+)?", data_match.group(1))]
    width = len(datasets)
    if width == 0 or len(values) < width:
        return []
    records: list[dict[str, float]] = []
    for offset in range(0, len(values) - width + 1, width):
        chunk = values[offset:offset + width]
        if len(chunk) == width:
            records.append(dict(zip(datasets, chunk)))
    return records


def convert_sentaurus_plt_record(record: dict[str, float]) -> dict[str, Any]:
    contacts: dict[str, dict[str, float]] = {}
    for dataset, value in record.items():
        if dataset == "time":
            continue
        parts = dataset.split(" ", 1)
        if len(parts) != 2:
            continue
        contact, field = parts
        contact_report = contacts.setdefault(contact, {})
        if field in SENTAURUS_CURRENT_FIELDS:
            contact_report[SENTAURUS_CURRENT_FIELDS[field]] = value
        elif field == "OuterVoltage":
            contact_report["voltage"] = value
    return contacts


def parse_sentaurus_current_plt(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    records = sentaurus_plt_records(path)
    current_blocks: list[dict[str, Any]] = []
    for index, record in enumerate(records):
        contacts = convert_sentaurus_plt_record(record)
        if contacts:
            current_blocks.append({
                "source": "plt",
                "record_index": index,
                "contacts": contacts,
            })
    if not current_blocks:
        return None
    return {
        "path": str(path),
        "initial_poisson": current_blocks[0],
        "final_coupled": current_blocks[-1],
        "blocks": current_blocks,
    }


def sentaurus_current_reference(reference_root: Path) -> dict[str, Any]:
    repo = Path(__file__).resolve().parents[1]
    source_dirs = [
        reference_root / "source",
        repo / "reference_tcad" / "pn2d_sentaurus2018" / "source",
    ]
    plt_reference = None
    log_reference = None
    for source_dir in source_dirs:
        if plt_reference is None:
            plt_reference = parse_sentaurus_current_plt(source_dir / "pn2d_0v.plt")
        if log_reference is None:
            log_reference = parse_sentaurus_current_log(source_dir / "pn2d_0v.log_des.log")
    selected = plt_reference or log_reference
    if not selected:
        return {
            "status": "missing",
            "source_preference": ["plt", "log"],
            "searched_source_dirs": [str(path) for path in source_dirs],
            "initial_poisson": None,
            "final_coupled": None,
        }
    selected["status"] = "present"
    selected["source_preference"] = ["plt", "log"]
    selected["searched_source_dirs"] = [str(path) for path in source_dirs]
    if log_reference and plt_reference:
        selected["log_crosscheck"] = {
            "path": log_reference.get("path"),
            "final_coupled": log_reference.get("final_coupled"),
        }
    return selected


def signed_relation(a: float, b: float) -> str:
    if a == 0.0 or b == 0.0:
        return "zero"
    return "same_sign" if math.copysign(1.0, a) == math.copysign(1.0, b) else "opposite_sign"


def current_component_presence(
    terminal: dict[str, Any],
    sentaurus_final: dict[str, Any] | None,
    zero_gate: float,
) -> dict[str, Any]:
    sentaurus_contacts = (sentaurus_final or {}).get("contacts", {})
    missing_components: list[str] = []
    by_contact: dict[str, Any] = {}
    for contact in DEFAULT_CONTACTS:
        vela = terminal.get("contacts", {}).get(contact, {})
        sentaurus = sentaurus_contacts.get(contact, {})
        contact_components: dict[str, Any] = {}
        for label, vela_key, sentaurus_key in [
            ("electron", "current_electron_A_per_m", "electron_current"),
            ("hole", "current_hole_A_per_m", "hole_current"),
        ]:
            vela_value = float(vela.get(vela_key) or 0.0)
            sentaurus_value = float(sentaurus.get(sentaurus_key) or 0.0)
            vela_zero = abs(vela_value) <= zero_gate
            sentaurus_nonzero = abs(sentaurus_value) > zero_gate
            if vela_zero and sentaurus_nonzero:
                missing_components.append(f"{contact}:{label}")
            contact_components[label] = {
                "vela_A_per_m": vela_value,
                "sentaurus_current": sentaurus_value,
                "vela_zero": vela_zero,
                "sentaurus_nonzero": sentaurus_nonzero,
                "status": "missing_in_vela" if vela_zero and sentaurus_nonzero else "present_or_both_zero",
            }
        by_contact[contact] = contact_components
    return {
        "zero_gate": zero_gate,
        "by_contact": by_contact,
        "missing_components": missing_components,
        "status": "missing_components" if missing_components else "ok",
    }


def sentaurus_current_parity(
    terminal: dict[str, Any],
    reference: dict[str, Any],
    ratio_gate: float = 10.0,
    zero_gate: float = 1.0e-40,
) -> dict[str, Any]:
    final = reference.get("final_coupled")
    if not final:
        return {
            "status": "missing_reference",
            "classification": "sentaurus_current_reference_missing",
        }
    sentaurus_contacts = final.get("contacts", {})
    vela_minus = terminal.get("conventions", {}).get("electron_minus_hole", {}).get("by_contact_A_per_um", {})
    by_contact: dict[str, Any] = {}
    ratios: list[float] = []
    signs: dict[str, str] = {}
    for contact in DEFAULT_CONTACTS:
        vela_current = float(vela_minus.get(contact) or 0.0)
        sentaurus_total = float(sentaurus_contacts.get(contact, {}).get("total_current") or 0.0)
        ratio = math.inf if abs(vela_current) <= zero_gate and abs(sentaurus_total) > zero_gate else (
            abs(sentaurus_total) / abs(vela_current) if vela_current != 0.0 else 0.0
        )
        if math.isfinite(ratio):
            ratios.append(ratio)
        signs[contact] = signed_relation(vela_current, sentaurus_total)
        by_contact[contact] = {
            "vela_current_total_A_per_um": vela_current,
            "sentaurus_total_current": sentaurus_total,
            "abs_ratio_sentaurus_to_vela": ratio,
            "sign_relation": signs[contact],
        }
    component_presence = current_component_presence(terminal, final, zero_gate)
    finite_mismatch = any(ratio > ratio_gate or ratio < 1.0 / ratio_gate for ratio in ratios)
    infinite_mismatch = any(
        math.isinf(item["abs_ratio_sentaurus_to_vela"]) for item in by_contact.values()
    )
    sign_mismatch = any(relation == "opposite_sign" for relation in signs.values())
    status = "mismatch" if finite_mismatch or infinite_mismatch or sign_mismatch else "match"
    classification = "sentaurus_current_parity_match"
    if status == "mismatch":
        classification = "sentaurus_abs_current_mismatch"
    return {
        "status": status,
        "classification": classification,
        "unit_hypothesis": "vela_current_total_A_per_um_vs_sentaurus_current",
        "ratio_gate": ratio_gate,
        "zero_gate": zero_gate,
        "by_contact": by_contact,
        "component_presence": component_presence,
        "minority_component_zero_contacts": sorted({
            entry.split(":", 1)[0] for entry in component_presence["missing_components"]
        }),
        "notes": [
            "This compares Vela electron_minus_hole_A_per_um with Sentaurus TotalCurrent as printed in pn2d_0v.plt/log.",
            "A mismatch here does not imply two-terminal imbalance; it flags absolute-current or definition parity.",
        ],
    }


def contact_metadata_indicates_weak_solve(terminal: dict[str, Any]) -> bool:
    for row in terminal.get("contacts", {}).values():
        converged = str(row.get("converged", "")).strip().lower()
        if converged in {"0", "false", "no"}:
            return True
        handoff = str(row.get("handoff_stage", "")).strip().lower()
        if "failed" in handoff or handoff in {"", "non_convergence"}:
            return True
        try:
            iterations = int(float(row.get("newton_iterations") or 0))
        except (TypeError, ValueError):
            iterations = 0
        if iterations <= 0:
            return True
    return False


def classify(terminal: dict[str, Any], edges: dict[str, Any], vtk: dict[str, Any]) -> tuple[str, list[str], dict[str, bool]]:
    minus = terminal["conventions"]["electron_minus_hole"]
    plus = terminal["conventions"]["electron_plus_hole"]
    flags = {
        "two_probe_artifact": False,
        "total_current_sign_convention": False,
        "contact_edge_coverage": False,
        "contact_current_aggregation": False,
        "contact_flux_formula": False,
        "newton_residual_not_tight_enough": False,
        "contact_boundary_qf_state": False,
    }
    reasons: list[str] = []

    if terminal["missing_contacts"]:
        flags["contact_edge_coverage"] = True
        reasons.append(f"missing terminal balance contacts: {terminal['missing_contacts']}")
    if edges["missing_contacts"]:
        flags["contact_edge_coverage"] = True
        reasons.append(f"missing contact-edge rows: {edges['missing_contacts']}")
    if flags["contact_edge_coverage"]:
        return "contact_edge_coverage", reasons, flags

    if edges["aggregation_status"] == "fail":
        flags["contact_current_aggregation"] = True
        reasons.append("edge current sum does not match terminal electron_minus_hole row")
        return "contact_current_aggregation", reasons, flags

    if minus["status"] == "fail" and plus["status"] == "pass":
        flags["total_current_sign_convention"] = True
        reasons.append("electron_plus_hole balances but electron_minus_hole does not")
        return "total_current_sign_convention", reasons, flags

    if minus["status"] == "fail" and contact_metadata_indicates_weak_solve(terminal):
        flags["newton_residual_not_tight_enough"] = True
        reasons.append("terminal current is above the absolute floor and solver metadata indicates weak/non-Newton convergence")
        return "newton_residual_not_tight_enough", reasons, flags

    if vtk["missing_fields"]:
        flags["contact_boundary_qf_state"] = True
        reasons.append(f"missing VTK fields for QF equilibrium check: {vtk['missing_fields']}")
        return "contact_boundary_qf_state", reasons, flags

    if not vtk["near_equilibrium_qf"]:
        flags["contact_boundary_qf_state"] = True
        reasons.append(f"quasi-Fermi span is {vtk['qf_max_span_V']} V at 0 V")
        return "contact_boundary_qf_state", reasons, flags

    if minus["status"] == "fail":
        flags["contact_flux_formula"] = True
        reasons.append("same-solution terminal current does not balance under either total-current convention")
        return "contact_flux_formula", reasons, flags

    reasons.append("same-solution electron_minus_hole terminal current balances")
    return "balanced", reasons, flags


def markdown_report(report: dict[str, Any]) -> str:
    lines = [
        "# PN2D 0V Current Balance Diagnosis",
        "",
        f"Classification: {report.get('classification')}",
        "",
        "## Terminal Balance",
        "",
        "| Convention | Status | Sum A/um | Relative balance |",
        "| --- | --- | ---: | ---: |",
    ]
    for name, stats in report.get("terminal_balance", {}).get("conventions", {}).items():
        lines.append(
            f"| {name} | {stats.get('status')} | {stats.get('sum_A_per_um')} | "
            f"{stats.get('pair_balance_relative')} |"
        )
    summary = report.get("conservation_summary", {})
    if summary:
        lines.extend([
            "",
            "## Conservation Gate",
            "",
            f"- electron_minus_hole sum: {summary.get('electron_minus_hole_sum_A_per_um')} A/um",
            f"- electron_plus_hole sum: {summary.get('electron_plus_hole_sum_A_per_um')} A/um",
            f"- electron-only sum: {summary.get('electron_sum_A_per_um')} A/um",
            f"- hole-only sum: {summary.get('hole_sum_A_per_um')} A/um",
            f"- max abs terminal current: {summary.get('max_abs_terminal_current_A_per_um')} A/um",
            f"- abs pair sum: {summary.get('abs_pair_sum_A_per_um')} A/um",
            f"- absolute floor gate pass: {summary.get('absolute_floor_gate_pass')}",
        ])
    parity = report.get("sentaurus_current_parity", {})
    if parity:
        lines.extend([
            "",
            "## Sentaurus Current Parity",
            "",
            f"- status: {parity.get('status')}",
            f"- classification: {parity.get('classification')}",
            f"- unit hypothesis: {parity.get('unit_hypothesis')}",
            f"- minority component zero contacts: {parity.get('minority_component_zero_contacts')}",
            "",
            "| Contact | Vela total A/um | Sentaurus total | Abs ratio | Sign |",
            "| --- | ---: | ---: | ---: | --- |",
        ])
        for contact, row in parity.get("by_contact", {}).items():
            lines.append(
                f"| {contact} | {row.get('vela_current_total_A_per_um')} | "
                f"{row.get('sentaurus_total_current')} | "
                f"{row.get('abs_ratio_sentaurus_to_vela')} | {row.get('sign_relation')} |"
            )
    lines.extend(["", "## Edge Coverage", ""])
    edges = report.get("contact_edges", {})
    mesh = report.get("mesh_reference", {})
    boundary_counts = mesh.get("sentaurus_contact_boundary_elements", {})
    for contact, count in edges.get("flux_link_count_by_contact", {}).items():
        boundary_count = boundary_counts.get(contact)
        lines.append(
            f"- {contact}: {count} Vela flux links, {boundary_count} Sentaurus boundary elements, "
            f"sum {edges.get('flux_link_sum_A_per_um_by_contact', {}).get(contact)} A/um"
        )
    if mesh:
        lines.extend([
            "",
            f"Sentaurus total elements including contacts: {mesh.get('sentaurus_total_elements_including_contacts')}",
            f"Sentaurus bulk triangle elements: {mesh.get('sentaurus_bulk_triangle_elements')}",
            "",
        ])
    lines.extend(["", "## Reasons", ""])
    for reason in report.get("classification_reasons", []):
        lines.append(f"- {reason}")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--existing-terminal-balance-csv", type=Path)
    parser.add_argument("--existing-contact-edge-csv", type=Path)
    parser.add_argument("--existing-vtk", type=Path)
    parser.add_argument("--pair-balance-relative-gate", "--current-pair-relative-gate",
                        dest="pair_balance_relative_gate", type=float, default=5.0e-2)
    parser.add_argument("--abs-current-gate-a-per-um", "--current-pair-abs-gate-A-per-um",
                        dest="abs_current_gate_a_per_um", type=float, default=1.0e-24)
    parser.add_argument("--require-balanced", action="store_true")
    parser.add_argument("--top-edges", type=int, default=10)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.output_dir / REPORT_NAME
    md_path = args.output_dir / MARKDOWN_NAME

    try:
        if args.existing_terminal_balance_csv and args.existing_contact_edge_csv and args.existing_vtk:
            terminal_csv = args.existing_terminal_balance_csv
            edge_csv = args.existing_contact_edge_csv
            vtk_path = args.existing_vtk
        else:
            if not args.runner:
                raise RuntimeError(
                    "--runner is required unless existing terminal, edge, and VTK files are supplied"
                )
            vtk_path, terminal_csv, edge_csv = run_probe(args.reference_root, args.output_dir, args.runner)

        terminal = terminal_report(
            terminal_csv,
            args.pair_balance_relative_gate,
            args.abs_current_gate_a_per_um,
        )
        edges = edge_report(edge_csv, terminal, args.top_edges)
        vtk = vtk_report(vtk_path)
        mesh_reference = mesh_reference_report(args.reference_root)
        classification, reasons, flags = classify(terminal, edges, vtk)
        sentaurus_current = sentaurus_current_reference(args.reference_root)
        sentaurus_parity = sentaurus_current_parity(terminal, sentaurus_current)
        report = {
            "status": "pass" if classification == "balanced" else "diagnostic_fail",
            "classification": classification,
            "classification_reasons": reasons,
            "root_cause_flags": flags,
            "reference_root": str(args.reference_root),
            "sentaurus_log_reference": sentaurus_log_reference(args.reference_root),
            "sentaurus_current_reference": sentaurus_current,
            "sentaurus_current_parity": sentaurus_parity,
            "mesh_reference": mesh_reference,
            "terminal_balance": terminal,
            "conservation_summary": conservation_summary(terminal),
            "contact_edges": edges,
            "vtk_fields": vtk,
            "thresholds": {
                "pair_balance_relative_gate": args.pair_balance_relative_gate,
                "abs_current_gate_a_per_um": args.abs_current_gate_a_per_um,
            },
        }
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1 if args.require_balanced and classification != "balanced" else 0
    except Exception as exc:
        report = {
            "status": "error",
            "classification": "input_error",
            "classification_reasons": [str(exc)],
            "root_cause_flags": {},
            "reference_root": str(args.reference_root),
        }
        write_json(report_path, report)
        md_path.write_text(markdown_report(report))
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
