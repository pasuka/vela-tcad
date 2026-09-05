#!/usr/bin/env python3
"""Audit and normalize the Genius NPN BJT Sentaurus WP1-WP2 oracle."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sentaurus_import import parse_quoted_list, parse_values_block  # noqa: E402


TERMINALS = ("collector", "base", "emitter")
TARGET_VCE = tuple(index / 10.0 for index in range(31))


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def outside_distance(value: float, lo: float, hi: float) -> float:
    return max(abs(value - 0.5 * (lo + hi)) - 0.5 * (hi - lo), 0.0)


def analytical_profile(
    x_um: float,
    y_um: float,
    peak_cm3: float,
    window_um: tuple[float, float, float, float],
    char_um: tuple[float, float],
) -> float:
    xmin, ymin, xmax, ymax = window_um
    xchar, ychar = char_um
    dx = outside_distance(x_um, xmin, xmax)
    dy = outside_distance(y_um, ymin, ymax)
    return peak_cm3 * math.exp(-((dx / xchar) ** 2) - ((dy / ychar) ** 2))


def expected_doping(x_um: float, y_um: float) -> tuple[float, float]:
    donors = 5.0e15
    donors += analytical_profile(
        x_um, y_um, 7.0e19, (2.75, 0.0, 4.25, 0.0), (0.1275, 0.17)
    )
    donors += analytical_profile(
        x_um, y_um, 1.0e19, (0.0, 2.0, 6.0, 2.0), (0.27, 0.27)
    )
    acceptors = analytical_profile(
        x_um, y_um, 6.0e17, (1.25, 0.0, 4.75, 0.35), (0.12, 0.16)
    )
    acceptors += analytical_profile(
        x_um, y_um, 4.0e18, (1.25, 0.0, 4.75, 0.0), (0.12, 0.16)
    )
    return donors, acceptors


def _relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(expected), 1.0)


def audit_structure(neutral_dir: Path) -> dict[str, object]:
    nodes = _read_csv(neutral_dir / "nodes.csv")
    doping = _read_csv(neutral_dir / "doping.csv")
    contacts = _read_csv(neutral_dir / "contacts.csv")
    node_by_id = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"])) for row in nodes
    }
    doping_by_id = {int(row["node_id"]): row for row in doping}

    donor_errors: list[float] = []
    acceptor_errors: list[float] = []
    worst: dict[str, object] = {}
    for node_id, (x_um, y_um) in node_by_id.items():
        row = doping_by_id[node_id]
        actual_donor = float(row["donors_cm3"])
        actual_acceptor = float(row["acceptors_cm3"])
        expected_donor, expected_acceptor = expected_doping(x_um, y_um)
        donor_error = _relative_error(actual_donor, expected_donor)
        acceptor_error = _relative_error(actual_acceptor, expected_acceptor)
        donor_errors.append(donor_error)
        acceptor_errors.append(acceptor_error)
        score = max(donor_error, acceptor_error)
        if not worst or score > float(worst["relative_error"]):
            worst = {
                "node_id": node_id,
                "x_um": x_um,
                "y_um": y_um,
                "relative_error": score,
                "actual_donor_cm3": actual_donor,
                "expected_donor_cm3": expected_donor,
                "actual_acceptor_cm3": actual_acceptor,
                "expected_acceptor_cm3": expected_acceptor,
            }

    expected_contacts = {
        "base": (0.0, 1.25, 2.0),
        "emitter": (0.0, 2.75, 4.25),
        "collector": (2.0, 0.0, 6.0),
    }
    contact_audit: dict[str, object] = {}
    contact_pass = True
    by_name = {row["name"].lower(): row for row in contacts}
    for name, (expected_y, expected_xmin, expected_xmax) in expected_contacts.items():
        row = by_name.get(name)
        if row is None:
            contact_audit[name] = {"present": False}
            contact_pass = False
            continue
        ids = [int(value) for value in row["node_ids"].split(";") if value]
        points = [node_by_id[node_id] for node_id in ids]
        xs = [point[0] for point in points]
        ys = [point[1] for point in points]
        max_coordinate_error = max(
            abs(min(xs) - expected_xmin),
            abs(max(xs) - expected_xmax),
            max(abs(value - expected_y) for value in ys),
        )
        passed = len(ids) >= 2 and max_coordinate_error <= 1.0e-8
        contact_pass = contact_pass and passed
        contact_audit[name] = {
            "present": True,
            "node_count": len(ids),
            "xmin_um": min(xs),
            "xmax_um": max(xs),
            "ymin_um": min(ys),
            "ymax_um": max(ys),
            "maximum_coordinate_error_um": max_coordinate_error,
            "pass": passed,
        }

    max_donor_error = max(donor_errors, default=math.inf)
    max_acceptor_error = max(acceptor_errors, default=math.inf)
    doping_pass = max(max_donor_error, max_acceptor_error) <= 1.0e-6
    return {
        "mesh": {
            "node_count": len(nodes),
            "bounds_um": {
                "xmin": min(point[0] for point in node_by_id.values()),
                "ymin": min(point[1] for point in node_by_id.values()),
                "xmax": max(point[0] for point in node_by_id.values()),
                "ymax": max(point[1] for point in node_by_id.values()),
            },
            "minimum_node_count": 1000,
            "pass": len(nodes) >= 1000,
        },
        "contacts": contact_audit,
        "contact_pass": contact_pass,
        "doping": {
            "maximum_donor_relative_error": max_donor_error,
            "maximum_acceptor_relative_error": max_acceptor_error,
            "worst_node": worst,
            "tolerance": 1.0e-6,
            "pass": doping_pass,
        },
        "pass": len(nodes) >= 1000 and contact_pass and doping_pass,
    }


def _find_dataset(datasets: Iterable[str], terminal: str, quantity: str) -> str:
    names = list(datasets)
    terminal_matches = [name for name in names if terminal in name.lower()]
    if quantity == "voltage":
        preferred = [name for name in terminal_matches if "outervoltage" in name.lower()]
        fallback = [name for name in terminal_matches if "voltage" in name.lower()]
    else:
        preferred = [name for name in terminal_matches if "totalcurrent" in name.lower()]
        fallback = [name for name in terminal_matches if "current" in name.lower()]
    matches = preferred or fallback
    if not matches:
        raise ValueError(f"cannot find {terminal} {quantity} in PLT datasets: {names}")
    return matches[0]


def normalize_plt(plt_path: Path, model: str) -> tuple[list[dict[str, float | str]], dict[str, object]]:
    text = plt_path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if not datasets:
        raise ValueError(f"{plt_path} does not declare datasets")
    rows = parse_values_block(text, len(datasets))
    indices: dict[str, int] = {}
    selected_datasets: dict[str, str] = {}
    for terminal in TERMINALS:
        for quantity in ("voltage", "current"):
            key = f"{terminal}_{quantity}"
            name = _find_dataset(datasets, terminal, quantity)
            selected_datasets[key] = name
            indices[key] = datasets.index(name)

    normalized: list[dict[str, float | str]] = []
    maximum_voltage_error = 0.0
    for target_vce in TARGET_VCE:
        row = min(rows, key=lambda values: abs(values[indices["collector_voltage"]] - target_vce))
        vce = row[indices["collector_voltage"]]
        voltage_error = abs(vce - target_vce)
        maximum_voltage_error = max(maximum_voltage_error, voltage_error)
        currents = {terminal: row[indices[f"{terminal}_current"]] for terminal in TERMINALS}
        kcl = sum(currents.values())
        current_scale = max((abs(value) for value in currents.values()), default=0.0)
        relative_kcl = abs(kcl) / max(current_scale, 1.0e-30)
        ib = currents["base"]
        ie = currents["emitter"]
        normalized.append(
            {
                "model": model,
                "VBE_V": row[indices["base_voltage"]] - row[indices["emitter_voltage"]],
                "VCE_V": vce - row[indices["emitter_voltage"]],
                "Ic_A_per_um": currents["collector"],
                "Ib_A_per_um": currents["base"],
                "Ie_A_per_um": currents["emitter"],
                "beta_abs": abs(currents["collector"] / ib) if ib else math.inf,
                "alpha_abs": abs(currents["collector"] / ie) if ie else math.inf,
                "kcl_abs_A_per_um": abs(kcl),
                "kcl_relative": relative_kcl,
            }
        )

    max_kcl_abs = max(float(row["kcl_abs_A_per_um"]) for row in normalized)
    max_kcl_relative = max(float(row["kcl_relative"]) for row in normalized)
    kcl_pass = max_kcl_abs <= 1.0e-18 or max_kcl_relative <= 1.0e-6
    audit = {
        "input": str(plt_path),
        "raw_row_count": len(rows),
        "normalized_row_count": len(normalized),
        "selected_datasets": selected_datasets,
        "maximum_collector_voltage_selection_error_V": maximum_voltage_error,
        "maximum_absolute_kcl_residual_A_per_um": max_kcl_abs,
        "maximum_relative_kcl_residual": max_kcl_relative,
        "pass": len(normalized) == 31 and maximum_voltage_error <= 1.0e-6 and kcl_pass,
    }
    return normalized, audit


def write_curve(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--neutral-dir", type=Path, required=True)
    parser.add_argument("--m0-plt", type=Path, required=True)
    parser.add_argument("--m1-plt", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    structure = audit_structure(args.neutral_dir)
    m0_rows, m0_audit = normalize_plt(args.m0_plt, "M0")
    m1_rows, m1_audit = normalize_plt(args.m1_plt, "M1")
    write_curve(args.output_root / "reference_curves" / "bjt_m0_output.csv", m0_rows)
    write_curve(args.output_root / "reference_curves" / "bjt_m1_output.csv", m1_rows)

    reports = args.output_root / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    (reports / "structure_audit.json").write_text(
        json.dumps(structure, indent=2) + "\n", encoding="utf-8"
    )
    validation = {
        "schema_version": 1,
        "structure": structure,
        "M0": m0_audit,
        "M1": m1_audit,
        "pass": bool(structure["pass"] and m0_audit["pass"] and m1_audit["pass"]),
    }
    (reports / "wp0_wp2_validation.json").write_text(
        json.dumps(validation, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(validation, indent=2))
    return 0 if validation["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

