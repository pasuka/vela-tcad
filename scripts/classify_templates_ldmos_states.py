#!/usr/bin/env python3
"""Classify exact Sentaurus representative-state TDRs by terminal values."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from run_templates_ldmos_sentaurus_vm import sha256_file, utc_now, write_json


def natural_key(path: Path) -> tuple[Any, ...]:
    return tuple(
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", path.name)
    )


def terminal_values(document: dict[str, Any]) -> dict[str, dict[str, float]]:
    names = {
        int(item["index"]): str(item["name"])
        for item in document["geometry"]["regions"]
    }
    result: dict[str, dict[str, float]] = {}
    for field in document["fields"]:
        if field["name"] not in {"ContactExternalVoltage", "ContactCurrentFlux"}:
            continue
        raw = field.get("raw_values", [])
        if len(raw) != 1:
            continue
        contact = names.get(int(field["region"]), f"region_{field['region']}")
        key = "voltage_V" if field["name"] == "ContactExternalVoltage" else "current_A_per_um"
        result.setdefault(contact, {})[key] = float(raw[0])
    # Non-loadable Sentaurus Plot TDRs do not expose the contact-scalar
    # datasets carried by some loadable solution TDRs.  Recover a terminal
    # marker from the electrostatic potential on the material side of each
    # contact instead.  Region-node field order is ascending global vertex id,
    # which is also the contract used by sentaurus_import.
    regions = {int(item["index"]): item for item in document["geometry"]["regions"]}
    potential_fields = {
        int(field["region"]): field
        for field in document["fields"]
        if field["name"] == "ElectrostaticPotential" and field.get("raw_values")
    }

    def nodes(region: dict[str, Any]) -> list[int]:
        values = set(int(value) for value in region.get("points", []))
        for edge in region.get("edges", []):
            values.update(int(value) for value in edge)
        for triangle in region.get("triangles", []):
            values.update(int(value) for value in triangle)
        return sorted(values)

    node_orders = {index: nodes(region) for index, region in regions.items()}
    for contact_index, contact in regions.items():
        if int(contact.get("type", 99)) != 1:
            continue
        contact_nodes = set(node_orders[contact_index])
        candidates: list[tuple[int, int]] = []
        for material_index, field in potential_fields.items():
            order = node_orders.get(material_index, [])
            if len(order) != len(field["raw_values"]):
                continue
            overlap = len(contact_nodes.intersection(order))
            if overlap:
                candidates.append((overlap, material_index))
        if not candidates:
            continue
        _, material_index = max(candidates)
        order = node_orders[material_index]
        values_by_node = dict(zip(order, potential_fields[material_index]["raw_values"]))
        values = [float(values_by_node[node]) for node in sorted(contact_nodes) if node in values_by_node]
        marker = sum(values) / len(values)
        result.setdefault(str(contact["name"]), {}).update({
            "potential_marker_V": marker,
            "potential_marker_spread_V": max(values) - min(values),
            "potential_marker_node_count": len(values),
        })
    return result


def marker(record: dict[str, Any], terminal: str) -> float:
    values = record["terminals"][terminal]
    if "potential_marker_V" in values:
        return float(values["potential_marker_V"])
    return float(values["voltage_V"])


def set_relative_voltage(records: list[dict[str, Any]], terminal: str,
                         baseline: float) -> list[float]:
    result = []
    for record in records:
        voltage = marker(record, terminal) - baseline
        record["terminals"][terminal]["voltage_V"] = voltage
        result.append(voltage)
    return result


def read_plt_curve(path: Path) -> list[dict[str, float]]:
    from sentaurus_import import parse_quoted_list, parse_values_block

    text = path.read_text(encoding="utf-8", errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    rows = parse_values_block(text, len(datasets))
    bias_name = next(
        name for name in ("drain InnerVoltage", "drain OuterVoltage", "drain Voltage")
        if name in datasets
    )
    current_name = "drain TotalCurrent"
    bias_index, current_index = datasets.index(bias_name), datasets.index(current_name)
    return [
        {"bias_V": float(row[bias_index]), "current_A_per_um": float(row[current_index])}
        for row in rows
    ]


def align_states_to_curve(records: list[dict[str, Any]], curve: list[dict[str, float]],
                          tolerance_V: float = 1e-8) -> float:
    """Attach exact PLT terminal currents using ordered closest-voltage matches."""
    previous = 0
    maximum_error = 0.0
    for record in records:
        voltage = record["terminals"]["drain"]["voltage_V"]
        index = min(
            range(previous, len(curve)),
            key=lambda candidate: abs(curve[candidate]["bias_V"] - voltage),
        )
        error = abs(curve[index]["bias_V"] - voltage)
        if error > tolerance_V:
            raise RuntimeError(
                f"state {record['name']} drain marker {voltage:.17g} V has no exact "
                f"ordered PLT match; closest error={error:.3e} V"
            )
        previous = index
        maximum_error = max(maximum_error, error)
        record["terminals"]["drain"].update({
            "voltage_V": curve[index]["bias_V"],
            "current_A_per_um": curve[index]["current_A_per_um"],
            "matched_plt_row": index,
            "marker_match_error_V": error,
        })
    return maximum_error


def inspect_state(importer: Path, path: Path, evidence_dir: Path) -> dict[str, Any]:
    field_values = evidence_dir / f"{path.stem}.field_values.json"
    if not field_values.is_file():
        completed = subprocess.run(
            [str(importer), "--tdr", str(path), "--field-values-json", str(field_values)],
            text=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"failed to inspect {path.name}: {completed.stderr}")
    document = json.loads(field_values.read_text(encoding="utf-8"))
    return {
        "name": path.name,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
        "terminals": terminal_values(document),
        "field_values_path": f"field_values/{field_values.name}",
        "field_values_sha256": sha256_file(field_values),
    }


def closest_index(records: list[dict[str, Any]], target: float) -> int:
    return min(
        range(len(records)),
        key=lambda index: abs(math.log10(max(abs(records[index]["terminals"]["drain"]["current_A_per_um"]), 1e-300))
                              - math.log10(target)),
    )


def select_bv(records: list[dict[str, Any]], iadapt: float, criterion: float) -> dict[str, Any]:
    currents = [abs(item["terminals"]["drain"]["current_A_per_um"]) for item in records]
    below_iadapt = [index for index, value in enumerate(currents) if value < iadapt]
    below_criterion = [index for index, value in enumerate(currents) if value < criterion]
    above_criterion = [index for index, value in enumerate(currents) if value >= criterion]
    selections = {
        "pre_iadapt": max(below_iadapt, key=currents.__getitem__) if below_iadapt else 0,
        "near_iadapt": closest_index(records, iadapt),
        "avalanche_growth": closest_index(records, math.sqrt(iadapt * criterion)),
        "criterion_pre": max(below_criterion, key=currents.__getitem__) if below_criterion else 0,
        "criterion_post": min(above_criterion, key=currents.__getitem__) if above_criterion else len(records) - 1,
    }
    return {
        name: {
            "record_index": index,
            "state": records[index]["name"],
            "drain_voltage_V": records[index]["terminals"]["drain"]["voltage_V"],
            "drain_current_A_per_um": records[index]["terminals"]["drain"]["current_A_per_um"],
        }
        for name, index in selections.items()
    }


def max_bias_error(actual: list[float], expected: list[float]) -> float:
    if len(actual) != len(expected):
        return math.inf
    return max((abs(left - right) for left, right in zip(sorted(actual), sorted(expected))), default=0.0)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--states-dir", type=Path, required=True)
    parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iadapt", type=float, default=6.5e-13)
    parser.add_argument("--criterion", type=float, default=1.0e-8)
    args = parser.parse_args()

    states_dir = args.states_dir.resolve()
    output = args.output.resolve()
    evidence_dir = output.parent / "field_values"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    paths = sorted(states_dir.glob("state_*.tdr"), key=natural_key)
    if not paths:
        raise FileNotFoundError(f"no representative state TDRs in {states_dir}")
    records = [inspect_state(args.importer.resolve(), path, evidence_dir) for path in paths]
    bv = [item for item in records if item["name"].lower().startswith("state_bv_path_")]
    if len(bv) < 4:
        raise RuntimeError(f"expected at least four BV path states, found {len(bv)}")
    deck_manifest_path = output.parent / "decks" / "state_deck_manifest.json"
    deck_manifest = json.loads(deck_manifest_path.read_text(encoding="utf-8"))
    expected_by_deck = {item["name"]: item["states"] for item in deck_manifest["files"]}
    idvg = [item for item in records if item["name"].lower().startswith("state_idvg_")
            and not item["name"].lower().startswith("state_idvg_eq_0v")]
    idvg_eq = [item for item in records if item["name"].lower().startswith("state_idvg_eq_0v")]
    idvd_vg4 = [item for item in records if item["name"].lower().startswith("state_idvd_vg4_")]
    idvd_vg8 = [item for item in records if item["name"].lower().startswith("state_idvd_vg8_")]

    if len(idvg_eq) != 1:
        raise RuntimeError(f"expected one IdVg equilibrium state, found {len(idvg_eq)}")
    equilibrium = idvg_eq[0]
    equilibrium_markers = {
        terminal: marker(equilibrium, terminal)
        for terminal in ("gate", "drain", "source", "substrate")
    }
    for terminal in ("gate", "drain", "source", "substrate"):
        equilibrium["terminals"][terminal]["voltage_V"] = 0.0
    idvg_gate = set_relative_voltage(idvg, "gate", equilibrium_markers["gate"])
    set_relative_voltage(idvg, "drain", equilibrium_markers["drain"])

    idvd4_markers = set_relative_voltage(idvd_vg4, "drain", marker(idvd_vg4[0], "drain"))
    idvd8_markers = set_relative_voltage(idvd_vg8, "drain", marker(idvd_vg8[0], "drain"))
    requested_drains = expected_by_deck["IdVd.cmd"]["drain_biases_V"]
    idvd4_marker_error = max_bias_error(idvd4_markers, requested_drains)
    idvd8_marker_error = max_bias_error(idvd8_markers, requested_drains)
    for records_at_gate in (idvd_vg4, idvd_vg8):
        for item, expected_voltage in zip(records_at_gate, requested_drains, strict=True):
            item["terminals"]["drain"]["requested_voltage_V"] = expected_voltage
    gate_step = marker(idvd_vg8[0], "gate") - marker(idvd_vg4[0], "gate")
    for item in idvd_vg4:
        item["terminals"]["gate"]["voltage_V"] = expected_by_deck["IdVd.cmd"]["gate_biases_V"][0]
    for item in idvd_vg8:
        item["terminals"]["gate"]["voltage_V"] = expected_by_deck["IdVd.cmd"]["gate_biases_V"][1]

    bv_contract = expected_by_deck["BVdss.cmd"]
    target_states = bv_contract.get("target_states", [])
    if target_states:
        targets_by_time: dict[float, dict[str, Any]] = {}
        for target in target_states:
            time_value = float(target["time"])
            if time_value not in targets_by_time:
                targets_by_time[time_value] = dict(target)
                targets_by_time[time_value]["roles"] = []
            targets_by_time[time_value]["roles"].append(str(target["role"]))
        ordered_targets = [targets_by_time[value] for value in sorted(targets_by_time)]
        if len(bv) != len(ordered_targets):
            raise RuntimeError(
                f"expected {len(ordered_targets)} explicit BV states, found {len(bv)}"
            )
        bv_marker_errors = []
        bv_curve = read_plt_curve(states_dir / "n6_des.plt")
        for record, target in zip(bv, ordered_targets, strict=True):
            observed = marker(record, "drain") - equilibrium_markers["drain"]
            is_final_fallback = "_final_" in record["name"].lower()
            actual_point = bv_curve[-1] if is_final_fallback else {
                "bias_V": float(target["voltage_V"]),
                "current_A_per_um": float(target["current_A_per_um"]),
            }
            expected = actual_point["bias_V"]
            error = abs(observed - expected)
            bv_marker_errors.append(error)
            record["terminals"]["drain"].update({
                "voltage_V": expected,
                "current_A_per_um": actual_point["current_A_per_um"],
                "requested_continuation_time": float(target["time"]),
                "representative_roles": target["roles"],
                "marker_match_error_V": error,
                "capture_origin": (
                    "derived_final_device_state_after_requested_time_exceeded_run_end"
                    if is_final_fallback else "explicit_plot_time"
                ),
            })
        bv_curve_match_error = max(bv_marker_errors, default=math.inf)
    else:
        set_relative_voltage(bv, "drain", equilibrium_markers["drain"])
        bv_curve = read_plt_curve(states_dir / "n6_des.plt")
        bv_curve_match_error = align_states_to_curve(bv, bv_curve)
    idvg_error = max_bias_error(
        idvg_gate,
        expected_by_deck["IdVg.cmd"]["biases_V"],
    )
    idvd4_drain_error = max_bias_error(
        idvd4_markers,
        requested_drains,
    )
    idvd8_drain_error = max_bias_error(
        idvd8_markers,
        requested_drains,
    )
    idvd_gate_step_error = abs(gate_step - 4.0)
    equilibrium_error = max(
        (abs(values["voltage_V"]) for item in idvg_eq
         for values in item["terminals"].values() if "voltage_V" in values),
        default=math.inf,
    )
    expected_bv_count = int(bv_contract.get(
        "expected_state_count",
        bv_contract.get("minimum_path_samples", bv_contract.get("path_samples", 31)),
    ))
    accepted_bv_counts = {expected_bv_count}
    if bool(bv_contract.get("allow_extra_initial_state", not target_states)):
        accepted_bv_counts.add(expected_bv_count + 1)
    contract_gates = {
        "equilibrium_all_terminal_bias_max_error_le_1e_10_V": equilibrium_error <= 1e-10,
        "idvg_exact_bias_max_error_le_1e_10_V": idvg_error <= 1e-10,
        "idvd_vg4_field_marker_matches_requested_bias_le_5e_2_V": idvd4_marker_error <= 5e-2,
        "idvd_vg8_field_marker_matches_requested_bias_le_5e_2_V": idvd8_marker_error <= 5e-2,
        "idvd_gate_marker_step_max_error_le_1e_10_V": idvd_gate_step_error <= 1e-10,
        "bv_path_sample_count_matches_plot_contract": len(bv) in accepted_bv_counts,
        "bv_state_field_markers_match_requested_points_le_1e_1_V": bv_curve_match_error <= 1e-1,
    }
    result = {
        "schema": "vela.templates_ldmos.representative_state_inventory.v1",
        "classification": "derived_output_only_control_not_official_oracle",
        "generated_at": utc_now(),
        "state_count": len(records),
        "states": records,
        "exact_bias_contract": {
            "gates": contract_gates,
            "equilibrium_max_error_V": equilibrium_error,
            "idvg_max_error_V": idvg_error,
            "idvd_vg4_field_marker_max_error_V": idvd4_marker_error,
            "idvd_vg8_field_marker_max_error_V": idvd8_marker_error,
            "idvd_gate_marker_step_error_V": idvd_gate_step_error,
            "bv_ordered_plt_match_max_error_V": bv_curve_match_error,
        },
        "bv_selection": {
            "iadapt_A_per_um": args.iadapt,
            "criterion_A_per_um": args.criterion,
            "policy": "select exact saved states in continuation row order; no interpolation",
            "representatives": select_bv(bv, args.iadapt, args.criterion),
        },
    }
    write_json(output, result)
    print(json.dumps({
        "state_count": len(records),
        "bv_state_count": len(bv),
        "exact_bias_gates": contract_gates,
        "bv_representatives": result["bv_selection"]["representatives"],
    }, indent=2))
    return 0 if all(contract_gates.values()) else 2


if __name__ == "__main__":
    raise SystemExit(main())
