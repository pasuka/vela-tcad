#!/usr/bin/env python3
"""Run the SimpleMOS M20 drain-adjacent continuity-row audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
M14 = ROOT / "m14_adjoint_attribution"
M19 = ROOT / "m19_drain_cut_audit"
M19_PORTABLE = REPO / "reference_tcad/simplemos_sentaurus2022/drain_cut_audit"
OUTPUT = ROOT / "m20_drain_adjacent_row_audit"
PORTABLE = REPO / "reference_tcad/simplemos_sentaurus2022/drain_adjacent_row_audit"
RUNNER = REPO / "build-release/vela_example_runner.exe"
KEY_STATE = "n21_vd_0p05_vg_0p8"
PROBE_TYPES = {
    "terms": "newton_carrier_term_probe",
    "rows": "newton_carrier_row_probe",
    "jacobian": "transport_edge_jacobian_probe",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def run_config(config: Path, runner: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run([str(runner), "--config", str(config)],
                               cwd=REPO, text=True, capture_output=True,
                               env=env, check=False)
    config.with_suffix(".stdout.txt").write_text(completed.stdout, encoding="utf-8")
    config.with_suffix(".stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"M20 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def f(row: dict[str, str], column: str) -> float:
    return float(row[column])


def signed_edge_flux(row: dict[str, str], node: int) -> float:
    if int(row["node0"]) == node:
        return f(row, "electron_flux")
    if int(row["node1"]) == node:
        return -f(row, "electron_flux")
    raise ValueError(f"node {node} is not on edge {row['edge_id']}")


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(actual), abs(expected), 1.0e-300)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable", type=Path, default=PORTABLE)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    runner = args.runner.resolve()
    output = args.output.resolve()
    portable = args.portable.resolve()
    output.mkdir(parents=True, exist_ok=True)
    portable.mkdir(parents=True, exist_ok=True)

    state_source = read_csv(M19_PORTABLE / "m19_state_summary.csv")
    cases: list[dict[str, Any]] = []
    jobs: list[tuple[Path, Path, str]] = []
    for state_row in state_source:
        state = state_row["state"]
        case_dir = output / state
        case_dir.mkdir(parents=True, exist_ok=True)
        base_config = read_json(M14 / state / "baseline_functional.json")
        mesh_path = Path(base_config["mesh_file"])
        mesh = read_json(mesh_path)
        drain_contacts = [item for item in mesh["contacts"]
                          if item["name"].lower() == "drain"]
        if len(drain_contacts) != 1:
            raise ValueError(f"{state}: expected one drain contact")
        contact_nodes = {int(node) for node in drain_contacts[0]["node_ids"]}
        baseline_sg = M19 / state / "cut_mask_0_sg_edges.csv"
        sentaurus_sg = M19 / state / "full_sentaurus_sg_edges.csv"
        sg_rows = read_csv(baseline_sg)
        cut_edges = {
            int(row["edge_id"]) for row in sg_rows
            if (int(row["node0"]) in contact_nodes)
            != (int(row["node1"]) in contact_nodes)
        }
        adjacent_nodes = sorted({
            node for row in sg_rows if int(row["edge_id"]) in cut_edges
            for node in (int(row["node0"]), int(row["node1"]))
            if node not in contact_nodes
        })
        variants: dict[str, dict[str, Path]] = {}
        for variant, state_file, sg_file in (
                ("baseline", M14 / state / "baseline_state.csv", baseline_sg),
                ("sentaurus", M14 / state / "field_all_state.csv", sentaurus_sg)):
            paths: dict[str, Path] = {"state": state_file, "sg": sg_file}
            for label, simulation_type in PROBE_TYPES.items():
                config_path = case_dir / f"{variant}_{label}.json"
                csv_path = case_dir / f"{variant}_{label}.csv"
                deck = dict(base_config)
                deck.update({
                    "simulation_type": simulation_type,
                    "state_file": str(state_file.resolve()),
                    "output_csv": str(csv_path.resolve()),
                    "simplemos_m20": {
                        "read_only": True,
                        "variant": variant,
                        "target_nodes": adjacent_nodes,
                    },
                })
                deck.pop("contact", None)
                deck.pop("residual_output_csv", None)
                if label == "terms":
                    deck["carrier_term_probe"] = {"solved_equation_terms": False}
                if label == "jacobian":
                    deck["physical_finite_difference_step_V"] = 1.0e-7
                write_json(config_path, deck)
                paths[label] = csv_path
                if (args.force or not csv_path.is_file()
                        or not config_path.with_suffix(".stdout.txt").is_file()):
                    jobs.append((config_path, csv_path, f"{state}:{variant}:{label}"))
            variants[variant] = paths
        cases.append({
            "state": state,
            "state_row": state_row,
            "mesh": mesh,
            "contact_nodes": contact_nodes,
            "cut_edges": cut_edges,
            "adjacent_nodes": adjacent_nodes,
            "variants": variants,
        })

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_config, config, runner): label
                   for config, _expected, label in jobs}
        for future in as_completed(futures):
            label = futures[future]
            future.result()
            print(f"completed {label}", flush=True)

    node_rows: list[dict[str, Any]] = []
    jacobian_rows: list[dict[str, Any]] = []
    delta_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    max_flux_closure = 0.0
    max_term_closure = 0.0
    max_jacobian_partition_closure = 0.0
    for case in cases:
        state = case["state"]
        variant_nodes: dict[str, dict[int, dict[str, Any]]] = {}
        for variant, paths in case["variants"].items():
            terms = {int(row["node_id"]): row for row in read_csv(paths["terms"])}
            rows = {int(row["node_id"]): row for row in read_csv(paths["rows"])}
            sg = read_csv(paths["sg"])
            jac = [row for row in read_csv(paths["jacobian"])
                   if row["carrier"] == "electron"]
            by_node: dict[int, dict[str, Any]] = {}
            for node in case["adjacent_nodes"]:
                incident = [row for row in sg
                            if node in (int(row["node0"]), int(row["node1"]))]
                contact_edges = [row for row in incident
                                 if int(row["edge_id"]) in case["cut_edges"]]
                internal_edges = [row for row in incident
                                  if int(row["edge_id"]) not in case["cut_edges"]]
                contact_flux = sum(signed_edge_flux(row, node)
                                   for row in contact_edges)
                internal_flux = sum(signed_edge_flux(row, node)
                                    for row in internal_edges)
                contact_flux_abs = sum(abs(signed_edge_flux(row, node))
                                       for row in contact_edges)
                internal_flux_abs = sum(abs(signed_edge_flux(row, node))
                                        for row in internal_edges)
                term = terms[node]
                row = rows[node]
                source = sum(f(term, column) for column in (
                    "electron_recombination", "electron_impact",
                    "electron_gauge", "electron_boundary"))
                flux_closure = contact_flux + internal_flux - f(term, "electron_flux")
                term_closure = contact_flux + internal_flux + source - f(
                    term, "electron_residual")
                max_flux_closure = max(max_flux_closure, abs(flux_closure))
                max_term_closure = max(max_term_closure, abs(term_closure))

                node_jac = [item for item in jac if int(item["row_node"]) == node]
                if not node_jac:
                    raise AssertionError(f"{state}:{variant}: no Jacobian edges at {node}")
                weights = {f(item, "continuity_row_weight") for item in node_jac}
                if len(weights) != 1:
                    raise AssertionError(f"{state}:{variant}: row-weight drift at {node}")
                weight = next(iter(weights))
                transport_by_column: dict[int, float] = defaultdict(float)
                for item in node_jac:
                    transport_by_column[int(item["column_node"])] += f(
                        item, "solver_production_edge_derivative")
                exact_partition_sum = sum(f(row, column) for column in (
                    "electron_psi_contact_column_abs_sum",
                    "electron_psi_free_column_abs_sum",
                    "electron_phin_contact_column_abs_sum",
                    "electron_phin_free_column_abs_sum",
                    "electron_phip_contact_column_abs_sum",
                    "electron_phip_free_column_abs_sum"))
                jacobian_partition_closure = relative_error(
                    exact_partition_sum, f(row, "electron_row_abs_sum"))
                max_jacobian_partition_closure = max(
                    max_jacobian_partition_closure, jacobian_partition_closure)
                transport_contact_abs = sum(
                    abs(value) for column, value in transport_by_column.items()
                    if column in case["contact_nodes"])
                transport_free_abs = sum(
                    abs(value) for column, value in transport_by_column.items()
                    if column not in case["contact_nodes"])
                node_record = {
                    "state": state,
                    "device": case["state_row"]["device"],
                    "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                    "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                    "variant": variant,
                    "node_id": node,
                    "x_um": term["x"],
                    "y_um": term["y"],
                    "contact_edge_count": len(contact_edges),
                    "internal_edge_count": len(internal_edges),
                    "raw_contact_edge_flux": contact_flux,
                    "raw_internal_edge_flux": internal_flux,
                    "raw_total_edge_flux": contact_flux + internal_flux,
                    "raw_contact_edge_flux_abs_sum": contact_flux_abs,
                    "raw_internal_edge_flux_abs_sum": internal_flux_abs,
                    "raw_flux_condition_number": (
                        (contact_flux_abs + internal_flux_abs)
                        / max(abs(contact_flux + internal_flux), 1.0e-300)),
                    "raw_srh_recombination": f(term, "electron_recombination"),
                    "raw_impact": f(term, "electron_impact"),
                    "raw_gauge": f(term, "electron_gauge"),
                    "raw_boundary": f(term, "electron_boundary"),
                    "raw_source_sum": source,
                    "raw_electron_residual": f(term, "electron_residual"),
                    "raw_edge_flux_closure": flux_closure,
                    "raw_term_sum_closure": term_closure,
                    "continuity_row_weight": weight,
                    "production_electron_residual": f(row, "electron_residual"),
                    "production_electron_diagonal": f(row, "electron_diagonal"),
                    "production_electron_row_abs_sum": f(row, "electron_row_abs_sum"),
                    "production_psi_column_abs_sum": f(
                        row, "electron_psi_column_abs_sum"),
                    "production_phin_column_abs_sum": f(
                        row, "electron_phin_column_abs_sum"),
                    "production_phip_column_abs_sum": f(
                        row, "electron_phip_column_abs_sum"),
                    "production_contact_column_abs_sum": f(
                        row, "electron_contact_column_abs_sum"),
                    "production_free_column_abs_sum": f(
                        row, "electron_free_column_abs_sum"),
                    "production_psi_contact_column_abs_sum": f(
                        row, "electron_psi_contact_column_abs_sum"),
                    "production_psi_free_column_abs_sum": f(
                        row, "electron_psi_free_column_abs_sum"),
                    "production_phin_contact_column_abs_sum": f(
                        row, "electron_phin_contact_column_abs_sum"),
                    "production_phin_free_column_abs_sum": f(
                        row, "electron_phin_free_column_abs_sum"),
                    "production_phip_contact_column_abs_sum": f(
                        row, "electron_phip_contact_column_abs_sum"),
                    "production_phip_free_column_abs_sum": f(
                        row, "electron_phip_free_column_abs_sum"),
                    "exact_jacobian_partition_sum": exact_partition_sum,
                    "exact_jacobian_partition_relative_error":
                        jacobian_partition_closure,
                    "edge_projection_transport_phin_column_abs_sum": sum(
                        abs(value) for value in transport_by_column.values()),
                    "edge_projection_transport_contact_phin_column_abs_sum":
                        transport_contact_abs,
                    "edge_projection_transport_free_phin_column_abs_sum":
                        transport_free_abs,
                }
                node_rows.append(node_record)
                by_node[node] = node_record
                for item in node_jac:
                    edge_id = int(item["edge_id"])
                    edge_class = "contact_cut" if edge_id in case["cut_edges"] else "internal"
                    jacobian_rows.append({
                        "state": state,
                        "variant": variant,
                        "target_row_node": node,
                        "edge_class": edge_class,
                        **item,
                    })
            variant_nodes[variant] = by_node

        component_abs = defaultdict(float)
        component_signed = defaultdict(float)
        primary_abs = defaultdict(float)
        for node in case["adjacent_nodes"]:
            base = variant_nodes["baseline"][node]
            sent = variant_nodes["sentaurus"][node]
            delta = {
                "contact_flux": sent["raw_contact_edge_flux"] - base["raw_contact_edge_flux"],
                "internal_flux": sent["raw_internal_edge_flux"] - base["raw_internal_edge_flux"],
                "srh": sent["raw_srh_recombination"] - base["raw_srh_recombination"],
                "other_source": ((sent["raw_impact"] + sent["raw_gauge"] + sent["raw_boundary"])
                                 - (base["raw_impact"] + base["raw_gauge"] + base["raw_boundary"])),
            }
            residual_delta = (sent["raw_electron_residual"]
                              - base["raw_electron_residual"])
            for component, value in delta.items():
                component_abs[component] += abs(value)
                component_signed[component] += value
                if node in {990, 991, 992}:
                    primary_abs[component] += abs(value)
            delta_rows.append({
                "state": state,
                "device": case["state_row"]["device"],
                "drain_voltage_V": case["state_row"]["drain_voltage_V"],
                "gate_voltage_V": case["state_row"]["gate_voltage_V"],
                "node_id": node,
                "delta_contact_flux": delta["contact_flux"],
                "delta_internal_flux": delta["internal_flux"],
                "delta_srh": delta["srh"],
                "delta_other_source": delta["other_source"],
                "delta_component_sum": sum(delta.values()),
                "delta_raw_residual": residual_delta,
                "delta_closure": sum(delta.values()) - residual_delta,
                "baseline_raw_residual": base["raw_electron_residual"],
                "sentaurus_state_raw_residual": sent["raw_electron_residual"],
                "baseline_production_residual": base["production_electron_residual"],
                "sentaurus_state_production_residual": sent["production_electron_residual"],
                "baseline_flux_condition_number": base["raw_flux_condition_number"],
                "sentaurus_flux_condition_number": sent["raw_flux_condition_number"],
                "sentaurus_psi_jacobian_share": (
                    sent["production_psi_column_abs_sum"]
                    / max(sent["production_electron_row_abs_sum"], 1.0e-300)),
                "sentaurus_phin_jacobian_share": (
                    sent["production_phin_column_abs_sum"]
                    / max(sent["production_electron_row_abs_sum"], 1.0e-300)),
                "sentaurus_phip_jacobian_share": (
                    sent["production_phip_column_abs_sum"]
                    / max(sent["production_electron_row_abs_sum"], 1.0e-300)),
                "sentaurus_contact_column_share": (
                    sent["production_contact_column_abs_sum"]
                    / max(sent["production_electron_row_abs_sum"], 1.0e-300)),
                "sentaurus_phin_contact_column_share": (
                    sent["production_phin_contact_column_abs_sum"]
                    / max(sent["production_electron_row_abs_sum"], 1.0e-300)),
                "sentaurus_phin_free_column_share": (
                    sent["production_phin_free_column_abs_sum"]
                    / max(sent["production_electron_row_abs_sum"], 1.0e-300)),
            })
        dominant = max(component_abs, key=component_abs.get)
        total_abs = sum(component_abs.values())
        state_rows.append({
            "state": state,
            "device": case["state_row"]["device"],
            "drain_voltage_V": case["state_row"]["drain_voltage_V"],
            "gate_voltage_V": case["state_row"]["gate_voltage_V"],
            "adjacent_node_count": len(case["adjacent_nodes"]),
            "adjacent_nodes": ";".join(str(node) for node in case["adjacent_nodes"]),
            "dominant_absolute_residual_delta_component": dominant,
            **{f"absolute_{key}_delta_sum": value for key, value in component_abs.items()},
            **{f"signed_{key}_delta_sum": value for key, value in component_signed.items()},
            **{f"absolute_{key}_delta_share": value / max(total_abs, 1.0e-300)
               for key, value in component_abs.items()},
            **{f"primary_990_992_absolute_{key}_delta_sum": value
               for key, value in primary_abs.items()},
        })

    write_csv(portable / "m20_node_balance.csv", node_rows)
    write_csv(portable / "m20_node_delta.csv", delta_rows)
    write_csv(portable / "m20_incident_edge_jacobian.csv", jacobian_rows)
    write_csv(portable / "m20_state_summary.csv", state_rows)
    key_nodes = [row for row in delta_rows if row["state"] == KEY_STATE]
    write_csv(portable / "m20_key_state_node_ledger.csv", key_nodes)

    key_state = next(row for row in state_rows if row["state"] == KEY_STATE)
    key_deltas = [row for row in delta_rows if row["state"] == KEY_STATE]
    primary = [row for row in key_deltas if int(row["node_id"]) in {990, 991, 992}]
    report = {
        "schema": "vela.simplemos.sdevice.m20_drain_adjacent_row_audit_report.v1",
        "execution": {
            "state_count": len(cases),
            "variant_count": 2,
            "probe_types_per_variant": len(PROBE_TYPES),
            "probe_count": len(cases) * 2 * len(PROBE_TYPES),
            "new_sentaurus_execution": False,
            "diagnostic_cpp_changed": True,
            "default_model_changed": False,
        },
        "closure": {
            "maximum_raw_edge_flux_closure": max_flux_closure,
            "maximum_raw_term_sum_closure": max_term_closure,
            "maximum_exact_jacobian_partition_relative_error":
                max_jacobian_partition_closure,
        },
        "findings": {
            "key_state": key_state,
            "key_state_primary_nodes": primary,
            "dominant_component_counts": dict(sorted({
                component: sum(row["dominant_absolute_residual_delta_component"] == component
                               for row in state_rows)
                for component in ("contact_flux", "internal_flux", "srh", "other_source")
            }.items())),
            "all_drain_adjacent_nodes": {
                row["state"]: [int(node) for node in row["adjacent_nodes"].split(";")]
                for row in state_rows
            },
        },
        "interpretation_limits": read_json(
            REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m20_drain_adjacent_row_audit_contract_v1.json"
        )["interpretation_limits"],
    }
    write_json(portable / "m20_drain_adjacent_row_audit_report.json", report)
    print(json.dumps({
        "status": "complete",
        "state_count": len(cases),
        "probe_count": len(cases) * 2 * len(PROBE_TYPES),
        "key_state_dominant_component": key_state[
            "dominant_absolute_residual_delta_component"],
        "maximum_raw_edge_flux_closure": max_flux_closure,
        "maximum_exact_jacobian_partition_relative_error":
            max_jacobian_partition_closure,
    }))


if __name__ == "__main__":
    main()
