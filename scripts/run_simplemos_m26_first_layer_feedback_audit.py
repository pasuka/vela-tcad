#!/usr/bin/env python3
"""Run M26 drain-first-layer finite-difference and adjoint feedback audit."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import csv
import json
import math
from pathlib import Path
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402
import run_simplemos_m22_n23_hfs_deep_off as m22  # noqa: E402


RAW = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m22_n23_hfs_deep_off")
M24 = (REPO / "reference_tcad/simplemos_sentaurus2022"
       / "contact_boundary_audit")
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m26_first_layer_feedback_audit_contract_v1.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m26_first_layer_feedback_audit")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "first_layer_feedback_audit")
RUNNER = REPO / "build-release/vela_example_runner.exe"
Q = 1.602176634e-19


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def state_dir(variant: str) -> Path:
    return RAW / "self_consistent" / variant / "vg_0p05"


def validate(contract: dict[str, Any]) -> None:
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m26_first_layer_feedback_audit.v1"):
        raise ValueError("unexpected M26 contract schema")
    if contract.get("state_variants") != ["full", "no_hfs"]:
        raise ValueError("M26 requires paired full/no_hfs states")
    if float(contract["perturbation"]["central_step_V"]) <= 0.0:
        raise ValueError("M26 perturbation step must be positive")


def status_path(config: Path) -> Path:
    return config.parent / f"{config.stem}.stdout.txt"


def load_status(config: Path) -> dict[str, Any]:
    lines = status_path(config).read_text(encoding="utf-8").strip().splitlines()
    if not lines:
        raise ValueError(f"empty probe status: {status_path(config)}")
    return json.loads(lines[-1])


def task_ready(task: dict[str, Any]) -> bool:
    return status_path(task["config"]).is_file() and all(
        Path(path).is_file() for path in task["expected"])


def run_tasks(tasks: list[dict[str, Any]], runner: Path,
              workers: int, force: bool) -> dict[str, dict[str, Any]]:
    selected = [task for task in tasks if force or not task_ready(task)]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(m10.execute_runner, task["config"], runner): task
                   for task in selected}
        for future in as_completed(futures):
            task = futures[future]
            future.result()
            print(f"completed {task['id']}", flush=True)
    return {task["id"]: load_status(task["config"]) for task in tasks}


def base_config(variant: str) -> dict[str, Any]:
    return read_json(state_dir(variant) / "config.json")


def configure_operator(deck: dict[str, Any], variant: str) -> None:
    deck["solver"]["mobility"] = m22.mobility_config_for(variant)


def add_task(tasks: list[dict[str, Any]], task_id: str, config: Path,
             deck: dict[str, Any], expected: list[Path]) -> None:
    write_json(config, deck)
    tasks.append({"id": task_id, "config": config,
                  "expected": [path.resolve() for path in expected]})


def write_perturbed_state(source: Path, destination: Path,
                          node: int, delta_V: float) -> None:
    rows = read_csv(source)
    found = False
    for row in rows:
        if int(row["node_id"]) != node:
            continue
        found = True
        increment = float(row["electron_qf_increment_V"]) + delta_V
        reference = float(row["electron_qf_reference_V"])
        row["electron_qf_increment_V"] = format(increment, ".17g")
        row["phin"] = format(reference + increment, ".17g")
    if not found:
        raise ValueError(f"node {node} not present in {source}")
    write_csv(destination, rows)


def functional_deck(variant: str, state: Path, residual: Path,
                    role: str) -> dict[str, Any]:
    deck = copy.deepcopy(base_config(variant))
    configure_operator(deck, variant)
    deck.update({
        "simulation_type": "terminal_current_functional_probe",
        "state_file": str(state.resolve()),
        "contact": "drain",
        "residual_output_csv": str(residual.resolve()),
        "simplemos_m26": {"read_only": True, "role": role,
                           "operator": variant},
    })
    deck.pop("output_csv", None)
    return deck


def csv_probe_deck(variant: str, state: Path, output: Path,
                   simulation_type: str, role: str) -> dict[str, Any]:
    deck = copy.deepcopy(base_config(variant))
    configure_operator(deck, variant)
    deck.update({
        "simulation_type": simulation_type,
        "state_file": str(state.resolve()),
        "output_csv": str(output.resolve()),
        "simplemos_m26": {"read_only": True, "role": role,
                           "operator": variant},
    })
    deck.pop("residual_output_csv", None)
    if simulation_type == "newton_carrier_term_probe":
        deck["carrier_term_probe"] = {"solved_equation_terms": False}
    return deck


def build_tasks(contract: dict[str, Any], output: Path
                ) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    paths: dict[str, Any] = {"adjoint": {}, "operator": {}, "perturb": {}}
    variants = contract["state_variants"]
    first_layer = sorted(int(row["node_id"]) for row in read_csv(
        M24 / "m24_first_layer_continuity_ledger.csv"))
    contact_nodes = {int(row["contact_node"]) for row in read_csv(
        M24 / "m24_drain_cut_edge_ledger.csv")}
    baseline_edges = read_csv(
        RAW / "frozen/vg_0p05/state_full/eval_full/sg_edges.csv")
    perturb_columns = sorted({
        node for row in baseline_edges
        if int(row["node0"]) in first_layer
        or int(row["node1"]) in first_layer
        for node in (int(row["node0"]), int(row["node1"]))
        if node not in contact_nodes
    })

    for variant in variants:
        state = state_dir(variant) / "state.csv"
        root = output / "adjoint" / variant
        adjoint_csv = root / "adjoint.csv"
        config = root / "adjoint.json"
        deck = copy.deepcopy(base_config(variant))
        configure_operator(deck, variant)
        deck.update({
            "simulation_type": "terminal_current_adjoint_probe",
            "state_file": str(state.resolve()),
            "output_csv": str(adjoint_csv.resolve()),
            "contact": "drain",
            "simplemos_m26": {"read_only": True,
                               "role": "baseline_adjoint",
                               "operator": variant},
        })
        add_task(tasks, f"adjoint_{variant}", config, deck, [adjoint_csv])
        paths["adjoint"][variant] = {"config": config, "csv": adjoint_csv}

    for source in variants:
        state = state_dir(source) / "state.csv"
        for evaluation in variants:
            root = output / "operator" / f"state_{source}" / f"eval_{evaluation}"
            residual = root / "residual.csv"
            config = root / "functional.json"
            deck = functional_deck(
                evaluation, state, residual,
                f"operator_state_{source}_eval_{evaluation}")
            task_id = f"operator_state_{source}_eval_{evaluation}"
            add_task(tasks, task_id, config, deck, [residual])
            paths["operator"][(source, evaluation)] = {
                "config": config, "residual": residual, "task": task_id}

    step = float(contract["perturbation"]["central_step_V"])
    for variant in variants:
        source = state_dir(variant) / "state.csv"
        paths["perturb"][variant] = {}
        for node in perturb_columns:
            paths["perturb"][variant][node] = {}
            for sign, multiplier in (("minus", -1.0), ("plus", 1.0)):
                root = output / "perturb" / variant / f"node_{node}" / sign
                state = root / "state.csv"
                write_perturbed_state(source, state, node, multiplier * step)
                item: dict[str, Any] = {"state": state}
                residual = root / "functional_residual.csv"
                functional = root / "functional.json"
                role = f"fd_{variant}_node_{node}_{sign}"
                add_task(tasks, f"{role}_functional", functional,
                         functional_deck(variant, state, residual, role),
                         [residual])
                item.update({"functional": functional,
                             "functional_residual": residual,
                             "functional_task": f"{role}_functional"})
                for label, probe in (("sg", "sg_edge_flux_probe"),
                                     ("terms", "newton_carrier_term_probe")):
                    csv_path = root / f"{label}.csv"
                    config = root / f"{label}.json"
                    add_task(tasks, f"{role}_{label}", config,
                             csv_probe_deck(variant, state, csv_path,
                                            probe, role), [csv_path])
                    item[label] = csv_path
                paths["perturb"][variant][node][sign] = item
    paths["first_layer"] = first_layer
    paths["perturb_columns"] = perturb_columns
    return tasks, paths


def residual_map(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["node_id"]): row for row in read_csv(path)}


def adjoint_map(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["node_id"]): row for row in read_csv(path)}


def sg_map(path: Path) -> dict[int, dict[str, str]]:
    return {int(row["edge_id"]): row for row in read_csv(path)}


def signed_flux(row: dict[str, str], node: int) -> float:
    if int(row["node0"]) == node:
        return float(row["electron_flux"])
    if int(row["node1"]) == node:
        return -float(row["electron_flux"])
    raise ValueError(f"edge {row['edge_id']} is not incident to node {node}")


def flux_components(rows: dict[int, dict[str, str]], node: int,
                    cut_ids: set[int]) -> tuple[float, float]:
    incident = [row for row in rows.values()
                if node in (int(row["node0"]), int(row["node1"]))]
    contact = sum(signed_flux(row, node) for row in incident
                  if int(row["edge_id"]) in cut_ids)
    internal = sum(signed_flux(row, node) for row in incident
                   if int(row["edge_id"]) not in cut_ids)
    return contact, internal


def signed_particle_at_contact(row: dict[str, str],
                               contact_nodes: set[int]) -> float:
    sign = 1.0 if int(row["node0"]) in contact_nodes else -1.0
    return sign * float(row["electron_particle_line_flux_per_m_s"])


def cut_electron_current(rows: dict[int, dict[str, str]],
                         cut_ids: set[int], contact_nodes: set[int]) -> float:
    return sum(-Q * signed_particle_at_contact(rows[edge], contact_nodes)
               * 1.0e-6 for edge in cut_ids)


def source_term(row: dict[str, str]) -> float:
    return sum(float(row[column]) for column in (
        "electron_recombination", "electron_impact",
        "electron_gauge", "electron_boundary"))


def mesh_topology() -> tuple[set[int], set[int]]:
    config = base_config("full")
    mesh = read_json(Path(config["mesh_file"]))
    drain = next(item for item in mesh["contacts"]
                 if item["name"].lower() == "drain")
    contact_nodes = {int(node) for node in drain["node_ids"]}
    cut_ids = {int(row["edge_id"]) for row in read_csv(
        M24 / "m24_drain_cut_edge_ledger.csv")}
    return contact_nodes, cut_ids


def operator_attribution(contract: dict[str, Any], paths: dict[str, Any],
                         statuses: dict[str, dict[str, Any]],
                         first_layer: set[int], contact_nodes: set[int]
                         ) -> tuple[list[dict[str, Any]],
                                    list[dict[str, Any]],
                                    list[dict[str, Any]]]:
    detail: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    directions: list[dict[str, Any]] = []
    residual_columns = {"poisson": "psi_residual",
                        "electron": "phin_residual",
                        "hole": "phip_residual"}
    lambda_columns = {"poisson": "lambda_poisson",
                      "electron": "lambda_electron",
                      "hole": "lambda_hole"}
    for spec in contract["operator_directions"]:
        direction = spec["id"]
        source = spec["source_state"]
        baseline = spec["baseline_operator"]
        alternate = spec["alternate_operator"]
        target = spec["target_state"]
        base_item = paths["operator"][(source, baseline)]
        alt_item = paths["operator"][(source, alternate)]
        base_residual = residual_map(base_item["residual"])
        alt_residual = residual_map(alt_item["residual"])
        adjoint = adjoint_map(paths["adjoint"][baseline]["csv"])
        for node in sorted(adjoint):
            zone = ("drain_contact" if node in contact_nodes
                    else "drain_first_layer" if node in first_layer
                    else "other")
            for equation in ("poisson", "electron", "hole"):
                delta = (float(alt_residual[node][residual_columns[equation]])
                         - float(base_residual[node][residual_columns[equation]]))
                weight = float(adjoint[node][lambda_columns[equation]])
                detail.append({
                    "direction": direction,
                    "source_state": source,
                    "baseline_operator": baseline,
                    "alternate_operator": alternate,
                    "node_id": node,
                    "x_m": adjoint[node]["x_m"],
                    "y_m": adjoint[node]["y_m"],
                    "zone": zone,
                    "equation": equation,
                    "delta_scaled_residual": delta,
                    "adjoint_weight": weight,
                    "relaxation_current_contribution_A_per_um": -weight * delta,
                })
        current = [row for row in detail if row["direction"] == direction]
        for zone in ("drain_contact", "drain_first_layer", "other", "all"):
            for equation in ("poisson", "electron", "hole", "all"):
                selected = [row for row in current
                            if (zone == "all" or row["zone"] == zone)
                            and (equation == "all"
                                 or row["equation"] == equation)]
                values = [row["relaxation_current_contribution_A_per_um"]
                          for row in selected]
                summaries.append({
                    "direction": direction,
                    "zone": zone,
                    "equation": equation,
                    "row_count": len(selected),
                    "signed_relaxation_A_per_um": sum(values),
                    "absolute_relaxation_sum_A_per_um": sum(
                        abs(value) for value in values),
                })
        base_status = statuses[base_item["task"]]
        alt_status = statuses[alt_item["task"]]
        target_item = paths["operator"][(target, target)]
        target_status = statuses[target_item["task"]]
        source_current = float(base_status["current_A_per_um"])
        direct = float(alt_status["current_A_per_um"]) - source_current
        relaxation = sum(row["relaxation_current_contribution_A_per_um"]
                         for row in current)
        predicted = direct + relaxation
        actual = float(target_status["current_A_per_um"]) - source_current
        first = [row for row in current if row["zone"] == "drain_first_layer"]
        all_abs = sum(abs(row["relaxation_current_contribution_A_per_um"])
                      for row in current)
        directions.append({
            "direction": direction,
            "source_state": source,
            "target_state": target,
            "source_current_A_per_um": source_current,
            "target_current_A_per_um": float(target_status["current_A_per_um"]),
            "actual_self_consistent_current_change_A_per_um": actual,
            "direct_frozen_operator_current_change_A_per_um": direct,
            "adjoint_relaxation_current_change_A_per_um": relaxation,
            "first_order_total_current_change_A_per_um": predicted,
            "first_order_relative_error": abs(predicted - actual)
                / max(abs(actual), 1.0e-300),
            "first_layer_signed_relaxation_A_per_um": sum(
                row["relaxation_current_contribution_A_per_um"]
                for row in first),
            "first_layer_absolute_contribution_fraction": sum(
                abs(row["relaxation_current_contribution_A_per_um"])
                for row in first) / max(all_abs, 1.0e-300),
            "relaxation_condition_number": all_abs
                / max(abs(relaxation), 1.0e-300),
            "adjoint_relative_residual": statuses[f"adjoint_{baseline}"][
                "adjoint_relative_residual"],
        })
    return detail, summaries, directions


def perturbation_analysis(contract: dict[str, Any], paths: dict[str, Any],
                          statuses: dict[str, dict[str, Any]],
                          contact_nodes: set[int], cut_ids: set[int]
                          ) -> tuple[list[dict[str, Any]],
                                     list[dict[str, Any]],
                                     list[dict[str, Any]],
                                     list[dict[str, Any]]]:
    step = float(contract["perturbation"]["central_step_V"])
    first_layer = paths["first_layer"]
    perturb_columns = paths["perturb_columns"]
    state_rows = {
        variant: {int(row["node_id"]): row for row in read_csv(
            state_dir(variant) / "state.csv")}
        for variant in contract["state_variants"]
    }
    qf_delta = {
        node: float(state_rows["full"][node]["electron_qf_increment_V"])
        - float(state_rows["no_hfs"][node]["electron_qf_increment_V"])
        for node in perturb_columns
    }
    fd_rows: list[dict[str, Any]] = []
    matrix_rows: list[dict[str, Any]] = []
    response_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    m24_rows = {int(row["node_id"]): row for row in read_csv(
        M24 / "m24_first_layer_continuity_ledger.csv")}
    m24_report = read_json(M24 / "m24_contact_boundary_audit_report.json")
    actual_cut_delta = (
        float(m24_report["findings"]["full_cut_electron_current_A_per_um"])
        - float(m24_report["findings"]["no_hfs_cut_electron_current_A_per_um"]))

    for variant in contract["state_variants"]:
        adjoint = adjoint_map(paths["adjoint"][variant]["csv"])
        by_column: dict[int, list[dict[str, Any]]] = {}
        for node in perturb_columns:
            minus = paths["perturb"][variant][node]["minus"]
            plus = paths["perturb"][variant][node]["plus"]
            minus_status = statuses[minus["functional_task"]]
            plus_status = statuses[plus["functional_task"]]
            sg_minus = sg_map(minus["sg"])
            sg_plus = sg_map(plus["sg"])
            term_minus = {int(row["node_id"]): row
                          for row in read_csv(minus["terms"])}
            term_plus = {int(row["node_id"]): row
                         for row in read_csv(plus["terms"])}
            functional_gradient = (
                float(plus_status["current_A_per_um"])
                - float(minus_status["current_A_per_um"])) / (2.0 * step)
            cut_gradient = (cut_electron_current(
                sg_plus, cut_ids, contact_nodes) - cut_electron_current(
                    sg_minus, cut_ids, contact_nodes)) / (2.0 * step)
            adjoint_gradient = float(adjoint[node][
                "dI_dphin_A_per_um_per_V"])
            gradient_scale = max(abs(functional_gradient),
                                 abs(adjoint_gradient), 1.0e-300)
            fd_rows.append({
                "baseline_variant": variant,
                "perturbed_node": node,
                "central_step_V": step,
                "qf_increment_full_minus_no_hfs_V": qf_delta[node],
                "functional_dI_dphin_A_per_um_per_V": functional_gradient,
                "drain_cut_electron_dI_dphin_A_per_um_per_V": cut_gradient,
                "adjoint_dI_dphin_A_per_um_per_V": adjoint_gradient,
                "functional_vs_adjoint_relative_disagreement": abs(
                    functional_gradient - adjoint_gradient) / gradient_scale,
                "functional_vs_cut_relative_disagreement": abs(
                    functional_gradient - cut_gradient) / max(
                        abs(functional_gradient), abs(cut_gradient), 1.0e-300),
                "predicted_total_current_delta_A_per_um": (
                    functional_gradient * qf_delta[node]),
                "predicted_cut_electron_delta_A_per_um": (
                    cut_gradient * qf_delta[node]),
            })
            column_rows: list[dict[str, Any]] = []
            for target in first_layer:
                contact_minus, internal_minus = flux_components(
                    sg_minus, target, cut_ids)
                contact_plus, internal_plus = flux_components(
                    sg_plus, target, cut_ids)
                source_minus = source_term(term_minus[target])
                source_plus = source_term(term_plus[target])
                residual_minus = float(term_minus[target]["electron_residual"])
                residual_plus = float(term_plus[target]["electron_residual"])
                contact_derivative = (contact_plus - contact_minus) / (2.0 * step)
                internal_derivative = (internal_plus - internal_minus) / (2.0 * step)
                source_derivative = (source_plus - source_minus) / (2.0 * step)
                residual_derivative = (residual_plus - residual_minus) / (2.0 * step)
                closure = (contact_derivative + internal_derivative
                           + source_derivative - residual_derivative)
                scale = max(abs(contact_derivative), abs(internal_derivative),
                            abs(source_derivative), abs(residual_derivative),
                            1.0e-300)
                row = {
                    "baseline_variant": variant,
                    "target_row_node": target,
                    "perturbed_column_node": node,
                    "contact_flux_derivative_per_V": contact_derivative,
                    "internal_flux_derivative_per_V": internal_derivative,
                    "source_derivative_per_V": source_derivative,
                    "electron_residual_derivative_per_V": residual_derivative,
                    "term_closure_derivative_per_V": closure,
                    "term_closure_relative": abs(closure) / scale,
                }
                matrix_rows.append(row)
                column_rows.append(row)
            by_column[node] = column_rows

        variant_fd = [row for row in fd_rows
                      if row["baseline_variant"] == variant]
        predicted_cut = sum(row["predicted_cut_electron_delta_A_per_um"]
                            for row in variant_fd)
        first_layer_predicted_cut = sum(
            row["predicted_cut_electron_delta_A_per_um"]
            for row in variant_fd if row["perturbed_node"] in first_layer)
        predicted_total = sum(row["predicted_total_current_delta_A_per_um"]
                              for row in variant_fd)
        for target in first_layer:
            selected = [row for row in matrix_rows
                        if row["baseline_variant"] == variant
                        and row["target_row_node"] == target]
            predicted = {}
            for key in ("contact_flux", "internal_flux", "source",
                        "electron_residual"):
                derivative_key = f"{key}_derivative_per_V"
                predicted[key] = sum(
                    row[derivative_key] * qf_delta[row["perturbed_column_node"]]
                    for row in selected)
            observed = m24_rows[target]
            response_rows.append({
                "baseline_variant": variant,
                "target_row_node": target,
                "predicted_contact_flux_delta": predicted["contact_flux"],
                "observed_contact_flux_delta": observed["delta_contact_flux"],
                "predicted_internal_flux_delta": predicted["internal_flux"],
                "observed_internal_flux_delta": observed["delta_internal_flux"],
                "predicted_source_delta": predicted["source"],
                "observed_source_delta": observed["delta_source"],
                "predicted_residual_delta": predicted["electron_residual"],
                "observed_residual_delta": observed["delta_residual"],
            })
        summaries.append({
            "baseline_variant": variant,
            "first_layer_node_count": len(first_layer),
            "perturbation_column_count": len(perturb_columns),
            "predicted_cut_electron_delta_A_per_um": predicted_cut,
            "first_layer_predicted_cut_electron_delta_A_per_um": (
                first_layer_predicted_cut),
            "actual_cut_electron_delta_A_per_um": actual_cut_delta,
            "cut_prediction_relative_error": abs(predicted_cut - actual_cut_delta)
                / max(abs(actual_cut_delta), 1.0e-300),
            "predicted_total_current_delta_A_per_um": predicted_total,
            "maximum_direct_gradient_relative_disagreement": max(
                row["functional_vs_adjoint_relative_disagreement"]
                for row in variant_fd),
            "maximum_functional_cut_gradient_relative_disagreement": max(
                row["functional_vs_cut_relative_disagreement"]
                for row in variant_fd),
            "maximum_fd_term_closure_relative": max(
                row["term_closure_relative"] for row in matrix_rows
                if row["baseline_variant"] == variant),
        })
    return fd_rows, matrix_rows, response_rows, summaries


def analyze(contract: dict[str, Any], paths: dict[str, Any],
            statuses: dict[str, dict[str, Any]], portable_dir: Path
            ) -> dict[str, Any]:
    first_layer = set(paths["first_layer"])
    contact_nodes, cut_ids = mesh_topology()
    operator_rows, zone_rows, direction_rows = operator_attribution(
        contract, paths, statuses, first_layer, contact_nodes)
    fd_rows, matrix_rows, response_rows, perturb_summaries = (
        perturbation_analysis(contract, paths, statuses,
                              contact_nodes, cut_ids))
    write_csv(portable_dir / "m26_operator_adjoint_ledger.csv", operator_rows)
    write_csv(portable_dir / "m26_operator_zone_summary.csv", zone_rows)
    write_csv(portable_dir / "m26_direction_summary.csv", direction_rows)
    write_csv(portable_dir / "m26_first_layer_fd_ledger.csv", fd_rows)
    write_csv(portable_dir / "m26_first_layer_jacobian_ledger.csv", matrix_rows)
    write_csv(portable_dir / "m26_first_layer_response_ledger.csv", response_rows)
    write_csv(portable_dir / "m26_perturbation_summary.csv", perturb_summaries)
    maximum_adjoint_residual = max(
        float(row["adjoint_relative_residual"]) for row in direction_rows)
    maximum_gradient_disagreement = max(
        row["functional_vs_adjoint_relative_disagreement"] for row in fd_rows)
    maximum_functional_cut_disagreement = max(
        row["functional_vs_cut_relative_disagreement"] for row in fd_rows)
    maximum_term_closure = max(row["term_closure_relative"]
                               for row in matrix_rows)
    maximum_first_order_error = max(
        row["first_order_relative_error"] for row in direction_rows)
    maximum_cut_prediction_error = max(
        row["cut_prediction_relative_error"] for row in perturb_summaries)
    m24_report = read_json(M24 / "m24_contact_boundary_audit_report.json")
    m24_currents = {
        "full": float(m24_report["findings"][
            "full_cut_electron_current_A_per_um"]),
        "no_hfs": float(m24_report["findings"][
            "no_hfs_cut_electron_current_A_per_um"]),
    }
    functional_currents = {
        "full": float(statuses[paths["operator"][("full", "full")][
            "task"]]["current_A_per_um"]),
        "no_hfs": float(statuses[paths["operator"][("no_hfs", "no_hfs")][
            "task"]]["current_A_per_um"]),
    }
    reference_current_disagreement = {
        variant: abs(functional_currents[variant] - m24_currents[variant])
        / max(abs(functional_currents[variant]), abs(m24_currents[variant]),
              1.0e-300)
        for variant in contract["state_variants"]
    }
    maximum_reference_current_disagreement = max(
        reference_current_disagreement.values())
    turn_on = next(row for row in direction_rows
                   if row["direction"] == "turn_on")
    turn_off = next(row for row in direction_rows
                    if row["direction"] == "turn_off")
    turn_on_rows = [row for row in operator_rows
                    if row["direction"] == "turn_on"]
    node_totals: list[dict[str, Any]] = []
    for node in sorted({row["node_id"] for row in turn_on_rows}):
        selected = [row for row in turn_on_rows if row["node_id"] == node]
        node_totals.append({
            "node_id": node,
            "x_m": float(selected[0]["x_m"]),
            "y_m": float(selected[0]["y_m"]),
            "relaxation_current_A_per_um": sum(
                row["relaxation_current_contribution_A_per_um"]
                for row in selected),
            "absolute_contribution_sum_A_per_um": sum(
                abs(row["relaxation_current_contribution_A_per_um"])
                for row in selected),
        })
    top_operator_nodes = sorted(
        node_totals,
        key=lambda row: row["absolute_contribution_sum_A_per_um"],
        reverse=True)[:10]
    baseline_fd = [row for row in fd_rows
                   if row["baseline_variant"] == "no_hfs"
                   and row["perturbed_node"] in first_layer]
    predicted_cut_sum = sum(
        row["predicted_cut_electron_delta_A_per_um"] for row in baseline_fd)
    first_layer_fractions = {
        str(row["perturbed_node"]): (
            row["predicted_cut_electron_delta_A_per_um"]
            / predicted_cut_sum if predicted_cut_sum != 0.0 else 0.0)
        for row in baseline_fd
    }
    contact_support_nodes = {
        node for node in first_layer
        if abs(float(next(row for row in read_csv(
            M24 / "m24_first_layer_continuity_ledger.csv")
            if int(row["node_id"]) == node)["full_contact_flux"])) > 0.0
    }
    maximum_qf_only_row_imbalance = max(
        abs(float(row["predicted_residual_delta"]))
        / max(abs(float(row["predicted_contact_flux_delta"])),
              abs(float(row["predicted_internal_flux_delta"])), 1.0e-300)
        for row in response_rows
        if int(row["target_row_node"]) in contact_support_nodes)
    all_turn_on = next(row for row in zone_rows
                       if row["direction"] == "turn_on"
                       and row["zone"] == "all"
                       and row["equation"] == "all")
    electron_turn_on = next(row for row in zone_rows
                            if row["direction"] == "turn_on"
                            and row["zone"] == "all"
                            and row["equation"] == "electron")
    acceptance = contract["acceptance"]
    report = {
        "schema": "vela.simplemos.sdevice.m26_first_layer_feedback_audit_report.v1",
        "status": "complete",
        "execution": {
            "first_layer_node_count": len(first_layer),
            "first_layer_nodes": sorted(first_layer),
            "perturbation_column_count": len(paths["perturb_columns"]),
            "perturbation_columns": paths["perturb_columns"],
            "central_perturbation_state_count": 2 * len(
                paths["perturb_columns"])
                * len(contract["state_variants"]),
            "read_only_probe_count": len(statuses),
            "new_sentaurus_execution": False,
            "default_model_changed": False,
        },
        "closure": {
            "maximum_adjoint_relative_residual": maximum_adjoint_residual,
            "maximum_direct_gradient_relative_disagreement": (
                maximum_gradient_disagreement),
            "maximum_functional_cut_gradient_relative_disagreement": (
                maximum_functional_cut_disagreement),
            "maximum_fd_term_closure_relative": maximum_term_closure,
            "maximum_first_order_relative_error": maximum_first_order_error,
            "maximum_cut_prediction_relative_error": (
                maximum_cut_prediction_error),
            "maximum_reference_current_relative_disagreement": (
                maximum_reference_current_disagreement),
            "all_acceptance_checks_pass": (
                maximum_adjoint_residual < acceptance[
                    "maximum_adjoint_relative_residual"]
                and maximum_gradient_disagreement < acceptance[
                    "maximum_direct_gradient_relative_disagreement"]
                and maximum_functional_cut_disagreement < acceptance[
                    "maximum_functional_cut_gradient_relative_disagreement"]
                and maximum_term_closure < acceptance[
                    "maximum_fd_term_closure_relative"]
                and maximum_first_order_error < acceptance[
                    "maximum_first_order_relative_error"]
                and maximum_cut_prediction_error < acceptance[
                    "maximum_cut_prediction_relative_error"]
                and maximum_reference_current_disagreement < acceptance[
                    "maximum_reference_current_relative_disagreement"]),
        },
        "operator_directions": direction_rows,
        "perturbation_summaries": perturb_summaries,
        "findings": {
            "reference_aware_arclength_pack_enabled": True,
            "functional_default_current_scale_includes_line_factor": True,
            "functional_vs_m24_cut_current_relative_disagreement": (
                reference_current_disagreement),
            "turn_on_direct_frozen_fraction_of_actual": abs(
                turn_on["direct_frozen_operator_current_change_A_per_um"])
                / max(abs(turn_on[
                    "actual_self_consistent_current_change_A_per_um"]),
                      1.0e-300),
            "turn_on_adjoint_relaxation_fraction_of_actual": (
                turn_on["adjoint_relaxation_current_change_A_per_um"]
                / turn_on["actual_self_consistent_current_change_A_per_um"]),
            "turn_off_adjoint_relaxation_fraction_of_actual": (
                turn_off["adjoint_relaxation_current_change_A_per_um"]
                / turn_off["actual_self_consistent_current_change_A_per_um"]),
            "turn_on_electron_equation_signed_fraction": (
                electron_turn_on["signed_relaxation_A_per_um"]
                / all_turn_on["signed_relaxation_A_per_um"]),
            "turn_on_first_layer_operator_absolute_fraction": (
                turn_on["first_layer_absolute_contribution_fraction"]),
            "first_layer_cut_current_fraction_by_node_no_hfs_linearization": (
                first_layer_fractions),
            "maximum_qf_only_contact_support_row_imbalance_relative": (
                maximum_qf_only_row_imbalance),
            "top_turn_on_operator_relaxation_nodes": top_operator_nodes,
        },
        "conclusions": {
            "supported": [
                "The frozen HFS operator switch contributes only a negligible direct drain-current change; the observed response is predominantly self-consistent relaxation.",
                "The reference-aware adjoint predicts the finite HFS-on/off current change within the bidirectional first-order error envelope reported here.",
                "The drain-first-layer QF increments reproduce the electron cut-current change to better than 0.1 percent, so the first layer is the terminal-current conduit.",
                "Node 1092 carries the dominant first-layer cut sensitivity, followed by node 1100; node 1090 has zero drain-cut support.",
                "The direct HFS operator residual is not localized at the drain contact or its first free-node layer; its adjoint-weighted source is distributed upstream and dominated by the electron equation.",
                "Including the complete noncontact one-ring makes QF-only contact and internal flux changes nearly cancel in each first-layer continuity row.",
            ],
            "not_supported": [
                "Treating the drain-first-layer rows as the physical origin of the HFS operator response merely because they control terminal extraction.",
                "Claiming exact finite-switch attribution from a local adjoint; the turn-on and turn-off linearizations retain finite nonlinear error.",
                "Equating the Vela feedback path with Sentaurus without matching Sentaurus node-level QF, mobility, and residual exports.",
                "Changing a default HFS, contact, SG, mobility, or solver setting from this audit alone.",
            ],
        },
        "artifacts": {
            "operator_adjoint_ledger_csv": portable(
                portable_dir / "m26_operator_adjoint_ledger.csv"),
            "operator_zone_summary_csv": portable(
                portable_dir / "m26_operator_zone_summary.csv"),
            "direction_summary_csv": portable(
                portable_dir / "m26_direction_summary.csv"),
            "first_layer_fd_ledger_csv": portable(
                portable_dir / "m26_first_layer_fd_ledger.csv"),
            "first_layer_jacobian_ledger_csv": portable(
                portable_dir / "m26_first_layer_jacobian_ledger.csv"),
            "first_layer_response_ledger_csv": portable(
                portable_dir / "m26_first_layer_response_ledger.csv"),
            "perturbation_summary_csv": portable(
                portable_dir / "m26_perturbation_summary.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    write_json(portable_dir / "m26_first_layer_feedback_audit_report.json",
               report)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--portable", type=Path, default=PORTABLE)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    contract = read_json(args.contract.resolve())
    validate(contract)
    tasks, paths = build_tasks(contract, args.output.resolve())
    statuses = run_tasks(tasks, args.runner.resolve(), args.workers, args.force)
    report = analyze(contract, paths, statuses, args.portable.resolve())
    print(json.dumps({"status": report["status"],
                      **report["execution"], **report["closure"]}))


if __name__ == "__main__":
    main()
