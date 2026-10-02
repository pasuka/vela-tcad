#!/usr/bin/env python3
"""Run the SimpleMOS M30 double-off causal closure audit."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import run_simplemos_m10_fixed_state_replay as m10


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m30_double_off_causal_closure_contract_v1.json")
M29 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m29_bgn_srh_factorial")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m30_double_off_causal_closure")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "double_off_causal_closure")
RUNNER = REPO / "build-release/vela_example_runner.exe"
Q = 1.602176634e-19
CELLS = ("bgn_on_srh_on", "bgn_on_srh_off",
         "bgn_off_srh_on", "bgn_off_srh_off")
PROBES = {
    "sg": "sg_edge_flux_probe",
    "terms": "newton_carrier_term_probe",
    "rows": "newton_carrier_row_probe",
}


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


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


def run_config(config: Path, runner: Path) -> dict[str, Any]:
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run([str(runner), "--config", str(config)],
                               cwd=REPO, text=True, capture_output=True,
                               env=env, check=False)
    config.with_suffix(".stdout.txt").write_text(
        completed.stdout, encoding="utf-8", newline="\n")
    config.with_suffix(".stderr.txt").write_text(
        completed.stderr, encoding="utf-8", newline="\n")
    if completed.returncode:
        raise RuntimeError(f"M30 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    status = json.loads(completed.stdout.strip().splitlines()[-1])
    write_json(config.with_suffix(".result.json"), status)
    return status


def set_biases(config: dict[str, Any]) -> None:
    for contact in config["contacts"]:
        name = contact["name"].lower()
        if name == "drain":
            contact["bias"] = 0.05
        elif name == "gate":
            contact["bias"] = 0.05


def prepare(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    sent_manifest = read_json(M29 / "sentaurus_export_manifest.json")
    vela_manifest = read_json(M29 / "vela_manifest.json")
    sent = {row["cell"]: row for row in sent_manifest["states"]}
    vela = {row["cell"]: row for row in vela_manifest["states"]}
    state_rows = []
    for cell in CELLS:
        vela_state = REPO / vela[cell]["state"]
        sent_state = OUTPUT / "states" / f"{cell}_sentaurus.csv"
        if force or not sent_state.is_file():
            m10.make_sentaurus_state(REPO / sent[cell]["export_dir"],
                                      vela_state, sent_state)
        state_rows.extend([
            {"solver": "vela", "cell": cell,
             "state_file": portable(vela_state),
             "native_current_A_per_um": vela[cell]["current_A_per_um"]},
            {"solver": "sentaurus", "cell": cell,
             "state_file": portable(sent_state),
             "native_current_A_per_um":
                 sent[cell]["terminal"]["drain_total_current_A_per_um"]},
        ])
    manifest = {
        "schema": "vela.simplemos.sdevice.m30_state_manifest.v1",
        "status": "prepared", "states": state_rows,
    }
    write_json(OUTPUT / "state_manifest.json", manifest)
    return manifest


def run_probes(contract: dict[str, Any], states: dict[str, Any],
               runner: Path, workers: int, force: bool) -> dict[str, Any]:
    state_map = {(row["solver"], row["cell"]): row
                 for row in states["states"]}
    jobs: list[tuple[Path, str]] = []
    pairs = []
    for state_solver in ("vela", "sentaurus"):
        for state_cell in CELLS:
            state = REPO / state_map[(state_solver, state_cell)]["state_file"]
            for operator_cell in CELLS:
                pair_id = f"{state_solver}_{state_cell}__op_{operator_cell}"
                pair_dir = OUTPUT / "probes" / pair_id
                pair_dir.mkdir(parents=True, exist_ok=True)
                base = read_json(M29 / "vela" / operator_cell
                                 / "20_gate/config.json")
                paths: dict[str, str] = {}
                for label, simulation_type in PROBES.items():
                    config_path = pair_dir / f"{label}.json"
                    output_path = pair_dir / f"{label}.csv"
                    deck = copy.deepcopy(base)
                    deck.update({
                        "simulation_type": simulation_type,
                        "state_file": str(state.resolve()),
                        "output_csv": str(output_path.resolve()),
                        "simplemos_m30": {
                            "read_only": True,
                            "state_solver": state_solver,
                            "state_cell": state_cell,
                            "operator_cell": operator_cell,
                        },
                    })
                    set_biases(deck)
                    deck.pop("sweep", None)
                    deck.pop("log_file", None)
                    if label == "terms":
                        deck["carrier_term_probe"] = {
                            "solved_equation_terms": False}
                    write_json(config_path, deck)
                    paths[label] = portable(output_path)
                    result_path = config_path.with_suffix(".result.json")
                    if force or not output_path.is_file() or not result_path.is_file():
                        jobs.append((config_path, f"{pair_id}:{label}"))
                pairs.append({
                    "pair_id": pair_id,
                    "state_solver": state_solver,
                    "state_cell": state_cell,
                    "operator_cell": operator_cell,
                    "state_file": portable(state),
                    "native_current_A_per_um":
                        state_map[(state_solver, state_cell)]["native_current_A_per_um"],
                    "outputs": paths,
                })
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_config, config, runner): label
                   for config, label in jobs}
        for future in as_completed(futures):
            label = futures[future]
            future.result()
            print(f"completed {label}", flush=True)
    result = {
        "schema": "vela.simplemos.sdevice.m30_probe_manifest.v1",
        "status": "complete", "probe_run_count": len(pairs) * len(PROBES),
        "pairs": pairs,
    }
    write_json(OUTPUT / "probe_manifest.json", result)
    return result


def contact_nodes(mesh: dict[str, Any], name: str) -> set[int]:
    matches = [row for row in mesh["contacts"]
               if row["name"].lower() == name]
    if len(matches) != 1:
        raise ValueError(f"expected one {name} contact")
    return {int(node) for node in matches[0]["node_ids"]}


def cut_metrics(edges: list[dict[str, str]], nodes: set[int]) -> dict[str, float]:
    electron_edges, hole_edges = [], []
    for row in edges:
        n0, n1 = int(row["node0"]), int(row["node1"])
        if (n0 in nodes) == (n1 in nodes):
            continue
        outward = 1.0 if n0 in nodes else -1.0
        electron_edges.append(-Q * outward
                              * float(row["electron_particle_line_flux_per_m_s"])
                              * 1.0e-6)
        hole_edges.append(-Q * outward
                          * float(row["hole_particle_line_flux_per_m_s"])
                          * 1.0e-6)
    electron = sum(electron_edges)
    hole = sum(hole_edges)
    total = electron - hole
    return {
        "edge_count": len(electron_edges),
        "electron_A_per_um": electron,
        "hole_A_per_um": hole,
        "total_A_per_um": total,
        "electron_condition": sum(map(abs, electron_edges))
            / max(abs(electron), 1.0e-300),
        "hole_condition": sum(map(abs, hole_edges))
            / max(abs(hole), 1.0e-300),
        "total_condition": (sum(map(abs, electron_edges))
                            + sum(map(abs, hole_edges)))
            / max(abs(total), 1.0e-300),
    }


def run_solver_controls(contract: dict[str, Any], runner: Path,
                        force: bool) -> dict[str, Any]:
    source_config = M29 / "vela/bgn_off_srh_off/20_gate/config.json"
    reclosure_state = M29 / "vela/bgn_off_srh_off/20_gate/state.csv"
    path_state = M29 / "vela/bgn_off_srh_off/10_drain/state.csv"
    results = []
    for control in contract["solver_controls"]:
        run_dir = OUTPUT / "solver_controls" / control
        run_dir.mkdir(parents=True, exist_ok=True)
        config_path = run_dir / "config.json"
        curve = run_dir / "curve.csv"
        state = run_dir / "state.csv"
        deck = read_json(source_config)
        deck["output_csv"] = str(curve.resolve())
        deck["log_file"] = str((run_dir / "run.log").resolve())
        path_repeat = control.endswith("_path_repeat")
        deck["sweep"].update({
            "start": 0.0 if path_repeat else 0.05,
            "stop": 0.05,
            "step": 0.05,
            "bias_points": [0.0, 0.05] if path_repeat else [0.05],
            "max_retries": 0,
            "initial_state_file": str(
                (path_state if path_repeat else reclosure_state).resolve()),
            "write_state_file": str(state.resolve()),
            "write_state_every_point_prefix": str((run_dir / "accepted_state").resolve()),
        })
        diagnostics = deck["sweep"].setdefault("diagnostics", {})
        diagnostics["terminal_balance"] = {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": str((run_dir / "terminal_balance.csv").resolve()),
        }
        diagnostics["newton_history"] = {
            "enabled": True,
            "csv_file": str((run_dir / "newton_history.csv").resolve()),
        }
        row_control = deck["solver"].setdefault("carrier_row_convergence", {})
        row_control.update({
            "mode": "report",
            "eps_row": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale_fraction": 0.0,
            "min_source_scale": 1.0e-18,
            "diagnostic_csv": str((run_dir / "carrier_row_violations.csv").resolve()),
        })
        deck["solver"]["global_continuity_closure"] = {
            "mode": "report", "tolerance": 0.01, "source_floor": 1.0e-18}
        if control != "baseline_report_reclosure":
            deck["solver"].update({
                "reltol": 1.0e-10, "abstol": 1.0e-15, "max_iter": 200})
        if control in ("strict_enforce_reclosure", "strict_row_scaled_reclosure",
                       "strict_path_repeat"):
            row_control["mode"] = "enforce"
            deck["solver"]["global_continuity_closure"]["mode"] = "enforce"
            deck["solver"]["carrier_row_qualified_stall_acceptance"] = True
        if control == "strict_row_scaled_reclosure":
            deck["solver"]["continuity_row_scaling"] = {
                "enabled": True,
                "flux_fraction": 0.0,
                "scale_floor": 1.0e-30,
                "min_source_scale": 1.0e-18,
                "min_weight": 1.0e-12,
                "max_weight": 1.0e18,
            }
        deck["simplemos_m30_control"] = {
            "read_only": False, "control": control,
            "same_bias_reclosure": not path_repeat,
            "continuation_path_repeat": path_repeat,
        }
        write_json(config_path, deck)
        stdout = config_path.with_suffix(".stdout.txt")
        stderr = config_path.with_suffix(".stderr.txt")
        if force or not stdout.is_file():
            env = os.environ.copy()
            env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
            completed = subprocess.run([str(runner), "--config", str(config_path)],
                                       cwd=REPO, text=True, capture_output=True,
                                       env=env, check=False)
            stdout.write_text(completed.stdout, encoding="utf-8", newline="\n")
            stderr.write_text(completed.stderr, encoding="utf-8", newline="\n")
            returncode = completed.returncode
        else:
            returncode = 0 if state.is_file() and curve.is_file() else 1
        current = None
        if curve.is_file():
            curve_rows = read_csv(curve)
            if curve_rows:
                current = abs(float(curve_rows[-1]["current_total_A_per_um"]))
        sg_path = run_dir / "sg.csv"
        if returncode == 0 and state.is_file():
            sg_config = copy.deepcopy(deck)
            sg_config.update({
                "simulation_type": "sg_edge_flux_probe",
                "state_file": str(state.resolve()),
                "output_csv": str(sg_path.resolve()),
            })
            sg_config.pop("sweep", None)
            sg_config.pop("log_file", None)
            sg_config_path = run_dir / "sg.json"
            write_json(sg_config_path, sg_config)
            if force or not sg_path.is_file():
                run_config(sg_config_path, runner)
        failure_path = curve.with_name(
            curve.stem + "_newton_failure_diagnostics.json")
        failure = None
        if failure_path.is_file():
            failure_rows = read_json(failure_path)
            if failure_rows:
                failure = failure_rows[-1]
        block = failure.get("block_residuals", {}) if failure else {}
        results.append({
            "control": control,
            "returncode": returncode,
            "converged": returncode == 0 and state.is_file(),
            "current_A_per_um": current,
            "same_bias_reclosure": not path_repeat,
            "continuation_path_repeat": path_repeat,
            "failure_reason": failure.get("failure_reason") if failure else None,
            "failed_iteration": failure.get("failed_iteration") if failure else None,
            "failure_residual_norm": failure.get("residual_norm") if failure else None,
            "failure_psi_residual": block.get("psi"),
            "failure_phin_residual": block.get("phin"),
            "failure_phip_residual": block.get("phip"),
            "config": portable(config_path),
            "state": portable(state) if state.is_file() else None,
            "sg": portable(sg_path) if sg_path.is_file() else None,
            "stdout": portable(stdout),
            "stderr": portable(stderr),
            "failure_diagnostics": portable(failure_path)
                if failure_path.is_file() else None,
        })
    manifest = {
        "schema": "vela.simplemos.sdevice.m30_solver_controls.v1",
        "status": "complete", "controls": results,
    }
    write_json(OUTPUT / "solver_control_manifest.json", manifest)
    return manifest


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def state_fields(path: Path) -> dict[str, dict[int, float]]:
    rows = read_csv(path)
    result: dict[str, dict[int, float]] = {}
    for field in ("psi", "phin", "phip", "electrons_m3", "holes_m3"):
        result[field] = {int(row["node_id"]): float(row[field]) for row in rows}
    result["logn"] = {node: math.log10(max(value, 1e-300))
                      for node, value in result["electrons_m3"].items()}
    result["logp"] = {node: math.log10(max(value, 1e-300))
                      for node, value in result["holes_m3"].items()}
    return result


def interaction(fields: dict[str, dict[str, dict[int, float]]], field: str,
                node: int) -> float:
    return (fields["bgn_on_srh_on"][field][node]
            - fields["bgn_on_srh_off"][field][node]
            - fields["bgn_off_srh_on"][field][node]
            + fields["bgn_off_srh_off"][field][node])


def percentile(values: list[float], fraction: float) -> float:
    return float(m10.percentile(values, fraction))


def pearson(left: list[float], right: list[float]) -> float:
    lm, rm = sum(left) / len(left), sum(right) / len(right)
    numerator = sum((a - lm) * (b - rm) for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - lm) ** 2 for a in left)
                            * sum((b - rm) ** 2 for b in right))
    return numerator / denominator if denominator else 0.0


def log_ratio(left: float, right: float) -> float:
    return math.log10(max(abs(left), 1.0e-300)
                      / max(abs(right), 1.0e-300))


def analyze(contract: dict[str, Any], states: dict[str, Any],
            probes: dict[str, Any], controls: dict[str, Any]) -> dict[str, Any]:
    state_map = {(row["solver"], row["cell"]): row
                 for row in states["states"]}
    mesh = read_json(M29 / "vela" / CELLS[0] / "20_gate/config.json")
    mesh = read_json(Path(mesh["mesh_file"]))
    source_nodes = contact_nodes(mesh, "source")
    drain_nodes = contact_nodes(mesh, "drain")
    substrate_nodes = contact_nodes(mesh, "substrate")
    gate_nodes = contact_nodes(mesh, "gate")
    all_contact_nodes = {int(node) for contact in mesh["contacts"]
                         for node in contact["node_ids"]}
    pair_rows = []
    for pair in probes["pairs"]:
        edges = read_csv(REPO / pair["outputs"]["sg"])
        terms = read_csv(REPO / pair["outputs"]["terms"])
        rows = read_csv(REPO / pair["outputs"]["rows"])
        source = cut_metrics(edges, source_nodes)
        drain = cut_metrics(edges, drain_nodes)
        substrate = cut_metrics(edges, substrate_nodes)
        gate = cut_metrics(edges, gate_nodes)
        all_terminal_total = sum(item["total_A_per_um"] for item in (
            source, drain, substrate, gate))
        terminal_balance_condition = sum(abs(item["total_A_per_um"])
                                         for item in (source, drain, substrate, gate)) \
            / max(abs(all_terminal_total), 1.0e-300)
        free_terms = [row for row in terms
                      if int(row["node_id"]) not in all_contact_nodes]
        free_rows = [row for row in rows
                     if int(row["node_id"]) not in all_contact_nodes]
        electron_residuals = [float(row["electron_residual"])
                              for row in free_rows]
        hole_residuals = [float(row["hole_residual"])
                          for row in free_rows]
        native = float(pair["native_current_A_per_um"])
        pair_rows.append({
            "pair_id": pair["pair_id"],
            "state_solver": pair["state_solver"],
            "state_cell": pair["state_cell"],
            "operator_cell": pair["operator_cell"],
            "is_diagonal": pair["state_cell"] == pair["operator_cell"],
            "native_current_A_per_um": native,
            "drain_replayed_current_A_per_um": drain["total_A_per_um"],
            "drain_replay_vs_native_dex": log_ratio(drain["total_A_per_um"], native),
            "source_total_current_A_per_um": source["total_A_per_um"],
            "substrate_total_current_A_per_um": substrate["total_A_per_um"],
            "gate_total_current_A_per_um": gate["total_A_per_um"],
            "drain_electron_current_A_per_um": drain["electron_A_per_um"],
            "drain_hole_current_A_per_um": drain["hole_A_per_um"],
            "source_electron_current_A_per_um": source["electron_A_per_um"],
            "source_hole_current_A_per_um": source["hole_A_per_um"],
            "source_drain_total_imbalance_A_per_um":
                source["total_A_per_um"] + drain["total_A_per_um"],
            "all_terminal_total_imbalance_A_per_um": all_terminal_total,
            "all_terminal_balance_condition": terminal_balance_condition,
            "source_electron_condition": source["electron_condition"],
            "source_hole_condition": source["hole_condition"],
            "source_total_condition": source["total_condition"],
            "drain_electron_condition": drain["electron_condition"],
            "drain_hole_condition": drain["hole_condition"],
            "drain_total_condition": drain["total_condition"],
            "source_cut_edge_count": int(source["edge_count"]),
            "drain_cut_edge_count": int(drain["edge_count"]),
            "substrate_cut_edge_count": int(substrate["edge_count"]),
            "gate_cut_edge_count": int(gate["edge_count"]),
            "free_electron_recombination_sum": sum(
                float(row["electron_recombination"]) for row in free_terms),
            "free_hole_recombination_sum": sum(
                float(row["hole_recombination"]) for row in free_terms),
            "free_electron_residual_l2": l2(electron_residuals),
            "free_hole_residual_l2": l2(hole_residuals),
            "free_electron_residual_max_abs": max(map(abs, electron_residuals)),
            "free_hole_residual_max_abs": max(map(abs, hole_residuals)),
        })
    write_csv(PORTABLE / "m30_state_operator_matrix.csv", pair_rows)
    indexed = {(row["state_solver"], row["state_cell"], row["operator_cell"]): row
               for row in pair_rows}
    decomposition_rows = []
    for cell in CELLS:
        vela_row = indexed[("vela", cell, cell)]
        sent_row = indexed[("sentaurus", cell, cell)]
        native_gap = log_ratio(vela_row["native_current_A_per_um"],
                               sent_row["native_current_A_per_um"])
        state_contribution = log_ratio(
            vela_row["drain_replayed_current_A_per_um"],
            sent_row["drain_replayed_current_A_per_um"])
        sent_operator = log_ratio(
            sent_row["drain_replayed_current_A_per_um"],
            sent_row["native_current_A_per_um"])
        vela_closure = log_ratio(
            vela_row["drain_replayed_current_A_per_um"],
            vela_row["native_current_A_per_um"])
        closure = native_gap - (state_contribution + sent_operator - vela_closure)
        decomposition_rows.append({
            "cell": cell,
            "native_gap_dex": native_gap,
            "common_operator_state_contribution_dex": state_contribution,
            "sentaurus_state_vela_operator_contribution_dex": sent_operator,
            "vela_cut_closure_dex": vela_closure,
            "reconstructed_gap_dex": state_contribution + sent_operator - vela_closure,
            "decomposition_closure_dex": closure,
            "vela_drain_condition": vela_row["drain_total_condition"],
            "sentaurus_state_drain_condition": sent_row["drain_total_condition"],
            "vela_source_condition": vela_row["source_total_condition"],
            "sentaurus_state_source_condition": sent_row["source_total_condition"],
        })
    write_csv(PORTABLE / "m30_native_gap_decomposition.csv", decomposition_rows)

    interaction_rows = []
    interaction_summary = []
    solver_fields: dict[str, dict[str, dict[str, dict[int, float]]]] = {}
    for solver in ("vela", "sentaurus"):
        solver_fields[solver] = {
            cell: state_fields(REPO / state_map[(solver, cell)]["state_file"])
            for cell in CELLS}
    sent_manifest = read_json(M29 / "sentaurus_export_manifest.json")
    sent_exports = {row["cell"]: REPO / row["export_dir"]
                    for row in sent_manifest["states"]}
    silicon = set(m10.scalar_field(sent_exports[CELLS[0]], "eDensity"))
    common = sorted(silicon & set.intersection(*(
        set(solver_fields[solver][cell]["psi"])
        for solver in ("vela", "sentaurus") for cell in CELLS)))
    for field in ("psi", "phin", "phip", "logn", "logp"):
        sent_values, vela_values = [], []
        for node in common:
            sent_value = interaction(solver_fields["sentaurus"], field, node)
            vela_value = interaction(solver_fields["vela"], field, node)
            sent_values.append(sent_value)
            vela_values.append(vela_value)
            interaction_rows.append({
                "field": field, "node_id": node,
                "sentaurus_interaction": sent_value,
                "vela_interaction": vela_value,
                "vela_minus_sentaurus": vela_value - sent_value,
            })
        differences = [abs(a - b) for a, b in zip(sent_values, vela_values)]
        interaction_summary.append({
            "field": field, "node_count": len(common),
            "pearson": pearson(sent_values, vela_values),
            "median_absolute_difference": percentile(differences, 0.5),
            "p95_absolute_difference": percentile(differences, 0.95),
            "maximum_absolute_difference": max(differences),
        })
    write_csv(PORTABLE / "m30_nodal_interaction.csv", interaction_rows)
    write_csv(PORTABLE / "m30_nodal_interaction_summary.csv", interaction_summary)

    double_off = next(row for row in decomposition_rows
                      if row["cell"] == "bgn_off_srh_off")
    control_rows = []
    for item in controls["controls"]:
        row = dict(item)
        if item["sg"]:
            edges = read_csv(REPO / item["sg"])
            source = cut_metrics(edges, source_nodes)
            drain = cut_metrics(edges, drain_nodes)
            substrate = cut_metrics(edges, substrate_nodes)
            gate = cut_metrics(edges, gate_nodes)
            all_total = sum(part["total_A_per_um"]
                            for part in (source, drain, substrate, gate))
            row.update({
                "source_current_A_per_um": source["total_A_per_um"],
                "drain_current_A_per_um": drain["total_A_per_um"],
                "substrate_current_A_per_um": substrate["total_A_per_um"],
                "gate_current_A_per_um": gate["total_A_per_um"],
                "all_terminal_imbalance_A_per_um": all_total,
                "drain_condition": drain["total_condition"],
                "source_condition": source["total_condition"],
            })
        control_rows.append(row)
    write_csv(PORTABLE / "m30_solver_control_summary.csv", control_rows)
    numeric_values = [float(value) for row in pair_rows for key, value in row.items()
                      if key not in ("pair_id", "state_solver", "state_cell",
                                     "operator_cell", "is_diagonal")]
    acceptance = {
        "state_count": len(states["states"]),
        "state_operator_pair_count": len(pair_rows),
        "probe_run_count": probes["probe_run_count"],
        "solver_control_count": len(control_rows),
        "common_silicon_node_count": len(common),
        "maximum_decomposition_closure_dex": max(
            abs(float(row["decomposition_closure_dex"]))
            for row in decomposition_rows),
        "nonfinite_fraction": sum(not math.isfinite(value)
                                  for value in numeric_values) / len(numeric_values),
    }
    expected = contract["acceptance"]
    checks = {
        "state_count": acceptance["state_count"] == expected["state_count"],
        "pair_count": acceptance["state_operator_pair_count"]
            == expected["state_operator_pair_count"],
        "probe_count": acceptance["probe_run_count"]
            == expected["probe_run_count"],
        "solver_control_count": acceptance["solver_control_count"]
            == expected["solver_control_count"],
        "common_nodes": acceptance["common_silicon_node_count"]
            >= expected["minimum_common_silicon_nodes"],
        "decomposition_closure": acceptance["maximum_decomposition_closure_dex"]
            <= expected["maximum_decomposition_closure_dex"],
        "finite": acceptance["nonfinite_fraction"] == 0.0,
    }
    acceptance["checks"] = checks
    acceptance["all_checks_pass"] = all(checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m30_double_off_causal_closure_report.v1",
        "status": "complete",
        "execution": {
            "new_sentaurus_execution": False,
            "cpp_changed": False,
            "default_model_changed": False,
        },
        "acceptance": acceptance,
        "double_off_decomposition": double_off,
        "all_cell_decomposition": decomposition_rows,
        "nodal_interaction_summary": interaction_summary,
        "solver_controls": control_rows,
        "artifacts": {
            "state_operator_matrix": portable(PORTABLE / "m30_state_operator_matrix.csv"),
            "native_gap_decomposition": portable(PORTABLE / "m30_native_gap_decomposition.csv"),
            "nodal_interaction": portable(PORTABLE / "m30_nodal_interaction.csv"),
            "nodal_interaction_summary": portable(
                PORTABLE / "m30_nodal_interaction_summary.csv"),
            "solver_control_summary": portable(
                PORTABLE / "m30_solver_control_summary.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    write_json(PORTABLE / "m30_double_off_causal_closure_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--run-controls", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    states = prepare(contract, args.force) if args.prepare else read_json(
        OUTPUT / "state_manifest.json")
    probes = run_probes(contract, states, args.runner.resolve(), args.workers,
                        args.force) if args.run else read_json(
                            OUTPUT / "probe_manifest.json")
    controls = run_solver_controls(
        contract, args.runner.resolve(), args.force) if args.run_controls else read_json(
            OUTPUT / "solver_control_manifest.json")
    if args.analyze:
        report = analyze(contract, states, probes, controls)
        print(json.dumps({
            "status": report["status"],
            "all_checks_pass": report["acceptance"]["all_checks_pass"],
            "double_off": report["double_off_decomposition"],
        }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
