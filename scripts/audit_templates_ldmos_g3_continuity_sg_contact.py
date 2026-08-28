#!/usr/bin/env python3
"""Audit finite-bias SG cancellation and the drain-contact branch for LDMOS G3.

The audit replays immutable Vela and Sentaurus states through production Vela
operators.  It checks three identities independently:

* SG edge divergence versus the Newton carrier flux term;
* SG drain-cut integration versus ContactCurrent;
* stable/long-double contact current versus the diagnostic drift+diffusion split.

Generated decks and large CSV files belong in an ignored staging directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics
import subprocess
from collections import defaultdict, deque
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
Q_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sentaurus_plt_current(path: Path, bias: float) -> dict[str, float | int | str]:
    """Read the highest-time DF-ISE row at the requested drain voltage."""
    text = path.read_text(encoding="utf-8", errors="replace")
    datasets_match = re.search(r"datasets\s*=\s*\[(.*?)\]", text, re.DOTALL)
    data_match = re.search(r"Data\s*\{(.*)\}\s*$", text, re.DOTALL)
    if datasets_match is None or data_match is None:
        raise ValueError(f"Not a DF-ISE xyplot file: {path}")
    datasets = re.findall(r'"([^"]+)"', datasets_match.group(1))
    values = [
        float(token)
        for token in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?",
            data_match.group(1),
        )
    ]
    if not datasets or len(values) % len(datasets) != 0:
        raise ValueError(
            f"DF-ISE data length {len(values)} is not divisible by {len(datasets)} datasets"
        )
    rows = [
        dict(zip(datasets, values[index : index + len(datasets)]))
        for index in range(0, len(values), len(datasets))
    ]
    closest_error = min(abs(row["drain OuterVoltage"] - bias) for row in rows)
    candidates = [
        row
        for row in rows
        if abs(row["drain OuterVoltage"] - bias) <= closest_error + 1.0e-15
    ]
    selected = max(candidates, key=lambda row: row["time"])
    return {
        "source": str(path.resolve()),
        "requested_bias_V": bias,
        "matched_bias_V": selected["drain OuterVoltage"],
        "time": selected["time"],
        "drain_electron_A_per_um": selected["drain eCurrent"],
        "drain_hole_A_per_um": selected["drain hCurrent"],
        "drain_total_A_per_um": selected["drain TotalCurrent"],
        "matching_rows": len(candidates),
    }


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def finite_distribution(values: Iterable[float]) -> dict[str, float | int]:
    samples = sorted(value for value in values if math.isfinite(value))
    if not samples:
        return {"count": 0, "median": math.nan, "p95": math.nan, "maximum": math.nan}
    p95_index = min(len(samples) - 1, math.ceil(0.95 * len(samples)) - 1)
    return {
        "count": len(samples),
        "median": statistics.median(samples),
        "p95": samples[p95_index],
        "maximum": samples[-1],
    }


def safe_condition(term0: float, term1: float, result: float) -> float:
    scale = abs(term0) + abs(term1)
    if scale == 0.0:
        return 0.0
    if result == 0.0:
        return math.inf
    return scale / abs(result)


def relative_error(candidate: float, reference: float, floor: float = 1.0e-300) -> float:
    return abs(candidate - reference) / max(abs(reference), floor)


def runner_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["VELA_LINEAR_SOLVER"] = "sparselu"
    return environment


def run_config(
    runner: Path,
    config: dict[str, Any],
    run_dir: Path,
    stem: str,
    allow_failure: bool = False,
) -> dict[str, Any]:
    run_dir.mkdir(parents=True, exist_ok=True)
    config_path = run_dir / f"{stem}.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [
            str(runner),
            "--config",
            str(config_path.resolve()),
            "--log",
            str((run_dir / f"{stem}.log").resolve()),
        ],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (run_dir / f"{stem}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_dir / f"{stem}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode and not allow_failure:
        raise RuntimeError(f"{stem} failed: {completed.stderr or completed.stdout}")
    status = (
        json.loads(completed.stdout.strip().splitlines()[-1])
        if completed.stdout.strip()
        else {}
    )
    status["return_code"] = completed.returncode
    return status


def probe_config(
    baseline: dict[str, Any], state: Path, output_csv: Path, kind: str, bias: float
) -> dict[str, Any]:
    config = deepcopy(baseline)
    config["simulation_type"] = kind
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output_csv.resolve())
    config.pop("sweep", None)
    config.pop("log_file", None)
    for contact in config["contacts"]:
        contact["bias"] = bias if contact["name"].lower() == "drain" else 0.0
    config["_comment"] = "Immutable finite-bias state replay for SG/contact audit."
    return config


def frozen_config(
    baseline: dict[str, Any], state: Path, output_dir: Path, bias: float
) -> dict[str, Any]:
    config = deepcopy(baseline)
    contacts = [contact["name"] for contact in config["contacts"]]
    config["output_csv"] = str((output_dir / "curve.csv").resolve())
    config["solver"]["method"] = "frozen_state"
    config["solver"].pop("local_update_diagnostics", None)
    config["sweep"] = {
        "mode": "iv",
        "contact": "drain",
        "current_contact": "drain",
        "start": bias,
        "stop": bias,
        "step": 1.0,
        "bias_points": [bias],
        "initial_state_file": str(state.resolve()),
        "write_state_file": str((output_dir / "replayed_state.csv").resolve()),
        "frozen_state_compute_current": True,
        "write_vtk": False,
        "diagnostics": {
            "terminal_balance": {
                "enabled": True,
                "contacts": contacts,
                "csv_file": str((output_dir / "terminal_balance.csv").resolve()),
            },
            "contact_edge": {
                "enabled": True,
                "contacts": contacts,
                "csv_file": str((output_dir / "contact_edges.csv").resolve()),
            },
            "continuity_balance": {
                "enabled": True,
                "contacts": contacts,
                "csv_file": str((output_dir / "continuity_balance.csv").resolve()),
            },
            "transport": {"enabled": True},
        },
    }
    config["_comment"] = "Frozen finite-bias state replay; physics and discretization unchanged."
    return config


def strict_reclose_config(
    baseline: dict[str, Any], state: Path, output_dir: Path, bias: float
) -> dict[str, Any]:
    """Single-factor qualification: tighten convergence, preserve all physics."""
    config = deepcopy(baseline)
    config["output_csv"] = str((output_dir / "curve.csv").resolve())
    config["solver"]["method"] = "newton"
    config["solver"]["reltol"] = 1.0e-10
    config["solver"]["abstol"] = 1.0e-14
    config["solver"]["stall_residual_floor"] = 1.0e-12
    config["sweep"] = {
        "mode": "iv",
        "contact": "drain",
        "current_contact": "drain",
        "start": bias,
        "stop": bias,
        "step": 1.0,
        "bias_points": [bias],
        "initial_state_file": str(state.resolve()),
        "write_state_file": str((output_dir / "state.csv").resolve()),
        "write_vtk": False,
        "diagnostics": {
            "terminal_balance": {
                "enabled": True,
                "contacts": [contact["name"] for contact in config["contacts"]],
                "csv_file": str((output_dir / "terminal_balance.csv").resolve()),
            },
            "contact_edge": {
                "enabled": True,
                "contacts": [contact["name"] for contact in config["contacts"]],
                "csv_file": str((output_dir / "contact_edges.csv").resolve()),
            },
            "newton_history": {
                "enabled": True,
                "csv_file": str((output_dir / "newton_history.csv").resolve()),
                "attempts_csv_file": str((output_dir / "newton_attempts.csv").resolve()),
                "iterations_csv_file": str((output_dir / "newton_iterations.csv").resolve()),
                "rejected_state_directory": str((output_dir / "rejected_states").resolve()),
            },
            "transport": {"enabled": True},
        },
    }
    config["_comment"] = (
        "G3 single-factor strict reclose: unchanged revision-4 physics, "
        "reltol=1e-10, abstol=1e-14, no predictor."
    )
    return config


def mesh_topology(mesh_path: Path) -> dict[str, Any]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = {int(node["id"]): (float(node["x"]), float(node["y"])) for node in mesh["nodes"]}
    contacts = {
        contact["name"]: {int(node) for node in contact.get("node_ids", [])}
        for contact in mesh["contacts"]
    }
    return {"nodes": nodes, "contacts": contacts}


def graph_rings(
    edge_rows: list[dict[str, str]], seeds: set[int], max_ring: int = 3
) -> dict[int, int]:
    graph: dict[int, set[int]] = defaultdict(set)
    for row in edge_rows:
        node0, node1 = int(row["node0"]), int(row["node1"])
        graph[node0].add(node1)
        graph[node1].add(node0)
    rings = {node: 0 for node in seeds}
    queue = deque(seeds)
    while queue:
        node = queue.popleft()
        if rings[node] >= max_ring:
            continue
        for neighbor in graph[node]:
            if neighbor not in rings:
                rings[neighbor] = rings[node] + 1
                queue.append(neighbor)
    return rings


def reconstruct_divergence(
    edge_rows: list[dict[str, str]], carrier: str
) -> dict[int, float]:
    divergence: dict[int, float] = defaultdict(float)
    column = f"{carrier}_flux"
    for row in edge_rows:
        value = f(row, column)
        divergence[int(row["node0"])] += value
        divergence[int(row["node1"])] -= value
    return dict(divergence)


def top_term_rows(
    rows: list[dict[str, str]], rings: dict[int, int], column: str, limit: int = 12
) -> list[dict[str, Any]]:
    ranked = sorted(rows, key=lambda row: abs(f(row, column)), reverse=True)[:limit]
    return [
        {
            "node_id": int(row["node_id"]),
            "x_um": f(row, "x"),
            "y_um": f(row, "y"),
            "drain_ring": rings.get(int(row["node_id"])),
            column: f(row, column),
            "electron_flux_abs_sum": f(row, "electron_flux_abs_sum"),
            "flux_cancellation_condition": safe_condition(
                f(row, "electron_flux_abs_sum"), 0.0, f(row, "electron_flux")
            ),
        }
        for row in ranked
    ]


def contact_summary(
    rows: list[dict[str, str]], contact: str, contact_nodes: set[int]
) -> dict[str, Any]:
    selected = [row for row in rows if row["current_contact"] == contact]
    electron = sum(f(row, "current_electron") for row in selected)
    electron_long = sum(f(row, "current_electron_long_double_reference") for row in selected)
    electron_drift = sum(f(row, "current_electron_drift") for row in selected)
    electron_diffusion = sum(f(row, "current_electron_diffusion") for row in selected)
    hole = sum(f(row, "current_hole") for row in selected)
    total = sum(f(row, "current_total") for row in selected)
    edge_records = []
    for row in selected:
        node0, node1 = int(row["node0"]), int(row["node1"])
        interior = node1 if node0 in contact_nodes else node0
        stable = f(row, "current_electron")
        drift = f(row, "current_electron_drift")
        diffusion = f(row, "current_electron_diffusion")
        edge_records.append(
            {
                "edge_id": int(row["edge_id"]),
                "contact_node": node0 if node0 in contact_nodes else node1,
                "interior_node": interior,
                "stable_electron_A_per_um": stable / 1.0e6,
                "long_double_electron_A_per_um": f(
                    row, "current_electron_long_double_reference"
                ) / 1.0e6,
                "drift_A_per_um": drift / 1.0e6,
                "diffusion_A_per_um": diffusion / 1.0e6,
                "drift_diffusion_sum_A_per_um": (drift + diffusion) / 1.0e6,
                "cancellation_condition": safe_condition(drift, diffusion, stable),
                "stable_vs_long_double_relative_error": relative_error(
                    stable, f(row, "current_electron_long_double_reference")
                ),
                "electron_qf_drop_V": abs(f(row, "phin1") - f(row, "phin0")),
                "psi_drop_V": abs(f(row, "psi1") - f(row, "psi0")),
                "electron_density_interior_m3": f(
                    row, "n1" if node0 in contact_nodes else "n0"
                ),
            }
        )
    top_stable = sorted(
        edge_records, key=lambda row: abs(row["stable_electron_A_per_um"]), reverse=True
    )[:12]
    top_cancellation = sorted(
        edge_records, key=lambda row: row["cancellation_condition"], reverse=True
    )[:12]
    return {
        "edge_count": len(selected),
        "electron_A_per_um": electron / 1.0e6,
        "electron_long_double_A_per_um": electron_long / 1.0e6,
        "electron_drift_A_per_um": electron_drift / 1.0e6,
        "electron_diffusion_A_per_um": electron_diffusion / 1.0e6,
        "electron_drift_diffusion_sum_A_per_um": (electron_drift + electron_diffusion) / 1.0e6,
        "hole_A_per_um": hole / 1.0e6,
        "total_A_per_um": total / 1.0e6,
        "stable_vs_long_double_relative_error": relative_error(electron, electron_long),
        "drift_diffusion_cancellation_condition": safe_condition(
            electron_drift, electron_diffusion, electron
        ),
        "edge_cancellation_condition": finite_distribution(
            record["cancellation_condition"] for record in edge_records
        ),
        "electron_qf_drop_V": finite_distribution(
            record["electron_qf_drop_V"] for record in edge_records
        ),
        "top_edges_by_stable_electron_current": top_stable,
        "top_edges_by_cancellation_condition": top_cancellation,
    }


def sg_contact_cut(
    rows: list[dict[str, str]], contact_nodes: set[int]
) -> dict[str, float | int]:
    electron = 0.0
    hole = 0.0
    count = 0
    for row in rows:
        node0, node1 = int(row["node0"]), int(row["node1"])
        n0_contact, n1_contact = node0 in contact_nodes, node1 in contact_nodes
        if n0_contact == n1_contact:
            continue
        outward_sign = 1.0 if n0_contact else -1.0
        electron += -Q_C * outward_sign * f(row, "electron_particle_line_flux_per_m_s")
        hole += -Q_C * outward_sign * f(row, "hole_particle_line_flux_per_m_s")
        count += 1
    return {
        "crossing_edge_count": count,
        "electron_A_per_um": electron / 1.0e6,
        "hole_A_per_um": hole / 1.0e6,
        "total_A_per_um": (electron - hole) / 1.0e6,
    }


def state_summary(
    label: str,
    run_dir: Path,
    topology: dict[str, Any],
) -> dict[str, Any]:
    sg_rows = read_csv(run_dir / "sg_edge_flux.csv")
    term_rows = read_csv(run_dir / "newton_carrier_term.csv")
    contact_rows = read_csv(run_dir / "contact_edges.csv")
    curve = read_csv(run_dir / "curve.csv")[-1]
    term_by_node = {int(row["node_id"]): row for row in term_rows}
    all_contact_nodes = set().union(*topology["contacts"].values())
    drain_nodes = topology["contacts"]["drain"]
    rings = graph_rings(sg_rows, drain_nodes)
    electron_divergence = reconstruct_divergence(sg_rows, "electron")
    free_nodes = sorted(set(term_by_node) - all_contact_nodes)
    divergence_errors = [
        abs(electron_divergence[node] - f(term_by_node[node], "electron_flux"))
        for node in free_nodes
    ]
    drain_contact = contact_summary(contact_rows, "drain", drain_nodes)
    drain_cut = sg_contact_cut(sg_rows, drain_nodes)
    sg_by_edge = {int(row["edge_id"]): row for row in sg_rows}
    density_errors_dex = []
    for row in contact_rows:
        if row["current_contact"] != "drain":
            continue
        sg = sg_by_edge[int(row["edge_id"])]
        for endpoint in ("0", "1"):
            raw_m3 = f(row, f"n{endpoint}") * 1.0e6
            reconstructed_m3 = f(sg, f"electron_density{endpoint}_m3")
            if raw_m3 <= 0.0 and reconstructed_m3 <= 0.0:
                continue
            density_errors_dex.append(
                abs(
                    math.log10(max(raw_m3, 1.0e-300))
                    - math.log10(max(reconstructed_m3, 1.0e-300))
                )
            )
    current_A_per_um = f(curve, "current_total_A_per_um")
    electron_residual_l2 = math.sqrt(
        sum(f(row, "electron_residual") ** 2 for row in term_rows)
    )
    ring_metrics = {}
    for ring in range(4):
        selected = [row for row in term_rows if rings.get(int(row["node_id"])) == ring]
        ring_metrics[str(ring)] = {
            "node_count": len(selected),
            "electron_flux_l2": math.sqrt(sum(f(row, "electron_flux") ** 2 for row in selected)),
            "electron_residual_l2": math.sqrt(
                sum(f(row, "electron_residual") ** 2 for row in selected)
            ),
            "flux_cancellation_condition": finite_distribution(
                safe_condition(
                    f(row, "electron_flux_abs_sum"), 0.0, f(row, "electron_flux")
                )
                for row in selected
            ),
        }
    return {
        "label": label,
        "current_A_per_um": current_A_per_um,
        "electron_continuity_residual_l2": electron_residual_l2,
        "drain_contact": drain_contact,
        "sg_drain_cut": drain_cut,
        "drain_raw_vs_operator_electron_density_abs_error_dex": finite_distribution(
            density_errors_dex
        ),
        "closure": {
            "sg_cut_vs_contact_total_relative_error": relative_error(
                float(drain_cut["total_A_per_um"]),
                float(drain_contact["total_A_per_um"]),
            ),
            "contact_total_vs_curve_relative_error": relative_error(
                float(drain_contact["total_A_per_um"]), current_A_per_um
            ),
            "free_node_sg_divergence_vs_carrier_flux_max_abs": max(
                divergence_errors, default=0.0
            ),
            "free_node_sg_divergence_vs_carrier_flux_l2": math.sqrt(
                sum(error * error for error in divergence_errors)
            ),
        },
        "drain_rings": ring_metrics,
        "top_electron_residual_nodes": top_term_rows(
            term_rows, rings, "electron_residual"
        ),
        "top_electron_flux_nodes": top_term_rows(term_rows, rings, "electron_flux"),
        "artifacts": {
            "curve": str((run_dir / "curve.csv").resolve()),
            "contact_edges": str((run_dir / "contact_edges.csv").resolve()),
            "continuity_balance": str((run_dir / "continuity_balance.csv").resolve()),
            "sg_edge_flux": str((run_dir / "sg_edge_flux.csv").resolve()),
            "newton_carrier_term": str((run_dir / "newton_carrier_term.csv").resolve()),
        },
    }


def classify(vela: dict[str, Any], sentaurus: dict[str, Any]) -> dict[str, Any]:
    vela_current = abs(float(vela["current_A_per_um"]))
    sentaurus_current = abs(float(sentaurus["current_A_per_um"]))
    current_gap_dex = abs(
        math.log10(max(vela_current, 1.0e-300))
        - math.log10(max(sentaurus_current, 1.0e-300))
    )
    vela_operator_closed = all(
        vela["closure"][key] <= 1.0e-9
        for key in (
            "sg_cut_vs_contact_total_relative_error",
            "contact_total_vs_curve_relative_error",
        )
    ) and vela["closure"]["free_node_sg_divergence_vs_carrier_flux_max_abs"] <= 1.0e-18
    divergence_closed = all(
        state["closure"]["free_node_sg_divergence_vs_carrier_flux_max_abs"] <= 1.0e-18
        for state in (vela, sentaurus)
    )
    sentaurus_external_state_qualified = (
        sentaurus["closure"]["sg_cut_vs_contact_total_relative_error"] <= 1.0e-3
        and sentaurus["closure"]["contact_total_vs_curve_relative_error"] <= 1.0e-9
        and sentaurus["drain_raw_vs_operator_electron_density_abs_error_dex"]["maximum"]
        > 0.0
    )
    precision_closed = all(
        state["drain_contact"]["stable_vs_long_double_relative_error"] <= 1.0e-9
        for state in (vela, sentaurus)
    )
    return {
        "same_operator_current_gap_dex": current_gap_dex,
        "vela_production_operator_closure_pass": vela_operator_closed,
        "sg_divergence_closure_pass": divergence_closed,
        "sentaurus_external_state_replay_qualified": sentaurus_external_state_qualified,
        "sentaurus_contact_mismatch_explanation": (
            "ContactCurrent consumes supplied density while the Newton SG operator "
            "reconstructs density from psi/QF; the external Sentaurus state is not "
            "exactly self-consistent under the Vela mapping."
        ),
        "stable_vs_long_double_precision_pass": precision_closed,
        "drift_diffusion_split_is_diagnostic_only": any(
            state["drain_contact"]["drift_diffusion_cancellation_condition"] >= 1.0e8
            for state in (vela, sentaurus)
        ),
        "root_cause_branch": (
            "self_consistent_contact_neighborhood_state_branch"
            if vela_operator_closed
            and divergence_closed
            and sentaurus_external_state_qualified
            and precision_closed
            and current_gap_dex >= 1.0
            else "operator_or_precision_not_closed"
        ),
        "ialmob_authorized": False,
    }


def replay_state(
    runner: Path,
    baseline: dict[str, Any],
    state: Path,
    run_dir: Path,
    bias: float,
    topology: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    statuses = {
        "frozen": run_config(
            runner,
            frozen_config(baseline, state, run_dir, bias),
            run_dir,
            "frozen_replay",
        ),
        "sg": run_config(
            runner,
            probe_config(
                baseline, state, run_dir / "sg_edge_flux.csv", "sg_edge_flux_probe", bias
            ),
            run_dir,
            "sg_edge_flux",
        ),
        "carrier_terms": run_config(
            runner,
            probe_config(
                baseline,
                state,
                run_dir / "newton_carrier_term.csv",
                "newton_carrier_term_probe",
                bias,
            ),
            run_dir,
            "newton_carrier_term",
        ),
    }
    return statuses, state_summary(run_dir.name, run_dir, topology)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--vela-state", type=Path, required=True)
    parser.add_argument("--sentaurus-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, default=0.0231559774221138)
    parser.add_argument("--run-strict-reclose", action="store_true")
    parser.add_argument("--sentaurus-reference-plt", type=Path)
    args = parser.parse_args()

    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    topology = mesh_topology(Path(baseline["mesh_file"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    statuses = {}
    for label, state in (("vela", args.vela_state), ("sentaurus", args.sentaurus_state)):
        run_dir = args.output_dir / label
        run_dir.mkdir(parents=True, exist_ok=True)
        statuses[label], _ = replay_state(
            args.runner.resolve(), baseline, state, run_dir, args.bias, topology
        )

    vela = state_summary("vela", args.output_dir / "vela", topology)
    sentaurus = state_summary("sentaurus", args.output_dir / "sentaurus", topology)
    controls: dict[str, Any] = {}
    if args.run_strict_reclose:
        reclose_dir = args.output_dir / "vela_strict_reclose"
        reclose_status = run_config(
            args.runner.resolve(),
            strict_reclose_config(baseline, args.vela_state, reclose_dir, args.bias),
            reclose_dir,
            "strict_reclose",
            allow_failure=True,
        )
        controls["strict_reclose"] = {
            "status": reclose_status,
            "config": str((reclose_dir / "strict_reclose.json").resolve()),
            "curve": str((reclose_dir / "curve.csv").resolve()),
            "state": str((reclose_dir / "state.csv").resolve()),
        }
        attempts_path = reclose_dir / "newton_attempts.csv"
        if attempts_path.is_file():
            attempts = read_csv(attempts_path)
            rejected_final = Path(attempts[-1].get("rejected_final_state_file", ""))
            if rejected_final.is_file():
                final_dir = reclose_dir / "final_state_replay"
                final_statuses, _ = replay_state(
                    args.runner.resolve(),
                    baseline,
                    rejected_final,
                    final_dir,
                    args.bias,
                    topology,
                )
                final_summary = state_summary("strict_reclose_final", final_dir, topology)
                controls["strict_reclose"]["rejected_final_state"] = str(
                    rejected_final.resolve()
                )
                controls["strict_reclose"]["final_state_replay_statuses"] = final_statuses
                controls["strict_reclose"]["final_state_replay"] = final_summary

    sentaurus_reference = (
        sentaurus_plt_current(args.sentaurus_reference_plt, args.bias)
        if args.sentaurus_reference_plt is not None
        else None
    )
    if sentaurus_reference is not None and "strict_reclose" in controls:
        final_replay = controls["strict_reclose"].get("final_state_replay")
        if final_replay is not None:
            candidate = abs(float(final_replay["current_A_per_um"]))
            reference = abs(float(sentaurus_reference["drain_total_A_per_um"]))
            controls["strict_reclose"]["final_vs_sentaurus_plt"] = {
                "candidate_A_per_um": candidate,
                "reference_A_per_um": reference,
                "signed_ratio": candidate / reference,
                "relative_error": relative_error(candidate, reference),
                "magnitude_error_dex": abs(math.log10(candidate) - math.log10(reference)),
                "electron_residual_improvement": (
                    vela["electron_continuity_residual_l2"]
                    / final_replay["electron_continuity_residual_l2"]
                ),
            }

    report = {
        "schema": "vela.templates_ldmos_g3_continuity_sg_contact_audit.v1",
        "status": "complete",
        "bias_V": args.bias,
        "contracts": {
            "physics": "revision-4 Fermi+OldSlotboom+constant-field; immutable",
            "state": "same-bias frozen replay; no predictor and no nonlinear update",
            "sg": "production CoupledDDAssembler edge flux",
            "contact": "production ContactCurrent SG cut, 1 um depth",
            "drift_diffusion": "diagnostic algebraic split only",
        },
        "statuses": statuses,
        "vela_state": vela,
        "sentaurus_state": sentaurus,
        "sentaurus_plt_reference": sentaurus_reference,
        "controls": controls,
        "decision": classify(vela, sentaurus),
    }
    summary = args.output_dir / "summary.json"
    summary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": str(summary.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
