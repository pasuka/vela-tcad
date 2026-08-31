#!/usr/bin/env python3
"""Audit the two-endpoint G3 gm reclose and localized electron-QF Jacobian."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any


def runner_environment() -> dict[str, str]:
    env = dict(os.environ)
    if os.name == "nt":
        env["PATH"] = os.pathsep.join(
            [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin", env.get("PATH", "")]
        )
    env["VELA_LINEAR_SOLVER"] = "sparselu"
    return env


def run_probe(runner: Path, config: dict[str, Any], path: Path) -> dict[str, Any]:
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(path.resolve())],
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (path.parent / f"{path.stem}.stdout.txt").write_text(
        completed.stdout, encoding="utf-8"
    )
    (path.parent / f"{path.stem}.stderr.txt").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{config['simulation_type']} failed: "
            f"{completed.stderr or completed.stdout}"
        )
    return json.loads(completed.stdout)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def read_state(path: Path) -> dict[int, dict[str, float]]:
    result: dict[int, dict[str, float]] = {}
    for row in read_rows(path):
        result[int(row["node_id"])] = {
            name: float(row[name]) for name in ("psi", "phin", "phip")
        }
    return result


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def vector_comparison(actual: list[float], target: list[float]) -> dict[str, float]:
    actual_norm = l2(actual)
    target_norm = l2(target)
    difference_norm = l2([a - b for a, b in zip(actual, target, strict=True)])
    dot = sum(a * b for a, b in zip(actual, target, strict=True))
    return {
        "actual_l2_V": actual_norm,
        "target_l2_V": target_norm,
        "actual_over_target_l2": actual_norm / target_norm if target_norm else 0.0,
        "relative_difference_l2": difference_norm / target_norm if target_norm else 0.0,
        "cosine_similarity": (
            dot / (actual_norm * target_norm) if actual_norm and target_norm else 0.0
        ),
    }


def focus_nodes(feedback_documents: list[dict[str, Any]], limit: int = 12) -> list[int]:
    nodes: set[int] = set()
    for document in feedback_documents:
        edges = document["points"][0]["edge_feedback"]["top_flux_feedback_edges"]
        for edge in edges[:limit]:
            nodes.update((int(edge["node0"]), int(edge["node1"])))
    return sorted(nodes)


def contact_nodes(mesh: dict[str, Any]) -> set[int]:
    return {
        int(node)
        for contact in mesh["contacts"]
        for node in contact.get("node_ids", [])
    }


def probe_config(
    baseline: dict[str, Any], state: Path, bias: float, output: Path,
    simulation_type: str,
) -> dict[str, Any]:
    config = deepcopy(baseline)
    config.pop("sweep", None)
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output.resolve())
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    return config


def endpoint_audit(
    runner: Path,
    baseline: dict[str, Any],
    bias: float,
    feedback: dict[str, Any],
    sentaurus_state_path: Path,
    reclosed_state_path: Path,
    reclose: dict[str, Any],
    focus: list[int],
    free: list[int],
    output: Path,
) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    step_csv = output / "newton_step.csv"
    step_config = probe_config(
        baseline, sentaurus_state_path, bias, step_csv, "newton_step_probe"
    )
    step_status = run_probe(runner, step_config, output / "newton_step.json")

    block_step_csv = output / "newton_block_step.csv"
    block_step_config = probe_config(
        baseline, sentaurus_state_path, bias, block_step_csv,
        "newton_block_step_probe",
    )
    block_step_config["block_modes"] = ["carrier_only"]
    block_step_status = run_probe(
        runner, block_step_config, output / "newton_block_step.json"
    )

    terms_csv = output / "carrier_terms.csv"
    terms_config = probe_config(
        baseline, sentaurus_state_path, bias, terms_csv, "newton_carrier_term_probe"
    )
    terms_status = run_probe(runner, terms_config, output / "carrier_terms.json")

    jvp_csv = output / "jvp.csv"
    jvp_config = probe_config(
        baseline, sentaurus_state_path, bias, jvp_csv, "newton_jvp_probe"
    )
    jvp_config["directions"] = [
        {
            "name": "phin_flux_feedback_nodes",
            "mode": "phin",
            "node_ids": focus,
            "exclude_contacts": False,
            "amplitude_V": 1.0e-6,
        },
        {
            "name": "phin_flux_feedback_patch_ring1",
            "mode": "phin",
            "node_ids": focus,
            "adjacent_cell_rings": 1,
            "exclude_contacts": False,
            "amplitude_V": 1.0e-6,
        },
    ] + [
        {
            "name": f"phin_node_{node}",
            "mode": "phin",
            "node_ids": [node],
            "exclude_contacts": False,
            "amplitude_V": 1.0e-6,
        }
        for node in focus[:4]
    ]
    jvp_status = run_probe(runner, jvp_config, output / "jvp.json")

    sentaurus = read_state(sentaurus_state_path)
    reclosed = read_state(reclosed_state_path)
    step_rows = {int(row["node_id"]): row for row in read_rows(step_csv)}
    block_step_rows = {
        int(row["node_id"]): row
        for row in read_rows(block_step_csv)
        if row["mode"] == "carrier_only"
    }
    term_rows = {int(row["node_id"]): row for row in read_rows(terms_csv)}
    total_delta = [reclosed[node]["phin"] - sentaurus[node]["phin"] for node in free]
    first_delta = [float(step_rows[node]["delta_phin_V"]) for node in free]
    carrier_only_delta = [
        float(block_step_rows[node]["delta_phin_V"]) for node in free
    ]
    focus_free = [node for node in focus if node in set(free)]
    focus_total = [reclosed[node]["phin"] - sentaurus[node]["phin"] for node in focus_free]
    focus_first = [float(step_rows[node]["delta_phin_V"]) for node in focus_free]
    all_residual = [float(term_rows[node]["electron_residual"]) for node in free]
    focus_residual = [float(term_rows[node]["electron_residual"]) for node in focus_free]
    all_total_norm = l2(total_delta)
    all_first_norm = l2(first_delta)
    all_residual_norm = l2(all_residual)

    top_update_nodes = sorted(
        (
            {
                "node_id": node,
                "x_um": float(step_rows[node]["x"]),
                "y_um": float(step_rows[node]["y"]),
                "total_delta_phin_V": reclosed[node]["phin"] - sentaurus[node]["phin"],
                "first_step_delta_phin_V": float(step_rows[node]["delta_phin_V"]),
                "electron_residual": float(term_rows[node]["electron_residual"]),
                "in_flux_feedback_focus": node in set(focus),
            }
            for node in free
        ),
        key=lambda row: abs(row["total_delta_phin_V"]),
        reverse=True,
    )[:20]
    top_residual_nodes = sorted(
        (
            {
                "node_id": node,
                "x_um": float(step_rows[node]["x"]),
                "y_um": float(step_rows[node]["y"]),
                "electron_residual": float(term_rows[node]["electron_residual"]),
                "total_delta_phin_V": reclosed[node]["phin"] - sentaurus[node]["phin"],
                "in_flux_feedback_focus": node in set(focus),
            }
            for node in free
        ),
        key=lambda row: abs(row["electron_residual"]),
        reverse=True,
    )[:20]
    jvp_rows = read_rows(jvp_csv)
    point = feedback["points"][0]
    return {
        "bias_V": bias,
        "reclose": reclose,
        "terminal_currents_A_per_um": {
            "sentaurus": point["sentaurus_terminal_A_per_um"],
            "frozen_operator_on_sentaurus_state": point["variants"]["SSS"]["current_A_per_um"],
            "self_consistent_vela": reclose["drain_current_A_per_um"],
        },
        "first_newton_step_vs_full_reclose_phin": vector_comparison(
            first_delta, total_delta
        ),
        "carrier_only_step_vs_full_reclose_phin": vector_comparison(
            carrier_only_delta, total_delta
        ),
        "focus_patch": {
            "direct_node_count": len(focus),
            "free_node_count": len(focus_free),
            "full_reclose_phin_l2_fraction": l2(focus_total) / all_total_norm if all_total_norm else 0.0,
            "first_step_phin_l2_fraction": l2(focus_first) / all_first_norm if all_first_norm else 0.0,
            "electron_residual_l2_fraction": l2(focus_residual) / all_residual_norm if all_residual_norm else 0.0,
        },
        "top_full_reclose_phin_update_nodes": top_update_nodes,
        "top_initial_electron_residual_nodes": top_residual_nodes,
        "jacobian_directional_derivative": {
            "max_relative_error": max(float(row["relative_error"]) for row in jvp_rows),
            "directions": jvp_rows,
        },
        "probe_status": {
            "newton_step": step_status,
            "newton_block_step": block_step_status,
            "carrier_terms": terms_status,
            "jvp": jvp_status,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--bias", type=float, action="append", required=True)
    parser.add_argument("--feedback-summary", type=Path, action="append", required=True)
    parser.add_argument("--sentaurus-state", type=Path, action="append", required=True)
    parser.add_argument("--reclosed-state", type=Path, action="append", required=True)
    parser.add_argument("--reclose-summary", type=Path, action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    counts = {
        len(args.bias), len(args.feedback_summary), len(args.sentaurus_state),
        len(args.reclosed_state), len(args.reclose_summary),
    }
    if len(counts) != 1 or len(args.bias) != 2:
        raise ValueError("exactly two aligned endpoint inputs are required")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    feedbacks = [json.loads(path.read_text(encoding="utf-8")) for path in args.feedback_summary]
    recloses = [json.loads(path.read_text(encoding="utf-8")) for path in args.reclose_summary]
    mesh = json.loads(Path(baseline["mesh_file"]).read_text(encoding="utf-8"))
    contacts = contact_nodes(mesh)
    free = [node for node in range(len(mesh["nodes"])) if node not in contacts]
    focus = focus_nodes(feedbacks)

    endpoints = []
    for index, bias in enumerate(args.bias):
        endpoints.append(endpoint_audit(
            args.runner.resolve(), baseline, bias, feedbacks[index],
            args.sentaurus_state[index], args.reclosed_state[index], recloses[index],
            focus, free, output / f"vg_{bias:.6f}".replace(".", "p"),
        ))

    delta_bias = args.bias[1] - args.bias[0]
    currents = [endpoint["terminal_currents_A_per_um"] for endpoint in endpoints]
    slopes = {
        name: (currents[1][name] - currents[0][name]) / delta_bias
        for name in currents[0]
    }
    sentaurus_gm = slopes["sentaurus"]
    summary = {
        "schema": "vela.templates_ldmos.g3_gm_reclose_jacobian.v1",
        "contract": {
            "operator": "qualified_contact_hfs_external_averagebox",
            "contact_boundary_reconstruction": "legacy_node_local",
            "initial_state": "sentaurus_endpoint_state",
            "reclose": "same_bias_full_newton",
            "jacobian_check": "double_symmetric_finite_difference",
            "production_defaults_changed": False,
        },
        "focus_nodes": focus,
        "endpoints": endpoints,
        "gm_A_per_um_V": slopes,
        "gm_ratio_to_sentaurus": {
            name: value / sentaurus_gm for name, value in slopes.items()
        },
        "classification": (
            "complete_sentaurus_state_requires_coupled_reclose"
            if abs(slopes["frozen_operator_on_sentaurus_state"] / sentaurus_gm - 1.0) < 0.01
            and abs(slopes["self_consistent_vela"] / sentaurus_gm - 1.0) > 0.2
            else "not_closed"
        ),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "classification": summary["classification"],
        "focus_node_count": len(focus),
        "gm_ratio_to_sentaurus": summary["gm_ratio_to_sentaurus"],
        "reclose_iterations": [endpoint["reclose"]["iterations"] for endpoint in endpoints],
        "max_jvp_relative_error": max(
            endpoint["jacobian_directional_derivative"]["max_relative_error"]
            for endpoint in endpoints
        ),
        "summary": str((output / "summary.json").resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
