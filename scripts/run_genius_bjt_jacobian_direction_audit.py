#!/usr/bin/env python3
"""Audit selected Genius BJT Newton Jacobian rows against centered differences."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
FIXTURE = REPO / "reference_tcad" / "genius_bjt_sentaurus2022"
BUILD_ROOT = REPO / "build-release" / "reference_tcad" / "genius_bjt_sentaurus2022"
DEFAULT_OUTPUT = BUILD_ROOT / "jacobian_direction_audit"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
INTERNAL_ROWS = (4438, 4440)
FD_STEPS_V = (1.0e-5, 1.0e-6, 1.0e-7)
BLOCKS = ("psi", "phin", "phip")


def absolute(path: Path) -> str:
    return str(path.resolve())


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def one_ring(mesh: dict, seeds: tuple[int, ...]) -> list[int]:
    seed_set = set(seeds)
    selected = set(seeds)
    for triangle in mesh["triangles"]:
        nodes = triangle.get("node_ids", triangle.get("nodes"))
        if seed_set.intersection(nodes):
            selected.update(nodes)
    return sorted(selected)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summarize-only", action="store_true")
    return parser.parse_args()


def build_config(output: Path) -> tuple[dict, list[int], list[int]]:
    base_path = FIXTURE / "vela" / "configs" / "m1_spatial_vce3.json"
    mesh_path = FIXTURE / "vela" / "input" / "mesh.json"
    mesh = load_json(mesh_path)
    collector = next(c for c in mesh["contacts"] if c["name"] == "collector")
    collector_nodes = sorted(collector["node_ids"])
    local_nodes = one_ring(mesh, INTERNAL_ROWS)

    config = load_json(base_path)
    config.pop("sweep", None)
    config["simulation_type"] = "newton_jvp_probe"
    config["mesh_file"] = absolute(mesh_path)
    config["node_doping_file"] = absolute(FIXTURE / "vela" / "input" / "doping.csv")
    config["materials_file"] = absolute(FIXTURE / "vela" / "materials_sentaurus2022.json")
    config["state_file"] = absolute(BUILD_ROOT / "vela_wp3_wp5" / "m1_vce300_state.csv")
    config["output_csv"] = absolute(output / "direction_summary.csv")
    config["row_output_csv"] = absolute(output / "row_samples.csv")

    directions = []
    # One complete local-stencil basis at 1e-6 V identifies an offending column.
    for block in BLOCKS:
        for node in local_nodes:
            directions.append({
                "name": f"basis_{block}_{node}_h1e-6",
                "mode": block,
                "node_ids": [node],
                "exclude_contacts": False,
                "amplitude_V": 1.0e-6,
            })
    # Multi-step aggregate probes expose finite-difference truncation/noise trends.
    for step in FD_STEPS_V:
        token = f"{step:.0e}".replace("-", "m")
        for block in BLOCKS:
            directions.append({
                "name": f"local_{block}_h{token}",
                "mode": block,
                "node_ids": local_nodes,
                "exclude_contacts": False,
                "amplitude_V": step,
            })
            directions.append({
                "name": f"collector_{block}_h{token}",
                "mode": block,
                "node_ids": collector_nodes,
                "exclude_contacts": False,
                "amplitude_V": step,
            })
    config["directions"] = directions
    config["sample_rows"] = [
        {"name": f"internal_phin_{node}", "block": "phin", "node_id": node}
        for node in INTERNAL_ROWS
    ] + [
        {"name": f"collector_{block}_{node}", "block": block, "node_id": node}
        for block in BLOCKS for node in collector_nodes
    ]
    return config, local_nodes, collector_nodes


def run_probe(runner: Path, config_path: Path, output: Path) -> None:
    process = subprocess.run(
        [absolute(runner), "--config", absolute(config_path)],
        cwd=REPO,
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        check=False,
    )
    (output / "runner.stdout.log").write_text(process.stdout, encoding="utf-8")
    (output / "runner.stderr.log").write_text(process.stderr, encoding="utf-8")
    write_json(output / "runner.status.json", {"return_code": process.returncode})
    if process.returncode:
        raise RuntimeError(f"JVP probe failed; see {output / 'runner.stderr.log'}")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def summarize_internal(rows: list[dict[str, str]]) -> list[dict]:
    result = []
    for row_node in INTERNAL_ROWS:
        chosen = [
            row for row in rows
            if int(row["row_node"]) == row_node
            and row["row_block"] == "phin"
            and row["direction"].startswith("basis_")
        ]
        analytic = [float(row["analytic_derivative"]) for row in chosen]
        finite = [float(row["finite_difference_derivative"]) for row in chosen]
        errors = [a - f for a, f in zip(analytic, finite)]
        worst = max(chosen, key=lambda row: float(row["absolute_derivative_error"]))
        scale = max(norm(analytic), norm(finite), 1.0e-300)
        result.append({
            "row_node": row_node,
            "basis_columns": len(chosen),
            "analytic_row_l2": norm(analytic),
            "finite_difference_row_l2": norm(finite),
            "difference_l2": norm(errors),
            "relative_row_l2_error": norm(errors) / scale,
            "worst_direction": worst["direction"],
            "worst_absolute_derivative_error": float(worst["absolute_derivative_error"]),
            "worst_analytic_derivative": float(worst["analytic_derivative"]),
            "worst_finite_difference_derivative": float(worst["finite_difference_derivative"]),
        })
    return result


def summarize_steps(rows: list[dict[str, str]]) -> list[dict]:
    grouped: dict[tuple[str, float, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        direction = row["direction"]
        if not direction.startswith("local_"):
            continue
        grouped[(row["direction_mode"], float(row["amplitude_V"]), int(row["row_node"]))].append(row)
    result = []
    for (block, step, node), chosen in sorted(grouped.items()):
        if node not in INTERNAL_ROWS:
            continue
        row = chosen[0]
        analytic = float(row["analytic_derivative"])
        finite = float(row["finite_difference_derivative"])
        result.append({
            "row_node": node,
            "direction_block": block,
            "step_V": step,
            "analytic_derivative": analytic,
            "finite_difference_derivative": finite,
            "absolute_error": abs(analytic - finite),
            "relative_error": abs(analytic - finite) / max(abs(analytic), abs(finite), 1.0e-300),
        })
    return result


def summarize_contacts(rows: list[dict[str, str]], collector_nodes: list[int]) -> list[dict]:
    result = []
    for step in FD_STEPS_V:
        for direction_block in BLOCKS:
            prefix = f"collector_{direction_block}_"
            chosen = [
                row for row in rows
                if row["direction"].startswith(prefix)
                and math.isclose(float(row["amplitude_V"]), step)
                and int(row["row_node"]) in collector_nodes
            ]
            analytic_fd = max(float(row["absolute_derivative_error"]) for row in chosen)
            own = [
                float(row["analytic_derivative"])
                for row in chosen if row["row_block"] == direction_block
            ]
            cross = [
                abs(float(row["analytic_derivative"]))
                for row in chosen if row["row_block"] != direction_block
            ]
            result.append({
                "direction_block": direction_block,
                "step_V": step,
                "sampled_contact_rows": len(chosen),
                "max_analytic_vs_fd_error": analytic_fd,
                "own_block_derivative_min": min(own),
                "own_block_derivative_max": max(own),
                "own_block_derivative_spread": max(own) - min(own),
                "max_cross_block_derivative": max(cross),
            })
    return result


def summarize_direction_evidence() -> dict:
    diagnosis = BUILD_ROOT / "hole_density_diagnosis"
    variants = []
    for name, path in (
        ("raw_uncapped", diagnosis / "ab_step" / "cap_disabled" / "one_step.csv"),
        ("capped_0.1V", diagnosis / "ab_step" / "cap_0.1V" / "one_step.csv"),
    ):
        rows = read_rows(path)
        residual = []
        trial = []
        delta = []
        for row in rows:
            residual.extend(float(row[key]) for key in ("psi_residual", "phin_residual", "phip_residual"))
            trial.extend(float(row[key]) for key in ("trial_psi_residual", "trial_phin_residual", "trial_phip_residual"))
            delta.extend(float(row[key]) for key in ("delta_psi_V", "delta_phin_V", "delta_phip_V"))
        variants.append({
            "name": name,
            "node_count": len(rows),
            "direction_l2_V": norm(delta),
            "direction_max_abs_V": max(abs(value) for value in delta),
            "initial_residual_l2": norm(residual),
            "trial_residual_l2": norm(trial),
            "initial_merit_half_l2_squared": 0.5 * norm(residual) ** 2,
            "trial_merit_half_l2_squared": 0.5 * norm(trial) ** 2,
        })

    trace_rows = read_rows(diagnosis / "enforce" / "enforced_carrier_trace.csv")
    state_residual: dict[int, float] = {}
    for row in trace_rows:
        if row["event"] == "initial":
            state_residual[0] = float(row["residual_norm"])
        elif row["event"] == "iteration":
            state_residual[int(row["iteration"])] = float(row["residual_norm"])

    local_rows = read_rows(diagnosis / "enforce" / "enforced_local_updates.csv")
    first_by_iteration: dict[int, dict[str, str]] = {}
    for row in local_rows:
        first_by_iteration.setdefault(int(row["iteration"]), row)
    iterations = []
    for iteration, row in sorted(first_by_iteration.items()):
        residual_l2 = state_residual.get(iteration - 1)
        if residual_l2 is None:
            continue
        item = {
            "iteration": iteration,
            "state_residual_l2": residual_l2,
            "raw_Jd_plus_r_l2": float(row["raw_linear_residual_l2"]),
            "capped_Jd_plus_r_l2": float(row["capped_linear_residual_l2"]),
            "selected_damping": float(row["selected_damping"]),
            "line_search_attempts": int(row["line_search_attempts"]),
        }
        for prefix in ("raw", "capped"):
            closure = item[f"{prefix}_Jd_plus_r_l2"]
            r2 = residual_l2 * residual_l2
            radius = residual_l2 * closure
            item[f"{prefix}_half_l2_merit_directional_derivative_lower_bound"] = -r2 - radius
            item[f"{prefix}_half_l2_merit_directional_derivative_upper_bound"] = -r2 + radius
            item[f"{prefix}_descent_guaranteed_by_norm_bound"] = closure < residual_l2
        iterations.append(item)
    return {
        "one_step_nonlinear_merit": variants,
        "enforced_newton_direction_bounds": iterations,
        "interpretation": (
            "The bound uses r dot Jd = -||r||^2 + r dot (Jd+r). "
            "A negative upper bound proves descent for the unweighted half-L2-squared merit; "
            "a nonnegative upper bound is inconclusive, not proof of ascent."
        ),
    }


def write_summary(output: Path, config_path: Path, local_nodes: list[int], collector_nodes: list[int]) -> dict:
    row_path = output / "row_samples.csv"
    rows = read_rows(row_path)
    internal = summarize_internal(rows)
    steps = summarize_steps(rows)
    contacts = summarize_contacts(rows, collector_nodes)
    evidence = {
        "schema_version": 1,
        "scope": "read_only_jacobian_direction_audit",
        "state_contract": "accepted collector-sweep terminal state at VBE=0.70 V, VCE=3.00 V",
        "internal_rows": list(INTERNAL_ROWS),
        "local_stencil_nodes": local_nodes,
        "collector_node_count": len(collector_nodes),
        "internal_basis_audit": internal,
        "internal_aggregate_step_sensitivity": steps,
        "collector_dirichlet_audit": contacts,
        "newton_direction_evidence": summarize_direction_evidence(),
        "source_sha256": {
            "config": file_sha256(config_path),
            "state": file_sha256(BUILD_ROOT / "vela_wp3_wp5" / "m1_vce300_state.csv"),
            "row_samples": file_sha256(row_path),
            "direction_summary": file_sha256(output / "direction_summary.csv"),
        },
    }
    write_json(output / "audit_summary.json", evidence)
    lines = [
        "# Genius BJT Jacobian direction audit",
        "",
        "This is a read-only diagnostic at VBE=0.70 V and VCE=3.00 V.",
        "",
        "## Internal electron rows 4438/4440",
        "",
        "| row | basis columns | relative row L2 error | worst direction | worst absolute error |",
        "|---:|---:|---:|---|---:|",
    ]
    for item in internal:
        lines.append(
            f"| {item['row_node']} | {item['basis_columns']} | "
            f"{item['relative_row_l2_error']:.6e} | {item['worst_direction']} | "
            f"{item['worst_absolute_derivative_error']:.6e} |"
        )
    lines += [
        "",
        "## Collector Dirichlet rows",
        "",
        "| direction block | FD step (V) | sampled rows | max analytic-vs-FD error | own derivative range | max cross derivative |",
        "|---|---:|---:|---:|---:|",
    ]
    for item in contacts:
        lines.append(
            f"| {item['direction_block']} | {item['step_V']:.1e} | "
            f"{item['sampled_contact_rows']} | {item['max_analytic_vs_fd_error']:.6e} | "
            f"{item['own_block_derivative_min']:.6e} .. {item['own_block_derivative_max']:.6e} | "
            f"{item['max_cross_block_derivative']:.6e} |"
        )
    (output / "audit_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return evidence


def main() -> int:
    args = parse_args()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    config, local_nodes, collector_nodes = build_config(output)
    config_path = output / "audit_config.json"
    write_json(config_path, config)
    if not args.summarize_only:
        run_probe(args.runner, config_path, output)
    elif not (output / "row_samples.csv").exists():
        raise FileNotFoundError("--summarize-only requires an existing row_samples.csv")
    evidence = write_summary(output, config_path, local_nodes, collector_nodes)
    print(json.dumps({
        "status": "complete",
        "output": absolute(output),
        "internal_basis_audit": evidence["internal_basis_audit"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
