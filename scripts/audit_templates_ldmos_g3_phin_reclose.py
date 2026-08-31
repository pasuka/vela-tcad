#!/usr/bin/env python3
"""Audit the bias-dependent electron-QF feedback on the G3 maximum-gm segment.

The production physics contract is immutable.  At Vg=1.0 and 7/6 V the tool
replays the existing SSS, VSV and VVV hybrid states, checks localized analytic
Jacobian-vector products against finite differences, and performs a diagnostic
carrier-block reclose with electrostatic potential held fixed.  The reclose is
not a production solver mode and does not authorize parameter calibration.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import subprocess
from collections import defaultdict
from copy import deepcopy
from pathlib import Path
from typing import Any

try:
    from scripts.audit_templates_ldmos_g3_idvg_shift_kcl import run_sg_probe
except ModuleNotFoundError:
    from audit_templates_ldmos_g3_idvg_shift_kcl import run_sg_probe


REPO = Path(__file__).resolve().parents[1]
VARIANTS = ("SSS", "VSV", "VVV")
DAMPINGS = (1.0, 0.5, 0.25, 0.125, 0.0625)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_state(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                field: row[field] if field == "node_id" else format(float(row[field]), ".17g")
                for field in fields
            })


def state_rows(path: Path) -> list[dict[str, float | int]]:
    return [
        {
            "node_id": int(row["node_id"]),
            "psi": float(row["psi"]),
            "phin": float(row["phin"]),
            "phip": float(row["phip"]),
            "electrons_m3": float(row["electrons_m3"]),
            "holes_m3": float(row["holes_m3"]),
        }
        for row in read_csv(path)
    ]


def runner_environment() -> dict[str, str]:
    environment = dict(os.environ)
    if os.name == "nt":
        environment["PATH"] = os.pathsep.join([
            r"D:\msys64\ucrt64\bin",
            r"D:\msys64\usr\bin",
            environment.get("PATH", ""),
        ])
    environment["VELA_LINEAR_SOLVER"] = "sparselu"
    return environment


def run_config(
    runner: Path, config: dict[str, Any], directory: Path, stem: str,
) -> dict[str, Any]:
    directory.mkdir(parents=True, exist_ok=True)
    config_path = directory / f"{stem}.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve()),
         "--log", str((directory / f"{stem}.log").resolve())],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=runner_environment(),
        check=False,
    )
    (directory / f"{stem}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (directory / f"{stem}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"{stem} failed: {completed.stderr or completed.stdout}")
    return json.loads(completed.stdout)


def probe_config(
    baseline: dict[str, Any], state: Path, bias: float,
    simulation_type: str, output: Path,
) -> dict[str, Any]:
    config = deepcopy(baseline)
    config.pop("sweep", None)
    config["simulation_type"] = simulation_type
    config["state_file"] = str(state.resolve())
    config["output_csv"] = str(output.resolve())
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = bias
    config["_comment"] = (
        "Templates/LDMOS G3 maximum-gm phin audit; IALMob and predictor off."
    )
    return config


def l2(rows: list[dict[str, str]], field: str) -> float:
    return math.sqrt(sum(float(row[field]) ** 2 for row in rows))


def l1(rows: list[dict[str, str]], field: str) -> float:
    return sum(abs(float(row[field])) for row in rows)


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def cosine(left: list[float], right: list[float]) -> float:
    numerator = sum(a * b for a, b in zip(left, right))
    left_norm = math.sqrt(sum(value * value for value in left))
    right_norm = math.sqrt(sum(value * value for value in right))
    denominator = left_norm * right_norm
    return numerator / denominator if denominator else 0.0


def mesh_node_classes(mesh_path: Path) -> dict[str, Any]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    material_by_region = {
        int(region["id"]): region.get("material", "").lower()
        for region in mesh["regions"]
    }
    silicon_regions = {
        region_id
        for region_id, material in material_by_region.items()
        if material in {"si", "silicon"}
    }
    node_regions: dict[int, set[int]] = defaultdict(set)
    node_neighbors: dict[int, set[int]] = defaultdict(set)
    silicon_nodes: set[int] = set()
    for cell in mesh["triangles"]:
        region_id = int(cell["region_id"])
        nodes = [int(node) for node in cell["node_ids"]]
        for node in nodes:
            node_regions[node].add(region_id)
            if region_id in silicon_regions:
                silicon_nodes.add(node)
            node_neighbors[node].update(other for other in nodes if other != node)

    contact_names_by_node: dict[int, set[str]] = defaultdict(set)
    for contact in mesh["contacts"]:
        for node in contact.get("node_ids", []):
            contact_names_by_node[int(node)].add(str(contact["name"]))
    contact_nodes = set(contact_names_by_node)

    interface_nodes = {
        node
        for node in silicon_nodes
        if any(region not in silicon_regions for region in node_regions[node])
    }
    interface_one_ring = set(interface_nodes)
    for node in interface_nodes:
        interface_one_ring.update(node_neighbors[node] & silicon_nodes)

    contact_one_ring_names: dict[int, set[str]] = defaultdict(set)
    for node, names in contact_names_by_node.items():
        for neighbor in {node} | node_neighbors[node]:
            if neighbor in silicon_nodes:
                contact_one_ring_names[neighbor].update(names)

    return {
        "free_nodes": silicon_nodes - contact_nodes,
        "interface_nodes": interface_nodes,
        "interface_one_ring": interface_one_ring,
        "contact_names_by_node": contact_names_by_node,
        "contact_one_ring_names": contact_one_ring_names,
    }


def term_summary(
    rows: list[dict[str, str]], node_classes: dict[str, Any], limit: int = 12
) -> dict[str, Any]:
    ranked = sorted(rows, key=lambda row: abs(float(row["electron_residual"])), reverse=True)
    residual_energy = sum(float(row["electron_residual"]) ** 2 for row in rows)

    def energy_fraction(nodes: set[int]) -> float:
        selected = sum(
            float(row["electron_residual"]) ** 2
            for row in rows
            if int(row["node_id"]) in nodes
        )
        return selected / residual_energy if residual_energy else 0.0

    return {
        "electron_residual_l2": l2(rows, "electron_residual"),
        "electron_residual_l1": l1(rows, "electron_residual"),
        "electron_flux_l1": l1(rows, "electron_flux"),
        "electron_recombination_l1": l1(rows, "electron_recombination"),
        "hole_residual_l2": l2(rows, "hole_residual"),
        "electron_residual_energy_fraction": {
            "silicon_dielectric_interface": energy_fraction(
                node_classes["interface_nodes"]
            ),
            "silicon_dielectric_interface_one_ring": energy_fraction(
                node_classes["interface_one_ring"]
            ),
            "contact_one_ring": energy_fraction(
                set(node_classes["contact_one_ring_names"])
            ),
        },
        "top_electron_residual_nodes": [
            {
                "node_id": int(row["node_id"]),
                "x_um": float(row["x"]),
                "y_um": float(row["y"]),
                "electron_residual": float(row["electron_residual"]),
                "electron_flux_abs_sum": float(row["electron_flux_abs_sum"]),
                "silicon_dielectric_interface": int(row["node_id"])
                in node_classes["interface_nodes"],
                "silicon_dielectric_interface_one_ring": int(row["node_id"])
                in node_classes["interface_one_ring"],
                "contact_names": sorted(
                    node_classes["contact_names_by_node"].get(
                        int(row["node_id"]), set()
                    )
                ),
                "contact_one_ring_names": sorted(
                    node_classes["contact_one_ring_names"].get(
                        int(row["node_id"]), set()
                    )
                ),
            }
            for row in ranked[:limit]
        ],
    }


def target_step_metrics(
    start: list[dict[str, float | int]],
    target: list[dict[str, float | int]],
    step_rows: list[dict[str, str]],
    free_nodes: set[int],
) -> dict[str, float]:
    start_by_node = {int(row["node_id"]): row for row in start}
    target_by_node = {int(row["node_id"]): row for row in target}
    carrier = {
        int(row["node_id"]): row
        for row in step_rows
        if row["mode"] == "carrier_only"
    }
    ordered = sorted(free_nodes & set(start_by_node) & set(target_by_node) & set(carrier))
    desired = [
        float(target_by_node[node]["phin"]) - float(start_by_node[node]["phin"])
        for node in ordered
    ]
    update = [float(carrier[node]["delta_phin_V"]) for node in ordered]
    before = math.sqrt(sum(value * value for value in desired) / max(len(desired), 1))
    after_error = [desired_value - update_value for desired_value, update_value in zip(desired, update)]
    after = math.sqrt(sum(value * value for value in after_error) / max(len(after_error), 1))
    improvement = 0.0 if before <= 1.0e-300 else 1.0 - after / before
    return {
        "free_silicon_nodes": len(ordered),
        "phin_update_target_cosine": cosine(update, desired),
        "phin_target_rmse_before_V": before,
        "phin_target_rmse_after_full_step_V": after,
        "phin_target_rmse_improvement_fraction": improvement,
        "phin_update_l2_V": math.sqrt(sum(value * value for value in update)),
        "phin_target_l2_V": math.sqrt(sum(value * value for value in desired)),
    }


def row_hotspot_summary(
    row_path: Path, hotspot_nodes: list[int], limit: int = 12,
) -> list[dict[str, float | int]]:
    wanted = set(hotspot_nodes)
    selected = [row for row in read_csv(row_path) if int(row["node_id"]) in wanted]
    selected.sort(key=lambda row: abs(float(row["electron_residual"])), reverse=True)
    result = []
    for row in selected[:limit]:
        diagonal = abs(float(row["electron_diagonal"]))
        row_sum = float(row["electron_row_abs_sum"])
        result.append({
            "node_id": int(row["node_id"]),
            "x_um": float(row["x"]),
            "y_um": float(row["y"]),
            "electron_residual": float(row["electron_residual"]),
            "electron_diagonal_fraction": diagonal / max(row_sum, 1.0e-300),
            "electron_offdiag_over_diagonal": (
                float(row["electron_offdiag_abs_sum"]) / max(diagonal, 1.0e-300)
            ),
            "raw_delta_phin_V": float(row["raw_delta_phin_V"]),
            "capped_delta_phin_V": float(row["capped_delta_phin_V"]),
        })
    return result


def point_directory(root: Path, bias: float) -> Path:
    return root / f"vg_{bias:.6f}".replace(".", "p")


def interpolate_trial_state(
    current: list[dict[str, float | int]],
    step_rows: list[dict[str, str]],
    damping: float,
) -> list[dict[str, float | int]]:
    carrier = {
        int(row["node_id"]): row
        for row in step_rows
        if row["mode"] == "carrier_only"
    }
    result = []
    for row in current:
        node = int(row["node_id"])
        trial = carrier[node]
        result.append({
            **row,
            "phin": float(row["phin"]) + damping * (
                float(trial["trial_phin"]) - float(row["phin"])
            ),
            "phip": float(row["phip"]) + damping * (
                float(trial["trial_phip"]) - float(row["phip"])
            ),
        })
    return result


def run_terms(
    runner: Path, baseline: dict[str, Any], state: Path, bias: float,
    directory: Path, stem: str,
) -> tuple[dict[str, Any], list[dict[str, str]], Path]:
    output = directory / f"{stem}.csv"
    status = run_config(
        runner,
        probe_config(baseline, state, bias, "newton_carrier_term_probe", output),
        directory,
        stem,
    )
    return status, read_csv(output), output


def reclose_carriers(
    runner: Path,
    baseline: dict[str, Any],
    initial_state: Path,
    target_state: Path,
    bias: float,
    free_nodes: set[int],
    directory: Path,
    max_iterations: int,
) -> dict[str, Any]:
    current_path = directory / "iteration_00_state.csv"
    current_rows = state_rows(initial_state)
    target_rows = state_rows(target_state)
    write_state(current_path, current_rows)
    _, current_terms, _ = run_terms(
        runner, baseline, current_path, bias, directory / "iteration_00", "terms"
    )
    current_electron = l2(current_terms, "electron_residual")
    history: list[dict[str, Any]] = [{
        "iteration": 0,
        "electron_residual_l2": current_electron,
        "hole_residual_l2": l2(current_terms, "hole_residual"),
        "accepted_damping": 0.0,
        "state": str(current_path.resolve()),
    }]
    stop_reason = "max_iterations"
    for iteration in range(1, max_iterations + 1):
        iteration_dir = directory / f"iteration_{iteration:02d}"
        block_csv = iteration_dir / "carrier_step.csv"
        block = probe_config(
            baseline, current_path, bias, "newton_block_step_probe", block_csv
        )
        block["block_modes"] = ["carrier_only"]
        block_status = run_config(runner, block, iteration_dir, "carrier_step")
        step_rows = read_csv(block_csv)
        candidates = []
        for damping in DAMPINGS:
            label = f"damping_{damping:.4f}".replace(".", "p")
            candidate_path = iteration_dir / label / "state.csv"
            write_state(
                candidate_path,
                interpolate_trial_state(current_rows, step_rows, damping),
            )
            _, candidate_terms, _ = run_terms(
                runner, baseline, candidate_path, bias,
                iteration_dir / label, "terms",
            )
            candidates.append({
                "damping": damping,
                "state": candidate_path,
                "rows": state_rows(candidate_path),
                "electron_residual_l2": l2(candidate_terms, "electron_residual"),
                "hole_residual_l2": l2(candidate_terms, "hole_residual"),
            })
        best = min(candidates, key=lambda item: item["electron_residual_l2"])
        if float(best["electron_residual_l2"]) >= current_electron * (1.0 - 1.0e-10):
            stop_reason = "no_electron_residual_decrease"
            break
        current_path = Path(best["state"])
        current_rows = list(best["rows"])
        current_electron = float(best["electron_residual_l2"])
        history.append({
            "iteration": iteration,
            "electron_residual_l2": current_electron,
            "hole_residual_l2": float(best["hole_residual_l2"]),
            "accepted_damping": float(best["damping"]),
            "raw_step_norm": block_status["block_steps"][0]["raw_step_norm"],
            "capped_step_norm": block_status["block_steps"][0]["step_norm"],
            "state": str(current_path.resolve()),
        })
        if current_electron <= 1.0e-9:
            stop_reason = "electron_residual_target"
            break

    final_probe = run_sg_probe(
        runner, baseline, current_path, bias, directory / "final_sg"
    )
    target = {int(row["node_id"]): row for row in target_rows}
    current = {int(row["node_id"]): row for row in current_rows}
    phin_error = [
        float(current[node]["phin"]) - float(target[node]["phin"])
        for node in sorted(free_nodes & set(current) & set(target))
    ]
    initial = {int(row["node_id"]): row for row in state_rows(initial_state)}
    phip_delta = [
        float(current[node]["phip"]) - float(initial[node]["phip"])
        for node in sorted(free_nodes & set(current) & set(initial))
    ]
    return {
        "stop_reason": stop_reason,
        "iterations_accepted": len(history) - 1,
        "history": history,
        "final_state": str(current_path.resolve()),
        "final_current_A_per_um": float(final_probe["drain_cut"]["total_A_per_um"]),
        "final_phin_error_vs_vvv": {
            "rmse_V": math.sqrt(
                sum(value * value for value in phin_error) / max(len(phin_error), 1)
            ),
            "median_V": statistics.median(phin_error),
            "p95_abs_V": percentile([abs(value) for value in phin_error], 0.95),
            "maximum_abs_V": max((abs(value) for value in phin_error), default=0.0),
        },
        "phip_change_from_initial": {
            "rmse_V": math.sqrt(
                sum(value * value for value in phip_delta) / max(len(phip_delta), 1)
            ),
            "maximum_abs_V": max((abs(value) for value in phip_delta), default=0.0),
        },
        "artifacts": final_probe["artifacts"],
    }


def audit_point(
    runner: Path,
    baseline: dict[str, Any],
    source_root: Path,
    point: dict[str, Any],
    output: Path,
    node_classes: dict[str, Any],
    max_reclose_iterations: int,
) -> dict[str, Any]:
    free_nodes = node_classes["free_nodes"]
    bias = float(point["bias_V"])
    source = point_directory(source_root, bias)
    destination = output / f"vg_{bias:.6f}".replace(".", "p")
    paths = {name: source / name / "hybrid_state.csv" for name in VARIANTS}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    variant_results: dict[str, Any] = {}
    term_rows_by_variant: dict[str, list[dict[str, str]]] = {}
    for name in VARIANTS:
        case = destination / name
        _, term_rows, term_path = run_terms(
            runner, baseline, paths[name], bias, case, "carrier_terms"
        )
        term_rows_by_variant[name] = term_rows
        row_path = case / "carrier_rows.csv"
        row_status = run_config(
            runner,
            probe_config(
                baseline, paths[name], bias, "newton_carrier_row_probe", row_path
            ),
            case,
            "carrier_rows",
        )
        step_path = case / "carrier_step.csv"
        step_config = probe_config(
            baseline, paths[name], bias, "newton_block_step_probe", step_path
        )
        step_config["block_modes"] = ["carrier_only"]
        step_status = run_config(runner, step_config, case, "carrier_step")
        step_rows = read_csv(step_path)
        terms = term_summary(term_rows, node_classes)
        hotspots = [row["node_id"] for row in terms["top_electron_residual_nodes"]]
        variant_results[name] = {
            "terms": terms,
            "carrier_step": step_status["block_steps"][0],
            "step_toward_vvv": target_step_metrics(
                state_rows(paths[name]), state_rows(paths["VVV"]), step_rows, free_nodes
            ),
            "hotspot_rows": row_hotspot_summary(row_path, hotspots),
            "artifacts": {
                "state": str(paths[name].resolve()),
                "carrier_terms": str(term_path.resolve()),
                "carrier_rows": str(row_path.resolve()),
                "carrier_step": str(step_path.resolve()),
            },
        }

    vsv_hotspots = [
        row["node_id"]
        for row in variant_results["VSV"]["terms"]["top_electron_residual_nodes"][:8]
    ]
    sss_hotspots = [
        row["node_id"]
        for row in variant_results["SSS"]["terms"]["top_electron_residual_nodes"][:8]
    ]
    jvp_path = destination / "jvp" / "localized_jvp.csv"
    jvp_config = probe_config(
        baseline, paths["VSV"], bias, "newton_jvp_probe", jvp_path
    )
    jvp_config["directions"] = [
        {
            "name": "vsv_hotspot_phin",
            "mode": "phin",
            "node_ids": vsv_hotspots,
            "adjacent_cell_rings": 1,
            "exclude_contacts": True,
            "amplitude_V": 1.0e-6,
        },
        {
            "name": "sss_hotspot_phin",
            "mode": "phin",
            "node_ids": sss_hotspots,
            "adjacent_cell_rings": 1,
            "exclude_contacts": True,
            "amplitude_V": 1.0e-6,
        },
        {
            "name": "vsv_hotspot_psi_minus_phin",
            "mode": "psi_minus_phin",
            "node_ids": vsv_hotspots,
            "adjacent_cell_rings": 1,
            "exclude_contacts": True,
            "amplitude_V": 1.0e-6,
        },
    ]
    jvp_status = run_config(runner, jvp_config, destination / "jvp", "localized_jvp")
    jvp_rows = read_csv(jvp_path)

    reclose = reclose_carriers(
        runner,
        baseline,
        paths["VSV"],
        paths["VVV"],
        bias,
        free_nodes,
        destination / "carrier_reclose",
        max_reclose_iterations,
    )
    reclose["current_ratio_to_sentaurus"] = (
        reclose["final_current_A_per_um"]
        / float(point["sentaurus_terminal_A_per_um"])
    )
    reclose["current_ratio_to_vvv"] = (
        reclose["final_current_A_per_um"]
        / float(point["vela_curve_A_per_um"])
    )
    return {
        "bias_V": bias,
        "sentaurus_terminal_A_per_um": float(point["sentaurus_terminal_A_per_um"]),
        "vela_curve_A_per_um": float(point["vela_curve_A_per_um"]),
        "variants": variant_results,
        "localized_jvp": {
            "max_relative_error": max(float(row["relative_error"]) for row in jvp_rows),
            "max_phin_relative_error": max(
                float(row["phin_relative_error"]) for row in jvp_rows
            ),
            "directions": jvp_status["directions"],
            "artifact": str(jvp_path.resolve()),
        },
        "fixed_psi_carrier_reclose_from_vsv": reclose,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument(
        "--state-feedback-root", type=Path, action="append", required=True
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-reclose-iterations", type=int, default=6)
    args = parser.parse_args()
    if args.max_reclose_iterations < 1:
        raise ValueError("max reclose iterations must be positive")

    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8"))
    if "ialmob" in json.dumps(baseline["solver"].get("mobility", {})).lower():
        raise ValueError("phin reclose audit requires IALMob disabled")
    if "predictor" in json.dumps(baseline.get("sweep", {})).lower():
        raise ValueError("phin reclose audit requires predictor disabled")
    if baseline["solver"].get("contact_boundary_reconstruction") != "legacy_node_local":
        raise ValueError("phin reclose audit requires legacy_node_local")
    if baseline.get("mesh_geometry", {}).get("carrier_transport_couple_profile") != (
        "templates_ldmos_external_averagebox"
    ):
        raise ValueError("phin reclose audit requires external AverageBox")

    sources = []
    for root in args.state_feedback_root:
        summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
        if len(summary["points"]) != 1:
            raise ValueError(f"expected one point in {root / 'summary.json'}")
        sources.append((root, summary["points"][0]))
    sources.sort(key=lambda item: float(item[1]["bias_V"]))
    if [round(float(item[1]["bias_V"]), 12) for item in sources] != [1.0, round(7.0 / 6.0, 12)]:
        raise ValueError("audit requires exact Vg=1.0 and 7/6 V endpoints")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    node_classes = mesh_node_classes(Path(baseline["mesh_file"]))
    points = [
        audit_point(
            args.runner, baseline, root, point, output, node_classes,
            args.max_reclose_iterations,
        )
        for root, point in sources
    ]
    low, high = points
    growth = {
        variant: (
            high["variants"][variant]["terms"]["electron_residual_l2"]
            / low["variants"][variant]["terms"]["electron_residual_l2"]
        )
        for variant in VARIANTS
    }
    low_reclose = low["fixed_psi_carrier_reclose_from_vsv"]
    high_reclose = high["fixed_psi_carrier_reclose_from_vsv"]
    delta_bias = high["bias_V"] - low["bias_V"]
    reclose_gm = (
        high_reclose["final_current_A_per_um"]
        - low_reclose["final_current_A_per_um"]
    ) / delta_bias
    sentaurus_gm = (
        high["sentaurus_terminal_A_per_um"] - low["sentaurus_terminal_A_per_um"]
    ) / delta_bias
    vela_gm = (
        high["vela_curve_A_per_um"] - low["vela_curve_A_per_um"]
    ) / delta_bias
    report = {
        "schema": "vela.templates_ldmos.g3_phin_reclose_audit.v1",
        "contracts": {
            "ialmob": "disabled",
            "predictor": "disabled",
            "physics_parameters_changed": False,
            "contact_boundary_reconstruction": "legacy_node_local",
            "carrier_transport_couple_profile": "templates_ldmos_external_averagebox",
            "reclose": "diagnostic carrier-block Newton with psi fixed",
        },
        "points": points,
        "aggregate": {
            "electron_residual_l2_growth_high_over_low": growth,
            "sentaurus_gm_A_per_um_V": sentaurus_gm,
            "vela_gm_A_per_um_V": vela_gm,
            "fixed_psi_carrier_reclose_gm_A_per_um_V": reclose_gm,
            "fixed_psi_carrier_reclose_gm_ratio_to_sentaurus": reclose_gm / sentaurus_gm,
            "fixed_psi_carrier_reclose_gm_ratio_to_vela": reclose_gm / vela_gm,
            "maximum_localized_jvp_relative_error": max(
                point["localized_jvp"]["max_relative_error"] for point in points
            ),
            "classification": "phin_continuity_fixed_point_feedback",
            "ialmob_authorized": False,
            "parameter_calibration_authorized": False,
            "ledger_status": "draft",
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"summary": str(summary_path), **report["aggregate"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
