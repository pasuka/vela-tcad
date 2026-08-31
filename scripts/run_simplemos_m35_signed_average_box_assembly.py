#!/usr/bin/env python3
"""Run the M35 signed-AverageBox region-resolved assembly factorial."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts import run_simplemos_m33_region_resolved_interface_assembly as m33


CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m35_signed_average_box_assembly_contract_v1.json")
SOURCE = m33.SOURCE
MESH = m33.MESH
REFERENCE = m33.REFERENCE
M34_NODE_ORACLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                   / "sentaurus_interface_box_probe/m34_interface_node_measures.csv")
M33_ROOT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "region_resolved_interface_assembly")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m35_signed_average_box_assembly")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "signed_average_box_assembly")
RUNNER = REPO / "build-release/vela_example_runner.exe"

FACTORS = (
    "poisson_edge_coupling",
    "transport_edge_coupling",
    "transport_signed_average_box_node_volume",
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def variants() -> list[tuple[str, dict[str, bool]]]:
    rows = []
    for mask in range(8):
        bits = tuple(bool(mask & (1 << bit)) for bit in range(3))
        name = "p{}_t{}_v{}".format(*(int(value) for value in bits))
        rows.append((name, dict(zip(FACTORS, bits))))
    return rows


def stage_config(stage: str, name: str, flags: dict[str, bool]) -> Path:
    run_dir = OUTPUT / name
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = read_json(SOURCE / f"{stage}.json")
    cfg["solver"]["region_resolved_interface_assembly"] = {
        "enabled": False,
        "transport_node_volume": False,
        **flags,
    }
    cfg["output_csv"] = str((run_dir / f"{stage}.csv").resolve())
    cfg["log_file"] = str((run_dir / f"{stage}.log").resolve())
    cfg["sweep"]["write_state_file"] = str(
        (run_dir / f"{stage}_accepted_state.csv").resolve())
    if stage == "10_drain_ramp":
        cfg["sweep"]["initial_state_file"] = str(
            (run_dir / "00_equilibrium_accepted_state.csv").resolve())
    elif stage == "20_gate_sweep":
        cfg["sweep"]["initial_state_file"] = str(
            (run_dir / "10_drain_ramp_accepted_state.csv").resolve())
        cfg["sweep"].update({
            "start": 0.0,
            "stop": 0.05,
            "step": 0.05,
            "min_step": min(float(cfg["sweep"].get("min_step", 0.05)), 0.05),
            "max_step": 0.05,
            "bias_points": [0.0, 0.05],
        })
    cfg["simplemos_m35"] = {
        "variant": name,
        "region_resolved_interface_assembly": flags,
        "sentaurus_geometry_oracle": portable(M34_NODE_ORACLE),
        "default_model_changed": False,
    }
    path = run_dir / f"{stage}.json"
    write_json(path, cfg)
    return path


def execute(config: Path, runner: Path) -> None:
    result = subprocess.run(
        [str(runner), "--config", str(config)], cwd=REPO,
        text=True, capture_output=True, encoding="utf-8", errors="replace")
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {config}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")


def run_variant(name: str, flags: dict[str, bool], runner: Path,
                run: bool) -> dict[str, Any]:
    stages = ["00_equilibrium", "10_drain_ramp", "20_gate_sweep"]
    configs = [stage_config(stage, name, flags) for stage in stages]
    if run:
        for config in configs:
            execute(config, runner)
    curve = OUTPUT / name / "20_gate_sweep.csv"
    state = OUTPUT / name / "20_gate_sweep_accepted_state.csv"
    if not curve.exists() or not state.exists():
        raise FileNotFoundError(f"missing M35 output for {name}; rerun with --run")
    rows = read_csv(curve)
    final = min(rows, key=lambda row: abs(float(row["bias_V"]) - 0.05))
    if not math.isclose(float(final["bias_V"]), 0.05, abs_tol=1.0e-12):
        raise RuntimeError(f"{name} did not reach the exact gate bias")
    return {
        "variant": name,
        **flags,
        "converged": str(final["converged"]).lower() in ("1", "true"),
        "iterations": int(final["iterations"]),
        "current_A_per_um": abs(float(final["current_total_A_per_um"])),
        "state": state,
        "curve": curve,
    }


def triangle_signed_measures(points: list[tuple[float, float]]) -> list[float]:
    def distance_squared(a: int, b: int) -> float:
        return sum((points[b][axis] - points[a][axis]) ** 2 for axis in (0, 1))

    def cotangent(vertex: int, a: int, b: int) -> float:
        ux, uy = (points[a][axis] - points[vertex][axis] for axis in (0, 1))
        vx, vy = (points[b][axis] - points[vertex][axis] for axis in (0, 1))
        return (ux * vx + uy * vy) / abs(ux * vy - uy * vx)

    result = []
    for index in range(3):
        j, k = (index + 1) % 3, (index + 2) % 3
        result.append(0.125 * (
            distance_squared(index, k) * cotangent(j, index, k)
            + distance_squared(index, j) * cotangent(k, index, j)))
    return result


def sentaurus_measure_oracle() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    mesh = read_json(MESH)
    nodes = {int(node["id"]): node for node in mesh["nodes"]}
    regions = {int(region["id"]): region for region in mesh["regions"]}
    totals = {node: 0.0 for node in nodes}
    si_nodes: set[int] = set()
    max_cell_closure = 0.0
    negative_cell_share_count = 0
    for triangle in mesh["triangles"]:
        if regions[int(triangle["region_id"])]["material"] != "Si":
            continue
        ids = list(map(int, triangle["node_ids"]))
        si_nodes.update(ids)
        points = [(float(nodes[node]["x"]), float(nodes[node]["y"])) for node in ids]
        measures = triangle_signed_measures(points)
        area = m33.triangle_area(nodes, triangle)
        max_cell_closure = max(max_cell_closure, abs(sum(measures) - area))
        negative_cell_share_count += sum(value < 0.0 for value in measures)
        for node, measure in zip(ids, measures):
            totals[node] += measure

    rows = []
    for oracle in read_csv(M34_NODE_ORACLE):
        node = int(oracle["node"])
        sentaurus = float(oracle["sentaurus_si_measure_um2"])
        calculated = totals[node]
        absolute = abs(calculated - sentaurus)
        relative = absolute / max(abs(sentaurus), 1e-300)
        rows.append({
            "node": node,
            "x_um": oracle["x_um"],
            "y_um": oracle["y_um"],
            "sentaurus_si_measure_um2": sentaurus,
            "vela_signed_average_box_si_measure_um2": calculated,
            "absolute_error_um2": absolute,
            "relative_error": relative,
            "barycentric_si_measure_um2": oracle["barycentric_si_measure_um2"],
        })
    summary = {
        "sentaurus_interface_node_count": len(rows),
        "maximum_sentaurus_measure_absolute_error_um2": max(
            row["absolute_error_um2"] for row in rows),
        "maximum_sentaurus_measure_relative_error": max(
            row["relative_error"] for row in rows),
        "maximum_si_cell_area_closure_error_um2": max_cell_closure,
        "negative_si_cell_share_count": negative_cell_share_count,
        "nonpositive_assembled_si_node_measure_count": sum(
            totals[node] <= 0.0 for node in si_nodes
        ),
        "minimum_assembled_si_node_measure_um2": min(
            totals[node] for node in si_nodes
        ),
    }
    return rows, summary


def factorial_effects(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    effects = []
    for factor in FACTORS:
        enabled = [row["log_current_shift_vs_baseline_dex"]
                   for row in rows if row[factor]]
        disabled = [row["log_current_shift_vs_baseline_dex"]
                    for row in rows if not row[factor]]
        effects.append({
            "factor": factor,
            "mean_enabled_log_current_shift_dex": sum(enabled) / len(enabled),
            "mean_disabled_log_current_shift_dex": sum(disabled) / len(disabled),
            "main_effect_log_current_dex": (
                sum(enabled) / len(enabled) - sum(disabled) / len(disabled)),
        })
    return effects


def m33_comparison(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    old = {row["variant"]: row for row in read_csv(
        M33_ROOT / "m33_factorial_response.csv")}
    comparison = []
    for row in rows:
        previous = old[row["variant"]]
        comparison.append({
            "variant": row["variant"],
            "m33_barycentric_current_A_per_um": previous["current_A_per_um"],
            "m35_signed_average_box_current_A_per_um": row["current_A_per_um"],
            "m33_barycentric_error_dex": previous["absolute_error_dex"],
            "m35_signed_average_box_error_dex": row["absolute_error_dex"],
            "signed_minus_barycentric_log_current_dex": math.log10(
                row["current_A_per_um"] / float(previous["current_A_per_um"])),
            "error_improvement_vs_m33_dex": (
                float(previous["absolute_error_dex"]) - row["absolute_error_dex"]),
        })
    return comparison


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    if contract["schema"] != "vela.simplemos.sdevice.m35_signed_average_box_assembly.v1":
        raise ValueError("unexpected M35 contract schema")

    results = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(run_variant, name, flags, args.runner, args.run): name
                   for name, flags in variants()}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"completed {result['variant']}", flush=True)
    results.sort(key=lambda row: row["variant"])

    sentaurus_current, historical_vela_current = m33.reference_currents()
    baseline = next(row for row in results if row["variant"] == "p0_t0_v0")
    baseline_error = abs(math.log10(baseline["current_A_per_um"] / sentaurus_current))
    for row in results:
        row["absolute_error_A_per_um"] = abs(row["current_A_per_um"] - sentaurus_current)
        row["absolute_error_dex"] = abs(math.log10(row["current_A_per_um"] / sentaurus_current))
        row["error_improvement_vs_baseline_dex"] = baseline_error - row["absolute_error_dex"]
        row["log_current_shift_vs_baseline_dex"] = math.log10(
            row["current_A_per_um"] / baseline["current_A_per_um"])
        row["curve"] = portable(Path(row["curve"]))
        row["state"] = portable(Path(row["state"]))

    geometry_rows, geometry_summary, interface_nodes = m33.geometry_audit()
    state_response = m33.interface_state_response(results, interface_nodes)
    state_by_variant = {row["variant"]: row for row in state_response}
    for row in results:
        row.update(state_by_variant[row["variant"]])
    measure_rows, measure_summary = sentaurus_measure_oracle()
    effects = factorial_effects(results)
    comparison = m33_comparison(results)

    limits = contract["acceptance"]
    baseline_relative_error = abs(
        baseline["current_A_per_um"] - historical_vela_current) / historical_vela_current
    acceptance = {
        "variant_count": len(results),
        "converged_count": sum(row["converged"] for row in results),
        "finite_result_count": sum(math.isfinite(row["current_A_per_um"]) for row in results),
        "baseline_replay_relative_error": baseline_relative_error,
        "si_sio2_interface_edge_count": geometry_summary["si_sio2_interface_edge_count"],
        **measure_summary,
    }
    checks = {
        "variants": acceptance["variant_count"] == limits["variant_count"],
        "converged": acceptance["converged_count"] == limits["variant_count"],
        "finite": acceptance["finite_result_count"] == limits["variant_count"],
        "baseline_replay": baseline_relative_error <= limits["maximum_baseline_replay_relative_error"],
        "interface_edges": acceptance["si_sio2_interface_edge_count"] >= limits["minimum_si_oxide_interface_edges"],
        "sentaurus_nodes": measure_summary["sentaurus_interface_node_count"] >= limits["minimum_sentaurus_interface_nodes"],
        "sentaurus_measure": measure_summary["maximum_sentaurus_measure_relative_error"] <= limits["maximum_sentaurus_measure_relative_error"],
    }
    acceptance["checks"] = checks
    acceptance["all_checks_pass"] = all(checks.values())

    best = min(results, key=lambda row: row["absolute_error_dex"])
    all_on = next(row for row in results if row["variant"] == "p1_t1_v1")
    report = {
        "schema": "vela.simplemos.sdevice.m35_signed_average_box_assembly_report.v1",
        "status": "complete" if acceptance["all_checks_pass"] else "failed_acceptance",
        "execution": {
            "new_sentaurus_execution": False,
            "cpp_changed": True,
            "default_model_changed": False,
            "self_consistent_workflows": len(results),
        },
        "reference": {
            "sentaurus_current_A_per_um": sentaurus_current,
            "historical_vela_current_A_per_um": historical_vela_current,
            "m34_measure_oracle": portable(M34_NODE_ORACLE),
        },
        "acceptance": acceptance,
        "geometry": geometry_summary,
        "sentaurus_measure_oracle": measure_summary,
        "factorial_response": results,
        "factorial_main_effects": effects,
        "m33_barycentric_comparison": comparison,
        "best_variant": best,
        "all_region_resolved_variant": all_on,
        "gap_closure": {
            "all_region_resolved_fraction": all_on["error_improvement_vs_baseline_dex"] / baseline_error,
            "best_factorial_fraction": best["error_improvement_vs_baseline_dex"] / baseline_error,
        },
        "claim_policy": contract["claim_policy"],
        "artifacts": {
            "factorial_response": portable(PORTABLE / "m35_factorial_response.csv"),
            "factorial_main_effects": portable(PORTABLE / "m35_factorial_main_effects.csv"),
            "interface_state_response": portable(PORTABLE / "m35_interface_state_response.csv"),
            "sentaurus_measure_oracle": portable(PORTABLE / "m35_sentaurus_measure_oracle.csv"),
            "m33_barycentric_comparison": portable(PORTABLE / "m35_m33_barycentric_comparison.csv"),
        },
    }
    PORTABLE.mkdir(parents=True, exist_ok=True)
    write_csv(PORTABLE / "m35_factorial_response.csv", results)
    write_csv(PORTABLE / "m35_factorial_main_effects.csv", effects)
    write_csv(PORTABLE / "m35_interface_state_response.csv", state_response)
    write_csv(PORTABLE / "m35_sentaurus_measure_oracle.csv", measure_rows)
    write_csv(PORTABLE / "m35_m33_barycentric_comparison.csv", comparison)
    write_json(PORTABLE / "m35_signed_average_box_assembly_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "maximum_measure_relative_error": measure_summary["maximum_sentaurus_measure_relative_error"],
        "baseline_error_dex": baseline_error,
        "all_on_error_dex": all_on["absolute_error_dex"],
        "all_on_improvement_dex": all_on["error_improvement_vs_baseline_dex"],
        "best_variant": best["variant"],
        "best_error_dex": best["absolute_error_dex"],
    }))
    return 0 if acceptance["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
