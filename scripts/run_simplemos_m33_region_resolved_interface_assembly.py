#!/usr/bin/env python3
"""Run the M33 SimpleMOS shared-node interface-assembly factorial audit."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import math
from pathlib import Path
import subprocess
from typing import Any


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m33_region_resolved_interface_assembly_contract_v1.json")
SOURCE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m8a_model_ablation/vela/full/workflow/vd_0p05")
MESH = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
        / "m8a_model_ablation/vela/mesh.json")
MATERIALS = REPO / "reference_tcad/transportmodels_sentaurus2022/vela/materials_sentaurus2022.json"
REFERENCE = (REPO / "reference_tcad/simplemos_sentaurus2022/model_ablation"
             / "first_round/comparisons/full_vd_0p05_comparison.csv")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m33_region_resolved_interface_assembly")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "region_resolved_interface_assembly")
RUNNER = REPO / "build-release/vela_example_runner.exe"

FACTORS = (
    "poisson_edge_coupling",
    "transport_edge_coupling",
    "transport_node_volume",
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
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def variant_name(bits: tuple[bool, bool, bool]) -> str:
    return "p{}_t{}_v{}".format(*(int(value) for value in bits))


def variants() -> list[tuple[str, dict[str, bool]]]:
    result = []
    for mask in range(8):
        bits = tuple(bool(mask & (1 << bit)) for bit in range(3))
        result.append((variant_name(bits), dict(zip(FACTORS, bits))))
    return result


def execute(config: Path, runner: Path) -> None:
    result = subprocess.run(
        [str(runner), "--config", str(config)], cwd=REPO,
        text=True, capture_output=True, encoding="utf-8", errors="replace")
    if result.returncode != 0:
        raise RuntimeError(
            f"runner failed for {config}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")


def stage_config(stage: str, name: str, flags: dict[str, bool]) -> Path:
    run_dir = OUTPUT / name
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg = read_json(SOURCE / f"{stage}.json")
    cfg["solver"]["region_resolved_interface_assembly"] = {
        "enabled": False,
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
    cfg["simplemos_m33"] = {
        "variant": name,
        "region_resolved_interface_assembly": flags,
        "default_model_changed": False,
    }
    path = run_dir / f"{stage}.json"
    write_json(path, cfg)
    return path


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
        raise FileNotFoundError(f"missing M33 output for {name}; rerun with --run")
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


def triangle_area(nodes: dict[int, dict[str, Any]], tri: dict[str, Any]) -> float:
    a, b, c = (nodes[int(node)] for node in tri["node_ids"])
    return 0.5 * abs((b["x"] - a["x"]) * (c["y"] - a["y"])
                     - (c["x"] - a["x"]) * (b["y"] - a["y"]))


def local_couple(nodes: dict[int, dict[str, Any]], tri: dict[str, Any],
                 edge: tuple[int, int]) -> float:
    a, b = (nodes[node] for node in edge)
    opposite_id = next(node for node in map(int, tri["node_ids"])
                       if node not in edge)
    o = nodes[opposite_id]
    ux, uy = a["x"] - o["x"], a["y"] - o["y"]
    vx, vy = b["x"] - o["x"], b["y"] - o["y"]
    cross = abs(ux * vy - uy * vx)
    length = math.hypot(b["x"] - a["x"], b["y"] - a["y"])
    if cross <= 1e-300 or length <= 1e-30:
        return 0.0
    cotangent = (ux * vx + uy * vy) / cross
    if cotangent < 0.0:
        return triangle_area(nodes, tri) / (3.0 * length)
    return 0.5 * cotangent * length


def geometry_audit() -> tuple[list[dict[str, Any]], dict[str, Any], set[int]]:
    mesh = read_json(MESH)
    nodes = {int(node["id"]): node for node in mesh["nodes"]}
    regions = {int(region["id"]): region for region in mesh["regions"]}
    triangles = {int(tri["id"]): tri for tri in mesh["triangles"]}
    edge_cells: dict[tuple[int, int], list[int]] = {}
    node_total = {node: 0.0 for node in nodes}
    node_si = {node: 0.0 for node in nodes}
    for tri in triangles.values():
        ids = list(map(int, tri["node_ids"]))
        area = triangle_area(nodes, tri)
        for node in ids:
            node_total[node] += area / 3.0
            if regions[int(tri["region_id"])]["material"] == "Si":
                node_si[node] += area / 3.0
        for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            edge_cells.setdefault(tuple(sorted((a, b))), []).append(int(tri["id"]))

    material_rows = read_json(MATERIALS)["materials"]
    eps = {row["name"]: float(row["eps_r"]) for row in material_rows}
    rows: list[dict[str, Any]] = []
    interface_nodes: set[int] = set()
    for edge, cell_ids in edge_cells.items():
        materials = {regions[int(triangles[cell]["region_id"])]["material"]
                     for cell in cell_ids}
        if materials != {"Si", "SiO2"}:
            continue
        couples: dict[str, float] = {"Si": 0.0, "SiO2": 0.0}
        for cell in cell_ids:
            tri = triangles[cell]
            material = regions[int(tri["region_id"])]["material"]
            couples[material] += local_couple(nodes, tri, edge)
        total = couples["Si"] + couples["SiO2"]
        legacy_poisson = 0.5 * (eps["Si"] + eps["SiO2"]) * total
        resolved_poisson = (eps["Si"] * couples["Si"]
                            + eps["SiO2"] * couples["SiO2"])
        interface_nodes.update(edge)
        rows.append({
            "node0": edge[0], "node1": edge[1],
            "si_local_couple": couples["Si"],
            "sio2_local_couple": couples["SiO2"],
            "legacy_total_couple": total,
            "legacy_over_resolved_transport_couple": total / couples["Si"],
            "legacy_over_resolved_poisson_coefficient": (
                legacy_poisson / resolved_poisson),
        })
    volume_ratios = [node_total[node] / node_si[node]
                     for node in interface_nodes if node_si[node] > 0.0]
    ratios = sorted(row["legacy_over_resolved_transport_couple"] for row in rows)
    poisson = sorted(row["legacy_over_resolved_poisson_coefficient"] for row in rows)
    volumes = sorted(volume_ratios)
    summary = {
        "mesh_node_count": len(nodes),
        "si_sio2_interface_edge_count": len(rows),
        "si_sio2_interface_node_count": len(interface_nodes),
        "transport_couple_ratio_min": min(ratios),
        "transport_couple_ratio_median": ratios[len(ratios) // 2],
        "transport_couple_ratio_max": max(ratios),
        "poisson_coefficient_ratio_min": min(poisson),
        "poisson_coefficient_ratio_median": poisson[len(poisson) // 2],
        "poisson_coefficient_ratio_max": max(poisson),
        "total_over_si_node_volume_min": min(volumes),
        "total_over_si_node_volume_median": volumes[len(volumes) // 2],
        "total_over_si_node_volume_max": max(volumes),
    }
    return rows, summary, interface_nodes


def interface_state_response(results: list[dict[str, Any]], interface_nodes: set[int]
                             ) -> list[dict[str, Any]]:
    states = {row["variant"]: {int(item["node_id"]): item
                               for item in read_csv(Path(row["state"]))}
              for row in results}
    baseline = states["p0_t0_v0"]
    response = []
    for result in results:
        state = states[result["variant"]]
        delta_psi = [float(state[node]["psi"]) - float(baseline[node]["psi"])
                     for node in interface_nodes]
        delta_n_dex = [math.log10(max(float(state[node]["electrons_m3"]), 1e-300)
                                  / max(float(baseline[node]["electrons_m3"]), 1e-300))
                       for node in interface_nodes]
        response.append({
            "variant": result["variant"],
            "max_abs_interface_delta_psi_V": max(map(abs, delta_psi)),
            "mean_abs_interface_delta_psi_V": sum(map(abs, delta_psi)) / len(delta_psi),
            "max_abs_interface_electron_density_shift_dex": max(map(abs, delta_n_dex)),
        })
    return response


def reference_currents() -> tuple[float, float]:
    rows = read_csv(REFERENCE)
    row = min(rows, key=lambda item: abs(float(item["gate_voltage_V"]) - 0.05))
    return (abs(float(row["sentaurus_current_A_per_um"])),
            abs(float(row["vela_current_A_per_um"])))


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


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    if contract["schema"] != "vela.simplemos.sdevice.m33_region_resolved_interface_assembly.v1":
        raise ValueError("unexpected M33 contract schema")

    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futures = {pool.submit(run_variant, name, flags, args.runner, args.run): name
                   for name, flags in variants()}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"completed {result['variant']}", flush=True)
    results.sort(key=lambda row: row["variant"])

    sentaurus_current, historical_vela_current = reference_currents()
    baseline = next(row for row in results if row["variant"] == "p0_t0_v0")
    baseline_error = abs(math.log10(baseline["current_A_per_um"] / sentaurus_current))
    for row in results:
        row["absolute_error_A_per_um"] = abs(row["current_A_per_um"] - sentaurus_current)
        row["absolute_error_dex"] = abs(math.log10(
            row["current_A_per_um"] / sentaurus_current))
        row["error_improvement_vs_baseline_dex"] = (
            baseline_error - row["absolute_error_dex"])
        row["log_current_shift_vs_baseline_dex"] = math.log10(
            row["current_A_per_um"] / baseline["current_A_per_um"])
        row["curve"] = portable(Path(row["curve"]))
        row["state"] = portable(Path(row["state"]))

    geometry_rows, geometry_summary, interface_nodes = geometry_audit()
    state_response = interface_state_response(results, interface_nodes)
    response_by_variant = {row["variant"]: row for row in state_response}
    for row in results:
        row.update(response_by_variant[row["variant"]])

    baseline_relative_error = (
        abs(baseline["current_A_per_um"] - historical_vela_current)
        / historical_vela_current)
    acceptance = {
        "variant_count": len(results),
        "converged_count": sum(row["converged"] for row in results),
        "baseline_replay_relative_error": baseline_relative_error,
        "finite_result_count": sum(math.isfinite(row["current_A_per_um"])
                                   for row in results),
        "si_sio2_interface_edge_count": geometry_summary[
            "si_sio2_interface_edge_count"],
    }
    limits = contract["acceptance"]
    checks = {
        "variants": acceptance["variant_count"] == limits["variant_count"],
        "converged": acceptance["converged_count"] == limits["variant_count"],
        "baseline_replay": baseline_relative_error <= limits[
            "maximum_baseline_replay_relative_error"],
        "finite": acceptance["finite_result_count"] == limits["variant_count"],
        "interface_edges": acceptance["si_sio2_interface_edge_count"] >= limits[
            "minimum_si_oxide_interface_edges"],
    }
    acceptance["checks"] = checks
    acceptance["all_checks_pass"] = all(checks.values())

    best = min(results, key=lambda row: row["absolute_error_dex"])
    all_on = next(row for row in results if row["variant"] == "p1_t1_v1")
    effects = factorial_effects(results)
    baseline_gap = baseline["absolute_error_dex"]
    report = {
        "schema": "vela.simplemos.sdevice.m33_region_resolved_interface_assembly_report.v1",
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
            "historical_cross_solver_error_dex": abs(math.log10(
                historical_vela_current / sentaurus_current)),
        },
        "acceptance": acceptance,
        "geometry": geometry_summary,
        "factorial_response": results,
        "factorial_main_effects": effects,
        "best_variant": best,
        "all_region_resolved_variant": all_on,
        "gap_closure": {
            "all_region_resolved_fraction": (
                all_on["error_improvement_vs_baseline_dex"] / baseline_gap),
            "best_incomplete_factorial_fraction": (
                best["error_improvement_vs_baseline_dex"] / baseline_gap),
        },
        "claim_policy": contract["claim_policy"],
        "artifacts": {
            "factorial_response": portable(PORTABLE / "m33_factorial_response.csv"),
            "interface_geometry": portable(PORTABLE / "m33_interface_geometry.csv"),
            "interface_state_response": portable(PORTABLE / "m33_interface_state_response.csv"),
            "factorial_main_effects": portable(PORTABLE / "m33_factorial_main_effects.csv"),
        },
    }
    PORTABLE.mkdir(parents=True, exist_ok=True)
    write_csv(PORTABLE / "m33_factorial_response.csv", results)
    write_csv(PORTABLE / "m33_interface_geometry.csv", geometry_rows)
    write_csv(PORTABLE / "m33_interface_state_response.csv", state_response)
    write_csv(PORTABLE / "m33_factorial_main_effects.csv", effects)
    write_json(PORTABLE / "m33_region_resolved_interface_assembly_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "baseline_error_dex": baseline_error,
        "all_on_error_dex": all_on["absolute_error_dex"],
        "all_on_improvement_dex": all_on["error_improvement_vs_baseline_dex"],
        "best_variant": best["variant"],
        "best_error_dex": best["absolute_error_dex"],
    }))
    return 0 if acceptance["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
