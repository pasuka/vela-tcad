#!/usr/bin/env python3
"""Run and analyze SimpleMOS M13 spatial-attribution diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
M10_ROOT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
            / "m10_fixed_state_replay")
M11_ROOT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
            / "m11_mobility_factorial/frozen_vela")
M8_ROOT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
           / "m8_original_physics/vela")
M8_DEEP = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
           / "m8_deep_off_diagnostics/vela")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m13_spatial_attribution")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "spatial_attribution")
M12_SPECTRUM = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "terminal_sensitivity/m12_error_spectrum_points.csv")
RUNNER = REPO / "build-release/vela_example_runner.exe"
VT_300 = 8.617333262145e-5 * 300.0


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
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def percentile(values: Iterable[float], fraction: float) -> float:
    data = sorted(value for value in values if math.isfinite(value))
    if not data:
        return math.nan
    position = (len(data) - 1) * fraction
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return data[lower]
    weight = position - lower
    return data[lower] * (1.0 - weight) + data[upper] * weight


def stats(values: Iterable[float]) -> dict[str, float | int | None]:
    data = [value for value in values if math.isfinite(value)]
    if not data:
        return {"count": 0, "median": None, "p95": None, "maximum": None}
    return {"count": len(data), "median": percentile(data, 0.5),
            "p95": percentile(data, 0.95), "maximum": max(data)}


def signed_stats(values: Iterable[float]) -> dict[str, float | int | None]:
    data = [value for value in values if math.isfinite(value)]
    result = stats(abs(value) for value in data)
    result["signed_mean"] = sum(data) / len(data) if data else None
    return result


def log_error(reference: float, candidate: float) -> float:
    if reference <= 0.0 or candidate <= 0.0:
        return math.nan
    return abs(math.log10(candidate / reference))


def state_tag(device: str, drain: float, gate: float) -> str:
    drain_tag = "0p05" if drain == 0.05 else "1"
    gate_tags = {0.0: "0", 0.05: "0p05", 0.8: "0p8", 2.5: "2p5"}
    return f"{device}_vd_{drain_tag}_vg_{gate_tags[gate]}"


def drain_tag(drain: float) -> str:
    return "0p05" if drain == 0.05 else "1"


def gate_tag(gate: float) -> str:
    return {0.0: "0", 0.05: "0p05", 0.8: "0p8", 2.5: "2p5"}[gate]


def prepare_projection(case: dict[str, Any], root: Path) -> tuple[Path, Path]:
    state = case["state"]
    source = M10_ROOT / "replay" / state / "vela_drive_edge_mobility.json"
    target_dir = root / "edge_projection" / state
    config_path = target_dir / "edge_projection.json"
    output = target_dir / "edge_projection.csv"
    drive_override = target_dir / "reference_aware_edge_projection_drive.csv"
    sg_rows = read_csv(M10_ROOT / "replay" / state / "vela_drive_sg_edges.csv")
    drive_rows = []
    for row in sg_rows:
        length = float(row["length_m"])
        electron_drive = abs(float(row["phin1_V"]) - float(row["phin0_V"])) / length
        hole_drive = abs(float(row["phip1_V"]) - float(row["phip0_V"])) / length
        drive_rows.append({
            "edge_id": int(row["edge_id"]),
            "electron_drive_V_m": electron_drive,
            "hole_drive_V_m": hole_drive,
        })
    write_csv(drive_override, drive_rows)
    config = read_json(source)
    config["output_csv"] = str(output.resolve())
    config["solver"]["mobility"]["high_field_gradient_discretization"] = "edge_projection"
    config["mobility_drive_override_csv"] = str(drive_override.resolve())
    config["mobility_drive_override_provenance"] = (
        "reference_aware_edge_projection_from_m10_sg_endpoints")
    config["m13"] = {"diagnostic_only": True,
                     "variant": "reference_aware_edge_projection_on_m10_sentaurus_state",
                     "contact_basin_reference_normalized": True}
    write_json(config_path, config)
    return config_path, output


def prepare_vg08(device: str, drain: float, root: Path) -> tuple[Path, Path]:
    tag = drain_tag(drain)
    source_dir = M8_ROOT / device / "workflow" / f"vd_{tag}"
    source = source_dir / "20_gate_sweep.json"
    target_dir = root / "self_consistent" / f"{device}_vd_{tag}_vg_0p8"
    config_path = target_dir / "simulation.json"
    output_state = target_dir / "accepted_state.csv"
    config = read_json(source)
    config["output_csv"] = str((target_dir / "iv.csv").resolve())
    config["log_file"] = str((target_dir / "run.log").resolve())
    config["sweep"]["stop"] = 0.8
    config["sweep"]["bias_points"] = [
        value for value in config["sweep"]["bias_points"] if value <= 0.8 + 1e-12]
    config["sweep"]["write_state_file"] = str(output_state.resolve())
    config["m13"] = {"diagnostic_only": True,
                     "purpose": "exact-lattice self-consistent Vg=0.8 state"}
    write_json(config_path, config)
    return config_path, output_state


def run_config(config: Path, runner: Path) -> dict[str, Any]:
    environment = os.environ.copy()
    environment["PATH"] = (r"D:\msys64\ucrt64\bin" + os.pathsep
                           + environment.get("PATH", ""))
    completed = subprocess.run([str(runner), "--config", str(config)],
                               cwd=REPO, text=True, capture_output=True,
                               env=environment, check=False)
    (config.parent / f"{config.stem}.stdout.txt").write_text(
        completed.stdout, encoding="utf-8")
    (config.parent / f"{config.stem}.stderr.txt").write_text(
        completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(f"Vela M13 run failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    final_line = completed.stdout.strip().splitlines()[-1]
    return json.loads(final_line)


def execute_missing(jobs: list[tuple[Path, Path]], runner: Path,
                    force: bool, workers: int) -> list[dict[str, Any]]:
    selected = [(config, output) for config, output in jobs
                if force or not output.exists()]
    results: list[dict[str, Any]] = []
    if not selected:
        return results
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(run_config, config, runner): (config, output)
                   for config, output in selected}
        for future in as_completed(futures):
            config, output = futures[future]
            result = future.result()
            if not output.exists():
                raise FileNotFoundError(f"M13 run did not create {output}")
            results.append({"config": portable(config), "output": portable(output),
                            "summary": result})
    return results


def scalar_field(export_dir: Path, name: str) -> dict[int, float]:
    path = export_dir / "fields" / f"{name}_region0.csv"
    return {int(row["node_id"]): float(row["component0"])
            for row in read_csv(path)}


def node_coordinates(export_dir: Path) -> dict[int, tuple[float, float]]:
    return {int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
            for row in read_csv(export_dir / "nodes.csv")}


def interface_support(export_dir: Path) -> tuple[set[tuple[int, int]], set[int]]:
    region_edges: dict[str, set[tuple[int, int]]] = {"Silicon_1": set(), "Oxide_1": set()}
    for row in read_csv(export_dir / "elements.csv"):
        region = row["region"]
        if region not in region_edges:
            continue
        nodes = [int(row["node0"]), int(row["node1"]), int(row["node2"])]
        for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
            region_edges[region].add(tuple(sorted((a, b))))
    interface = region_edges["Silicon_1"] & region_edges["Oxide_1"]
    nodes = {node for edge in interface for node in edge}
    return interface, nodes


def state_file(device: str, drain: float, gate: float, root: Path) -> Path:
    tag = drain_tag(drain)
    if gate == 0.0:
        return M8_ROOT / device / "workflow" / f"vd_{tag}" / "10_drain_ramp_accepted_state.csv"
    if gate == 0.05:
        return M8_DEEP / f"{device}_vd_{tag}" / "vg_0p05" / "accepted_state.csv"
    if gate == 0.8:
        return root / "self_consistent" / f"{device}_vd_{tag}_vg_0p8" / "accepted_state.csv"
    if gate == 2.5:
        return M8_ROOT / device / "workflow" / f"vd_{tag}" / "20_gate_sweep_accepted_state.csv"
    raise ValueError(f"unsupported gate voltage {gate}")


def read_state(path: Path) -> dict[int, dict[str, float]]:
    return {int(row["node_id"]): {
                "psi": float(row["psi"]), "phin": float(row["phin"]),
                "phip": float(row["phip"]),
                "electrons_m3": float(row["electrons_m3"]),
                "holes_m3": float(row["holes_m3"])}
            for row in read_csv(path)}


def old_slotboom(total_impurity_cm3: float) -> float:
    if total_impurity_cm3 <= 0.0:
        return 0.0
    x = math.log(total_impurity_cm3 / 1.0e17)
    return max(0.0, 0.009 * (x + math.sqrt(x * x + 0.5)))


def profile_analysis(case: dict[str, Any], root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    state = case["state"]
    device = case["device"]
    drain = float(case["drain_voltage_V"])
    gate = float(case["gate_voltage_V"])
    export_dir = M10_ROOT / "sentaurus_exports" / state
    coordinates = node_coordinates(export_dir)
    _, interface_nodes = interface_support(export_dir)
    channel_nodes = sorted(
        node for node in interface_nodes
        if -0.3 - 1e-12 <= coordinates[node][1] <= 0.3 + 1e-12)
    source_nodes = [node for node in channel_nodes if coordinates[node][1] <= 0.0]
    sent = {
        "psi": scalar_field(export_dir, "ElectrostaticPotential"),
        "phin": scalar_field(export_dir, "eQuasiFermiPotential"),
        "phip": scalar_field(export_dir, "hQuasiFermiPotential"),
        "electrons": scalar_field(export_dir, "eDensity"),
        "holes": scalar_field(export_dir, "hDensity"),
        "bgn": scalar_field(export_dir, "BandgapNarrowing"),
        "donors": scalar_field(export_dir, "DonorConcentration"),
        "acceptors": scalar_field(export_dir, "AcceptorConcentration"),
        "conduction_band": scalar_field(export_dir, "ConductionBandEnergy"),
    }
    vela_path = state_file(device, drain, gate, root)
    vela = read_state(vela_path)
    rows: list[dict[str, Any]] = []
    for node in channel_nodes:
        depth, lateral = coordinates[node]
        sent_bgn = sent["bgn"][node]
        total_impurity = sent["donors"][node] + sent["acceptors"][node]
        vela_bgn = old_slotboom(total_impurity)
        rows.append({
            "state": state, "device": device, "drain_voltage_V": drain,
            "gate_voltage_V": gate, "node_id": node,
            "depth_um": depth, "lateral_um": lateral,
            "sentaurus_psi_V": sent["psi"][node],
            "vela_psi_V": vela[node]["psi"],
            "psi_difference_mV": 1.0e3 * (vela[node]["psi"] - sent["psi"][node]),
            "sentaurus_phin_V": sent["phin"][node],
            "vela_phin_V": vela[node]["phin"],
            "phin_difference_mV": 1.0e3 * (vela[node]["phin"] - sent["phin"][node]),
            "sentaurus_barrier_proxy_V": sent["phin"][node] - sent["psi"][node],
            "vela_barrier_proxy_V": vela[node]["phin"] - vela[node]["psi"],
            "sentaurus_eDensity_cm_minus3": sent["electrons"][node],
            "vela_eDensity_cm_minus3": vela[node]["electrons_m3"] / 1.0e6,
            "electron_density_log10_ratio_dex": math.log10(
                max(vela[node]["electrons_m3"] / 1.0e6, 1e-300)
                / max(sent["electrons"][node], 1e-300)),
            "sentaurus_bgn_eV": sent_bgn,
            "total_ionized_impurity_cm_minus3": total_impurity,
            "vela_old_slotboom_total_impurity_eV": vela_bgn,
            "old_slotboom_difference_meV": 1.0e3 * (vela_bgn - sent_bgn),
            "sentaurus_ni_eff_multiplier": math.exp(sent_bgn / (2.0 * VT_300)),
            "vela_old_slotboom_ni_eff_multiplier": math.exp(vela_bgn / (2.0 * VT_300)),
            "sentaurus_conduction_band_eV": sent["conduction_band"][node],
        })
    by_node = {int(row["node_id"]): row for row in rows}
    sent_barrier_node = max(source_nodes,
                            key=lambda node: float(by_node[node]["sentaurus_barrier_proxy_V"]))
    vela_barrier_node = max(source_nodes,
                            key=lambda node: float(by_node[node]["vela_barrier_proxy_V"]))
    sent_barrier = float(by_node[sent_barrier_node]["sentaurus_barrier_proxy_V"])
    vela_barrier = float(by_node[vela_barrier_node]["vela_barrier_proxy_V"])
    sent_min_n_node = min(source_nodes, key=lambda node: sent["electrons"][node])
    vela_min_n_node = min(source_nodes,
                          key=lambda node: vela[node]["electrons_m3"])
    metrics = {
        "interface_node_count": len(interface_nodes),
        "channel_profile_node_count": len(channel_nodes),
        "vela_state": portable(vela_path),
        "vela_state_sha256": sha256(vela_path),
        "psi_difference_mV": signed_stats(
            float(row["psi_difference_mV"]) for row in rows),
        "phin_difference_mV": signed_stats(
            float(row["phin_difference_mV"]) for row in rows),
        "electron_density_log10_ratio_dex": signed_stats(
            float(row["electron_density_log10_ratio_dex"]) for row in rows),
        "old_slotboom_difference_meV": signed_stats(
            float(row["old_slotboom_difference_meV"]) for row in rows),
        "source_barrier": {
            "sentaurus_node": sent_barrier_node,
            "sentaurus_lateral_um": coordinates[sent_barrier_node][1],
            "sentaurus_proxy_V": sent_barrier,
            "vela_node": vela_barrier_node,
            "vela_lateral_um": coordinates[vela_barrier_node][1],
            "vela_proxy_V": vela_barrier,
            "vela_minus_sentaurus_mV": 1.0e3 * (vela_barrier - sent_barrier),
            "sentaurus_min_density_node": sent_min_n_node,
            "vela_min_density_node": vela_min_n_node,
            "minimum_density_log10_ratio_dex": math.log10(
                max(vela[vela_min_n_node]["electrons_m3"] / 1.0e6, 1e-300)
                / max(sent["electrons"][sent_min_n_node], 1e-300)),
        },
    }
    return metrics, rows


def edge_analysis(case: dict[str, Any], root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    state = case["state"]
    m10_dir = M10_ROOT / "replay" / state
    comparison = read_csv(m10_dir / "edge_replay_comparison.csv")
    projection = {int(row["edge_id"]): row for row in read_csv(
        root / "edge_projection" / state / "edge_projection.csv")}
    sent_drive = {int(row["edge_id"]): row for row in read_csv(
        m10_dir / "sentaurus_mean_endpoint_magnitude_edge_mobility.csv")}
    m11_dir = M11_ROOT / state
    phumob = {int(row["edge_id"]): row for row in read_csv(
        m11_dir / "p1_e0_h0_edge_mobility.csv")}
    enormal = {int(row["edge_id"]): row for row in read_csv(
        m11_dir / "p1_e1_h0_edge_mobility.csv")}
    full = {int(row["edge_id"]): row for row in read_csv(
        m11_dir / "p1_e1_h1_edge_mobility.csv")}
    export_dir = M10_ROOT / "sentaurus_exports" / state
    interface_edges, _ = interface_support(export_dir)
    rows: list[dict[str, Any]] = []
    for source in comparison:
        edge = int(source["edge_id"])
        if not all(edge in mapping for mapping in
                   (projection, sent_drive, phumob, enormal, full)):
            continue
        n0, n1 = int(source["node0"]), int(source["node1"])
        sent_mu = float(source["sentaurus_eMobility_cm2_V_s"])
        # M10's comparison records the mobility used by the production SG
        # contact-current path.  Use it instead of re-evaluating raw external
        # quasi-Fermi values in the standalone probe, because production first
        # normalizes contact-basin references.
        transport_mu = float(source["vela_eMobility_vela_drive_cm2_V_s"])
        projection_mu = float(projection[edge]["electron_final_mobility_m2_V_s"]) * 1e4
        sent_drive_mu = float(sent_drive[edge]["electron_final_mobility_m2_V_s"]) * 1e4
        sent_grad = float(source["sentaurus_eGradQF_mean_magnitude_V_cm"])
        rows.append({
            "state": state, "edge_id": edge, "node0": n0, "node1": n1,
            "x_mid_um": float(source["x_mid_um"]),
            "y_mid_um": float(source["y_mid_um"]),
            "active_current_edge": source["active_current_edge"].lower() == "true",
            "silicon_oxide_interface_edge": tuple(sorted((n0, n1))) in interface_edges,
            "sentaurus_gradqf_V_cm": sent_grad,
            "transport_cell_vector_drive_V_cm": float(
                source["vela_eGradQF_drive_V_cm"]),
            "edge_projection_drive_V_cm": float(
                projection[edge]["electron_mobility_field_V_m"]) / 100.0,
            "sentaurus_external_drive_V_cm": float(
                sent_drive[edge]["electron_mobility_field_V_m"]) / 100.0,
            "sentaurus_final_mobility_cm2_V_s": sent_mu,
            "transport_final_mobility_cm2_V_s": transport_mu,
            "edge_projection_final_mobility_cm2_V_s": projection_mu,
            "sentaurus_drive_final_mobility_cm2_V_s": sent_drive_mu,
            "phumob_only_mobility_cm2_V_s": float(
                phumob[edge]["electron_final_mobility_m2_V_s"]) * 1e4,
            "phumob_enormal_mobility_cm2_V_s": float(
                enormal[edge]["electron_final_mobility_m2_V_s"]) * 1e4,
            "phumob_enormal_hfs_mobility_cm2_V_s": float(
                full[edge]["electron_final_mobility_m2_V_s"]) * 1e4,
            "transport_mobility_error_dex": log_error(sent_mu, transport_mu),
            "edge_projection_mobility_error_dex": log_error(sent_mu, projection_mu),
            "sentaurus_drive_mobility_error_dex": log_error(sent_mu, sent_drive_mu),
            "transport_drive_error_dex": log_error(sent_grad, float(
                source["vela_eGradQF_drive_V_cm"])),
            "edge_projection_drive_error_dex": log_error(sent_grad, float(
                projection[edge]["electron_mobility_field_V_m"]) / 100.0),
        })
    active = [row for row in rows if row["active_current_edge"]]
    interface = [row for row in rows if row["silicon_oxide_interface_edge"]]
    def errors(column: str, selected: list[dict[str, Any]]) -> dict[str, Any]:
        return stats(float(row[column]) for row in selected)
    metrics = {
        "silicon_edge_count": len(rows),
        "active_edge_count": len(active),
        "interface_edge_count": len(interface),
        "active_mobility_error_dex": {
            "transport_cell_vector": errors("transport_mobility_error_dex", active),
            "edge_projection": errors("edge_projection_mobility_error_dex", active),
            "sentaurus_drive": errors("sentaurus_drive_mobility_error_dex", active),
        },
        "active_drive_error_dex": {
            "transport_cell_vector": errors("transport_drive_error_dex", active),
            "edge_projection": errors("edge_projection_drive_error_dex", active),
        },
        "interface_mobility_error_dex": {
            "transport_cell_vector": errors("transport_mobility_error_dex", interface),
            "edge_projection": errors("edge_projection_mobility_error_dex", interface),
            "sentaurus_drive": errors("sentaurus_drive_mobility_error_dex", interface),
        },
        "edge_projection_p95_change_vs_transport_dex": (
            float(errors("edge_projection_mobility_error_dex", active)["p95"])
            - float(errors("transport_mobility_error_dex", active)["p95"])),
        "mobility_stage_active_median_cm2_V_s": {
            "phumob": percentile((float(row["phumob_only_mobility_cm2_V_s"])
                                   for row in active), 0.5),
            "phumob_enormal": percentile((float(row["phumob_enormal_mobility_cm2_V_s"])
                                           for row in active), 0.5),
            "phumob_enormal_hfs": percentile((float(row["phumob_enormal_hfs_mobility_cm2_V_s"])
                                               for row in active), 0.5),
            "sentaurus_drive_hfs": percentile((float(row["sentaurus_drive_final_mobility_cm2_V_s"])
                                                for row in active), 0.5),
            "sentaurus_final": percentile((float(row["sentaurus_final_mobility_cm2_V_s"])
                                            for row in active), 0.5),
        },
    }
    write_csv(root / "edge_comparisons" / f"{state}.csv", rows)
    return metrics, rows


def flatten_state(case: dict[str, Any]) -> dict[str, Any]:
    edge = case["edge"]
    profile = case["profile"]
    return {
        "state": case["state"], "device": case["device"],
        "drain_voltage_V": case["drain_voltage_V"],
        "gate_voltage_V": case["gate_voltage_V"],
        "transport_active_p95_mobility_error_dex": edge[
            "active_mobility_error_dex"]["transport_cell_vector"]["p95"],
        "edge_projection_active_p95_mobility_error_dex": edge[
            "active_mobility_error_dex"]["edge_projection"]["p95"],
        "sentaurus_drive_active_p95_mobility_error_dex": edge[
            "active_mobility_error_dex"]["sentaurus_drive"]["p95"],
        "edge_projection_p95_change_vs_transport_dex": edge[
            "edge_projection_p95_change_vs_transport_dex"],
        "psi_p95_difference_mV": profile["psi_difference_mV"]["p95"],
        "phin_p95_difference_mV": profile["phin_difference_mV"]["p95"],
        "electron_density_p95_error_dex": profile[
            "electron_density_log10_ratio_dex"]["p95"],
        "source_barrier_difference_mV": profile[
            "source_barrier"]["vela_minus_sentaurus_mV"],
        "terminal_signed_log10_error_dex": case["terminal_signed_log10_error_dex"],
        "barrier_predicted_log10_shift_dex": case["barrier_predicted_log10_shift_dex"],
        "barrier_prediction_residual_dex": case["barrier_prediction_residual_dex"],
        "old_slotboom_p95_difference_meV": profile[
            "old_slotboom_difference_meV"]["p95"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--portable-output-dir", type=Path, default=PORTABLE)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    root = args.output_dir.resolve()
    portable_root = args.portable_output_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    portable_root.mkdir(parents=True, exist_ok=True)
    m10_report_path = M10_ROOT / "fixed_state_replay_report.json"
    m10_report = read_json(m10_report_path)
    cases_source = m10_report["cases"]
    spectrum = {
        (row["device"], round(float(row["drain_voltage_V"]), 8),
         round(float(row["gate_voltage_V"]), 8)):
            float(row["signed_log10_ratio_dex"])
        for row in read_csv(M12_SPECTRUM)
    }
    projection_jobs = [prepare_projection(case, root) for case in cases_source]
    vg08_jobs = [prepare_vg08(device, drain, root)
                 for device in ("n17", "n21") for drain in (0.05, 1.0)]
    execution: list[dict[str, Any]] = []
    if args.execute:
        execution = execute_missing(projection_jobs + vg08_jobs,
                                    args.runner.resolve(), args.force,
                                    max(1, args.workers))
    required = [output for _, output in projection_jobs + vg08_jobs]
    missing = [path for path in required if not path.exists()]
    if missing:
        write_json(root / "prepared_manifest.json", {
            "status": "prepared", "missing_outputs": [portable(path) for path in missing],
            "configs": [portable(config) for config, _ in projection_jobs + vg08_jobs]})
        print(json.dumps({"status": "prepared", "missing_outputs": len(missing)}))
        return 0

    cases: list[dict[str, Any]] = []
    profiles: list[dict[str, Any]] = []
    for source in cases_source:
        edge_metrics, _ = edge_analysis(source, root)
        profile_metrics, profile_rows = profile_analysis(source, root)
        profiles.extend(profile_rows)
        key = (source["device"], round(float(source["drain_voltage_V"]), 8),
               round(float(source["gate_voltage_V"]), 8))
        terminal_error = spectrum[key]
        barrier_shift = -1.0e-3 * float(
            profile_metrics["source_barrier"]["vela_minus_sentaurus_mV"]
        ) / (VT_300 * math.log(10.0))
        cases.append({"state": source["state"], "device": source["device"],
                      "drain_voltage_V": source["drain_voltage_V"],
                      "gate_voltage_V": source["gate_voltage_V"],
                      "terminal_signed_log10_error_dex": terminal_error,
                      "barrier_predicted_log10_shift_dex": barrier_shift,
                      "barrier_prediction_residual_dex": terminal_error - barrier_shift,
                      "edge": edge_metrics, "profile": profile_metrics})
    state_rows = [flatten_state(case) for case in cases]
    state_summary = portable_root / "m13_state_summary.csv"
    surface_profiles = portable_root / "m13_surface_profiles.csv"
    write_csv(state_summary, state_rows)
    write_csv(surface_profiles, profiles)
    projection_improves = sum(
        float(case["edge"]["edge_projection_p95_change_vs_transport_dex"]) < 0.0
        for case in cases)
    low_drain_cases = [case for case in cases
                       if float(case["drain_voltage_V"]) == 0.05]
    high_drain_cases = [case for case in cases
                        if float(case["drain_voltage_V"]) == 1.0]
    transport_p95 = [float(case["edge"]["active_mobility_error_dex"]
                           ["transport_cell_vector"]["p95"]) for case in cases]
    sentaurus_drive_p95 = [float(case["edge"]["active_mobility_error_dex"]
                                ["sentaurus_drive"]["p95"]) for case in cases]
    barrier_differences = [float(case["profile"]["source_barrier"]
                                 ["vela_minus_sentaurus_mV"]) for case in cases]
    weak_region = [case for case in cases
                   if case["device"] == "n21"
                   and float(case["drain_voltage_V"]) == 0.05
                   and float(case["gate_voltage_V"]) <= 0.8]
    report = {
        "schema": "vela.simplemos.sdevice.m13_spatial_attribution_report.v1",
        "status": "complete",
        "scope": "SDevice-only spatial attribution on 16 exact frozen/self-consistent states",
        "execution": {"state_count": len(cases),
                      "edge_projection_probe_count": len(projection_jobs),
                      "new_self_consistent_vg0p8_state_count": len(vg08_jobs),
                      "required_local_outputs_verified": len(required),
                      "new_sentaurus_execution": False,
                      "cpp_changed": False,
                      "default_model_changed": False,
                      "runs_executed_this_invocation": len(execution)},
        "findings": {
            "edge_projection_improves_active_p95_state_count": projection_improves,
            "transport_cell_vector_better_or_equal_state_count": len(cases) - projection_improves,
            "edge_projection_low_drain_improves_state_count": sum(
                case["edge"]["edge_projection_p95_change_vs_transport_dex"] < 0.0
                for case in low_drain_cases),
            "edge_projection_high_drain_worsens_state_count": sum(
                case["edge"]["edge_projection_p95_change_vs_transport_dex"] > 0.0
                for case in high_drain_cases),
            "mean_active_p95_mobility_error_dex": {
                "transport_cell_vector": sum(transport_p95) / len(transport_p95),
                "edge_projection": sum(
                    float(case["edge"]["active_mobility_error_dex"]
                          ["edge_projection"]["p95"]) for case in cases) / len(cases),
                "sentaurus_drive": sum(sentaurus_drive_p95) / len(sentaurus_drive_p95),
            },
            "sentaurus_drive_improves_active_p95_state_count": sum(
                sent < transport
                for sent, transport in zip(sentaurus_drive_p95, transport_p95)),
            "maximum_abs_edge_projection_p95_change_dex": max(
                abs(float(case["edge"]["edge_projection_p95_change_vs_transport_dex"]))
                for case in cases),
            "source_barrier_difference_mV": {
                "minimum": min(barrier_differences),
                "maximum": max(barrier_differences),
                "mean": sum(barrier_differences) / len(barrier_differences),
                "same_positive_sign_count": sum(value > 0.0 for value in barrier_differences),
            },
            "high_nwell_low_drain_weak_region_barrier": {
                "representative_state_count": len(weak_region),
                "sign_consistent_with_terminal_error_count": sum(
                    case["terminal_signed_log10_error_dex"]
                    * case["barrier_predicted_log10_shift_dex"] > 0.0
                    for case in weak_region),
                "mean_absolute_prediction_residual_dex": sum(
                    abs(case["barrier_prediction_residual_dex"])
                    for case in weak_region) / len(weak_region),
            },
            "maximum_surface_psi_p95_difference_mV": max(
                float(case["profile"]["psi_difference_mV"]["p95"])
                for case in cases),
            "maximum_surface_phin_p95_difference_mV": max(
                float(case["profile"]["phin_difference_mV"]["p95"])
                for case in cases),
            "maximum_old_slotboom_p95_difference_meV": max(
                float(case["profile"]["old_slotboom_difference_meV"]["p95"])
                for case in cases),
        },
        "cases": cases,
        "interpretation_guards": [
            "The edge_projection comparison changes only Vela HFS drive discretization on the same imported Sentaurus state.",
            "The source-barrier proxy is phin-psi on exact Silicon/Oxide interface nodes; it is diagnostic, not a fitted threshold voltage.",
            "OldSlotboom production replay uses total ionized impurity only; carrier populations do not enter buildEffectiveNodeNi.",
            "The barrier-predicted current shift is a Boltzmann proxy, -Delta(phin-psi)/(Vt ln(10)), not a solved terminal current.",
            "Sentaurus low-field mobility is not separately exported, so PhuMob and Enormal stage values are Vela-only diagnostics."
        ],
        "inputs": {
            "m10_report": portable(m10_report_path),
            "m10_report_sha256": sha256(m10_report_path),
            "m12_evidence": "reference_tcad/simplemos_sentaurus2022/simplemos_m12_terminal_sensitivity_evidence.json",
            "m12_evidence_sha256": sha256(REPO / "reference_tcad/simplemos_sentaurus2022"
                                           / "simplemos_m12_terminal_sensitivity_evidence.json"),
        },
        "artifacts": {
            "state_summary": {"path": portable(state_summary),
                              "sha256": sha256(state_summary)},
            "surface_profiles": {"path": portable(surface_profiles),
                                 "sha256": sha256(surface_profiles)},
        },
    }
    report_path = portable_root / "m13_spatial_attribution_report.json"
    write_json(report_path, report)
    write_json(root / "m13_spatial_attribution_report.json", report)
    print(json.dumps({"status": "complete", "states": len(cases),
                      "projection_improves": projection_improves,
                      "report": portable(report_path)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
