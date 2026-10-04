#!/usr/bin/env python3
"""Run the SimpleMOS M32 scaled electron quasi-Fermi perturbation audit."""

from __future__ import annotations

import argparse
import copy
import json
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import run_simplemos_m31_minority_poisson_perturbation as m31


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m32_electron_qf_scaled_perturbation_contract_v1.json")
M29 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m29_bgn_srh_factorial")
M30 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m30_double_off_causal_closure")
M31 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m31_minority_poisson_perturbation")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m32_electron_qf_scaled_perturbation")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "electron_qf_scaled_perturbation")
RUNNER = REPO / "build-release/vela_example_runner.exe"
Q = 1.602176634e-19


def label(fraction: float) -> str:
    return f"f{fraction:.2f}".replace(".", "p")


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def scaled_state(fields: list[str], baseline: dict[int, dict[str, str]],
                 sentaurus: dict[int, dict[str, str]], fraction: float,
                 output: Path) -> None:
    rows = {node: dict(row) for node, row in baseline.items()}
    for node, row in rows.items():
        base = float(baseline[node]["phin"])
        target = float(sentaurus[node]["phin"])
        row["phin"] = str(base + fraction * (target - base))
    m31.write_state(output, fields, rows)


def drain_nodes(mesh: dict[str, Any]) -> set[int]:
    contacts = [contact for contact in mesh["contacts"]
                if contact["name"].lower() == "drain"]
    if len(contacts) != 1:
        raise ValueError("M32 requires exactly one drain contact")
    return {int(node) for node in contacts[0]["node_ids"]}


def cut_rows(path: Path, contact_nodes: set[int]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for row in m31.read_csv(path):
        node0, node1 = int(row["node0"]), int(row["node1"])
        node0_contact = node0 in contact_nodes
        node1_contact = node1 in contact_nodes
        if node0_contact == node1_contact:
            continue
        outward = 1.0 if node0_contact else -1.0
        electron = (-Q * outward
                    * float(row["electron_particle_line_flux_per_m_s"])
                    * 1.0e-6)
        hole = (-Q * outward
                * float(row["hole_particle_line_flux_per_m_s"])
                * 1.0e-6)
        result.append({
            **row,
            "outward_sign": outward,
            "electron_current_A_per_um": electron,
            "hole_current_A_per_um": hole,
            "total_current_A_per_um": electron - hole,
        })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = m31.read_json(CONTRACT)
    fractions = [float(value) for value in contract["perturbation"]["fractions"]]
    base = m31.read_json(M29 / "vela/bgn_off_srh_off/20_gate/config.json")
    baseline_path = M29 / "vela/bgn_off_srh_off/20_gate/state.csv"
    sentaurus_path = M30 / "states/bgn_off_srh_off_sentaurus.csv"
    original_fields, baseline = m31.state_map(baseline_path)
    sentaurus_fields, sentaurus = m31.state_map(sentaurus_path)
    core_fields = ["node_id", "psi", "phin", "phip",
                   "electrons_m3", "holes_m3"]
    required = set(core_fields)
    if (not required.issubset(original_fields)
            or not required.issubset(sentaurus_fields)
            or set(baseline) != set(sentaurus)):
        raise ValueError("M32 baseline and mapped Sentaurus states do not align")

    mesh = m31.read_json(Path(base["mesh_file"]))
    contact_nodes = drain_nodes(mesh)
    sent_manifest = m31.read_json(M29 / "sentaurus_export_manifest.json")
    sent_export = REPO / next(
        row["export_dir"] for row in sent_manifest["states"]
        if row["cell"] == "bgn_off_srh_off")
    silicon_nodes = set(m31.m10.scalar_field(sent_export, "eDensity"))
    free_silicon = silicon_nodes - {
        int(node) for contact in mesh["contacts"] for node in contact["node_ids"]}

    states: dict[float, Path] = {}
    functionals: dict[float, dict[str, Any]] = {}
    sg_paths: dict[float, Path] = {}
    ordinary_jobs: list[tuple[float, str, Path, list[Path]]] = []
    for fraction in fractions:
        name = label(fraction)
        run_dir = OUTPUT / name
        state = run_dir / "state.csv"
        scaled_state(core_fields, baseline, sentaurus, fraction, state)
        states[fraction] = state

        functional_config = run_dir / "functional.json"
        residual_csv = run_dir / "residual.csv"
        deck = m31.probe_deck(base, "terminal_current_functional_probe", state)
        deck.update({"contact": "drain",
                     "residual_output_csv": str(residual_csv.resolve()),
                     "simplemos_m32": {"fraction": fraction, "read_only": True}})
        m31.write_json(functional_config, deck)
        ordinary_jobs.append((fraction, "functional", functional_config,
                              [residual_csv]))

        sg_config = run_dir / "sg_edges.json"
        sg_csv = run_dir / "sg_edges.csv"
        deck = m31.probe_deck(base, "sg_edge_flux_probe", state)
        deck.update({"output_csv": str(sg_csv.resolve()),
                     "simplemos_m32": {"fraction": fraction, "read_only": True}})
        m31.write_json(sg_config, deck)
        ordinary_jobs.append((fraction, "sg", sg_config, [sg_csv]))
        sg_paths[fraction] = sg_csv

    ordinary_status: dict[tuple[float, str], dict[str, Any]] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {
            pool.submit(m31.run_config, config, args.runner.resolve(),
                        args.force if args.run else False, required_paths):
                (fraction, probe)
            for fraction, probe, config, required_paths in ordinary_jobs
        }
        for future in as_completed(futures):
            key = futures[future]
            ordinary_status[key] = future.result()
            print(f"completed {label(key[0])}:{key[1]}", flush=True)
    functionals = {fraction: ordinary_status[(fraction, "functional")]
                   for fraction in fractions}

    cross_results: dict[float, dict[str, Any]] = {}
    cross_jobs: list[tuple[float, Path, list[Path]]] = []
    for fraction in fractions:
        run_dir = OUTPUT / label(fraction) / "cross_block"
        feedback_dir = run_dir / "feedback_fields"
        _, replacement = m31.state_map(states[fraction])
        m31.feedback_fields(feedback_dir, baseline, replacement, True, False)
        config = run_dir / "cross_block.json"
        output_csv = run_dir / "cross_block.csv"
        blocks_csv = run_dir / "jacobian_blocks.csv"
        loop_csv = run_dir / "schur_loop.csv"
        deck = m31.probe_deck(
            base, "newton_poisson_qfp_cross_block_probe", baseline_path)
        deck.update({
            "output_csv": str(output_csv.resolve()),
            "feedback_state_fields_dir": str(feedback_dir.resolve()),
            "jacobian_blocks_csv": str(blocks_csv.resolve()),
            "schur_loop_csv": str(loop_csv.resolve()),
            "compute_condition_estimates": False,
            "simplemos_m32": {"fraction": fraction, "read_only": True},
        })
        m31.write_json(config, deck)
        cross_jobs.append((fraction, config, [output_csv, blocks_csv, loop_csv]))
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = {
            pool.submit(m31.run_config, config, args.runner.resolve(),
                        args.force if args.run else False, required_paths): fraction
            for fraction, config, required_paths in cross_jobs
        }
        for future in as_completed(futures):
            fraction = futures[future]
            cross_results[fraction] = future.result()
            print(f"completed {label(fraction)}:cross", flush=True)

    if not args.analyze:
        return 0

    adjoint = {int(row["node_id"]): row for row in m31.read_csv(
        M31 / "baseline/adjoint.csv")}
    tangent = sum(
        float(adjoint[node]["dI_dphin_A_per_um_per_V"])
        * (float(sentaurus[node]["phin"]) - float(baseline[node]["phin"]))
        for node in free_silicon)
    baseline_current = float(functionals[0.0]["current_A_per_um"])
    baseline_cut = cut_rows(sg_paths[0.0], contact_nodes)
    baseline_edges = {int(row["edge_id"]): row for row in baseline_cut}
    if not baseline_edges:
        raise AssertionError("M32 drain contact cut is empty")

    fraction_rows: list[dict[str, Any]] = []
    all_edge_rows: list[dict[str, Any]] = []
    cut_edge_rows: list[dict[str, Any]] = []
    cross_rows: list[dict[str, Any]] = []
    cross_nodes_by_fraction: dict[float, dict[int, dict[str, str]]] = {}
    baseline_cross_nodes: dict[int, dict[str, str]] | None = None
    for fraction in fractions:
        status = functionals[fraction]
        current = float(status["current_A_per_um"])
        exact_delta = current - baseline_current
        predicted = fraction * tangent
        nonlinear_remainder = exact_delta - predicted
        tangent_error = (abs(nonlinear_remainder)
                         / max(abs(exact_delta), abs(predicted), 1.0e-300))
        cut = cut_rows(sg_paths[fraction], contact_nodes)
        cut_by_edge = {int(row["edge_id"]): row for row in cut}
        if set(cut_by_edge) != set(baseline_edges):
            raise AssertionError("M32 drain-cut support changed across fractions")
        electron_cut = sum(float(row["electron_current_A_per_um"])
                           for row in cut)
        hole_cut = sum(float(row["hole_current_A_per_um"]) for row in cut)
        total_cut = electron_cut - hole_cut
        condition = sum(abs(float(row["total_current_A_per_um"]))
                        for row in cut) / max(abs(total_cut), 1.0e-300)
        cross_status = cross_results[fraction]
        rows = m31.read_csv(
            OUTPUT / label(fraction) / "cross_block/cross_block.csv")
        current_cross_nodes = {int(row["node_id"]): row for row in rows}
        cross_nodes_by_fraction[fraction] = current_cross_nodes
        if baseline_cross_nodes is None:
            baseline_cross_nodes = current_cross_nodes
        induced_psi = [
            float(current_cross_nodes[node]["full_raw_delta_psi_V"])
            - float(baseline_cross_nodes[node]["full_raw_delta_psi_V"])
            for node in current_cross_nodes]
        induced_psi_qfp = [
            float(current_cross_nodes[node]["psi_qfp_product"])
            - float(baseline_cross_nodes[node]["psi_qfp_product"])
            for node in current_cross_nodes]
        induced_phin = [
            float(current_cross_nodes[node]["full_raw_delta_phin_V"])
            - float(baseline_cross_nodes[node]["full_raw_delta_phin_V"])
            for node in current_cross_nodes]
        cross_summary = {
            "target_delta_phin_l2_V": l2([
                float(row["target_delta_phin_V"]) for row in rows]),
            "psi_qfp_product_l2": l2([
                float(row["psi_qfp_product"]) for row in rows]),
            "induced_psi_qfp_product_l2": l2(induced_psi_qfp),
            "full_raw_delta_psi_l2_V": l2([
                float(row["full_raw_delta_psi_V"]) for row in rows]),
            "induced_full_raw_delta_psi_l2_V": l2(induced_psi),
            "full_raw_delta_phin_l2_V": l2([
                float(row["full_raw_delta_phin_V"]) for row in rows]),
            "induced_full_raw_delta_phin_l2_V": l2(induced_phin),
            "schur_relative_closure": cross_status["schur_relative_closure"],
            "condition_estimates_computed": bool(
                cross_status["condition_estimates_computed"]),
        }
        cross_rows.append({"fraction": fraction, **cross_summary})
        fraction_rows.append({
            "fraction": fraction,
            "current_A_per_um": current,
            "exact_delta_current_A_per_um": exact_delta,
            "tangent_predicted_delta_current_A_per_um": predicted,
            "nonlinear_remainder_A_per_um": nonlinear_remainder,
            "tangent_relative_error": tangent_error,
            "absolute_log10_current_shift_dex": abs(math.log10(
                max(abs(current), 1.0e-300)
                / max(abs(baseline_current), 1.0e-300))),
            "psi_residual_norm": status["block_residuals"]["psi"],
            "phin_residual_norm": status["block_residuals"]["phin"],
            "phip_residual_norm": status["block_residuals"]["phip"],
            "drain_cut_electron_A_per_um": electron_cut,
            "drain_cut_hole_A_per_um": hole_cut,
            "drain_cut_total_A_per_um": total_cut,
            "drain_cut_condition": condition,
            "functional_minus_cut_A_per_um": current - total_cut,
            **cross_summary,
        })
        for row in cut:
            edge_id = int(row["edge_id"])
            base_row = baseline_edges[edge_id]
            cut_edge_rows.append({
                "fraction": fraction,
                "edge_id": edge_id,
                "node0": int(row["node0"]),
                "node1": int(row["node1"]),
                "x0_m": float(row["x0"]),
                "y0_m": float(row["y0"]),
                "x1_m": float(row["x1"]),
                "y1_m": float(row["y1"]),
                "electron_current_A_per_um": float(
                    row["electron_current_A_per_um"]),
                "delta_electron_current_A_per_um": float(
                    row["electron_current_A_per_um"])
                    - float(base_row["electron_current_A_per_um"]),
                "hole_current_A_per_um": float(row["hole_current_A_per_um"]),
                "total_current_A_per_um": float(row["total_current_A_per_um"]),
                "phin0_V": float(row["phin0_V"]),
                "phin1_V": float(row["phin1_V"]),
                "electron_sg_cancellation_condition": float(
                    row["electron_sg_cancellation_condition"]),
            })

        sg_rows = m31.read_csv(sg_paths[fraction])
        base_sg = {int(row["edge_id"]): row for row in m31.read_csv(sg_paths[0.0])}
        for row in sg_rows:
            edge_id = int(row["edge_id"])
            base_row = base_sg[edge_id]
            particle_current = (-Q
                * float(row["electron_particle_line_flux_per_m_s"]) * 1.0e-6)
            base_particle_current = (-Q
                * float(base_row["electron_particle_line_flux_per_m_s"]) * 1.0e-6)
            all_edge_rows.append({
                "fraction": fraction,
                "edge_id": edge_id,
                "node0": int(row["node0"]),
                "node1": int(row["node1"]),
                "x_mid_m": 0.5 * (float(row["x0"]) + float(row["x1"])),
                "y_mid_m": 0.5 * (float(row["y0"]) + float(row["y1"])),
                "is_drain_cut": edge_id in baseline_edges,
                "electron_line_current_A_per_um": particle_current,
                "delta_electron_line_current_A_per_um": (
                    particle_current - base_particle_current),
                "electron_mobility_drive_V_m": float(
                    row["electron_mobility_drive_V_m"]),
                "electron_mobility_m2_V_s": float(
                    row["electron_mobility_m2_V_s"]),
                "electron_sg_cancellation_condition": float(
                    row["electron_sg_cancellation_condition"]),
            })

    full_cross = cross_nodes_by_fraction[1.0]
    top_nodes = sorted(
        full_cross,
        key=lambda node: abs(float(full_cross[node]["psi_qfp_product"])),
        reverse=True)[:20]
    node_rows: list[dict[str, Any]] = []
    for fraction in fractions[1:]:
        rows = cross_nodes_by_fraction[fraction]
        for rank, node in enumerate(top_nodes, 1):
            row = rows[node]
            node_rows.append({
                "fraction": fraction,
                "rank_at_full_fraction": rank,
                "node_id": node,
                "x": float(row["x"]),
                "y": float(row["y"]),
                "is_contact": int(row["is_contact"]),
                "target_delta_phin_V": float(row["target_delta_phin_V"]),
                "psi_qfp_product": float(row["psi_qfp_product"]),
                "induced_psi_qfp_product": (
                    float(row["psi_qfp_product"])
                    - float(cross_nodes_by_fraction[0.0][node][
                        "psi_qfp_product"])),
                "full_raw_delta_psi_V": float(row["full_raw_delta_psi_V"]),
                "induced_full_raw_delta_psi_V": (
                    float(row["full_raw_delta_psi_V"])
                    - float(cross_nodes_by_fraction[0.0][node][
                        "full_raw_delta_psi_V"])),
                "full_raw_delta_phin_V": float(row["full_raw_delta_phin_V"]),
                "induced_full_raw_delta_phin_V": (
                    float(row["full_raw_delta_phin_V"])
                    - float(cross_nodes_by_fraction[0.0][node][
                        "full_raw_delta_phin_V"])),
                "qfp_psi_electron_product": float(
                    row["qfp_psi_electron_product"]),
            })

    full_edge_rows = [row for row in all_edge_rows if row["fraction"] == 1.0]
    top_edges = sorted(
        full_edge_rows,
        key=lambda row: abs(row["delta_electron_line_current_A_per_um"]),
        reverse=True)[:20]

    numeric_rows = fraction_rows + cross_rows + cut_edge_rows + node_rows + top_edges
    nonfinite = sum(
        not math.isfinite(float(value))
        for row in numeric_rows for key, value in row.items()
        if key not in ("is_drain_cut",))
    numeric_count = sum(len(row) for row in numeric_rows)
    expected = contract["acceptance"]
    acceptance = {
        "fraction_count": len(fraction_rows),
        "functional_probe_count": len(functionals),
        "sg_probe_count": len(sg_paths),
        "cross_block_probe_count": len(cross_rows),
        "common_silicon_node_count": len(silicon_nodes),
        "drain_cut_edge_count": len(baseline_edges),
        "condition_estimates_skipped_count": sum(
            not row["condition_estimates_computed"] for row in cross_rows),
        "maximum_schur_relative_closure": max(
            row["schur_relative_closure"] for row in cross_rows),
        "nonfinite_fraction": nonfinite / max(numeric_count, 1),
    }
    checks = {
        "fractions": acceptance["fraction_count"] == expected["fraction_count"],
        "functional_probes": acceptance["functional_probe_count"]
            == expected["functional_probe_count"],
        "sg_probes": acceptance["sg_probe_count"] == expected["sg_probe_count"],
        "cross_probes": acceptance["cross_block_probe_count"]
            == expected["cross_block_probe_count"],
        "common_nodes": acceptance["common_silicon_node_count"]
            >= expected["minimum_common_silicon_nodes"],
        "drain_cut": acceptance["drain_cut_edge_count"]
            >= expected["minimum_drain_cut_edges"],
        "schur_closure": acceptance["maximum_schur_relative_closure"]
            <= expected["maximum_schur_relative_closure"],
        "condition_estimates_skipped": (
            not expected["require_condition_estimates_skipped"]
            or acceptance["condition_estimates_skipped_count"]
            == len(cross_rows)),
        "finite": acceptance["nonfinite_fraction"] == 0.0,
    }
    acceptance["checks"] = checks
    acceptance["all_checks_pass"] = all(checks.values())

    m31.write_csv(PORTABLE / "m32_fraction_response.csv", fraction_rows)
    m31.write_csv(PORTABLE / "m32_cross_block_scaling.csv", cross_rows)
    m31.write_csv(PORTABLE / "m32_drain_cut_edge_response.csv", cut_edge_rows)
    m31.write_csv(PORTABLE / "m32_cross_block_node_localization.csv", node_rows)
    m31.write_csv(PORTABLE / "m32_top_transport_edges.csv", top_edges)
    valid = [row["fraction"] for row in fraction_rows
             if row["fraction"] > 0.0 and row["tangent_relative_error"] <= 0.1]
    report = {
        "schema": "vela.simplemos.sdevice.m32_electron_qf_scaled_perturbation_report.v1",
        "status": "complete",
        "execution": {
            "new_sentaurus_execution": False,
            "cpp_changed": False,
            "default_model_changed": False,
        },
        "acceptance": acceptance,
        "linearization": {
            "baseline_current_A_per_um": baseline_current,
            "full_direction_tangent_A_per_um": tangent,
            "ten_percent_relative_error_threshold": 0.1,
            "largest_tested_valid_fraction": max(valid) if valid else None,
        },
        "fraction_response": fraction_rows,
        "cross_block_scaling": cross_rows,
        "top_cross_block_nodes_at_full_fraction": top_nodes,
        "top_transport_edges_at_full_fraction": top_edges,
        "artifacts": {
            "fraction_response": m31.portable(
                PORTABLE / "m32_fraction_response.csv"),
            "cross_block_scaling": m31.portable(
                PORTABLE / "m32_cross_block_scaling.csv"),
            "drain_cut_edge_response": m31.portable(
                PORTABLE / "m32_drain_cut_edge_response.csv"),
            "cross_block_node_localization": m31.portable(
                PORTABLE / "m32_cross_block_node_localization.csv"),
            "top_transport_edges": m31.portable(
                PORTABLE / "m32_top_transport_edges.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    m31.write_json(
        PORTABLE / "m32_electron_qf_scaled_perturbation_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "all_checks_pass": acceptance["all_checks_pass"],
        "linearization": report["linearization"],
        "drain_cut_edge_count": acceptance["drain_cut_edge_count"],
        "top_cross_block_nodes": top_nodes[:10],
        "top_transport_edge_ids": [row["edge_id"] for row in top_edges[:10]],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
