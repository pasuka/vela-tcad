#!/usr/bin/env python3
"""Run the SimpleMOS M31 minority-hole and Poisson-floor perturbation audit."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import run_simplemos_m10_fixed_state_replay as m10


REPO = Path(__file__).resolve().parents[1]
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m31_minority_poisson_perturbation_contract_v1.json")
M29 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m29_bgn_srh_factorial")
M30 = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
       / "m30_double_off_causal_closure")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m31_minority_poisson_perturbation")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "minority_poisson_perturbation")
RUNNER = REPO / "build-release/vela_example_runner.exe"


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


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def run_config(config: Path, runner: Path, force: bool,
               required: list[Path]) -> dict[str, Any]:
    result_path = config.with_suffix(".result.json")
    if not force and result_path.is_file() and all(path.is_file() for path in required):
        return read_json(result_path)
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
        raise RuntimeError(f"M31 probe failed for {config}: "
                           f"{completed.stderr or completed.stdout}")
    status = json.loads(completed.stdout.strip().splitlines()[-1])
    write_json(result_path, status)
    return status


def set_biases(deck: dict[str, Any]) -> None:
    for contact in deck["contacts"]:
        name = contact["name"].lower()
        if name in ("drain", "gate"):
            contact["bias"] = 0.05


def probe_deck(base: dict[str, Any], simulation_type: str,
               state: Path) -> dict[str, Any]:
    deck = copy.deepcopy(base)
    deck.update({
        "simulation_type": simulation_type,
        "state_file": str(state.resolve()),
        "simplemos_m31": {"read_only": True},
    })
    set_biases(deck)
    deck.pop("sweep", None)
    deck.pop("log_file", None)
    deck.pop("output_csv", None)
    return deck


def state_map(path: Path) -> tuple[list[str], dict[int, dict[str, str]]]:
    rows = read_csv(path)
    return list(rows[0]), {int(row["node_id"]): row for row in rows}


def write_state(path: Path, fields: list[str],
                rows: dict[int, dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        for node in sorted(rows):
            writer.writerow(rows[node])


def hybrid_state(fields: list[str], baseline: dict[int, dict[str, str]],
                 replacement: dict[int, dict[str, str]],
                 replace_fields: tuple[str, ...], nodes: set[int] | None,
                 output: Path) -> None:
    result = {node: dict(row) for node, row in baseline.items()}
    chosen = set(result) if nodes is None else nodes
    for node in chosen:
        for field in replace_fields:
            result[node][field] = replacement[node][field]
    write_state(output, fields, result)


def feedback_fields(path: Path, baseline: dict[int, dict[str, str]],
                    replacement: dict[int, dict[str, str]],
                    electron: bool, hole: bool) -> None:
    path.mkdir(parents=True, exist_ok=True)
    specs = {
        "eQuasiFermiPotential": ("phin", electron),
        "hQuasiFermiPotential": ("phip", hole),
        "eDensity_m3": ("electrons_m3", electron),
        "hDensity_m3": ("holes_m3", hole),
    }
    for name, (field, use_replacement) in specs.items():
        output = path / f"{name}_region0.csv"
        with output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(["node_id", "component0"])
            source = replacement if use_replacement else baseline
            for node in sorted(source):
                writer.writerow([node, source[node][field]])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--runner", type=Path, default=RUNNER)
    args = parser.parse_args()
    contract = read_json(CONTRACT)
    base_config_path = M29 / "vela/bgn_off_srh_off/20_gate/config.json"
    base = read_json(base_config_path)
    baseline_path = M29 / "vela/bgn_off_srh_off/20_gate/state.csv"
    sentaurus_path = M30 / "states/bgn_off_srh_off_sentaurus.csv"
    fields, baseline = state_map(baseline_path)
    sent_fields, sentaurus = state_map(sentaurus_path)
    required = {"node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"}
    if (not required.issubset(fields) or not required.issubset(sent_fields)
            or set(baseline) != set(sentaurus)):
        raise ValueError("M31 baseline and Sentaurus mapped-state core schemas differ")
    core_fields = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]

    sent_manifest = read_json(M29 / "sentaurus_export_manifest.json")
    sent_export = REPO / next(row["export_dir"] for row in sent_manifest["states"]
                              if row["cell"] == "bgn_off_srh_off")
    silicon = set(m10.scalar_field(sent_export, "eDensity"))
    mesh = read_json(Path(base["mesh_file"]))
    contacts = {int(node) for contact in mesh["contacts"]
                for node in contact["node_ids"]}
    free_silicon = silicon - contacts
    hole_rank = sorted(
        free_silicon,
        key=lambda node: abs(float(sentaurus[node]["phip"])
                             - float(baseline[node]["phip"])),
        reverse=True)
    hole_top = set(hole_rank[:contract["selection"]["minority_hole_nodes"]])

    variants_dir = OUTPUT / "states"
    variants = {"vela_baseline": baseline_path}
    specifications = {
        "hole_top20_replace": (("phip", "holes_m3"), hole_top),
        "hole_full_replace": (("phip", "holes_m3"), None),
        "electron_full_replace": (("phin", "electrons_m3"), None),
        "qf_full_replace": (("phin", "electrons_m3", "phip", "holes_m3"), None),
        "sentaurus_full_replace": (("psi", "phin", "phip", "electrons_m3", "holes_m3"), None),
    }
    for name, (replace_fields, nodes) in specifications.items():
        path = variants_dir / f"{name}.csv"
        hybrid_state(core_fields, baseline, sentaurus, replace_fields, nodes, path)
        variants[name] = path

    residual_dir = OUTPUT / "baseline"
    residual_dir.mkdir(parents=True, exist_ok=True)
    residual_config = residual_dir / "residual.json"
    residual_csv = residual_dir / "residual.csv"
    deck = probe_deck(base, "newton_residual_probe", baseline_path)
    deck["output_csv"] = str(residual_csv.resolve())
    write_json(residual_config, deck)
    residual_status = run_config(residual_config, args.runner.resolve(),
                                 args.force if args.run else False,
                                 [residual_csv])

    adjoint_config = residual_dir / "adjoint.json"
    adjoint_csv = residual_dir / "adjoint.csv"
    deck = probe_deck(base, "terminal_current_adjoint_probe", baseline_path)
    deck.update({"output_csv": str(adjoint_csv.resolve()), "contact": "drain"})
    write_json(adjoint_config, deck)
    adjoint_status = run_config(adjoint_config, args.runner.resolve(),
                                args.force if args.run else False,
                                [adjoint_csv])

    block_config = residual_dir / "block_step.json"
    block_csv = residual_dir / "block_step.csv"
    deck = probe_deck(base, "newton_block_step_probe", baseline_path)
    deck.update({"output_csv": str(block_csv.resolve()),
                 "block_modes": ["poisson_only", "carrier_only"]})
    write_json(block_config, deck)
    block_status = run_config(block_config, args.runner.resolve(),
                              args.force if args.run else False,
                              [block_csv])

    residual_rows = {int(row["node_id"]): row for row in read_csv(residual_csv)}
    poisson_rank = sorted(
        free_silicon,
        key=lambda node: float(residual_rows[node]["abs_psi_residual"]),
        reverse=True)
    poisson_top = set(poisson_rank[:contract["selection"]["poisson_residual_nodes"]])
    block_rows = [row for row in read_csv(block_csv) if row["mode"] == "poisson_only"]
    block_map = {int(row["node_id"]): row for row in block_rows}
    for name, selected in (("poisson_top10_step", poisson_top),
                           ("poisson_full_step", set(baseline))):
        result = {node: dict(row) for node, row in baseline.items()}
        for node in selected:
            result[node]["psi"] = str(float(result[node]["psi"])
                                      + float(block_map[node]["delta_psi_V"]))
        path = variants_dir / f"{name}.csv"
        write_state(path, core_fields, result)
        variants[name] = path

    functional_results: dict[str, dict[str, Any]] = {}
    for name in contract["state_variants"]:
        run_dir = OUTPUT / "functionals" / name
        run_dir.mkdir(parents=True, exist_ok=True)
        config = run_dir / "functional.json"
        residual_output = run_dir / "residual.csv"
        deck = probe_deck(base, "terminal_current_functional_probe", variants[name])
        deck.update({"contact": "drain",
                     "residual_output_csv": str(residual_output.resolve())})
        write_json(config, deck)
        functional_results[name] = run_config(
            config, args.runner.resolve(), args.force if args.run else False,
            [residual_output])

    cross_results: dict[str, dict[str, Any]] = {}
    cross_jobs: list[tuple[str, Path, list[Path]]] = []
    for name, electron, hole in (
            ("hole_full_replace", False, True),
            ("electron_full_replace", True, False),
            ("qf_full_replace", True, True)):
        run_dir = OUTPUT / "cross_block" / name
        feedback_dir = run_dir / "feedback_fields"
        feedback_fields(feedback_dir, baseline, sentaurus, electron, hole)
        config = run_dir / "cross_block.json"
        output_csv = run_dir / "cross_block.csv"
        blocks_csv = run_dir / "jacobian_blocks.csv"
        loop_csv = run_dir / "schur_loop.csv"
        deck = probe_deck(base, "newton_poisson_qfp_cross_block_probe",
                          baseline_path)
        deck.update({
            "output_csv": str(output_csv.resolve()),
            "feedback_state_fields_dir": str(feedback_dir.resolve()),
            "jacobian_blocks_csv": str(blocks_csv.resolve()),
            "schur_loop_csv": str(loop_csv.resolve()),
            "compute_condition_estimates": False,
        })
        write_json(config, deck)
        cross_jobs.append((name, config, [output_csv, blocks_csv, loop_csv]))
    with ThreadPoolExecutor(max_workers=len(cross_jobs)) as pool:
        futures = {
            name: pool.submit(run_config, config, args.runner.resolve(),
                              args.force if args.run else False, required)
            for name, config, required in cross_jobs
        }
        for name, future in futures.items():
            cross_results[name] = future.result()
            print(f"completed cross-block {name}", flush=True)

    if not args.analyze:
        return 0

    adjoint = {int(row["node_id"]): row for row in read_csv(adjoint_csv)}
    baseline_current = float(functional_results["vela_baseline"]["current_A_per_um"])
    response_rows = []
    for name in contract["state_variants"]:
        _, state = state_map(variants[name])
        status = functional_results[name]
        current = float(status["current_A_per_um"])
        prediction = 0.0
        for node in free_silicon:
            prediction += (
                float(adjoint[node]["dI_dpsi_A_per_um_per_V"])
                * (float(state[node]["psi"]) - float(baseline[node]["psi"]))
                + float(adjoint[node]["dI_dphin_A_per_um_per_V"])
                * (float(state[node]["phin"]) - float(baseline[node]["phin"]))
                + float(adjoint[node]["dI_dphip_A_per_um_per_V"])
                * (float(state[node]["phip"]) - float(baseline[node]["phip"])))
        variant_residual = {int(row["node_id"]): row for row in read_csv(
            OUTPUT / "functionals" / name / "residual.csv")}
        response_rows.append({
            "variant": name,
            "current_A_per_um": current,
            "exact_delta_current_A_per_um": current - baseline_current,
            "absolute_log10_current_shift_dex": abs(math.log10(
                max(abs(current), 1e-300) / max(abs(baseline_current), 1e-300))),
            "adjoint_predicted_delta_current_A_per_um": prediction,
            "prediction_error_A_per_um": prediction - (current - baseline_current),
            "psi_residual_norm": status["block_residuals"]["psi"],
            "phin_residual_norm": status["block_residuals"]["phin"],
            "phip_residual_norm": status["block_residuals"]["phip"],
            "top10_psi_residual_l2": l2([
                float(variant_residual[node]["psi_residual"])
                for node in poisson_top]),
        })
    write_csv(PORTABLE / "m31_state_variant_response.csv", response_rows)

    all_psi_l2 = l2([float(residual_rows[node]["psi_residual"])
                     for node in free_silicon])
    top_psi_l2 = l2([float(residual_rows[node]["psi_residual"])
                     for node in poisson_top])
    node_rows = []
    for rank, node in enumerate(poisson_rank[:10], 1):
        row = residual_rows[node]
        delta = float(block_map[node]["delta_psi_V"])
        dcurrent = float(adjoint[node]["dI_dpsi_A_per_um_per_V"])
        node_rows.append({
            "rank": rank,
            "node_id": node,
            "x": float(row["x"]),
            "y": float(row["y"]),
            "psi_residual": float(row["psi_residual"]),
            "abs_psi_residual": float(row["abs_psi_residual"]),
            "poisson_only_delta_psi_V": delta,
            "dI_dpsi_A_per_um_per_V": dcurrent,
            "adjoint_current_contribution_A_per_um": dcurrent * delta,
            "phip_difference_V": float(sentaurus[node]["phip"])
                - float(baseline[node]["phip"]),
        })
    write_csv(PORTABLE / "m31_poisson_node_audit.csv", node_rows)

    cross_rows = []
    for name, status in cross_results.items():
        rows = read_csv(OUTPUT / "cross_block" / name / "cross_block.csv")
        cross_rows.append({
            "variant": name,
            "condition_estimates_computed": bool(
                status["condition_estimates_computed"]),
            "target_delta_phin_l2_V": l2([float(row["target_delta_phin_V"])
                                           for row in rows]),
            "target_delta_phip_l2_V": l2([float(row["target_delta_phip_V"])
                                           for row in rows]),
            "psi_qfp_product_l2": l2([float(row["psi_qfp_product"])
                                       for row in rows]),
            "full_raw_delta_psi_l2_V": l2([float(row["full_raw_delta_psi_V"])
                                            for row in rows]),
            "full_raw_delta_phin_l2_V": l2([float(row["full_raw_delta_phin_V"])
                                             for row in rows]),
            "full_raw_delta_phip_l2_V": l2([float(row["full_raw_delta_phip_V"])
                                             for row in rows]),
            "schur_relative_closure": status["schur_relative_closure"],
            "j_psi_qfp_fd_relative_error": status[
                "directional_derivative_check"]["J_psi_qfp_relative_error"],
            "j_qfp_psi_fd_relative_error": status[
                "directional_derivative_check"]["J_qfp_psi_relative_error"],
            "transport_boundary_loop_norm": status["loop_component_norms"][
                "transport_boundary"]["C_Ainv_B_norm"],
            "srh_auger_loop_norm": status["loop_component_norms"][
                "srh_auger"]["C_Ainv_B_norm"],
            "sg_avalanche_loop_norm": status["loop_component_norms"][
                "sg_avalanche"]["C_Ainv_B_norm"],
        })
    write_csv(PORTABLE / "m31_cross_block_summary.csv", cross_rows)

    numeric = [float(value) for row in response_rows + node_rows + cross_rows
               for key, value in row.items()
               if key not in ("variant", "condition_estimates_computed")]
    acceptance = {
        "state_variant_count": len(response_rows),
        "cross_block_variant_count": len(cross_rows),
        "common_silicon_node_count": len(silicon),
        "free_silicon_node_count": len(free_silicon),
        "poisson_top10_l2_share": top_psi_l2 / max(all_psi_l2, 1e-300),
        "adjoint_relative_residual": adjoint_status["adjoint_relative_residual"],
        "condition_estimates_skipped_count": sum(
            not row["condition_estimates_computed"] for row in cross_rows),
        "nonfinite_fraction": sum(not math.isfinite(value) for value in numeric)
            / len(numeric),
    }
    expected = contract["acceptance"]
    checks = {
        "state_variants": acceptance["state_variant_count"]
            == expected["state_variant_count"],
        "cross_block_variants": acceptance["cross_block_variant_count"]
            == expected["cross_block_variant_count"],
        "common_nodes": acceptance["common_silicon_node_count"]
            >= expected["minimum_common_silicon_nodes"],
        "poisson_concentration": acceptance["poisson_top10_l2_share"]
            >= expected["minimum_poisson_top10_l2_share"],
        "adjoint": acceptance["adjoint_relative_residual"]
            <= expected["maximum_adjoint_relative_residual"],
        "condition_estimates_skipped": (
            not expected["require_condition_estimates_skipped"]
            or acceptance["condition_estimates_skipped_count"]
            == len(cross_rows)),
        "finite": acceptance["nonfinite_fraction"] == 0.0,
    }
    acceptance["checks"] = checks
    acceptance["all_checks_pass"] = all(checks.values())
    by_variant = {row["variant"]: row for row in response_rows}
    report = {
        "schema": "vela.simplemos.sdevice.m31_minority_poisson_perturbation_report.v1",
        "status": "complete",
        "execution": {
            "new_sentaurus_execution": False,
            "cpp_changed": True,
            "diagnostic_only_cpp_changed": True,
            "default_model_changed": False,
        },
        "acceptance": acceptance,
        "minority_hole": {
            "selected_node_count": len(hole_top),
            "top20_exact_delta_current_A_per_um": by_variant[
                "hole_top20_replace"]["exact_delta_current_A_per_um"],
            "full_exact_delta_current_A_per_um": by_variant[
                "hole_full_replace"]["exact_delta_current_A_per_um"],
            "full_current_shift_dex": by_variant[
                "hole_full_replace"]["absolute_log10_current_shift_dex"],
            "full_psi_residual_norm": by_variant[
                "hole_full_replace"]["psi_residual_norm"],
        },
        "poisson_floor": {
            "top_node_ids": poisson_rank[:10],
            "top10_l2_share": acceptance["poisson_top10_l2_share"],
            "baseline_psi_residual_norm": residual_status[
                "block_residuals"]["psi"],
            "top10_step_psi_residual_norm": by_variant[
                "poisson_top10_step"]["psi_residual_norm"],
            "full_step_psi_residual_norm": by_variant[
                "poisson_full_step"]["psi_residual_norm"],
            "top10_step_current_shift_dex": by_variant[
                "poisson_top10_step"]["absolute_log10_current_shift_dex"],
            "full_step_current_shift_dex": by_variant[
                "poisson_full_step"]["absolute_log10_current_shift_dex"],
        },
        "state_variants": response_rows,
        "cross_block": cross_rows,
        "block_step_status": block_status,
        "artifacts": {
            "state_variant_response": portable(
                PORTABLE / "m31_state_variant_response.csv"),
            "poisson_node_audit": portable(
                PORTABLE / "m31_poisson_node_audit.csv"),
            "cross_block_summary": portable(
                PORTABLE / "m31_cross_block_summary.csv"),
        },
        "claim_policy": contract["claim_policy"],
    }
    write_json(PORTABLE / "m31_minority_poisson_perturbation_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "all_checks_pass": acceptance["all_checks_pass"],
        "minority_hole": report["minority_hole"],
        "poisson_floor": report["poisson_floor"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
