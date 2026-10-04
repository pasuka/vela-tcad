#!/usr/bin/env python3
"""Run SimpleMOS M41 equal-ni SG numerical-path ablation."""

from __future__ import annotations

import argparse
import copy
import math
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import run_simplemos_m10_fixed_state_replay as m10
import run_simplemos_m40_true_no_bgn_factorial as m40


REPO = Path(__file__).resolve().parents[1]
RUNNER = REPO / "build-release/vela_example_runner.exe"
CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m41_equal_ni_flux_ablation_contract_v1.json")
M40_OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
              / "m40_true_no_bgn_factorial")
M40_REPORT = (REPO / "reference_tcad/simplemos_sentaurus2022"
              / "true_no_bgn_factorial/m40_true_no_bgn_factorial_report.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m41_equal_ni_flux_ablation")
PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "equal_ni_flux_ablation")
CELLS = {
    "bgn_off_srh_on": True,
    "bgn_off_srh_off": False,
}
MODES = {
    "legacy": "legacy_factor_difference",
    "compensated": "compensated_log_expm1",
}
Q = 1.602176634e-19


def configure(config: dict[str, Any], srh: bool, evaluation: str) -> None:
    m40.configure_vela(config, "none", srh)
    config["solver"]["bandgap_narrowing"][
        "equal_ni_flux_evaluation"] = evaluation


def workflow(cell: str, srh: bool, mode: str, force: bool) -> dict[str, Any]:
    root = OUTPUT / "self_consistent" / cell / mode
    final_curve = root / "20_gate/curve.csv"
    final_state = root / "20_gate/state.csv"
    final_config = root / "20_gate/config.json"
    if not force and final_curve.is_file() and final_state.is_file():
        rows = m40.read_csv(final_curve)
        row = min(rows, key=lambda item: abs(float(item["bias_V"]) - 0.05))
        return {
            "cell": cell,
            "srh": srh,
            "mode": mode,
            "equal_ni_flux_evaluation": MODES[mode],
            "current_A_per_um": abs(float(row["current_total_A_per_um"])),
            "state": m40.portable(final_state),
            "config": m40.portable(final_config),
            "converged": True,
            "reused": True,
        }
    if force and root.exists():
        shutil.rmtree(root)
    previous: Path | None = None
    phases = (("00_equilibrium", "00_equilibrium.json"),
              ("10_drain", "10_drain_ramp.json"),
              ("20_gate", "20_gate_sweep.json"))
    statuses: dict[str, Any] = {}
    for label, source_name in phases:
        phase = root / label
        phase.mkdir(parents=True, exist_ok=True)
        config = copy.deepcopy(m40.read_json(m40.VELA_BASE / source_name))
        configure(config, srh, MODES[mode])
        curve = phase / "curve.csv"
        state = phase / "state.csv"
        config["output_csv"] = str(curve.resolve())
        config["log_file"] = str((phase / "run.log").resolve())
        config["sweep"]["write_state_file"] = str(state.resolve())
        if previous is not None:
            config["sweep"]["initial_state_file"] = str(previous.resolve())
        if label == "20_gate":
            config["sweep"].update({
                "start": 0.0, "stop": 0.05, "step": 0.05,
                "bias_points": [0.0, 0.05], "max_step": 0.05,
            })
        config["simplemos_m41"] = {
            "cell": cell,
            "srh": srh,
            "equal_ni_flux_evaluation": MODES[mode],
            "default_model_changed": False,
        }
        config_path = phase / "config.json"
        m40.write_json(config_path, config)
        status = m10.execute_runner(config_path, RUNNER)
        if not status.get("converged"):
            raise RuntimeError(f"M41 workflow failed: {cell} {mode} {label}")
        statuses[label] = status
        previous = state
    rows = m40.read_csv(final_curve)
    row = min(rows, key=lambda item: abs(float(item["bias_V"]) - 0.05))
    return {
        "cell": cell,
        "srh": srh,
        "mode": mode,
        "equal_ni_flux_evaluation": MODES[mode],
        "current_A_per_um": abs(float(row["current_total_A_per_um"])),
        "state": m40.portable(root / "20_gate/state.csv"),
        "config": m40.portable(root / "20_gate/config.json"),
        "converged": all(value.get("converged") for value in statuses.values()),
    }


def drain_nodes(mesh_path: Path) -> set[int]:
    mesh = m40.read_json(mesh_path)
    for contact in mesh["contacts"]:
        if contact["name"] == "drain":
            return {int(node) for node in contact["node_ids"]}
    raise ValueError("drain contact not found")


def frozen_probe(cell: str, mode: str) -> tuple[dict[str, Any], list[dict[str, str]]]:
    source_root = M40_OUTPUT / "vela" / cell / "20_gate"
    source_config = source_root / "config.json"
    source_state = source_root / "state.csv"
    root = OUTPUT / "frozen" / cell / mode
    root.mkdir(parents=True, exist_ok=True)
    config = m10.probe_config(
        source_config, source_state, root / "edges.csv", 0.05, 0.05,
        "sg_edge_flux_probe")
    config["solver"]["bandgap_narrowing"][
        "equal_ni_flux_evaluation"] = MODES[mode]
    config["simplemos_m41"] = {
        "cell": cell,
        "frozen_state": m40.portable(source_state),
        "equal_ni_flux_evaluation": MODES[mode],
    }
    config_path = root / "config.json"
    m40.write_json(config_path, config)
    status = m10.execute_runner(config_path, RUNNER)
    rows = m40.read_csv(root / "edges.csv")
    nodes = drain_nodes(Path(config["mesh_file"]))
    terminal = m10.drain_cut_current(rows, nodes)
    return ({
        "cell": cell,
        "mode": mode,
        "equal_ni_flux_evaluation": MODES[mode],
        "edge_count": len(rows),
        "converged": bool(status.get("converged")),
        "drain_cut_total_A_per_um": terminal["total_A_per_um"],
        "drain_cut_electron_A_per_um": terminal["electron_A_per_um"],
        "config": m40.portable(config_path),
        "edges": m40.portable(root / "edges.csv"),
    }, rows)


def finite_stats(values: list[float]) -> dict[str, float | int]:
    data = sorted(value for value in values if math.isfinite(value))
    return {
        "count": len(data),
        "median": m10.percentile(data, 0.5),
        "p95": m10.percentile(data, 0.95),
        "maximum": max(data, default=math.nan),
    }


def analyze_frozen(cell: str, legacy: list[dict[str, str]],
                   compensated: list[dict[str, str]]) -> dict[str, Any]:
    legacy_by_edge = {int(row["edge_id"]): row for row in legacy}
    compensated_by_edge = {int(row["edge_id"]): row for row in compensated}
    ledger: list[dict[str, Any]] = []
    legacy_errors: list[float] = []
    compensated_errors: list[float] = []
    for edge in sorted(set(legacy_by_edge) & set(compensated_by_edge)):
        old = legacy_by_edge[edge]
        new = compensated_by_edge[edge]
        reference = float(new["electron_sg_high_precision_reference_flux"])
        legacy_flux = float(old["electron_flux"])
        compensated_flux = float(new["electron_flux"])
        scale = max(abs(reference),
                    float(new["electron_sg_high_precision_reference_term_scale"]),
                    1.0e-300)
        legacy_error = abs(legacy_flux - reference) / scale
        compensated_error = abs(compensated_flux - reference) / scale
        legacy_errors.append(legacy_error)
        compensated_errors.append(compensated_error)
        ledger.append({
            "edge_id": edge,
            "node0": old["node0"],
            "node1": old["node1"],
            "legacy_electron_flux": legacy_flux,
            "compensated_electron_flux": compensated_flux,
            "high_precision_reference_flux": reference,
            "legacy_normalized_error": legacy_error,
            "compensated_normalized_error": compensated_error,
            "cancellation_condition": float(
                new["electron_sg_cancellation_condition"]),
            "delta_phin_over_vt": abs(
                float(new["electron_sg_log_left_over_right"])),
        })
    m40.write_csv(PORTABLE / f"m41_{cell}_frozen_edge_ledger.csv", ledger)
    return {
        "cell": cell,
        "common_edge_count": len(ledger),
        "legacy_error": finite_stats(legacy_errors),
        "compensated_error": finite_stats(compensated_errors),
        "edges_delta_phin_over_vt_below_1e_minus_12": sum(
            row["delta_phin_over_vt"] < 1.0e-12 for row in ledger),
        "maximum_cancellation_condition": max(
            (row["cancellation_condition"] for row in ledger), default=0.0),
        "ledger": m40.portable(
            PORTABLE / f"m41_{cell}_frozen_edge_ledger.csv"),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    if not RUNNER.is_file() or not M40_REPORT.is_file() or not CONTRACT.is_file():
        raise FileNotFoundError("M41 requires the runner and frozen M40 evidence")
    contract = m40.read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m41_equal_ni_flux_ablation.v1":
        raise ValueError("unexpected M41 contract schema")

    probes: list[dict[str, Any]] = []
    probe_rows: dict[tuple[str, str], list[dict[str, str]]] = {}
    for cell in CELLS:
        for mode in MODES:
            summary, rows = frozen_probe(cell, mode)
            probes.append(summary)
            probe_rows[(cell, mode)] = rows
    frozen = [analyze_frozen(
        cell, probe_rows[(cell, "legacy")],
        probe_rows[(cell, "compensated")]) for cell in CELLS]

    jobs = [(cell, srh, mode, args.force)
            for cell, srh in CELLS.items() for mode in MODES]
    with ThreadPoolExecutor(max_workers=min(args.jobs, len(jobs))) as pool:
        workflows = list(pool.map(lambda value: workflow(*value), jobs))

    m40_report = m40.read_json(M40_REPORT)
    sentaurus = {row["cell"]: float(row["sentaurus_current_A_per_um"])
                 for row in m40_report["matrix"] if row["cell"] in CELLS}
    m40_vela = {row["cell"]: float(row["vela_current_A_per_um"])
                for row in m40_report["matrix"] if row["cell"] in CELLS}
    probe_by_cell_mode = {(row["cell"], row["mode"]): row for row in probes}
    frozen_cut_reconciliation = []
    for cell in CELLS:
        legacy_cut = abs(probe_by_cell_mode[(cell, "legacy")][
            "drain_cut_total_A_per_um"])
        compensated_cut = abs(probe_by_cell_mode[(cell, "compensated")][
            "drain_cut_total_A_per_um"])
        frozen_cut_reconciliation.append({
            "cell": cell,
            "m40_reported_vela_current_A_per_um": m40_vela[cell],
            "legacy_continuity_cut_A_per_um": legacy_cut,
            "compensated_contact_semantics_cut_A_per_um": compensated_cut,
            "legacy_cut_vs_reported_dex": abs(math.log10(
                legacy_cut / m40_vela[cell])),
            "compensated_cut_vs_reported_dex": abs(math.log10(
                compensated_cut / m40_vela[cell])),
        })
    self_consistent = []
    by_cell_mode = {(row["cell"], row["mode"]): row for row in workflows}
    for cell in CELLS:
        legacy = by_cell_mode[(cell, "legacy")]
        compensated = by_cell_mode[(cell, "compensated")]
        reference = sentaurus[cell]
        legacy_gap = abs(math.log10(legacy["current_A_per_um"] / reference))
        compensated_gap = abs(math.log10(
            compensated["current_A_per_um"] / reference))
        self_consistent.append({
            "cell": cell,
            "srh": CELLS[cell],
            "sentaurus_current_A_per_um": reference,
            "legacy_current_A_per_um": legacy["current_A_per_um"],
            "compensated_current_A_per_um": compensated["current_A_per_um"],
            "legacy_gap_dex": legacy_gap,
            "compensated_gap_dex": compensated_gap,
            "gap_improvement_dex": legacy_gap - compensated_gap,
            "compensated_log_current_shift_dex": math.log10(
                compensated["current_A_per_um"] /
                legacy["current_A_per_um"]),
        })
    m40.write_csv(PORTABLE / "m41_self_consistent_response.csv", self_consistent)

    off = next(row for row in self_consistent if not row["srh"])
    on = next(row for row in self_consistent if row["srh"])
    checks = {
        "frozen_probe_count": len(probes) == 4,
        "self_consistent_workflow_count": len(workflows) == 4,
        "all_workflows_converged": all(row["converged"] for row in workflows),
        "compensated_frozen_cut_replays_m40_terminal": max(
            row["compensated_cut_vs_reported_dex"]
            for row in frozen_cut_reconciliation) < 0.01,
        "historical_legacy_cut_differs_from_reported_terminal": min(
            row["legacy_cut_vs_reported_dex"]
            for row in frozen_cut_reconciliation) > 1.0,
        "legacy_m40_reproduced": max(
            abs(math.log10(by_cell_mode[(cell, "legacy")]["current_A_per_um"] /
                           next(row["vela_current_A_per_um"]
                                for row in m40_report["matrix"]
                                if row["cell"] == cell))) for cell in CELLS) < 1e-8,
    }
    report = {
        "schema": "vela.simplemos.sdevice.m41_equal_ni_flux_ablation_report.v1",
        "status": "complete" if all(checks.values()) else "failed",
        "execution": {
            "default_model_changed": False,
            "new_sentaurus_execution": False,
            "frozen_state_count": 2,
            "self_consistent_workflow_count": len(workflows),
            "contract": m40.portable(CONTRACT),
            "contract_sha256": m10.sha256(CONTRACT),
        },
        "numerical_contract": {
            "legacy": MODES["legacy"],
            "compensated": MODES["compensated"],
            "physical_equivalence": (
                "For ni0 == ni1, the VariableNi logarithmic density-gradient "
                "term is zero; only floating-point evaluation changes."),
            "historical_mixed_semantics": (
                "M40 continuity residuals used legacy factor subtraction "
                "while ContactCurrent used the compensated VariableNi kernel."),
        },
        "frozen_probes": probes,
        "frozen_cut_reconciliation": frozen_cut_reconciliation,
        "frozen_edge_analysis": frozen,
        "self_consistent_response": self_consistent,
        "causal_result": {
            "srh_off_legacy_gap_dex": off["legacy_gap_dex"],
            "srh_off_compensated_gap_dex": off["compensated_gap_dex"],
            "srh_off_gap_improvement_dex": off["gap_improvement_dex"],
            "srh_on_legacy_gap_dex": on["legacy_gap_dex"],
            "srh_on_compensated_gap_dex": on["compensated_gap_dex"],
            "srh_on_gap_improvement_dex": on["gap_improvement_dex"],
            "hypothesis_supported": (
                off["gap_improvement_dex"] > 0.5 and
                abs(on["compensated_log_current_shift_dex"]) < 0.01),
            "m40_no_bgn_factorial_clean": False,
            "m40_srh_on_close_parity_survives_consistent_flux": (
                on["compensated_gap_dex"] < 0.05),
            "srh_off_anomaly_partially_closed": (
                off["gap_improvement_dex"] > 0.5),
        },
        "acceptance": {**checks, "all_checks_pass": all(checks.values())},
        "artifacts": {
            "self_consistent_response": m40.portable(
                PORTABLE / "m41_self_consistent_response.csv"),
        },
        "claim_policy": [
            "M41 changes no production default.",
            "The frozen comparison changes only equal-ni floating-point evaluation.",
            "A causal root claim requires both frozen edge evidence and self-consistent closure.",
            "M40 no-BGN causal fractions are superseded because continuity and terminal extraction used different equal-ni numerical paths.",
        ],
    }
    m40.write_json(PORTABLE / "m41_equal_ni_flux_ablation_report.json", report)
    artifact_paths = [
        PORTABLE / "m41_equal_ni_flux_ablation_report.json",
        PORTABLE / "m41_self_consistent_response.csv",
        PORTABLE / "m41_bgn_off_srh_on_frozen_edge_ledger.csv",
        PORTABLE / "m41_bgn_off_srh_off_frozen_edge_ledger.csv",
    ]
    source_paths = [
        CONTRACT,
        Path(__file__).resolve(),
        M40_REPORT,
        REPO / "include/vela/physics/BandgapNarrowing.h",
        REPO / "include/vela/equation/CoupledDDAssembler.h",
        REPO / "src/equation/CoupledDDAssembler.cpp",
        REPO / "src/solver/NewtonSolver.cpp",
        REPO / "src/solver/GummelSolver.cpp",
    ]
    evidence = {
        "schema": "vela.simplemos.sdevice.m41_equal_ni_flux_ablation_evidence.v1",
        "status": "frozen" if report["status"] == "complete" else "failed",
        "artifacts": [
            {"path": m40.portable(path), "sha256": m10.sha256(path)}
            for path in artifact_paths
        ],
        "source_hashes": {
            m40.portable(path): m10.sha256(path) for path in source_paths
        },
        "default_model_changed": False,
        "acceptance": report["acceptance"],
    }
    m40.write_json(REPO / "reference_tcad/simplemos_sentaurus2022"
                   / "simplemos_m41_equal_ni_flux_ablation_evidence.json",
                   evidence)
    print(f"wrote {PORTABLE / 'm41_equal_ni_flux_ablation_report.json'}")
    print(report["causal_result"])
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
