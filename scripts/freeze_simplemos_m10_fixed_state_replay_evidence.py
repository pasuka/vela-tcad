#!/usr/bin/env python3
"""Freeze compact, portable evidence for the SimpleMOS M10 fixed-state replay."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import statistics
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                  / "m10_fixed_state_replay")
DEFAULT_DESTINATION = (REPO / "reference_tcad/simplemos_sentaurus2022"
                       / "fixed_state_replay")
DEFAULT_EVIDENCE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "simplemos_m10_fixed_state_replay_evidence.json")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def aggregate(report: dict[str, Any]) -> dict[str, Any]:
    cases = report["cases"]
    mobility = {}
    for stage in ("vela_drive", "sentaurus_drive",
                  "sentaurus_drive_sensitivity"):
        mobility[stage] = {
            metric: statistics.fmean(float(case["electron_mobility"][stage]
                                                ["active_edges_abs_error_dex"]
                                                [metric]) for case in cases)
            for metric in ("median", "p90", "p95", "maximum")
        }
    terminal = {}
    for stage in ("vela_drive_vela_hfs", "sentaurus_drive_vela_hfs",
                  "sentaurus_drive_vela_hfs_sensitivity",
                  "sentaurus_final_mobility"):
        values = [float(case["terminal_current_replay"][stage]
                             ["absolute_log10_error_dex"]) for case in cases]
        terminal[stage] = {
            "mean_dex": statistics.fmean(values),
            "median_dex": statistics.median(values),
            "maximum_dex": max(values),
            "maximum_state": cases[values.index(max(values))]["state"],
        }
    on_cases = [case for case in cases if float(case["gate_voltage_V"]) >= 0.8]
    on_terminal = {
        stage: {
            "mean_dex": statistics.fmean(float(case["terminal_current_replay"]
                                                   [stage]
                                                   ["absolute_log10_error_dex"])
                                         for case in on_cases),
            "maximum_dex": max(float(case["terminal_current_replay"][stage]
                                      ["absolute_log10_error_dex"])
                               for case in on_cases),
        }
        for stage in terminal
    }
    base_p95 = mobility["vela_drive"]["p95"]
    replay_p95 = mobility["sentaurus_drive"]["p95"]
    base_median = mobility["vela_drive"]["median"]
    replay_median = mobility["sentaurus_drive"]["median"]
    return {
        "case_count": len(cases),
        "mobility_active_edge_mean_statistics_dex": mobility,
        "mobility_improvement": {
            "median_reduction_dex": base_median - replay_median,
            "median_reduction_fraction": 1.0 - replay_median / base_median,
            "p95_reduction_dex": base_p95 - replay_p95,
            "p95_reduction_fraction": 1.0 - replay_p95 / base_p95,
            "sentaurus_drive_improves_all_states_at_p95": all(
                float(case["electron_mobility"]["sentaurus_drive"]
                           ["active_edges_abs_error_dex"]["p95"])
                <= float(case["electron_mobility"]["vela_drive"]
                            ["active_edges_abs_error_dex"]["p95"])
                for case in cases),
        },
        "terminal_current_all_states": terminal,
        "terminal_current_strong_inversion": on_terminal,
        "mapping_sensitivity_max_p95_difference_dex": max(abs(
            float(case["electron_mobility"]["sentaurus_drive"]
                       ["active_edges_abs_error_dex"]["p95"])
            - float(case["electron_mobility"]["sentaurus_drive_sensitivity"]
                    ["active_edges_abs_error_dex"]["p95"])) for case in cases),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination", type=Path, default=DEFAULT_DESTINATION)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    args = parser.parse_args()
    source = args.source.resolve()
    destination = args.destination.resolve()
    destination.mkdir(parents=True, exist_ok=True)

    report = read_json(source / "fixed_state_replay_report.json")
    summary = aggregate(report)
    report_target = destination / "fixed_state_replay_report.json"
    csv_target = destination / "fixed_state_replay_summary.csv"
    aggregate_target = destination / "aggregate_summary.json"
    shutil.copy2(source / "fixed_state_replay_report.json", report_target)
    csv_target.write_text(
        (source / "fixed_state_replay_summary.csv").read_text(
            encoding="utf-8-sig"),
        encoding="utf-8", newline="\n")
    write_json(aggregate_target, summary)

    contract = (REPO / "reference_tcad/simplemos_sentaurus2022"
                / "simplemos_m10_fixed_state_replay_contract_v1.json")
    implementations = [
        contract,
        REPO / "src/tools/vela_example_runner.cpp",
        REPO / "scripts/run_simplemos_m10_fixed_state_replay.py",
        REPO / "scripts/plot_simplemos_m10_fixed_state_replay.py",
        REPO / "scripts/freeze_simplemos_m10_fixed_state_replay_evidence.py",
        REPO / "scripts/build_simplemos_m10_report_artifact.py",
        REPO / "tests/regression/test_reference_tcad_tools.py",
        REPO / "tests/regression/test_simplemos_m10_fixed_state_replay.py",
    ]
    figures = sorted((REPO / "docs/validation/figures/simplemos_m10").glob("*.png"))
    evidence = {
        "schema": "vela.simplemos.sdevice.m10_fixed_state_replay_evidence.v1",
        "status": "complete",
        "scope": "SDevice-only 16-state intermediate-quantity cross replay",
        "sentaurus_release": report["sentaurus_release"],
        "execution": {
            "devices": ["n17", "n21"],
            "drain_voltages_V": [0.05, 1.0],
            "gate_voltages_V": [0.0, 0.05, 0.8, 2.5],
            "state_count": report["case_count"],
            "interpolation": "forbidden",
            "default_model_changed": report["default_model_changed"],
        },
        "baseline_guard": read_json(contract)["default_model_policy"],
        "findings": summary,
        "interpretation": {
            "supported": [
                "Sentaurus GradQuasiFermi drive substitution reduces Vela-to-Sentaurus active-edge electron-mobility residuals in all 16 states at P95",
                "The two documented node-to-edge drive mappings have negligible sensitivity at the aggregate P95 level",
                "Strong-inversion terminal-current replay is nearly identical after the frozen Sentaurus state is mapped into Vela",
            ],
            "not_supported": [
                "Attributing all remaining mobility residual to the HFS formula alone because Sentaurus low-field mobility is not separately exported",
                "Using deep-off-state terminal-current replay for HFS attribution because cancellation and sub-fA current magnify mapping/discretization residuals",
                "Treating projected Sentaurus nodal current as an exact finite-volume edge-current identity",
            ],
        },
        "artifacts": {
            "report": {"path": portable(report_target), "sha256": sha256(report_target)},
            "summary": {"path": portable(csv_target), "sha256": sha256(csv_target)},
            "aggregate": {"path": portable(aggregate_target), "sha256": sha256(aggregate_target)},
            "figures": [{"path": portable(path), "sha256": sha256(path)}
                        for path in figures],
        },
        "source_integrity": {
            "sentaurus_manifest_sha256": sha256(source / "sentaurus_manifest.json"),
            "sentaurus_export_manifest_sha256": sha256(source / "sentaurus_export_manifest.json"),
            "raw_report_sha256": sha256(source / "fixed_state_replay_report.json"),
        },
        "implementation_sha256": {portable(path): sha256(path)
                                  for path in implementations},
    }
    write_json(args.evidence.resolve(), evidence)
    print(json.dumps({"status": "frozen", "states": report["case_count"],
                      "evidence": portable(args.evidence.resolve())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
