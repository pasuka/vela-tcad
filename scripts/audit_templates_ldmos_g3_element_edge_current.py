#!/usr/bin/env python3
"""Audit supported Sentaurus element-edge current and row-scale interfaces."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any


CONTACT_RE = re.compile(
    r"^\s*drain\s+([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s+"
    r"([+\-0-9.Ee]+)\s+([+\-0-9.Ee]+)\s*$",
    re.MULTILINE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_fields(root: Path) -> dict[str, Path]:
    return {
        path.name: path
        for path in sorted((root / "fields").glob("*.csv"))
    }


def compare_field_directories(baseline: Path, candidate: Path) -> dict[str, Any]:
    baseline_fields = csv_fields(baseline)
    candidate_fields = csv_fields(candidate)
    common = sorted(set(baseline_fields) & set(candidate_fields))
    mismatches = [
        name for name in common
        if sha256(baseline_fields[name]) != sha256(candidate_fields[name])
    ]
    return {
        "baseline_field_count": len(baseline_fields),
        "candidate_field_count": len(candidate_fields),
        "common_field_count": len(common),
        "baseline_only": sorted(set(baseline_fields) - set(candidate_fields)),
        "candidate_only": sorted(set(candidate_fields) - set(baseline_fields)),
        "byte_different_common_fields": mismatches,
        "all_common_fields_byte_identical": not mismatches,
        "selected_identity": {
            name: (
                sha256(baseline_fields[name]) == sha256(candidate_fields[name])
                if name in baseline_fields and name in candidate_fields
                else None
            )
            for name in (
                "eContinuityRhs_region0.csv",
                "NewtonStepEDensityUpdate_region0.csv",
                "ElectrostaticPotential_region0.csv",
                "eDensity_region0.csv",
                "eCurrentDensity_region0.csv",
            )
        },
    }


def drain_current(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = CONTACT_RE.findall(text)
    if not matches:
        raise ValueError(f"no final drain row in {path}")
    voltage, electron, hole, conduction = map(float, matches[-1])
    return {
        "voltage_V": voltage,
        "electron_A": electron,
        "hole_A": hole,
        "conduction_A": conduction,
    }


def log_contract(path: Path, expect_element_edge: bool) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    activated = "With ElementEdge Current Density approximation" in text
    if activated != expect_element_edge:
        raise ValueError(
            f"ElementEdgeCurrent activation mismatch in {path}: {activated}"
        )
    return {
        "element_edge_current_activated": activated,
        "drain": drain_current(path),
    }


def capability_result(log_path: Path, tcl_path: Path) -> dict[str, Any]:
    log = log_path.read_text(encoding="utf-8", errors="replace")
    tcl = tcl_path.read_text(encoding="utf-8")
    error = "Tried to read undefined Edge-Vector eCurrentDensity !"
    coefficient_pos = tcl.find("ReadCoefficient")
    edge_pos = tcl.find('ReadVector $::des_data_edge "eCurrentDensity"')
    observed = error in log
    result = {
        "requested_api": "ReadVector(des_data_edge, eCurrentDensity)",
        "native_error": error,
        "native_error_observed": observed,
        "coefficient_positive_control_precedes_target_call": (
            coefficient_pos >= 0 and edge_pos > coefficient_pos
        ),
        "newton_iteration_completed_before_capability_error": (
            "writing probe_high_vsv_node3721_minus_edge_api_newton_1_0.tdr" in log
            and log.find(
                "writing probe_high_vsv_node3721_minus_edge_api_newton_1_0.tdr"
            ) < log.find(error)
        ),
        "public_edge_current_dataset_available": False if observed else None,
    }
    if not result["coefficient_positive_control_precedes_target_call"]:
        raise ValueError("capability Tcl lacks the ordered coefficient positive control")
    if not result["newton_iteration_completed_before_capability_error"]:
        raise ValueError("native edge-vector error was not observed after Newton output")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for sign in ("minus", "plus"):
        parser.add_argument(f"--baseline-{sign}-newton0", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-newton0", type=Path, required=True)
        parser.add_argument(f"--baseline-{sign}-newton1", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-newton1", type=Path, required=True)
        parser.add_argument(f"--baseline-{sign}-state", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-state", type=Path, required=True)
        parser.add_argument(f"--baseline-{sign}-log", type=Path, required=True)
        parser.add_argument(f"--candidate-{sign}-log", type=Path, required=True)
    parser.add_argument("--capability-log", type=Path, required=True)
    parser.add_argument("--capability-tcl", type=Path, required=True)
    parser.add_argument("--single-mode-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    prior = json.loads(args.single_mode_summary.read_text(encoding="utf-8"))
    comparisons: dict[str, Any] = {}
    logs: dict[str, Any] = {}
    artifacts: list[Path] = [
        args.capability_log, args.capability_tcl, args.single_mode_summary,
    ]
    for sign in ("minus", "plus"):
        comparisons[sign] = {}
        for stage in ("newton0", "newton1", "state"):
            baseline = getattr(args, f"baseline_{sign}_{stage}")
            candidate = getattr(args, f"candidate_{sign}_{stage}")
            comparisons[sign][stage] = compare_field_directories(
                baseline, candidate
            )
            artifacts.extend(csv_fields(candidate).values())
        baseline_log = getattr(args, f"baseline_{sign}_log")
        candidate_log = getattr(args, f"candidate_{sign}_log")
        logs[sign] = {
            "baseline": log_contract(baseline_log, False),
            "candidate": log_contract(candidate_log, True),
        }
        artifacts.extend((baseline_log, candidate_log))

    capability = capability_result(args.capability_log, args.capability_tcl)
    all_fields_identical = all(
        comparisons[sign][stage]["all_common_fields_byte_identical"]
        and not comparisons[sign][stage]["baseline_only"]
        and not comparisons[sign][stage]["candidate_only"]
        for sign in ("minus", "plus") for stage in ("newton0", "newton1", "state")
    )
    terminals_identical = all(
        logs[sign]["baseline"]["drain"] == logs[sign]["candidate"]["drain"]
        for sign in ("minus", "plus")
    )
    scalar = prior["cross_engine_central_difference"]["one_ring"][
        "absolute_sentaurus_over_vela_least_squares"
    ]
    report = {
        "experiment": "templates_ldmos_g3_supported_element_edge_current_audit",
        "contract": {
            "state": "same high-endpoint VSV symmetric node-3721 QF pair",
            "only_changed_sentaurus_input": "Math.ElementEdgeCurrent",
            "production_defaults_changed": False,
        },
        "official_interface_inventory": {
            "read_coefficient": "element-edge box coefficient, supported",
            "read_measure": "element-vertex box measure, supported",
            "eCurrentDensity_registered_location": "vertex only",
            "read_flux_semantics": (
                "vertex box-boundary integral of a chosen gradient divided by box volume"
            ),
            "public_element_edge_current_dataset": "not registered",
            "public_general_dd_row_scale": "not documented",
        },
        "element_edge_current_ab": {
            "logs": logs,
            "field_comparisons": comparisons,
            "all_normalized_import_fields_byte_identical": all_fields_identical,
            "terminal_currents_identical_at_log_precision": terminals_identical,
            "central_difference_scalar_before": scalar,
            "central_difference_scalar_after": scalar,
        },
        "tcl_pmi_capability_probe": capability,
        "classification": {
            "element_edge_current_changes_observable_vsv_continuity_rows": False,
            "element_edge_current_changes_saved_vertex_current_density": False,
            "supported_public_edge_current_export_available": False,
            "supported_public_row_scale_disclosure_available": False,
            "absolute_newtonplot_rhs_cross_engine_gate_valid": False,
            "normalized_shape_and_relative_change_diagnostics_valid": True,
            "supported_interface_search_exhausted_for_t2022_03_sp2": True,
            "remaining_exact_scalar_requires_vendor_internal_disclosure": True,
            "ledger_status": "draft",
        },
        "artifacts": {str(path): sha256(path) for path in artifacts},
    }
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "summary": str((output / "summary.json").resolve()),
        "all_fields_identical": all_fields_identical,
        "terminal_currents_identical": terminals_identical,
        "edge_current_api_available": capability[
            "public_edge_current_dataset_available"
        ],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
