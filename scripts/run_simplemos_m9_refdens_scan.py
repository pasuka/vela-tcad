#!/usr/bin/env python3
"""Run the frozen SimpleMOS M9 HFS reference-density scan."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m8a_confirmation as confirmation  # noqa: E402
import run_simplemos_m8b_hfs_controls as m8b  # noqa: E402


DEFAULT_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m9_hfs_diagnostics_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m9_hfs_diagnostics"
)
DEFAULT_BASELINE = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_confirmation"
)
DEFAULT_RAW = REPO / "build-release/m9_hfs_raw"
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m9_hfs_diagnostics_20260828_v1"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n")


def validate_contract(contract: dict[str, Any]) -> None:
    if contract.get("schema") != "vela.simplemos.sdevice.m9_hfs_diagnostics.v1":
        raise ValueError("unexpected M9 contract schema")
    if contract["default_model_policy"]["modify_vela_default_model"]:
        raise ValueError("M9 must not change the Vela default mobility model")
    if [item["id"] for item in contract["devices"]] != ["n17", "n21"]:
        raise ValueError("M9 diagnostic devices must be n17 and n21")
    expected_refdens = [0.0, 1e2, 1e4, 1e6, 1e8, 1e10]
    actual_refdens = [float(item["refdens_cm3"])
                      for item in contract["variants"]]
    if actual_refdens != expected_refdens:
        raise ValueError("M9 RefDens scan must use the frozen logarithmic grid")
    if contract["variants"][0]["id"] != "explicit_gradqf":
        raise ValueError("M9 first variant must be the explicit GradQF null")
    lattice = [float(value) for value in contract["bias_matrix"][
        "gate_lattice"]["values_V"]]
    if len(lattice) != 51 or any(not math.isclose(
            value, 0.05 * index, abs_tol=1e-12)
            for index, value in enumerate(lattice)):
        raise ValueError("M9 gate lattice must be 0:0.05:2.5 V")
    if contract["comparison"]["interpolation"] != "forbidden":
        raise ValueError("M9 forbids interpolation")


def validate_baseline_guard(contract: dict[str, Any]) -> dict[str, Any]:
    policy = contract["default_model_policy"]
    evidence_path = REPO / policy["m8_baseline_evidence"]
    report_path = REPO / policy["m8_comparison_report"]
    if m8b.sha256(evidence_path) != policy["m8_baseline_evidence_sha256"]:
        raise ValueError("M8 evidence changed after the M9 baseline was frozen")
    if m8b.sha256(report_path) != policy["m8_comparison_report_sha256"]:
        raise ValueError("M8 comparison report changed after the M9 baseline was frozen")
    evidence = read_json(evidence_path)
    report = read_json(report_path)
    if evidence["comparison"]["passed_cases"] != policy["required_case_count"]:
        raise ValueError("M8 baseline no longer contains sixteen passing cases")
    if evidence["reference_qualification"]["direct_point_count"] != policy[
            "required_direct_bias_points"]:
        raise ValueError("M8 baseline no longer contains 816 direct points")
    if len(report["cases"]) != policy["required_case_count"]:
        raise ValueError("M8 comparison case count changed")
    if not all(item["status"] == "pass" for item in report["cases"]):
        raise ValueError("M8 baseline contains a non-passing case")
    return {
        "status": "unchanged",
        "evidence_sha256": m8b.sha256(evidence_path),
        "comparison_report_sha256": m8b.sha256(report_path),
        "passing_cases": len(report["cases"]),
        "direct_bias_points": evidence["reference_qualification"][
            "direct_point_count"],
    }


def prepare(contract: dict[str, Any], contract_path: Path,
            output: Path) -> dict[str, Any]:
    tdrs = {key: value.resolve() for key, value in confirmation.DEFAULT_TDRS.items()}
    manifest = m8b.prepare_sentaurus(
        contract, contract_path, tdrs, output)
    manifest["schema"] = "vela.simplemos.sdevice.m9_refdens_sentaurus.v1"
    manifest["scan_values_cm3"] = [
        float(item["refdens_cm3"]) for item in contract["variants"]]
    m8b.write_json(output / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def extract(contract: dict[str, Any], manifest: dict[str, Any], output: Path,
            raw: Path, banner: str) -> dict[str, Any]:
    report = m8b.extract_references(contract, manifest, output, raw, banner)
    report["schema"] = "vela.simplemos.sdevice.m9_refdens_reference.v1"
    m8b.write_json(output / "sentaurus_reference/reference_manifest.json", report)
    return report


def write_summary(path: Path, report: dict[str, Any],
                  variants: dict[str, dict[str, Any]]) -> None:
    rows = []
    for item in report["cases"]:
        variant = variants[item["variant"]]
        rows.append({
            "device": item["device"],
            "drain_voltage_V": item["drain_voltage_V"],
            "variant": item["variant"],
            "refdens_cm3": variant["refdens_cm3"],
            "maximum_sentaurus_response_dex": item[
                "maximum_sentaurus_response_dex"],
            "maximum_absolute_log10_ratio_above_floor": item[
                "maximum_absolute_log10_ratio_above_floor"],
            "maximum_error_change_dex": item["effect_vs_default"][
                "maximum_log10_error_change"],
            "median_error_change_dex": item["effect_vs_default"][
                "median_log10_error_change"],
        })
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def compare(contract: dict[str, Any], output: Path,
            baseline: Path) -> dict[str, Any]:
    report = m8b.compare(contract, output, baseline)
    variants = {item["id"]: item for item in contract["variants"]}
    for item in report["cases"]:
        item["refdens_cm3"] = float(variants[item["variant"]]["refdens_cm3"])
    for item in report["variant_summaries"]:
        item["refdens_cm3"] = float(variants[item["variant"]]["refdens_cm3"])
        item["mean_maximum_error_change_dex"] = statistics.fmean(
            float(value) for value in item["maximum_error_changes_dex"])
    positive = [item for item in report["variant_summaries"]
                if float(item["refdens_cm3"]) > 0.0]
    best = min(positive, key=lambda item: item[
        "mean_maximum_error_change_dex"])
    report.update({
        "schema": "vela.simplemos.sdevice.m9_refdens_comparison.v1",
        "status": "complete",
        "scan_values_cm3": [float(item["refdens_cm3"])
                            for item in contract["variants"]],
        "best_positive_refdens_by_mean_error_change_cm3": float(
            best["refdens_cm3"]),
        "best_positive_mean_error_change_dex": float(
            best["mean_maximum_error_change_dex"]),
        "default_model_changed": False,
    })
    m8b.write_json(output / "comparisons/comparison_report.json", report)
    write_summary(
        output / "comparisons/refdens_scan_summary.csv", report, variants)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--reuse-prepared", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-existing", action="store_true")
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--sentaurus-jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=m8b.m8a.executable("ssh"))
    parser.add_argument("--scp-bin", default=m8b.m8a.executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    args = parser.parse_args()

    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    baseline_guard = validate_baseline_guard(contract)
    tdrs = {key: value.resolve() for key, value in confirmation.DEFAULT_TDRS.items()}
    baseline = args.baseline_dir.resolve()
    m8b.validate_inputs(contract, tdrs, baseline)
    output = args.output_dir.resolve()
    raw = args.raw_dir.resolve()
    manifest = (read_json(output / "sentaurus_matrix_manifest.json")
                if args.reuse_prepared else prepare(contract, contract_path, output))

    reference = None
    if args.live_sentaurus:
        banner = m8b.run_sentaurus(
            manifest, output, args.ssh_target, args.ssh_bin, args.scp_bin,
            args.remote_root, args.sentaurus_jobs, raw)
        reference = extract(contract, manifest, output, raw, banner)
    elif args.extract_existing:
        banner = (output / "sentaurus_banner.txt").read_text(
            encoding="utf-8").strip()
        archive = (output / "sentaurus_raw"
                   / "simplemos_m8b_hfs_controls_results.tgz")
        m8b.extract_archive(archive, raw)
        reference = extract(contract, manifest, output, raw, banner)

    comparison = compare(contract, output, baseline) if args.compare else None
    write_json(output / "m9_execution_summary.json", {
        "schema": "vela.simplemos.sdevice.m9_execution.v1",
        "status": "complete" if comparison else "prepared",
        "baseline_guard": baseline_guard,
        "sentaurus_cases": len(manifest["cases"]),
        "qualified_references": len((reference or {}).get("artifacts", [])),
        "comparison_status": (comparison or {}).get("status"),
        "default_model_changed": False,
    })
    print(json.dumps({
        "sentaurus_cases": len(manifest["cases"]),
        "baseline_guard": baseline_guard["status"],
        "comparison_status": (comparison or {}).get("status"),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
