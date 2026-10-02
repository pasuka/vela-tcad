#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M53 DirectCurrent full-matrix attribution."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m4_controlled_matrix as m4  # noqa: E402
import sentaurus_import  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m53_direct_current_full_matrix_contract_v1.json"
FREEZE = ROOT / "simplemos_m53_direct_current_full_matrix_contract_freeze.json"
M8_CONTRACT = ROOT / "simplemos_m8_original_physics_contract_v1.json"
M46_EVIDENCE = ROOT / "simplemos_m46_full_matrix_requalification_evidence.json"
M46_REPORT = ROOT / "full_matrix_requalification/m46_full_matrix_requalification_report.json"
M46_COMPARISONS = ROOT / "full_matrix_requalification/m46_current_comparisons"
M52_EVIDENCE = ROOT / "simplemos_m52_direct_current_attribution_evidence.json"
M52_REPORT = ROOT / "direct_current_attribution/m52_direct_current_attribution_report.json"
M52_TERMINALS = ROOT / "direct_current_attribution/m52_terminal_component_ledger.csv"
M8_OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m53_direct_current_full_matrix"
PORTABLE = ROOT / "direct_current_full_matrix"
REFERENCES = PORTABLE / "m53_direct_references"
REPORT = PORTABLE / "m53_direct_current_full_matrix_report.json"
CASES = PORTABLE / "m53_case_summary.csv"
POINTS = PORTABLE / "m53_pointwise_comparison.csv"
NWELL = PORTABLE / "m53_nwell_pair_summary.csv"
TERMINALS = PORTABLE / "m53_terminal_replay_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m53_direct_current_full_matrix_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m53/artifact.json"
EVIDENCE = ROOT / "simplemos_m53_direct_current_full_matrix_evidence.json"
TEST = REPO / "tests/regression/test_simplemos_m53_direct_current_full_matrix.py"
ARCHIVE = "simplemos_m53_direct_current_full_matrix_results.tgz"
REMOTE_ROOT = "/tmp/vela_simplemos_m53_direct_current_full_matrix"
CONTACTS = ("source", "drain", "gate", "substrate")
COMPONENTS = ("electron", "hole", "total")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(resolved)


def voltage_tag(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p")


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def close_enough(left: float, right: float, abs_tol: float,
                 rel_tol: float) -> bool:
    difference = abs(left - right)
    scale = max(abs(left), abs(right), 1e-300)
    return difference <= abs_tol or difference / scale <= rel_tol


def trend(values: list[float]) -> str:
    diffs = [right - left for left, right in zip(values, values[1:])]
    tolerance = max(max((abs(value) for value in values), default=0.0) * 1e-12,
                    1e-300)
    if all(value >= -tolerance for value in diffs):
        return "nondecreasing"
    if all(value <= tolerance for value in diffs):
        return "nonincreasing"
    return "mixed"


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m53_direct_current_full_matrix_contract.v1"):
        raise ValueError("unexpected M53 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M53 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M53 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M53 upstream evidence changed: {relative}")
    m46_evidence = read_json(M46_EVIDENCE)
    m52_evidence = read_json(M52_EVIDENCE)
    if m46_evidence.get("status") != "frozen" or m52_evidence.get("status") != "frozen":
        raise ValueError("M53 requires frozen M46 and M52 evidence")
    m46 = read_json(M46_REPORT)
    m52 = read_json(M52_REPORT)
    upstream = contract["upstream"]
    checks = (
        (m46["findings"]["curve_count"], upstream["required_m46_curve_count"]),
        (m46["findings"]["direct_bias_point_count"],
         upstream["required_m46_point_count"]),
    )
    if any(int(observed) != int(expected) for observed, expected in checks):
        raise ValueError("M46 matrix anchor changed")
    if not math.isclose(
            float(m46["findings"]["maximum_absolute_log10_ratio_dex"]),
            float(upstream["required_m46_maximum_absolute_log10_ratio_dex"]),
            rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("M46 log-error anchor changed")
    if m52.get("classification") != upstream["required_m52_classification"]:
        raise ValueError("M52 classification anchor changed")
    devices = tuple(contract["matrix"]["devices"])
    drains = tuple(map(float, contract["matrix"]["drain_voltages_V"]))
    if devices != tuple(f"n{index}" for index in range(17, 25)):
        raise ValueError("M53 device matrix changed")
    if drains != (0.05, 1.0):
        raise ValueError("M53 drain matrix changed")
    m8 = read_json(M8_CONTRACT)
    if tuple(item["id"] for item in m8["devices"]) != devices:
        raise ValueError("M8 device order changed")
    if len(m8["bias_matrix"]["gate_lattice"]["values_V"]) != 51:
        raise ValueError("M8 gate lattice changed")
    return contract


def prepare(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force:
        for path in (bundle, OUTPUT / "sentaurus_raw"):
            if path.exists():
                shutil.rmtree(path)
    bundle.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    marker = "Math { Extrapolate Iterations=20 ExitOnFailure }"
    for device in contract["matrix"]["devices"]:
        source_root = M8_OUTPUT / "sentaurus_bundle" / device
        source_tdr = source_root / "input_fps.tdr"
        if not source_tdr.is_file():
            raise FileNotFoundError(source_tdr)
        target_root = bundle / device
        target_root.mkdir(parents=True, exist_ok=True)
        target_tdr = target_root / "input_fps.tdr"
        shutil.copy2(source_tdr, target_tdr)
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            case = f"{device}_vd_{voltage_tag(drain)}"
            source_deck = source_root / f"{case}_des.cmd"
            baseline = source_deck.read_text(encoding="utf-8")
            if baseline.count(marker) != 1:
                raise RuntimeError(f"M8 Math insertion point changed: {source_deck}")
            if "DirectCurrent" in baseline or "CurrentWeighting" in baseline:
                raise RuntimeError(f"M8 deck contains a current algorithm flag: {source_deck}")
            direct = baseline.replace(
                marker, "Math { Extrapolate Iterations=20 ExitOnFailure DirectCurrent }")
            if direct.replace(" DirectCurrent", "") != baseline:
                raise RuntimeError(f"M53 deck differs by more than DirectCurrent: {case}")
            target_deck = target_root / f"{case}_des.cmd"
            target_deck.write_text(direct, encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "device": device,
                "drain_voltage_V": drain,
                "deck": portable(target_deck),
                "deck_sha256": sha256(target_deck),
                "baseline_deck": portable(source_deck),
                "baseline_deck_sha256": sha256(source_deck),
                "input_tdr": portable(target_tdr),
                "input_tdr_sha256": sha256(target_tdr),
                "m8_input_tdr_sha256": sha256(source_tdr),
                "expected_plot": f"IdVg_{case}_des.plt",
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m53_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "intervention": "Math DirectCurrent flag only",
        "cases": cases,
    }
    write_json(OUTPUT / "sentaurus_manifest.json", manifest)
    return manifest


def run_sentaurus(manifest: dict[str, Any], ssh_target: str, ssh_bin: str,
                   scp_bin: str, remote_root: str, jobs: int) -> str:
    banner = run([ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
                 capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target,
         f"set -eu; test ! -e {remote_root}; mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(OUTPUT / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])
    by_device: dict[str, list[dict[str, Any]]] = {}
    for item in manifest["cases"]:
        by_device.setdefault(item["device"], []).append(item)

    def execute_device(item: tuple[str, list[dict[str, Any]]]) -> None:
        device, cases = item
        root = f"{remote_root}/sentaurus_bundle/{device}"
        commands = [f"cd {root}"]
        for case in cases:
            name = case["case"]
            commands.append(f"sdevice {name}_des.cmd > {name}.console.log 2>&1")
        run([ssh_bin, ssh_target, "set -eu; " + "; ".join(commands)])

    items = list(by_device.items())
    with ThreadPoolExecutor(max_workers=min(max(jobs, 1), len(items))) as pool:
        list(pool.map(execute_device, items))
    run([ssh_bin, ssh_target,
         f"cd {remote_root}; tar -czf {ARCHIVE} sentaurus_bundle"])
    raw = OUTPUT / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / ARCHIVE
    run([scp_bin, f"{ssh_target}:{remote_root}/{ARCHIVE}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (OUTPUT / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def exact_terminal_curve(path: Path, gates: list[float], tolerance: float
                         ) -> list[dict[str, float]]:
    text = path.read_text(errors="ignore")
    datasets = sentaurus_import.parse_quoted_list(text, "datasets")
    required = ["gate OuterVoltage", "drain OuterVoltage"]
    for contact in CONTACTS:
        required.extend((f"{contact} eCurrent", f"{contact} hCurrent",
                         f"{contact} TotalCurrent"))
    missing = set(required) - set(datasets)
    if missing:
        raise RuntimeError(f"missing M53 datasets in {path}: {sorted(missing)}")
    source = sentaurus_import.parse_values_block(text, len(datasets))
    rows: list[dict[str, float]] = []
    for gate in gates:
        matches = [values for values in source
                   if abs(float(values[datasets.index("gate OuterVoltage")]) - gate)
                   <= tolerance]
        if len(matches) != 1:
            raise RuntimeError(
                f"{path}: expected one direct row at Vg={gate:g}; got {len(matches)}")
        raw = dict(zip(datasets, matches[0], strict=True))
        row = {
            "gate_voltage_V": gate,
            "drain_voltage_V": float(raw["drain OuterVoltage"]),
        }
        for contact in CONTACTS:
            row[f"{contact}_electron_A_per_um"] = float(raw[f"{contact} eCurrent"])
            row[f"{contact}_hole_A_per_um"] = float(raw[f"{contact} hCurrent"])
            row[f"{contact}_total_A_per_um"] = float(raw[f"{contact} TotalCurrent"])
        if not all(math.isfinite(value) for value in row.values()):
            raise RuntimeError(f"non-finite M53 terminal value in {path}")
        rows.append(row)
    return rows


def log_error(candidate: float, reference: float, floor: float
              ) -> tuple[float, float, bool]:
    candidate_magnitude = abs(candidate)
    reference_magnitude = abs(reference)
    above_floor = reference_magnitude >= floor and candidate_magnitude > 0.0
    if not above_floor:
        return math.nan, math.nan, False
    return (abs(math.log10(candidate_magnitude / reference_magnitude)),
            abs(candidate_magnitude - reference_magnitude) / reference_magnitude,
            True)


def analyze(contract: dict[str, Any], manifest: dict[str, Any], banner: str
            ) -> dict[str, Any]:
    gates = [0.05 * index for index in range(51)]
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    floor = float(contract["analysis"]["current_floor_A_per_um"])
    raw_root = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    parsed: dict[str, list[dict[str, float]]] = {}
    REFERENCES.mkdir(parents=True, exist_ok=True)
    for item in manifest["cases"]:
        plot = raw_root / item["device"] / item["expected_plot"]
        if not plot.is_file():
            raise FileNotFoundError(plot)
        rows = exact_terminal_curve(plot, gates, tolerance)
        parsed[item["case"]] = rows
        write_csv(REFERENCES / f"{item['case']}_reference.csv", [
            {"gate_voltage_V": row["gate_voltage_V"],
             "drain_total_current_A_per_um": row["drain_total_A_per_um"]}
            for row in rows])

    point_rows: list[dict[str, Any]] = []
    case_rows: list[dict[str, Any]] = []
    worst: dict[str, Any] | None = None
    parity = contract["acceptance"]
    for item in manifest["cases"]:
        case = item["case"]
        baseline_rows = read_csv(M46_COMPARISONS / f"{case}_comparison.csv")
        if len(baseline_rows) != 51:
            raise RuntimeError(f"M46 baseline is not 51 points: {case}")
        direct_values: list[float] = []
        vela_values: list[float] = []
        case_points: list[dict[str, Any]] = []
        for direct, old in zip(parsed[case], baseline_rows, strict=True):
            gate = float(direct["gate_voltage_V"])
            baseline_gate = float(old["gate_voltage_V"])
            if not math.isclose(gate, baseline_gate, rel_tol=0.0,
                                abs_tol=tolerance):
                raise RuntimeError(
                    f"M46 exact gate mismatch for {case}: "
                    f"M53={gate}, M46={baseline_gate}")
            default_current = float(old["sentaurus_current_A_per_um"])
            vela_current = float(old["vela_current_A_per_um"])
            direct_current = float(direct["drain_total_A_per_um"])
            direct_log, direct_relative, direct_above = log_error(
                vela_current, direct_current, floor)
            default_log, default_relative, default_above = log_error(
                vela_current, default_current, floor)
            default_direct_log = (
                abs(math.log10(abs(direct_current) / abs(default_current)))
                if direct_current != 0.0 and default_current != 0.0 else math.nan)
            row = {
                "case": case,
                "device": item["device"],
                "drain_voltage_V": item["drain_voltage_V"],
                "gate_voltage_V": gate,
                "default_sentaurus_current_A_per_um": default_current,
                "direct_sentaurus_current_A_per_um": direct_current,
                "vela_current_A_per_um": vela_current,
                "direct_minus_default_A_per_um": direct_current - default_current,
                "absolute_default_to_direct_log10_ratio_dex": default_direct_log,
                "default_vela_above_floor": int(default_above),
                "default_vela_absolute_log10_ratio_dex": default_log,
                "default_vela_relative_error": default_relative,
                "direct_vela_above_floor": int(direct_above),
                "direct_vela_absolute_log10_ratio_dex": direct_log,
                "direct_vela_relative_error": direct_relative,
            }
            point_rows.append(row)
            case_points.append(row)
            direct_values.append(abs(direct_current))
            vela_values.append(abs(vela_current))
            if direct_above and (worst is None or direct_log >
                                 float(worst["direct_vela_absolute_log10_ratio_dex"])):
                worst = row
        direct_logs = [float(row["direct_vela_absolute_log10_ratio_dex"])
                       for row in case_points if row["direct_vela_above_floor"]]
        direct_relatives = [float(row["direct_vela_relative_error"])
                            for row in case_points if row["direct_vela_above_floor"]]
        default_logs = [float(row["default_vela_absolute_log10_ratio_dex"])
                        for row in case_points if row["default_vela_above_floor"]]
        direct_trend = trend(direct_values)
        vela_trend = trend(vela_values)
        max_direct = max(direct_logs)
        max_relative = max(direct_relatives)
        passes = (
            max_direct <= float(parity["original_parity_maximum_absolute_log10_ratio_dex"])
            and max_relative <= float(parity["original_parity_maximum_relative_error"])
            and (not parity["trend_match_required_for_original_parity"]
                 or direct_trend == vela_trend))
        case_rows.append({
            "case": case,
            "device": item["device"],
            "drain_voltage_V": item["drain_voltage_V"],
            "point_count": len(case_points),
            "default_maximum_absolute_log10_ratio_dex": max(default_logs),
            "direct_maximum_absolute_log10_ratio_dex": max_direct,
            "direct_maximum_relative_error": max_relative,
            "maximum_error_change_direct_minus_default_dex":
                max_direct - max(default_logs),
            "direct_reference_trend": direct_trend,
            "vela_trend": vela_trend,
            "trend_match": direct_trend == vela_trend,
            "original_parity_status": "pass" if passes else "fail",
        })

    case_index = {row["case"]: row for row in case_rows}
    nwell_rows: list[dict[str, Any]] = []
    pair_tolerance = float(
        contract["acceptance"]["global_error_comparison_tolerance_dex"])
    for pair in contract["matrix"]["matched_nwell_pairs"]:
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            tag = voltage_tag(drain)
            low = case_index[f"{pair['low']}_vd_{tag}"]
            high = case_index[f"{pair['high']}_vd_{tag}"]
            default_amplification = (
                float(high["default_maximum_absolute_log10_ratio_dex"])
                - float(low["default_maximum_absolute_log10_ratio_dex"]))
            direct_amplification = (
                float(high["direct_maximum_absolute_log10_ratio_dex"])
                - float(low["direct_maximum_absolute_log10_ratio_dex"]))
            nwell_rows.append({
                "low_nwell_device": pair["low"],
                "high_nwell_device": pair["high"],
                "drain_voltage_V": drain,
                "low_default_max_error_dex":
                    low["default_maximum_absolute_log10_ratio_dex"],
                "high_default_max_error_dex":
                    high["default_maximum_absolute_log10_ratio_dex"],
                "default_high_minus_low_amplification_dex": default_amplification,
                "low_direct_max_error_dex":
                    low["direct_maximum_absolute_log10_ratio_dex"],
                "high_direct_max_error_dex":
                    high["direct_maximum_absolute_log10_ratio_dex"],
                "direct_high_minus_low_amplification_dex": direct_amplification,
                "direct_minus_default_amplification_change_dex":
                    direct_amplification - default_amplification,
                "direct_did_not_increase_nwell_amplification":
                    direct_amplification <= default_amplification + pair_tolerance,
            })

    m52_rows = read_csv(M52_TERMINALS)
    terminal_rows: list[dict[str, Any]] = []
    terminal_abs = float(
        contract["acceptance"]["terminal_replay_absolute_tolerance_A_per_um"])
    terminal_rel = float(
        contract["acceptance"]["terminal_replay_relative_tolerance"])
    for anchor in m52_rows:
        device = anchor["device"]
        gate = float(anchor["gate_voltage_V"])
        case = f"{device}_vd_0p05"
        direct = next(row for row in parsed[case]
                      if float(row["gate_voltage_V"]) == gate)
        key = f"{anchor['contact']}_{anchor['component']}_A_per_um"
        observed = float(direct[key])
        expected = float(anchor["m52_direct_A_per_um"])
        difference = observed - expected
        relative = abs(difference) / max(abs(observed), abs(expected), 1e-300)
        terminal_rows.append({
            "state": anchor["state"],
            "device": device,
            "gate_voltage_V": gate,
            "contact": anchor["contact"],
            "component": anchor["component"],
            "m52_direct_A_per_um": expected,
            "m53_direct_A_per_um": observed,
            "m53_minus_m52_A_per_um": difference,
            "symmetric_relative_difference": relative,
            "within_tolerance": close_enough(
                observed, expected, terminal_abs, terminal_rel),
        })

    required_curves = int(contract["acceptance"]["required_curve_count"])
    required_points = int(contract["acceptance"]["required_point_count"])
    required_terminals = int(
        contract["acceptance"]["required_terminal_replay_rows"])
    complete = (len(case_rows) == required_curves
                and len(point_rows) == required_points
                and all(int(row["point_count"]) == 51 for row in case_rows)
                and len(terminal_rows) == required_terminals)
    terminal_replay = all(bool(row["within_tolerance"]) for row in terminal_rows)
    default_global = max(float(row["default_maximum_absolute_log10_ratio_dex"])
                         for row in case_rows)
    direct_global = max(float(row["direct_maximum_absolute_log10_ratio_dex"])
                        for row in case_rows)
    if not complete:
        classification = "replay_incomplete"
    elif not terminal_replay:
        classification = "terminal_replay_mismatch"
    elif direct_global < default_global - pair_tolerance:
        classification = (
            "direct_reduces_global_and_nwell_amplification"
            if all(bool(row["direct_did_not_increase_nwell_amplification"])
                   for row in nwell_rows)
            else "direct_reduces_global_partial_nwell")
    elif abs(direct_global - default_global) <= pair_tolerance:
        classification = "direct_global_error_unchanged"
    else:
        classification = "direct_increases_global_error"

    all_finite = all(math.isfinite(float(row[key])) for row in point_rows
                     for key in ("default_sentaurus_current_A_per_um",
                                 "direct_sentaurus_current_A_per_um",
                                 "vela_current_A_per_um"))
    only_direct_intervention = True
    for item in manifest["cases"]:
        direct_text = (REPO / item["deck"]).read_text(encoding="utf-8")
        baseline_text = (REPO / item["baseline_deck"]).read_text(encoding="utf-8")
        only_direct_intervention = only_direct_intervention and (
            direct_text.count(" DirectCurrent") == 1
            and "CurrentWeighting" not in direct_text
            and direct_text.replace(" DirectCurrent", "") == baseline_text
            and item["input_tdr_sha256"] == item["m8_input_tdr_sha256"])
    acceptance = {
        "contract_hash_unchanged": sha256(CONTRACT) ==
            read_json(FREEZE)["contract_sha256"],
        "curve_count": len(case_rows) == required_curves,
        "point_count": len(point_rows) == required_points,
        "points_per_curve": all(int(row["point_count"]) == 51
                                for row in case_rows),
        "terminal_replay_row_count": len(terminal_rows) == required_terminals,
        "terminal_replay_within_tolerance": terminal_replay,
        "all_values_finite": all_finite,
        "classification_declared": classification in
            contract["analysis"]["classifications"],
        "only_direct_current_intervention": only_direct_intervention,
        "closed_topics_not_reopened": True,
        "hypothesis_rejection_is_valid_completion": True,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    target = next(row for row in point_rows
                  if row["case"] == "n23_vd_0p05"
                  and float(row["gate_voltage_V"]) == 0.05)
    max_terminal_abs = max(abs(float(row["m53_minus_m52_A_per_um"]))
                           for row in terminal_rows)
    max_terminal_rel = max(float(row["symmetric_relative_difference"])
                           for row in terminal_rows)
    report = {
        "schema": "vela.simplemos.sdevice.m53_direct_current_full_matrix_report.v1",
        "status": "complete" if acceptance["all_checks_pass"] else "failed",
        "classification": classification,
        "execution": {
            "sentaurus_release": "T-2022.03-SP2",
            "sentaurus_banner": banner,
            "new_sentaurus_execution": True,
            "new_vela_execution": False,
            "intervention": "Math DirectCurrent flag only",
            "production_defaults_changed": False,
            "historical_artifacts_rewritten": False,
            "contract_sha256_before_and_after": sha256(CONTRACT),
        },
        "matrix": {
            "device_count": 8,
            "curve_count": len(case_rows),
            "direct_bias_point_count": len(point_rows),
            "points_per_curve": 51,
        },
        "global_findings": {
            "default_maximum_absolute_log10_ratio_dex": default_global,
            "direct_maximum_absolute_log10_ratio_dex": direct_global,
            "direct_minus_default_global_maximum_error_dex":
                direct_global - default_global,
            "direct_original_parity_passing_curve_count": sum(
                row["original_parity_status"] == "pass" for row in case_rows),
            "direct_original_parity_failing_curve_count": sum(
                row["original_parity_status"] == "fail" for row in case_rows),
            "worst_direct_point": worst,
        },
        "target": target,
        "nwell_findings": {
            "pair_drain_row_count": len(nwell_rows),
            "rows_without_increased_amplification": sum(
                bool(row["direct_did_not_increase_nwell_amplification"])
                for row in nwell_rows),
            "all_rows_without_increased_amplification": all(
                bool(row["direct_did_not_increase_nwell_amplification"])
                for row in nwell_rows),
        },
        "m52_terminal_replay": {
            "row_count": len(terminal_rows),
            "all_within_tolerance": terminal_replay,
            "maximum_absolute_difference_A_per_um": max_terminal_abs,
            "maximum_symmetric_relative_difference": max_terminal_rel,
        },
        "case_results": case_rows,
        "nwell_pairs": nwell_rows,
        "acceptance": acceptance,
        "claim_guard": contract["analysis"]["claim_guard"],
    }
    write_csv(CASES, case_rows)
    write_csv(POINTS, point_rows)
    write_csv(NWELL, nwell_rows)
    write_csv(TERMINALS, terminal_rows)
    write_json(REPORT, report)
    return report


def freeze_artifacts(report: dict[str, Any], manifest: dict[str, Any]) -> None:
    findings = report["global_findings"]
    target = report["target"]
    worst = findings["worst_direct_point"]
    nwell_table = "\n".join(
        f"| {row['low_nwell_device']}/{row['high_nwell_device']} | "
        f"{float(row['drain_voltage_V']):.2f} | "
        f"{float(row['default_high_minus_low_amplification_dex']):.6f} | "
        f"{float(row['direct_high_minus_low_amplification_dex']):.6f} | "
        f"{row['direct_did_not_increase_nwell_amplification']} |"
        for row in report["nwell_pairs"])
    doc = f'''# SimpleMOS M53 DirectCurrent 完整矩阵重资格

## 结论

M53 分类为 `{report['classification']}`。在冻结 M8/M46 TDR、物理、网格、偏压路径和 Vela 曲线后，唯一干预是在 16 个 Sentaurus 原始 deck 的全局 `Math` 中加入 `DirectCurrent`。

完整执行 8 个器件、2 个漏压、16 条 Id-Vg 曲线和 816 个精确栅压点。默认端口观测下的全局最大 Vela-Sentaurus 误差为 `{float(findings['default_maximum_absolute_log10_ratio_dex']):.12g}` dex；`DirectCurrent` 下为 `{float(findings['direct_maximum_absolute_log10_ratio_dex']):.12g}` dex，变化 `{float(findings['direct_minus_default_global_maximum_error_dex']):+.12g}` dex。按原 M8 数值/趋势门槛，`{findings['direct_original_parity_passing_curve_count']}` 条曲线通过、`{findings['direct_original_parity_failing_curve_count']}` 条失败。

最差 DirectCurrent 点为 `{worst['case']}`、Vg=`{float(worst['gate_voltage_V']):.2f}` V，误差 `{float(worst['direct_vela_absolute_log10_ratio_dex']):.12g}` dex。原 M46 目标 n23、Vd=0.05 V、Vg=0.05 V 的默认、DirectCurrent 和 Vela 漏电流分别为 `{float(target['default_sentaurus_current_A_per_um']):.12e}`、`{float(target['direct_sentaurus_current_A_per_um']):.12e}`、`{float(target['vela_current_A_per_um']):.12e}` A/um。

## M52 锚点回放

六个共享状态的四端电子、空穴和总电流共 `{report['m52_terminal_replay']['row_count']}` 行全部回放：`{report['m52_terminal_replay']['all_within_tolerance']}`。最大绝对差为 `{float(report['m52_terminal_replay']['maximum_absolute_difference_A_per_um']):.12e}` A/um，最大对称相对差为 `{float(report['m52_terminal_replay']['maximum_symmetric_relative_difference']):.12e}`。

## NWell 配对

| 低/高 NWell 器件 | Vd (V) | 默认高减低最大误差 (dex) | Direct 高减低最大误差 (dex) | 未增加放大 |
|---|---:|---:|---:|---|
{nwell_table}

`DirectCurrent` 在 `{report['nwell_findings']['rows_without_increased_amplification']}/{report['nwell_findings']['pair_drain_row_count']}` 个 NWell-漏压配对中没有增加高 NWell 相对低 NWell 的最大误差放大。

## 边界

- M53 检验的是 Sentaurus 端口观测算法，不修改任何求解状态物理或生产默认值。
- M8/M46 官方原始 deck 默认结果仍是主要兼容性基线；`DirectCurrent` 结果不会静默替代它。
- 假设被否证也是合同允许的完成结果；未针对运行结果修改阈值。
- 未重新排查 HFS、SG、准费米、BGN、SRH、迁移率、接触模型或网格。

机器报告：`{portable(REPORT)}`。
'''
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc, encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.simplemos.sdevice.m53_artifact.v1",
        "title": "SimpleMOS M53 DirectCurrent full-matrix attribution",
        "status": report["status"],
        "classification": report["classification"],
        "report": portable(REPORT),
        "ledgers": [portable(CASES), portable(POINTS), portable(NWELL),
                    portable(TERMINALS)],
        "direct_references": portable(REFERENCES),
        "document": portable(DOC),
        "global_findings": report["global_findings"],
        "acceptance": report["acceptance"],
    })
    artifacts = [REPORT, CASES, POINTS, NWELL, TERMINALS, DOC, ARTIFACT,
                 *sorted(REFERENCES.glob("*_reference.csv"))]
    sources = [CONTRACT, FREEZE, Path(__file__).resolve(), TEST,
               M46_EVIDENCE, M52_EVIDENCE, M46_REPORT, M52_REPORT,
               M52_TERMINALS, OUTPUT / "sentaurus_manifest.json",
               OUTPUT / "sentaurus_banner.txt"]
    for item in manifest["cases"]:
        sources.append(REPO / item["deck"])
        sources.append(OUTPUT / "sentaurus_raw/sentaurus_bundle" /
                       item["device"] / item["expected_plot"])
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m53_direct_current_full_matrix_evidence.v1",
        "status": "frozen",
        "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "new_sentaurus_execution": True,
        "new_vela_execution": False,
        "only_intervention": "Math DirectCurrent flag",
        "production_defaults_changed": False,
        "historical_artifacts_rewritten": False,
        "closed_topics_reinvestigated": False,
        "classification": report["classification"],
        "acceptance": report["acceptance"],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-sentaurus", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=r"C:\Windows\System32\OpenSSH\ssh.exe")
    parser.add_argument("--scp-bin", default=r"C:\Windows\System32\OpenSSH\scp.exe")
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    contract = validate_contract()
    do_prepare = args.prepare or args.all
    do_run = args.run_sentaurus or args.all
    do_analyze = args.analyze or args.all
    manifest_path = OUTPUT / "sentaurus_manifest.json"
    manifest = prepare(contract, args.force) if do_prepare else read_json(manifest_path)
    if not do_run and not do_analyze:
        print(json.dumps({"status": "prepared", "cases": len(manifest["cases"]),
                          "manifest": portable(manifest_path)}, indent=2))
        return
    banner_path = OUTPUT / "sentaurus_banner.txt"
    banner = (run_sentaurus(manifest, args.ssh_target, args.ssh_bin,
                            args.scp_bin, args.remote_root, args.jobs)
              if do_run else banner_path.read_text(encoding="utf-8").strip())
    if do_analyze:
        report = analyze(contract, manifest, banner)
        freeze_artifacts(report, manifest)
        print(json.dumps({
            "status": report["status"],
            "classification": report["classification"],
            "all_checks_pass": report["acceptance"]["all_checks_pass"],
            "report": portable(REPORT),
        }, indent=2))


if __name__ == "__main__":
    main()
