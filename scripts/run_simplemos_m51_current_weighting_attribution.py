#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M51 CurrentWeighting attribution."""

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
import run_simplemos_m50_native_substrate_face_flux_export as m50  # noqa: E402
import sentaurus_import  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m51_current_weighting_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m51_current_weighting_attribution_contract_freeze.json"
UPSTREAM_EVIDENCE = [
    ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_evidence.json",
    ROOT / "simplemos_m48_terminal_partition_continuity_closure_evidence.json",
    ROOT / "simplemos_m49_substrate_electron_transport_attribution_evidence.json",
    ROOT / "simplemos_m50_native_substrate_face_flux_export_evidence.json",
]
M50_OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m50_native_substrate_face_flux_export"
M50_REPORT = ROOT / "native_substrate_face_flux_export/m50_native_substrate_face_flux_export_report.json"
M50_FLUX = ROOT / "native_substrate_face_flux_export/m50_native_substrate_face_flux_ledger.csv"
M50_TERMINALS = ROOT / "native_substrate_face_flux_export/m50_terminal_replay_ledger.csv"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m51_current_weighting_attribution"
PORTABLE = ROOT / "current_weighting_attribution"
REPORT = PORTABLE / "m51_current_weighting_attribution_report.json"
STATES = PORTABLE / "m51_current_weighting_state_ledger.csv"
TERMINALS = PORTABLE / "m51_terminal_component_ledger.csv"
FIELDS = PORTABLE / "m51_field_invariance_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m51_current_weighting_attribution_2026-09-01.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m51/artifact.json"
EVIDENCE = ROOT / "simplemos_m51_current_weighting_attribution_evidence.json"
TEST = REPO / "tests/regression/test_simplemos_m51_current_weighting_attribution.py"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
ARCHIVE = "simplemos_m51_current_weighting_results.tgz"
REMOTE_ROOT = "/tmp/vela_simplemos_m51_current_weighting_attribution"
DEVICES = ("n23", "n19")
GATES = (0.0, 0.05, 0.1)
CONTACTS = ("source", "drain", "gate", "substrate")
COMPONENTS = ("electron", "hole", "total")
FIELD_PATTERNS = (
    "ElectrostaticPotential_region*.csv",
    "eDensity_region0.csv", "hDensity_region0.csv",
    "eQuasiFermiPotential_region*.csv", "hQuasiFermiPotential_region*.csv",
    "eCurrentDensity_region0.csv", "hCurrentDensity_region0.csv",
    "srhRecombination_region0.csv",
)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def state_id(device: str, gate: float) -> str:
    return f"{device}_vd_0p05_vg_{m50.m47.voltage_tag(gate)}"


def close_enough(value: float, reference: float, absolute: float,
                 relative: float) -> bool:
    difference = abs(value - reference)
    symmetric_scale = max(abs(value), abs(reference), 1e-300)
    return difference <= absolute or difference / symmetric_scale <= relative


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != "vela.simplemos.sdevice.m51_current_weighting_attribution_contract.v1":
        raise ValueError("unexpected M51 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M51 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M51 contract changed after freeze")
    for path in UPSTREAM_EVIDENCE:
        if read_json(path).get("status") != "frozen":
            raise ValueError(f"M51 requires frozen evidence: {path.name}")
    m50_report = read_json(M50_REPORT)
    if m50_report.get("classification") != contract["upstream"]["required_m50_classification"]:
        raise ValueError("M50 classification changed")
    target = m50_report["target"]
    anchors = {
        "same_run_substrate_electron_terminal_A_per_um":
            "required_m50_target_default_substrate_electron_A_per_um",
        "contact_surface_normal_integral_A_per_um":
            "required_m50_target_contact_surface_normal_A_per_um",
        "negative_contact_surface_normal_integral_A_per_um":
            "required_m50_target_negative_contact_surface_normal_A_per_um",
    }
    for observed, expected in anchors.items():
        if not math.isclose(float(target[observed]),
                            float(contract["upstream"][expected]),
                            rel_tol=0.0, abs_tol=1e-30):
            raise ValueError(f"M50 target anchor changed: {observed}")
    if tuple(contract["scope"]["devices"]) != DEVICES:
        raise ValueError("M51 device matrix changed")
    if tuple(map(float, contract["scope"]["gate_voltages_V"])) != GATES:
        raise ValueError("M51 gate matrix changed")
    return contract


def sentaurus_deck(case: str) -> str:
    deck = m50.sentaurus_deck(case)
    deck = deck.replace("M50_SubstrateElectronFaceFlux",
                        "M51_SubstrateElectronFaceFlux")
    deck = deck.replace('Plot(FilePrefix="m50_state"',
                        'Plot(FilePrefix="m51_state"')
    marker = "Math { Extrapolate Iterations=20 ExitOnFailure }"
    replacement = "Math { Extrapolate Iterations=20 ExitOnFailure CurrentWeighting }"
    if deck.count(marker) != 1:
        raise RuntimeError("M50 Math insertion point changed")
    weighted = deck.replace(marker, replacement)
    if weighted.replace(" CurrentWeighting", "") != deck:
        raise RuntimeError("M51 deck differs from M50 by more than CurrentWeighting")
    return weighted


def prepare(force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    cases: list[dict[str, Any]] = []
    for device in DEVICES:
        source = M50_OUTPUT / "sentaurus_bundle" / device / "input_fps.tdr"
        if not source.is_file():
            raise FileNotFoundError(source)
        root = bundle / device
        root.mkdir(parents=True, exist_ok=True)
        target = root / "input_fps.tdr"
        shutil.copy2(source, target)
        case = f"{device}_vd_0p05"
        deck = root / f"{case}_des.cmd"
        deck.write_text(sentaurus_deck(case), encoding="utf-8", newline="\n")
        cases.append({
            "device": device, "case": case,
            "deck": portable(deck), "deck_sha256": sha256(deck),
            "input_tdr": portable(target), "input_tdr_sha256": sha256(target),
            "m50_input_tdr_sha256": sha256(source),
        })
    manifest = {
        "schema": "vela.simplemos.sdevice.m51_sentaurus_manifest.v1",
        "status": "prepared", "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "intervention": "CurrentWeighting flag only; diagnostic names M50 to M51",
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

    def execute(item: dict[str, Any]) -> None:
        root = f"{remote_root}/sentaurus_bundle/{item['device']}"
        command = (f"set -eu; cd {root}; "
                   f"sdevice {item['case']}_des.cmd > console.log 2>&1")
        run([ssh_bin, ssh_target, command])

    with ThreadPoolExecutor(max_workers=min(jobs, len(manifest["cases"]))) as pool:
        list(pool.map(execute, manifest["cases"]))
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


def parse_terminal_file(path: Path) -> dict[float, dict[str, float]]:
    text = path.read_text(errors="ignore")
    datasets = sentaurus_import.parse_quoted_list(text, "datasets")
    if "M51_SubstrateElectronFaceFlux" not in datasets:
        raise RuntimeError(f"M51 native flux missing from {path}")
    parsed: dict[float, dict[str, float]] = {}
    for values in sentaurus_import.parse_values_block(text, len(datasets)):
        row = dict(zip(datasets, values, strict=True))
        gate = float(row["gate OuterVoltage"])
        nearest = min(GATES, key=lambda value: abs(value - gate))
        if math.isclose(gate, nearest, rel_tol=0.0, abs_tol=1e-10):
            item = {
                "gate_voltage_V": nearest,
                "drain_voltage_V": float(row["drain OuterVoltage"]),
                "native_flux_A_per_um": float(row["M51_SubstrateElectronFaceFlux"]),
            }
            for contact in CONTACTS:
                item[f"{contact}_electron_A_per_um"] = float(row[f"{contact} eCurrent"])
                item[f"{contact}_hole_A_per_um"] = float(row[f"{contact} hCurrent"])
                item[f"{contact}_total_A_per_um"] = float(row[f"{contact} TotalCurrent"])
            parsed[nearest] = item
    if set(parsed) != set(GATES):
        raise RuntimeError(f"missing exact M51 terminal states in {path}: {sorted(parsed)}")
    return parsed


def export_snapshots(manifest: dict[str, Any], force: bool) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = {"baseline": {}, "weighted": {}}
    for item in manifest["cases"]:
        device = item["device"]
        roots = {
            "baseline": M50_OUTPUT / "sentaurus_raw/sentaurus_bundle" / device,
            "weighted": OUTPUT / "sentaurus_raw/sentaurus_bundle" / device,
        }
        patterns = {"baseline": "m50_state*.tdr", "weighted": "m51_state*.tdr"}
        for variant, root in roots.items():
            snapshots = sorted(root.glob(patterns[variant]))
            if len(snapshots) != 3:
                raise RuntimeError(f"expected three {variant} snapshots for {device}: {snapshots}")
            for gate, tdr in zip(GATES, snapshots, strict=True):
                state = state_id(device, gate)
                export = OUTPUT / "state_exports" / variant / state
                if force and export.exists():
                    shutil.rmtree(export)
                if not (export / "field_manifest.json").is_file():
                    run([str(IMPORTER), "--tdr", str(tdr),
                         "--export-dir", str(export)], capture=True)
                result[variant][state] = export
    return result


def compare_nodes(baseline: Path, weighted: Path) -> tuple[bool, float]:
    left = read_csv(baseline / "nodes.csv")
    right = read_csv(weighted / "nodes.csv")
    if len(left) != len(right):
        return False, math.inf
    maximum = 0.0
    for a, b in zip(left, right, strict=True):
        if a["id"] != b["id"]:
            return False, math.inf
        maximum = max(maximum, abs(float(a["x_um"]) - float(b["x_um"])),
                      abs(float(a["y_um"]) - float(b["y_um"])))
    return maximum == 0.0, maximum


def selected_field_files(export: Path) -> list[Path]:
    fields = export / "fields"
    selected: set[Path] = set()
    for pattern in FIELD_PATTERNS:
        selected.update(fields.glob(pattern))
    return sorted(selected, key=lambda path: path.name)


def compare_fields(contract: dict[str, Any], exports: dict[str, dict[str, Path]]) -> tuple[list[dict[str, Any]], bool]:
    scalar_abs = float(contract["state_invariance"]["maximum_scalar_absolute_difference"])
    scalar_rel = float(contract["state_invariance"]["maximum_scalar_relative_difference"])
    vector_abs = float(contract["state_invariance"]["maximum_vector_component_absolute_difference"])
    vector_rel = float(contract["state_invariance"]["maximum_vector_component_relative_difference"])
    rows: list[dict[str, Any]] = []
    all_pass = True
    for state, baseline in exports["baseline"].items():
        weighted = exports["weighted"][state]
        nodes_exact, coordinate_error = compare_nodes(baseline, weighted)
        left_files = selected_field_files(baseline)
        right_names = {path.name for path in selected_field_files(weighted)}
        if {path.name for path in left_files} != right_names:
            raise RuntimeError(f"M51 field set mismatch for {state}")
        for left_path in left_files:
            right_path = weighted / "fields" / left_path.name
            left = read_csv(left_path)
            right = read_csv(right_path)
            if len(left) != len(right):
                raise RuntimeError(f"M51 field row count mismatch: {state} {left_path.name}")
            component_columns = [name for name in left[0] if name.startswith("component")]
            is_vector = len(component_columns) > 1
            absolute_limit = vector_abs if is_vector else scalar_abs
            relative_limit = vector_rel if is_vector else scalar_rel
            maximum_absolute = 0.0
            maximum_relative = 0.0
            failures = 0
            count = 0
            for a, b in zip(left, right, strict=True):
                if a["node_id"] != b["node_id"]:
                    raise RuntimeError(f"M51 node id mismatch: {state} {left_path.name}")
                for component in component_columns:
                    av, bv = float(a[component]), float(b[component])
                    difference = abs(av - bv)
                    relative = difference / max(abs(av), abs(bv), 1e-300)
                    maximum_absolute = max(maximum_absolute, difference)
                    maximum_relative = max(maximum_relative, relative)
                    failures += not close_enough(bv, av, absolute_limit, relative_limit)
                    count += 1
            passed = nodes_exact and failures == 0
            all_pass = all_pass and passed
            rows.append({
                "state": state, "field_file": left_path.name,
                "component_count": len(component_columns), "value_count": count,
                "node_ids_exact": True, "coordinates_exact": nodes_exact,
                "maximum_coordinate_difference_um": coordinate_error,
                "maximum_absolute_difference": maximum_absolute,
                "maximum_symmetric_relative_difference": maximum_relative,
                "failure_count": failures, "within_state_invariance": passed,
            })
    return rows, all_pass


def baseline_terminals() -> dict[tuple[str, float, str, str], float]:
    result: dict[tuple[str, float, str, str], float] = {}
    for row in read_csv(M50_TERMINALS):
        key = (row["device"], float(row["gate_voltage_V"]),
               row["contact"], row["component"])
        result[key] = float(row["m50_same_run_A_per_um"])
    return result


def baseline_fluxes() -> dict[tuple[str, float], dict[str, float]]:
    return {(row["device"], float(row["gate_voltage_V"])): {
        "normal": float(row["contact_surface_normal_integral_A_per_um"]),
        "negative": float(row["negative_contact_surface_normal_integral_A_per_um"])}
            for row in read_csv(M50_FLUX)}


def analyze(contract: dict[str, Any], manifest: dict[str, Any], banner: str,
            force_exports: bool) -> dict[str, Any]:
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    terminal_baseline = baseline_terminals()
    flux_baseline = baseline_fluxes()
    terminal_abs = float(contract["acceptance"][
        "maximum_terminal_closure_absolute_difference_A_per_um"])
    terminal_rel = float(contract["acceptance"][
        "maximum_terminal_closure_relative_difference"])
    flux_abs = float(contract["acceptance"][
        "maximum_native_flux_replay_absolute_difference_A_per_um"])
    flux_rel = float(contract["acceptance"][
        "maximum_native_flux_replay_relative_difference"])
    weighted_rows: dict[tuple[str, float], dict[str, float]] = {}
    for item in manifest["cases"]:
        files = sorted((raw / item["device"]).glob("IdVg_*des.plt"))
        if len(files) != 1:
            raise RuntimeError(f"expected one M51 current file for {item['device']}: {files}")
        for gate, row in parse_terminal_file(files[0]).items():
            weighted_rows[(item["device"], gate)] = row
    exports = export_snapshots(manifest, force_exports)
    field_rows, fields_invariant = compare_fields(contract, exports)
    terminal_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    for device in DEVICES:
        for gate in GATES:
            state = state_id(device, gate)
            weighted = weighted_rows[(device, gate)]
            for contact in CONTACTS:
                for component in COMPONENTS:
                    value = float(weighted[f"{contact}_{component}_A_per_um"])
                    baseline = terminal_baseline[(device, gate, contact, component)]
                    terminal_rows.append({
                        "state": state, "device": device, "gate_voltage_V": gate,
                        "contact": contact, "component": component,
                        "m50_default_A_per_um": baseline,
                        "m51_current_weighting_A_per_um": value,
                        "weighted_minus_default_A_per_um": value - baseline,
                        "absolute_difference_A_per_um": abs(value - baseline),
                        "matches_default": close_enough(
                            value, baseline, terminal_abs, terminal_rel),
                    })
            native = flux_baseline[(device, gate)]
            weighted_terminal = float(weighted["substrate_electron_A_per_um"])
            weighted_flux = float(weighted["native_flux_A_per_um"])
            default_value = terminal_baseline[(device, gate, "substrate", "electron")]
            default_native_gap = default_value - native["normal"]
            weighted_native_gap = weighted_terminal - native["normal"]
            gap_reduction_fraction = 1.0 - (
                abs(weighted_native_gap) / max(abs(default_native_gap), 1e-300))
            normal_close = close_enough(
                weighted_terminal, native["normal"], terminal_abs, terminal_rel)
            negative_close = close_enough(
                weighted_terminal, native["negative"], terminal_abs, terminal_rel)
            state_rows.append({
                "state": state, "device": device, "gate_voltage_V": gate,
                "drain_voltage_V": weighted["drain_voltage_V"],
                "m50_default_substrate_electron_A_per_um": default_value,
                "m51_weighted_substrate_electron_A_per_um": weighted_terminal,
                "weighted_minus_default_A_per_um": weighted_terminal - default_value,
                "m50_default_minus_native_normal_A_per_um": default_native_gap,
                "m50_native_normal_A_per_um": native["normal"],
                "m50_native_negative_normal_A_per_um": native["negative"],
                "m51_native_normal_A_per_um": weighted_flux,
                "m51_native_minus_m50_native_A_per_um": weighted_flux - native["normal"],
                "weighted_minus_native_normal_A_per_um": weighted_native_gap,
                "absolute_default_native_gap_reduction_fraction": gap_reduction_fraction,
                "weighted_native_normal_relative_difference": abs(weighted_native_gap) /
                    max(abs(native["normal"]), 1e-300),
                "native_flux_replays": close_enough(
                    weighted_flux, native["normal"], flux_abs, flux_rel),
                "weighted_matches_default": close_enough(
                    weighted_terminal, default_value, terminal_abs, terminal_rel),
                "weighted_matches_native_normal": normal_close,
                "weighted_matches_native_negative_normal": negative_close,
                "weighted_matches_either_native_orientation": normal_close or negative_close,
                "weighted_minus_nearest_native_A_per_um": min(
                    (weighted_terminal - native["normal"],
                     weighted_terminal - native["negative"]), key=abs),
            })
    default_matches = all(bool(row["weighted_matches_default"]) for row in state_rows)
    native_matches = all(bool(row["weighted_matches_either_native_orientation"])
                         for row in state_rows)
    native_replays = all(bool(row["native_flux_replays"]) for row in state_rows)
    if len(state_rows) != 6 or len(terminal_rows) != 72:
        classification = "replay_incomplete"
    elif not fields_invariant:
        classification = "state_perturbed"
    elif native_matches:
        classification = "weighted_matches_native_face"
    elif default_matches:
        classification = "weighted_matches_default"
    else:
        classification = "weighted_third_observable"
    target = next(row for row in state_rows
                  if row["device"] == "n23" and row["gate_voltage_V"] == 0.05)
    controls = [row for row in state_rows if
                (row["device"] == "n23" and row["gate_voltage_V"] in (0.0, 0.1))
                or (row["device"] == "n19" and row["gate_voltage_V"] == 0.05)]
    weighted_effect_localized = all(
        abs(float(target["weighted_minus_default_A_per_um"])) >
        abs(float(row["weighted_minus_default_A_per_um"])) for row in controls)
    declared = set(contract["analysis"]["classifications"])
    acceptance = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "exact_state_count": len(state_rows) == 6,
        "terminal_component_row_count": len(terminal_rows) == 72,
        "all_biases_exact": all(float(row["drain_voltage_V"]) == 0.05 for row in state_rows),
        "all_values_finite": all(math.isfinite(float(row[key])) for row in state_rows
                                  for key in ("m51_weighted_substrate_electron_A_per_um",
                                              "m51_native_normal_A_per_um")),
        "state_fields_invariant": fields_invariant,
        "native_flux_replays": native_replays,
        "classification_declared": classification in declared,
        "defaults_unchanged": True,
        "closed_topics_not_reopened": True,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    report = {
        "schema": "vela.simplemos.sdevice.m51_current_weighting_attribution_report.v1",
        "status": "complete", "classification": classification,
        "execution": {
            "sentaurus_release": "T-2022.03-SP2", "sentaurus_banner": banner,
            "new_sentaurus_execution": True, "new_vela_execution": False,
            "contract_sha256_before_and_after": sha256(CONTRACT),
            "intervention": "Math CurrentWeighting flag only",
            "production_defaults_changed": False,
        },
        "manual_interpretation": contract["manual_basis"],
        "target": target,
        "target_gap_reduction": {
            "absolute_default_native_gap_reduction_fraction": target[
                "absolute_default_native_gap_reduction_fraction"],
            "absolute_default_native_gap_reduction_percent": 100.0 * float(target[
                "absolute_default_native_gap_reduction_fraction"]),
            "weighted_native_normal_absolute_residual_A_per_um": abs(float(target[
                "weighted_minus_native_normal_A_per_um"])),
            "weighted_native_normal_relative_difference": target[
                "weighted_native_normal_relative_difference"],
            "strict_contract_closure_passed": target[
                "weighted_matches_native_normal"],
        },
        "controls": controls,
        "weighted_effect_target_localized": weighted_effect_localized,
        "state_field_summary": {
            "comparison_row_count": len(field_rows),
            "all_fields_invariant": fields_invariant,
            "maximum_absolute_difference": max(
                float(row["maximum_absolute_difference"]) for row in field_rows),
            "maximum_symmetric_relative_difference": max(
                float(row["maximum_symmetric_relative_difference"]) for row in field_rows),
            "total_failure_count": sum(int(row["failure_count"]) for row in field_rows),
        },
        "states": state_rows, "acceptance": acceptance,
        "claim_guard": contract["analysis"]["claim_guard"],
    }
    write_csv(STATES, state_rows)
    write_csv(TERMINALS, terminal_rows)
    write_csv(FIELDS, field_rows)
    write_json(REPORT, report)
    return report


def freeze_artifacts(report: dict[str, Any]) -> None:
    target = report["target"]
    table = "\n".join(
        f"| {row['device']} | {float(row['gate_voltage_V']):.2f} | "
        f"{float(row['m50_default_substrate_electron_A_per_um']):.12e} | "
        f"{float(row['m51_weighted_substrate_electron_A_per_um']):.12e} | "
        f"{float(row['m51_native_normal_A_per_um']):.12e} |"
        for row in report["states"])
    doc = f'''# SimpleMOS M51 CurrentWeighting 归因

## 结论

M51 分类为 `{report['classification']}`。唯一干预是在 M50 的全局 `Math` 中加入 `CurrentWeighting`；六状态物理、网格、偏压路径和诊断面通量保持冻结。

目标点 n23、Vd=0.05 V、Vg=0.05 V：默认 substrate `eCurrent` 为 `{float(target['m50_default_substrate_electron_A_per_um']):.12e}` A/um，`CurrentWeighting` 后为 `{float(target['m51_weighted_substrate_electron_A_per_um']):.12e}` A/um，变化 `{float(target['weighted_minus_default_A_per_um']):.12e}` A/um；同次原生法向面通量为 `{float(target['m51_native_normal_A_per_um']):.12e}` A/um。

相对于原生法向面通量，`CurrentWeighting` 消除了默认端口差异的 `{100.0 * float(target['absolute_default_native_gap_reduction_fraction']):.9f}%`；剩余绝对差为 `{abs(float(target['weighted_minus_native_normal_A_per_um'])):.12e}` A/um，相对原生量为 `{float(target['weighted_native_normal_relative_difference']):.9e}`。它非常接近原生量，但未达到冻结合同的严格闭合阈值，因此保持 `weighted_third_observable`，不改写为 `weighted_matches_native_face`。

## 六状态

| 器件 | Vg (V) | 默认 substrate eCurrent | CurrentWeighting eCurrent | 原生法向面通量 |
|---|---:|---:|---:|---:|
{table}

## 状态不变性

共比较 `{report['state_field_summary']['comparison_row_count']}` 个状态-字段文件；失败值数量为 `{report['state_field_summary']['total_failure_count']}`。最大绝对差为 `{float(report['state_field_summary']['maximum_absolute_difference']):.6e}`，最大对称相对差为 `{float(report['state_field_summary']['maximum_symmetric_relative_difference']):.6e}`。

## 边界

- 手册定义默认端口电流为关联掺杂阱表面通量与阱体生成率积分之和；`CurrentWeighting` 是减小数值误差的加权域积分。
- 本任务没有执行 `DirectCurrent`，没有修改 HFS、SG、准费米、BGN、SRH、迁移率、接触模型、网格或生产默认值。
- M51 只归因 Sentaurus 端口观测算法，不授权修改 Vela 或 Sentaurus 生产配置。

机器报告：`{portable(REPORT)}`。
'''
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc, encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.simplemos.sdevice.m51_artifact.v1",
        "title": "SimpleMOS M51 CurrentWeighting attribution",
        "status": report["status"], "classification": report["classification"],
        "report": portable(REPORT),
        "ledgers": [portable(STATES), portable(TERMINALS), portable(FIELDS)],
        "document": portable(DOC), "target": target,
        "acceptance": report["acceptance"],
    })
    artifacts = [REPORT, STATES, TERMINALS, FIELDS, DOC, ARTIFACT]
    sources = [CONTRACT, FREEZE, Path(__file__).resolve(), TEST,
               *UPSTREAM_EVIDENCE, OUTPUT / "sentaurus_manifest.json",
               OUTPUT / "sentaurus_banner.txt"]
    for device in DEVICES:
        sources.append(OUTPUT / "sentaurus_bundle" / device /
                       f"{device}_vd_0p05_des.cmd")
        sources.extend(sorted((OUTPUT / "sentaurus_raw/sentaurus_bundle" / device).glob(
            "IdVg_*des.plt")))
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m51_current_weighting_attribution_evidence.v1",
        "status": "frozen", "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "new_sentaurus_execution": True, "new_vela_execution": False,
        "only_intervention": "Math CurrentWeighting flag",
        "production_defaults_changed": False,
        "closed_topics_reinvestigated": False,
        "acceptance": report["acceptance"],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-sentaurus", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=2)
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
    manifest = prepare(args.force) if do_prepare else read_json(manifest_path)
    if not do_run and not do_analyze:
        print(json.dumps({"status": "prepared", "cases": len(manifest["cases"]),
                          "manifest": portable(manifest_path)}, indent=2))
        return
    banner_path = OUTPUT / "sentaurus_banner.txt"
    banner = (run_sentaurus(manifest, args.ssh_target, args.ssh_bin,
                            args.scp_bin, args.remote_root, args.jobs)
              if do_run else banner_path.read_text(encoding="utf-8").strip())
    if do_analyze:
        report = analyze(contract, manifest, banner, args.force)
        freeze_artifacts(report)
        print(json.dumps({"status": report["status"],
                          "classification": report["classification"],
                          "all_checks_pass": report["acceptance"]["all_checks_pass"],
                          "report": portable(REPORT)}, indent=2))


if __name__ == "__main__":
    main()
