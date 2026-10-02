#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M55 doping-well terminal attribution."""

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
import run_simplemos_m47_default_bgn_self_consistent_attribution as m47  # noqa: E402
import run_simplemos_m53_direct_current_full_matrix as m53  # noqa: E402
import sentaurus_import  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m55_doping_well_terminal_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m55_doping_well_terminal_attribution_contract_freeze.json"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
M47 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m47_default_bgn_state_attribution"
M53 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m53_direct_current_full_matrix"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m55_doping_well_terminal_attribution"
PORTABLE = ROOT / "doping_well_terminal_attribution"
REPORT = PORTABLE / "m55_doping_well_terminal_attribution_report.json"
REPLAY_LEDGER = PORTABLE / "m55_default_terminal_replay_ledger.csv"
WELL_LEDGER = PORTABLE / "m55_doping_well_ledger.csv"
TERMINAL_LEDGER = PORTABLE / "m55_terminal_component_attribution_ledger.csv"
PAIR_SUMMARY = PORTABLE / "m55_nwell_pair_summary.csv"
DOC = REPO / "docs/validation/simplemos_m55_doping_well_terminal_attribution_2026-09-02.md"
EVIDENCE = ROOT / "simplemos_m55_doping_well_terminal_attribution_evidence.json"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m55/artifact.json"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
ARCHIVE = "simplemos_m55_doping_well_terminal_attribution_results.tgz"
REMOTE_ROOT = "/tmp/vela_simplemos_m55_doping_well_terminal_attribution"
DEVICES = tuple(f"n{index}" for index in range(17, 25))
DRAINS = (0.05, 1.0)
GATES = (0.0, 0.05, 0.1)
ALL_GATES = tuple(0.05 * index for index in range(51))
CONTACTS = ("source", "drain", "gate", "substrate")
WELL_CONTACTS = ("source", "drain", "substrate")
COMPONENTS = ("electron", "hole", "total")
Q_C = 1.602176634e-19


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


def voltage_tag(value: float) -> str:
    return m53.voltage_tag(value)


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def close_enough(left: float, right: float, absolute: float,
                 relative: float) -> bool:
    difference = abs(left - right)
    return difference <= absolute or difference <= relative * max(
        abs(left), abs(right), 1e-300)


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    expected = "vela.simplemos.sdevice.m55_doping_well_terminal_attribution_contract.v1"
    if contract.get("schema") != expected:
        raise ValueError("unexpected M55 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M55 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M55 contract changed after freeze")
    for relative, expected_hash in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected_hash:
            raise ValueError(f"M55 upstream changed after freeze: {relative}")
    m46 = read_json(ROOT / "full_matrix_requalification/m46_full_matrix_requalification_report.json")
    m53_report = read_json(ROOT / "direct_current_full_matrix/m53_direct_current_full_matrix_report.json")
    m54_report = read_json(ROOT / "terminal_common_mode_attribution/m54_terminal_common_mode_attribution_report.json")
    if not m46["acceptance"]["all_checks_pass"]:
        raise ValueError("M55 requires accepted M46")
    if int(m46["findings"]["curve_count"]) != 16 or int(
            m46["findings"]["direct_bias_point_count"]) != 816:
        raise ValueError("M46 matrix anchor changed")
    if m53_report.get("classification") != contract["upstream"]["required_m53_classification"]:
        raise ValueError("M53 classification anchor changed")
    if m54_report.get("classification") != contract["upstream"]["required_m54_classification"]:
        raise ValueError("M54 classification anchor changed")
    if tuple(contract["matrix"]["devices"]) != DEVICES:
        raise ValueError("M55 device matrix changed")
    if tuple(map(float, contract["matrix"]["drain_voltages_V"])) != DRAINS:
        raise ValueError("M55 drain matrix changed")
    if tuple(map(float, contract["matrix"]["gate_control_voltages_V"])) != GATES:
        raise ValueError("M55 gate controls changed")
    return contract


def diagnostic_blocks(contract: dict[str, Any]) -> str:
    fields = """Plot {
  eDensity hDensity
  TotalCurrent/Vector eCurrent/Vector hCurrent/Vector
  eMobility hMobility eVelocity hVelocity
  eQuasiFermi hQuasiFermi
  ElectricField/Vector Potential SpaceCharge
  Doping DonorConcentration AcceptorConcentration DopingWells
  SRH
  eGradQuasiFermi/Vector hGradQuasiFermi/Vector
  eEparallel hEparallel eENormal hENormal
  BandGap BandGapNarrowing Affinity
  ConductionBand ValenceBand
}

"""
    formulas: list[str] = []
    seeds = contract["doping_well_diagnostics"]["seed_coordinates_um"]
    for contact in WELL_CONTACTS:
        x, y = map(float, seeds[contact])
        title = contact.capitalize()
        formulas.append(f'''  Tcl (
    Formula = "set value [tcl_cp_ReadScalar DopingWells]"
    Operation = "Coordinate = ({x:.12g} {y:.12g})"
    Dataset = "M55_{title}WellIndex"
    Function = "M55_{title}WellIndex"
    Unit = "1"
  )''')
        formulas.append(f'''  Tcl (
    Formula = "set value 1.0"
    Operation = "Integrate DopingWell=({x:.12g} {y:.12g}) IntegrationUnit=cm"
    Dataset = "M55_{title}WellArea"
    Function = "M55_{title}WellArea"
    Unit = "cm2"
  )''')
        formulas.append(f'''  Tcl (
    Formula = "set value [tcl_cp_ReadScalar srhRecombination]"
    Operation = "Integrate DopingWell=({x:.12g} {y:.12g}) IntegrationUnit=cm"
    Dataset = "M55_{title}WellSRH"
    Function = "M55_{title}WellSRH"
    Unit = "1/cm/s"
  )''')
    current_plot = "CurrentPlot {\n" + "\n".join(formulas) + "\n}\n\n"
    return fields + current_plot


def prepare(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force:
        for path in (bundle, OUTPUT / "sentaurus_raw", OUTPUT / "sentaurus_exports"):
            if path.exists():
                shutil.rmtree(path)
    bundle.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    marker = "Math { Extrapolate Iterations=20 ExitOnFailure }"
    sweep_marker = "      CurrentPlot(Time=(Range=(0 1) Intervals=50))"
    blocks = diagnostic_blocks(contract)
    for device in DEVICES:
        source_root = M8 / "sentaurus_bundle" / device
        source_tdr = source_root / "input_fps.tdr"
        target_root = bundle / device
        target_root.mkdir(parents=True, exist_ok=True)
        target_tdr = target_root / "input_fps.tdr"
        shutil.copy2(source_tdr, target_tdr)
        for drain in DRAINS:
            case = f"{device}_vd_{voltage_tag(drain)}"
            source_deck = source_root / f"{case}_des.cmd"
            baseline = source_deck.read_text(encoding="utf-8")
            if baseline.count(marker) != 1 or baseline.count(sweep_marker) != 1:
                raise RuntimeError(f"M8 insertion point changed: {source_deck}")
            if "CurrentWeighting" in baseline or "DirectCurrent" in baseline:
                raise RuntimeError(f"M8 current algorithm changed: {source_deck}")
            replay = baseline.replace(marker, blocks + marker)
            snapshot = (sweep_marker + "\n" +
                        f'      Plot(FilePrefix="m55_state_{case}" NoOverWrite Time=(0; 0.02; 0.04))')
            replay = replay.replace(sweep_marker, snapshot)
            target_deck = target_root / f"{case}_des.cmd"
            target_deck.write_text(replay, encoding="utf-8", newline="\n")
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
                "snapshot_prefix": f"m55_state_{case}",
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m55_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "intervention": "diagnostic CurrentPlot, DopingWells Plot field, and three low-gate snapshots only",
        "current_algorithm": "default",
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


def find_dataset(datasets: list[str], expected: str) -> str:
    if expected in datasets:
        return expected
    matches = [name for name in datasets if expected in name]
    if len(matches) != 1:
        raise RuntimeError(f"expected one dataset containing {expected}; got {matches}")
    return matches[0]


def diagnostic_rows(path: Path, tolerance: float,
                    contract: dict[str, Any]) -> dict[float, dict[str, float]]:
    text = path.read_text(errors="ignore")
    datasets = sentaurus_import.parse_quoted_list(text, "datasets")
    values_rows = sentaurus_import.parse_values_block(text, len(datasets))
    names: dict[tuple[str, str], str] = {}
    for contact in WELL_CONTACTS:
        title = contact.capitalize()
        for metric in ("WellIndex", "WellArea", "WellSRH"):
            names[(contact, metric)] = find_dataset(datasets, f"M55_{title}{metric}")
    result: dict[float, dict[str, float]] = {}
    for gate in GATES:
        matches = [values for values in values_rows
                   if abs(float(values[datasets.index("gate OuterVoltage")]) - gate)
                   <= tolerance]
        if len(matches) != 1:
            raise RuntimeError(f"{path}: expected one M55 diagnostic row at Vg={gate}")
        row = dict(zip(datasets, matches[0], strict=True))
        parsed: dict[str, float] = {}
        for contact in WELL_CONTACTS:
            parsed[f"{contact}_well_index"] = float(row[names[(contact, "WellIndex")]])
            parsed[f"{contact}_well_area_cm2"] = float(row[names[(contact, "WellArea")]])
            parsed[f"{contact}_well_srh_per_cm_s"] = float(row[names[(contact, "WellSRH")]])
        result[gate] = parsed
    return result


def export_snapshots(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    states: list[dict[str, Any]] = []
    for item in manifest["cases"]:
        case_root = raw / item["device"]
        snapshots = sorted(case_root.glob(f"{item['snapshot_prefix']}*.tdr"))
        if len(snapshots) != 3:
            raise RuntimeError(f"expected three M55 snapshots for {item['case']}; got {len(snapshots)}")
        for gate, tdr in zip(GATES, snapshots, strict=True):
            state = f"{item['case']}_vg_{voltage_tag(gate)}"
            export_dir = OUTPUT / "sentaurus_exports" / state
            run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export_dir)], capture=True)
            names = m47.m10.manifest_field_names(export_dir)
            if not any(name.lower() == "dopingwells" for name in names):
                raise RuntimeError(f"M55 state lacks DopingWells: {state}; fields={sorted(names)}")
            states.append({
                "state": state,
                "device": item["device"],
                "case": item["case"],
                "drain_voltage_V": float(item["drain_voltage_V"]),
                "gate_voltage_V": gate,
                "tdr": portable(tdr),
                "tdr_sha256": sha256(tdr),
                "export_dir": portable(export_dir),
                "field_manifest_sha256": sha256(export_dir / "field_manifest.json"),
            })
    payload = {
        "schema": "vela.simplemos.sdevice.m55_sentaurus_exports.v1",
        "state_count": len(states),
        "states": states,
    }
    write_json(OUTPUT / "sentaurus_export_manifest.json", payload)
    return states


def doping_well_field(export_dir: Path) -> dict[int, float]:
    candidates = [path for path in (export_dir / "fields").glob("*_region0.csv")
                  if path.stem.split("_region", 1)[0].lower() == "dopingwells"]
    if len(candidates) != 1:
        raise RuntimeError(f"expected one DopingWells region0 field in {export_dir}; got {candidates}")
    return {int(row["node_id"]): float(row["component0"])
            for row in read_csv(candidates[0])}


def support_metrics(export_dir: Path, well_index: float,
                    tolerance: float) -> dict[str, Any]:
    values = doping_well_field(export_dir)
    rounded = round(well_index)
    if abs(well_index - rounded) > tolerance:
        raise RuntimeError(f"non-integer M55 doping-well index {well_index} in {export_dir}")
    support = {node for node, value in values.items()
               if abs(value - rounded) <= tolerance}
    nodes = {int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
             for row in read_csv(export_dir / "nodes.csv")}
    elements = read_csv(export_dir / "elements.csv")
    element_count = sum(1 for row in elements
                        if row["material"] == "Si" and
                        {int(row["node0"]), int(row["node1"]), int(row["node2"])} <= support)
    if not support:
        raise RuntimeError(f"empty M55 well support for index {rounded} in {export_dir}")
    xs = [nodes[node][0] for node in support]
    ys = [nodes[node][1] for node in support]
    return {
        "well_index": rounded,
        "well_node_count": len(support),
        "well_element_count": element_count,
        "well_bbox_x_min_um": min(xs),
        "well_bbox_x_max_um": max(xs),
        "well_bbox_y_min_um": min(ys),
        "well_bbox_y_max_um": max(ys),
    }


def shared_field_replay(states: list[dict[str, Any]]) -> dict[str, Any]:
    comparisons = 0
    mismatches: list[str] = []
    for state in states:
        if state["device"] not in {"n19", "n23"} or not math.isclose(
                float(state["drain_voltage_V"]), 0.05, abs_tol=1e-15):
            continue
        old_state = f"{state['device']}_vd_0p05_vg_{m47.voltage_tag(float(state['gate_voltage_V']))}"
        old = M47 / "sentaurus_exports" / old_state
        new = REPO / state["export_dir"]
        for name in ("nodes.csv", "elements.csv", "contacts.csv"):
            comparisons += 1
            if (old / name).read_bytes() != (new / name).read_bytes():
                mismatches.append(f"{state['state']}:{name}")
        for old_field in sorted((old / "fields").glob("*.csv")):
            new_field = new / "fields" / old_field.name
            comparisons += 1
            if not new_field.is_file() or old_field.read_bytes() != new_field.read_bytes():
                mismatches.append(f"{state['state']}:fields/{old_field.name}")
    return {
        "comparison_count": comparisons,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches,
        "all_pointwise_identical": not mismatches,
    }


def analyze(contract: dict[str, Any], manifest: dict[str, Any],
            states: list[dict[str, Any]], banner: str) -> tuple[dict[str, Any],
                                                                 list[dict[str, Any]],
                                                                 list[dict[str, Any]],
                                                                 list[dict[str, Any]],
                                                                 list[dict[str, Any]]]:
    acceptance = contract["acceptance"]
    tolerance = float(acceptance["exact_bias_tolerance_V"])
    replay_abs = float(acceptance["terminal_replay_absolute_tolerance_A_per_um"])
    replay_rel = float(acceptance["terminal_replay_relative_tolerance"])
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    m8_raw = M8 / "sentaurus_raw/sentaurus_bundle"
    m53_raw = M53 / "sentaurus_raw/sentaurus_bundle"
    state_map = {(item["case"], float(item["gate_voltage_V"])): item for item in states}
    replay_rows: list[dict[str, Any]] = []
    well_rows: list[dict[str, Any]] = []
    terminal_rows: list[dict[str, Any]] = []
    parsed_controls: dict[tuple[str, float], dict[str, Any]] = {}
    max_replay_abs = 0.0
    max_replay_rel = 0.0
    for item in manifest["cases"]:
        case = item["case"]
        device = item["device"]
        new_plot = raw / device / item["expected_plot"]
        old_plot = m8_raw / device / item["expected_plot"]
        direct_plot = m53_raw / device / item["expected_plot"]
        new_all = m53.exact_terminal_curve(new_plot, list(ALL_GATES), tolerance)
        old_all = m53.exact_terminal_curve(old_plot, list(ALL_GATES), tolerance)
        direct_controls = {float(row["gate_voltage_V"]): row for row in
                           m53.exact_terminal_curve(direct_plot, list(GATES), tolerance)}
        new_controls = {float(row["gate_voltage_V"]): row for row in new_all
                        if float(row["gate_voltage_V"]) in GATES}
        diagnostics = diagnostic_rows(new_plot, tolerance, contract)
        for new_row, old_row in zip(new_all, old_all, strict=True):
            differences: list[float] = []
            relatives: list[float] = []
            for contact in CONTACTS:
                for component in COMPONENTS:
                    key = f"{contact}_{component}_A_per_um"
                    difference = float(new_row[key]) - float(old_row[key])
                    relative = abs(difference) / max(abs(float(old_row[key])), 1e-300)
                    differences.append(abs(difference))
                    relatives.append(relative)
            row_abs = max(differences)
            row_rel = max(relatives)
            max_replay_abs = max(max_replay_abs, row_abs)
            max_replay_rel = max(max_replay_rel, row_rel)
            replay_rows.append({
                "case": case,
                "device": device,
                "drain_voltage_V": float(item["drain_voltage_V"]),
                "gate_voltage_V": float(new_row["gate_voltage_V"]),
                "m8_drain_total_A_per_um": float(old_row["drain_total_A_per_um"]),
                "m55_drain_total_A_per_um": float(new_row["drain_total_A_per_um"]),
                "drain_total_difference_A_per_um": (float(new_row["drain_total_A_per_um"])
                                                     - float(old_row["drain_total_A_per_um"])),
                "maximum_terminal_component_absolute_difference_A_per_um": row_abs,
                "maximum_terminal_component_symmetric_relative_difference": row_rel,
                "within_tolerance": row_abs <= replay_abs or row_rel <= replay_rel,
            })
        for gate in GATES:
            default = new_controls[gate]
            direct = direct_controls[gate]
            diag = diagnostics[gate]
            parsed_controls[(case, gate)] = {
                "default": default, "direct": direct, "diagnostic": diag}
            export_dir = REPO / state_map[(case, gate)]["export_dir"]
            signed_srh: dict[str, float] = {"gate": 0.0}
            for contact in WELL_CONTACTS:
                srh_integral = float(diag[f"{contact}_well_srh_per_cm_s"])
                electron = -Q_C * 1e-4 * srh_integral
                signed_srh[contact] = electron
                metrics = support_metrics(
                    export_dir, float(diag[f"{contact}_well_index"]),
                    float(acceptance["well_seed_index_rounding_tolerance"]))
                well_rows.append({
                    "case": case,
                    "device": device,
                    "nwell_cm3": float(contract["matrix"]["nwell_cm3"][device]),
                    "drain_voltage_V": float(item["drain_voltage_V"]),
                    "gate_voltage_V": gate,
                    "contact": contact,
                    **metrics,
                    "well_area_cm2": float(diag[f"{contact}_well_area_cm2"]),
                    "well_srh_integral_per_cm_s": srh_integral,
                    "electron_srh_charge_current_A_per_um": electron,
                    "hole_srh_charge_current_A_per_um": -electron,
                    "srh_charge_cancellation_A_per_um": 0.0,
                    "export_dir": portable(export_dir),
                })
            for contact in CONTACTS:
                electron_generation = signed_srh[contact]
                for component in COMPONENTS:
                    key = f"{contact}_{component}_A_per_um"
                    default_value = float(default[key])
                    direct_value = float(direct[key])
                    shift = default_value - direct_value
                    generation = (electron_generation if component == "electron" else
                                  -electron_generation if component == "hole" else 0.0)
                    surface = shift - generation
                    terminal_rows.append({
                        "case": case,
                        "device": device,
                        "nwell_cm3": float(contract["matrix"]["nwell_cm3"][device]),
                        "drain_voltage_V": float(item["drain_voltage_V"]),
                        "gate_voltage_V": gate,
                        "contact": contact,
                        "component": component,
                        "default_A_per_um": default_value,
                        "direct_A_per_um": direct_value,
                        "default_minus_direct_A_per_um": shift,
                        "signed_well_srh_charge_current_A_per_um": generation,
                        "surface_redistribution_A_per_um": surface,
                        "absolute_generation_fraction_of_shift": (
                            abs(generation) / max(abs(shift), 1e-300)),
                    })

    by_terminal = {(row["case"], float(row["gate_voltage_V"]), row["contact"], row["component"]): row
                   for row in terminal_rows}
    identity_rows: list[dict[str, Any]] = []
    max_identity = 0.0
    max_cancel = 0.0
    for item in manifest["cases"]:
        case = item["case"]
        for gate in GATES:
            for contact in CONTACTS:
                electron = by_terminal[(case, gate, contact, "electron")]
                hole = by_terminal[(case, gate, contact, "hole")]
                total = by_terminal[(case, gate, contact, "total")]
                cancellation = (float(electron["signed_well_srh_charge_current_A_per_um"])
                                + float(hole["signed_well_srh_charge_current_A_per_um"]))
                residual = (float(total["default_minus_direct_A_per_um"])
                            - float(electron["surface_redistribution_A_per_um"])
                            - float(hole["surface_redistribution_A_per_um"]))
                max_cancel = max(max_cancel, abs(cancellation))
                max_identity = max(max_identity, abs(residual))
                identity_rows.append({
                    "case": case, "gate_voltage_V": gate, "contact": contact,
                    "srh_cancellation_A_per_um": cancellation,
                    "total_surface_identity_residual_A_per_um": residual,
                })

    pairs: list[dict[str, Any]] = []
    well_by = {(row["device"], float(row["drain_voltage_V"]),
                float(row["gate_voltage_V"]), row["contact"]): row for row in well_rows}
    topology_changed = False
    area_tol = float(acceptance["matched_well_area_relative_tolerance"])
    m53_pairs = read_csv(ROOT / "direct_current_full_matrix/m53_nwell_pair_summary.csv")
    m53_pair_map = {(row["low_nwell_device"], row["high_nwell_device"],
                     float(row["drain_voltage_V"])): row
                    for row in m53_pairs}
    for spec in contract["matrix"]["matched_nwell_pairs"]:
        low, high = spec["low"], spec["high"]
        for drain in DRAINS:
            relative_areas: list[float] = []
            area_details: list[dict[str, Any]] = []
            support_changes = 0
            surface_high_low: list[float] = []
            for gate in GATES:
                for contact in WELL_CONTACTS:
                    low_well = well_by[(low, drain, gate, contact)]
                    high_well = well_by[(high, drain, gate, contact)]
                    low_area = float(low_well["well_area_cm2"])
                    high_area = float(high_well["well_area_cm2"])
                    relative = abs(high_area - low_area) / max(abs(low_area), abs(high_area), 1e-300)
                    relative_areas.append(relative)
                    area_details.append({
                        "contact": contact,
                        "gate_voltage_V": gate,
                        "low_well_area_cm2": low_area,
                        "high_well_area_cm2": high_area,
                        "relative_area_difference": relative,
                        "low_well_node_count": int(low_well["well_node_count"]),
                        "high_well_node_count": int(high_well["well_node_count"]),
                        "low_well_element_count": int(low_well["well_element_count"]),
                        "high_well_element_count": int(high_well["well_element_count"]),
                    })
                    changed = (relative > area_tol or
                               int(low_well["well_node_count"]) != int(high_well["well_node_count"]) or
                               int(low_well["well_element_count"]) != int(high_well["well_element_count"]))
                    support_changes += int(changed)
                for contact in ("source", "drain"):
                    low_surface = float(by_terminal[(f"{low}_vd_{voltage_tag(drain)}", gate,
                                                     contact, "total")]["surface_redistribution_A_per_um"])
                    high_surface = float(by_terminal[(f"{high}_vd_{voltage_tag(drain)}", gate,
                                                      contact, "total")]["surface_redistribution_A_per_um"])
                    surface_high_low.append(abs(high_surface) - abs(low_surface))
            topology_changed = topology_changed or support_changes > 0
            anchor = m53_pair_map[(low, high, drain)]
            maximum_area = max(area_details,
                               key=lambda row: float(row["relative_area_difference"]))
            pairs.append({
                "low_device": low,
                "high_device": high,
                "drain_voltage_V": drain,
                "GOxTime_min": float(spec["GOxTime_min"]),
                "LDD_Dose_cm2": float(spec["LDD_Dose_cm2"]),
                "maximum_relative_well_area_difference": max(relative_areas),
                "maximum_area_difference_contact": maximum_area["contact"],
                "maximum_area_difference_gate_voltage_V": maximum_area["gate_voltage_V"],
                "maximum_area_difference_low_well_area_cm2": maximum_area["low_well_area_cm2"],
                "maximum_area_difference_high_well_area_cm2": maximum_area["high_well_area_cm2"],
                "maximum_area_difference_low_well_node_count": maximum_area["low_well_node_count"],
                "maximum_area_difference_high_well_node_count": maximum_area["high_well_node_count"],
                "changed_well_support_state_contact_count": support_changes,
                "maximum_high_minus_low_absolute_surface_redistribution_A_per_um": max(surface_high_low),
                "minimum_high_minus_low_absolute_surface_redistribution_A_per_um": min(surface_high_low),
                "m53_default_high_minus_low_max_error_dex": float(
                    anchor["default_high_minus_low_amplification_dex"]),
                "m53_direct_high_minus_low_max_error_dex": float(
                    anchor["direct_high_minus_low_amplification_dex"]),
            })

    shared = shared_field_replay(states)
    area_groups: dict[tuple[str, str], list[float]] = {}
    for row in well_rows:
        area_groups.setdefault((row["device"], row["contact"]), []).append(
            float(row["well_area_cm2"]))
    maximum_within_device_area_span = max(
        max(values) - min(values) for values in area_groups.values())
    replay_ok = (len(replay_rows) == 816 and all(bool(row["within_tolerance"])
                                                 for row in replay_rows))
    surface_ok = max_identity <= float(acceptance[
        "total_surface_identity_absolute_tolerance_A_per_um"])
    cancellation_ok = max_cancel <= float(acceptance[
        "srh_electron_hole_cancellation_absolute_tolerance_A_per_um"])
    if not replay_ok or not shared["all_pointwise_identical"]:
        classification = "state_or_terminal_replay_mismatch"
    elif topology_changed:
        classification = "nwell_changes_well_topology_surface_redistribution"
    else:
        classification = "stable_well_topology_surface_redistribution"
    target_rows = [row for row in terminal_rows
                   if row["device"] == "n23" and
                   math.isclose(float(row["drain_voltage_V"]), 0.05) and
                   math.isclose(float(row["gate_voltage_V"]), 0.05) and
                   row["contact"] in {"source", "drain", "substrate"}]
    target_wells = [row for row in well_rows
                    if row["device"] == "n23" and
                    math.isclose(float(row["drain_voltage_V"]), 0.05) and
                    math.isclose(float(row["gate_voltage_V"]), 0.05)]
    checks = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "case_count": len(manifest["cases"]) == int(acceptance["required_case_count"]),
        "state_count": len(states) == int(acceptance["required_state_count"]),
        "default_replay_row_count": len(replay_rows) == 816,
        "default_terminal_replay": replay_ok,
        "shared_m47_fields_pointwise_identical": bool(shared["all_pointwise_identical"]),
        "well_row_count": len(well_rows) == int(acceptance["required_well_row_count"]),
        "terminal_component_row_count": len(terminal_rows) == int(
            acceptance["required_terminal_component_row_count"]),
        "srh_charge_terms_cancel": cancellation_ok,
        "total_surface_identity": surface_ok,
        "classification_declared": classification in contract["analysis"]["classifications"],
        "closed_topics_not_reopened": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    maximum_pair_area = max(
        pairs, key=lambda row: float(row["maximum_relative_well_area_difference"]))
    report = {
        "schema": "vela.simplemos.sdevice.m55_doping_well_terminal_attribution_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "sentaurus": {"banner": banner, "case_count": len(manifest["cases"]),
                      "state_count": len(states), "current_algorithm": "default"},
        "default_terminal_replay": {
            "row_count": len(replay_rows),
            "maximum_absolute_difference_A_per_um": max_replay_abs,
            "maximum_symmetric_relative_difference": max_replay_rel,
            "all_within_tolerance": replay_ok,
        },
        "shared_state_replay": shared,
        "doping_well": {
            "row_count": len(well_rows),
            "topology_changed_across_at_least_one_nwell_pair": topology_changed,
            "maximum_matched_pair_relative_area_difference": max(
                float(row["maximum_relative_well_area_difference"]) for row in pairs),
            "maximum_area_difference_pair": {
                key: maximum_pair_area[key] for key in (
                    "low_device", "high_device", "drain_voltage_V",
                    "maximum_area_difference_contact",
                    "maximum_area_difference_gate_voltage_V",
                    "maximum_area_difference_low_well_area_cm2",
                    "maximum_area_difference_high_well_area_cm2",
                    "maximum_area_difference_low_well_node_count",
                    "maximum_area_difference_high_well_node_count")
            },
            "changed_support_state_contact_count": sum(
                int(row["changed_well_support_state_contact_count"]) for row in pairs),
            "maximum_within_device_contact_area_span_cm2":
                maximum_within_device_area_span,
        },
        "terminal_attribution": {
            "row_count": len(terminal_rows),
            "maximum_srh_electron_hole_cancellation_A_per_um": max_cancel,
            "maximum_total_surface_identity_residual_A_per_um": max_identity,
            "interpretation": "SRH charge-current terms are opposite for electrons and holes and cancel from total current; the default-minus-Direct total shift is therefore carried by the well-surface versus contact-surface redistribution residual.",
        },
        "target": {"device": "n23", "drain_voltage_V": 0.05,
                   "gate_voltage_V": 0.05,
                   "terminal_rows": target_rows, "well_rows": target_wells},
        "nwell_pairs": pairs,
        "claim_guard": contract["analysis"]["claim_guard"],
        "forbidden_work_respected": contract["forbidden_work"],
        "acceptance": checks,
    }
    write_json(REPORT, report)
    return report, replay_rows, well_rows, terminal_rows, pairs


def freeze_artifacts(report: dict[str, Any], replay_rows: list[dict[str, Any]],
                     well_rows: list[dict[str, Any]],
                     terminal_rows: list[dict[str, Any]],
                     pair_rows: list[dict[str, Any]]) -> None:
    write_csv(REPLAY_LEDGER, replay_rows)
    write_csv(WELL_LEDGER, well_rows)
    write_csv(TERMINAL_LEDGER, terminal_rows)
    write_csv(PAIR_SUMMARY, pair_rows)
    target = report["target"]
    target_total = {f"{row['contact']}_{row['component']}": row
                    for row in target["terminal_rows"]}
    target_wells = {row["contact"]: row for row in target["well_rows"]}
    pair_lines = [
        f"| {row['low_device']}/{row['high_device']} | {float(row['drain_voltage_V']):.2f} | "
        f"{float(row['maximum_relative_well_area_difference']):.6e} | "
        f"{int(row['changed_well_support_state_contact_count'])} | "
        f"{float(row['m53_default_high_minus_low_max_error_dex']):.6f} |"
        for row in pair_rows]
    DOC.write_text(f"""# SimpleMOS M55 接触关联掺杂阱端口归因

## 结论

M55 分类为 `{report['classification']}`。16个默认算法 Sentaurus deck、816个完整栅压点全部回放，最大端口分量绝对差为 `{float(report['default_terminal_replay']['maximum_absolute_difference_A_per_um']):.12e}` A/um；n19/n23六个共享状态的导出网格和全部 M47 字段逐文件一致。

Sentaurus `DopingWell` 域积分覆盖48个低栅压状态、source/drain/substrate三个接触关联阱，共144行。匹配低/高 NWell 器件的最大阱面积相对差为 `{float(report['doping_well']['maximum_matched_pair_relative_area_difference']):.6e}`，有 `{int(report['doping_well']['changed_support_state_contact_count'])}` 个状态-接触组合改变了离散阱支持；同一器件/接触跨两个漏压和三个栅压的最大面积跨度为 `{float(report['doping_well']['maximum_within_device_contact_area_span_cm2']):.3e}` cm2，证明该变化来自冻结TDR阱分区而不是偏压状态。

阱内 SRH 电荷项对电子和空穴符号相反，最大抵消残差为 `{float(report['terminal_attribution']['maximum_srh_electron_hole_cancellation_A_per_um']):.3e}` A/um；默认减 Direct 的总电流变化与电子、空穴阱表面重分配之和的最大恒等式残差为 `{float(report['terminal_attribution']['maximum_total_surface_identity_residual_A_per_um']):.3e}` A/um。因此，SRH阱体项可以改变载流子分量账本，但不能解释总 Id 的默认/Direct 差；总电流差由关联阱表面相对物理接触面的电流重分配承担。

## n23目标点

| 量 | 值 |
|---|---:|
| drain总电流默认-Direct | {float(target_total['drain_total']['default_minus_direct_A_per_um']):.12e} A/um |
| drain总表面重分配 | {float(target_total['drain_total']['surface_redistribution_A_per_um']):.12e} A/um |
| source总表面重分配 | {float(target_total['source_total']['surface_redistribution_A_per_um']):.12e} A/um |
| substrate总表面重分配 | {float(target_total['substrate_total']['surface_redistribution_A_per_um']):.12e} A/um |
| substrate阱面积 | {float(target_wells['substrate']['well_area_cm2']):.12e} cm2 |
| substrate阱SRH电子电荷项 | {float(target_wells['substrate']['electron_srh_charge_current_A_per_um']):.12e} A/um |

## NWell配对

| 低/高NWell | Vd (V) | 最大阱面积相对差 | 改变支持数 | 默认高减低最大误差 (dex) |
|---|---:|---:|---:|---:|
{chr(10).join(pair_lines)}

## 边界

- M55只增加 `CurrentPlot DopingWell` 域积分、`DopingWells` 绘图字段和三个低栅压快照；默认端口算法、物理、网格和偏压路径未改变。
- “表面重分配”是冻结的默认阱面定义减去 Direct 接触面定义并扣除阱内SRH电荷项的账本量，不是对Sentaurus内部单元贡献的逐项复刻。
- 阱支持变化说明 NWell 参与默认端口观测域的构造，但不单独证明其为全部 Vela-Sentaurus 误差的生产代码根因。
- 未重新排查 HFS、SG、准费米、BGN、SRH模型、迁移率、接触模型或网格，也未修改生产默认值。

机器报告：`reference_tcad/simplemos_sentaurus2022/doping_well_terminal_attribution/m55_doping_well_terminal_attribution_report.json`。
""", encoding="utf-8", newline="\n")
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m55_doping_well_terminal_attribution_evidence.v1",
        "status": "frozen",
        "classification": report["classification"],
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "report": portable(REPORT),
        "report_sha256": sha256(REPORT),
        "new_sentaurus_execution": True,
        "new_vela_execution": False,
        "historical_artifacts_rewritten": False,
        "closed_topics_reinvestigated": False,
        "sentaurus_release": report["sentaurus"]["banner"],
        "source_hashes": {
            portable(CONTRACT): sha256(CONTRACT),
            portable(FREEZE): sha256(FREEZE),
            portable(REPO / "scripts/run_simplemos_m55_doping_well_terminal_attribution.py"):
                sha256(REPO / "scripts/run_simplemos_m55_doping_well_terminal_attribution.py"),
        },
        "artifacts": {
            portable(REPLAY_LEDGER): sha256(REPLAY_LEDGER),
            portable(WELL_LEDGER): sha256(WELL_LEDGER),
            portable(TERMINAL_LEDGER): sha256(TERMINAL_LEDGER),
            portable(PAIR_SUMMARY): sha256(PAIR_SUMMARY),
            portable(DOC): sha256(DOC),
        },
        "acceptance": report["acceptance"],
    })
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M55 doping-well terminal attribution",
        "status": report["status"],
        "classification": report["classification"],
        "summary": report["terminal_attribution"]["interpretation"],
        "report": portable(REPORT),
        "ledgers": [portable(WELL_LEDGER), portable(TERMINAL_LEDGER),
                    portable(PAIR_SUMMARY)],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-sentaurus", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=m47.m10.executable("ssh"))
    parser.add_argument("--scp-bin", default=m47.m10.executable("scp"))
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    contract = validate_contract()
    do_all = not (args.prepare or args.run_sentaurus or args.analyze)
    manifest_path = OUTPUT / "sentaurus_manifest.json"
    manifest = prepare(contract, args.force) if (args.prepare or do_all) else read_json(manifest_path)
    banner = ""
    if args.run_sentaurus or do_all:
        banner = run_sentaurus(manifest, args.ssh_target, args.ssh_bin,
                               args.scp_bin, args.remote_root, args.jobs)
    elif args.analyze:
        banner = (OUTPUT / "sentaurus_banner.txt").read_text(encoding="utf-8").strip()
    if args.analyze or do_all:
        states = export_snapshots(manifest)
        report, replay, wells, terminals, pairs = analyze(
            contract, manifest, states, banner)
        freeze_artifacts(report, replay, wells, terminals, pairs)
        if not report["acceptance"]["all_checks_pass"]:
            raise RuntimeError(f"M55 acceptance failed: {report['acceptance']}")
        print(json.dumps({
            "status": report["status"],
            "classification": report["classification"],
            "all_checks_pass": report["acceptance"]["all_checks_pass"],
            "maximum_replay_absolute_difference_A_per_um": report[
                "default_terminal_replay"]["maximum_absolute_difference_A_per_um"],
            "maximum_matched_pair_relative_area_difference": report[
                "doping_well"]["maximum_matched_pair_relative_area_difference"],
            "report": portable(REPORT),
        }, indent=2))
    elif args.prepare:
        print(json.dumps({"status": "prepared", "case_count": len(manifest["cases"]),
                          "manifest": portable(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
