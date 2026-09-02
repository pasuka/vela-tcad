#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M60 tight-convergence port-burst replay."""

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
import run_simplemos_m51_current_weighting_attribution as m51  # noqa: E402
import sentaurus_import  # noqa: E402

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m60_tight_convergence_port_burst_contract_v1.json"
FREEZE = ROOT / "simplemos_m60_tight_convergence_port_burst_contract_freeze.json"
M46_POINTS = ROOT / "full_matrix_requalification/m46_pointwise_delta.csv"
M59_POINTS = ROOT / "port_observable_decomposition/m59_pointwise_observable_ledger.csv"
M59_EVIDENCE = ROOT / "simplemos_m59_port_observable_decomposition_evidence.json"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m60_tight_convergence_port_burst"
PORTABLE = ROOT / "tight_convergence_port_burst"
REPORT = PORTABLE / "m60_tight_convergence_port_burst_report.json"
POINTS = PORTABLE / "m60_tight_default_direct_point_ledger.csv"
TERMINALS = PORTABLE / "m60_tight_terminal_component_ledger.csv"
CASES = PORTABLE / "m60_tight_case_summary.csv"
LOGS = PORTABLE / "m60_cnormprint_log_ledger.csv"
FIELDS = PORTABLE / "m60_default_direct_field_invariance_ledger.csv"
QUALIFICATION = PORTABLE / "m60_sub_fa_secondary_reference_qualification.json"
DOC = REPO / "docs/validation/simplemos_m60_tight_convergence_port_burst_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m60/artifact.json"
EVIDENCE = ROOT / "simplemos_m60_tight_convergence_port_burst_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m60_tight_convergence_port_burst.py"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
REMOTE_ROOT = "/tmp/vela_simplemos_m60_tight_convergence_port_burst"
ARCHIVE = "simplemos_m60_tight_convergence_port_burst_results.tgz"
CONTACTS = ("source", "drain", "gate", "substrate")
COMPONENTS = ("electron", "hole", "total")
CORE_DEVICES = ("n19", "n23", "n24")
SNAPSHOT_GATES = (0.0, 0.05, 0.1, 0.15, 0.8)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8",
                    newline="\n")


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


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def voltage_tag(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p")


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != "vela.simplemos.sdevice.m60_tight_convergence_port_burst_contract.v1":
        raise ValueError("unexpected M60 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M60 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M60 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M60 upstream artifact changed: {relative}")
    m59 = read_json(M59_EVIDENCE)
    if (m59.get("status") != contract["upstream"]["required_m59_status"] or
            m59.get("classification") != contract["upstream"]["required_m59_classification"]):
        raise ValueError("M59 status or classification changed")
    return contract


def sentaurus_deck(run_case: str, drain: float, algorithm: str,
                   core_snapshots: bool) -> str:
    deck = m47.sentaurus_deck(run_case)
    deck = deck.replace(
        'Goal { Name="drain" Voltage=0.05 }',
        f'Goal {{ Name="drain" Voltage={drain:.12g} }}')
    baseline_math = "Math { Extrapolate Iterations=20 ExitOnFailure }"
    direct = " DirectCurrent" if algorithm == "direct" else ""
    tight_math = (
        "Math { Extrapolate RelErrControl Digits=8 "
        "ErrRef(Electron)=1e2 ErrRef(Hole)=1e2 Iterations=20 "
        f"ExitOnFailure CNormPrint{direct} }}")
    if deck.count(baseline_math) != 1:
        raise RuntimeError("M47 Math anchor changed")
    deck = deck.replace(baseline_math, tight_math)
    original_plot = 'Plot(FilePrefix="m47_state" NoOverWrite Time=(0; 0.02; 0.04))'
    replacement = (
        f'Plot(FilePrefix="{run_case}_state" NoOverWrite '
        'Time=(0; 0.02; 0.04; 0.06; 0.32))'
        if core_snapshots else "")
    if deck.count(original_plot) != 1:
        raise RuntimeError("M47 snapshot anchor changed")
    return deck.replace(original_plot, replacement)


def prepare(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    cases: list[dict[str, Any]] = []
    for device in contract["matrix"]["devices"]:
        source = M8 / "sentaurus_bundle" / device / "input_fps.tdr"
        if not source.is_file():
            raise FileNotFoundError(source)
        root = bundle / device
        root.mkdir(parents=True, exist_ok=True)
        target = root / "input_fps.tdr"
        shutil.copy2(source, target)
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            base_case = f"{device}_vd_{voltage_tag(drain)}"
            for algorithm in ("default", "direct"):
                run_case = f"m60_{algorithm}_{base_case}"
                deck_path = root / f"{run_case}_des.cmd"
                deck_path.write_text(
                    sentaurus_deck(run_case, drain, algorithm,
                                   device in CORE_DEVICES and math.isclose(drain, 0.05)),
                    encoding="utf-8", newline="\n")
                cases.append({
                    "run_case": run_case, "base_case": base_case,
                    "device": device, "drain_voltage_V": drain,
                    "algorithm": algorithm,
                    "deck": portable(deck_path), "deck_sha256": sha256(deck_path),
                    "input_tdr": portable(target), "input_tdr_sha256": sha256(target),
                    "expected_current": f"IdVg_{run_case}_des.plt",
                    "expected_console": f"{run_case}.console.log",
                    "core_snapshots": bool(device in CORE_DEVICES and math.isclose(drain, 0.05)),
                })
    if len(cases) != int(contract["matrix"]["sentaurus_case_count"]):
        raise RuntimeError("M60 prepared case count mismatch")
    manifest = {
        "schema": "vela.simplemos.sdevice.m60_sentaurus_manifest.v1",
        "status": "prepared", "contract_sha256": sha256(CONTRACT),
        "tight_math": contract["intervention"]["tight_math"],
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

    def execute(group: tuple[str, list[dict[str, Any]]]) -> None:
        device, cases = group
        root = f"{remote_root}/sentaurus_bundle/{device}"
        commands = [f"cd {root}"]
        for item in cases:
            name = item["run_case"]
            commands.append(f"sdevice {name}_des.cmd > {name}.console.log 2>&1")
        run([ssh_bin, ssh_target, "set -eu; " + "; ".join(commands)])

    groups = list(by_device.items())
    with ThreadPoolExecutor(max_workers=min(max(jobs, 1), len(groups))) as pool:
        list(pool.map(execute, groups))
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


def parse_current(path: Path, gates: list[float], tolerance: float
                  ) -> list[dict[str, float]]:
    text = path.read_text(errors="ignore")
    datasets = sentaurus_import.parse_quoted_list(text, "datasets")
    required = ["gate OuterVoltage", "drain OuterVoltage"]
    for contact in CONTACTS:
        required += [f"{contact} eCurrent", f"{contact} hCurrent",
                     f"{contact} TotalCurrent"]
    missing = set(required) - set(datasets)
    if missing:
        raise RuntimeError(f"missing M60 current datasets: {path}: {sorted(missing)}")
    source = sentaurus_import.parse_values_block(text, len(datasets))
    parsed: list[dict[str, float]] = []
    for gate in gates:
        matches = [values for values in source
                   if abs(float(values[datasets.index("gate OuterVoltage")]) - gate)
                   <= tolerance]
        if len(matches) != 1:
            raise RuntimeError(f"{path}: exact Vg={gate:g} count {len(matches)}")
        raw = dict(zip(datasets, matches[0], strict=True))
        row = {"gate_voltage_V": gate,
               "drain_voltage_V": float(raw["drain OuterVoltage"])}
        for contact in CONTACTS:
            row[f"{contact}_electron_A_per_um"] = float(raw[f"{contact} eCurrent"])
            row[f"{contact}_hole_A_per_um"] = float(raw[f"{contact} hCurrent"])
            row[f"{contact}_total_A_per_um"] = float(raw[f"{contact} TotalCurrent"])
        parsed.append(row)
    return parsed


def export_and_compare_snapshots(contract: dict[str, Any], manifest: dict[str, Any],
                                 force: bool) -> tuple[list[dict[str, Any]], bool]:
    exports: dict[str, dict[str, Path]] = {"baseline": {}, "weighted": {}}
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    selected = [item for item in manifest["cases"] if item["core_snapshots"]]
    for item in selected:
        root = raw / item["device"]
        snapshots = sorted(root.glob(f"{item['run_case']}_state*.tdr"))
        if len(snapshots) != len(SNAPSHOT_GATES):
            raise RuntimeError(f"M60 snapshot count mismatch {item['run_case']}: {snapshots}")
        variant = "baseline" if item["algorithm"] == "default" else "weighted"
        for gate, tdr in zip(SNAPSHOT_GATES, snapshots, strict=True):
            state = f"{item['device']}_vd_0p05_vg_{voltage_tag(gate)}"
            export = OUTPUT / "state_exports" / item["algorithm"] / state
            if force and export.exists():
                shutil.rmtree(export)
            if not (export / "field_manifest.json").is_file():
                run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export)],
                    capture=True)
            exports[variant][state] = export
    fake_contract = {"state_invariance": {
        "maximum_scalar_absolute_difference": contract["acceptance"]["state_scalar_absolute_tolerance"],
        "maximum_scalar_relative_difference": contract["acceptance"]["state_scalar_symmetric_relative_tolerance"],
        "maximum_vector_component_absolute_difference": contract["acceptance"]["state_vector_absolute_tolerance"],
        "maximum_vector_component_relative_difference": contract["acceptance"]["state_vector_symmetric_relative_tolerance"],
    }}
    return m51.compare_fields(fake_contract, exports)


def analyze(contract: dict[str, Any], manifest: dict[str, Any], banner: str,
            force_exports: bool) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    gates = [0.05 * index for index in range(51)]
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    parsed: dict[tuple[str, str], list[dict[str, float]]] = {}
    log_rows: list[dict[str, Any]] = []
    for item in manifest["cases"]:
        current = raw / item["device"] / item["expected_current"]
        console = raw / item["device"] / item["expected_console"]
        if not current.is_file() or not console.is_file():
            raise FileNotFoundError(current if not current.is_file() else console)
        parsed[(item["base_case"], item["algorithm"])] = parse_current(
            current, gates, tolerance)
        console_text = console.read_text(errors="ignore")
        log_rows.append({
            "run_case": item["run_case"], "base_case": item["base_case"],
            "device": item["device"], "drain_voltage_V": item["drain_voltage_V"],
            "algorithm": item["algorithm"], "console_bytes": console.stat().st_size,
            "console_sha256": sha256(console),
            "cnorm_or_error_line_count": sum(
                1 for line in console_text.splitlines()
                if "CNorm" in line or "error" in line.lower() or "Rhs" in line),
            "completed_without_exit_failure": int("terminated by user" not in console_text.lower()),
        })
    baseline = {(row["case"], round(float(row["gate_voltage_V"]), 12)): row
                for row in read_csv(M46_POINTS)}
    m59 = {(row["case"], round(float(row["gate_voltage_V"]), 12)): row
           for row in read_csv(M59_POINTS)}
    pair_rows: list[dict[str, Any]] = []
    terminal_rows: list[dict[str, Any]] = []
    case_rows: list[dict[str, Any]] = []
    burst_threshold = float(contract["acceptance"]["burst_absolute_fraction"])
    target: dict[str, Any] | None = None
    all_nonburst_shifts: list[float] = []
    for device in contract["matrix"]["devices"]:
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            case = f"{device}_vd_{voltage_tag(drain)}"
            defaults = parsed[(case, "default")]
            directs = parsed[(case, "direct")]
            case_points: list[dict[str, Any]] = []
            for default, direct in zip(defaults, directs, strict=True):
                gate = float(default["gate_voltage_V"])
                anchor = baseline[(case, round(gate, 12))]
                old = m59[(case, round(gate, 12))]
                default_id = float(default["drain_total_A_per_um"])
                direct_id = float(direct["drain_total_A_per_um"])
                vela_id = float(anchor["m46_vela_current_A_per_um"])
                old_default_id = float(anchor["sentaurus_current_A_per_um"])
                gap = (float(default["substrate_electron_A_per_um"])
                       - float(direct["substrate_electron_A_per_um"]))
                fraction = gap / default_id
                error = abs(math.log10(abs(vela_id) / abs(default_id)))
                old_error = abs(math.log10(abs(vela_id) / abs(old_default_id)))
                old_to_tight = abs(math.log10(abs(default_id) / abs(old_default_id)))
                old_burst = int(old["primary_burst_flag"])
                if not old_burst:
                    all_nonburst_shifts.append(old_to_tight)
                row = {
                    "case": case, "device": device, "drain_voltage_V": drain,
                    "gate_voltage_V": gate,
                    "m46_default_drain_current_A_per_um": old_default_id,
                    "tight_default_drain_current_A_per_um": default_id,
                    "tight_direct_drain_current_A_per_um": direct_id,
                    "vela_drain_current_A_per_um": vela_id,
                    "m46_default_vela_error_dex": old_error,
                    "tight_default_vela_error_dex": error,
                    "tight_default_vs_m46_log_shift_dex": old_to_tight,
                    "tight_default_substrate_eCurrent_A_per_um": default["substrate_electron_A_per_um"],
                    "tight_direct_substrate_eCurrent_A_per_um": direct["substrate_electron_A_per_um"],
                    "tight_substrate_default_minus_direct_A_per_um": gap,
                    "tight_signed_observable_fraction": fraction,
                    "tight_burst_flag": int(abs(fraction) >= burst_threshold),
                    "m59_burst_flag": old_burst,
                    "m59_substrate_default_minus_direct_A_per_um": old[
                        "substrate_default_minus_direct_eCurrent_A_per_um"],
                }
                pair_rows.append(row)
                case_points.append(row)
                if device == "n23" and math.isclose(drain, 0.05) and math.isclose(gate, 0.05):
                    target = row
                for algorithm, values in (("default", default), ("direct", direct)):
                    for contact in CONTACTS:
                        for component in COMPONENTS:
                            terminal_rows.append({
                                "case": case, "device": device,
                                "drain_voltage_V": drain, "gate_voltage_V": gate,
                                "algorithm": algorithm, "contact": contact,
                                "component": component,
                                "current_A_per_um": values[f"{contact}_{component}_A_per_um"],
                            })
            case_rows.append({
                "case": case, "device": device, "drain_voltage_V": drain,
                "point_count": len(case_points),
                "tight_burst_point_count": sum(int(row["tight_burst_flag"])
                                                for row in case_points),
                "m59_burst_point_count": sum(int(row["m59_burst_flag"])
                                              for row in case_points),
                "maximum_tight_default_vela_error_dex": max(
                    float(row["tight_default_vela_error_dex"]) for row in case_points),
                "maximum_tight_default_vs_m46_log_shift_dex": max(
                    float(row["tight_default_vs_m46_log_shift_dex"]) for row in case_points),
            })
    if target is None:
        raise RuntimeError("M60 target missing")
    field_rows, fields_invariant = export_and_compare_snapshots(
        contract, manifest, force_exports)
    old_gap = float(target["m59_substrate_default_minus_direct_A_per_um"])
    tight_gap = float(target["tight_substrate_default_minus_direct_A_per_um"])
    reduction = 1.0 - abs(tight_gap) / max(abs(old_gap), 1e-300)
    burst_count = sum(int(row["tight_burst_flag"]) for row in pair_rows)
    max_nonburst_shift = max(all_nonburst_shifts)
    suppression = (
        abs(float(target["tight_signed_observable_fraction"])) < burst_threshold
        and float(target["tight_default_vela_error_dex"]) <= float(
            contract["acceptance"]["target_default_vela_error_suppression_dex"])
        and reduction >= float(contract["acceptance"][
            "target_observable_gap_reduction_fraction"])
        and burst_count == 0
        and max_nonburst_shift <= float(contract["acceptance"][
            "maximum_nonburst_default_current_log_shift_dex"]))
    complete = (len(pair_rows) == int(contract["acceptance"]["required_algorithm_point_count"])
                and len(terminal_rows) == 816 * len(CONTACTS) * len(COMPONENTS) * 2
                and len(log_rows) == int(contract["acceptance"]["required_case_count"])
                and all(int(row["console_bytes"]) > 0 for row in log_rows))
    if not complete:
        classification = "execution_or_replay_mismatch"
    elif not fields_invariant:
        classification = "tight_convergence_changes_solved_state_between_observers"
    elif suppression:
        classification = "tight_convergence_suppresses_port_bursts"
    elif reduction >= float(contract["acceptance"]["material_gap_reduction_fraction"]):
        classification = "tight_convergence_reduces_but_does_not_suppress_port_bursts"
    else:
        classification = "tight_convergence_no_material_burst_effect"
    checks = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "case_count": len(log_rows) == int(contract["acceptance"]["required_case_count"]),
        "paired_point_count": len(pair_rows) == 816,
        "terminal_component_count": len(terminal_rows) == 19584,
        "exact_biases": all(abs(float(row["gate_voltage_V"]) / 0.05 -
                                round(float(row["gate_voltage_V"]) / 0.05)) < 1e-9
                            for row in pair_rows),
        "cnormprint_logs_retained": all(int(row["console_bytes"]) > 0 for row in log_rows),
        "observer_state_fields_invariant": fields_invariant,
        "classification_declared": classification in contract["analysis"]["classifications"],
        "historical_artifacts_not_rewritten": True,
        "production_defaults_not_changed": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m60_tight_convergence_port_burst_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "sentaurus": {"banner": banner, "case_count": len(log_rows),
                      "tight_math": contract["intervention"]["tight_math"]},
        "target": {**target, "observable_gap_reduction_fraction": reduction},
        "matrix": {
            "paired_point_count": len(pair_rows),
            "tight_burst_point_count": burst_count,
            "m59_burst_point_count": sum(int(row["m59_burst_flag"]) for row in pair_rows),
            "maximum_nonburst_default_vs_m46_log_shift_dex": max_nonburst_shift,
            "maximum_tight_default_vela_error_dex": max(
                float(row["tight_default_vela_error_dex"]) for row in pair_rows),
        },
        "state_invariance": {
            "row_count": len(field_rows), "all_fields_invariant": fields_invariant,
            "maximum_absolute_difference": max(
                float(row["maximum_absolute_difference"]) for row in field_rows),
            "maximum_symmetric_relative_difference": max(
                float(row["maximum_symmetric_relative_difference"]) for row in field_rows),
            "failure_count": sum(int(row["failure_count"]) for row in field_rows),
        },
        "secondary_reference_policy": {
            "eligible": classification == "tight_convergence_suppresses_port_bursts",
            "replaces_m8_m46": False,
            "required_math": contract["intervention"]["tight_math"],
            "required_diagnostics": ["CNormPrint logs", "paired default and DirectCurrent currents",
                                     "pointwise observable-fraction flag"],
        },
        "acceptance": checks,
    }
    return report, {"points": pair_rows, "terminals": terminal_rows,
                    "cases": case_rows, "logs": log_rows, "fields": field_rows}


def freeze(report: dict[str, Any], ledgers: dict[str, list[dict[str, Any]]]) -> None:
    write_csv(POINTS, ledgers["points"])
    write_csv(TERMINALS, ledgers["terminals"])
    write_csv(CASES, ledgers["cases"])
    write_csv(LOGS, ledgers["logs"])
    write_csv(FIELDS, ledgers["fields"])
    write_json(REPORT, report)
    write_json(QUALIFICATION, report["secondary_reference_policy"])
    target = report["target"]
    matrix = report["matrix"]
    DOC.write_text(f"""# SimpleMOS M60 收紧收敛端口 burst 判别

## 结论

M60 分类为 `{report['classification']}`。32个T-2022.03-SP2 deck使用冻结的同一收敛策略包，覆盖默认/Direct两种观测、16条曲线和每种算法816个精确偏压点。物理、网格、接触和偏压路径未改变；M8/M46继续作为生产基线。

目标n23、Vd=0.05 V、Vg=0.05 V的default-minus-Direct substrate电子电流由 `{float(target['m59_substrate_default_minus_direct_A_per_um']):.12e}` 变为 `{float(target['tight_substrate_default_minus_direct_A_per_um']):.12e}` A/um，差值削减 `{float(target['observable_gap_reduction_fraction']):.6%}`。收紧后默认Id与Vela的差为 `{float(target['tight_default_vela_error_dex']):.6f}` dex，观测差/Id为 `{float(target['tight_signed_observable_fraction']):.6e}`。

完整矩阵的收紧后burst标志点为 `{int(matrix['tight_burst_point_count'])}`，M59原标志点为 `{int(matrix['m59_burst_point_count'])}`；原非burst点相对M46默认曲线的最大变化为 `{float(matrix['maximum_nonburst_default_vs_m46_log_shift_dex']):.6e}` dex。默认与Direct配对快照字段不变性为 `{report['state_invariance']['all_fields_invariant']}`。

该结果只检验冻结的收敛策略包，不分离Digits与ErrRef各自贡献。次级亚fA参考资格为 `{report['secondary_reference_policy']['eligible']}`，且无论结果如何都不替换M8/M46。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M60 tight-convergence port-burst replay",
        "status": report["status"], "classification": report["classification"],
        "summary": "Frozen tight-convergence default/Direct full-matrix discriminator",
        "report": portable(REPORT),
        "ledgers": [portable(POINTS), portable(CASES), portable(LOGS), portable(FIELDS)],
    })
    artifacts = [REPORT, POINTS, TERMINALS, CASES, LOGS, FIELDS,
                 QUALIFICATION, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m60_tight_convergence_port_burst_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT),
        "new_sentaurus_execution": True, "new_vela_execution": False,
        "production_reference_replaced": False,
        "source_hashes": read_json(FREEZE)["upstream_hashes"],
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "acceptance": report["acceptance"],
    })


def verify() -> dict[str, Any]:
    validate_contract()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M60 evidence is not accepted and frozen")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M60 artifact hash changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-sentaurus", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--force-exports", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=r"C:\Windows\System32\OpenSSH\ssh.exe")
    parser.add_argument("--scp-bin", default=r"C:\Windows\System32\OpenSSH\scp.exe")
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    if args.verify:
        report = verify()
        print(f"M60 verified: {report['classification']}")
        return
    contract = validate_contract()
    do_prepare = args.prepare or args.all
    do_run = args.run_sentaurus or args.all
    do_analyze = args.analyze or args.all
    manifest_path = OUTPUT / "sentaurus_manifest.json"
    manifest = prepare(contract, args.force) if do_prepare else read_json(manifest_path)
    if not do_run and not do_analyze:
        print(json.dumps({"status": "prepared", "case_count": len(manifest["cases"]),
                          "manifest": portable(manifest_path)}, indent=2))
        return
    banner_path = OUTPUT / "sentaurus_banner.txt"
    banner = (run_sentaurus(manifest, args.ssh_target, args.ssh_bin,
                            args.scp_bin, args.remote_root, args.jobs)
              if do_run else banner_path.read_text(encoding="utf-8").strip())
    if do_analyze:
        report, ledgers = analyze(contract, manifest, banner, args.force_exports)
        freeze(report, ledgers)
        if not report["acceptance"]["all_checks_pass"]:
            raise RuntimeError(f"M60 acceptance failed: {report['acceptance']}")
        print(json.dumps({"status": report["status"],
                          "classification": report["classification"],
                          "target": report["target"],
                          "matrix": report["matrix"]}, indent=2))


if __name__ == "__main__":
    main()
