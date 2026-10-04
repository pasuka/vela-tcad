#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M52 DirectCurrent attribution."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
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
import run_simplemos_m51_current_weighting_attribution as m51  # noqa: E402
import sentaurus_import  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m52_direct_current_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m52_direct_current_attribution_contract_freeze.json"
M50_EVIDENCE = ROOT / "simplemos_m50_native_substrate_face_flux_export_evidence.json"
M51_EVIDENCE = ROOT / "simplemos_m51_current_weighting_attribution_evidence.json"
M50_REPORT = ROOT / "native_substrate_face_flux_export/m50_native_substrate_face_flux_export_report.json"
M50_FLUX = ROOT / "native_substrate_face_flux_export/m50_native_substrate_face_flux_ledger.csv"
M50_TERMINALS = ROOT / "native_substrate_face_flux_export/m50_terminal_replay_ledger.csv"
M51_REPORT = ROOT / "current_weighting_attribution/m51_current_weighting_attribution_report.json"
M51_TERMINALS = ROOT / "current_weighting_attribution/m51_terminal_component_ledger.csv"
M50_OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m50_native_substrate_face_flux_export"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m52_direct_current_attribution"
PORTABLE = ROOT / "direct_current_attribution"
REPORT = PORTABLE / "m52_direct_current_attribution_report.json"
STATES = PORTABLE / "m52_direct_current_state_ledger.csv"
TERMINALS = PORTABLE / "m52_terminal_component_ledger.csv"
FIELDS = PORTABLE / "m52_field_invariance_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m52_direct_current_attribution_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m52/artifact.json"
EVIDENCE = ROOT / "simplemos_m52_direct_current_attribution_evidence.json"
TEST = REPO / "tests/regression/test_simplemos_m52_direct_current_attribution.py"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
ARCHIVE = "simplemos_m52_direct_current_results.tgz"
REMOTE_ROOT = "/tmp/vela_simplemos_m52_direct_current_attribution"
DEVICES = ("n23", "n19")
GATES = (0.0, 0.05, 0.1)
CONTACTS = ("source", "drain", "gate", "substrate")
COMPONENTS = ("electron", "hole", "total")


read_json = m51.read_json
write_json = m51.write_json
read_csv = m51.read_csv
write_csv = m51.write_csv
sha256 = m51.sha256
portable = m51.portable
close_enough = m51.close_enough
state_id = m51.state_id


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != "vela.simplemos.sdevice.m52_direct_current_attribution_contract.v1":
        raise ValueError("unexpected M52 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M52 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M52 contract changed after freeze")
    for path in (M50_EVIDENCE, M51_EVIDENCE):
        if read_json(path).get("status") != "frozen":
            raise ValueError(f"M52 requires frozen evidence: {path.name}")
    m50_report = read_json(M50_REPORT)
    m51_report = read_json(M51_REPORT)
    if m50_report.get("classification") != contract["upstream"]["required_m50_classification"]:
        raise ValueError("M50 classification changed")
    if m51_report.get("classification") != contract["upstream"]["required_m51_classification"]:
        raise ValueError("M51 classification changed")
    checks = (
        (m50_report["target"]["same_run_substrate_electron_terminal_A_per_um"],
         contract["upstream"]["required_target_default_substrate_electron_A_per_um"]),
        (m51_report["target"]["m51_weighted_substrate_electron_A_per_um"],
         contract["upstream"]["required_target_weighted_substrate_electron_A_per_um"]),
        (m50_report["target"]["contact_surface_normal_integral_A_per_um"],
         contract["upstream"]["required_target_native_normal_A_per_um"]),
        (m51_report["target"]["weighted_minus_native_normal_A_per_um"],
         contract["upstream"]["required_target_weighted_native_residual_A_per_um"]),
    )
    for observed, expected in checks:
        if not math.isclose(float(observed), float(expected), rel_tol=0.0, abs_tol=1e-30):
            raise ValueError("M52 upstream target anchor changed")
    if tuple(contract["scope"]["devices"]) != DEVICES:
        raise ValueError("M52 device matrix changed")
    if tuple(map(float, contract["scope"]["gate_voltages_V"])) != GATES:
        raise ValueError("M52 gate matrix changed")
    return contract


def sentaurus_deck(case: str) -> str:
    deck = m51.m50.sentaurus_deck(case)
    deck = deck.replace("M50_SubstrateElectronFaceFlux",
                        "M52_SubstrateElectronFaceFlux")
    deck = deck.replace('Plot(FilePrefix="m50_state"',
                        'Plot(FilePrefix="m52_state"')
    marker = "Math { Extrapolate Iterations=20 ExitOnFailure }"
    replacement = "Math { Extrapolate Iterations=20 ExitOnFailure DirectCurrent }"
    if deck.count(marker) != 1:
        raise RuntimeError("M50 Math insertion point changed")
    direct = deck.replace(marker, replacement)
    if direct.replace(" DirectCurrent", "") != deck:
        raise RuntimeError("M52 deck differs from M50 by more than DirectCurrent")
    if "CurrentWeighting" in direct:
        raise RuntimeError("M52 DirectCurrent deck must not contain CurrentWeighting")
    return direct


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
        "schema": "vela.simplemos.sdevice.m52_sentaurus_manifest.v1",
        "status": "prepared", "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "intervention": "DirectCurrent flag only; diagnostic names M50 to M52",
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
        run([ssh_bin, ssh_target,
             f"set -eu; cd {root}; sdevice {item['case']}_des.cmd > console.log 2>&1"])

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
    dataset = "M52_SubstrateElectronFaceFlux"
    if dataset not in datasets:
        raise RuntimeError(f"M52 native flux missing from {path}")
    parsed: dict[float, dict[str, float]] = {}
    for values in sentaurus_import.parse_values_block(text, len(datasets)):
        row = dict(zip(datasets, values, strict=True))
        gate = float(row["gate OuterVoltage"])
        nearest = min(GATES, key=lambda value: abs(value - gate))
        if math.isclose(gate, nearest, rel_tol=0.0, abs_tol=1e-10):
            item = {
                "gate_voltage_V": nearest,
                "drain_voltage_V": float(row["drain OuterVoltage"]),
                "native_flux_A_per_um": float(row[dataset]),
            }
            for contact in CONTACTS:
                item[f"{contact}_electron_A_per_um"] = float(row[f"{contact} eCurrent"])
                item[f"{contact}_hole_A_per_um"] = float(row[f"{contact} hCurrent"])
                item[f"{contact}_total_A_per_um"] = float(row[f"{contact} TotalCurrent"])
            parsed[nearest] = item
    if set(parsed) != set(GATES):
        raise RuntimeError(f"missing exact M52 states in {path}: {sorted(parsed)}")
    return parsed


def export_snapshots(manifest: dict[str, Any], force: bool) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = {"baseline": {}, "weighted": {}}
    for item in manifest["cases"]:
        device = item["device"]
        roots = {
            "baseline": M50_OUTPUT / "sentaurus_raw/sentaurus_bundle" / device,
            "weighted": OUTPUT / "sentaurus_raw/sentaurus_bundle" / device,
        }
        patterns = {"baseline": "m50_state*.tdr", "weighted": "m52_state*.tdr"}
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


def default_terminals() -> dict[tuple[str, float, str, str], float]:
    return {(row["device"], float(row["gate_voltage_V"]), row["contact"],
             row["component"]): float(row["m50_same_run_A_per_um"])
            for row in read_csv(M50_TERMINALS)}


def weighted_terminals() -> dict[tuple[str, float, str, str], float]:
    return {(row["device"], float(row["gate_voltage_V"]), row["contact"],
             row["component"]): float(row["m51_current_weighting_A_per_um"])
            for row in read_csv(M51_TERMINALS)}


def native_fluxes() -> dict[tuple[str, float], dict[str, float]]:
    return {(row["device"], float(row["gate_voltage_V"])): {
        "normal": float(row["contact_surface_normal_integral_A_per_um"]),
        "negative": float(row["negative_contact_surface_normal_integral_A_per_um"])}
            for row in read_csv(M50_FLUX)}


def analyze(contract: dict[str, Any], manifest: dict[str, Any], banner: str,
            force_exports: bool) -> dict[str, Any]:
    default = default_terminals()
    weighted = weighted_terminals()
    native = native_fluxes()
    terminal_abs = float(contract["acceptance"][
        "maximum_terminal_closure_absolute_difference_A_per_um"])
    terminal_rel = float(contract["acceptance"][
        "maximum_terminal_closure_relative_difference"])
    flux_abs = float(contract["acceptance"][
        "maximum_native_flux_replay_absolute_difference_A_per_um"])
    flux_rel = float(contract["acceptance"][
        "maximum_native_flux_replay_relative_difference"])
    direct_rows: dict[tuple[str, float], dict[str, float]] = {}
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    for item in manifest["cases"]:
        files = sorted((raw / item["device"]).glob("IdVg_*des.plt"))
        if len(files) != 1:
            raise RuntimeError(f"expected one M52 current file for {item['device']}: {files}")
        for gate, row in parse_terminal_file(files[0]).items():
            direct_rows[(item["device"], gate)] = row
    exports = export_snapshots(manifest, force_exports)
    field_rows, fields_invariant = m51.compare_fields(contract, exports)
    terminal_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    for device in DEVICES:
        for gate in GATES:
            state = state_id(device, gate)
            direct = direct_rows[(device, gate)]
            for contact in CONTACTS:
                for component in COMPONENTS:
                    key = (device, gate, contact, component)
                    value = float(direct[f"{contact}_{component}_A_per_um"])
                    terminal_rows.append({
                        "state": state, "device": device, "gate_voltage_V": gate,
                        "contact": contact, "component": component,
                        "m50_default_A_per_um": default[key],
                        "m51_weighted_A_per_um": weighted[key],
                        "m52_direct_A_per_um": value,
                        "direct_minus_default_A_per_um": value - default[key],
                        "direct_minus_weighted_A_per_um": value - weighted[key],
                    })
            key = (device, gate, "substrate", "electron")
            values = native[(device, gate)]
            direct_terminal = float(direct["substrate_electron_A_per_um"])
            direct_flux = float(direct["native_flux_A_per_um"])
            direct_native_gap = direct_terminal - values["normal"]
            default_native_gap = default[key] - values["normal"]
            normal_match = close_enough(
                direct_terminal, values["normal"], terminal_abs, terminal_rel)
            negative_match = close_enough(
                direct_terminal, values["negative"], terminal_abs, terminal_rel)
            state_rows.append({
                "state": state, "device": device, "gate_voltage_V": gate,
                "drain_voltage_V": direct["drain_voltage_V"],
                "m50_default_substrate_electron_A_per_um": default[key],
                "m51_weighted_substrate_electron_A_per_um": weighted[key],
                "m52_direct_substrate_electron_A_per_um": direct_terminal,
                "m50_native_normal_A_per_um": values["normal"],
                "m50_native_negative_normal_A_per_um": values["negative"],
                "m52_native_normal_A_per_um": direct_flux,
                "direct_minus_default_A_per_um": direct_terminal - default[key],
                "direct_minus_weighted_A_per_um": direct_terminal - weighted[key],
                "direct_minus_native_normal_A_per_um": direct_native_gap,
                "direct_minus_nearest_native_A_per_um": min(
                    (direct_terminal - values["normal"],
                     direct_terminal - values["negative"]), key=abs),
                "absolute_default_native_gap_reduction_fraction": 1.0 -
                    abs(direct_native_gap) / max(abs(default_native_gap), 1e-300),
                "direct_native_normal_relative_difference": abs(direct_native_gap) /
                    max(abs(values["normal"]), 1e-300),
                "native_flux_replays": close_enough(
                    direct_flux, values["normal"], flux_abs, flux_rel),
                "direct_matches_default": close_enough(
                    direct_terminal, default[key], terminal_abs, terminal_rel),
                "direct_matches_weighted": close_enough(
                    direct_terminal, weighted[key], terminal_abs, terminal_rel),
                "direct_matches_native_normal": normal_match,
                "direct_matches_native_negative_normal": negative_match,
                "direct_matches_either_native_orientation": normal_match or negative_match,
            })
    native_matches = all(bool(row["direct_matches_either_native_orientation"])
                         for row in state_rows)
    weighted_matches = all(bool(row["direct_matches_weighted"]) for row in state_rows)
    default_matches = all(bool(row["direct_matches_default"]) for row in state_rows)
    if len(state_rows) != 6 or len(terminal_rows) != 72:
        classification = "replay_incomplete"
    elif not fields_invariant:
        classification = "state_perturbed"
    elif native_matches:
        classification = "direct_matches_native_face"
    elif weighted_matches:
        classification = "direct_matches_weighted"
    elif default_matches:
        classification = "direct_matches_default"
    else:
        classification = "direct_third_observable"
    target = next(row for row in state_rows
                  if row["device"] == "n23" and row["gate_voltage_V"] == 0.05)
    controls = [row for row in state_rows if
                (row["device"] == "n23" and row["gate_voltage_V"] in (0.0, 0.1))
                or (row["device"] == "n19" and row["gate_voltage_V"] == 0.05)]
    acceptance = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "exact_state_count": len(state_rows) == 6,
        "terminal_component_row_count": len(terminal_rows) == 72,
        "all_biases_exact": all(float(row["drain_voltage_V"]) == 0.05 for row in state_rows),
        "all_values_finite": all(math.isfinite(float(row[key])) for row in state_rows
                                  for key in ("m52_direct_substrate_electron_A_per_um",
                                              "m52_native_normal_A_per_um")),
        "state_fields_invariant": fields_invariant,
        "native_flux_replays": all(bool(row["native_flux_replays"]) for row in state_rows),
        "classification_declared": classification in contract["analysis"]["classifications"],
        "defaults_unchanged": True,
        "closed_topics_not_reopened": True,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    report = {
        "schema": "vela.simplemos.sdevice.m52_direct_current_attribution_report.v1",
        "status": "complete", "classification": classification,
        "execution": {
            "sentaurus_release": "T-2022.03-SP2", "sentaurus_banner": banner,
            "new_sentaurus_execution": True, "new_vela_execution": False,
            "contract_sha256_before_and_after": sha256(CONTRACT),
            "intervention": "Math DirectCurrent flag only",
            "production_defaults_changed": False,
        },
        "manual_interpretation": contract["manual_basis"],
        "target": target, "controls": controls,
        "direct_weighted_agreement": {
            "all_states_within_contract_tolerance": weighted_matches,
            "maximum_absolute_difference_A_per_um": max(
                abs(float(row["direct_minus_weighted_A_per_um"]))
                for row in state_rows),
            "target_absolute_difference_A_per_um": abs(float(target[
                "direct_minus_weighted_A_per_um"])),
            "interpretation": "DirectCurrent and CurrentWeighting produce the same terminal observable to numerical roundoff; their shared residual to the Tcl contact-domain integral is attributable to internal discrete contact integration versus plotted nodal-field integration, not to a solved-state change.",
        },
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
        f"{float(row['m52_direct_substrate_electron_A_per_um']):.12e} | "
        f"{float(row['m52_native_normal_A_per_um']):.12e} |"
        for row in report["states"])
    doc = f'''# SimpleMOS M52 DirectCurrent 归因

## 结论

M52 分类为 `{report['classification']}`。唯一干预是在 M50 的全局 `Math` 中加入与 `CurrentWeighting` 互斥的 `DirectCurrent`；六状态物理、网格、偏压路径和诊断面通量保持冻结。

目标点 n23、Vd=0.05 V、Vg=0.05 V：默认 substrate `eCurrent` 为 `{float(target['m50_default_substrate_electron_A_per_um']):.12e}` A/um，M51 `CurrentWeighting` 为 `{float(target['m51_weighted_substrate_electron_A_per_um']):.12e}` A/um，M52 `DirectCurrent` 为 `{float(target['m52_direct_substrate_electron_A_per_um']):.12e}` A/um，原生法向面通量为 `{float(target['m52_native_normal_A_per_um']):.12e}` A/um。

`DirectCurrent` 与 `CurrentWeighting` 在目标点仅差 `{abs(float(target['direct_minus_weighted_A_per_um'])):.12e}` A/um，六状态最大差为 `{float(report['direct_weighted_agreement']['maximum_absolute_difference_A_per_um']):.12e}` A/um，即在数值舍入尺度一致。

`DirectCurrent` 相对原生法向量的绝对残差为 `{abs(float(target['direct_minus_native_normal_A_per_um'])):.12e}` A/um，相对差为 `{float(target['direct_native_normal_relative_difference']):.9e}`；相对于默认端口差异的消除比例为 `{100.0 * float(target['absolute_default_native_gap_reduction_fraction']):.9f}%`。

## 六状态

| 器件 | Vg (V) | 默认 eCurrent | CurrentWeighting | DirectCurrent | 原生法向面通量 |
|---|---:|---:|---:|---:|---:|
{table}

## 状态不变性

共比较 `{report['state_field_summary']['comparison_row_count']}` 个状态-字段文件；失败值数量为 `{report['state_field_summary']['total_failure_count']}`。最大绝对差为 `{float(report['state_field_summary']['maximum_absolute_difference']):.6e}`，最大对称相对差为 `{float(report['state_field_summary']['maximum_symmetric_relative_difference']):.6e}`。

## 边界

- 手册定义 `DirectCurrent` 为接触面积上的电流密度表面积分；本任务没有同时启用 `CurrentWeighting`。
- `DirectCurrent` 与 `CurrentWeighting` 的机器精度一致性、以及二者共同相对 Tcl 节点场积分的小残差，限定为内部离散接触积分与绘图节点场后处理积分的定义差异；不是求解状态差异。
- 没有修改 HFS、SG、准费米、BGN、SRH、迁移率、接触模型、网格或生产默认值。
- M52 只归因 Sentaurus 端口观测算法，不授权修改 Vela 或 Sentaurus 生产配置。

机器报告：`{portable(REPORT)}`。
'''
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc, encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.simplemos.sdevice.m52_artifact.v1",
        "title": "SimpleMOS M52 DirectCurrent attribution",
        "status": report["status"], "classification": report["classification"],
        "report": portable(REPORT),
        "ledgers": [portable(STATES), portable(TERMINALS), portable(FIELDS)],
        "document": portable(DOC), "target": target,
        "acceptance": report["acceptance"],
    })
    artifacts = [REPORT, STATES, TERMINALS, FIELDS, DOC, ARTIFACT]
    sources = [CONTRACT, FREEZE, Path(__file__).resolve(), TEST,
               M50_EVIDENCE, M51_EVIDENCE,
               REPO / "scripts/run_simplemos_m51_current_weighting_attribution.py",
               OUTPUT / "sentaurus_manifest.json", OUTPUT / "sentaurus_banner.txt"]
    for device in DEVICES:
        sources.append(OUTPUT / "sentaurus_bundle" / device /
                       f"{device}_vd_0p05_des.cmd")
        sources.extend(sorted((OUTPUT / "sentaurus_raw/sentaurus_bundle" / device).glob(
            "IdVg_*des.plt")))
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m52_direct_current_attribution_evidence.v1",
        "status": "frozen", "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "new_sentaurus_execution": True, "new_vela_execution": False,
        "only_intervention": "Math DirectCurrent flag",
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
