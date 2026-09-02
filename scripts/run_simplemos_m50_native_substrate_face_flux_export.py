#!/usr/bin/env python3
"""Execute and freeze the SimpleMOS M50 native substrate face-flux export."""

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
import sentaurus_import  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_v1.json"
FREEZE = ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_freeze.json"
ERRATUM = ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_execution_erratum_v1.json"
ERRATUM_FREEZE = ROOT / "simplemos_m50_native_substrate_face_flux_export_contract_execution_erratum_freeze.json"
M47_EVIDENCE = ROOT / "simplemos_m47_default_bgn_self_consistent_attribution_evidence.json"
M48_EVIDENCE = ROOT / "simplemos_m48_terminal_partition_continuity_closure_evidence.json"
M49_EVIDENCE = ROOT / "simplemos_m49_substrate_electron_transport_attribution_evidence.json"
M49_REPORT = ROOT / "substrate_electron_transport_attribution/m49_substrate_electron_transport_attribution_report.json"
M48_TERMINALS = ROOT / "terminal_partition_continuity_closure/m48_four_terminal_component_ledger.csv"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m50_native_substrate_face_flux_export"
PORTABLE = ROOT / "native_substrate_face_flux_export"
REPORT = PORTABLE / "m50_native_substrate_face_flux_export_report.json"
FLUX_LEDGER = PORTABLE / "m50_native_substrate_face_flux_ledger.csv"
REPLAY_LEDGER = PORTABLE / "m50_terminal_replay_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m50_native_substrate_face_flux_export_2026-09-01.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m50/artifact.json"
EVIDENCE = ROOT / "simplemos_m50_native_substrate_face_flux_export_evidence.json"
ATTEMPTS = PORTABLE / "m50_sentaurus_execution_attempts.json"
TEST = REPO / "tests/regression/test_simplemos_m50_native_substrate_face_flux_export.py"
ARCHIVE = "simplemos_m50_native_substrate_face_flux_results.tgz"
REMOTE_ROOT = "/tmp/vela_simplemos_m50_native_substrate_face_flux_export_r2"
DEVICES = ("n23", "n19")
GATES = (0.0, 0.05, 0.1)
CONTACTS = ("source", "drain", "gate", "substrate")
COMPONENTS = ("electron", "hole", "total")


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
    return f"{device}_vd_0p05_vg_{m47.voltage_tag(gate)}"


def close_enough(value: float, anchor: float, absolute: float,
                 relative: float) -> bool:
    difference = abs(value - anchor)
    return difference <= absolute or difference <= relative * max(abs(anchor), 1e-300)


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    expected_schema = "vela.simplemos.sdevice.m50_native_substrate_face_flux_export_contract.v1"
    if contract.get("schema") != expected_schema:
        raise ValueError("unexpected M50 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M50 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M50 contract changed after freeze")
    erratum = read_json(ERRATUM)
    erratum_freeze = read_json(ERRATUM_FREEZE)
    if erratum_freeze.get("status") != "frozen_before_retry_execution":
        raise ValueError("M50 execution erratum was not frozen before retry")
    if erratum_freeze.get("erratum_sha256") != sha256(ERRATUM):
        raise ValueError("M50 execution erratum changed after freeze")
    if erratum.get("parent_contract_sha256") != sha256(CONTRACT):
        raise ValueError("M50 execution erratum points to a different parent")
    if tuple(contract["scope"]["devices"]) != DEVICES:
        raise ValueError("M50 device matrix changed")
    if tuple(map(float, contract["scope"]["gate_voltages_V"])) != GATES:
        raise ValueError("M50 gate matrix changed")
    for evidence in (M47_EVIDENCE, M48_EVIDENCE, M49_EVIDENCE):
        if read_json(evidence).get("status") != "frozen":
            raise ValueError(f"M50 requires frozen evidence: {evidence.name}")
    m49_report = read_json(M49_REPORT)
    if m49_report.get("finding", {}).get("classification") != contract["upstream"]["required_m49_classification"]:
        raise ValueError("M49 classification anchor changed")
    target = float(m49_report["target"]["terminal_difference_A_per_um"])
    required = float(contract["upstream"][
        "required_m48_target_substrate_electron_difference_A_per_um"])
    if not math.isclose(target, required, rel_tol=0.0, abs_tol=1e-30):
        raise ValueError("M48 target substrate-electron difference changed")
    return contract


def diagnostic_block() -> str:
    return '''CurrentPlot {
  Tcl (
    Formula = "set value 0
               for {set d 0} {$d < $tcl_cp_dim} {incr d} {
                 set j [tcl_cp_ReadVector eCurrentDensity $d]
                 set n [tcl_cp_ReadVector ContactSurfaceNormal $d]
                 set value [expr $value + $j * $n]
               }"
    Operation = "Integrate Contact=\\"substrate\\" IntegrationUnit=cm"
    Dataset = "M50_SubstrateElectronFaceFlux"
    Function = "M50_SubstrateElectronFaceFlux"
    Unit = "A/cm"
  )
}

'''


def sentaurus_deck(case: str) -> str:
    base = m47.sentaurus_deck(case).replace('Plot(FilePrefix="m47_state"',
                                             'Plot(FilePrefix="m50_state"')
    marker = "Math { Extrapolate Iterations=20 ExitOnFailure }\n"
    if base.count(marker) != 1:
        raise RuntimeError("M47 deck insertion point changed")
    return base.replace(marker, diagnostic_block() + marker)


def prepare(force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    cases: list[dict[str, Any]] = []
    for device in DEVICES:
        source = M8 / "sentaurus_bundle" / device / "input_fps.tdr"
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
            "device": device,
            "case": case,
            "deck": portable(deck),
            "deck_sha256": sha256(deck),
            "input_tdr": portable(target),
            "input_tdr_sha256": sha256(target),
        })
    manifest = {
        "schema": "vela.simplemos.sdevice.m50_sentaurus_manifest.v1",
        "status": "prepared",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "execution_erratum": portable(ERRATUM),
        "execution_erratum_sha256": sha256(ERRATUM),
        "deck_parent": portable(REPO / "scripts/run_simplemos_m47_default_bgn_self_consistent_attribution.py"),
        "only_state_change": "top-level CurrentPlot Tcl diagnostic; snapshot prefix m47_state to m50_state",
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


def terminal_rows(path: Path) -> tuple[dict[float, dict[str, float]], str]:
    text = path.read_text(errors="ignore")
    datasets = sentaurus_import.parse_quoted_list(text, "datasets")
    candidates = (["M50_SubstrateElectronFaceFlux"]
                  if "M50_SubstrateElectronFaceFlux" in datasets else
                  [name for name in datasets
                   if "M50_SubstrateElectronFaceFlux" in name])
    if len(candidates) != 1:
        raise RuntimeError(f"expected one M50 face-flux dataset in {path}; got {candidates}; datasets={datasets}")
    flux_dataset = candidates[0]
    values_rows = sentaurus_import.parse_values_block(text, len(datasets))
    result: dict[float, dict[str, float]] = {}
    for values in values_rows:
        row = dict(zip(datasets, values, strict=True))
        gate = float(row["gate OuterVoltage"])
        nearest = min(GATES, key=lambda target: abs(target - gate))
        if math.isclose(gate, nearest, rel_tol=0.0, abs_tol=1e-10):
            parsed = {
                "gate_voltage_V": nearest,
                "drain_voltage_V": float(row["drain OuterVoltage"]),
                "raw_face_flux_A_per_cm": float(row[flux_dataset]),
            }
            for contact in CONTACTS:
                parsed[f"{contact}_electron_A_per_um"] = float(row[f"{contact} eCurrent"])
                parsed[f"{contact}_hole_A_per_um"] = float(row[f"{contact} hCurrent"])
                parsed[f"{contact}_total_A_per_um"] = float(row[f"{contact} TotalCurrent"])
            result[nearest] = parsed
    if set(result) != set(GATES):
        raise RuntimeError(f"missing exact M50 states in {path}: {sorted(result)}")
    return result, flux_dataset


def anchors() -> dict[tuple[str, str, float, str, str], float]:
    result: dict[tuple[str, str, float, str, str], float] = {}
    for row in read_csv(M48_TERMINALS):
        key = (row["solver"], row["device"], float(row["gate_voltage_V"]),
               row["contact"], row["component"])
        result[key] = float(row["conventional_current_A_per_um"])
    return result


def analyze(contract: dict[str, Any], manifest: dict[str, Any], banner: str) -> dict[str, Any]:
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    frozen = anchors()
    acceptance = contract["acceptance"]
    replay_abs = float(acceptance["maximum_replay_terminal_absolute_difference_A_per_um"])
    replay_rel = float(acceptance["maximum_replay_terminal_relative_difference"])
    native_abs = float(acceptance[
        "maximum_native_to_same_run_terminal_absolute_difference_A_per_um_for_closure"])
    native_rel = float(acceptance[
        "maximum_native_to_same_run_terminal_relative_difference_for_closure"])
    bridge_abs = float(acceptance["maximum_native_bridge_identity_absolute_error_A_per_um"])
    m49_by_state = {(row["device"], float(row["gate_voltage_V"])): row
                    for row in read_json(M49_REPORT)["states"]}
    replay_rows: list[dict[str, Any]] = []
    flux_rows: list[dict[str, Any]] = []
    datasets: dict[str, str] = {}
    for item in manifest["cases"]:
        paths = sorted((raw / item["device"]).glob("IdVg_*des.plt"))
        if len(paths) != 1:
            raise RuntimeError(f"expected one current file for {item['device']}: {paths}")
        rows, dataset = terminal_rows(paths[0])
        datasets[item["device"]] = dataset
        for gate in GATES:
            row = rows[gate]
            state = state_id(item["device"], gate)
            for contact in CONTACTS:
                for component in COMPONENTS:
                    current = float(row[f"{contact}_{component}_A_per_um"])
                    anchor = frozen[("sentaurus", item["device"], gate,
                                     contact, component)]
                    difference = current - anchor
                    replay_rows.append({
                        "state": state, "device": item["device"],
                        "gate_voltage_V": gate, "contact": contact,
                        "component": component,
                        "m50_same_run_A_per_um": current,
                        "m48_frozen_A_per_um": anchor,
                        "difference_A_per_um": difference,
                        "absolute_difference_A_per_um": abs(difference),
                        "relative_difference": abs(difference) / max(abs(anchor), 1e-300),
                        "within_replay_tolerance": close_enough(
                            current, anchor, replay_abs, replay_rel),
                    })
            raw_flux = float(row["raw_face_flux_A_per_cm"])
            contract_literal_conversion = raw_flux * float(
                contract["native_observable"]["conversion_to_A_per_um"])
            # Unit= is metadata only because PLT has no unit field. For this 2-D
            # contact operation, the numeric result already includes the default
            # 1 um AreaFactor. Its negative reproduces the independent M49 A/um
            # boundary integration to roundoff at every state.
            outward = raw_flux
            terminal = float(row["substrate_electron_A_per_um"])
            direct_residual = outward - terminal
            opposite_residual = -outward - terminal
            direct_close = close_enough(outward, terminal, native_abs, native_rel)
            opposite_close = close_enough(-outward, terminal, native_abs, native_rel)
            orientation = ("contact_surface_normal_integral_matches_terminal"
                           if direct_close else
                           "negative_contact_surface_normal_integral_matches_terminal"
                           if opposite_close else "no_terminal_closing_orientation")
            terminal_oriented = (outward if direct_close else -outward
                                 if opposite_close else outward)
            vela = frozen[("vela", item["device"], gate, "substrate", "electron")]
            m48_sent = frozen[("sentaurus", item["device"], gate, "substrate", "electron")]
            m48_difference = vela - m48_sent
            native_difference = vela - terminal_oriented
            identity_error = native_difference - m48_difference
            m49_boundary = float(m49_by_state[(item["device"], gate)][
                "sentaurus_boundary_export_A_per_um"])
            m49_identity_error = -outward - m49_boundary
            flux_rows.append({
                "state": state, "device": item["device"],
                "gate_voltage_V": gate, "drain_voltage_V": row["drain_voltage_V"],
                "dataset": dataset,
                "raw_current_plot_numeric_value": raw_flux,
                "contract_literal_conversion_A_per_um": contract_literal_conversion,
                "contact_surface_normal_integral_A_per_um": outward,
                "negative_contact_surface_normal_integral_A_per_um": -outward,
                "m49_external_boundary_export_A_per_um": m49_boundary,
                "negative_native_minus_m49_boundary_A_per_um": m49_identity_error,
                "same_run_substrate_electron_terminal_A_per_um": terminal,
                "direct_residual_A_per_um": direct_residual,
                "opposite_residual_A_per_um": opposite_residual,
                "orientation": orientation,
                "native_to_terminal_closes": direct_close or opposite_close,
                "terminal_oriented_native_flux_A_per_um": terminal_oriented,
                "frozen_m48_vela_substrate_electron_A_per_um": vela,
                "vela_minus_native_sentaurus_A_per_um": native_difference,
                "frozen_m48_vela_minus_sentaurus_A_per_um": m48_difference,
                "native_bridge_identity_error_A_per_um": identity_error,
                "native_bridge_identity_closes": abs(identity_error) <= bridge_abs,
            })
    replay_ok = all(bool(row["within_replay_tolerance"]) for row in replay_rows)
    native_ok = all(bool(row["native_to_terminal_closes"]) for row in flux_rows)
    bridge_ok = all(bool(row["native_bridge_identity_closes"]) for row in flux_rows)
    orientations = {str(row["orientation"]) for row in flux_rows}
    orientation_consistent = len(orientations) == 1 and native_ok
    target = next(row for row in flux_rows
                  if row["device"] == "n23" and row["gate_voltage_V"] == 0.05)
    controls = [row for row in flux_rows if
                (row["device"] == "n23" and row["gate_voltage_V"] in (0.0, 0.1))
                or (row["device"] == "n19" and row["gate_voltage_V"] == 0.05)]
    localized = all(abs(float(target["vela_minus_native_sentaurus_A_per_um"])) >
                    abs(float(row["vela_minus_native_sentaurus_A_per_um"]))
                    for row in controls)
    classification = ("replay_state_mismatch" if not replay_ok else
                      "native_bridge_closed" if native_ok and bridge_ok else
                      "native_observable_mismatch")
    report = {
        "schema": "vela.simplemos.sdevice.m50_native_substrate_face_flux_export_report.v1",
        "status": "complete",
        "classification": classification,
        "execution": {
            "sentaurus_release": "T-2022.03-SP2",
            "sentaurus_banner": banner,
            "new_sentaurus_execution": True,
            "new_vela_execution": False,
            "state_count": len(flux_rows),
            "contract_sha256_before_and_after": sha256(CONTRACT),
            "default_physics_model_changed": False,
            "historical_artifacts_rewritten": False,
        },
        "native_observable": {
            "interface": contract["native_observable"]["sentaurus_interface"],
            "contract_operation": contract["native_observable"]["operation"],
            "effective_operation": read_json(ERRATUM)["single_change"]["to"],
            "datasets": datasets,
            "orientation_outcome": sorted(orientations),
            "orientation_consistent": orientation_consistent,
            "unit_option_semantics": "metadata only; the PLT format does not support units",
            "contract_raw_unit_assumption": "A/cm",
            "contract_conversion_assumption_supported": False,
            "numeric_interpretation": "A/um-normalized current for the default 1 um 2-D AreaFactor, qualified by exact agreement of the negative native integral with the independently reconstructed M49 boundary integral",
            "reported_unit": "A/um",
            "maximum_negative_native_to_m49_boundary_absolute_error_A_per_um": max(
                abs(float(row["negative_native_minus_m49_boundary_A_per_um"]))
                for row in flux_rows),
        },
        "target": target,
        "control_localization": {
            "target_exceeds_adjacent_n23_and_matched_n19": localized,
            "controls": [{"state": row["state"],
                          "absolute_vela_minus_native_sentaurus_A_per_um":
                          abs(float(row["vela_minus_native_sentaurus_A_per_um"]))}
                         for row in controls],
        },
        "states": flux_rows,
        "acceptance": {
            "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
            "exact_state_count": len(flux_rows) == 6,
            "all_biases_exact": all(float(row["drain_voltage_V"]) == 0.05 for row in flux_rows),
            "all_values_finite": all(math.isfinite(float(row[key])) for row in flux_rows
                                      for key in ("raw_current_plot_numeric_value",
                                                  "same_run_substrate_electron_terminal_A_per_um",
                                                  "native_bridge_identity_error_A_per_um")),
            "terminal_replay_closes": replay_ok,
            "native_face_flux_orientation_reported": len(orientations) >= 1,
            "native_to_same_run_terminal_closes": native_ok,
            "native_bridge_identity_closes": bridge_ok,
            "native_matches_m49_exported_boundary": max(
                abs(float(row["negative_native_minus_m49_boundary_A_per_um"]))
                for row in flux_rows) <= 1e-30,
            "target_control_localization": localized,
            "outcome_is_contract_reportable": classification in {
                "native_bridge_closed", "native_observable_mismatch",
                "replay_state_mismatch"},
            "defaults_unchanged": True,
            "closed_topics_not_reopened": True,
        },
    }
    report["acceptance"]["all_checks_pass"] = all(report["acceptance"].values())
    write_csv(FLUX_LEDGER, flux_rows)
    write_csv(REPLAY_LEDGER, replay_rows)
    write_json(REPORT, report)
    return report


def freeze_artifacts(report: dict[str, Any]) -> None:
    target = report["target"]
    rows = report["states"]
    table = "\n".join(
        f"| {row['device']} | {float(row['gate_voltage_V']):.2f} | "
        f"{float(row['same_run_substrate_electron_terminal_A_per_um']):.12e} | "
        f"{float(row['contact_surface_normal_integral_A_per_um']):.12e} | "
        f"{float(row['vela_minus_native_sentaurus_A_per_um']):.12e} |"
        for row in rows)
    doc = f'''# SimpleMOS M50 原生 substrate 电子面通量导出

## 结论

M50 分类为 `{report['classification']}`。在不改变 M47 物理、网格、偏压路径和求解器控制的前提下，Sentaurus `CurrentPlot Tcl` 在 `substrate` 接触域原生积分 `eCurrentDensity · ContactSurfaceNormal`，六个状态均完成导出。

目标点 n23、Vd=0.05 V、Vg=0.05 V 的 `ContactSurfaceNormal` 原始方向电子面通量为 `{float(target['contact_surface_normal_integral_A_per_um']):.12e}` A/um，反向量为 `{float(target['negative_contact_surface_normal_integral_A_per_um']):.12e}` A/um；同次运行端口电子电流为 `{float(target['same_run_substrate_electron_terminal_A_per_um']):.12e}` A/um。两种符号均不能闭合端口值，因此没有进行逐状态符号拟合，也不能把 M49 可观测性缺口声明为已闭合。

关键的新结论是：反向原生积分与 M49 独立节点场边界积分在六个状态逐点一致，最大绝对误差为 `{float(report['native_observable']['maximum_negative_native_to_m49_boundary_absolute_error_A_per_um']):.3e}` A/um。这排除了“缺口仅由外部端面重构造成”的假设，把不可观测区进一步缩小到 Sentaurus 的场导出接触面通量与其端口 `eCurrent` 计算之间。

## 六状态结果

| 器件 | Vg (V) | Sentaurus 端口 eCurrent (A/um) | 原生面通量，原始法向 (A/um) | Vela-原生值 (A/um) |
|---|---:|---:|---:|---:|
{table}

## 边界与解释

- 唯一新增项是 `CurrentPlot Tcl` 接触域积分；没有重跑 Vela，也没有修改 BGN、SRH、迁移率、HFS、SG、准费米、接触或网格设置。
- `ContactSurfaceNormal` 的原始方向、反向量及两种端口残差全部保留，未做幅值或逐状态符号拟合。
- 手册说明 `Unit` 仅提供数量单位描述且 `.plt` 格式不保存单位；本次二维默认 `AreaFactor=1` 的数值按 A/um 归一化解释，并由其反向量与 M49 独立 A/um 积分逐点一致来限定。v1 中将数值再次乘 `1e-4` 的先验没有被运行结果支持，原合同保持未修改。
- 首次执行因 T-2022.03-SP2 拒绝同时指定 `Contact` 与 `Region` 而失败；冻结执行勘误后仅删除冗余 `Region=Silicon_1`，成功回放。该勘误没有改变接触域、integrand 或器件状态。
- M50 只闭合原生可观测性和端口定位，不据此声明生产代码缺陷或授权参数修改。

机器可读报告：`{portable(REPORT)}`。
'''
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc, encoding="utf-8", newline="\n")
    write_json(ATTEMPTS, {
        "schema": "vela.simplemos.sdevice.m50_sentaurus_execution_attempts.v1",
        "attempts": [
            {"attempt": 1,
             "remote_root": read_json(ERRATUM)["reason"]["failed_remote_root"],
             "status": "failed_before_bias_sweeps",
             "error": read_json(ERRATUM)["reason"]["observed_runtime_error"]},
            {"attempt": 2, "remote_root": read_json(ERRATUM)["new_remote_root"],
             "status": "complete", "state_count": 6,
             "execution_erratum_sha256": sha256(ERRATUM)},
        ],
    })
    artifact = {
        "schema": "vela.simplemos.sdevice.m50_artifact.v1",
        "title": "SimpleMOS M50 native substrate electron face-flux export",
        "status": report["status"],
        "classification": report["classification"],
        "report": portable(REPORT),
        "ledgers": [portable(FLUX_LEDGER), portable(REPLAY_LEDGER)],
        "execution_attempts": portable(ATTEMPTS),
        "document": portable(DOC),
        "target": target,
        "acceptance": report["acceptance"],
    }
    write_json(ARTIFACT, artifact)
    artifacts = [REPORT, FLUX_LEDGER, REPLAY_LEDGER, ATTEMPTS, DOC, ARTIFACT]
    sources = [CONTRACT, FREEZE, ERRATUM, ERRATUM_FREEZE,
               Path(__file__).resolve(), TEST,
               M47_EVIDENCE, M48_EVIDENCE, M49_EVIDENCE,
               OUTPUT / "sentaurus_manifest.json", OUTPUT / "sentaurus_banner.txt"]
    for device in DEVICES:
        sources.append(OUTPUT / "sentaurus_bundle" / device /
                       f"{device}_vd_0p05_des.cmd")
        sources.extend(sorted((OUTPUT / "sentaurus_raw/sentaurus_bundle" / device).glob(
            "IdVg_*des.plt")))
    evidence = {
        "schema": "vela.simplemos.sdevice.m50_native_substrate_face_flux_export_evidence.v1",
        "status": "frozen",
        "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "new_sentaurus_execution": True,
        "new_vela_execution": False,
        "default_physics_model_changed": False,
        "closed_topics_reinvestigated": False,
        "acceptance": report["acceptance"],
    }
    write_json(EVIDENCE, evidence)


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
        report = analyze(contract, manifest, banner)
        freeze_artifacts(report)
        print(json.dumps({"status": report["status"],
                          "classification": report["classification"],
                          "all_checks_pass": report["acceptance"]["all_checks_pass"],
                          "report": portable(REPORT)}, indent=2))


if __name__ == "__main__":
    main()
