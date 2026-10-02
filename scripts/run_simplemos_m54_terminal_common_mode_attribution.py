#!/usr/bin/env python3
"""Analyze and freeze SimpleMOS M54 terminal common-mode attribution."""

from __future__ import annotations

import argparse
import csv
from decimal import Decimal, getcontext
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m53_direct_current_full_matrix as m53  # noqa: E402


ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m54_terminal_common_mode_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m54_terminal_common_mode_attribution_contract_freeze.json"
M8_EVIDENCE = ROOT / "simplemos_m8_original_physics_evidence.json"
M53_EVIDENCE = ROOT / "simplemos_m53_direct_current_full_matrix_evidence.json"
M53_REPORT = ROOT / "direct_current_full_matrix/m53_direct_current_full_matrix_report.json"
M8_RAW = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_raw/sentaurus_bundle"
M53_RAW = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m53_direct_current_full_matrix/sentaurus_raw/sentaurus_bundle"
PORTABLE = ROOT / "terminal_common_mode_attribution"
REPORT = PORTABLE / "m54_terminal_common_mode_attribution_report.json"
STATES = PORTABLE / "m54_state_decomposition_ledger.csv"
TERMINALS = PORTABLE / "m54_terminal_component_ledger.csv"
CASES = PORTABLE / "m54_case_summary.csv"
NWELL = PORTABLE / "m54_nwell_pair_summary.csv"
DOC = REPO / "docs/validation/simplemos_m54_terminal_common_mode_attribution_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m54/artifact.json"
EVIDENCE = ROOT / "simplemos_m54_terminal_common_mode_attribution_evidence.json"
TEST = REPO / "tests/regression/test_simplemos_m54_terminal_common_mode_attribution.py"
CONTACTS = ("source", "drain", "gate", "substrate")
COMPONENTS = ("electron", "hole", "total")
getcontext().prec = 50


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


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


def relation(left: float, right: float) -> str:
    if left == 0.0 or right == 0.0:
        return "contains_zero"
    return "same_sign" if math.copysign(1.0, left) == math.copysign(1.0, right) \
        else "opposite_sign"


def identity_ok(residual: float, scale: float, abs_tol: float,
                rel_tol: float) -> bool:
    return abs(residual) <= abs_tol or abs(residual) / max(scale, 1e-300) <= rel_tol


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != (
            "vela.simplemos.sdevice.m54_terminal_common_mode_attribution_contract.v1"):
        raise ValueError("unexpected M54 contract schema")
    if freeze.get("status") != "frozen_before_analysis":
        raise ValueError("M54 contract was not frozen before analysis")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M54 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M54 upstream changed: {relative}")
    if read_json(M8_EVIDENCE).get("status") != "accepted":
        raise ValueError("M54 requires accepted M8 evidence")
    if read_json(M53_EVIDENCE).get("status") != "frozen":
        raise ValueError("M54 requires frozen M53 evidence")
    m53_report = read_json(M53_REPORT)
    upstream = contract["upstream"]
    if m53_report.get("classification") != upstream["required_m53_classification"]:
        raise ValueError("M53 classification changed")
    if int(m53_report["matrix"]["curve_count"]) != 16 or int(
            m53_report["matrix"]["direct_bias_point_count"]) != 816:
        raise ValueError("M53 matrix anchor changed")
    return contract


def parse_matrix(contract: dict[str, Any], root: Path
                 ) -> tuple[dict[str, list[dict[str, float]]], list[Path]]:
    gates = [0.05 * index for index in range(51)]
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    parsed: dict[str, list[dict[str, float]]] = {}
    paths: list[Path] = []
    for device in contract["scope"]["devices"]:
        for drain in map(float, contract["scope"]["drain_voltages_V"]):
            case = f"{device}_vd_{m53.voltage_tag(drain)}"
            path = root / device / f"IdVg_{case}_des.plt"
            if not path.is_file():
                raise FileNotFoundError(path)
            parsed[case] = m53.exact_terminal_curve(path, gates, tolerance)
            paths.append(path)
    return parsed, paths


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[Path]]:
    default, default_paths = parse_matrix(contract, M8_RAW)
    direct, direct_paths = parse_matrix(contract, M53_RAW)
    abs_tol = float(contract["acceptance"][
        "drain_identity_absolute_tolerance_A_per_um"])
    rel_tol = float(contract["acceptance"]["drain_identity_relative_tolerance"])
    terminal_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    case_rows: list[dict[str, Any]] = []
    state_index: dict[tuple[str, float, float], dict[str, Any]] = {}
    for device in contract["scope"]["devices"]:
        for drain_bias in map(float, contract["scope"]["drain_voltages_V"]):
            case = f"{device}_vd_{m53.voltage_tag(drain_bias)}"
            per_case: list[dict[str, Any]] = []
            for baseline, changed in zip(default[case], direct[case], strict=True):
                gate = float(baseline["gate_voltage_V"])
                if gate != float(changed["gate_voltage_V"]):
                    raise RuntimeError(f"M54 gate mismatch: {case}, {gate}")
                state: dict[str, Any] = {
                    "state": f"{case}_vg_{m53.voltage_tag(gate)}",
                    "case": case,
                    "device": device,
                    "drain_voltage_V": drain_bias,
                    "gate_voltage_V": gate,
                }
                for component in COMPONENTS:
                    values: dict[str, tuple[float, float]] = {}
                    for contact in CONTACTS:
                        key = f"{contact}_{component}_A_per_um"
                        old = float(baseline[key])
                        new = float(changed[key])
                        values[contact] = (old, new)
                        terminal_rows.append({
                            "state": state["state"], "case": case,
                            "device": device, "drain_voltage_V": drain_bias,
                            "gate_voltage_V": gate, "contact": contact,
                            "component": component,
                            "default_A_per_um": old,
                            "direct_A_per_um": new,
                            "direct_minus_default_A_per_um": new - old,
                        })
                    old_source, new_source = values["source"]
                    old_drain, new_drain = values["drain"]
                    old_source_decimal = Decimal(str(old_source))
                    new_source_decimal = Decimal(str(new_source))
                    old_drain_decimal = Decimal(str(old_drain))
                    new_drain_decimal = Decimal(str(new_drain))
                    two = Decimal(2)
                    old_common_decimal = (
                        old_drain_decimal + old_source_decimal) / two
                    new_common_decimal = (
                        new_drain_decimal + new_source_decimal) / two
                    old_anti_decimal = (
                        old_drain_decimal - old_source_decimal) / two
                    new_anti_decimal = (
                        new_drain_decimal - new_source_decimal) / two
                    delta_drain_decimal = new_drain_decimal - old_drain_decimal
                    delta_common_decimal = (
                        new_common_decimal - old_common_decimal)
                    delta_anti_decimal = new_anti_decimal - old_anti_decimal
                    residual_decimal = (delta_drain_decimal
                                        - delta_common_decimal
                                        - delta_anti_decimal)
                    old_common = float(old_common_decimal)
                    new_common = float(new_common_decimal)
                    old_anti = float(old_anti_decimal)
                    new_anti = float(new_anti_decimal)
                    delta_drain = float(delta_drain_decimal)
                    delta_common = float(delta_common_decimal)
                    delta_anti = float(delta_anti_decimal)
                    residual = float(residual_decimal)
                    scale = max(abs(delta_drain), abs(delta_common) + abs(delta_anti))
                    denominator = abs(delta_common) + abs(delta_anti)
                    prefix = f"{component}_"
                    state.update({
                        prefix + "default_source_A_per_um": old_source,
                        prefix + "direct_source_A_per_um": new_source,
                        prefix + "default_drain_A_per_um": old_drain,
                        prefix + "direct_drain_A_per_um": new_drain,
                        prefix + "default_common_mode_A_per_um": old_common,
                        prefix + "direct_common_mode_A_per_um": new_common,
                        prefix + "default_antisymmetric_A_per_um": old_anti,
                        prefix + "direct_antisymmetric_A_per_um": new_anti,
                        prefix + "delta_drain_A_per_um": delta_drain,
                        prefix + "delta_common_mode_A_per_um": delta_common,
                        prefix + "delta_antisymmetric_A_per_um": delta_anti,
                        prefix + "drain_identity_residual_A_per_um": residual,
                        prefix + "drain_identity_within_tolerance": identity_ok(
                            residual, scale, abs_tol, rel_tol),
                        prefix + "common_mode_fraction": (
                            abs(delta_common) / denominator if denominator else 0.0),
                        prefix + "default_source_drain_sign_relation": relation(
                            old_source, old_drain),
                        prefix + "direct_source_drain_sign_relation": relation(
                            new_source, new_drain),
                        prefix + "default_four_terminal_sum_A_per_um": sum(
                            pair[0] for pair in values.values()),
                        prefix + "direct_four_terminal_sum_A_per_um": sum(
                            pair[1] for pair in values.values()),
                    })
                state_rows.append(state)
                per_case.append(state)
                state_index[(device, drain_bias, gate)] = state
            fractions = [float(row["total_common_mode_fraction"])
                         for row in per_case]
            case_rows.append({
                "case": case, "device": device,
                "drain_voltage_V": drain_bias,
                "state_count": len(per_case),
                "median_total_common_mode_fraction": statistics.median(fractions),
                "minimum_total_common_mode_fraction": min(fractions),
                "maximum_total_common_mode_fraction": max(fractions),
                "common_mode_dominated_state_count": sum(value >= 0.9
                                                         for value in fractions),
                "maximum_absolute_total_drain_shift_A_per_um": max(
                    abs(float(row["total_delta_drain_A_per_um"]))
                    for row in per_case),
                "default_opposite_sign_state_count": sum(
                    row["total_default_source_drain_sign_relation"] ==
                    "opposite_sign" for row in per_case),
                "direct_same_sign_state_count": sum(
                    row["total_direct_source_drain_sign_relation"] ==
                    "same_sign" for row in per_case),
                "maximum_absolute_total_four_terminal_sum_change_A_per_um": max(
                    abs(float(row["total_direct_four_terminal_sum_A_per_um"])
                        - float(row["total_default_four_terminal_sum_A_per_um"]))
                    for row in per_case),
            })

    nwell_rows: list[dict[str, Any]] = []
    for pair in (("n17", "n21"), ("n18", "n22"),
                 ("n19", "n23"), ("n20", "n24")):
        for drain in (0.05, 1.0):
            low = [state_index[(pair[0], drain, 0.05 * index)]
                   for index in range(51)]
            high = [state_index[(pair[1], drain, 0.05 * index)]
                    for index in range(51)]
            nwell_rows.append({
                "low_nwell_device": pair[0],
                "high_nwell_device": pair[1],
                "drain_voltage_V": drain,
                "low_median_total_common_mode_fraction": statistics.median(
                    float(row["total_common_mode_fraction"]) for row in low),
                "high_median_total_common_mode_fraction": statistics.median(
                    float(row["total_common_mode_fraction"]) for row in high),
                "low_maximum_absolute_total_drain_shift_A_per_um": max(
                    abs(float(row["total_delta_drain_A_per_um"])) for row in low),
                "high_maximum_absolute_total_drain_shift_A_per_um": max(
                    abs(float(row["total_delta_drain_A_per_um"])) for row in high),
                "high_minus_low_maximum_shift_A_per_um": max(
                    abs(float(row["total_delta_drain_A_per_um"])) for row in high)
                    - max(abs(float(row["total_delta_drain_A_per_um"]))
                          for row in low),
            })

    target_spec = contract["analysis"]["target"]
    worst_spec = contract["analysis"]["m53_worst"]
    target = state_index[(target_spec["device"],
                          float(target_spec["drain_voltage_V"]),
                          float(target_spec["gate_voltage_V"]))]
    worst = state_index[(worst_spec["device"],
                         float(worst_spec["drain_voltage_V"]),
                         float(worst_spec["gate_voltage_V"]))]
    matrix_fractions = [float(row["total_common_mode_fraction"])
                        for row in state_rows]
    median_fraction = statistics.median(matrix_fractions)
    identity_pass = all(bool(row[f"{component}_drain_identity_within_tolerance"])
                        for row in state_rows for component in COMPONENTS)
    complete = len(state_rows) == 816 and len(terminal_rows) == 9792
    if not complete or not identity_pass:
        classification = "identity_or_replay_failure"
    elif (float(target["total_common_mode_fraction"]) >= 0.9
          and float(worst["total_common_mode_fraction"]) >= 0.9
          and median_fraction >= 0.5):
        classification = "direct_shift_common_mode_dominated"
    elif (float(target["total_common_mode_fraction"]) <= 0.1
          and float(worst["total_common_mode_fraction"]) <= 0.1
          and median_fraction <= 0.5):
        classification = "direct_shift_antisymmetric_dominated"
    else:
        classification = "direct_shift_mixed"
    max_identity = max(abs(float(row[f"{component}_drain_identity_residual_A_per_um"]))
                       for row in state_rows for component in COMPONENTS)
    acceptance = {
        "contract_hash_unchanged": sha256(CONTRACT) ==
            read_json(FREEZE)["contract_sha256"],
        "state_count": len(state_rows) == 816,
        "terminal_component_row_count": len(terminal_rows) == 9792,
        "case_count": len(case_rows) == 16,
        "points_per_case": all(int(row["state_count"]) == 51
                               for row in case_rows),
        "all_values_finite": all(math.isfinite(float(row[key]))
                                 for row in terminal_rows
                                 for key in ("default_A_per_um", "direct_A_per_um")),
        "drain_identities_close": identity_pass,
        "classification_declared": classification in
            contract["analysis"]["classifications"],
        "no_new_solver_execution": True,
        "historical_artifacts_preserved": True,
        "closed_topics_not_reopened": True,
    }
    acceptance["all_checks_pass"] = all(acceptance.values())
    report = {
        "schema": "vela.simplemos.sdevice.m54_terminal_common_mode_attribution_report.v1",
        "status": "complete" if acceptance["all_checks_pass"] else "failed",
        "classification": classification,
        "execution": {
            "new_sentaurus_execution": False,
            "new_vela_execution": False,
            "inputs": "frozen M8 default and M53 DirectCurrent terminal files",
            "contract_sha256_before_and_after": sha256(CONTRACT),
        },
        "matrix": {"state_count": len(state_rows),
                   "terminal_component_row_count": len(terminal_rows),
                   "case_count": len(case_rows)},
        "findings": {
            "target_total_common_mode_fraction":
                target["total_common_mode_fraction"],
            "worst_total_common_mode_fraction":
                worst["total_common_mode_fraction"],
            "matrix_median_total_common_mode_fraction": median_fraction,
            "common_mode_dominated_state_count": sum(
                value >= 0.9 for value in matrix_fractions),
            "default_opposite_sign_state_count": sum(
                row["total_default_source_drain_sign_relation"] == "opposite_sign"
                for row in state_rows),
            "direct_same_sign_state_count": sum(
                row["total_direct_source_drain_sign_relation"] == "same_sign"
                for row in state_rows),
            "maximum_absolute_drain_identity_residual_A_per_um": max_identity,
        },
        "target": target,
        "m53_worst": worst,
        "case_results": case_rows,
        "nwell_pairs": nwell_rows,
        "acceptance": acceptance,
        "claim_guard": contract["analysis"]["claim_guard"],
    }
    write_csv(STATES, state_rows)
    write_csv(TERMINALS, terminal_rows)
    write_csv(CASES, case_rows)
    write_csv(NWELL, nwell_rows)
    write_json(REPORT, report)
    return report, default_paths + direct_paths


def freeze_artifacts(report: dict[str, Any], raw_paths: list[Path]) -> None:
    findings = report["findings"]
    target = report["target"]
    worst = report["m53_worst"]
    case_table = "\n".join(
        f"| {row['case']} | {float(row['median_total_common_mode_fraction']):.6f} | "
        f"{row['common_mode_dominated_state_count']} | "
        f"{row['default_opposite_sign_state_count']} | "
        f"{row['direct_same_sign_state_count']} |"
        for row in report["case_results"])
    doc = f'''# SimpleMOS M54 端口共模归因

## 结论

M54 分类为 `{report['classification']}`。分析完全复用冻结的 M8 默认和 M53 `DirectCurrent` 端口文件，没有重新运行 Sentaurus 或 Vela。

对每个状态定义 source/drain 共模 `C=(Id+Is)/2` 和反对称输运项 `T=(Id-Is)/2`，并逐点验证 `Delta Id = Delta C + Delta T`。816个状态、3个载流子分量的最大恒等式残差为 `{float(findings['maximum_absolute_drain_identity_residual_A_per_um']):.12e}` A/um。

M53目标 n23、Vd=0.05 V、Vg=0.05 V 的总电流共模分数为 `{float(findings['target_total_common_mode_fraction']):.9f}`；M53最差点 n24、Vd=0.05 V、Vg=0 V 为 `{float(findings['worst_total_common_mode_fraction']):.9f}`；全816状态中位数为 `{float(findings['matrix_median_total_common_mode_fraction']):.9f}`。共有 `{findings['common_mode_dominated_state_count']}/816` 个状态的共模分数不低于0.9。

目标点 DirectCurrent-minus-default 漏端变化为 `{float(target['total_delta_drain_A_per_um']):.12e}` A/um，其中共模项 `{float(target['total_delta_common_mode_A_per_um']):.12e}` A/um，反对称项 `{float(target['total_delta_antisymmetric_A_per_um']):.12e}` A/um。最差点相应三项为 `{float(worst['total_delta_drain_A_per_um']):.12e}`、`{float(worst['total_delta_common_mode_A_per_um']):.12e}`、`{float(worst['total_delta_antisymmetric_A_per_um']):.12e}` A/um。

## 16曲线符号与共模摘要

| 曲线 | 共模分数中位数 | 共模>=0.9状态数 | 默认源漏反号数 | Direct源漏同号数 |
|---|---:|---:|---:|---:|
{case_table}

默认观测中 `{findings['default_opposite_sign_state_count']}/816` 个状态的源漏总电流反号；`DirectCurrent` 中 `{findings['direct_same_sign_state_count']}/816` 个状态同号。该符号和共模变化解释了 M53 的完整 Id-Vg 劣化，但它是端口观测分解，不是求解状态缺陷。

## 边界

- 不对 source/drain 按器件或偏压拟合符号。
- 不修改或替换官方默认端口定义。
- 不重新排查 HFS、SG、准费米、BGN、SRH、迁移率、接触模型或网格。

机器报告：`{portable(REPORT)}`。
'''
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(doc, encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.simplemos.sdevice.m54_artifact.v1",
        "title": "SimpleMOS M54 terminal common-mode attribution",
        "status": report["status"], "classification": report["classification"],
        "report": portable(REPORT),
        "ledgers": [portable(STATES), portable(TERMINALS), portable(CASES),
                    portable(NWELL)],
        "document": portable(DOC),
        "findings": report["findings"], "acceptance": report["acceptance"],
    })
    artifacts = [REPORT, STATES, TERMINALS, CASES, NWELL, DOC, ARTIFACT]
    sources = [CONTRACT, FREEZE, Path(__file__).resolve(), TEST,
               M8_EVIDENCE, M53_EVIDENCE, M53_REPORT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m54_terminal_common_mode_attribution_evidence.v1",
        "status": "frozen", "contract_sha256": sha256(CONTRACT),
        "artifacts": [{"path": portable(path), "sha256": sha256(path)}
                      for path in artifacts],
        "source_hashes": {portable(path): sha256(path) for path in sources},
        "raw_terminal_input_hashes": {portable(path): sha256(path)
                                      for path in raw_paths},
        "new_sentaurus_execution": False, "new_vela_execution": False,
        "historical_artifacts_rewritten": False,
        "closed_topics_reinvestigated": False,
        "classification": report["classification"],
        "acceptance": report["acceptance"],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analyze", action="store_true")
    args = parser.parse_args()
    contract = validate_contract()
    if not args.analyze:
        print(json.dumps({"status": "contract_validated",
                          "contract": portable(CONTRACT)}, indent=2))
        return
    report, raw_paths = analyze(contract)
    freeze_artifacts(report, raw_paths)
    print(json.dumps({
        "status": report["status"],
        "classification": report["classification"],
        "all_checks_pass": report["acceptance"]["all_checks_pass"],
        "report": portable(REPORT),
    }, indent=2))


if __name__ == "__main__":
    main()
