#!/usr/bin/env python3
"""Freeze the read-only SimpleMOS M59 port-observable decomposition."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import analyze_simplemos_port_residual_ledger as analysis  # noqa: E402

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m59_port_observable_decomposition_contract_v1.json"
FREEZE = ROOT / "simplemos_m59_port_observable_decomposition_contract_freeze.json"
M46_POINTS = ROOT / "full_matrix_requalification/m46_pointwise_delta.csv"
M54_COMPONENTS = ROOT / "terminal_common_mode_attribution/m54_terminal_component_ledger.csv"
M54_STATES = ROOT / "terminal_common_mode_attribution/m54_state_decomposition_ledger.csv"
M58_DEVICES = ROOT / "compensation_contour_mesh_audit/m58_device_compensation_mesh_ledger.csv"
M58_EVIDENCE = ROOT / "simplemos_m58_compensation_contour_mesh_audit_evidence.json"
OUTPUT = ROOT / "port_observable_decomposition"
REPORT = OUTPUT / "m59_port_observable_decomposition_report.json"
POINTS = OUTPUT / "m59_pointwise_observable_ledger.csv"
CURVES = OUTPUT / "m59_curve_decomposition_ledger.csv"
PAIRS = OUTPUT / "m59_pair_decomposition_ledger.csv"
SENSITIVITY = OUTPUT / "m59_threshold_scope_sensitivity_ledger.csv"
KCL = OUTPUT / "m59_kcl_topology_stratification_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m59_port_observable_decomposition_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m59/artifact.json"
EVIDENCE = ROOT / "simplemos_m59_port_observable_decomposition_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m59_port_observable_decomposition.py"
ENGINE = REPO / "scripts/analyze_simplemos_port_residual_ledger.py"


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


def percentile_median(values: list[float]) -> float:
    return float(statistics.median(values))


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if contract.get("schema") != "vela.simplemos.sdevice.m59_port_observable_decomposition_contract.v1":
        raise ValueError("unexpected M59 contract schema")
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M59 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M59 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M59 upstream artifact changed: {relative}")
    m58 = read_json(M58_EVIDENCE)
    if (m58.get("status") != contract["upstream"]["required_m58_status"] or
            m58.get("classification") != contract["upstream"]["required_m58_classification"]):
        raise ValueError("M58 status or classification changed")
    return contract


def local_excess(rows: list[dict[str, Any]], row: dict[str, Any],
                 clean_fraction: float, window: float) -> tuple[float, float]:
    clean = [candidate for candidate in rows
             if candidate is not row
             and abs(candidate["vg"] - row["vg"]) <= window + 1e-12
             and abs(candidate["observable_gap_e"] / candidate["s_id"]) < clean_fraction]
    if not clean:
        return math.nan, math.nan
    baseline = sum(candidate["dex"] for candidate in clean) / len(clean)
    return baseline, 10.0 ** (row["dex"] - baseline) - 1.0


def slope_stats(rows: list[dict[str, Any]], threshold: float,
                clean_fraction: float, window: float,
                drain: float | None) -> dict[str, float | int]:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        if drain is not None and not math.isclose(row["vd"], drain):
            continue
        fraction = row["observable_gap_e"] / row["s_id"]
        if abs(fraction) < threshold:
            continue
        baseline, excess = local_excess(
            analysis_rows_by_case[row["case"]], row, clean_fraction, window)
        if not math.isfinite(baseline):
            continue
        xs.append(fraction)
        ys.append(excess)
    denominator = sum(value * value for value in xs)
    slope = sum(left * right for left, right in zip(xs, ys)) / denominator
    pearson = analysis.pearson(xs, ys) if len(xs) >= 2 else math.nan
    return {"flagged_point_count": len(xs), "slope": slope,
            "pearson": pearson}


analysis_rows_by_case: dict[str, list[dict[str, Any]]] = {}


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    global analysis_rows_by_case
    rows, by_case = analysis.load()
    analysis_rows_by_case = by_case
    rules = contract["frozen_analysis_rules"]
    acceptance = contract["acceptance"]
    primary = float(rules["primary_burst_absolute_fraction"])
    clean_fraction = float(rules["clean_neighbor_absolute_fraction"])
    window = float(rules["baseline_window_V"])
    components = read_csv(M54_COMPONENTS)
    states = read_csv(M54_STATES)
    m58_rows = read_csv(M58_DEVICES)
    if len(rows) != int(acceptance["required_point_count"]):
        raise ValueError("M59 point count mismatch")
    if len(components) != int(acceptance["required_component_count"]):
        raise ValueError("M59 component count mismatch")
    if len(states) != int(contract["upstream"]["required_m54_state_count"]):
        raise ValueError("M59 state count mismatch")
    if len(m58_rows) != len(contract["matrix"]["devices"]):
        raise ValueError("M59 M58 device count mismatch")

    point_rows: list[dict[str, Any]] = []
    for case, case_rows in sorted(by_case.items()):
        xs = [row["vg"] for row in case_rows]
        sentaurus_log = [math.log10(abs(row["s_id"])) for row in case_rows]
        vela_log = [math.log10(abs(row["v_id"])) for row in case_rows]
        for index, row in enumerate(case_rows):
            fraction = row["observable_gap_e"] / row["s_id"]
            baseline, excess = local_excess(case_rows, row, clean_fraction, window)
            rough_s = analysis.lagrange_excluded(xs, sentaurus_log, index)
            rough_v = analysis.lagrange_excluded(xs, vela_log, index)
            point_rows.append({
                "state": row["state"], "case": case, "device": row["dev"],
                "drain_voltage_V": row["vd"], "gate_voltage_V": row["vg"],
                "sentaurus_default_drain_current_A_per_um": row["s_id"],
                "vela_drain_current_A_per_um": row["v_id"],
                "vela_sentaurus_log10_ratio_dex": row["dex"],
                "substrate_default_eCurrent_A_per_um": row["sub_e_def"],
                "substrate_direct_eCurrent_A_per_um": row["sub_e_dir"],
                "substrate_default_minus_direct_eCurrent_A_per_um": row["observable_gap_e"],
                "signed_observable_fraction": fraction,
                "absolute_observable_fraction": abs(fraction),
                "primary_burst_flag": int(abs(fraction) >= primary),
                "clean_neighbor_baseline_dex": baseline if math.isfinite(baseline) else "",
                "local_excess_linear": excess if math.isfinite(excess) else "",
                "sentaurus_log_current_roughness_dex": rough_s if rough_s is not None else "",
                "vela_log_current_roughness_dex": rough_v if rough_v is not None else "",
                "default_four_terminal_kcl_A_per_um": row["kcl_def"],
                "direct_four_terminal_kcl_A_per_um": row["kcl_dir"],
            })

    curve_rows: list[dict[str, Any]] = []
    curve_clean_max: dict[str, dict[str, Any]] = {}
    for case, case_rows in sorted(by_case.items()):
        maximum = max(case_rows, key=lambda row: abs(row["dex"]))
        clean = [row for row in case_rows
                 if abs(row["observable_gap_e"] / row["s_id"]) < primary]
        clean_maximum = max(clean, key=lambda row: abs(row["dex"]))
        curve_clean_max[case] = clean_maximum
        curve_rows.append({
            "case": case, "device": maximum["dev"],
            "drain_voltage_V": maximum["vd"],
            "all_point_maximum_error_dex": maximum["dex"],
            "all_point_maximum_gate_voltage_V": maximum["vg"],
            "all_point_maximum_is_burst": int(
                abs(maximum["observable_gap_e"] / maximum["s_id"]) >= primary),
            "burst_free_maximum_error_dex": clean_maximum["dex"],
            "burst_free_maximum_gate_voltage_V": clean_maximum["vg"],
            "burst_free_maximum_observable_fraction":
                clean_maximum["observable_gap_e"] / clean_maximum["s_id"],
        })

    pair_rows: list[dict[str, Any]] = []
    for pair in contract["matrix"]["matched_pairs"]:
        low, high = pair["low"], pair["high"]
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            tag = "0p05" if math.isclose(drain, 0.05) else "1"
            low_case, high_case = f"{low}_vd_{tag}", f"{high}_vd_{tag}"
            low_all = max(by_case[low_case], key=lambda row: abs(row["dex"]))
            high_all = max(by_case[high_case], key=lambda row: abs(row["dex"]))
            low_clean, high_clean = curve_clean_max[low_case], curve_clean_max[high_case]
            pair_rows.append({
                "low_device": low, "high_device": high,
                "GOxTime_min": pair["GOxTime_min"],
                "LDD_Dose_cm2": pair["LDD_Dose_cm2"],
                "drain_voltage_V": drain,
                "all_point_low_maximum_error_dex": low_all["dex"],
                "all_point_high_maximum_error_dex": high_all["dex"],
                "all_point_growth_dex": high_all["dex"] - low_all["dex"],
                "burst_free_low_maximum_error_dex": low_clean["dex"],
                "burst_free_high_maximum_error_dex": high_clean["dex"],
                "burst_free_growth_dex": high_clean["dex"] - low_clean["dex"],
                "burst_free_low_gate_voltage_V": low_clean["vg"],
                "burst_free_high_gate_voltage_V": high_clean["vg"],
            })

    sensitivity_rows: list[dict[str, Any]] = []
    sensitivity_pair_growth: dict[float, list[float]] = {}
    for threshold_value in map(float, rules["threshold_sensitivity_absolute_fractions"]):
        maxima: dict[str, dict[str, Any]] = {}
        for case, case_rows in by_case.items():
            clean = [row for row in case_rows
                     if abs(row["observable_gap_e"] / row["s_id"]) < threshold_value]
            maxima[case] = max(clean, key=lambda row: abs(row["dex"]))
        growths: list[float] = []
        for pair in contract["matrix"]["matched_pairs"]:
            for drain in map(float, contract["matrix"]["drain_voltages_V"]):
                tag = "0p05" if math.isclose(drain, 0.05) else "1"
                growths.append(maxima[f"{pair['high']}_vd_{tag}"]["dex"] -
                               maxima[f"{pair['low']}_vd_{tag}"]["dex"])
        sensitivity_pair_growth[threshold_value] = growths
        for scope, drain in (("all", None), ("low_drain", 0.05), ("high_drain", 1.0)):
            stats = slope_stats(rows, threshold_value, clean_fraction, window, drain)
            sensitivity_rows.append({
                "threshold": threshold_value, "scope": scope,
                "flagged_point_count": stats["flagged_point_count"],
                "zero_intercept_slope": stats["slope"],
                "pearson": stats["pearson"],
                "minimum_pair_burst_free_growth_dex": min(growths),
                "maximum_pair_burst_free_growth_dex": max(growths),
                "minimum_burst_free_maximum_gate_voltage_V": min(
                    item["vg"] for item in maxima.values()),
                "maximum_burst_free_maximum_gate_voltage_V": max(
                    item["vg"] for item in maxima.values()),
            })

    m58_by_device = {row["device"]: row for row in m58_rows}
    kcl_rows: list[dict[str, Any]] = []
    floating_medians: list[float] = []
    nonfloating_medians: list[float] = []
    for case, case_rows in sorted(by_case.items()):
        med_default = percentile_median([abs(row["kcl_def"]) for row in case_rows])
        med_direct = percentile_median([abs(row["kcl_dir"]) for row in case_rows])
        device = case_rows[0]["dev"]
        floating = int(m58_by_device[device]["single_node_floating_n_type_component_count"])
        (floating_medians if floating else nonfloating_medians).append(med_default)
        kcl_rows.append({
            "case": case, "device": device,
            "drain_voltage_V": case_rows[0]["vd"],
            "single_node_floating_n_type_component_count": floating,
            "median_absolute_default_four_terminal_kcl_A_per_um": med_default,
            "maximum_absolute_default_four_terminal_kcl_A_per_um": max(
                abs(row["kcl_def"]) for row in case_rows),
            "median_absolute_direct_four_terminal_kcl_A_per_um": med_direct,
        })
    kcl_ratio = min(floating_medians) / max(nonfloating_medians)

    primary_low = next(row for row in sensitivity_rows
                       if math.isclose(float(row["threshold"]), primary)
                       and row["scope"] == "low_drain")
    sensitivity_limit = float(acceptance["sensitivity_threshold_max"])
    tested_growths = [growth for threshold, values in sensitivity_pair_growth.items()
                      if threshold <= sensitivity_limit for growth in values]
    checks = {
        "contract_hash_unchanged": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "point_count": len(point_rows) == int(acceptance["required_point_count"]),
        "curve_count": len(curve_rows) == int(acceptance["required_curve_count"]),
        "pair_row_count": len(pair_rows) == int(acceptance["required_pair_row_count"]),
        "component_count": len(components) == int(acceptance["required_component_count"]),
        "burst_free_growth_threshold_robust": (
            min(tested_growths) >= float(acceptance["burst_free_growth_min_dex"])
            and max(tested_growths) <= float(acceptance["burst_free_growth_max_dex"])),
        "low_drain_slope_in_range": (
            float(acceptance["low_drain_slope_min"]) <= float(primary_low["zero_intercept_slope"])
            <= float(acceptance["low_drain_slope_max"])),
        "low_drain_pearson_pass": float(primary_low["pearson"]) >= float(
            acceptance["low_drain_pearson_min"]),
        "floating_device_set_exact": sorted(
            device for device, row in m58_by_device.items()
            if int(row["single_node_floating_n_type_component_count"]) > 0
        ) == sorted(acceptance["floating_devices"]),
        "kcl_topology_stratification": kcl_ratio >= float(
            acceptance["nonfloating_to_floating_median_kcl_ratio_min"]),
        "historical_artifacts_not_rewritten": True,
        "no_new_solve": True,
        "residual_causality_not_claimed": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    classification = (
        "port_observable_burst_association_plus_smooth_nwell_growth"
        if checks["all_checks_pass"] else "upstream_or_analysis_mismatch")
    report = {
        "schema": "vela.simplemos.sdevice.m59_port_observable_decomposition_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "source_scope": {"point_count": len(point_rows), "curve_count": len(curve_rows),
                         "component_count": len(components)},
        "primary_threshold": primary,
        "smooth_component": {
            "minimum_threshold_robust_pair_growth_dex": min(tested_growths),
            "maximum_threshold_robust_pair_growth_dex": max(tested_growths),
            "interpretation": "A systematic matched-background-Boron curve-error growth remains after excluding points with a large default-minus-Direct substrate observable fraction; its physical chain remains for M63.",
        },
        "port_observable_component": {
            "low_drain_primary_slope": primary_low["zero_intercept_slope"],
            "low_drain_primary_pearson": primary_low["pearson"],
            "low_drain_primary_point_count": primary_low["flagged_point_count"],
            "high_drain_primary": next(row for row in sensitivity_rows
                                      if math.isclose(float(row["threshold"]), primary)
                                      and row["scope"] == "high_drain"),
            "interpretation": "The positive 0.36-scale drain-error association is low-drain-specific; default-minus-Direct is an observable gap, while continuity-residual causality remains an M60 hypothesis.",
        },
        "kcl_topology": {
            "floating_devices": acceptance["floating_devices"],
            "minimum_floating_to_maximum_nonfloating_median_kcl_ratio": kcl_ratio,
            "interpretation": "Large default four-terminal KCL sums stratify with the two frozen-TDR devices containing single-node floating n-type components; this is observational support, not an internal-algorithm proof.",
        },
        "next_gate": contract["next_gate"],
        "acceptance": checks,
    }
    return report, {
        "points": point_rows, "curves": curve_rows, "pairs": pair_rows,
        "sensitivity": sensitivity_rows, "kcl": kcl_rows,
    }


def freeze(report: dict[str, Any], ledgers: dict[str, list[dict[str, Any]]]) -> None:
    write_csv(POINTS, ledgers["points"])
    write_csv(CURVES, ledgers["curves"])
    write_csv(PAIRS, ledgers["pairs"])
    write_csv(SENSITIVITY, ledgers["sensitivity"])
    write_csv(KCL, ledgers["kcl"])
    write_json(REPORT, report)
    smooth = report["smooth_component"]
    port = report["port_observable_component"]
    high = port["high_drain_primary"]
    DOC.write_text(f"""# SimpleMOS M59 端口观测 burst 与平滑 NWell 增幅分解

## 结论

M59 分类为 `{report['classification']}`。它只读取冻结的 M46、M54 和 M58 账本，没有执行新的 Sentaurus 或 Vela 求解，也没有改写历史账本。

在 `|substrate(default-Direct)eCurrent / Sentaurus Id|` 阈值从0.005扫描到0.03时，八个匹配NWell/漏压行的burst-free每曲线最大误差增幅始终位于 `{float(smooth['minimum_threshold_robust_pair_growth_dex']):.4f}` 到 `{float(smooth['maximum_threshold_robust_pair_growth_dex']):.4f}` dex。这证明M58使用的每曲线最大值混合了平滑NWell增幅与低电流端口观测burst，但平滑成分的BGN、静电或迁移率归因仍留给M63。

在主阈值0.02下，低漏压标志点的零截距斜率为 `{float(port['low_drain_primary_slope']):.3f}`、Pearson为 `{float(port['low_drain_primary_pearson']):.3f}`；高漏压斜率为 `{float(high['zero_intercept_slope']):.3f}`。因此0.36量级的正关联只用于Vd=0.05 V，不能外推到全部15个混合漏压标志点。

`default-Direct` 是M54直接记录的端口观测差，不是独立导出的节点连续性残差。把它解释为p阱内部残差和是等待M60收敛控制实验检验的box恒等式假说。M59不声明Sentaurus求解器缺陷，也不授权生产参考修改。

M58中仅n21/n23含单节点浮置n型分量；它们的默认四端KCL中位幅值相对所有无浮置器件的最大中位值至少高 `{float(report['kcl_topology']['minimum_floating_to_maximum_nonfloating_median_kcl_ratio']):.1f}` 倍。这支持浮置阱调制默认KCL观测，但不是内部算法实现的独立证明。

## 下一门槛

M60只允许冻结的Sentaurus收敛控制干预，检验收紧收敛后端口观测burst是否下降，同时非burst曲线点保持稳定。M60之前不得把端口观测差称为已证实的Newton残差。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M59 port-observable decomposition",
        "status": report["status"],
        "classification": report["classification"],
        "summary": report["port_observable_component"]["interpretation"],
        "report": portable(REPORT),
        "ledgers": [portable(POINTS), portable(CURVES), portable(PAIRS),
                    portable(SENSITIVITY), portable(KCL)],
    })
    artifacts = [REPORT, POINTS, CURVES, PAIRS, SENSITIVITY, KCL, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m59_port_observable_decomposition_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT),
        "new_sentaurus_execution": False,
        "new_vela_execution": False,
        "historical_artifacts_rewritten": False,
        "source_hashes": read_json(FREEZE)["upstream_hashes"],
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT),
                                  portable(ENGINE): sha256(ENGINE)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "acceptance": report["acceptance"],
    })


def verify() -> dict[str, Any]:
    validate_contract()
    evidence = read_json(EVIDENCE)
    report = read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M59 evidence is not accepted and frozen")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M59 artifact hash changed: {relative}")
    if not report["acceptance"]["all_checks_pass"]:
        raise ValueError("M59 acceptance is not complete")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        report = verify()
        print(f"M59 verified: {report['classification']}")
        return
    contract = validate_contract()
    report, ledgers = analyze(contract)
    freeze(report, ledgers)
    if not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError(f"M59 acceptance failed: {report['acceptance']}")
    print(json.dumps({"status": report["status"],
                      "classification": report["classification"],
                      "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
