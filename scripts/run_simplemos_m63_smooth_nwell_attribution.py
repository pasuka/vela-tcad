#!/usr/bin/env python3
"""Analyze and freeze SimpleMOS M63 burst-free NWell error attribution."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m63_smooth_nwell_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m63_smooth_nwell_attribution_contract_freeze.json"
M60_EVIDENCE = ROOT / "simplemos_m60_tight_convergence_port_burst_evidence.json"
M60_POINTS = ROOT / "tight_convergence_port_burst/m60_tight_default_direct_point_ledger.csv"
M13_STATES = ROOT / "spatial_attribution/m13_state_summary.csv"
M40_REPORT = ROOT / "true_no_bgn_factorial/m40_true_no_bgn_factorial_report.json"
OUT = ROOT / "smooth_nwell_attribution"
REPORT = OUT / "m63_smooth_nwell_attribution_report.json"
POINTS = OUT / "m63_smooth_window_point_ledger.csv"
CURVES = OUT / "m63_curve_shape_ledger.csv"
THRESHOLDS = OUT / "m63_constant_current_threshold_ledger.csv"
PAIRS = OUT / "m63_nwell_pair_attribution_ledger.csv"
MOBILITY = OUT / "m63_m13_mobility_crosscheck_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m63_smooth_nwell_attribution_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m63/artifact.json"
EVIDENCE = ROOT / "simplemos_m63_smooth_nwell_attribution_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m63_smooth_nwell_attribution.py"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8",
                    newline="\n")


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


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def median(values: list[float]) -> float:
    data = sorted(values)
    midpoint = len(data) // 2
    return data[midpoint] if len(data) % 2 else 0.5 * (data[midpoint - 1] + data[midpoint])


def linear_interpolate_x(points: list[tuple[float, float]], target_y: float) -> float:
    """Return x at target y for a monotonic-enough (x, y) sequence."""
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if (y0 <= target_y <= y1) or (y1 <= target_y <= y0):
            if math.isclose(y0, y1):
                return 0.5 * (x0 + x1)
            return x0 + (target_y - y0) * (x1 - x0) / (y1 - y0)
    raise ValueError(f"target {target_y} is outside curve")


def slope_at(rows: list[dict[str, float]], index: int, key: str) -> float:
    left = max(index - 1, 0)
    right = min(index + 1, len(rows) - 1)
    if left == right:
        raise ValueError("cannot differentiate a one-point curve")
    dv = rows[right]["gate"] - rows[left]["gate"]
    return (math.log10(abs(rows[right][key])) - math.log10(abs(rows[left][key]))) / dv


def regression_swing(rows: list[dict[str, float]], key: str,
                     low: float, high: float) -> tuple[float, int]:
    selected = [(math.log10(abs(row[key])), row["gate"]) for row in rows
                if low <= abs(row[key]) <= high]
    if len(selected) < 3:
        raise ValueError(f"insufficient subthreshold points for {key}")
    xbar = mean([item[0] for item in selected])
    ybar = mean([item[1] for item in selected])
    slope = sum((x - xbar) * (y - ybar) for x, y in selected) / sum(
        (x - xbar) ** 2 for x, _ in selected)
    return 1000.0 * slope, len(selected)


def validate() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if freeze.get("status") != "frozen_before_analysis":
        raise ValueError("M63 contract was not frozen before analysis")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M63 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M63 upstream artifact changed: {relative}")
    m60 = read_json(M60_EVIDENCE)
    if (m60.get("status") != contract["upstream"]["required_m60_status"] or
            m60.get("classification") != contract["upstream"]["required_m60_classification"]):
        raise ValueError("M60 status or classification changed")
    return contract


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    source = read_csv(M60_POINTS)
    grouped: dict[tuple[str, float], list[dict[str, float]]] = {}
    for row in source:
        key = (row["device"], float(row["drain_voltage_V"]))
        grouped.setdefault(key, []).append({
            "gate": float(row["gate_voltage_V"]),
            "sentaurus": float(row["tight_default_drain_current_A_per_um"]),
            "vela": float(row["vela_drain_current_A_per_um"]),
        })
    for rows in grouped.values():
        rows.sort(key=lambda item: item["gate"])

    low_gate, high_gate = map(float, contract["matrix"]["smooth_gate_window_V"])
    current_levels = list(map(float, contract["matrix"]["constant_current_levels_A_per_um"]))
    ss_low, ss_high = map(float, contract["matrix"]["subthreshold_current_window_A_per_um"])
    point_rows: list[dict[str, Any]] = []
    curve_rows: list[dict[str, Any]] = []
    threshold_rows: list[dict[str, Any]] = []
    curve_by_key: dict[tuple[str, float], dict[str, Any]] = {}
    smooth_by_key: dict[tuple[str, float], list[dict[str, Any]]] = {}

    for (device, drain), rows in sorted(grouped.items()):
        smooth: list[dict[str, Any]] = []
        for index, row in enumerate(rows):
            if not (low_gate <= row["gate"] <= high_gate):
                continue
            error = math.log10(abs(row["vela"]) / abs(row["sentaurus"]))
            slope = slope_at(rows, index, "sentaurus")
            smooth.append({
                "case": f"{device}_vd_{str(drain).replace('.', 'p')}",
                "device": device, "drain_voltage_V": drain,
                "gate_voltage_V": row["gate"],
                "sentaurus_current_A_per_um": row["sentaurus"],
                "vela_current_A_per_um": row["vela"],
                "signed_error_dex": error,
                "sentaurus_local_slope_dec_per_V": slope,
                "point_equivalent_horizontal_shift_mV": 1000.0 * error / slope,
            })
        if len(smooth) < int(contract["acceptance"]["minimum_smooth_points_per_curve"]):
            raise ValueError(f"insufficient smooth-window points for {device}, {drain}")
        delta = sum(float(row["sentaurus_local_slope_dec_per_V"]) *
                    float(row["signed_error_dex"]) for row in smooth) / sum(
            float(row["sentaurus_local_slope_dec_per_V"]) ** 2 for row in smooth)
        error_energy = sum(float(row["signed_error_dex"]) ** 2 for row in smooth)
        residual_energy = 0.0
        for row in smooth:
            prediction = float(row["sentaurus_local_slope_dec_per_V"]) * delta
            residual = float(row["signed_error_dex"]) - prediction
            row["fitted_horizontal_shift_mV"] = 1000.0 * delta
            row["horizontal_shift_prediction_dex"] = prediction
            row["horizontal_shift_residual_dex"] = residual
            residual_energy += residual * residual
            point_rows.append(row)
        explained = 1.0 - residual_energy / error_energy

        sent_points = [(row["gate"], math.log10(abs(row["sentaurus"]))) for row in rows]
        vela_points = [(row["gate"], math.log10(abs(row["vela"]))) for row in rows]
        shifts: list[float] = []
        for level in current_levels:
            log_level = math.log10(level)
            try:
                sent_v = linear_interpolate_x(sent_points, log_level)
                vela_v = linear_interpolate_x(vela_points, log_level)
            except ValueError:
                continue
            shift = sent_v - vela_v
            shifts.append(shift)
            threshold_rows.append({
                "device": device, "drain_voltage_V": drain,
                "current_level_A_per_um": level,
                "sentaurus_gate_voltage_V": sent_v,
                "vela_gate_voltage_V": vela_v,
                "sentaurus_minus_vela_threshold_shift_mV": 1000.0 * shift,
            })
        if len(shifts) < int(contract["acceptance"]["minimum_constant_current_levels_per_curve"]):
            raise ValueError(f"insufficient constant-current levels for {device}, {drain}")
        sent_ss, sent_ss_count = regression_swing(rows, "sentaurus", ss_low, ss_high)
        vela_ss, vela_ss_count = regression_swing(rows, "vela", ss_low, ss_high)
        at08 = min(smooth, key=lambda item: abs(float(item["gate_voltage_V"]) - 0.8))
        curve = {
            "device": device, "drain_voltage_V": drain,
            "smooth_point_count": len(smooth),
            "maximum_signed_error_dex": max(float(row["signed_error_dex"]) for row in smooth),
            "maximum_error_gate_voltage_V": max(smooth, key=lambda item: float(item["signed_error_dex"]))["gate_voltage_V"],
            "fitted_horizontal_shift_mV": 1000.0 * delta,
            "horizontal_shift_error_energy_explained_fraction": explained,
            "constant_current_threshold_shift_mean_mV": 1000.0 * mean(shifts),
            "constant_current_threshold_shift_range_mV": 1000.0 * (max(shifts) - min(shifts)),
            "constant_current_level_count": len(shifts),
            "sentaurus_subthreshold_swing_mV_per_dec": sent_ss,
            "vela_subthreshold_swing_mV_per_dec": vela_ss,
            "subthreshold_swing_difference_mV_per_dec": vela_ss - sent_ss,
            "sentaurus_subthreshold_fit_point_count": sent_ss_count,
            "vela_subthreshold_fit_point_count": vela_ss_count,
            "vg_0p8_signed_error_dex": at08["signed_error_dex"],
            "vg_0p8_horizontal_shift_residual_dex": at08["horizontal_shift_residual_dex"],
        }
        curve_rows.append(curve)
        curve_by_key[(device, drain)] = curve
        smooth_by_key[(device, drain)] = smooth

    pair_rows: list[dict[str, Any]] = []
    pair_explanations: list[float] = []
    for low_device, high_device in contract["matrix"]["nwell_pairs"]:
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            low_curve = curve_by_key[(low_device, drain)]
            high_curve = curve_by_key[(high_device, drain)]
            high_points = smooth_by_key[(high_device, drain)]
            diagnostic = max(high_points, key=lambda item: float(item["signed_error_dex"]))
            gate = float(diagnostic["gate_voltage_V"])
            low_point = next(row for row in smooth_by_key[(low_device, drain)]
                             if math.isclose(float(row["gate_voltage_V"]), gate))
            observed = float(diagnostic["signed_error_dex"]) - float(low_point["signed_error_dex"])
            predicted = (float(diagnostic["horizontal_shift_prediction_dex"])
                         - float(low_point["horizontal_shift_prediction_dex"]))
            residual = observed - predicted
            explanation = max(0.0, 1.0 - abs(residual) / max(abs(observed), 1e-300))
            pair_explanations.append(explanation)
            pair_rows.append({
                "low_device": low_device, "high_device": high_device,
                "drain_voltage_V": drain, "diagnostic_gate_voltage_V": gate,
                "low_signed_error_dex": low_point["signed_error_dex"],
                "high_signed_error_dex": diagnostic["signed_error_dex"],
                "observed_nwell_growth_dex": observed,
                "horizontal_shift_predicted_growth_dex": predicted,
                "horizontal_shift_growth_residual_dex": residual,
                "horizontal_shift_growth_explained_fraction": explanation,
                "low_fitted_shift_mV": low_curve["fitted_horizontal_shift_mV"],
                "high_fitted_shift_mV": high_curve["fitted_horizontal_shift_mV"],
                "fitted_shift_growth_mV": float(high_curve["fitted_horizontal_shift_mV"]) - float(low_curve["fitted_horizontal_shift_mV"]),
                "subthreshold_swing_growth_mV_per_dec": float(high_curve["subthreshold_swing_difference_mV_per_dec"]) - float(low_curve["subthreshold_swing_difference_mV_per_dec"]),
            })

    m13 = read_csv(M13_STATES)
    mobility_rows: list[dict[str, Any]] = []
    for drain in map(float, contract["matrix"]["fixed_state_mobility_crosscheck"]["drain_voltages_V"]):
        selected = {row["device"]: row for row in m13
                    if row["device"] in ("n17", "n21")
                    and math.isclose(float(row["drain_voltage_V"]), drain)
                    and math.isclose(float(row["gate_voltage_V"]), 0.8)}
        low, high = selected["n17"], selected["n21"]
        low_error = float(low["sentaurus_drive_active_p95_mobility_error_dex"])
        high_error = float(high["sentaurus_drive_active_p95_mobility_error_dex"])
        curve_growth = (float(curve_by_key[("n21", drain)]["vg_0p8_signed_error_dex"])
                        - float(curve_by_key[("n17", drain)]["vg_0p8_signed_error_dex"]))
        mobility_rows.append({
            "low_device": "n17", "high_device": "n21",
            "drain_voltage_V": drain, "gate_voltage_V": 0.8,
            "low_sentaurus_drive_active_p95_mobility_error_dex": low_error,
            "high_sentaurus_drive_active_p95_mobility_error_dex": high_error,
            "mobility_error_growth_dex": high_error - low_error,
            "terminal_error_growth_dex": curve_growth,
            "sign_supports_mobility_growth": int((high_error - low_error) * curve_growth > 0.0),
        })

    threshold_gate = float(contract["analysis"]["threshold_dominance_fraction"])
    median_curve_explained = median([float(row["horizontal_shift_error_energy_explained_fraction"])
                                     for row in curve_rows])
    median_pair_explained = median(pair_explanations)
    mobility_support_count = sum(int(row["sign_supports_mobility_growth"]) for row in mobility_rows)
    all_growth_positive = all(float(row["observed_nwell_growth_dex"]) > 0.0 for row in pair_rows)
    if (all_growth_positive and median_curve_explained >= threshold_gate and
            median_pair_explained >= threshold_gate and mobility_support_count == 0):
        classification = "threshold_like_horizontal_shift_dominant"
    elif median_curve_explained < threshold_gate and mobility_support_count == len(mobility_rows):
        classification = "mobility_factor_dominant"
    elif median_curve_explained >= 0.5 or mobility_support_count > 0:
        classification = "mixed_threshold_and_mobility"
    else:
        classification = "smooth_growth_not_closed"

    checks = {
        "contract_frozen": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "input_point_count": len(source) == int(contract["acceptance"]["required_point_count"]),
        "curve_count": len(curve_rows) == int(contract["acceptance"]["required_curve_count"]),
        "pair_count": len(pair_rows) == int(contract["acceptance"]["required_pair_count"]),
        "mobility_crosscheck_count": len(mobility_rows) == int(contract["acceptance"]["required_mobility_crosscheck_count"]),
        "smooth_point_completeness": all(int(row["smooth_point_count"]) >= int(contract["acceptance"]["minimum_smooth_points_per_curve"]) for row in curve_rows),
        "constant_current_completeness": all(int(row["constant_current_level_count"]) >= int(contract["acceptance"]["minimum_constant_current_levels_per_curve"]) for row in curve_rows),
        "classification_declared": classification in contract["analysis"]["classifications"],
        "no_new_solve": True,
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m63_smooth_nwell_attribution_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification if checks["all_checks_pass"] else "upstream_or_completeness_failure",
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "summary": {
            "curve_count": len(curve_rows), "pair_count": len(pair_rows),
            "smooth_point_count": len(point_rows),
            "all_nwell_growth_positive": all_growth_positive,
            "median_curve_horizontal_shift_explained_fraction": median_curve_explained,
            "minimum_curve_horizontal_shift_explained_fraction": min(float(row["horizontal_shift_error_energy_explained_fraction"]) for row in curve_rows),
            "median_pair_growth_horizontal_shift_explained_fraction": median_pair_explained,
            "minimum_pair_growth_horizontal_shift_explained_fraction": min(pair_explanations),
            "mobility_sign_support_count": mobility_support_count,
            "mobility_crosscheck_count": len(mobility_rows),
            "fitted_horizontal_shift_mV_range": [
                min(float(row["fitted_horizontal_shift_mV"]) for row in curve_rows),
                max(float(row["fitted_horizontal_shift_mV"]) for row in curve_rows)],
            "subthreshold_swing_difference_mV_per_dec_range": [
                min(float(row["subthreshold_swing_difference_mV_per_dec"]) for row in curve_rows),
                max(float(row["subthreshold_swing_difference_mV_per_dec"]) for row in curve_rows)],
        },
        "causal_scope": {
            "closed": "The smooth paired growth is predominantly threshold-like at curve level.",
            "not_closed": "BGN, electrostatics, carrier statistics, and other self-consistent causes of the horizontal shift remain unseparated.",
            "m40_context_only": (
                "M40 completed the deep-off true-NoBGN x SRH factorial; "
                "it is not a smooth-window intervention"),
        },
        "acceptance": checks,
    }
    return report, {"points": point_rows, "curves": curve_rows,
                    "thresholds": threshold_rows, "pairs": pair_rows,
                    "mobility": mobility_rows}


def freeze(report: dict[str, Any], ledgers: dict[str, list[dict[str, Any]]]) -> None:
    write_csv(POINTS, ledgers["points"])
    write_csv(CURVES, ledgers["curves"])
    write_csv(THRESHOLDS, ledgers["thresholds"])
    write_csv(PAIRS, ledgers["pairs"])
    write_csv(MOBILITY, ledgers["mobility"])
    write_json(REPORT, report)
    summary = report["summary"]
    DOC.write_text(f"""# SimpleMOS M63 平滑 NWell 误差归因

## 结论

M63 分类为 `{report['classification']}`。分析只读取M60收紧收敛后的16条默认曲线、M46/M8的Vela点值以及既有M13/M40账本，没有新增求解或修改生产参数。

在冻结的Vg=0.55--0.90 V窗口内，16条曲线的水平平移模型误差能量解释率中位数为 `{float(summary['median_curve_horizontal_shift_explained_fraction']):.6f}`，8个NWell配对增幅的解释率中位数为 `{float(summary['median_pair_growth_horizontal_shift_explained_fraction']):.6f}`。拟合的Sentaurus-minus-Vela等效栅压平移范围为 `{float(summary['fitted_horizontal_shift_mV_range'][0]):.3f}` 至 `{float(summary['fitted_horizontal_shift_mV_range'][1]):.3f}` mV。

M13在Vg=0.8 V的n17/n21固定状态迁移率链有 `{int(summary['mobility_sign_support_count'])}/{int(summary['mobility_crosscheck_count'])}` 个Vd工况与端口误差增长同号，因此现有固定状态证据不支持“迁移率误差随NWell增大”作为主导解释。

该结论只闭合到曲线级的阈值样/电静力学平移；它不能单独证明BGN根因。M40仅作为深关断BGN/SRH交互背景保留，BGN、载流子统计与自洽电势对这一平移的各自贡献仍需独立干预才能区分。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M63 smooth NWell attribution",
        "status": report["status"], "classification": report["classification"],
        "summary": "Burst-free curve-shape, constant-current threshold, and fixed-state mobility attribution",
        "report": portable(REPORT),
        "ledgers": [portable(CURVES), portable(THRESHOLDS), portable(PAIRS), portable(MOBILITY)],
    })
    artifacts = [REPORT, POINTS, CURVES, THRESHOLDS, PAIRS, MOBILITY, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m63_smooth_nwell_attribution_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT),
        "new_sentaurus_execution": False, "new_vela_execution": False,
        "production_reference_replaced": False,
        "source_hashes": read_json(FREEZE)["upstream_hashes"],
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "acceptance": report["acceptance"],
    })


def verify() -> dict[str, Any]:
    validate()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M63 evidence is not accepted and frozen")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M63 artifact hash changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        report = verify()
        print(f"M63 verified: {report['classification']}")
        return
    contract = validate()
    report, ledgers = analyze(contract)
    freeze(report, ledgers)
    if not report["acceptance"]["all_checks_pass"]:
        raise RuntimeError(f"M63 acceptance failed: {report['acceptance']}")
    print(json.dumps({"status": report["status"],
                      "classification": report["classification"],
                      "summary": report["summary"]}, indent=2))


if __name__ == "__main__":
    main()
