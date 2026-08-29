#!/usr/bin/env python3
"""Analyze SimpleMOS M12 frozen-state terminal-current sensitivity.

This is a read-only post-processing stage.  It reuses M10 edge probes and M8
comparison/SRH evidence; it does not run a nonlinear solve or change a Vela
physics default.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
Q = 1.602176634e-19
DEFAULT_M10_ROOT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                    / "m10_fixed_state_replay")
DEFAULT_COMPARISONS = (REPO / "reference_tcad/simplemos_sentaurus2022"
                       / "original_physics/comparisons")
DEFAULT_M8_DEEP_OFF = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                       / "m8_deep_off_diagnostics/vela")
DEFAULT_OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
                  / "m12_terminal_sensitivity")
DEFAULT_PORTABLE = (REPO / "reference_tcad/simplemos_sentaurus2022"
                    / "terminal_sensitivity")
STATE_RE = re.compile(
    r"^(n\d+)_vd_(0p05|1)_vg_(0|0p05|0p8|2p5)$")
CURVE_RE = re.compile(r"^(n\d+)_vd_(0p05|1)_comparison\.csv$")
REGIMES = {
    "deep_off": (0.0, 0.1),
    "subthreshold": (0.15, 0.6),
    "weak_inversion": (0.65, 0.95),
    "on_state": (1.0, 2.45),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def finite(value: float) -> float | None:
    return value if math.isfinite(value) else None


def weighted_mean(pairs: Iterable[tuple[float, float]]) -> float:
    numerator = 0.0
    denominator = 0.0
    for value, weight in pairs:
        if math.isfinite(value) and math.isfinite(weight) and weight > 0.0:
            numerator += value * weight
            denominator += weight
    return numerator / denominator if denominator > 0.0 else math.nan


def linear_fit(points: list[tuple[float, float]]) -> dict[str, float | int | None]:
    if len(points) < 2:
        return {"count": len(points), "slope_dex_per_V": None,
                "intercept_dex": None, "r_squared": None,
                "mean_signed_error_dex": None,
                "maximum_abs_error_dex": None}
    mean_x = sum(x for x, _ in points) / len(points)
    mean_y = sum(y for _, y in points) / len(points)
    xx = sum((x - mean_x) ** 2 for x, _ in points)
    xy = sum((x - mean_x) * (y - mean_y) for x, y in points)
    slope = xy / xx if xx > 0.0 else 0.0
    intercept = mean_y - slope * mean_x
    residual = sum((y - (intercept + slope * x)) ** 2 for x, y in points)
    total = sum((y - mean_y) ** 2 for _, y in points)
    r2 = 1.0 - residual / total if total > 0.0 else 1.0
    return {
        "count": len(points),
        "slope_dex_per_V": slope,
        "intercept_dex": intercept,
        "r_squared": r2,
        "mean_signed_error_dex": mean_y,
        "maximum_abs_error_dex": max(abs(y) for _, y in points),
    }


def error_spectrum(comparisons: Path) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    fits: list[dict[str, Any]] = []
    points_out: list[dict[str, Any]] = []
    curves = sorted(comparisons.glob("n*_vd_*_comparison.csv"))
    for path in curves:
        match = CURVE_RE.match(path.name)
        if not match:
            continue
        device, drain_tag = match.groups()
        drain = 0.05 if drain_tag == "0p05" else 1.0
        nwell = "1e17" if int(device[1:]) <= 20 else "2e17"
        points: list[tuple[float, float]] = []
        for row in read_csv(path):
            sentaurus = abs(float(row["sentaurus_current_A_per_um"]))
            vela = abs(float(row["vela_current_A_per_um"]))
            if sentaurus <= 0.0 or vela <= 0.0:
                continue
            gate = float(row["gate_voltage_V"])
            signed = math.log10(vela / sentaurus)
            points.append((gate, signed))
            points_out.append({
                "device": device,
                "NWell_cm_minus3": nwell,
                "drain_voltage_V": drain,
                "gate_voltage_V": gate,
                "signed_log10_ratio_dex": signed,
            })
        for regime, (lower, upper) in REGIMES.items():
            selected = [(x, y) for x, y in points
                        if lower - 1e-12 <= x <= upper + 1e-12]
            fit = linear_fit(selected)
            fits.append({
                "device": device,
                "NWell_cm_minus3": nwell,
                "drain_voltage_V": drain,
                "regime": regime,
                "gate_min_V": lower,
                "gate_max_V": upper,
                **fit,
            })
    return fits, points_out


def contact_nodes(path: Path) -> dict[str, set[int]]:
    result: dict[str, set[int]] = {}
    for row in read_csv(path):
        result[row["name"]] = {
            int(item) for item in row["node_ids"].split(";") if item}
    return result


def edge_contact(row: dict[str, str], contacts: dict[str, set[int]]) -> tuple[str, float] | None:
    n0, n1 = int(row["node0"]), int(row["node1"])
    matches: list[tuple[str, float]] = []
    for name in ("source", "drain", "substrate"):
        nodes = contacts[name]
        at0, at1 = n0 in nodes, n1 in nodes
        if at0 != at1:
            matches.append((name, 1.0 if at0 else -1.0))
    if len(matches) > 1:
        raise ValueError(f"edge {row['edge_id']} crosses multiple contact cuts")
    return matches[0] if matches else None


def condition_number(values: list[float]) -> float | None:
    denominator = abs(sum(values))
    numerator = sum(abs(value) for value in values)
    if denominator <= 1e-300:
        return None
    return finite(numerator / denominator)


def contact_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    base_e = [float(row["baseline_electron_A_per_um"]) for row in rows]
    base_h = [float(row["baseline_hole_A_per_um"]) for row in rows]
    alt_e = [float(row["sentaurus_drive_electron_A_per_um"]) for row in rows]
    alt_h = [float(row["sentaurus_drive_hole_A_per_um"]) for row in rows]
    base_total = sum(base_e) + sum(base_h)
    alt_total = sum(alt_e) + sum(alt_h)
    before = weighted_mean(
        (float(row["mobility_vela_drive_abs_error_dex"]),
         abs(float(row["sentaurus_projected_eLineCurrent_A_m"])))
        for row in rows)
    after = weighted_mean(
        (float(row["mobility_sentaurus_drive_abs_error_dex"]),
         abs(float(row["sentaurus_projected_eLineCurrent_A_m"])))
        for row in rows)
    direct_change = weighted_mean(
        (abs(float(row["electron_log10_mobility_scale"])),
         abs(float(row["baseline_electron_A_per_um"])))
        for row in rows)
    return {
        "edge_count": len(rows),
        "baseline_electron_A_per_um": sum(base_e),
        "baseline_hole_A_per_um": sum(base_h),
        "baseline_total_A_per_um": base_total,
        "sentaurus_drive_electron_A_per_um": sum(alt_e),
        "sentaurus_drive_hole_A_per_um": sum(alt_h),
        "sentaurus_drive_total_A_per_um": alt_total,
        "delta_electron_A_per_um": sum(alt_e) - sum(base_e),
        "delta_hole_A_per_um": sum(alt_h) - sum(base_h),
        "delta_total_A_per_um": alt_total - base_total,
        "kappa_electron": condition_number(base_e),
        "kappa_hole": condition_number(base_h),
        "kappa_total": condition_number(base_e + base_h),
        "current_weighted_abs_log10_mobility_change_dex": finite(direct_change),
        "projected_current_weighted_mobility_error_before_dex": finite(before),
        "projected_current_weighted_mobility_error_after_dex": finite(after),
        "projected_current_weighted_mobility_error_reduction_dex": finite(before - after),
    }


def domain_weighted_reduction(rows: Iterable[dict[str, str]]) -> dict[str, float | int | None]:
    selected = list(rows)
    before = weighted_mean(
        (float(row["mobility_vela_drive_abs_error_dex"]),
         abs(float(row["sentaurus_projected_eLineCurrent_A_m"])))
        for row in selected)
    after = weighted_mean(
        (float(row["mobility_sentaurus_drive_abs_error_dex"]),
         abs(float(row["sentaurus_projected_eLineCurrent_A_m"])))
        for row in selected)
    return {"edge_count": len(selected), "before_dex": finite(before),
            "after_dex": finite(after), "reduction_dex": finite(before - after)}


def analyze_state(case: dict[str, Any], m10_root: Path, output: Path,
                  threshold: float, identity_relative: float,
                  identity_absolute: float) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    state = case["state"]
    match = STATE_RE.match(state)
    if not match:
        raise ValueError(f"unexpected M10 state name: {state}")
    run_dir = m10_root / "replay" / state
    export_dir = m10_root / "sentaurus_exports" / state
    sg_rows = read_csv(run_dir / "vela_drive_sg_edges.csv")
    comparison_rows = read_csv(run_dir / "edge_replay_comparison.csv")
    comparison_by_edge = {int(row["edge_id"]): row for row in comparison_rows}
    primary = {int(row["edge_id"]): row for row in read_csv(
        run_dir / "sentaurus_mean_endpoint_magnitude_edge_mobility.csv")}
    nodes = contact_nodes(export_dir / "contacts.csv")

    contribution_rows: list[dict[str, Any]] = []
    contact_edges: set[int] = set()
    for row in sg_rows:
        cut = edge_contact(row, nodes)
        if cut is None:
            continue
        contact, sign = cut
        edge = int(row["edge_id"])
        if edge not in comparison_by_edge or edge not in primary:
            continue
        contact_edges.add(edge)
        spatial = comparison_by_edge[edge]
        base_mu_e = float(row["electron_mobility_m2_V_s"])
        base_mu_h = float(row["hole_mobility_m2_V_s"])
        alt_mu_e = float(primary[edge]["electron_final_mobility_m2_V_s"])
        alt_mu_h = float(primary[edge]["hole_final_mobility_m2_V_s"])
        scale_e = alt_mu_e / base_mu_e if base_mu_e > 0.0 else 1.0
        scale_h = alt_mu_h / base_mu_h if base_mu_h > 0.0 else 1.0
        base_e = (-Q * sign * float(row["electron_particle_line_flux_per_m_s"])
                  * 1.0e-6)
        base_h = (Q * sign * float(row["hole_particle_line_flux_per_m_s"])
                  * 1.0e-6)
        alt_e, alt_h = base_e * scale_e, base_h * scale_h
        contribution_rows.append({
            "state": state,
            "contact": contact,
            "edge_id": edge,
            "node0": int(row["node0"]),
            "node1": int(row["node1"]),
            "orientation_sign": sign,
            "x_mid_um": 0.5e6 * (float(row["x0"]) + float(row["x1"])),
            "y_mid_um": 0.5e6 * (float(row["y0"]) + float(row["y1"])),
            "baseline_electron_A_per_um": base_e,
            "baseline_hole_A_per_um": base_h,
            "sentaurus_drive_electron_A_per_um": alt_e,
            "sentaurus_drive_hole_A_per_um": alt_h,
            "delta_electron_A_per_um": alt_e - base_e,
            "delta_hole_A_per_um": alt_h - base_h,
            "baseline_electron_mobility_cm2_V_s": base_mu_e * 1.0e4,
            "sentaurus_drive_electron_mobility_cm2_V_s": alt_mu_e * 1.0e4,
            "electron_log10_mobility_scale": math.log10(scale_e) if scale_e > 0 else math.nan,
            "hole_log10_mobility_scale": math.log10(scale_h) if scale_h > 0 else math.nan,
            "sentaurus_projected_eLineCurrent_A_m": float(
                spatial["sentaurus_projected_eLineCurrent_A_m"]),
            "active_current_edge": spatial["active_current_edge"],
            "mobility_vela_drive_abs_error_dex": float(
                spatial["mobility_vela_drive_abs_error_dex"]),
            "mobility_sentaurus_drive_abs_error_dex": float(
                spatial["mobility_sentaurus_drive_abs_error_dex"]),
        })

    by_contact = {
        name: contact_metrics([row for row in contribution_rows
                               if row["contact"] == name])
        for name in ("source", "drain", "substrate")
    }
    replay = case["terminal_current_replay"]
    observed_delta = (
        float(replay["sentaurus_drive_vela_hfs"]["total_A_per_um"])
        - float(replay["vela_drive_vela_hfs"]["total_A_per_um"]))
    reconstructed = float(by_contact["drain"]["delta_total_A_per_um"])
    identity_abs = abs(reconstructed - observed_delta)
    identity_rel = identity_abs / max(abs(reconstructed), abs(observed_delta), 1e-300)
    identity_pass = identity_abs <= identity_absolute or identity_rel <= identity_relative

    active = [row for row in comparison_rows
              if row["active_current_edge"].lower() == "true"]
    internal_active = [row for row in active if int(row["edge_id"]) not in contact_edges]
    active_weighted = domain_weighted_reduction(active)
    internal_weighted = domain_weighted_reduction(internal_active)
    p95_before = float(case["electron_mobility"]["vela_drive"]
                       ["active_edges_abs_error_dex"]["p95"])
    p95_after = float(case["electron_mobility"]["sentaurus_drive"]
                      ["active_edges_abs_error_dex"]["p95"])
    p95_reduction = p95_before - p95_after
    decision_reductions = [
        abs(float(by_contact[name]
                  ["projected_current_weighted_mobility_error_reduction_dex"] or 0.0))
        for name in ("source", "drain")]
    cut_ratio = max(decision_reductions) / max(abs(p95_reduction), 1e-300)
    h1_pass = cut_ratio < threshold and identity_pass
    baseline_total = float(replay["vela_drive_vela_hfs"]["total_A_per_um"])
    alternate_total = float(replay["sentaurus_drive_vela_hfs"]["total_A_per_um"])
    replay_shift = (math.log10(abs(alternate_total) / abs(baseline_total))
                    if baseline_total != 0.0 and alternate_total != 0.0 else math.nan)
    frozen_contact_sum = sum(float(metrics["baseline_total_A_per_um"])
                             for metrics in by_contact.values())
    result = {
        "state": state,
        "device": case["device"],
        "drain_voltage_V": case["drain_voltage_V"],
        "gate_voltage_V": case["gate_voltage_V"],
        "active_edge_p95": {"before_dex": p95_before, "after_dex": p95_after,
                            "reduction_dex": p95_reduction},
        "active_edge_projected_current_weighted": active_weighted,
        "internal_active_edge_projected_current_weighted": internal_weighted,
        "contacts": by_contact,
        "drain_replay_identity": {
            "reconstructed_delta_A_per_um": reconstructed,
            "reported_delta_A_per_um": observed_delta,
            "absolute_error_A_per_um": identity_abs,
            "relative_error": identity_rel,
            "pass": identity_pass,
        },
        "terminal_current_drive_substitution_signed_shift_dex": replay_shift,
        "frozen_source_drain_substrate_contact_sum_A_per_um": frozen_contact_sum,
        "h1": {"contact_to_active_reduction_ratio": cut_ratio,
               "threshold": threshold, "pass": h1_pass},
    }
    write_csv(output / "edge_contributions" / f"{state}.csv", contribution_rows)
    return result, contribution_rows


def deep_off_context(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for device in ("n17", "n21"):
        for drain_tag, drain in (("0p05", 0.05), ("1", 1.0)):
            for gate_tag, gate in (("0", 0.0), ("0p05", 0.05)):
                path = root / f"{device}_vd_{drain_tag}" / f"vg_{gate_tag}" / "srh_balance.csv"
                if not path.exists():
                    continue
                source = read_csv(path)[-1]
                rows.append({
                    "state": f"{device}_vd_{drain_tag}_vg_{gate_tag}",
                    "device": device,
                    "drain_voltage_V": drain,
                    "gate_voltage_V": gate,
                    "srh_net_current_A_per_um": float(source["srh_net_current_A_per_um"]),
                    "srh_generation_current_A_per_um": float(
                        source["srh_generation_current_A_per_um"]),
                    "srh_recombination_current_A_per_um": float(
                        source["srh_recombination_current_A_per_um"]),
                    "four_terminal_kcl_residual_A_per_um": float(
                        source["four_terminal_kcl_residual_A_per_um"]),
                    "numerical_status": source["numerical_status"],
                    "scope": "M8 self-consistent Vela state; context only",
                })
    return rows


def flatten_state(case: dict[str, Any]) -> dict[str, Any]:
    source = case["contacts"]["source"]
    drain = case["contacts"]["drain"]
    substrate = case["contacts"]["substrate"]
    return {
        "state": case["state"],
        "device": case["device"],
        "drain_voltage_V": case["drain_voltage_V"],
        "gate_voltage_V": case["gate_voltage_V"],
        "active_p95_reduction_dex": case["active_edge_p95"]["reduction_dex"],
        "source_weighted_reduction_dex": source[
            "projected_current_weighted_mobility_error_reduction_dex"],
        "drain_weighted_reduction_dex": drain[
            "projected_current_weighted_mobility_error_reduction_dex"],
        "substrate_weighted_reduction_dex": substrate[
            "projected_current_weighted_mobility_error_reduction_dex"],
        "contact_to_active_reduction_ratio": case["h1"]
            ["contact_to_active_reduction_ratio"],
        "terminal_drive_shift_dex": case[
            "terminal_current_drive_substitution_signed_shift_dex"],
        "drain_identity_absolute_error_A_per_um": case[
            "drain_replay_identity"]["absolute_error_A_per_um"],
        "drain_identity_relative_error": case[
            "drain_replay_identity"]["relative_error"],
        "drain_kappa_electron": drain["kappa_electron"],
        "drain_kappa_hole": drain["kappa_hole"],
        "drain_kappa_total": drain["kappa_total"],
        "source_kappa_total": source["kappa_total"],
        "substrate_kappa_total": substrate["kappa_total"],
        "frozen_contact_sum_A_per_um": case[
            "frozen_source_drain_substrate_contact_sum_A_per_um"],
        "h1_pass": case["h1"]["pass"],
    }


def flatten_contacts(case: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, metrics in case["contacts"].items():
        rows.append({"state": case["state"], "device": case["device"],
                     "drain_voltage_V": case["drain_voltage_V"],
                     "gate_voltage_V": case["gate_voltage_V"],
                     "contact": name, **metrics})
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--m10-root", type=Path, default=DEFAULT_M10_ROOT)
    parser.add_argument("--comparisons", type=Path, default=DEFAULT_COMPARISONS)
    parser.add_argument("--m8-deep-off", type=Path, default=DEFAULT_M8_DEEP_OFF)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--portable-output-dir", type=Path, default=DEFAULT_PORTABLE)
    args = parser.parse_args()
    m10_root = args.m10_root.resolve()
    output = args.output_dir.resolve()
    portable_output = args.portable_output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    portable_output.mkdir(parents=True, exist_ok=True)

    contract = read_json(REPO / "reference_tcad/simplemos_sentaurus2022"
                         / "simplemos_m12_terminal_sensitivity_contract_v1.json")
    h1_contract = contract["h1"]
    m10_report_path = m10_root / "fixed_state_replay_report.json"
    m10_report = read_json(m10_report_path)
    spectrum_rows, spectrum_points = error_spectrum(args.comparisons.resolve())
    cases: list[dict[str, Any]] = []
    contacts_flat: list[dict[str, Any]] = []
    for source_case in m10_report["cases"]:
        case, _ = analyze_state(
            source_case, m10_root, output,
            float(h1_contract["contact_to_active_reduction_ratio_threshold"]),
            float(h1_contract["identity_relative_tolerance"]),
            float(h1_contract["identity_absolute_tolerance_A_per_um"]))
        cases.append(case)
        contacts_flat.extend(flatten_contacts(case))
    context = deep_off_context(args.m8_deep_off.resolve())

    state_rows = [flatten_state(case) for case in cases]
    artifacts = {
        "state_summary": portable_output / "m12_state_summary.csv",
        "contact_summary": portable_output / "m12_contact_summary.csv",
        "error_spectrum": portable_output / "m12_error_spectrum.csv",
        "error_spectrum_points": portable_output / "m12_error_spectrum_points.csv",
        "deep_off_context": portable_output / "m12_deep_off_context.csv",
    }
    write_csv(artifacts["state_summary"], state_rows)
    write_csv(artifacts["contact_summary"], contacts_flat)
    write_csv(artifacts["error_spectrum"], spectrum_rows)
    write_csv(artifacts["error_spectrum_points"], spectrum_points)
    write_csv(artifacts["deep_off_context"], context)

    identity_passes = sum(bool(case["drain_replay_identity"]["pass"])
                          for case in cases)
    h1_passes = sum(bool(case["h1"]["pass"]) for case in cases)
    deep_off_cases = [case for case in cases if case["gate_voltage_V"] <= 0.05]
    kappa_maxima = {
        name: {
            carrier: max(float(case["contacts"][name][carrier] or 0.0)
                         for case in deep_off_cases)
            for carrier in ("kappa_electron", "kappa_hole", "kappa_total")}
        for name in ("source", "drain", "substrate")}
    high_nwell_low_drain_weak = [
        float(row["signed_log10_ratio_dex"])
        for row in spectrum_points
        if row["NWell_cm_minus3"] == "2e17"
        and float(row["drain_voltage_V"]) == 0.05
        and float(row["gate_voltage_V"]) <= 0.95]
    report = {
        "schema": "vela.simplemos.sdevice.m12_terminal_sensitivity_report.v1",
        "status": "complete",
        "scope": contract["scope"],
        "execution": {
            "state_count": len(cases),
            "residual_curve_count": len({(row["device"], row["drain_voltage_V"])
                                          for row in spectrum_points}),
            "residual_fit_count": len(spectrum_rows),
            "cpp_changed": False,
            "default_model_changed": False,
        },
        "findings": {
            "drain_replay_identity_pass_count": identity_passes,
            "h1_state_pass_count": h1_passes,
            "h1_all_states_pass": h1_passes == len(cases),
            "maximum_contact_to_active_reduction_ratio": max(
                float(case["h1"]["contact_to_active_reduction_ratio"])
                for case in cases),
            "maximum_abs_terminal_drive_shift_dex": max(
                abs(float(case["terminal_current_drive_substitution_signed_shift_dex"]))
                for case in cases),
            "deep_off_contact_kappa_maxima": kappa_maxima,
            "deep_off_any_component_kappa_ge_100_count": sum(
                1 for case in deep_off_cases
                for name in ("source", "drain", "substrate")
                for carrier in ("kappa_electron", "kappa_hole", "kappa_total")
                if float(case["contacts"][name][carrier] or 0.0) >= 100.0),
            "high_nwell_low_drain_vg_le_0p95": {
                "point_count": len(high_nwell_low_drain_weak),
                "all_vela_above_sentaurus": all(
                    value > 0.0 for value in high_nwell_low_drain_weak),
                "mean_signed_error_dex": (
                    sum(high_nwell_low_drain_weak)
                    / len(high_nwell_low_drain_weak)),
                "maximum_signed_error_dex": max(high_nwell_low_drain_weak),
            },
        },
        "interpretation_guards": [
            "Contact-cut delta-current reconstruction is an exact frozen-state mobility scaling identity.",
            "Internal edges have no unique direct terminal-current contribution under a frozen state.",
            "The frozen source-drain-substrate contact sum is a mapped-state continuity residual, not an SRH closure.",
            "M8 self-consistent deep-off SRH rows are included only as separately labelled context.",
            "Large current condition number proves amplification potential, not complete attribution by itself."
        ],
        "cases": cases,
        "self_consistent_deep_off_context": context,
        "inputs": {
            "contract": portable(REPO / "reference_tcad/simplemos_sentaurus2022"
                                 / "simplemos_m12_terminal_sensitivity_contract_v1.json"),
            "m10_report": portable(m10_report_path),
            "m10_report_sha256": sha256(m10_report_path),
            "m8_comparison_report": portable(args.comparisons.resolve()
                / "comparison_report.json"),
            "m8_comparison_report_sha256": sha256(args.comparisons.resolve()
                / "comparison_report.json"),
        },
        "artifacts": {
            name: {"path": portable(path), "sha256": sha256(path)}
            for name, path in artifacts.items()
        },
    }
    report_path = portable_output / "m12_terminal_sensitivity_report.json"
    write_json(report_path, report)
    write_json(output / "m12_terminal_sensitivity_report.json", report)
    print(json.dumps({"status": "complete", "states": len(cases),
                      "h1_passes": h1_passes,
                      "identity_passes": identity_passes,
                      "report": portable(report_path)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
