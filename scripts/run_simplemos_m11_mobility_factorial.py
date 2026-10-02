#!/usr/bin/env python3
"""Run the complete SimpleMOS PhuMob x Enormal x HFS factorial analysis."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import statistics
import sys
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m8a_confirmation as m8c  # noqa: E402
import run_simplemos_m8a_model_ablation as m8a  # noqa: E402
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402


CONTRACT = (REPO / "reference_tcad/simplemos_sentaurus2022"
            / "simplemos_m11_mobility_factorial_contract_v1.json")
OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
          / "m11_mobility_factorial")
M10_OUTPUT = (REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
              / "m10_fixed_state_replay")
RUNNER = REPO / "build-release/vela_example_runner.exe"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/simplemos_m11_factorial_20260828_v1"
FACTORS = ("phumob", "enormal", "hfs")
TERMS = ("phumob", "enormal", "hfs", "phumob:enormal",
         "phumob:hfs", "enormal:hfs", "phumob:enormal:hfs")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8", newline="\n")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]),
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def validate_contract(contract: dict[str, Any]) -> None:
    if contract.get("schema") != "vela.simplemos.sdevice.m11_mobility_factorial.v1":
        raise ValueError("unexpected M11 contract schema")
    if tuple(contract["factors"]) != FACTORS or len(contract["variants"]) != 8:
        raise ValueError("M11 requires a complete three-factor design")
    combinations = {tuple(bool(item[factor]) for factor in FACTORS)
                    for item in contract["variants"]}
    if len(combinations) != 8:
        raise ValueError("M11 variants do not cover every 2^3 combination")
    if [item["id"] for item in contract["devices"]] != ["n17", "n21"]:
        raise ValueError("M11 devices must be n17 and n21")
    lattice = [float(value) for value in contract["bias_matrix"]
               ["gate_lattice"]["values_V"]]
    if len(lattice) != 51 or any(not math.isclose(
            value, index * 0.05, rel_tol=0.0, abs_tol=1e-12)
            for index, value in enumerate(lattice)):
        raise ValueError("M11 gate lattice must be exact 0:0.05:2.5 V")
    if contract["comparison"]["interpolation"] != "forbidden":
        raise ValueError("M11 forbids interpolation")


def coded(variant: dict[str, Any], term: str) -> int:
    value = 1
    for factor in term.split(":"):
        value *= 1 if variant[factor] else -1
    return value


def factorial_effects(rows: list[dict[str, Any]], response: str,
                      variants: dict[str, dict[str, Any]]) -> dict[str, float]:
    if {row["variant"] for row in rows} != set(variants):
        raise ValueError("factorial effect requires all eight variants")
    result = {}
    for term in TERMS:
        positive = [float(row[response]) for row in rows
                    if coded(variants[row["variant"]], term) > 0]
        negative = [float(row[response]) for row in rows
                    if coded(variants[row["variant"]], term) < 0]
        result[term] = statistics.fmean(positive) - statistics.fmean(negative)
    return result


def compare_factorial(contract: dict[str, Any], output: Path) -> dict[str, Any]:
    destination = output / "comparisons"
    destination.mkdir(parents=True, exist_ok=True)
    cases = []
    for device in contract["devices"]:
        for variant in contract["variants"]:
            for vd in contract["bias_matrix"]["drain_voltages_V"]:
                tag = m8a.voltage_tag(float(vd))
                case_id = f"{device['id']}_{variant['id']}_vd_{tag}"
                reference = output / "sentaurus_reference" / f"{case_id}_reference.csv"
                candidate = (output / "vela" / device["id"] / variant["id"]
                             / "workflow" / f"vd_{tag}" / "20_gate_sweep.csv")
                result = m8c.compare_m4.compare_case(reference, candidate, contract)
                rows = result.pop("rows")
                result["diagnostic_metrics"] = m8a.diagnostic_metrics(rows, contract)
                comparison_csv = destination / f"{case_id}_comparison.csv"
                m8c.compare_m4.write_case_csv(comparison_csv, rows)
                result.update({
                    "case": case_id, "device": device["id"],
                    "variant": variant["id"], "drain_voltage_V": float(vd),
                    "comparison_csv": m8a.portable_path(comparison_csv),
                    "comparison_csv_sha256": m8a.sha256(comparison_csv),
                })
                cases.append(result)
    for item in cases:
        baseline = next(case for case in cases
                        if case["device"] == item["device"]
                        and case["drain_voltage_V"] == item["drain_voltage_V"]
                        and case["variant"] == "p1_e1_h1")
        item["effect_vs_full"] = {
            "maximum_log10_ratio_change": (
                item["maximum_absolute_log10_ratio_above_floor"]
                - baseline["maximum_absolute_log10_ratio_above_floor"]),
            "median_log10_ratio_change": (
                item["median_absolute_log10_ratio_above_floor"]
                - baseline["median_absolute_log10_ratio_above_floor"]),
        }
    report = {
        "schema": "vela.simplemos.sdevice.m11_factorial_comparison.v1",
        "status": "complete", "interpolation": "forbidden",
        "full_variant": "p1_e1_h1", "cases": cases,
    }
    write_json(destination / "comparison_report.json", report)
    return report


def analyze_self_consistent(contract: dict[str, Any], output: Path) -> dict[str, Any]:
    comparison = compare_factorial(contract, output)
    variants = {item["id"]: item for item in contract["variants"]}
    responses = ("maximum_absolute_log10_ratio_above_floor",
                 "median_absolute_log10_ratio_above_floor")
    conditions = []
    for device in ("n17", "n21"):
        for vd in (0.05, 1.0):
            rows = [case for case in comparison["cases"]
                    if case["device"] == device
                    and math.isclose(float(case["drain_voltage_V"]), vd)]
            conditions.append({
                "device": device, "drain_voltage_V": vd,
                "effects": {response: factorial_effects(rows, response, variants)
                            for response in responses},
            })
    aggregate_rows = []
    for variant in variants:
        selected = [case for case in comparison["cases"]
                    if case["variant"] == variant]
        aggregate_rows.append({
            "variant": variant,
            **{response: statistics.fmean(float(case[response])
                                          for case in selected)
               for response in responses},
        })
    report = {
        "schema": "vela.simplemos.sdevice.m11_self_consistent_factorial.v1",
        "status": "complete",
        "curve_count": len(comparison["cases"]),
        "direct_bias_points_per_solver": 51 * len(comparison["cases"]),
        "effect_definition": contract["effect_definition"],
        "conditions": conditions,
        "aggregate_effects": {
            response: factorial_effects(aggregate_rows, response, variants)
            for response in responses},
        "variant_aggregate": aggregate_rows,
        "comparison_report": m8a.portable_path(
            output / "comparisons/comparison_report.json"),
        "interpretation": (
            "Paired self-consistent Sentaurus/Vela effects include both direct "
            "mobility-model response and nonlinear electrostatic/carrier-state feedback."),
    }
    write_json(output / "self_consistent_factorial_report.json", report)
    rows = []
    for condition in conditions:
        for response, effects in condition["effects"].items():
            for term, value in effects.items():
                rows.append({"device": condition["device"],
                             "drain_voltage_V": condition["drain_voltage_V"],
                             "response": response, "term": term,
                             "effect_dex": value})
    write_csv(output / "self_consistent_factorial_effects.csv", rows)
    return report


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return float("nan")
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def frozen_case(case: dict[str, Any], variants: list[dict[str, Any]],
                export_by_state: dict[str, dict[str, Any]], output: Path,
                runner: Path) -> list[dict[str, Any]]:
    state_id = case["state"]
    root = output / "frozen_vela" / state_id
    root.mkdir(parents=True, exist_ok=True)
    source_config = Path(case["probe_statuses"]["vela_drive_edge_mobility"]["config"])
    state_file = M10_OUTPUT / "replay" / state_id / "sentaurus_state_for_vela.csv"
    reference_rows = m10.read_csv(
        M10_OUTPUT / "replay" / state_id / "edge_replay_comparison.csv")
    active = {int(row["edge_id"]) for row in reference_rows
              if row["active_current_edge"].lower() == "true"}
    sent_mobility = {int(row["edge_id"]): float(row["sentaurus_eMobility_cm2_V_s"])
                     for row in reference_rows}
    export_dir = REPO / export_by_state[state_id]["export_dir"]
    drain_nodes = m10.contact_nodes(export_dir, "drain")
    sent_current = abs(float(case["sentaurus_terminal_current_A_per_um"]))
    results = []
    for variant in variants:
        variant_id = variant["id"]
        mobility_path = root / f"{variant_id}_edge_mobility.csv"
        mobility_config_path = root / f"{variant_id}_edge_mobility.json"
        mobility_config = m10.probe_config(
            source_config, state_file, mobility_path,
            float(case["drain_voltage_V"]), float(case["gate_voltage_V"]),
            "edge_mobility_probe")
        mobility_config["solver"]["mobility"] = m8a.mobility_config(variant)
        mobility_config["simplemos_m11"] = {
            "variant": variant_id, "frozen_state": state_id,
            "read_only_formula_evaluation": True}
        write_json(mobility_config_path, mobility_config)
        m10.execute_runner(mobility_config_path, runner)

        sg_path = root / f"{variant_id}_sg_edges.csv"
        sg_config_path = root / f"{variant_id}_sg_edges.json"
        sg_config = m10.probe_config(
            source_config, state_file, sg_path,
            float(case["drain_voltage_V"]), float(case["gate_voltage_V"]),
            "sg_edge_flux_probe")
        sg_config["solver"]["mobility"] = m8a.mobility_config(variant)
        sg_config["simplemos_m11"] = {
            "variant": variant_id, "frozen_state": state_id,
            "read_only_formula_evaluation": True}
        write_json(sg_config_path, sg_config)
        m10.execute_runner(sg_config_path, runner)

        mobility_rows = {int(row["edge_id"]): row
                         for row in m10.read_csv(mobility_path)}
        errors = []
        log_mobility = []
        for edge in active:
            candidate = float(mobility_rows[edge]
                              ["electron_final_mobility_m2_V_s"]) * 1.0e4
            reference = sent_mobility[edge]
            if candidate > 0.0 and reference > 0.0:
                errors.append(abs(math.log10(candidate / reference)))
                log_mobility.append(math.log10(candidate))
        terminal = m10.drain_cut_current(m10.read_csv(sg_path), drain_nodes)
        terminal_abs = abs(float(terminal["total_A_per_um"]))
        results.append({
            "state": state_id, "device": case["device"],
            "drain_voltage_V": case["drain_voltage_V"],
            "gate_voltage_V": case["gate_voltage_V"],
            "variant": variant_id,
            "active_edge_count": len(errors),
            "median_log10_electron_mobility_cm2_V_s": statistics.median(log_mobility),
            "mobility_active_median_error_dex": statistics.median(errors),
            "mobility_active_p95_error_dex": percentile(errors, 0.95),
            "terminal_current_A_per_um": terminal_abs,
            "log10_terminal_current_A_per_um": math.log10(max(terminal_abs, 1e-300)),
            "terminal_current_error_dex": abs(math.log10(
                max(terminal_abs, 1e-300) / max(sent_current, 1e-300))),
        })
    return results


def run_frozen(contract: dict[str, Any], output: Path, runner: Path,
               jobs: int) -> dict[str, Any]:
    m10_report = read_json(M10_OUTPUT / "fixed_state_replay_report.json")
    exports = read_json(M10_OUTPUT / "sentaurus_export_manifest.json")
    export_by_state = {item["state"]: item for item in exports["states"]}
    variants = contract["variants"]
    with ThreadPoolExecutor(max_workers=min(jobs, len(m10_report["cases"]))) as pool:
        nested = list(pool.map(
            lambda case: frozen_case(case, variants, export_by_state,
                                     output, runner),
            m10_report["cases"]))
    rows = [row for group in nested for row in group]
    write_csv(output / "frozen_vela/frozen_factorial_metrics.csv", rows)
    variant_map = {item["id"]: item for item in variants}
    responses = ("median_log10_electron_mobility_cm2_V_s",
                 "mobility_active_median_error_dex",
                 "mobility_active_p95_error_dex",
                 "log10_terminal_current_A_per_um",
                 "terminal_current_error_dex")
    states = []
    for state_id in [case["state"] for case in m10_report["cases"]]:
        selected = [row for row in rows if row["state"] == state_id]
        states.append({
            "state": state_id,
            "device": selected[0]["device"],
            "drain_voltage_V": selected[0]["drain_voltage_V"],
            "gate_voltage_V": selected[0]["gate_voltage_V"],
            "effects": {response: factorial_effects(selected, response, variant_map)
                        for response in responses},
        })
    aggregate_rows = []
    for variant in variant_map:
        selected = [row for row in rows if row["variant"] == variant]
        aggregate_rows.append({
            "variant": variant,
            **{response: statistics.fmean(float(row[response])
                                          for row in selected)
               for response in responses},
        })
    report = {
        "schema": "vela.simplemos.sdevice.m11_frozen_vela_factorial.v1",
        "status": "complete", "state_count": len(states),
        "variant_evaluations": len(rows),
        "effect_definition": contract["effect_definition"],
        "sentaurus_variant_re_evaluation": False,
        "states": states,
        "aggregate_effects": {
            response: factorial_effects(aggregate_rows, response, variant_map)
            for response in responses},
        "variant_aggregate": aggregate_rows,
        "interpretation": (
            "All eight Vela mobility compositions are evaluated on each same "
            "full-physics Sentaurus state. Effects isolate Vela formula response "
            "but are not frozen Sentaurus model re-evaluations."),
    }
    write_json(output / "frozen_vela/frozen_factorial_report.json", report)
    effect_rows = []
    for state in states:
        for response, effects in state["effects"].items():
            for term, value in effects.items():
                effect_rows.append({
                    "state": state["state"], "device": state["device"],
                    "drain_voltage_V": state["drain_voltage_V"],
                    "gate_voltage_V": state["gate_voltage_V"],
                    "response": response, "term": term, "effect": value})
    write_csv(output / "frozen_vela/frozen_factorial_effects.csv", effect_rows)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-sentaurus", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--analyze-self-consistent", action="store_true")
    parser.add_argument("--run-frozen-vela", action="store_true")
    parser.add_argument("--sentaurus-jobs", type=int, default=4)
    parser.add_argument("--vela-jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=m8a.executable("ssh"))
    parser.add_argument("--scp-bin", default=m8a.executable("scp"))
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    if not any((args.prepare, args.live_sentaurus, args.extract_sentaurus,
                args.run_vela, args.analyze_self_consistent,
                args.run_frozen_vela)):
        parser.error("select at least one action")
    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    output = args.output_dir.resolve()
    tdrs = {key: value.resolve() for key, value in m8c.DEFAULT_TDRS.items()}
    m8c.validate_tdrs(contract, tdrs)
    manifest_path = output / "sentaurus_matrix_manifest.json"
    if args.prepare or not manifest_path.is_file():
        manifest = m8c.prepare_sentaurus(contract, contract_path, tdrs, output)
        m8c.prepare_vela(contract, tdrs, output, m8a.DEFAULT_MATERIALS,
                         IMPORTER)
    else:
        manifest = read_json(manifest_path)
    vela_report = (m8c.execute_vela(output, RUNNER, args.vela_jobs)
                   if args.run_vela else None)
    if args.live_sentaurus:
        banner = m8c.run_sentaurus(
            manifest, output, args.ssh_target, args.ssh_bin, args.scp_bin,
            args.remote_root, args.sentaurus_jobs)
        m8c.extract_references(contract, manifest, output, banner)
    elif args.extract_sentaurus:
        banner = (output / "sentaurus_banner.txt").read_text(
            encoding="utf-8").strip()
        m8c.extract_references(contract, manifest, output, banner)
    self_report = (analyze_self_consistent(contract, output)
                   if args.analyze_self_consistent else None)
    frozen_report = (run_frozen(contract, output, RUNNER, args.vela_jobs)
                     if args.run_frozen_vela else None)
    print(json.dumps({
        "status": "complete",
        "sentaurus_cases": len(manifest["cases"]),
        "vela_status": (vela_report or {}).get("status"),
        "self_consistent": (self_report or {}).get("status"),
        "frozen_vela": (frozen_report or {}).get("status"),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
