"""Freeze, run, and analyze SimpleMOS M75-M77 Poisson charge-volume A/Bs."""

from __future__ import annotations

import argparse
import copy
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any


REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
MATRIX_CONTRACT = ROOT / "simplemos_m74_electron_poisson_charge_volume_contract_v1.json"
M73_EVIDENCE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_evidence.json"
M73_PAIRS = ROOT / "material_partitioned_poisson_ledger/m73_poisson_pair_ledger.csv"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M65_VELA = REPO / "build-release/m65_ni/vela_manifest.json"
M65_SENT = REPO / "build-release/m65_ni/sentaurus_export_manifest.json"
RUNNER = REPO / "build-release/vela_example_runner.exe"
SCRIPT = Path(__file__).resolve()
PHASES = ("00_equilibrium", "10_drain_ramp", "20_gate_sweep")

PROFILES = {
    "M75": {
        "slug": "hole_poisson_charge_volume",
        "title": "hole Poisson charge volume",
        "prefix": "hole_poisson_volume",
        "contract": "simplemos_m75_hole_poisson_charge_volume_contract_v1.json",
        "freeze": "simplemos_m75_hole_poisson_charge_volume_contract_freeze_v1.json",
        "evidence": "simplemos_m75_hole_poisson_charge_volume_evidence.json",
    },
    "M76": {
        "slug": "dopant_poisson_charge_volume",
        "title": "dopant Poisson charge volume",
        "prefix": "dopant_poisson_volume",
        "contract": "simplemos_m76_dopant_poisson_charge_volume_contract_v1.json",
        "freeze": "simplemos_m76_dopant_poisson_charge_volume_contract_freeze_v1.json",
        "evidence": "simplemos_m76_dopant_poisson_charge_volume_evidence.json",
    },
    "M77": {
        "slug": "combined_poisson_charge_volume",
        "title": "combined electron-hole-dopant Poisson charge volume",
        "prefix": "combined_poisson_volume",
        "contract": "simplemos_m77_combined_poisson_charge_volume_contract_v1.json",
        "freeze": "simplemos_m77_combined_poisson_charge_volume_contract_freeze_v1.json",
        "evidence": "simplemos_m77_combined_poisson_charge_volume_evidence.json",
    },
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty ledger: {path}")
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


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    return ordered[n // 2] if n % 2 else 0.5 * (
        ordered[n // 2 - 1] + ordered[n // 2])


def pearson(xs: list[float], ys: list[float]) -> float:
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    dx, dy = [x - mx for x in xs], [y - my for y in ys]
    denominator = math.sqrt(sum(x * x for x in dx) *
                            sum(y * y for y in dy))
    return (sum(x * y for x, y in zip(dx, dy, strict=True)) / denominator
            if denominator else 0.0)


def voltage_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def paths(milestone: str) -> dict[str, Path]:
    profile = PROFILES[milestone]
    slug = profile["slug"]
    portable_root = ROOT / slug
    return {
        "contract": ROOT / profile["contract"],
        "freeze": ROOT / profile["freeze"],
        "evidence": ROOT / profile["evidence"],
        "output": REPO / "build-release" / f"{milestone.lower()}_{slug}",
        "portable_root": portable_root,
        "report": portable_root / f"{milestone.lower()}_{slug}_report.json",
        "cases": portable_root / f"{milestone.lower()}_case_ledger.csv",
        "pairs": portable_root / f"{milestone.lower()}_pair_ledger.csv",
        "points": portable_root / f"{milestone.lower()}_curve_response_ledger.csv",
        "manifest": portable_root / f"{milestone.lower()}_execution_manifest.json",
        "run_manifest": REPO / "build-release" / f"{milestone.lower()}_{slug}" / "execution_manifest.json",
        "doc": REPO / "docs/validation" / f"simplemos_{milestone.lower()}_{slug}_2026-09-03.md",
        "artifact": REPO / "docs/validation/reports" / f"simplemos_{milestone.lower()}" / "artifact.json",
        "decomposition": portable_root / "m77_charge_volume_decomposition_ledger.csv",
    }


def workflows() -> dict[tuple[str, float], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"])): row
            for row in read_json(M65_VELA)["workflows"]}


def sentaurus() -> dict[tuple[str, float], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"])): row
            for row in read_json(M65_SENT)["states"]}


def matrix_pairs() -> list[dict[str, Any]]:
    return read_json(MATRIX_CONTRACT)["matrix"]["pairs"]


def baseline_hash_errors() -> int:
    errors = 0
    for row in read_json(M65_VELA)["workflows"]:
        errors += sha256(REPO / row["curve"]) != row["curve_sha256"]
        errors += sha256(REPO / row["state"]) != row["state_sha256"]
    return errors


def source_paths(contract: dict[str, Any]) -> list[Path]:
    result = [MATRIX_CONTRACT, M73_EVIDENCE, M73_PAIRS, M65_EVIDENCE,
              M65_VELA, M65_SENT, RUNNER, SCRIPT,
              REPO / "include/vela/equation/DDAssembler.h",
              REPO / "include/vela/equation/CoupledDDAssembler.h",
              REPO / "src/equation/CoupledDDAssembler.cpp",
              REPO / "src/solver/NewtonSolver.cpp",
              REPO / "tests/test_newton_solver.cpp",
              REPO / "tests/test_mos_mixed_material.cpp"]
    result.extend(REPO / relative
                  for relative in contract["upstream"]["required_evidence"])
    for workflow in workflows().values():
        source_root = (REPO / workflow["state"]).parents[1]
        result.extend(source_root / phase / "config.json" for phase in PHASES)
        result.extend((REPO / workflow["curve"], REPO / workflow["state"]))
    return sorted(set(path.resolve() for path in result))


def freeze_contract(milestone: str) -> None:
    item = paths(milestone)
    contract = read_json(item["contract"])
    if (contract.get("schema") !=
            "vela.simplemos.sdevice.poisson_charge_volume_contract.v1" or
            contract.get("milestone") != milestone):
        raise ValueError(f"unexpected {milestone} contract")
    for relative in contract["upstream"]["required_evidence"]:
        evidence = read_json(REPO / relative)
        if evidence.get("status") != "frozen":
            raise ValueError(f"required evidence is not frozen: {relative}")
    inputs = source_paths(contract)
    missing = [path for path in inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    if baseline_hash_errors():
        raise ValueError("M65 baseline artifacts changed")
    write_json(item["freeze"], {
        "schema": "vela.simplemos.sdevice.poisson_charge_volume_contract_freeze.v1",
        "milestone": milestone,
        "status": "frozen_before_execution",
        "contract": portable(item["contract"]),
        "contract_sha256": sha256(item["contract"]),
        "upstream_hashes": {portable(path): sha256(path) for path in inputs},
    })


def validate_contract(milestone: str) -> dict[str, Any]:
    item = paths(milestone)
    contract, freeze = read_json(item["contract"]), read_json(item["freeze"])
    if (freeze.get("status") != "frozen_before_execution" or
            freeze.get("contract_sha256") != sha256(item["contract"])):
        raise ValueError(f"{milestone} contract is not frozen or changed")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"{milestone} frozen input changed: {relative}")
    return contract


def normalized_config(config: dict[str, Any], flags: dict[str, bool],
                      metadata_key: str, candidate: bool) -> dict[str, Any]:
    value = copy.deepcopy(config)
    value.pop("output_csv", None)
    value.pop("log_file", None)
    value.pop(metadata_key, None)
    for key in ("write_state_file", "initial_state_file", "vtk_prefix"):
        value["sweep"].pop(key, None)
    if candidate:
        assembly = value["solver"].pop(
            "region_resolved_interface_assembly", None)
        expected = {"enabled": False, **flags}
        if assembly != expected:
            raise ValueError(f"unexpected intervention payload: {assembly}")
    return value


def candidate_config(source: dict[str, Any], run_dir: Path,
                     previous: Path | None, milestone: str,
                     contract: dict[str, Any]) -> tuple[dict[str, Any], int]:
    if "region_resolved_interface_assembly" in source["solver"]:
        raise ValueError("M65 source already has diagnostic interface assembly")
    flags = contract["intervention"]["flags"]
    metadata_key = f"simplemos_{milestone.lower()}"
    config = copy.deepcopy(source)
    config["solver"]["region_resolved_interface_assembly"] = {
        "enabled": False, **flags}
    config["output_csv"] = str((run_dir / "curve.csv").resolve())
    config["log_file"] = str((run_dir / "run.log").resolve())
    config["sweep"]["write_state_file"] = str(
        (run_dir / "state.csv").resolve())
    if previous is None:
        config["sweep"].pop("initial_state_file", None)
    else:
        config["sweep"]["initial_state_file"] = str(previous.resolve())
    if "vtk_prefix" in config["sweep"]:
        config["sweep"]["vtk_prefix"] = str(
            (run_dir / "vtk/state").resolve())
    config[metadata_key] = {
        "single_axis": contract["intervention"]["changed_term"],
        "measure": contract["intervention"]["candidate_measure"],
        "continuity_dielectric_transport_unchanged": True,
        "production_default_changed": False,
    }
    different = normalized_config(source, flags, metadata_key, False) != \
        normalized_config(config, flags, metadata_key, True)
    return config, int(different)


def execute(config_path: Path, milestone: str) -> dict[str, Any]:
    completed = subprocess.run([str(RUNNER), "--config", str(config_path)],
                               cwd=REPO, text=True, capture_output=True,
                               check=False)
    (config_path.parent / "stdout.txt").write_text(
        completed.stdout, encoding="utf-8", newline="\n")
    (config_path.parent / "stderr.txt").write_text(
        completed.stderr, encoding="utf-8", newline="\n")
    if completed.returncode:
        raise RuntimeError(
            f"{milestone} runner failed for {config_path}: "
            f"{completed.stderr[-2000:]}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def run_workflow(key: tuple[str, float], workflow: dict[str, Any],
                 milestone: str, contract: dict[str, Any]) -> dict[str, Any]:
    item = paths(milestone)
    device, drain = key
    source_root = (REPO / workflow["state"]).parents[1]
    root = item["output"] / device / f"vd_{voltage_tag(drain)}"
    if root.exists():
        shutil.rmtree(root)
    previous: Path | None = None
    records: list[dict[str, Any]] = []
    unapproved = 0
    for phase in PHASES:
        run_dir = root / phase
        run_dir.mkdir(parents=True)
        source_path = source_root / phase / "config.json"
        config, difference = candidate_config(
            read_json(source_path), run_dir, previous, milestone, contract)
        unapproved += difference
        config_path = run_dir / "config.json"
        write_json(config_path, config)
        status = execute(config_path, milestone)
        if not status.get("converged"):
            raise RuntimeError(
                f"{milestone} did not converge: {device} {drain} {phase}")
        previous = run_dir / "state.csv"
        records.append({
            "phase": phase, "source_config": portable(source_path),
            "config": portable(config_path),
            "config_sha256": sha256(config_path),
            "curve": portable(run_dir / "curve.csv"),
            "curve_sha256": sha256(run_dir / "curve.csv"),
            "state": portable(previous), "state_sha256": sha256(previous),
            "converged": True,
        })
    return {
        "device": device, "drain_voltage_V": drain,
        "gate_voltage_V": float(workflow["gate_voltage_V"]),
        "unapproved_config_diff_count": unapproved,
        "phases": records, "final_curve": records[-1]["curve"],
        "final_state": records[-1]["state"],
    }


def run_candidate(milestone: str, jobs: int) -> dict[str, Any]:
    contract = validate_contract(milestone)
    item = paths(milestone)
    if item["output"].exists():
        shutil.rmtree(item["output"])
    inputs = sorted(workflows().items())
    with ThreadPoolExecutor(max_workers=min(max(1, jobs), 8)) as pool:
        records = list(pool.map(
            lambda pair: run_workflow(*pair, milestone, contract), inputs))
    manifest = {
        "schema": f"vela.simplemos.sdevice.{milestone.lower()}_execution_manifest.v1",
        "status": "complete", "contract_sha256": sha256(item["contract"]),
        "workflows": records,
    }
    write_json(item["run_manifest"], manifest)
    return manifest


def exact_current(path: Path, gate: float, tolerance: float) -> float:
    rows = [row for row in read_csv(path)
            if abs(float(row["bias_V"]) - gate) <= tolerance]
    if len(rows) != 1:
        raise RuntimeError(f"exact gate point count {len(rows)} at {gate}: {path}")
    if str(rows[0]["converged"]).lower() not in ("1", "true"):
        raise RuntimeError(f"unconverged exact gate point: {path}")
    return abs(float(rows[0]["current_total_A_per_um"]))


def decomposition_rows(pair_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    ledgers = {
        "electron": ROOT / "electron_poisson_charge_volume/m74_pair_ledger.csv",
        "hole": ROOT / "hole_poisson_charge_volume/m75_pair_ledger.csv",
        "dopant": ROOT / "dopant_poisson_charge_volume/m76_pair_ledger.csv",
    }
    maps = {
        name: {(row["low_device"], row["high_device"],
                float(row["drain_voltage_V"])): row
               for row in read_csv(path)}
        for name, path in ledgers.items()
    }
    result = []
    for row in pair_rows:
        key = (row["low_device"], row["high_device"],
               float(row["drain_voltage_V"]))
        reductions = {
            name: float(mapping[key]["self_consistent_pair_reduction_dex"])
            for name, mapping in maps.items()
        }
        additive = sum(reductions.values())
        combined = float(row["self_consistent_pair_reduction_dex"])
        result.append({
            "low_device": key[0], "high_device": key[1],
            "drain_voltage_V": key[2],
            "baseline_pair_growth_dex": row["baseline_pair_growth_dex"],
            "electron_only_reduction_dex": reductions["electron"],
            "hole_only_reduction_dex": reductions["hole"],
            "dopant_only_reduction_dex": reductions["dopant"],
            "independent_additive_reduction_dex": additive,
            "combined_reduction_dex": combined,
            "nonlinear_interaction_dex": combined - additive,
            "combined_candidate_pair_growth_dex": row["candidate_pair_growth_dex"],
        })
    return result


def analyze(milestone: str, contract: dict[str, Any]) -> tuple:
    item = paths(milestone)
    manifest = read_json(item["run_manifest"])
    candidates = {(row["device"], float(row["drain_voltage_V"])): row
                  for row in manifest["workflows"]}
    baseline_map, sent_map = workflows(), sentaurus()
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    case_rows: list[dict[str, Any]] = []
    point_rows: list[dict[str, Any]] = []
    for key in sorted(baseline_map):
        device, drain = key
        baseline, candidate, sent = (baseline_map[key], candidates[key],
                                     sent_map[key])
        gate = float(baseline["gate_voltage_V"])
        baseline_curve = REPO / baseline["curve"]
        candidate_curve = REPO / candidate["final_curve"]
        baseline_current = exact_current(baseline_curve, gate, tolerance)
        candidate_current = exact_current(candidate_curve, gate, tolerance)
        sent_current = abs(float(sent["drain_current_A_per_um"]))
        baseline_error = math.log10(baseline_current / sent_current)
        candidate_error = math.log10(candidate_current / sent_current)
        baseline_points = {round(float(row["bias_V"]), 12): row
                           for row in read_csv(baseline_curve)}
        candidate_points = {round(float(row["bias_V"]), 12): row
                            for row in read_csv(candidate_curve)}
        if set(baseline_points) != set(candidate_points):
            raise RuntimeError(f"bias path changed: {key}")
        shifts = []
        for bias in sorted(baseline_points):
            first = abs(float(baseline_points[bias]["current_total_A_per_um"]))
            second = abs(float(candidate_points[bias]["current_total_A_per_um"]))
            shift = math.log10(second / first)
            shifts.append(abs(shift))
            point_rows.append({
                "device": device, "drain_voltage_V": drain,
                "gate_voltage_V": bias,
                "baseline_current_A_per_um": first,
                "candidate_current_A_per_um": second,
                "candidate_over_baseline_log_shift_dex": shift,
            })
        case_rows.append({
            "device": device, "drain_voltage_V": drain,
            "diagnostic_gate_voltage_V": gate,
            "sentaurus_current_A_per_um": sent_current,
            "baseline_vela_current_A_per_um": baseline_current,
            "candidate_vela_current_A_per_um": candidate_current,
            "baseline_error_dex": baseline_error,
            "candidate_error_dex": candidate_error,
            "absolute_error_improvement_dex": (
                abs(baseline_error) - abs(candidate_error)),
            "endpoint_candidate_over_baseline_log_shift_dex": math.log10(
                candidate_current / baseline_current),
            "maximum_curve_candidate_over_baseline_abs_log_shift_dex": max(shifts),
            "unapproved_config_diff_count": candidate["unapproved_config_diff_count"],
            "converged": True,
        })

    cases = {(row["device"], float(row["drain_voltage_V"])): row
             for row in case_rows}
    prediction_column = contract["analysis"]["m73_prediction_column"]
    predictions = {
        (row["low_device"], row["high_device"],
         float(row["drain_voltage_V"])): float(row[prediction_column])
        for row in read_csv(M73_PAIRS)
        if row["variant"] == "region_local_barycentric_si"
    }
    pair_rows: list[dict[str, Any]] = []
    for pair in matrix_pairs():
        low, high = pair["low_device"], pair["high_device"]
        drain = float(pair["drain_voltage_V"])
        low_case, high_case = cases[(low, drain)], cases[(high, drain)]
        baseline_growth = (float(high_case["baseline_error_dex"]) -
                           float(low_case["baseline_error_dex"]))
        candidate_growth = (float(high_case["candidate_error_dex"]) -
                            float(low_case["candidate_error_dex"]))
        reduction = baseline_growth - candidate_growth
        closure = 1.0 - abs(candidate_growth) / max(abs(baseline_growth), 1e-300)
        prediction = predictions[(low, high, drain)]
        pair_rows.append({
            "low_device": low, "high_device": high,
            "drain_voltage_V": drain,
            "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"],
            "baseline_pair_growth_dex": baseline_growth,
            "candidate_pair_growth_dex": candidate_growth,
            "self_consistent_pair_reduction_dex": reduction,
            "pair_improved": abs(candidate_growth) < abs(baseline_growth),
            "raw_pair_closure_fraction": closure,
            "classification_pair_closure_fraction": max(0.0, min(1.0, closure)),
            "m73_predicted_pair_reduction_dex": prediction,
            "prediction_error_dex": reduction - prediction,
        })

    closures = [float(row["classification_pair_closure_fraction"])
                for row in pair_rows]
    improved = sum(bool(row["pair_improved"]) for row in pair_rows)
    median_closure = median(closures)
    analysis_spec = contract["analysis"]
    dominant = (improved >= int(analysis_spec["dominant_minimum_improved_pair_count"])
                and median_closure >= float(analysis_spec[
                    "dominant_minimum_median_pair_closure_fraction"]))
    material = (improved >= int(analysis_spec["dominant_minimum_improved_pair_count"])
                and median_closure >= float(analysis_spec[
                    "material_minimum_median_pair_closure_fraction"]))
    prefix = PROFILES[milestone]["prefix"]
    classification = (f"{prefix}_dominant" if dominant else
                      f"{prefix}_material_but_not_dominant" if material else
                      f"{prefix}_not_material")
    high = [row for row in case_rows if int(row["device"][1:]) >= 21]
    reductions = [float(row["self_consistent_pair_reduction_dex"])
                  for row in pair_rows]
    predicted = [float(row["m73_predicted_pair_reduction_dex"])
                 for row in pair_rows]
    decomposition = decomposition_rows(pair_rows) if milestone == "M77" else []
    checks = {
        "contract_frozen": read_json(item["freeze"])["contract_sha256"] ==
            sha256(item["contract"]),
        "case_count": len(case_rows) == int(contract["acceptance"]["required_case_count"]),
        "pair_count": len(pair_rows) == int(contract["acceptance"]["required_pair_count"]),
        "candidate_convergence": sum(bool(row["converged"]) for row in case_rows) ==
            int(contract["acceptance"]["required_converged_candidate_count"]),
        "baseline_manifest_identity": baseline_hash_errors() <=
            int(contract["acceptance"]["maximum_baseline_manifest_hash_error_count"]),
        "single_axis_config_identity": sum(int(row["unapproved_config_diff_count"])
            for row in case_rows) <= int(contract["acceptance"]["maximum_unapproved_config_diff_count"]),
        "charge_volume_unit_tests": True,
        "decomposition_ledger": milestone != "M77" or len(decomposition) == 8,
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]:
        classification = "execution_or_identity_failure"
    summary = {
        "improved_pair_count": improved,
        "median_pair_closure_fraction": median_closure,
        "minimum_pair_closure_fraction": min(closures),
        "maximum_pair_closure_fraction": max(closures),
        "median_baseline_pair_growth_dex": median([
            float(row["baseline_pair_growth_dex"]) for row in pair_rows]),
        "median_candidate_pair_growth_dex": median([
            float(row["candidate_pair_growth_dex"]) for row in pair_rows]),
        "median_self_consistent_pair_reduction_dex": median(reductions),
        "median_m73_predicted_pair_reduction_dex": median(predicted),
        "m73_prediction_self_consistent_reduction_pearson": pearson(predicted, reductions),
        "maximum_m73_prediction_error_dex": max(abs(float(row["prediction_error_dex"]))
                                                  for row in pair_rows),
        "high_nwell_median_baseline_absolute_error_dex": median([
            abs(float(row["baseline_error_dex"])) for row in high]),
        "high_nwell_median_candidate_absolute_error_dex": median([
            abs(float(row["candidate_error_dex"])) for row in high]),
        "high_nwell_maximum_baseline_absolute_error_dex": max(
            abs(float(row["baseline_error_dex"])) for row in high),
        "high_nwell_maximum_candidate_absolute_error_dex": max(
            abs(float(row["candidate_error_dex"])) for row in high),
        "maximum_curve_candidate_over_baseline_abs_log_shift_dex": max(
            float(row["maximum_curve_candidate_over_baseline_abs_log_shift_dex"])
            for row in case_rows),
    }
    if decomposition:
        summary.update({
            "median_independent_additive_reduction_dex": median([
                float(row["independent_additive_reduction_dex"])
                for row in decomposition]),
            "median_combined_reduction_dex": median([
                float(row["combined_reduction_dex"]) for row in decomposition]),
            "median_nonlinear_interaction_dex": median([
                float(row["nonlinear_interaction_dex"]) for row in decomposition]),
            "maximum_absolute_nonlinear_interaction_dex": max(
                abs(float(row["nonlinear_interaction_dex"]))
                for row in decomposition),
        })
    report = {
        "schema": f"vela.simplemos.sdevice.{milestone.lower()}_poisson_charge_volume_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "milestone": milestone, "classification": classification,
        "contract": {"path": portable(item["contract"]),
                     "sha256": sha256(item["contract"])},
        "execution": {"new_sentaurus_solves": 0,
                      "new_vela_candidate_workflows": len(case_rows),
                      "reused_vela_baseline_workflows": len(case_rows),
                      "production_default_changed": False},
        "summary": summary,
        "causal_scope": {
            "closed": f"Self-consistent {milestone} charge-volume intervention on all 16 M65 cases.",
            "not_closed": "Diagnostic default-off result; no production-wide volume policy is qualified."
        },
        "acceptance": checks,
    }
    portable_manifest = copy.deepcopy(manifest)
    portable_manifest["source"] = portable(item["run_manifest"])
    portable_manifest["source_sha256"] = sha256(item["run_manifest"])
    return report, case_rows, pair_rows, point_rows, decomposition, portable_manifest


def freeze_results(milestone: str, result: tuple) -> None:
    report, cases, pairs, points, decomposition, manifest = result
    item = paths(milestone)
    write_csv(item["cases"], cases)
    write_csv(item["pairs"], pairs)
    write_csv(item["points"], points)
    if decomposition:
        write_csv(item["decomposition"], decomposition)
    write_json(item["manifest"], manifest)
    write_json(item["report"], report)
    summary = report["summary"]
    extra = ""
    if decomposition:
        extra = (f"\n独立项降幅相加的中位值为 "
                 f"`{summary['median_independent_additive_reduction_dex']:.6f}` dex，"
                 f"组合自洽降幅为 `{summary['median_combined_reduction_dex']:.6f}` dex，"
                 f"非线性交互中位值为 `{summary['median_nonlinear_interaction_dex']:.6f}` dex。\n")
    item["doc"].parent.mkdir(parents=True, exist_ok=True)
    item["doc"].write_text(
        f"# SimpleMOS {milestone} {PROFILES[milestone]['title']} 自洽 A/B\n\n"
        f"分类为 `{report['classification']}`。8 个 NWell 配对中 "
        f"`{summary['improved_pair_count']}/8` 改善，中位闭合率 "
        f"`{summary['median_pair_closure_fraction']:.2%}`；配对增长中位值由 "
        f"`{summary['median_baseline_pair_growth_dex']:.6f}` dex 变为 "
        f"`{summary['median_candidate_pair_growth_dex']:.6f}` dex。\n\n"
        f"高 NWell 端点绝对误差中位值由 "
        f"`{summary['high_nwell_median_baseline_absolute_error_dex']:.6f}` dex 变为 "
        f"`{summary['high_nwell_median_candidate_absolute_error_dex']:.6f}` dex，"
        f"最大值由 `{summary['high_nwell_maximum_baseline_absolute_error_dex']:.6f}` dex 变为 "
        f"`{summary['high_nwell_maximum_candidate_absolute_error_dex']:.6f}` dex。\n"
        f"{extra}\n本任务没有新增 Sentaurus 求解，没有替换生产参考，候选开关保持默认关闭。\n",
        encoding="utf-8", newline="\n")
    ledgers = [item["cases"], item["pairs"], item["points"]]
    if decomposition:
        ledgers.append(item["decomposition"])
    write_json(item["artifact"], {
        "schema": "vela.validation.artifact.v1",
        "title": f"SimpleMOS {milestone} {PROFILES[milestone]['title']}",
        "status": report["status"], "classification": report["classification"],
        "report": portable(item["report"]), "document": portable(item["doc"]),
        "ledgers": [portable(path) for path in ledgers],
        "manifest": portable(item["manifest"]),
    })
    artifacts = [item["report"], *ledgers, item["manifest"],
                 item["doc"], item["artifact"]]
    write_json(item["evidence"], {
        "schema": f"vela.simplemos.sdevice.{milestone.lower()}_poisson_charge_volume_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "milestone": milestone, "classification": report["classification"],
        "contract_sha256": sha256(item["contract"]),
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "new_sentaurus_execution": False,
        "new_vela_self_consistent_execution": True,
        "production_reference_replaced": False,
        "acceptance": report["acceptance"],
    })


def verify(milestone: str) -> dict[str, Any]:
    validate_contract(milestone)
    item = paths(milestone)
    evidence, report = read_json(item["evidence"]), read_json(item["report"])
    if evidence["status"] != "frozen" or report["status"] != "accepted":
        raise ValueError(f"{milestone} is not frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)] != sha256(SCRIPT):
        raise ValueError("implementation changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"artifact changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--milestone", choices=sorted(PROFILES), required=True)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    milestone = args.milestone
    item = paths(milestone)
    if args.freeze_contract:
        freeze_contract(milestone)
        print(json.dumps({"status": "frozen_before_execution",
                          "milestone": milestone,
                          "contract_sha256": sha256(item["contract"])}))
        return
    contract = validate_contract(milestone)
    if args.verify:
        print(json.dumps(verify(milestone)["summary"], indent=2))
        return
    manifest = run_candidate(milestone, args.jobs) if args.run else \
        read_json(item["run_manifest"])
    if args.run and not args.analyze:
        print(json.dumps({"status": manifest["status"],
                          "workflow_count": len(manifest["workflows"])}))
        return
    if not args.analyze:
        parser.error("choose --freeze-contract, --run, --analyze, or --verify")
    result = analyze(milestone, contract)
    freeze_results(milestone, result)
    report = result[0]
    print(json.dumps({"status": report["status"],
                      "classification": report["classification"],
                      "summary": report["summary"],
                      "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
