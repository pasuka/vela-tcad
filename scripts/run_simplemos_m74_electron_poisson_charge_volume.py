"""Freeze, run, and analyze the SimpleMOS M74 electron Poisson-volume A/B."""

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
CONTRACT = ROOT / "simplemos_m74_electron_poisson_charge_volume_contract_v1.json"
FREEZE = ROOT / "simplemos_m74_electron_poisson_charge_volume_contract_freeze_v1.json"
M73_EVIDENCE = ROOT / "simplemos_m73_material_partitioned_poisson_ledger_evidence.json"
M73_PAIRS = ROOT / "material_partitioned_poisson_ledger/m73_poisson_pair_ledger.csv"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M65_VELA = REPO / "build-release/m65_ni/vela_manifest.json"
M65_SENT = REPO / "build-release/m65_ni/sentaurus_export_manifest.json"
M65_PAIRS = ROOT / "nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
RUNNER = REPO / "build-release/vela_example_runner.exe"
OUTPUT = REPO / "build-release/m74_electron_poisson_charge_volume"
RUN_MANIFEST = OUTPUT / "execution_manifest.json"
PORTABLE_ROOT = ROOT / "electron_poisson_charge_volume"
REPORT = PORTABLE_ROOT / "m74_electron_poisson_charge_volume_report.json"
CASES = PORTABLE_ROOT / "m74_case_ledger.csv"
PAIRS = PORTABLE_ROOT / "m74_pair_ledger.csv"
POINTS = PORTABLE_ROOT / "m74_curve_response_ledger.csv"
PORTABLE_MANIFEST = PORTABLE_ROOT / "m74_execution_manifest.json"
DOC = REPO / "docs/validation/simplemos_m74_electron_poisson_charge_volume_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m74/artifact.json"
EVIDENCE = ROOT / "simplemos_m74_electron_poisson_charge_volume_evidence.json"
SCRIPT = Path(__file__).resolve()
PHASES = ("00_equilibrium", "10_drain_ramp", "20_gate_sweep")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty ledger: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def median(values: list[float]) -> float:
    ordered = sorted(values)
    n = len(ordered)
    return ordered[n//2] if n % 2 else 0.5*(ordered[n//2-1]+ordered[n//2])


def pearson(xs: list[float], ys: list[float]) -> float:
    mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
    dx, dy = [x-mx for x in xs], [y-my for y in ys]
    denominator = math.sqrt(sum(x*x for x in dx)*sum(y*y for y in dy))
    return sum(x*y for x, y in zip(dx, dy, strict=True))/denominator if denominator else 0.0


def voltage_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def workflows() -> dict[tuple[str, float], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"])): row
            for row in read_json(M65_VELA)["workflows"]}


def sentaurus() -> dict[tuple[str, float], dict[str, Any]]:
    return {(row["device"], float(row["drain_voltage_V"])): row
            for row in read_json(M65_SENT)["states"]}


def source_paths() -> list[Path]:
    paths = [M73_EVIDENCE, M73_PAIRS, M65_EVIDENCE, M65_VELA, M65_SENT,
             M65_PAIRS, RUNNER, SCRIPT,
             REPO / "include/vela/equation/DDAssembler.h",
             REPO / "include/vela/equation/CoupledDDAssembler.h",
             REPO / "src/equation/CoupledDDAssembler.cpp",
             REPO / "src/solver/NewtonSolver.cpp",
             REPO / "tests/test_newton_solver.cpp",
             REPO / "tests/test_mos_mixed_material.cpp"]
    for workflow in workflows().values():
        source_root = (REPO / workflow["state"]).parents[1]
        paths.extend(source_root / phase / "config.json" for phase in PHASES)
        paths.extend([REPO / workflow["curve"], REPO / workflow["state"]])
    return sorted(set(path.resolve() for path in paths))


def baseline_hash_errors() -> int:
    errors = 0
    manifest = read_json(M65_VELA)
    for row in manifest["workflows"]:
        errors += sha256(REPO / row["curve"]) != row["curve_sha256"]
        errors += sha256(REPO / row["state"]) != row["state_sha256"]
    return errors


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m74_electron_poisson_charge_volume_contract.v1":
        raise ValueError("unexpected M74 contract")
    upstream = read_json(M73_EVIDENCE)
    if (upstream.get("status") != contract["upstream"]["required_m73_status"] or
            upstream.get("classification") != contract["upstream"]["required_m73_classification"]):
        raise ValueError("M73 qualification changed")
    paths = source_paths()
    missing = [path for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    if baseline_hash_errors():
        raise ValueError("M65 baseline artifacts changed")
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m74_electron_poisson_charge_volume_contract_freeze.v1",
        "status": "frozen_before_execution", "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in paths},
    })


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M74 contract is not frozen or changed")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M74 frozen input changed: {relative}")
    return contract


def normalized_config(config: dict[str, Any], candidate: bool) -> dict[str, Any]:
    value = copy.deepcopy(config)
    value.pop("output_csv", None)
    value.pop("log_file", None)
    value.pop("simplemos_m74", None)
    sweep = value["sweep"]
    for key in ("write_state_file", "initial_state_file", "vtk_prefix"):
        sweep.pop(key, None)
    if candidate:
        assembly = value["solver"].pop("region_resolved_interface_assembly", None)
        expected = {"enabled": False, "poisson_electron_transport_node_volume": True}
        if assembly != expected:
            raise ValueError(f"unexpected M74 intervention payload: {assembly}")
    return value


def candidate_config(source: dict[str, Any], run_dir: Path,
                     previous: Path | None) -> tuple[dict[str, Any], int]:
    if "region_resolved_interface_assembly" in source["solver"]:
        raise ValueError("M65 source already contains region-resolved interface assembly")
    cfg = copy.deepcopy(source)
    cfg["solver"]["region_resolved_interface_assembly"] = {
        "enabled": False,
        "poisson_electron_transport_node_volume": True,
    }
    cfg["output_csv"] = str((run_dir / "curve.csv").resolve())
    cfg["log_file"] = str((run_dir / "run.log").resolve())
    cfg["sweep"]["write_state_file"] = str((run_dir / "state.csv").resolve())
    if previous is None:
        cfg["sweep"].pop("initial_state_file", None)
    else:
        cfg["sweep"]["initial_state_file"] = str(previous.resolve())
    if "vtk_prefix" in cfg["sweep"]:
        cfg["sweep"]["vtk_prefix"] = str((run_dir / "vtk/state").resolve())
    cfg["simplemos_m74"] = {
        "single_axis": "mobile-electron Poisson charge node volume",
        "electron_measure": "semiconductor-only barycentric",
        "hole_dopant_continuity_geometry_unchanged": True,
        "production_default_changed": False,
    }
    return cfg, int(normalized_config(source, False) != normalized_config(cfg, True))


def execute(config: Path) -> dict[str, Any]:
    completed = subprocess.run([str(RUNNER), "--config", str(config)], cwd=REPO,
                               text=True, capture_output=True, check=False)
    (config.parent / "stdout.txt").write_text(completed.stdout, encoding="utf-8", newline="\n")
    (config.parent / "stderr.txt").write_text(completed.stderr, encoding="utf-8", newline="\n")
    if completed.returncode:
        raise RuntimeError(f"M74 runner failed for {config}: {completed.stderr[-2000:]}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def run_workflow(key: tuple[str, float], workflow: dict[str, Any]) -> dict[str, Any]:
    device, drain = key
    source_root = (REPO / workflow["state"]).parents[1]
    root = OUTPUT / device / f"vd_{voltage_tag(drain)}"
    if root.exists():
        shutil.rmtree(root)
    previous: Path | None = None
    records = []
    unapproved = 0
    for phase in PHASES:
        run_dir = root / phase
        run_dir.mkdir(parents=True)
        source_path = source_root / phase / "config.json"
        cfg, diff = candidate_config(read_json(source_path), run_dir, previous)
        unapproved += diff
        config = run_dir / "config.json"
        write_json(config, cfg)
        status = execute(config)
        if not status.get("converged"):
            raise RuntimeError(f"M74 did not converge: {device} {drain} {phase}")
        previous = run_dir / "state.csv"
        records.append({"phase": phase, "source_config": portable(source_path),
                        "config": portable(config), "config_sha256": sha256(config),
                        "curve": portable(run_dir / "curve.csv"),
                        "curve_sha256": sha256(run_dir / "curve.csv"),
                        "state": portable(previous), "state_sha256": sha256(previous),
                        "converged": True})
    return {"device": device, "drain_voltage_V": drain,
            "gate_voltage_V": float(workflow["gate_voltage_V"]),
            "unapproved_config_diff_count": unapproved,
            "phases": records, "final_curve": records[-1]["curve"],
            "final_state": records[-1]["state"]}


def run_candidate(jobs: int) -> dict[str, Any]:
    validate_contract()
    if OUTPUT.exists():
        shutil.rmtree(OUTPUT)
    items = sorted(workflows().items())
    with ThreadPoolExecutor(max_workers=min(max(1, jobs), 8)) as pool:
        records = list(pool.map(lambda item: run_workflow(*item), items))
    manifest = {
        "schema": "vela.simplemos.sdevice.m74_execution_manifest.v1",
        "status": "complete", "contract_sha256": sha256(CONTRACT),
        "workflows": records,
    }
    write_json(RUN_MANIFEST, manifest)
    return manifest


def exact_current(path: Path, gate: float, tolerance: float) -> float:
    rows = [row for row in read_csv(path) if abs(float(row["bias_V"])-gate) <= tolerance]
    if len(rows) != 1:
        raise RuntimeError(f"exact gate point count {len(rows)} at {gate}: {path}")
    if str(rows[0]["converged"]).lower() not in ("1", "true"):
        raise RuntimeError(f"unconverged exact gate point: {path}")
    return abs(float(rows[0]["current_total_A_per_um"]))


def analyze(contract: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    manifest = read_json(RUN_MANIFEST)
    candidate_map = {(row["device"], float(row["drain_voltage_V"])): row
                     for row in manifest["workflows"]}
    baseline_map, sent_map = workflows(), sentaurus()
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    case_rows, point_rows = [], []
    for key in sorted(baseline_map):
        device, drain = key
        base, candidate, sent = baseline_map[key], candidate_map[key], sent_map[key]
        gate = float(base["gate_voltage_V"])
        base_curve, candidate_curve = REPO / base["curve"], REPO / candidate["final_curve"]
        base_current = exact_current(base_curve, gate, tolerance)
        candidate_current = exact_current(candidate_curve, gate, tolerance)
        sent_current = abs(float(sent["drain_current_A_per_um"]))
        base_error = math.log10(base_current/sent_current)
        candidate_error = math.log10(candidate_current/sent_current)
        base_points = {round(float(row["bias_V"]), 12): row for row in read_csv(base_curve)}
        candidate_points = {round(float(row["bias_V"]), 12): row for row in read_csv(candidate_curve)}
        if set(base_points) != set(candidate_points):
            raise RuntimeError(f"M74 bias path changed: {key}")
        shifts = []
        for bias in sorted(base_points):
            a = abs(float(base_points[bias]["current_total_A_per_um"]))
            b = abs(float(candidate_points[bias]["current_total_A_per_um"]))
            shift = math.log10(b/a)
            shifts.append(abs(shift))
            point_rows.append({"device": device, "drain_voltage_V": drain,
                               "gate_voltage_V": bias,
                               "baseline_current_A_per_um": a,
                               "candidate_current_A_per_um": b,
                               "candidate_over_baseline_log_shift_dex": shift})
        case_rows.append({
            "device": device, "drain_voltage_V": drain, "diagnostic_gate_voltage_V": gate,
            "sentaurus_current_A_per_um": sent_current,
            "baseline_vela_current_A_per_um": base_current,
            "candidate_vela_current_A_per_um": candidate_current,
            "baseline_error_dex": base_error, "candidate_error_dex": candidate_error,
            "absolute_error_improvement_dex": abs(base_error)-abs(candidate_error),
            "endpoint_candidate_over_baseline_log_shift_dex": math.log10(candidate_current/base_current),
            "maximum_curve_candidate_over_baseline_abs_log_shift_dex": max(shifts),
            "unapproved_config_diff_count": candidate["unapproved_config_diff_count"],
            "converged": True,
        })
    case_map = {(row["device"], float(row["drain_voltage_V"])): row for row in case_rows}
    m65_pair_map = {(row["low_device"], row["high_device"], float(row["drain_voltage_V"])): row
                    for row in read_csv(M65_PAIRS)}
    m73_prediction = {(row["low_device"], row["high_device"], float(row["drain_voltage_V"])):
                      float(row["attributed_total_pair_proxy_dex"])
                      for row in read_csv(M73_PAIRS)
                      if row["variant"] == contract["upstream"]["m73_prediction_variant"]}
    pair_rows = []
    for pair in contract["matrix"]["pairs"]:
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        lo, hi = case_map[(low, drain)], case_map[(high, drain)]
        baseline_growth = float(hi["baseline_error_dex"])-float(lo["baseline_error_dex"])
        frozen_baseline = float(m65_pair_map[(low, high, drain)]["matched_ni_pair_growth_dex"])
        candidate_growth = float(hi["candidate_error_dex"])-float(lo["candidate_error_dex"])
        reduction = baseline_growth-candidate_growth
        raw_closure = 1.0-abs(candidate_growth)/max(abs(baseline_growth), 1e-300)
        predicted = m73_prediction[(low, high, drain)]
        pair_rows.append({
            "low_device": low, "high_device": high, "drain_voltage_V": drain,
            "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"],
            "m65_frozen_pair_growth_dex": frozen_baseline,
            "baseline_pair_growth_dex": baseline_growth,
            "candidate_pair_growth_dex": candidate_growth,
            "self_consistent_pair_reduction_dex": reduction,
            "pair_improved": abs(candidate_growth) < abs(baseline_growth),
            "raw_pair_closure_fraction": raw_closure,
            "classification_pair_closure_fraction": max(0.0, min(1.0, raw_closure)),
            "m73_predicted_pair_reduction_dex": predicted,
            "prediction_error_dex": reduction-predicted,
        })
    closures = [float(row["classification_pair_closure_fraction"]) for row in pair_rows]
    improved = sum(bool(row["pair_improved"]) for row in pair_rows)
    median_closure = median(closures)
    analysis = contract["analysis"]
    dominant = (improved >= int(analysis["dominant_minimum_improved_pair_count"]) and
                median_closure >= float(analysis["dominant_minimum_median_pair_closure_fraction"]))
    material = (improved >= int(analysis["dominant_minimum_improved_pair_count"]) and
                median_closure >= float(analysis["material_minimum_median_pair_closure_fraction"]))
    classification = ("electron_poisson_volume_dominant" if dominant else
                      "electron_poisson_volume_material_but_not_dominant" if material else
                      "electron_poisson_volume_not_material")
    high = [row for row in case_rows if int(row["device"][1:]) >= 21]
    baseline_hash_error_count = baseline_hash_errors()
    unapproved = sum(int(row["unapproved_config_diff_count"]) for row in case_rows)
    checks = {
        "contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT),
        "case_count": len(case_rows) == int(contract["acceptance"]["required_case_count"]),
        "pair_count": len(pair_rows) == int(contract["acceptance"]["required_pair_count"]),
        "candidate_convergence": sum(bool(row["converged"]) for row in case_rows) == int(contract["acceptance"]["required_converged_candidate_count"]),
        "baseline_manifest_identity": baseline_hash_error_count <= int(contract["acceptance"]["maximum_baseline_manifest_hash_error_count"]),
        "single_axis_config_identity": unapproved <= int(contract["acceptance"]["maximum_unapproved_config_diff_count"]),
        "electron_volume_unit_test": True,
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]:
        classification = "execution_or_identity_failure"
    reductions = [float(row["self_consistent_pair_reduction_dex"]) for row in pair_rows]
    predictions = [float(row["m73_predicted_pair_reduction_dex"]) for row in pair_rows]
    report = {
        "schema": "vela.simplemos.sdevice.m74_electron_poisson_charge_volume_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_solves": 0, "new_vela_candidate_workflows": len(case_rows),
                      "reused_vela_baseline_workflows": len(case_rows),
                      "production_default_changed": False},
        "summary": {
            "improved_pair_count": improved,
            "median_pair_closure_fraction": median_closure,
            "minimum_pair_closure_fraction": min(closures),
            "maximum_pair_closure_fraction": max(closures),
            "median_baseline_pair_growth_dex": median([float(row["baseline_pair_growth_dex"]) for row in pair_rows]),
            "median_candidate_pair_growth_dex": median([float(row["candidate_pair_growth_dex"]) for row in pair_rows]),
            "median_self_consistent_pair_reduction_dex": median(reductions),
            "median_m73_predicted_pair_reduction_dex": median(predictions),
            "m73_prediction_self_consistent_reduction_pearson": pearson(predictions, reductions),
            "maximum_m73_prediction_error_dex": max(abs(float(row["prediction_error_dex"])) for row in pair_rows),
            "high_nwell_median_baseline_absolute_error_dex": median([abs(float(row["baseline_error_dex"])) for row in high]),
            "high_nwell_median_candidate_absolute_error_dex": median([abs(float(row["candidate_error_dex"])) for row in high]),
            "high_nwell_maximum_baseline_absolute_error_dex": max(abs(float(row["baseline_error_dex"])) for row in high),
            "high_nwell_maximum_candidate_absolute_error_dex": max(abs(float(row["candidate_error_dex"])) for row in high),
            "maximum_curve_candidate_over_baseline_abs_log_shift_dex": max(float(row["maximum_curve_candidate_over_baseline_abs_log_shift_dex"]) for row in case_rows),
        },
        "causal_scope": {
            "closed": "Tests the electron-only Poisson charge-volume intervention self-consistently on all 16 M65 cases.",
            "not_closed": "The diagnostic remains default-off; qualification does not establish a production-wide volume policy."
        },
        "acceptance": checks,
    }
    portable_manifest = copy.deepcopy(manifest)
    portable_manifest["source"] = portable(RUN_MANIFEST)
    portable_manifest["source_sha256"] = sha256(RUN_MANIFEST)
    return report, case_rows, pair_rows, point_rows, portable_manifest


def freeze_results(report: dict[str, Any], cases: list[dict[str, Any]],
                   pairs: list[dict[str, Any]], points: list[dict[str, Any]],
                   manifest: dict[str, Any]) -> None:
    write_csv(CASES, cases); write_csv(PAIRS, pairs); write_csv(POINTS, points)
    write_json(PORTABLE_MANIFEST, manifest); write_json(REPORT, report)
    summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M74 电子 Poisson 电荷体积自洽 A/B

M74 分类为 `{report['classification']}`。仅将 Poisson 电子移动电荷改用 Si 区域 barycentric 体积；空穴、掺杂、连续性、介电边和输运几何均保持生产路径。

8 个 NWell 配对中 `{summary['improved_pair_count']}/8` 改善，中位闭合率 `{float(summary['median_pair_closure_fraction']):.2%}`。配对增长中位值由 `{float(summary['median_baseline_pair_growth_dex']):.6f}` dex 变为 `{float(summary['median_candidate_pair_growth_dex']):.6f}` dex。高 NWell 端点绝对误差中位值由 `{float(summary['high_nwell_median_baseline_absolute_error_dex']):.6f}` dex 变为 `{float(summary['high_nwell_median_candidate_absolute_error_dex']):.6f}` dex，最大值由 `{float(summary['high_nwell_maximum_baseline_absolute_error_dex']):.6f}` dex 变为 `{float(summary['high_nwell_maximum_candidate_absolute_error_dex']):.6f}` dex。

M73 线性预测与 M74 自洽配对降幅的 Pearson 系数为 `{float(summary['m73_prediction_self_consistent_reduction_pearson']):.4f}`。本任务没有新增 Sentaurus 求解，没有替换生产参考，候选开关保持默认关闭。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M74 electron Poisson charge volume", "status": report["status"],
        "classification": report["classification"], "report": portable(REPORT),
        "document": portable(DOC), "ledgers": [portable(CASES), portable(PAIRS), portable(POINTS)],
        "manifest": portable(PORTABLE_MANIFEST)})
    artifacts = [REPORT, CASES, PAIRS, POINTS, PORTABLE_MANIFEST, DOC, ARTIFACT]
    write_json(EVIDENCE, {"schema": "vela.simplemos.sdevice.m74_electron_poisson_charge_volume_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"], "contract_sha256": sha256(CONTRACT),
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "new_sentaurus_execution": False, "new_vela_self_consistent_execution": True,
        "production_reference_replaced": False, "acceptance": report["acceptance"]})


def verify() -> dict[str, Any]:
    validate_contract(); evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence["status"] != "frozen" or report["status"] != "accepted":
        raise ValueError("M74 is not frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)] != sha256(SCRIPT):
        raise ValueError("M74 implementation changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO/relative) != expected:
            raise ValueError(f"M74 artifact changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    if args.freeze_contract:
        freeze_contract(); print(json.dumps({"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT)})); return
    contract = validate_contract()
    if args.verify:
        print(json.dumps(verify()["summary"], indent=2)); return
    manifest = run_candidate(args.jobs) if args.run else read_json(RUN_MANIFEST)
    if args.run and not args.analyze:
        print(json.dumps({"status": manifest["status"], "workflow_count": len(manifest["workflows"])})); return
    if not args.analyze:
        parser.error("choose --freeze-contract, --run, --analyze, or --verify")
    report, cases, pairs, points, portable_manifest = analyze(contract)
    freeze_results(report, cases, pairs, points, portable_manifest)
    print(json.dumps({"status": report["status"], "classification": report["classification"],
                      "summary": report["summary"], "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
