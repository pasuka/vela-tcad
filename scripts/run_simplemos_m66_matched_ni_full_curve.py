"""Execute and freeze SimpleMOS M66 matched-ni full-curve continuation."""

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
import sys
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m60_tight_convergence_port_burst as m60  # noqa: E402
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m66_matched_ni_full_curve_contract_v1.json"
FREEZE = ROOT / "simplemos_m66_matched_ni_full_curve_contract_freeze.json"
M65_EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
M64_EVIDENCE = ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_evidence.json"
M65_REPORT = ROOT / "nobgn_intrinsic_density_attribution/m65_nobgn_intrinsic_density_attribution_report.json"
M65_MANIFEST = REPO / "build-release/m65_ni/vela_manifest.json"
M65_MATERIALS = REPO / "build-release/m65_ni/vela/matched_materials.json"
M64_POINTS = ROOT / "bgn_smooth_nwell_intervention/m64_bgn_on_off_point_ledger.csv"
M64_CASES = ROOT / "bgn_smooth_nwell_intervention/m64_execution_case_ledger.csv"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
OUTPUT = REPO / "build-release/m66_full"
PORTABLE = ROOT / "matched_ni_full_curve"
REPORT = PORTABLE / "m66_matched_ni_full_curve_report.json"
POINTS = PORTABLE / "m66_matched_ni_full_curve_point_ledger.csv"
PAIRS = PORTABLE / "m66_matched_ni_full_curve_pair_ledger.csv"
CASES = PORTABLE / "m66_execution_case_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m66_matched_ni_full_curve_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m66/artifact.json"
EVIDENCE = ROOT / "simplemos_m66_matched_ni_full_curve_evidence.json"
SCRIPT = Path(__file__).resolve()
RUNNER = REPO / "build-release/vela_example_runner.exe"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                               allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"empty CSV: {path}")
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
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else 0.5 * (ordered[middle - 1] + ordered[middle])


def voltage_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def source_paths() -> list[Path]:
    manifest = read_json(M65_MANIFEST)
    paths = [M65_EVIDENCE, M65_REPORT, M65_MANIFEST, M65_MATERIALS, M64_POINTS, RUNNER]
    for row in manifest["workflows"]:
        paths.extend([REPO / row["curve"], REPO / row["state"]])
        tag = "vd_0p05" if math.isclose(float(row["drain_voltage_V"]), 0.05) else "vd_1"
        paths.append(M8 / "vela" / row["device"] / "workflow" / tag / "20_gate_sweep.json")
    return paths


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m66_matched_ni_full_curve_contract.v1":
        raise ValueError("unexpected M66 contract schema")
    upstream = read_json(M65_EVIDENCE)
    if (upstream.get("status") != contract["upstream"]["required_m65_status"] or
            upstream.get("classification") != contract["upstream"]["required_m65_classification"]):
        raise ValueError("M65 upstream qualification changed")
    sources = source_paths()
    missing = [path for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m66_matched_ni_full_curve_contract_freeze.v1",
        "status": "frozen_before_execution",
        "contract": portable(CONTRACT), "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in sources},
    })


def validate_contract() -> dict[str, Any]:
    contract, freeze = read_json(CONTRACT), read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M66 contract was not frozen before execution")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M66 frozen input changed: {relative}")
    return contract


def exact(rows: list[dict[str, str]], gate: float) -> dict[str, str]:
    matches = [row for row in rows if abs(float(row["bias_V"]) - gate) <= 1e-10]
    if len(matches) != 1:
        raise RuntimeError(f"exact Vg={gate} count {len(matches)}")
    return matches[0]


def continue_case(row: dict[str, Any], force: bool) -> dict[str, Any]:
    device, drain, gate = row["device"], float(row["drain_voltage_V"]), float(row["gate_voltage_V"])
    tag = "vd_0p05" if math.isclose(drain, 0.05) else "vd_1"
    root = OUTPUT / device / f"vd_{voltage_tag(drain)}"
    curve, state, stdout = root / "curve.csv", root / "state.csv", root / "stdout.txt"
    complete = curve.is_file() and state.is_file() and stdout.is_file() and '"converged":true' in stdout.read_text(errors="replace")
    if root.exists() and (force or not complete):
        shutil.rmtree(root)
    root.mkdir(parents=True, exist_ok=True)
    source = M8 / "vela" / device / "workflow" / tag / "20_gate_sweep.json"
    config = copy.deepcopy(read_json(source))
    config["materials_file"] = str(M65_MATERIALS.resolve())
    config["solver"]["bandgap_narrowing"] = {"model": "none", "fermi_statistics_correction": False}
    config["contacts"] = [{**item, "bias": gate if item["name"] == "gate" else item["bias"]}
                          for item in config["contacts"]]
    config["output_csv"] = str(curve.resolve())
    config["log_file"] = str((root / "run.log").resolve())
    config["sweep"]["start"] = gate
    config["sweep"]["stop"] = 2.5
    config["sweep"]["initial_state_file"] = str((REPO / row["state"]).resolve())
    config["sweep"]["write_state_file"] = str(state.resolve())
    config["sweep"]["write_vtk"] = False
    config["sweep"].pop("vtk_prefix", None)
    config["sweep"]["bias_points"] = [round(0.05 * index, 12) for index in range(51)
                                       if 0.05 * index >= gate - 1e-12]
    config["simplemos_m66"] = {"restart_gate_V": gate, "prefix_reused_from_m65": True,
                                "new_equilibrium_or_drain_solve": False,
                                "production_default_changed": False}
    config_path = root / "config.json"
    write_json(config_path, config)
    if not complete:
        completed = subprocess.run([str(RUNNER), "--config", str(config_path)], cwd=REPO,
                                   text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   encoding="utf-8", errors="replace", check=False)
        stdout.write_text(completed.stdout or "", encoding="utf-8", newline="\n")
        if completed.returncode or '"converged":true' not in (completed.stdout or ""):
            raise RuntimeError(f"M66 continuation failed: {device} Vd={drain}\n{completed.stdout}")
    continuation = read_csv(curve)
    prefix = read_csv(REPO / row["curve"])
    replay_old = abs(float(exact(prefix, gate)["current_total_A_per_um"]))
    newton_start = abs(float(exact(continuation, gate)["current_total_A_per_um"]))
    newton_reclosure_error = math.log10(newton_start / replay_old)

    combined_rows = [item for item in prefix if float(item["bias_V"]) <= gate + 1e-10]
    combined_rows.extend(item for item in continuation if float(item["bias_V"]) > gate + 1e-10)
    combined = root / "combined_curve.csv"
    write_csv(combined, combined_rows)
    combined_restart = abs(float(exact(combined_rows, gate)["current_total_A_per_um"]))
    replay_error = math.log10(combined_restart / replay_old)
    return {"device": device, "drain_voltage_V": drain, "restart_gate_voltage_V": gate,
            "m65_prefix_curve": row["curve"], "m65_restart_state": row["state"],
            "continuation_curve": portable(curve), "combined_curve": portable(combined),
            "final_state": portable(state), "restart_current_replay_error_dex": replay_error,
            "newton_start_reclosure_error_dex": newton_reclosure_error,
            "point_count": len(combined_rows), "combined_curve_sha256": sha256(combined)}


def run_vela(contract: dict[str, Any], jobs: int, force: bool) -> dict[str, Any]:
    manifest = read_json(M65_MANIFEST)
    with ThreadPoolExecutor(max_workers=min(max(1, jobs), 8)) as pool:
        cases = list(pool.map(lambda row: continue_case(row, force), manifest["workflows"]))
    result = {"schema": "vela.simplemos.sdevice.m66_vela_manifest.v1", "status": "complete",
              "contract_sha256": sha256(CONTRACT), "cases": cases}
    write_json(OUTPUT / "manifest.json", result)
    return result


def local_slope(rows: list[dict[str, Any]], index: int, key: str) -> float:
    left, right = max(0, index - 1), min(len(rows) - 1, index + 1)
    return ((math.log10(float(rows[right][key])) - math.log10(float(rows[left][key]))) /
            (float(rows[right]["gate_voltage_V"]) - float(rows[left]["gate_voltage_V"])))


def fit_shift(rows: list[dict[str, Any]], low: float, high: float) -> tuple[float, float]:
    selected = [(local_slope(rows, index, "sentaurus_no_bgn_A_per_um"), float(row["error_dex"]))
                for index, row in enumerate(rows) if low <= float(row["gate_voltage_V"]) <= high]
    shift = sum(slope * error for slope, error in selected) / sum(slope * slope for slope, _ in selected)
    energy = sum(error * error for _, error in selected)
    residual = sum((error - slope * shift) ** 2 for slope, error in selected)
    return 1000.0 * shift, 1.0 - residual / energy if energy else 1.0


def analyze(contract: dict[str, Any], manifest: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    # M64's portable point ledger intentionally contains only the smooth
    # 0.55--0.90 V window.  Its frozen execution ledger seals the complete
    # 51-point raw PLT paths and hashes, so parse those for M66 and cross-check
    # the overlapping smooth points below.
    m64_evidence = read_json(M64_EVIDENCE)
    if m64_evidence["artifacts"][portable(M64_CASES)] != sha256(M64_CASES):
        raise RuntimeError("M64 execution ledger is not sealed by M64 evidence")
    gates = [round(0.05 * index, 12) for index in range(51)]
    sent: dict[tuple[str, float, float], float] = {}
    for row in read_csv(M64_CASES):
        if row["solver"] != "sentaurus":
            continue
        curve = REPO / row["curve"]
        if sha256(curve) != row["curve_sha256"]:
            raise RuntimeError(f"M64 raw Sentaurus curve hash changed: {curve}")
        parsed = m60.parse_current(curve, gates, 1e-10)
        for point in parsed:
            key = (row["device"], round(float(row["drain_voltage_V"]), 12),
                   round(float(point["gate_voltage_V"]), 12))
            sent[key] = abs(float(point["drain_total_A_per_um"]))
    for row in read_csv(M64_POINTS):
        key = (row["device"], round(float(row["drain_voltage_V"]), 12),
               round(float(row["gate_voltage_V"]), 12))
        if not math.isclose(sent[key], float(row["sentaurus_off_A_per_um"]), rel_tol=1e-14):
            raise RuntimeError(f"M64 smooth point/raw PLT mismatch: {key}")
    point_rows: list[dict[str, Any]] = []
    by_curve: dict[tuple[str, float], list[dict[str, Any]]] = {}
    case_rows: list[dict[str, Any]] = []
    for item in manifest["cases"]:
        curve = read_csv(REPO / item["combined_curve"])
        if len(curve) != int(contract["acceptance"]["required_points_per_curve"]):
            raise RuntimeError(f"M66 curve point count: {item}")
        rows = []
        for source in curve:
            gate = float(source["bias_V"])
            vela = abs(float(source["current_total_A_per_um"]))
            sentaurus = sent[(item["device"], round(float(item["drain_voltage_V"]), 12), round(gate, 12))]
            row = {"device": item["device"], "drain_voltage_V": item["drain_voltage_V"],
                   "gate_voltage_V": gate, "sentaurus_no_bgn_A_per_um": sentaurus,
                   "vela_matched_ni_A_per_um": vela, "error_dex": math.log10(vela / sentaurus),
                   "absolute_relative_error_percent": 100.0 * abs(vela / sentaurus - 1.0)}
            rows.append(row)
            point_rows.append(row)
        by_curve[(item["device"], float(item["drain_voltage_V"]))] = rows
        shift, explained = fit_shift(rows, *contract["matrix"]["smooth_gate_window_V"])
        case_rows.append({"device": item["device"], "drain_voltage_V": item["drain_voltage_V"],
                          "point_count": len(rows), "restart_gate_voltage_V": item["restart_gate_voltage_V"],
                          "restart_current_replay_error_dex": item["restart_current_replay_error_dex"],
                          "newton_start_reclosure_error_dex": item["newton_start_reclosure_error_dex"],
                          "smooth_fit_horizontal_shift_mV": shift,
                          "smooth_fit_explained_fraction": explained,
                          "maximum_absolute_error_dex": max(abs(float(row["error_dex"])) for row in rows),
                          "maximum_relative_error_percent": max(float(row["absolute_relative_error_percent"]) for row in rows)})
    pair_rows: list[dict[str, Any]] = []
    low_gate, high_gate = contract["matrix"]["smooth_gate_window_V"]
    for pair in contract["matrix"]["pairs"]:
        for drain in contract["matrix"]["drain_voltages_V"]:
            low_rows, high_rows = by_curve[(pair["low_device"], float(drain))], by_curve[(pair["high_device"], float(drain))]
            growths = [float(high["error_dex"]) - float(low["error_dex"])
                       for low, high in zip(low_rows, high_rows)
                       if low_gate <= float(low["gate_voltage_V"]) <= high_gate]
            pair_rows.append({"low_device": pair["low_device"], "high_device": pair["high_device"],
                              "drain_voltage_V": drain, "smooth_point_count": len(growths),
                              "median_smooth_pair_growth_dex": median(growths),
                              "maximum_absolute_smooth_pair_growth_dex": max(abs(value) for value in growths),
                              "minimum_smooth_pair_growth_dex": min(growths),
                              "maximum_smooth_pair_growth_dex": max(growths)})
    metric = median([abs(float(row["median_smooth_pair_growth_dex"])) for row in pair_rows])
    analysis = contract["analysis"]
    if metric <= float(analysis["dominant_maximum_median_absolute_pair_growth_dex"]):
        classification = "matched_ni_full_curve_dominant"
    elif metric <= float(analysis["material_maximum_median_absolute_pair_growth_dex"]):
        classification = "matched_ni_full_curve_material_but_not_complete"
    else:
        classification = "matched_ni_full_curve_not_material"
    maximum_replay = max(abs(float(row["restart_current_replay_error_dex"])) for row in case_rows)
    checks = {"contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT),
              "curve_count": len(case_rows) == int(contract["acceptance"]["required_curve_count_per_solver"]),
              "point_count": len(point_rows) == int(contract["acceptance"]["required_m64_sentaurus_point_count"]),
              "restart_identity": maximum_replay <= float(contract["acceptance"]["maximum_restart_current_replay_error_dex"]),
              "classification_declared": classification in analysis["classifications"],
              "production_reference_not_replaced": True}
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]:
        classification = "execution_or_restart_identity_failure"
    high = [row for row in point_rows if row["device"] in {"n21", "n22", "n23", "n24"}
            and low_gate <= float(row["gate_voltage_V"]) <= high_gate]
    report = {"schema": "vela.simplemos.sdevice.m66_matched_ni_full_curve_report.v1",
              "status": "accepted" if checks["all_checks_pass"] else "failed",
              "classification": classification,
              "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
              "execution": {"new_sentaurus_solves": 0, "new_vela_equilibrium_or_drain_solves": 0,
                            "vela_continuations": len(case_rows), "production_default_changed": False},
              "summary": {"curve_count": len(case_rows), "point_count": len(point_rows),
                          "maximum_restart_current_replay_error_dex": maximum_replay,
                          "maximum_newton_start_reclosure_error_dex": max(abs(float(row["newton_start_reclosure_error_dex"])) for row in case_rows),
                          "median_absolute_smooth_pair_growth_dex": metric,
                          "maximum_absolute_smooth_pair_growth_dex": max(float(row["maximum_absolute_smooth_pair_growth_dex"]) for row in pair_rows),
                          "high_nwell_smooth_median_absolute_error_dex": median([abs(float(row["error_dex"])) for row in high]),
                          "high_nwell_smooth_maximum_absolute_error_dex": max(abs(float(row["error_dex"])) for row in high),
                          "high_nwell_smooth_maximum_relative_error_percent": max(float(row["absolute_relative_error_percent"]) for row in high),
                          "median_smooth_fit_explained_fraction": median([float(row["smooth_fit_explained_fraction"]) for row in case_rows])},
              "causal_scope": {"closed": "M65 matched-ni behavior is evaluated on every native Id-Vg point without a new Sentaurus solve.",
                               "not_closed": "M66 tests representativeness only; remaining state and transport differences are assigned to M67-M70."},
              "acceptance": checks}
    return report, {"points": point_rows, "pairs": pair_rows, "cases": case_rows}


def freeze_results(report: dict[str, Any], rows: dict[str, Any]) -> None:
    write_csv(POINTS, rows["points"]); write_csv(PAIRS, rows["pairs"]); write_csv(CASES, rows["cases"])
    write_json(REPORT, report)
    summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M66 matched-ni 完整曲线复核

M66 分类为 `{report['classification']}`。16 条 Vela matched-ni 曲线均从 M65 冻结诊断点续算，形成每条 51 点的 0--2.5 V 曲线；未启动新的 Sentaurus、平衡态或漏压态求解。组合曲线在检查点拼接处逐位复用 M65 电流行，最大身份误差为 `{float(summary['maximum_restart_current_replay_error_dex']):.3e} dex`；Newton 续算起点的再次闭合变化单独记录，最大为 `{float(summary['maximum_newton_start_reclosure_error_dex']):.3e} dex`，不被隐去或误作拼接身份误差。

在 0.55--0.90 V 平滑窗口，8 个低/高 NWell 配对的中位绝对误差增幅为 `{float(summary['median_absolute_smooth_pair_growth_dex']):.6f} dex`，高 NWell 工况逐点最大绝对误差为 `{float(summary['high_nwell_smooth_maximum_absolute_error_dex']):.6f} dex`（相对误差 `{float(summary['high_nwell_smooth_maximum_relative_error_percent']):.2f}%`）。这确认 M65 的基础 ni 差异在完整曲线上是实质成分，但没有闭合剩余平滑差异。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1", "title": "SimpleMOS M66 matched-ni full curve",
                          "status": report["status"], "classification": report["classification"],
                          "report": portable(REPORT), "ledgers": [portable(POINTS), portable(PAIRS), portable(CASES)]})
    artifacts = [REPORT, POINTS, PAIRS, CASES, DOC, ARTIFACT]
    write_json(EVIDENCE, {"schema": "vela.simplemos.sdevice.m66_matched_ni_full_curve_evidence.v1",
                          "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
                          "classification": report["classification"], "contract_sha256": sha256(CONTRACT),
                          "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
                          "artifacts": {portable(path): sha256(path) for path in artifacts},
                          "new_sentaurus_execution": False, "production_reference_replaced": False,
                          "acceptance": report["acceptance"]})


def verify() -> dict[str, Any]:
    validate_contract(); evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence["status"] != "frozen" or report["status"] != "accepted":
        raise ValueError("M66 is not accepted and frozen")
    if evidence["implementation_hashes"][portable(SCRIPT)] != sha256(SCRIPT):
        raise ValueError("M66 implementation changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M66 artifact changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    if args.freeze_contract:
        freeze_contract(); print(json.dumps({"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT)})); return
    contract = validate_contract()
    if args.verify:
        print(json.dumps(verify()["summary"], indent=2)); return
    manifest = run_vela(contract, args.jobs, args.force) if args.run_vela else read_json(OUTPUT / "manifest.json")
    if args.analyze:
        report, rows = analyze(contract, manifest); freeze_results(report, rows)
        print(json.dumps({"status": report["status"], "classification": report["classification"],
                          "summary": report["summary"], "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
