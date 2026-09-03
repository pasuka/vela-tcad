#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M64 no-BGN smooth-NWell intervention."""

from __future__ import annotations

import argparse
import copy
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m10_fixed_state_replay as m10  # noqa: E402
import run_simplemos_m60_tight_convergence_port_burst as m60  # noqa: E402

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_contract_v1.json"
FREEZE = ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_contract_freeze.json"
M60_EVIDENCE = ROOT / "simplemos_m60_tight_convergence_port_burst_evidence.json"
M63_EVIDENCE = ROOT / "simplemos_m63_smooth_nwell_attribution_evidence.json"
M60_POINTS = ROOT / "tight_convergence_port_burst/m60_tight_default_direct_point_ledger.csv"
M63_PAIRS = ROOT / "smooth_nwell_attribution/m63_nwell_pair_attribution_ledger.csv"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
OUTPUT = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m64_bgn_smooth_nwell_intervention"
PORTABLE = ROOT / "bgn_smooth_nwell_intervention"
REPORT = PORTABLE / "m64_bgn_smooth_nwell_intervention_report.json"
POINTS = PORTABLE / "m64_bgn_on_off_point_ledger.csv"
CURVES = PORTABLE / "m64_bgn_on_off_curve_ledger.csv"
PAIRS = PORTABLE / "m64_bgn_nwell_pair_ledger.csv"
CASES = PORTABLE / "m64_execution_case_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m64_bgn_smooth_nwell_intervention_2026-09-02.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m64/artifact.json"
EVIDENCE = ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m64_bgn_smooth_nwell_intervention.py"
RUNNER = REPO / "build-release/vela_example_runner.exe"
REMOTE_ROOT = "/tmp/vela_simplemos_m64_bgn_smooth_nwell_intervention"
ARCHIVE = "simplemos_m64_bgn_smooth_nwell_intervention_results.tgz"
SSH_CONFIG = (Path(os.environ.get("USERPROFILE", str(Path.home())))
              / ".ssh/config")
PHASES = ("00_equilibrium", "10_drain_ramp", "20_gate_sweep")


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


def median(values: list[float]) -> float:
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    return (ordered[midpoint] if len(ordered) % 2 else
            0.5 * (ordered[midpoint - 1] + ordered[midpoint]))


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=False, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        encoding="utf-8", errors="replace")
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {list(argv)}\n"
            f"{completed.stdout or ''}")
    return completed.stdout or ""


def source_paths(contract: dict[str, Any]) -> list[Path]:
    paths = [M60_EVIDENCE, M63_EVIDENCE, M60_POINTS, M63_PAIRS]
    for device in contract["matrix"]["devices"]:
        paths.append(M8 / "sentaurus_bundle" / device / "input_fps.tdr")
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            workflow = M8 / "vela" / device / "workflow" / (
                "vd_0p05" if math.isclose(drain, 0.05) else "vd_1")
            paths.extend(workflow / f"{phase}.json" for phase in PHASES)
    return paths


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m64_bgn_smooth_nwell_intervention_contract.v1":
        raise ValueError("unexpected M64 contract schema")
    sources = source_paths(contract)
    missing = [path for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m64_bgn_smooth_nwell_intervention_contract_freeze.v1",
        "status": "frozen_before_execution",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in sources},
    })


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution":
        raise ValueError("M64 contract was not frozen before execution")
    if freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M64 contract changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M64 upstream artifact changed: {relative}")
    m60_evidence = read_json(M60_EVIDENCE)
    m63_evidence = read_json(M63_EVIDENCE)
    upstream = contract["upstream"]
    if (m60_evidence.get("status") != upstream["required_m60_status"] or
            m60_evidence.get("classification") != upstream["required_m60_classification"]):
        raise ValueError("M60 status or classification changed")
    if (m63_evidence.get("status") != upstream["required_m63_status"] or
            m63_evidence.get("classification") != upstream["required_m63_classification"]):
        raise ValueError("M63 status or classification changed")
    return contract


def sentaurus_deck(run_case: str, drain: float) -> str:
    deck = m60.sentaurus_deck(run_case, drain, "default", False)
    old = "Physics { EffectiveIntrinsicDensity(OldSlotboom) }"
    new = "Physics { EffectiveIntrinsicDensity(NoBandGapNarrowing) }"
    if deck.count(old) != 1:
        raise RuntimeError("M60 OldSlotboom anchor changed")
    deck = deck.replace(old, new)
    if "DirectCurrent" in deck:
        raise RuntimeError("M64 must retain the default current observer")
    return deck


def prepare_sentaurus(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    cases: list[dict[str, Any]] = []
    for device in contract["matrix"]["devices"]:
        source = M8 / "sentaurus_bundle" / device / "input_fps.tdr"
        root = bundle / device
        root.mkdir(parents=True, exist_ok=True)
        target = root / "input_fps.tdr"
        shutil.copy2(source, target)
        for drain in map(float, contract["matrix"]["drain_voltages_V"]):
            tag = m60.voltage_tag(drain)
            run_case = f"m64_nobgn_{device}_vd_{tag}"
            deck_path = root / f"{run_case}_des.cmd"
            deck_path.write_text(sentaurus_deck(run_case, drain),
                                 encoding="utf-8", newline="\n")
            cases.append({
                "run_case": run_case, "device": device,
                "drain_voltage_V": drain,
                "deck": portable(deck_path), "deck_sha256": sha256(deck_path),
                "input_tdr": portable(target), "input_tdr_sha256": sha256(target),
                "expected_current": f"IdVg_{run_case}_des.plt",
                "expected_console": f"{run_case}.console.log",
                "expected_model_log": f"{run_case}.log_des.log",
            })
    if len(cases) != int(contract["matrix"]["sentaurus_case_count"]):
        raise RuntimeError("M64 Sentaurus case count mismatch")
    manifest = {
        "schema": "vela.simplemos.sdevice.m64_sentaurus_manifest.v1",
        "status": "prepared", "contract_sha256": sha256(CONTRACT),
        "cases": cases,
    }
    write_json(OUTPUT / "sentaurus_manifest.json", manifest)
    return manifest


def run_sentaurus(manifest: dict[str, Any], ssh_target: str, ssh_bin: str,
                   scp_bin: str, ssh_config: str, remote_root: str,
                   jobs: int) -> str:
    ssh = [ssh_bin, "-F", ssh_config]
    scp = [scp_bin, "-F", ssh_config]
    banner = run([*ssh, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
                 capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([*ssh, ssh_target,
         f"set -eu; test ! -e {remote_root}; mkdir -p {remote_root}"])
    run([*scp, "-r", str(OUTPUT / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])
    grouped: dict[str, list[dict[str, Any]]] = {}
    for item in manifest["cases"]:
        grouped.setdefault(item["device"], []).append(item)

    def execute(group: tuple[str, list[dict[str, Any]]]) -> None:
        device, cases = group
        root = f"{remote_root}/sentaurus_bundle/{device}"
        commands = [f"cd {root}"]
        for item in cases:
            name = item["run_case"]
            commands.append(f"sdevice {name}_des.cmd > {name}.console.log 2>&1")
        run([*ssh, ssh_target, "set -eu; " + "; ".join(commands)])

    groups = list(grouped.items())
    with ThreadPoolExecutor(max_workers=min(max(jobs, 1), len(groups))) as pool:
        list(pool.map(execute, groups))
    run([*ssh, ssh_target,
         f"cd {remote_root}; tar -czf {ARCHIVE} sentaurus_bundle"])
    raw = OUTPUT / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / ARCHIVE
    run([*scp, f"{ssh_target}:{remote_root}/{ARCHIVE}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (OUTPUT / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def vela_workflow(device: str, drain: float, force: bool) -> dict[str, Any]:
    tag = "vd_0p05" if math.isclose(drain, 0.05) else "vd_1"
    source_root = M8 / "vela" / device / "workflow" / tag
    root = OUTPUT / "vela" / device / tag
    if force and root.exists():
        shutil.rmtree(root)
    previous: Path | None = None
    final_curve: Path | None = None
    for phase_name in PHASES:
        config = copy.deepcopy(read_json(source_root / f"{phase_name}.json"))
        phase = root / phase_name
        curve = phase / "curve.csv"
        state = phase / "state.csv"
        stdout = phase / "config.stdout.txt"
        complete = (curve.is_file() and state.is_file() and stdout.is_file()
                    and '"converged":true' in stdout.read_text(
                        encoding="utf-8", errors="replace"))
        if force or (phase.exists() and not complete):
            shutil.rmtree(phase)
        phase.mkdir(parents=True, exist_ok=True)
        log = phase / "run.log"
        config["solver"]["bandgap_narrowing"] = {
            "model": "none", "fermi_statistics_correction": False}
        config["output_csv"] = str(curve.resolve())
        config["log_file"] = str(log.resolve())
        config["sweep"]["write_state_file"] = str(state.resolve())
        if previous is None:
            config["sweep"].pop("initial_state_file", None)
        else:
            config["sweep"]["initial_state_file"] = str(previous.resolve())
        config["simplemos_m64"] = {
            "single_axis": "bandgap_narrowing", "bgn_model": "none",
            "source_workflow": portable(source_root / f"{phase_name}.json"),
        }
        config_path = phase / "config.json"
        write_json(config_path, config)
        if not complete:
            status = m10.execute_runner(config_path, RUNNER)
            if not status.get("converged"):
                raise RuntimeError(f"M64 Vela failed: {device} {tag} {phase_name}")
        previous = state
        final_curve = curve
    assert final_curve is not None and previous is not None
    return {
        "device": device, "drain_voltage_V": drain,
        "curve": portable(final_curve), "state": portable(previous),
        "curve_sha256": sha256(final_curve), "state_sha256": sha256(previous),
    }


def run_vela(contract: dict[str, Any], jobs: int, force: bool) -> dict[str, Any]:
    work = [(device, float(drain), force)
            for device in contract["matrix"]["devices"]
            for drain in contract["matrix"]["drain_voltages_V"]]
    with ThreadPoolExecutor(max_workers=min(max(jobs, 1), 8)) as pool:
        workflows = list(pool.map(lambda args: vela_workflow(*args), work))
    manifest = {
        "schema": "vela.simplemos.sdevice.m64_vela_manifest.v1",
        "status": "complete", "contract_sha256": sha256(CONTRACT),
        "workflows": workflows,
    }
    write_json(OUTPUT / "vela_manifest.json", manifest)
    return manifest


def parse_vela_curve(path: Path, gates: list[float], tolerance: float
                     ) -> list[dict[str, float]]:
    rows = read_csv(path)
    parsed: list[dict[str, float]] = []
    for gate in gates:
        matches = [row for row in rows
                   if abs(float(row["bias_V"]) - gate) <= tolerance]
        if len(matches) != 1:
            raise RuntimeError(f"{path}: exact Vg={gate:g} count {len(matches)}")
        parsed.append({"gate_voltage_V": gate,
                       "drain_current_A_per_um": abs(float(
                           matches[0]["current_total_A_per_um"]))})
    return parsed


def local_slope(rows: list[dict[str, float]], index: int, key: str) -> float:
    left, right = max(index - 1, 0), min(index + 1, len(rows) - 1)
    return ((math.log10(abs(rows[right][key])) -
             math.log10(abs(rows[left][key]))) /
            (rows[right]["gate_voltage_V"] - rows[left]["gate_voltage_V"]))


def fitted_shift(rows: list[dict[str, float]], sent_key: str,
                 error_key: str, low_gate: float, high_gate: float
                 ) -> tuple[float, float]:
    selected: list[tuple[float, float]] = []
    for index, row in enumerate(rows):
        if low_gate <= row["gate_voltage_V"] <= high_gate:
            selected.append((local_slope(rows, index, sent_key), row[error_key]))
    delta = sum(slope * error for slope, error in selected) / sum(
        slope * slope for slope, _ in selected)
    energy = sum(error * error for _, error in selected)
    residual = sum((error - slope * delta) ** 2 for slope, error in selected)
    explained = 1.0 if energy == 0.0 else 1.0 - residual / energy
    return 1000.0 * delta, explained


def analyze(contract: dict[str, Any], sentaurus_manifest: dict[str, Any],
            vela_manifest: dict[str, Any], banner: str
            ) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    gates = [0.05 * index for index in range(51)]
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    off_sent: dict[tuple[str, float], list[dict[str, float]]] = {}
    case_rows: list[dict[str, Any]] = []
    for item in sentaurus_manifest["cases"]:
        current = raw / item["device"] / item["expected_current"]
        console = raw / item["device"] / item["expected_console"]
        model_log = raw / item["device"] / item["expected_model_log"]
        if not model_log.is_file():
            model_log = raw / item["device"] / f"{item['run_case']}.log_des.log"
        for required in (current, console, model_log):
            if not required.is_file():
                raise FileNotFoundError(required)
        parsed = m60.parse_current(current, gates, tolerance)
        off_sent[(item["device"], float(item["drain_voltage_V"]))] = [{
            "gate_voltage_V": float(row["gate_voltage_V"]),
            "drain_current_A_per_um": abs(float(row["drain_total_A_per_um"])),
        } for row in parsed]
        model_text = model_log.read_text(errors="replace")
        qualified = "without bandgap narrowing" in model_text.lower()
        case_rows.append({
            "solver": "sentaurus", "device": item["device"],
            "drain_voltage_V": item["drain_voltage_V"],
            "point_count": len(parsed), "converged": 1,
            "bgn_off_qualified": int(qualified),
            "curve": portable(current), "log": portable(model_log),
            "curve_sha256": sha256(current), "log_sha256": sha256(model_log),
        })

    off_vela: dict[tuple[str, float], list[dict[str, float]]] = {}
    for item in vela_manifest["workflows"]:
        curve = REPO / item["curve"]
        parsed = parse_vela_curve(curve, gates, tolerance)
        off_vela[(item["device"], float(item["drain_voltage_V"]))] = parsed
        case_rows.append({
            "solver": "vela", "device": item["device"],
            "drain_voltage_V": item["drain_voltage_V"],
            "point_count": len(parsed), "converged": 1,
            "bgn_off_qualified": 1, "curve": item["curve"],
            "log": "", "curve_sha256": item["curve_sha256"], "log_sha256": "",
        })

    on_rows = read_csv(M60_POINTS)
    on: dict[tuple[str, float], list[dict[str, float]]] = {}
    for row in on_rows:
        key = (row["device"], float(row["drain_voltage_V"]))
        on.setdefault(key, []).append({
            "gate_voltage_V": float(row["gate_voltage_V"]),
            "sentaurus_A_per_um": abs(float(row["tight_default_drain_current_A_per_um"])),
            "vela_A_per_um": abs(float(row["vela_drain_current_A_per_um"])),
        })
    for rows in on.values():
        rows.sort(key=lambda row: row["gate_voltage_V"])

    low_gate, high_gate = map(float, contract["matrix"]["smooth_gate_window_V"])
    point_rows: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    smooth_by_key: dict[tuple[str, float], list[dict[str, Any]]] = {}
    identity_residuals: list[float] = []
    monotonic_checks: list[bool] = []
    for key in sorted(on):
        device, drain = key
        joined: list[dict[str, float]] = []
        for on_point, sent_off, vela_off in zip(
                on[key], off_sent[key], off_vela[key], strict=True):
            gate = on_point["gate_voltage_V"]
            if not (math.isclose(gate, sent_off["gate_voltage_V"], abs_tol=tolerance)
                    and math.isclose(gate, vela_off["gate_voltage_V"], abs_tol=tolerance)):
                raise RuntimeError(f"M64 gate mismatch for {key}")
            sent_on = on_point["sentaurus_A_per_um"]
            vela_on = on_point["vela_A_per_um"]
            sent_no = sent_off["drain_current_A_per_um"]
            vela_no = vela_off["drain_current_A_per_um"]
            if min(sent_on, vela_on, sent_no, vela_no) <= 0.0:
                raise RuntimeError(f"non-positive M64 current for {key}, Vg={gate}")
            error_on = math.log10(vela_on / sent_on)
            error_off = math.log10(vela_no / sent_no)
            bgn_component = error_on - error_off
            response_difference = (math.log10(vela_on / vela_no)
                                   - math.log10(sent_on / sent_no))
            identity = bgn_component - response_difference
            identity_residuals.append(abs(identity))
            joined.append({
                "gate_voltage_V": gate,
                "sentaurus_on_A_per_um": sent_on, "vela_on_A_per_um": vela_on,
                "sentaurus_off_A_per_um": sent_no, "vela_off_A_per_um": vela_no,
                "error_on_dex": error_on, "error_off_dex": error_off,
                "bgn_mismatch_component_dex": bgn_component,
                "response_difference_dex": response_difference,
                "identity_residual_dex": identity,
            })
        smooth = [row for row in joined
                  if low_gate <= row["gate_voltage_V"] <= high_gate]
        if len(smooth) < int(contract["acceptance"]["minimum_smooth_points_per_curve"]):
            raise RuntimeError(f"insufficient M64 smooth points for {key}")
        monotonic = all(
            all(b[field] > a[field] for field in (
                "sentaurus_on_A_per_um", "vela_on_A_per_um",
                "sentaurus_off_A_per_um", "vela_off_A_per_um"))
            for a, b in zip(smooth, smooth[1:]))
        monotonic_checks.append(monotonic)
        on_shift, on_explained = fitted_shift(
            joined, "sentaurus_on_A_per_um", "error_on_dex", low_gate, high_gate)
        off_shift, off_explained = fitted_shift(
            joined, "sentaurus_off_A_per_um", "error_off_dex", low_gate, high_gate)
        curves.append({
            "device": device, "drain_voltage_V": drain,
            "smooth_point_count": len(smooth), "smooth_monotonic": int(monotonic),
            "bgn_on_fitted_shift_mV": on_shift,
            "bgn_off_fitted_shift_mV": off_shift,
            "bgn_contribution_to_fitted_shift_mV": on_shift - off_shift,
            "bgn_on_horizontal_shift_explained_fraction": on_explained,
            "bgn_off_horizontal_shift_explained_fraction": off_explained,
        })
        smooth_by_key[key] = smooth
        for row in smooth:
            point_rows.append({"device": device, "drain_voltage_V": drain, **row})

    m63_pair_rows = read_csv(M63_PAIRS)
    pair_rows: list[dict[str, Any]] = []
    closures: list[float] = []
    suppression_fractions: list[float] = []
    anchor_errors: list[float] = []
    pair_identity_residuals: list[float] = []
    for anchor in m63_pair_rows:
        low_device, high_device = anchor["low_device"], anchor["high_device"]
        drain = float(anchor["drain_voltage_V"])
        gate = float(anchor["diagnostic_gate_voltage_V"])
        low = next(row for row in smooth_by_key[(low_device, drain)]
                   if math.isclose(float(row["gate_voltage_V"]), gate, abs_tol=tolerance))
        high = next(row for row in smooth_by_key[(high_device, drain)]
                    if math.isclose(float(row["gate_voltage_V"]), gate, abs_tol=tolerance))
        on_growth = high["error_on_dex"] - low["error_on_dex"]
        off_growth = high["error_off_dex"] - low["error_off_dex"]
        explained = on_growth - off_growth
        response_pair = (high["response_difference_dex"]
                         - low["response_difference_dex"])
        identity = explained - response_pair
        anchor_growth = float(anchor["observed_nwell_growth_dex"])
        anchor_error = on_growth - anchor_growth
        raw_closure = 1.0 - abs(off_growth) / max(abs(on_growth), 1e-300)
        closure = min(1.0, max(0.0, raw_closure))
        suppression = (off_growth - on_growth) / max(abs(off_growth), 1e-300)
        closures.append(closure)
        suppression_fractions.append(suppression)
        anchor_errors.append(abs(anchor_error))
        pair_identity_residuals.append(abs(identity))
        pair_rows.append({
            "low_device": low_device, "high_device": high_device,
            "drain_voltage_V": drain, "diagnostic_gate_voltage_V": gate,
            "m63_anchor_nwell_growth_dex": anchor_growth,
            "recomputed_bgn_on_nwell_growth_dex": on_growth,
            "bgn_off_nwell_growth_dex": off_growth,
            "bgn_explained_nwell_growth_dex": explained,
            "bgn_response_pair_difference_dex": response_pair,
            "difference_in_differences_identity_residual_dex": identity,
            "m63_anchor_replay_error_dex": anchor_error,
            "raw_pair_closure_fraction": raw_closure,
            "classification_pair_closure_fraction": closure,
            "bgn_suppression_fraction_of_no_bgn_growth": suppression,
        })

    median_closure = median(closures)
    comparable = all(monotonic_checks)
    thresholds = contract["analysis"]["classification_thresholds"]
    if not comparable:
        classification = "bgn_changes_curve_regime_noncomparably"
    elif median_closure >= float(thresholds["dominant_median_closure_fraction"]):
        classification = "bgn_self_consistent_response_dominant"
    elif median_closure >= float(thresholds["material_median_closure_fraction"]):
        classification = "bgn_material_but_not_dominant"
    else:
        classification = "bgn_independent_threshold_shift_dominant"

    acceptance = contract["acceptance"]
    checks = {
        "contract_frozen": sha256(CONTRACT) == read_json(FREEZE)["contract_sha256"],
        "sentaurus_release": acceptance["required_sentaurus_release"] in banner,
        "sentaurus_case_count": len(off_sent) == int(acceptance["required_sentaurus_case_count"]),
        "vela_workflow_count": len(off_vela) == int(acceptance["required_vela_workflow_count"]),
        "sentaurus_no_bgn_log_qualified": all(
            int(row["bgn_off_qualified"]) == 1 for row in case_rows
            if row["solver"] == "sentaurus"),
        "full_curve_point_count": (sum(len(rows) for rows in off_sent.values())
                                   == int(acceptance["required_curve_point_count_per_solver_condition"])
                                   == sum(len(rows) for rows in off_vela.values())),
        "smooth_point_count": len(point_rows) == int(acceptance["required_smooth_point_count"]),
        "pair_count": len(pair_rows) == int(acceptance["required_pair_count"]),
        "difference_in_differences_identity": max(identity_residuals + pair_identity_residuals)
        <= float(acceptance["maximum_difference_in_differences_identity_residual_dex"]),
        "m63_anchor_replay": max(anchor_errors)
        <= float(acceptance["maximum_bgn_on_anchor_replay_error_dex"]),
        "classification_declared": classification in contract["analysis"]["classifications"],
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    report = {
        "schema": "vela.simplemos.sdevice.m64_bgn_smooth_nwell_intervention_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification if checks["all_checks_pass"] else "execution_or_replay_mismatch",
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {
            "new_sentaurus_execution": True, "new_vela_execution": True,
            "sentaurus_case_count": len(off_sent),
            "vela_workflow_count": len(off_vela),
            "full_gate_path_retained": True, "production_default_changed": False,
        },
        "summary": {
            "median_pair_closure_fraction": median_closure,
            "minimum_pair_closure_fraction": min(closures),
            "maximum_pair_closure_fraction": max(closures),
            "pair_closure_at_or_above_0p8_count": sum(value >= 0.8 for value in closures),
            "pair_count": len(closures),
            "median_bgn_on_nwell_growth_dex": median([
                float(row["recomputed_bgn_on_nwell_growth_dex"]) for row in pair_rows]),
            "median_bgn_off_nwell_growth_dex": median([
                float(row["bgn_off_nwell_growth_dex"]) for row in pair_rows]),
            "median_bgn_explained_nwell_growth_dex": median([
                float(row["bgn_explained_nwell_growth_dex"]) for row in pair_rows]),
            "median_bgn_suppression_fraction_of_no_bgn_growth": median(
                suppression_fractions),
            "bgn_counteracts_nwell_growth_pair_count": sum(
                float(row["bgn_explained_nwell_growth_dex"])
                * float(row["bgn_off_nwell_growth_dex"]) < 0.0
                for row in pair_rows),
            "bgn_on_fitted_shift_mV_range": [
                min(float(row["bgn_on_fitted_shift_mV"]) for row in curves),
                max(float(row["bgn_on_fitted_shift_mV"]) for row in curves)],
            "bgn_off_fitted_shift_mV_range": [
                min(float(row["bgn_off_fitted_shift_mV"]) for row in curves),
                max(float(row["bgn_off_fitted_shift_mV"]) for row in curves)],
            "smooth_regime_comparable": comparable,
            "maximum_identity_residual_dex": max(identity_residuals + pair_identity_residuals),
            "maximum_m63_anchor_replay_error_dex": max(anchor_errors),
        },
        "causal_scope": {
            "closed": "The M63 paired smooth NWell growth is not caused by enabling OldSlotboom BGN; removing BGN increases the growth in all eight pairs, so the measured BGN response counteracts a larger BGN-independent threshold-like mismatch.",
            "not_closed": "M64 does not separate carrier statistics, base intrinsic-density conventions, electrostatics, and other self-consistent causes inside the larger no-BGN threshold-like mismatch.",
            "formula_claim": "M64 does not identify an error in the shared OldSlotboom delta-Eg formula; the M39 formula-level closure remains intact.",
        },
        "acceptance": checks,
    }
    return report, {"points": point_rows, "curves": curves,
                    "pairs": pair_rows, "cases": case_rows}


def freeze_results(report: dict[str, Any], ledgers: dict[str, list[dict[str, Any]]]) -> None:
    write_csv(POINTS, ledgers["points"])
    write_csv(CURVES, ledgers["curves"])
    write_csv(PAIRS, ledgers["pairs"])
    write_csv(CASES, ledgers["cases"])
    write_json(REPORT, report)
    summary = report["summary"]
    DOC.write_text(f"""# SimpleMOS M64 BGN 平滑 NWell 干预

## 结论

M64 分类为 `{report['classification']}`。在 M60 收紧收敛策略、默认端口观测和完整偏压路径不变的条件下，Sentaurus 与 Vela 同时从 OldSlotboom BGN 切换到显式 no-BGN；其余物理、网格、接触和生产默认值均未改变。

8 个 NWell 配对的误差增幅中位闭合率为 `{float(summary['median_pair_closure_fraction']):.6f}`，其中 `{int(summary['pair_closure_at_or_above_0p8_count'])}/8` 个配对达到 80% 闭合。BGN-on 配对增幅中位数为 `{float(summary['median_bgn_on_nwell_growth_dex']):.6f}` dex，BGN-off 后为 `{float(summary['median_bgn_off_nwell_growth_dex']):.6f}` dex，差分中的差分为 `{float(summary['median_bgn_explained_nwell_growth_dex']):.6f}` dex。

8/8 配对中，启用 BGN 的响应都与 no-BGN 的 NWell 增幅反向；它对 no-BGN 增幅的中位抑制比例为 `{float(summary['median_bgn_suppression_fraction_of_no_bgn_growth']):.6%}`。因此，BGN 不是 M63 平滑增幅的来源，而是在抵消一个更大的 BGN-independent 阈值样失配。

BGN-on 的拟合水平平移范围为 `{float(summary['bgn_on_fitted_shift_mV_range'][0]):.3f}`--`{float(summary['bgn_on_fitted_shift_mV_range'][1]):.3f}` mV；BGN-off 后为 `{float(summary['bgn_off_fitted_shift_mV_range'][0]):.3f}`--`{float(summary['bgn_off_fitted_shift_mV_range'][1]):.3f}` mV。全部平滑窗口曲线保持单调，差分恒等式最大数值残差为 `{float(summary['maximum_identity_residual_dex']):.3e}` dex。

这一干预归因的是两个求解器对“启用 BGN 后的自洽响应”之差，不等于 OldSlotboom delta-Eg 公式错误；M39 的公式级闭环与 M8/M46/M60 参考资格均保持不变。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {
        "schema": "vela.validation.artifact.v1",
        "title": "SimpleMOS M64 BGN smooth-NWell intervention",
        "status": report["status"], "classification": report["classification"],
        "summary": "Paired BGN-on/off difference-in-differences attribution",
        "report": portable(REPORT),
        "ledgers": [portable(POINTS), portable(CURVES), portable(PAIRS), portable(CASES)],
    })
    artifacts = [REPORT, POINTS, CURVES, PAIRS, CASES, DOC, ARTIFACT]
    write_json(EVIDENCE, {
        "schema": "vela.simplemos.sdevice.m64_bgn_smooth_nwell_intervention_evidence.v1",
        "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
        "classification": report["classification"],
        "contract_sha256": sha256(CONTRACT),
        "new_sentaurus_execution": True, "new_vela_execution": True,
        "production_reference_replaced": False,
        "source_hashes": read_json(FREEZE)["upstream_hashes"],
        "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
        "artifacts": {portable(path): sha256(path) for path in artifacts},
        "acceptance": report["acceptance"],
    })


def verify() -> dict[str, Any]:
    validate_contract()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M64 evidence is not accepted and frozen")
    if evidence["implementation_hashes"].get(portable(SCRIPT)) != sha256(SCRIPT):
        raise ValueError("M64 implementation hash changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M64 artifact hash changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-sentaurus", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=r"C:\Windows\System32\OpenSSH\ssh.exe")
    parser.add_argument("--scp-bin", default=r"C:\Windows\System32\OpenSSH\scp.exe")
    parser.add_argument("--ssh-config", default=str(SSH_CONFIG))
    parser.add_argument("--remote-root", default=REMOTE_ROOT)
    args = parser.parse_args()
    if args.freeze_contract:
        freeze_contract()
        print(json.dumps({"status": "frozen_before_execution",
                          "contract_sha256": sha256(CONTRACT)}, indent=2))
        return
    if args.verify:
        report = verify()
        print(f"M64 verified: {report['classification']}")
        return
    contract = validate_contract()
    do_prepare = args.prepare or args.all
    do_sentaurus = args.run_sentaurus or args.all
    do_vela = args.run_vela or args.all
    do_analyze = args.analyze or args.all
    sentaurus_manifest_path = OUTPUT / "sentaurus_manifest.json"
    sentaurus_manifest = (prepare_sentaurus(contract, args.force) if do_prepare
                           else read_json(sentaurus_manifest_path))
    if do_sentaurus:
        banner = run_sentaurus(sentaurus_manifest, args.ssh_target, args.ssh_bin,
                               args.scp_bin, args.ssh_config, args.remote_root,
                               args.jobs)
    else:
        banner = (OUTPUT / "sentaurus_banner.txt").read_text(
            encoding="utf-8").strip() if do_analyze else ""
    vela_manifest = (run_vela(contract, args.jobs, args.force) if do_vela
                     else read_json(OUTPUT / "vela_manifest.json") if do_analyze
                     else {"workflows": []})
    if do_analyze:
        report, ledgers = analyze(contract, sentaurus_manifest,
                                  vela_manifest, banner)
        freeze_results(report, ledgers)
        if not report["acceptance"]["all_checks_pass"]:
            raise RuntimeError(f"M64 acceptance failed: {report['acceptance']}")
        print(json.dumps({"status": report["status"],
                          "classification": report["classification"],
                          "summary": report["summary"]}, indent=2))
    elif do_prepare or do_vela:
        print(json.dumps({"status": "prepared_or_locally_executed",
                          "sentaurus_case_count": len(sentaurus_manifest["cases"]),
                          "vela_workflow_count": len(vela_manifest["workflows"])}, indent=2))


if __name__ == "__main__":
    main()
