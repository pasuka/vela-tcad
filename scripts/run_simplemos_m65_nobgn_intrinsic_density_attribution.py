#!/usr/bin/env python3
"""Execute and freeze SimpleMOS M65 no-BGN intrinsic-density attribution."""

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
import run_simplemos_m47_default_bgn_self_consistent_attribution as m47  # noqa: E402
import run_simplemos_m60_tight_convergence_port_burst as m60  # noqa: E402

ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
CONTRACT = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_contract_freeze.json"
M64_EVIDENCE = ROOT / "simplemos_m64_bgn_smooth_nwell_intervention_evidence.json"
M64_REPORT = ROOT / "bgn_smooth_nwell_intervention/m64_bgn_smooth_nwell_intervention_report.json"
M64_POINTS = ROOT / "bgn_smooth_nwell_intervention/m64_bgn_on_off_point_ledger.csv"
M64_PAIRS = ROOT / "bgn_smooth_nwell_intervention/m64_bgn_nwell_pair_ledger.csv"
M40_NI = REPO / ("build-release/reference_tcad/simplemos_sentaurus2022/"
                 "m40_true_no_bgn_factorial/sentaurus_exports/bgn_off_srh_on/"
                 "fields/EffectiveIntrinsicDensity_region0.csv")
MATERIALS = REPO / "reference_tcad/transportmodels_sentaurus2022/vela/materials_sentaurus2022.json"
M8 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics"
# Keep this deliberately short: Sentaurus snapshot sidecars can otherwise cross
# the legacy Windows MAX_PATH boundary during tar extraction.
OUTPUT = REPO / "build-release/m65_ni"
PORTABLE = ROOT / "nobgn_intrinsic_density_attribution"
REPORT = PORTABLE / "m65_nobgn_intrinsic_density_attribution_report.json"
CURRENT = PORTABLE / "m65_current_intervention_ledger.csv"
STATES = PORTABLE / "m65_state_attribution_ledger.csv"
PAIRS = PORTABLE / "m65_nwell_pair_closure_ledger.csv"
CASES = PORTABLE / "m65_execution_case_ledger.csv"
DOC = REPO / "docs/validation/simplemos_m65_nobgn_intrinsic_density_attribution_2026-09-03.md"
ARTIFACT = REPO / "docs/validation/reports/simplemos_m65/artifact.json"
EVIDENCE = ROOT / "simplemos_m65_nobgn_intrinsic_density_attribution_evidence.json"
SCRIPT = REPO / "scripts/run_simplemos_m65_nobgn_intrinsic_density_attribution.py"
RUNNER = REPO / "build-release/vela_example_runner.exe"
IMPORTER = REPO / "build-release/sentaurus_import.exe"
REMOTE_ROOT = "/tmp/vela_simplemos_m65_nobgn_intrinsic_density_attribution"
ARCHIVE = "simplemos_m65_nobgn_intrinsic_density_attribution_results.tgz"
SSH_CONFIG = Path(os.environ.get("USERPROFILE", str(Path.home()))) / ".ssh/config"
PHASES = ("00_equilibrium", "10_drain_ramp", "20_gate_sweep")
VT = 8.617333262145e-5 * 300.0


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
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable(path: Path) -> str:
    return path.resolve().relative_to(REPO.resolve()).as_posix()


def median(values: list[float]) -> float:
    values = sorted(values)
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else 0.5 * (values[middle - 1] + values[middle])


def percentile(values: list[float], fraction: float) -> float:
    values = sorted(values)
    if not values:
        return 0.0
    position = fraction * (len(values) - 1)
    left = int(math.floor(position))
    right = int(math.ceil(position))
    if left == right:
        return values[left]
    return values[left] + (values[right] - values[left]) * (position - left)


def voltage_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=False, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
        encoding="utf-8", errors="replace")
    if completed.returncode:
        raise RuntimeError(f"command failed ({completed.returncode}): {list(argv)}\n{completed.stdout or ''}")
    return completed.stdout or ""


def cases(contract: dict[str, Any]) -> list[dict[str, Any]]:
    result: dict[tuple[str, float], dict[str, Any]] = {}
    for pair in contract["matrix"]["pairs"]:
        drain = float(pair["drain_voltage_V"])
        gate = float(pair["diagnostic_gate_voltage_V"])
        for device in (pair["low_device"], pair["high_device"]):
            key = (device, drain)
            if key in result and not math.isclose(result[key]["gate_voltage_V"], gate):
                raise ValueError(f"ambiguous M65 gate for {key}")
            result[key] = {"device": device, "drain_voltage_V": drain,
                           "gate_voltage_V": gate}
    return [result[key] for key in sorted(result)]


def source_paths(contract: dict[str, Any]) -> list[Path]:
    paths = [M64_EVIDENCE, M64_REPORT, M64_POINTS, M64_PAIRS, M40_NI, MATERIALS]
    for item in cases(contract):
        device = item["device"]
        drain = float(item["drain_voltage_V"])
        tag = "vd_0p05" if math.isclose(drain, 0.05) else "vd_1"
        paths.append(M8 / "sentaurus_bundle" / device / "input_fps.tdr")
        paths.extend(M8 / "vela" / device / "workflow" / tag / f"{phase}.json"
                     for phase in PHASES)
    return paths


def validate_ni_anchor(contract: dict[str, Any]) -> None:
    values = [float(row["component0"]) for row in read_csv(M40_NI)]
    expected = float(contract["intervention"]["matched_vela_ni_cm3"])
    if len(values) != 942 or min(values) != max(values) or not math.isclose(values[0], expected, rel_tol=1e-15):
        raise ValueError("M65 Sentaurus no-BGN intrinsic-density anchor changed")
    materials = read_json(MATERIALS)
    silicon = next(item for item in materials["materials"] if item["name"] == "Si")
    if not math.isclose(float(silicon["ni"]),
                        float(contract["intervention"]["baseline_vela_ni_cm3"]), rel_tol=1e-15):
        raise ValueError("M65 Vela baseline Si ni changed")


def freeze_contract() -> None:
    contract = read_json(CONTRACT)
    if contract.get("schema") != "vela.simplemos.sdevice.m65_nobgn_intrinsic_density_attribution_contract.v1":
        raise ValueError("unexpected M65 contract schema")
    validate_ni_anchor(contract)
    sources = source_paths(contract)
    missing = [path for path in sources if not path.is_file()]
    if missing:
        raise FileNotFoundError(missing[0])
    write_json(FREEZE, {
        "schema": "vela.simplemos.sdevice.m65_nobgn_intrinsic_density_attribution_contract_freeze.v1",
        "status": "frozen_before_execution",
        "contract": portable(CONTRACT),
        "contract_sha256": sha256(CONTRACT),
        "upstream_hashes": {portable(path): sha256(path) for path in sources},
    })


def validate_contract() -> dict[str, Any]:
    contract = read_json(CONTRACT)
    freeze = read_json(FREEZE)
    if freeze.get("status") != "frozen_before_execution" or freeze.get("contract_sha256") != sha256(CONTRACT):
        raise ValueError("M65 contract is not frozen or changed after freeze")
    for relative, expected in freeze["upstream_hashes"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M65 upstream artifact changed: {relative}")
    evidence = read_json(M64_EVIDENCE)
    required = contract["upstream"]
    if (evidence.get("status") != required["required_m64_status"] or
            evidence.get("classification") != required["required_m64_classification"]):
        raise ValueError("M64 status or classification changed")
    validate_ni_anchor(contract)
    return contract


def sentaurus_deck(run_case: str, drain: float, gate: float) -> str:
    deck = m60.sentaurus_deck(run_case, drain, "default", False)
    old = "Physics { EffectiveIntrinsicDensity(OldSlotboom) }"
    if deck.count(old) != 1:
        raise RuntimeError("M60 BGN anchor changed")
    deck = deck.replace(old, "Physics { EffectiveIntrinsicDensity(NoBandGapNarrowing) }")
    goal = 'Goal { Name="gate" Voltage=2.5 }'
    if deck.count(goal) != 1:
        raise RuntimeError("M60 gate goal anchor changed")
    deck = deck.replace(goal, f'Goal {{ Name="gate" Voltage={gate:.12g} }}')
    # Sentaurus Quasistationary steps are normalized to the Goal delta.  Scale
    # the normalized values so truncating the 2.5 V M64 path does not change
    # its physical 0.025 V initial, 2.5e-5 V minimum, or 0.125 V maximum step.
    gate_steps = "DoZero InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05"
    scaled_steps = (f"DoZero InitialStep={0.025 / gate:.17g} Increment=1.5 "
                    f"MinStep={2.5e-5 / gate:.17g} MaxStep={0.125 / gate:.17g}")
    if deck.count(gate_steps) != 1:
        raise RuntimeError("M60 normalized gate-step anchor changed")
    deck = deck.replace(gate_steps, scaled_steps)
    current_plot = "CurrentPlot(Time=(Range=(0 1) Intervals=50))"
    replacement = (f"CurrentPlot(Time=(1))\n    Plot(FilePrefix=\"{run_case}_state\" "
                   "NoOverWrite Time=(1))")
    if deck.count(current_plot) != 1:
        raise RuntimeError("M60 CurrentPlot anchor changed")
    deck = deck.replace(current_plot, replacement)
    deck = deck.replace("BandGap BandGapNarrowing Affinity",
                        "BandGap BandGapNarrowing Affinity EffectiveIntrinsicDensity")
    if "DirectCurrent" in deck:
        raise RuntimeError("M65 must retain default Sentaurus current observer")
    return deck


def prepare_sentaurus(contract: dict[str, Any], force: bool) -> dict[str, Any]:
    bundle = OUTPUT / "sentaurus_bundle"
    if force and bundle.exists():
        shutil.rmtree(bundle)
    prepared = []
    for item in cases(contract):
        device = item["device"]
        drain = float(item["drain_voltage_V"])
        gate = float(item["gate_voltage_V"])
        root = bundle / device / f"vd_{voltage_tag(drain)}"
        root.mkdir(parents=True, exist_ok=True)
        source = M8 / "sentaurus_bundle" / device / "input_fps.tdr"
        target = root / "input_fps.tdr"
        shutil.copy2(source, target)
        run_case = f"m65_nobgn_{device}_vd_{voltage_tag(drain)}_vg_{voltage_tag(gate)}"
        deck = root / f"{run_case}_des.cmd"
        deck.write_text(sentaurus_deck(run_case, drain, gate), encoding="utf-8", newline="\n")
        prepared.append({**item, "run_case": run_case, "deck": portable(deck),
                         "deck_sha256": sha256(deck), "input_tdr": portable(target),
                         "input_tdr_sha256": sha256(target)})
    if len(prepared) != int(contract["matrix"]["sentaurus_case_count"]):
        raise RuntimeError("M65 Sentaurus case count mismatch")
    manifest = {"schema": "vela.simplemos.sdevice.m65_sentaurus_manifest.v1",
                "status": "prepared", "contract_sha256": sha256(CONTRACT), "cases": prepared}
    write_json(OUTPUT / "sentaurus_manifest.json", manifest)
    return manifest


def run_sentaurus(manifest: dict[str, Any], ssh_target: str, ssh_bin: str,
                   scp_bin: str, ssh_config: str, remote_root: str, jobs: int) -> str:
    ssh = [ssh_bin, "-F", ssh_config]
    scp = [scp_bin, "-F", ssh_config]
    banner = run([*ssh, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"], capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([*ssh, ssh_target, f"set -eu; test ! -e {remote_root}; mkdir -p {remote_root}"])
    run([*scp, "-r", str(OUTPUT / "sentaurus_bundle"), f"{ssh_target}:{remote_root}/"])

    def execute(item: dict[str, Any]) -> None:
        relative = (REPO / item["deck"]).parent.relative_to(OUTPUT / "sentaurus_bundle").as_posix()
        name = item["run_case"]
        root = f"{remote_root}/sentaurus_bundle/{relative}"
        run([*ssh, ssh_target, f"set -eu; cd {root}; sdevice {name}_des.cmd > {name}.console.log 2>&1"])

    with ThreadPoolExecutor(max_workers=min(max(1, jobs), len(manifest["cases"]))) as pool:
        list(pool.map(execute, manifest["cases"]))
    run([*ssh, ssh_target, f"cd {remote_root}; tar -czf {ARCHIVE} sentaurus_bundle"])
    raw = OUTPUT / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / ARCHIVE
    run([*scp, f"{ssh_target}:{remote_root}/{ARCHIVE}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (OUTPUT / "sentaurus_banner.txt").write_text(banner + "\n", encoding="utf-8", newline="\n")
    return banner


def matched_materials(contract: dict[str, Any]) -> Path:
    target = OUTPUT / "vela/matched_materials.json"
    materials = copy.deepcopy(read_json(MATERIALS))
    changed = []
    for item in materials["materials"]:
        if item["name"] == "Si":
            item["ni"] = float(contract["intervention"]["matched_vela_ni_cm3"])
            changed.append(item["name"])
    if changed != ["Si"]:
        raise RuntimeError(f"unexpected M65 material changes: {changed}")
    materials["simplemos_m65"] = {"single_axis": "Si.ni",
                                  "source": portable(MATERIALS),
                                  "production_file_changed": False}
    write_json(target, materials)
    return target


def vela_workflow(item: dict[str, Any], contract: dict[str, Any], materials: Path,
                  force: bool) -> dict[str, Any]:
    device = item["device"]
    drain = float(item["drain_voltage_V"])
    gate = float(item["gate_voltage_V"])
    tag = "vd_0p05" if math.isclose(drain, 0.05) else "vd_1"
    source_root = M8 / "vela" / device / "workflow" / tag
    root = OUTPUT / "vela" / device / f"vd_{voltage_tag(drain)}"
    if force and root.exists():
        shutil.rmtree(root)
    previous: Path | None = None
    final_curve: Path | None = None
    final_state: Path | None = None
    for phase_name in PHASES:
        phase = root / phase_name
        config = copy.deepcopy(read_json(source_root / f"{phase_name}.json"))
        curve, state = phase / "curve.csv", phase / "state.csv"
        stdout = phase / "config.stdout.txt"
        complete = (curve.is_file() and state.is_file() and stdout.is_file() and
                    '"converged":true' in stdout.read_text(encoding="utf-8", errors="replace"))
        if phase.exists() and (force or not complete):
            shutil.rmtree(phase)
        phase.mkdir(parents=True, exist_ok=True)
        config["materials_file"] = str(materials.resolve())
        config["solver"]["bandgap_narrowing"] = {"model": "none", "fermi_statistics_correction": False}
        config["output_csv"] = str(curve.resolve())
        config["log_file"] = str((phase / "run.log").resolve())
        config["sweep"]["write_state_file"] = str(state.resolve())
        if previous is None:
            config["sweep"].pop("initial_state_file", None)
        else:
            config["sweep"]["initial_state_file"] = str(previous.resolve())
        if phase_name == "20_gate_sweep":
            bias_points = [float(value) for value in config["sweep"]["bias_points"]
                           if float(value) <= gate + 1e-12]
            if not bias_points or not math.isclose(bias_points[-1], gate, abs_tol=1e-12):
                raise RuntimeError(f"M65 diagnostic gate absent from M8 path: {item}")
            config["sweep"]["stop"] = gate
            config["sweep"]["bias_points"] = bias_points
            config["sweep"]["write_vtk"] = True
            config["sweep"]["vtk_prefix"] = str((phase / "vtk/state").resolve())
        config["simplemos_m65"] = {"single_axis": "Si.ni", "matched_ni_cm3":
                                    float(contract["intervention"]["matched_vela_ni_cm3"]),
                                    "production_material_changed": False}
        config_path = phase / "config.json"
        write_json(config_path, config)
        if not complete:
            status = m10.execute_runner(config_path, RUNNER)
            if not status.get("converged"):
                raise RuntimeError(f"M65 Vela failed: {device} {drain} {phase_name}")
        previous, final_curve, final_state = state, curve, state
    assert final_curve is not None and final_state is not None
    vtks = sorted((root / "20_gate_sweep/vtk").glob("state_*.vtk"))
    if not vtks:
        raise RuntimeError(f"M65 final VTK missing: {item}")
    return {**item, "curve": portable(final_curve), "state": portable(final_state),
            "vtk": portable(vtks[-1]), "curve_sha256": sha256(final_curve),
            "state_sha256": sha256(final_state), "vtk_sha256": sha256(vtks[-1])}


def run_vela(contract: dict[str, Any], jobs: int, force: bool) -> dict[str, Any]:
    materials = matched_materials(contract)
    work = cases(contract)
    with ThreadPoolExecutor(max_workers=min(max(1, jobs), 8)) as pool:
        workflows = list(pool.map(lambda item: vela_workflow(item, contract, materials, force), work))
    manifest = {"schema": "vela.simplemos.sdevice.m65_vela_manifest.v1",
                "status": "complete", "contract_sha256": sha256(CONTRACT),
                "materials": portable(materials), "materials_sha256": sha256(materials),
                "workflows": workflows}
    write_json(OUTPUT / "vela_manifest.json", manifest)
    return manifest


def export_sentaurus(contract: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    required = {"ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential",
                "eDensity", "hDensity", "srhRecombination", "BandgapNarrowing",
                "EffectiveIntrinsicDensity", "ConductionBandEnergy", "ValenceBandEnergy"}
    exported = []
    raw = OUTPUT / "sentaurus_raw/sentaurus_bundle"
    for item in manifest["cases"]:
        relative = (REPO / item["deck"]).parent.relative_to(OUTPUT / "sentaurus_bundle")
        root = raw / relative
        name = item["run_case"]
        tdrs = sorted(root.glob(f"{name}_state*.tdr"))
        plots = sorted(root.glob(f"IdVg_{name}*des.plt"))
        logs = sorted(root.glob("*.log"))
        if len(tdrs) != 1 or len(plots) != 1:
            raise RuntimeError(f"M65 expected one TDR/PLT for {name}: {len(tdrs)}/{len(plots)}")
        log_text = "\n".join(path.read_text(errors="replace") for path in logs)
        if "without bandgap narrowing" not in log_text.lower():
            raise RuntimeError(f"M65 no-BGN log qualification failed: {name}")
        state_id = f"{item['device']}_vd_{voltage_tag(float(item['drain_voltage_V']))}"
        export_dir = OUTPUT / "sentaurus_exports" / state_id
        run([str(IMPORTER), "--tdr", str(tdrs[0]), "--export-dir", str(export_dir)], capture=True)
        names = m10.manifest_field_names(export_dir)
        missing = required - names
        if missing:
            raise RuntimeError(f"M65 state missing fields {sorted(missing)}: {state_id}")
        parsed = m60.parse_current(plots[0], [float(item["gate_voltage_V"])],
                                   float(contract["acceptance"]["exact_bias_tolerance_V"]))
        exported.append({**item, "state_id": state_id, "tdr": portable(tdrs[0]),
                         "tdr_sha256": sha256(tdrs[0]), "current": portable(plots[0]),
                         "current_sha256": sha256(plots[0]), "export_dir": portable(export_dir),
                         "field_manifest_sha256": sha256(export_dir / "field_manifest.json"),
                         "drain_current_A_per_um": abs(float(parsed[0]["drain_total_A_per_um"])),
                         "no_bgn_log_qualified": True})
    result = {"schema": "vela.simplemos.sdevice.m65_sentaurus_exports.v1",
              "status": "complete", "states": exported}
    write_json(OUTPUT / "sentaurus_export_manifest.json", result)
    return result


def vela_state(path: Path) -> dict[str, dict[int, float]]:
    rows = read_csv(path)
    names = ("psi", "phin", "phip", "electrons_m3", "holes_m3")
    return {name: {int(row["node_id"]): float(row[name]) for row in rows} for name in names}


def exact_vela_current(path: Path, gate: float, tolerance: float) -> float:
    rows = [row for row in read_csv(path)
            if abs(float(row["bias_V"]) - gate) <= tolerance]
    if len(rows) != 1:
        raise RuntimeError(f"M65 exact Vela gate count {len(rows)} at {gate}: {path}")
    return abs(float(rows[0]["current_total_A_per_um"]))


def baseline_map() -> dict[tuple[str, float], dict[str, float]]:
    result = {}
    for row in read_csv(M64_POINTS):
        result[(row["device"], round(float(row["drain_voltage_V"]), 12),
                round(float(row["gate_voltage_V"]), 12))] = {
            "sentaurus": float(row["sentaurus_off_A_per_um"]),
            "vela": float(row["vela_off_A_per_um"]),
            "error": float(row["error_off_dex"])}
    return result


def analyze_state(item: dict[str, Any], workflow: dict[str, Any],
                  matched_ni: float) -> tuple[dict[str, Any], list[float]]:
    export_dir = REPO / item["export_dir"]
    sent = {name: m47.scalar_field(export_dir, field) for name, field in {
        "psi": "ElectrostaticPotential", "phin": "eQuasiFermiPotential",
        "phip": "hQuasiFermiPotential", "n": "eDensity", "p": "hDensity",
        "srh": "srhRecombination", "ni": "EffectiveIntrinsicDensity",
        "ec": "ConductionBandEnergy", "ev": "ValenceBandEnergy"}.items()}
    vela = vela_state(REPO / workflow["state"])
    common = sorted(set(sent["psi"]) & set(vela["psi"]))
    if common != sorted(sent["psi"]):
        raise RuntimeError(f"M65 Silicon node mismatch: {item['state_id']}")
    coords = m47.coordinates(export_dir)
    mesh = read_json(M8 / "vela" / item["device"] / "mesh.json")
    mesh_coords = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    coordinate_error = max(max(abs(coords[node][axis] - mesh_coords[node][axis])
                               for axis in (0, 1)) for node in common)
    _, source_support = m47.interface_nodes(export_dir)
    source_contact = m47.contact_nodes(export_dir, "source") & set(common)
    sent_ni_values = [sent["ni"][node] for node in common]
    identities, density_errors = [], []
    for node in common:
        n_vela_cm3 = vela["electrons_m3"][node] / 1e6
        n_sent_cm3 = sent["n"][node]
        if n_vela_cm3 <= 0.0 or n_sent_cm3 <= 0.0:
            continue
        observed = math.log10(n_vela_cm3 / n_sent_cm3)
        predicted = (math.log10(matched_ni / sent["ni"][node]) +
                     ((vela["psi"][node] - vela["phin"][node]) -
                      (sent["psi"][node] - sent["phin"][node])) / (VT * math.log(10.0)))
        identities.append(abs(observed - predicted))
        if node in source_support:
            density_errors.append(abs(observed))
    sent_qf_barrier = max(sent["phin"][node] - sent["psi"][node] for node in source_support)
    vela_qf_barrier = max(vela["phin"][node] - vela["psi"][node] for node in source_support)
    sent_source_psi = median([sent["psi"][node] for node in source_contact])
    vela_source_psi = median([vela["psi"][node] for node in source_contact])
    sent_electro_barrier = sent_source_psi - min(sent["psi"][node] for node in source_support)
    vela_electro_barrier = vela_source_psi - min(vela["psi"][node] for node in source_support)
    sent_ec_barrier = max(sent["ec"][node] for node in source_support) - median(
        [sent["ec"][node] for node in source_contact])
    vtks = REPO / workflow["vtk"]
    vela_srh = m47.vtk_scalar(vtks, "SRHRecombinationCm3PerS", len(mesh["nodes"]))
    srh_floor = max(max(abs(sent["srh"][node]) for node in common),
                    max(abs(vela_srh[node]) for node in common)) * 1e-14
    srh_logs = [abs(math.log10(abs(vela_srh[node]) / abs(sent["srh"][node])))
                for node in common if abs(vela_srh[node]) > srh_floor and abs(sent["srh"][node]) > srh_floor]
    row = {"device": item["device"], "drain_voltage_V": item["drain_voltage_V"],
           "gate_voltage_V": item["gate_voltage_V"], "silicon_node_count": len(common),
           "maximum_coordinate_error_um": coordinate_error,
           "sentaurus_ni_min_cm3": min(sent_ni_values), "sentaurus_ni_max_cm3": max(sent_ni_values),
           "sentaurus_ni_relative_spread": (max(sent_ni_values) - min(sent_ni_values)) / median(sent_ni_values),
           "vela_matched_ni_cm3": matched_ni,
           "ni_log10_ratio_dex": math.log10(matched_ni / median(sent_ni_values)),
           "boltzmann_identity_p95_residual_dex": percentile(identities, 0.95),
           "boltzmann_identity_max_residual_dex": max(identities),
           "source_electron_density_log_error_p95_dex": percentile(density_errors, 0.95),
           "sentaurus_source_qf_barrier_V": sent_qf_barrier,
           "vela_source_qf_barrier_V": vela_qf_barrier,
           "qf_barrier_delta_V": vela_qf_barrier - sent_qf_barrier,
           "qf_barrier_current_proxy_dex": -(vela_qf_barrier - sent_qf_barrier) / (VT * math.log(10.0)),
           "sentaurus_source_electrostatic_barrier_V": sent_electro_barrier,
           "vela_source_electrostatic_barrier_V": vela_electro_barrier,
           "electrostatic_barrier_delta_V": vela_electro_barrier - sent_electro_barrier,
           "sentaurus_conduction_band_barrier_eV": sent_ec_barrier,
           "srh_log_error_p95_dex": percentile(srh_logs, 0.95),
           "srh_qualified_node_count": len(srh_logs)}
    return row, identities


def analyze(contract: dict[str, Any], sentaurus: dict[str, Any], vela: dict[str, Any],
            banner: str) -> tuple[dict[str, Any], dict[str, list[dict[str, Any]]]]:
    tolerance = float(contract["acceptance"]["exact_bias_tolerance_V"])
    baseline = baseline_map()
    workflows = {(row["device"], float(row["drain_voltage_V"])): row for row in vela["workflows"]}
    state_rows, current_rows, case_rows, identities = [], [], [], []
    matched_ni = float(contract["intervention"]["matched_vela_ni_cm3"])
    anchor_errors = []
    for item in sentaurus["states"]:
        device, drain, gate = item["device"], float(item["drain_voltage_V"]), float(item["gate_voltage_V"])
        workflow = workflows[(device, drain)]
        vela_current = exact_vela_current(REPO / workflow["curve"], gate, tolerance)
        sent_current = float(item["drain_current_A_per_um"])
        anchor = baseline[(device, round(drain, 12), round(gate, 12))]
        anchor_error = math.log10(sent_current / anchor["sentaurus"])
        anchor_errors.append(abs(anchor_error))
        matched_error = math.log10(vela_current / sent_current)
        current_rows.append({"device": device, "drain_voltage_V": drain, "gate_voltage_V": gate,
                             "m64_sentaurus_no_bgn_A_per_um": anchor["sentaurus"],
                             "replayed_sentaurus_no_bgn_A_per_um": sent_current,
                             "sentaurus_anchor_replay_error_dex": anchor_error,
                             "m64_vela_baseline_no_bgn_A_per_um": anchor["vela"],
                             "matched_ni_vela_A_per_um": vela_current,
                             "baseline_error_dex": anchor["error"],
                             "matched_ni_error_dex": matched_error,
                             "error_change_dex": matched_error - anchor["error"]})
        state_row, residuals = analyze_state(item, workflow, matched_ni)
        state_rows.append(state_row)
        identities.extend(residuals)
        case_rows.extend([
            {"solver": "sentaurus", "device": device, "drain_voltage_V": drain,
             "gate_voltage_V": gate, "converged": 1, "state": item["tdr"],
             "current": item["current"], "no_bgn_qualified": int(item["no_bgn_log_qualified"])},
            {"solver": "vela_matched_ni", "device": device, "drain_voltage_V": drain,
             "gate_voltage_V": gate, "converged": 1, "state": workflow["state"],
             "current": workflow["curve"], "no_bgn_qualified": 1}])

    current_map = {(row["device"], float(row["drain_voltage_V"])): row for row in current_rows}
    state_map = {(row["device"], float(row["drain_voltage_V"])): row for row in state_rows}
    pair_rows, closures = [], []
    for pair in contract["matrix"]["pairs"]:
        low, high, drain = pair["low_device"], pair["high_device"], float(pair["drain_voltage_V"])
        low_current, high_current = current_map[(low, drain)], current_map[(high, drain)]
        baseline_growth = float(high_current["baseline_error_dex"]) - float(low_current["baseline_error_dex"])
        matched_growth = float(high_current["matched_ni_error_dex"]) - float(low_current["matched_ni_error_dex"])
        raw_closure = 1.0 - abs(matched_growth) / max(abs(baseline_growth), 1e-300)
        closure = min(1.0, max(0.0, raw_closure))
        closures.append(closure)
        low_state, high_state = state_map[(low, drain)], state_map[(high, drain)]
        qf_proxy_growth = (float(high_state["qf_barrier_current_proxy_dex"]) -
                           float(low_state["qf_barrier_current_proxy_dex"]))
        electro_proxy_growth = (-(float(high_state["electrostatic_barrier_delta_V"]) -
                                  float(low_state["electrostatic_barrier_delta_V"])) /
                                  (VT * math.log(10.0)))
        pair_rows.append({"low_device": low, "high_device": high, "drain_voltage_V": drain,
                          "diagnostic_gate_voltage_V": pair["diagnostic_gate_voltage_V"],
                          "m64_baseline_ni_pair_growth_dex": baseline_growth,
                          "matched_ni_pair_growth_dex": matched_growth,
                          "pair_growth_change_dex": matched_growth - baseline_growth,
                          "raw_pair_closure_fraction": raw_closure,
                          "classification_pair_closure_fraction": closure,
                          "matched_ni_qf_barrier_pair_proxy_dex": qf_proxy_growth,
                          "matched_ni_electrostatic_barrier_pair_proxy_dex": electro_proxy_growth,
                          "qf_proxy_minus_current_growth_dex": qf_proxy_growth - matched_growth})

    median_closure = median(closures)
    thresholds = contract["analysis"]["classification_thresholds"]
    if median_closure >= float(thresholds["base_ni_dominant_minimum_median_closure_fraction"]):
        classification = "base_intrinsic_density_convention_dominant"
    elif median_closure >= float(thresholds["base_ni_material_minimum_median_closure_fraction"]):
        classification = "base_intrinsic_density_convention_material_but_not_dominant"
    else:
        classification = "base_intrinsic_density_convention_not_dominant"
    acceptance = contract["acceptance"]
    ni_anchor = float(contract["intervention"]["matched_vela_ni_cm3"])
    checks = {
        "contract_frozen": read_json(FREEZE)["contract_sha256"] == sha256(CONTRACT),
        "sentaurus_release": acceptance["required_sentaurus_release"] in banner,
        "case_counts": len(current_rows) == int(acceptance["required_case_count_per_solver"])
                       and len(case_rows) == 2 * int(acceptance["required_case_count_per_solver"]),
        "pair_count": len(pair_rows) == int(acceptance["required_pair_count"]),
        "sentaurus_no_bgn_qualified": all(int(row["no_bgn_qualified"]) == 1 for row in case_rows),
        "m64_sentaurus_anchor_replay": max(anchor_errors) <= float(acceptance["maximum_m64_sentaurus_anchor_replay_error_dex"]),
        "sentaurus_ni_uniform": max(float(row["sentaurus_ni_relative_spread"]) for row in state_rows)
                                 <= float(acceptance["maximum_sentaurus_ni_relative_spread"]),
        "sentaurus_ni_anchor": max(abs(float(row["sentaurus_ni_min_cm3"]) / ni_anchor - 1.0)
                                    for row in state_rows)
                                <= float(acceptance["maximum_sentaurus_ni_relative_anchor_error"]),
        "state_coordinates": max(float(row["maximum_coordinate_error_um"]) for row in state_rows)
                             <= float(acceptance["maximum_state_coordinate_error_um"]),
        "boltzmann_identity": percentile(identities, 0.95)
                              <= float(acceptance["maximum_boltzmann_identity_p95_residual_dex"]),
        "classification_declared": classification in contract["analysis"]["classifications"],
        "production_reference_not_replaced": True,
    }
    checks["all_checks_pass"] = all(checks.values())
    if not checks["all_checks_pass"]:
        classification = "execution_or_state_identity_failure"
    report = {
        "schema": "vela.simplemos.sdevice.m65_nobgn_intrinsic_density_attribution_report.v1",
        "status": "accepted" if checks["all_checks_pass"] else "failed",
        "classification": classification,
        "contract": {"path": portable(CONTRACT), "sha256": sha256(CONTRACT)},
        "execution": {"new_sentaurus_execution": True, "new_vela_execution": True,
                      "sentaurus_case_count": len(current_rows), "vela_workflow_count": len(workflows),
                      "single_axis": "copied Vela Si.ni",
                      "sentaurus_normalized_steps_rescaled_to_preserve_m64_physical_path": True,
                      "production_default_changed": False},
        "summary": {
            "sentaurus_no_bgn_ni_cm3": ni_anchor,
            "vela_baseline_no_bgn_ni_cm3": float(contract["intervention"]["baseline_vela_ni_cm3"]),
            "ni_log10_ratio_dex": math.log10(ni_anchor / float(contract["intervention"]["baseline_vela_ni_cm3"])),
            "median_pair_closure_fraction": median_closure,
            "minimum_pair_closure_fraction": min(closures),
            "maximum_pair_closure_fraction": max(closures),
            "pair_closure_at_or_above_0p8_count": sum(value >= 0.8 for value in closures),
            "median_baseline_pair_growth_dex": median([float(row["m64_baseline_ni_pair_growth_dex"]) for row in pair_rows]),
            "median_matched_ni_pair_growth_dex": median([float(row["matched_ni_pair_growth_dex"]) for row in pair_rows]),
            "median_matched_ni_current_error_dex": median([float(row["matched_ni_error_dex"]) for row in current_rows]),
            "median_qf_barrier_pair_proxy_dex": median([float(row["matched_ni_qf_barrier_pair_proxy_dex"]) for row in pair_rows]),
            "maximum_m64_sentaurus_anchor_replay_error_dex": max(anchor_errors),
            "maximum_sentaurus_ni_relative_spread": max(float(row["sentaurus_ni_relative_spread"]) for row in state_rows),
            "boltzmann_identity_p95_residual_dex": percentile(identities, 0.95),
            "maximum_coordinate_error_um": max(float(row["maximum_coordinate_error_um"]) for row in state_rows)},
        "causal_scope": {
            "closed": "The only causal intervention is the copied Vela Silicon base intrinsic density under explicit no-BGN; pair closure measures how much of M64's larger BGN-independent NWell growth this convention explains.",
            "state_identity": "The matched-ni electron-density difference is independently checked against the electrostatic and electron quasi-Fermi potentials through the Boltzmann identity.",
            "not_closed": "Any residual matched-ni pair growth remains a self-consistent electrostatics, statistics, transport, or material-convention difference; M65 does not tune those terms."},
        "acceptance": checks}
    return report, {"current": current_rows, "states": state_rows, "pairs": pair_rows,
                    "cases": case_rows}


def freeze_results(report: dict[str, Any], ledgers: dict[str, list[dict[str, Any]]]) -> None:
    write_csv(CURRENT, ledgers["current"])
    write_csv(STATES, ledgers["states"])
    write_csv(PAIRS, ledgers["pairs"])
    write_csv(CASES, ledgers["cases"])
    write_json(REPORT, report)
    summary = report["summary"]
    DOC.parent.mkdir(parents=True, exist_ok=True)
    DOC.write_text(f"""# SimpleMOS M65 no-BGN 本征浓度归因

## 结论

M65 分类为 `{report['classification']}`。Sentaurus 显式 no-BGN 的硅区 `EffectiveIntrinsicDensity` 为 `{float(summary['sentaurus_no_bgn_ni_cm3']):.16g} cm^-3`，Vela 原基准为 `{float(summary['vela_baseline_no_bgn_ni_cm3']):.16g} cm^-3`；M65 只在复制材料文件中把 Si `ni` 匹配到前者，未修改生产材料、BGN、HFS、SG、接触、网格或收敛参数。

8 个冻结 NWell 配对的误差增幅由中位 `{float(summary['median_baseline_pair_growth_dex']):.6f} dex` 变为 `{float(summary['median_matched_ni_pair_growth_dex']):.6f} dex`，中位闭合率 `{float(summary['median_pair_closure_fraction']):.2%}`，其中 `{int(summary['pair_closure_at_or_above_0p8_count'])}/8` 个配对达到 80% 闭合。

状态账本覆盖 16 个固定诊断点。Sentaurus `ni` 最大相对空间展宽为 `{float(summary['maximum_sentaurus_ni_relative_spread']):.3e}`；匹配后的电子浓度差由 `ni + (psi-phin)` Boltzmann 恒等式闭合，p95 残差 `{float(summary['boltzmann_identity_p95_residual_dex']):.3e} dex`。这使电势/准费米势垒账本可用于解释匹配后剩余项，而不是把其误记为基础 `ni` 差异。

机器报告：`{portable(REPORT)}`。
""", encoding="utf-8", newline="\n")
    write_json(ARTIFACT, {"schema": "vela.validation.artifact.v1",
                          "title": "SimpleMOS M65 no-BGN intrinsic-density attribution",
                          "status": report["status"], "classification": report["classification"],
                          "summary": "Matched-Si-ni causal intervention at the eight frozen M63 pair gates",
                          "report": portable(REPORT),
                          "ledgers": [portable(CURRENT), portable(STATES), portable(PAIRS), portable(CASES)]})
    artifacts = [REPORT, CURRENT, STATES, PAIRS, CASES, DOC, ARTIFACT]
    write_json(EVIDENCE, {"schema": "vela.simplemos.sdevice.m65_nobgn_intrinsic_density_attribution_evidence.v1",
                          "status": "frozen" if report["acceptance"]["all_checks_pass"] else "failed",
                          "classification": report["classification"],
                          "contract_sha256": sha256(CONTRACT),
                          "new_sentaurus_execution": True, "new_vela_execution": True,
                          "production_reference_replaced": False,
                          "source_hashes": read_json(FREEZE)["upstream_hashes"],
                          "implementation_hashes": {portable(SCRIPT): sha256(SCRIPT)},
                          "artifacts": {portable(path): sha256(path) for path in artifacts},
                          "acceptance": report["acceptance"]})


def verify() -> dict[str, Any]:
    validate_contract()
    evidence, report = read_json(EVIDENCE), read_json(REPORT)
    if evidence.get("status") != "frozen" or report.get("status") != "accepted":
        raise ValueError("M65 evidence is not accepted and frozen")
    if evidence["implementation_hashes"].get(portable(SCRIPT)) != sha256(SCRIPT):
        raise ValueError("M65 implementation hash changed")
    for relative, expected in evidence["artifacts"].items():
        if sha256(REPO / relative) != expected:
            raise ValueError(f"M65 artifact hash changed: {relative}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-contract", action="store_true")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run-sentaurus", action="store_true")
    parser.add_argument("--run-vela", action="store_true")
    parser.add_argument("--export", action="store_true")
    parser.add_argument("--analyze", action="store_true")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--all", action="store_true")
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
        print(json.dumps({"status": "frozen_before_execution", "contract_sha256": sha256(CONTRACT)}))
        return
    contract = validate_contract()
    if args.verify:
        print(json.dumps(verify()["summary"], indent=2))
        return
    if args.prepare or args.all:
        smanifest = prepare_sentaurus(contract, args.force)
    else:
        smanifest = read_json(OUTPUT / "sentaurus_manifest.json")
    banner = ""
    if args.run_sentaurus or args.all:
        banner = run_sentaurus(smanifest, args.ssh_target, args.ssh_bin, args.scp_bin,
                               args.ssh_config, args.remote_root, args.jobs)
    elif (OUTPUT / "sentaurus_banner.txt").is_file():
        banner = (OUTPUT / "sentaurus_banner.txt").read_text(encoding="utf-8").strip()
    if args.run_vela or args.all:
        vmanifest = run_vela(contract, args.jobs, args.force)
    else:
        vmanifest = read_json(OUTPUT / "vela_manifest.json")
    if args.export or args.all:
        exports = export_sentaurus(contract, smanifest)
    else:
        exports = read_json(OUTPUT / "sentaurus_export_manifest.json")
    if args.analyze or args.all:
        report, ledgers = analyze(contract, exports, vmanifest, banner)
        freeze_results(report, ledgers)
        print(json.dumps({"status": report["status"], "classification": report["classification"],
                          "summary": report["summary"], "acceptance": report["acceptance"]}, indent=2))


if __name__ == "__main__":
    main()
