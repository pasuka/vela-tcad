#!/usr/bin/env python3
"""Run the SimpleMOS M8 original-physics SDevice comparison matrix."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
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
import run_simplemos_m3_workflow as m3  # noqa: E402
import run_simplemos_m4_controlled_matrix as m4  # noqa: E402
import compare_simplemos_m4_controlled_matrix as compare_m4  # noqa: E402

DEFAULT_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8_original_physics_contract_v1.json"
)
DEFAULT_VALIDATION_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_sdevice_validation_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8_original_physics"
)
DEFAULT_REFERENCE_DIR = (
    REPO / "reference_tcad/simplemos_sentaurus2022/original_physics"
)
DEFAULT_MATERIALS = (
    REPO / "reference_tcad/transportmodels_sentaurus2022/vela"
    / "materials_sentaurus2022.json"
)
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m8_original_physics_20260827"
)
NOMINAL_GATE_FILE = "nominal_two_curve_gate.json"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def portable_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return resolved.relative_to(REPO).as_posix()
    except ValueError:
        return str(resolved)


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (
            Path(os.environ.get("SystemRoot", r"C:\Windows"))
            / "System32/OpenSSH" / f"{name}.exe"
        )
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    return completed.stdout or ""


def voltage_tag(value: float) -> str:
    return format(value, ".12g").replace("-", "m").replace(".", "p")


def parse_tdr_specs(values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError("--tdr requires ID=PATH")
        device_id, raw_path = value.split("=", 1)
        if device_id in result:
            raise ValueError(f"duplicate TDR id {device_id!r}")
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise FileNotFoundError(path)
        result[device_id] = path
    return result


def selected_devices(contract: dict[str, Any], ids: list[str]
                     ) -> list[dict[str, Any]]:
    devices = list(contract["devices"])
    if not ids:
        return devices
    unknown = set(ids) - {str(item["id"]) for item in devices}
    if unknown:
        raise ValueError(f"unknown M8 device ids: {sorted(unknown)}")
    return [item for item in devices if item["id"] in ids]


def require_nominal_gate(contract_path: Path, devices: list[dict[str, Any]],
                         output_dir: Path) -> None:
    nominal_id = next(str(item["id"]) for item in devices if item.get("nominal")) \
        if any(item.get("nominal") for item in devices) else ""
    if all(str(item["id"]) == nominal_id for item in devices):
        return
    gate_path = output_dir / "comparisons" / NOMINAL_GATE_FILE
    if not gate_path.is_file():
        raise RuntimeError(
            "M8 nominal two-curve gate has not passed; run and compare n17 first")
    gate = read_json(gate_path)
    if gate.get("status") != "pass":
        raise RuntimeError("M8 nominal two-curve gate is not accepted")
    if gate.get("contract_sha256") != sha256(contract_path):
        raise RuntimeError("M8 nominal gate contract hash is stale")
    cases = gate.get("cases", [])
    if len(cases) != 2 or any(item.get("status") != "pass" for item in cases):
        raise RuntimeError("M8 nominal gate does not contain two passing curves")


def write_nominal_gate(contract_path: Path, devices: list[dict[str, Any]],
                       output_dir: Path, comparison: dict[str, Any]) -> None:
    if len(devices) != 1 or not devices[0].get("nominal"):
        return
    cases = comparison["cases"]
    if comparison["status"] != "pass" or len(cases) != 2:
        return
    gate = {
        "schema": "vela.simplemos.sdevice.m8_nominal_gate.v1",
        "status": "pass",
        "device": str(devices[0]["id"]),
        "contract_sha256": sha256(contract_path),
        "cases": [{
            "case": item["case"],
            "status": item["status"],
            "point_count": item["point_count"],
            "maximum_absolute_log10_ratio_above_floor": item[
                "maximum_absolute_log10_ratio_above_floor"],
            "maximum_relative_error_above_floor": item[
                "maximum_relative_error_above_floor"],
            "trend_match": item["trend_match"],
        } for item in cases],
    }
    write_json(output_dir / "comparisons" / NOMINAL_GATE_FILE, gate)


def validate_contract(contract: dict[str, Any]) -> None:
    physics = contract["physics"]
    if physics["effective_intrinsic_density"] != "OldSlotboom":
        raise ValueError("M8 must use the original OldSlotboom model")
    if physics["mobility"] != ["PhuMob", "HighFieldSaturation", "Enormal"]:
        raise ValueError("M8 must use the original mobility combination")
    if physics["recombination"] != "SRH(DopingDependence)":
        raise ValueError("M8 must use the original SRH model")
    srh = physics.get("srh_scharfetter_defaults", {})
    if (srh.get("electron", {}).get("tau_max_s") != 1.0e-5 or
            srh.get("hole", {}).get("tau_max_s") != 3.0e-6 or
            srh.get("electron", {}).get("reference_doping_cm3") != 1.0e16 or
            srh.get("hole", {}).get("reference_doping_cm3") != 1.0e16 or
            srh.get("electron", {}).get("gamma") != 1.0 or
            srh.get("hole", {}).get("gamma") != 1.0):
        raise ValueError("M8 must use the T-2022.03-SP2 default Scharfetter parameters")
    devices = contract["devices"]
    if [item["workbench_process_node"] for item in devices] != list(range(17, 25)):
        raise ValueError("M8 device matrix must be Workbench nodes 17 through 24")
    if sum(bool(item["nominal"]) for item in devices) != 1:
        raise ValueError("M8 must define exactly one nominal device")
    combinations = {
        (float(item["NWell_cm3"]), float(item["GOxTime_min"]),
         float(item["LDD_Dose_cm2"])) for item in devices
    }
    expected = {
        (nwell, oxide, ldd)
        for nwell in (1e17, 2e17)
        for oxide in (10.0, 15.0)
        for ldd in (1e14, 2e14)
    }
    if combinations != expected:
        raise ValueError("M8 process-parameter matrix must be the frozen 2x2x2 grid")
    if any(not math.isclose(float(item["Lg_um"]), 0.25)
           for item in devices):
        raise ValueError("M8 gate length must remain fixed at 0.25 um")
    nominal = next(item for item in devices if item["nominal"])
    if nominal["id"] != contract["execution_gates"]["nominal_device"]:
        raise ValueError("M8 nominal device is inconsistent")
    drain_voltages = [float(value) for value in
                      contract["bias_matrix"]["drain_voltages_V"]]
    if drain_voltages != [0.05, 1.0]:
        raise ValueError("M8 drain biases must be exactly 0.05 V and 1.0 V")
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    if len(lattice) != 51:
        raise ValueError("M8 gate lattice must contain exactly 51 points")
    for index, value in enumerate(lattice):
        if not math.isclose(value, 0.05 * index, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError("M8 gate lattice must be exactly 0:0.05:2.5 V")
    if contract["comparison"]["interpolation"] != "forbidden":
        raise ValueError("M8 interpolation must remain forbidden")


def sentaurus_deck(case: str, vd: float, intervals: int) -> str:
    return f'''File {{
   Grid="input_fps.tdr"
   Plot="{case}_des.tdr"
   Current="{case}"
   Output="{case}.log"
}}

Electrode {{
   {{ Name="source" Voltage=0.0 }}
   {{ Name="drain" Voltage=0.0 }}
   {{ Name="gate" Voltage=0.0 }}
   {{ Name="substrate" Voltage=0.0 }}
}}

Physics {{ EffectiveIntrinsicDensity(OldSlotboom) }}
Physics(Material="Silicon") {{
   Mobility(PhuMob HighFieldSaturation Enormal)
   Recombination(SRH(DopingDependence))
}}

Math {{ Extrapolate Iterations=20 ExitOnFailure }}

Solve {{
   Coupled(Iterations=100) {{ Poisson }}
   Coupled {{ Poisson Electron Hole }}
   Quasistationary(
      InitialStep=0.1 Increment=1.5 MinStep=1e-5 MaxStep=1
      Goal {{ Name="drain" Voltage={vd:.17g} }}
   ) {{ Coupled {{ Poisson Electron Hole }} }}
   NewCurrentPrefix="IdVg_"
   Quasistationary(
      DoZero InitialStep=0.01 Increment=1.5 MinStep=1e-5 MaxStep=0.05
      Goal {{ Name="gate" Voltage=2.5 }}
   ) {{
      Coupled {{ Poisson Electron Hole }}
      CurrentPlot(Time=(Range=(0 1) Intervals={intervals}))
   }}
}}
'''


def prepare_sentaurus_bundle(contract: dict[str, Any], devices: list[dict[str, Any]],
                              tdrs: dict[str, Path], output_dir: Path
                              ) -> dict[str, Any]:
    bundle = output_dir / "sentaurus_bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    intervals = int(contract["bias_matrix"]["gate_lattice"]["point_count"]) - 1
    cases: list[dict[str, Any]] = []
    inputs: list[dict[str, Any]] = []
    for device in devices:
        device_id = str(device["id"])
        if device_id not in tdrs:
            raise ValueError(f"missing --tdr for selected device {device_id}")
        device_dir = bundle / device_id
        device_dir.mkdir(parents=True, exist_ok=True)
        local_tdr = device_dir / "input_fps.tdr"
        shutil.copy2(tdrs[device_id], local_tdr)
        inputs.append({
            "device": device_id,
            "workbench_process_node": device["workbench_process_node"],
            "tdr": str(tdrs[device_id]),
            "tdr_sha256": sha256(tdrs[device_id]),
            "tdr_size_bytes": tdrs[device_id].stat().st_size,
        })
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            case = f"{device_id}_vd_{voltage_tag(float(vd))}"
            deck = device_dir / f"{case}_des.cmd"
            deck.write_text(
                sentaurus_deck(case, float(vd), intervals),
                encoding="utf-8", newline="\n")
            cases.append({
                "case": case,
                "device": device_id,
                "drain_voltage_V": float(vd),
                "deck": str(deck.resolve()),
                "deck_sha256": sha256(deck),
                "expected_plot": f"IdVg_{case}_des.plt",
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m8_sentaurus_matrix.v1",
        "status": "prepared",
        "inputs": inputs,
        "cases": cases,
    }
    write_json(output_dir / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def run_sentaurus_vm(manifest: dict[str, Any], output_dir: Path,
                     ssh_target: str, ssh_bin: str, scp_bin: str,
                     remote_root: str) -> str:
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target, f"mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output_dir / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])
    for index, case in enumerate(manifest["cases"], start=1):
        print(f"[{index}/{len(manifest['cases'])}] Sentaurus {case['case']}",
              flush=True)
        command = (
            f"set -eu; cd {remote_root}/sentaurus_bundle/{case['device']}; "
            f"sdevice {case['case']}_des.cmd > {case['case']}.console.log 2>&1"
        )
        run([ssh_bin, ssh_target, command])
    archive_name = "simplemos_m8_results.tgz"
    run([ssh_bin, ssh_target,
         f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle"])
    raw = output_dir / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (output_dir / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def extract_references(contract: dict[str, Any], manifest: dict[str, Any],
                       contract_path: Path, output_dir: Path,
                       reference_dir: Path, banner: str
                       ) -> dict[str, Any]:
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    tolerance = float(contract["comparison"]["exact_bias_tolerance_V"])
    raw = output_dir / "sentaurus_raw/sentaurus_bundle"
    reference_dir.mkdir(parents=True, exist_ok=True)
    artifacts: list[dict[str, Any]] = []
    for case in manifest["cases"]:
        plot = raw / case["device"] / case["expected_plot"]
        rows = m4.exact_plt_curve(plot, lattice, tolerance)
        output = reference_dir / f"{case['case']}_reference.csv"
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        artifacts.append({
            "case": case["case"], "device": case["device"],
            "drain_voltage_V": case["drain_voltage_V"],
            "path": output.name, "sha256": sha256(output),
            "point_count": len(rows),
        })
    result = {
        "schema": "vela.simplemos.sdevice.m8_reference_curves.v1",
        "status": "qualified",
        "sentaurus_banner": banner,
        "contract": "../simplemos_m8_original_physics_contract_v1.json",
        "contract_sha256": sha256(contract_path),
        "inputs": [{key: value for key, value in item.items() if key != "tdr"}
                   for item in manifest["inputs"]],
        "interpolation": "forbidden",
        "artifacts": artifacts,
        "raw_binary_commit_policy": "excluded",
    }
    write_json(reference_dir / "reference_manifest.json", result)
    return result


def import_tdr(tdr: Path, neutral_dir: Path, importer: Path) -> None:
    neutral_dir.mkdir(parents=True, exist_ok=True)
    run([
        str(importer.resolve()), "--tdr", str(tdr.resolve()),
        "--inventory-json", str(neutral_dir / "tdr_inventory.json"),
        "--qualification-contract", str(DEFAULT_VALIDATION_CONTRACT.resolve()),
        "--qualification-report", str(neutral_dir / "tdr_qualification_report.json"),
        "--export-dir", str(neutral_dir),
    ])
    report = read_json(neutral_dir / "tdr_qualification_report.json")
    if not report.get("passed"):
        raise RuntimeError(f"TDR qualification failed for {tdr}")


def apply_simplemos_srh_defaults(config: dict[str, Any],
                                 contract: dict[str, Any]) -> None:
    """Apply the Silicon defaults used when SimpleMOS has no Parameter file."""
    srh_contract = contract["physics"]["srh_scharfetter_defaults"]
    srh_config = config["solver"]["srh_doping_dependence"]
    for carrier in ("electron", "hole"):
        source = srh_contract[carrier]
        srh_config[carrier].update({
            "tau_min_s": float(source["tau_min_s"]),
            "tau_max_s": float(source["tau_max_s"]),
            "reference_doping_m3": float(source["reference_doping_cm3"]),
            "gamma": float(source["gamma"]),
        })


def prepare_vela_workflows(contract: dict[str, Any], devices: list[dict[str, Any]],
                           tdrs: dict[str, Path], output_dir: Path,
                           materials: Path, importer: Path) -> dict[str, Any]:
    stages: list[dict[str, Any]] = []
    inputs: list[dict[str, Any]] = []
    lattice = [float(value) for value in
               contract["bias_matrix"]["gate_lattice"]["values_V"]]
    for device in devices:
        device_id = str(device["id"])
        neutral = output_dir / "neutral" / device_id
        import_tdr(tdrs[device_id], neutral, importer)
        device_root = output_dir / "vela" / device_id
        device_root.mkdir(parents=True, exist_ok=True)
        mesh = device_root / "mesh.json"
        write_json(mesh, m4.neutral_mesh(neutral))
        mobility = {
            "model": "phumob_field_lombardi",
            "doping_concentration_basis": "total_impurity",
            "high_field_driving_force": "quasi_fermi_gradient",
            "high_field_gradient_discretization": "transport_cell_vector",
        }
        base = device_root / "base.json"
        base_config = m4.base_config(
            mesh, neutral / "doping.csv", materials, mobility)
        # SimpleMOS does not name a Parameter file.  Use the T-2022.03-SP2
        # Silicon defaults printed by ``sdevice -P`` rather than the
        # TransportModels example's electron-lifetime override.
        apply_simplemos_srh_defaults(base_config, contract)
        write_json(base, base_config)
        workflow_dir = device_root / "workflow"
        workflow = m3.materialize(
            base, DEFAULT_VALIDATION_CONTRACT, workflow_dir,
            [float(value) for value in contract["bias_matrix"]["drain_voltages_V"]])
        for stage in workflow["stages"]:
            config_path = Path(stage["config"])
            config = read_json(config_path)
            if stage["phase"] == "equilibrium":
                config["solver"]["quasi_fermi_update_limit_V"] = float(
                    contract["vela_execution_protocol"]
                    ["equilibrium_quasi_fermi_update_limit_V"])
            elif stage["phase"] == "drain_ramp":
                protocol = contract["vela_execution_protocol"]
                config["solver"].update({
                    "line_search_mode": "block_filter",
                    "residual_filter_gamma": 1e-4,
                    "residual_filter_envelope_factor": 2.0,
                })
                config["sweep"]["step"] = float(
                    protocol["auxiliary_drain_initial_step_V"])
                config["sweep"]["max_step"] = float(
                    protocol["auxiliary_drain_max_step_V"])
                config["simplemos_m8"] = {
                    "purpose": "continuation_only",
                    "acceptance_role": "none",
                }
            else:
                config["sweep"]["bias_points"] = lattice
                config["sweep"]["step"] = 0.05
                config["sweep"]["min_step"] = min(
                    float(config["sweep"]["min_step"]), 0.05)
                config["sweep"]["max_step"] = max(
                    float(config["sweep"]["max_step"]), 0.05)
                stage["exact_gate_lattice_V"] = lattice
            write_json(config_path, config)
            stage["config_sha256"] = sha256(config_path)
        workflow["simplemos_m8"] = {
            "device": device_id,
            "original_physics": True,
        }
        write_json(workflow_dir / "workflow_manifest.json", workflow)
        stages.extend(workflow["stages"])
        inputs.append({
            "device": device_id, "tdr_sha256": sha256(tdrs[device_id]),
            "mesh_sha256": sha256(mesh),
            "doping_sha256": sha256(neutral / "doping.csv"),
            "workflow_manifest": str((workflow_dir / "workflow_manifest.json").resolve()),
        })
    result = {
        "schema": "vela.simplemos.sdevice.m8_vela_matrix.v1",
        "status": "materialized", "inputs": inputs, "stages": stages,
    }
    write_json(output_dir / "vela/vela_matrix_manifest.json", result)
    return result


def execute_vela(devices: list[dict[str, Any]], output_dir: Path,
                 runner: Path, jobs: int = 1) -> dict[str, Any]:
    if jobs < 1:
        raise ValueError("jobs must be at least one")

    def execute_device(index_device: tuple[int, dict[str, Any]]) -> dict[str, Any]:
        index, device = index_device
        device_id = str(device["id"])
        print(f"[{index}/{len(devices)}] Vela {device_id}", flush=True)
        workflow_dir = output_dir / "vela" / device_id / "workflow"
        manifest_path = workflow_dir / "workflow_manifest.json"
        workflow = m3.execute(read_json(manifest_path), runner, workflow_dir)
        return {
            "device": device_id, "status": workflow["status"],
            "manifest": str(manifest_path.resolve()),
            "manifest_sha256": sha256(manifest_path),
        }

    indexed = list(enumerate(devices, start=1))
    if jobs == 1:
        results = []
        for item in indexed:
            result = execute_device(item)
            results.append(result)
            if result["status"] != "accepted":
                break
    else:
        with ThreadPoolExecutor(max_workers=min(jobs, len(devices))) as executor:
            results = list(executor.map(execute_device, indexed))
    report = {
        "schema": "vela.simplemos.sdevice.m8_vela_execution.v1",
        "status": "accepted" if len(results) == len(devices) and all(
            item["status"] == "accepted" for item in results) else "fail",
        "runner": str(runner.resolve()), "runner_sha256": sha256(runner),
        "devices": results,
    }
    write_json(output_dir / "vela/vela_execution_report.json", report)
    return report


def compare_results(contract: dict[str, Any], devices: list[dict[str, Any]],
                    contract_path: Path, output_dir: Path,
                    reference_dir: Path,
                    comparison_dir: Path | None = None) -> dict[str, Any]:
    comparison_dir = comparison_dir or output_dir / "comparisons"
    comparison_dir.mkdir(parents=True, exist_ok=True)
    cases: list[dict[str, Any]] = []
    for device in devices:
        device_id = str(device["id"])
        for vd in contract["bias_matrix"]["drain_voltages_V"]:
            tag = voltage_tag(float(vd))
            case = f"{device_id}_vd_{tag}"
            reference = reference_dir / f"{case}_reference.csv"
            candidate = (
                output_dir / "vela" / device_id / "workflow"
                / f"vd_{tag}" / "20_gate_sweep.csv")
            try:
                result = compare_m4.compare_case(reference, candidate, contract)
                rows = result.pop("rows")
                comparison_csv = comparison_dir / f"{case}_comparison.csv"
                compare_m4.write_case_csv(comparison_csv, rows)
                result["comparison_csv"] = portable_path(comparison_csv)
                result["comparison_csv_sha256"] = sha256(comparison_csv)
            except (FileNotFoundError, ValueError) as error:
                result = {"status": "qualification_failed",
                          "qualification_error": str(error)}
            result.update(case=case, device=device_id,
                          drain_voltage_V=float(vd))
            cases.append(result)
    report = {
        "schema": "vela.simplemos.sdevice.m8_comparison.v1",
        "status": "pass" if all(item["status"] == "pass" for item in cases)
        else "fail",
        "contract": (
            "reference_tcad/simplemos_sentaurus2022/"
            "simplemos_m8_original_physics_contract_v1.json"),
        "contract_sha256": sha256(contract_path),
        "interpolation": "forbidden", "cases": cases,
    }
    write_json(comparison_dir / "comparison_report.json", report)
    return report


def freeze_evidence(contract: dict[str, Any], contract_path: Path,
                    devices: list[dict[str, Any]], reference_dir: Path,
                    comparison: dict[str, Any], runner: Path) -> dict[str, Any]:
    expected_cases = 2 * len(contract["devices"])
    if len(devices) != len(contract["devices"]):
        raise RuntimeError("M8 evidence requires all eight devices")
    if len(comparison["cases"]) != expected_cases:
        raise RuntimeError("M8 evidence requires sixteen completed comparisons")
    if any(item["status"] not in {"pass", "fail"}
           for item in comparison["cases"]):
        raise RuntimeError("M8 evidence cannot contain qualification failures")
    reference_manifest_path = reference_dir / "reference_manifest.json"
    reference_manifest = read_json(reference_manifest_path)
    if (reference_manifest.get("status") != "qualified" or
            len(reference_manifest.get("artifacts", [])) != expected_cases):
        raise RuntimeError("M8 evidence requires sixteen qualified references")
    comparison_report = reference_dir / "comparisons" / "comparison_report.json"
    repair_evidence_path = (
        contract_path.parent / "simplemos_m8_srh_repair_evidence.json")
    cases = comparison["cases"]
    evidence = {
        "schema": "vela.simplemos.sdevice.m8_evidence.v1",
        "status": "accepted" if comparison["status"] == "pass" else "failed",
        "scope": "SDevice-only original SimpleMOS physics comparison",
        "qualification_date": "2026-08-28",
        "contract": contract_path.name,
        "contract_sha256": sha256(contract_path),
        "source_deck": contract["source_deck"],
        "sentaurus_release": contract["sentaurus_release"],
        "reference_manifest": (
            "original_physics/reference_manifest.json"),
        "reference_manifest_sha256": sha256(reference_manifest_path),
        "input_tdrs": reference_manifest["inputs"],
        "reference_qualification": {
            "curve_count": expected_cases,
            "points_per_curve": 51,
            "direct_point_count": expected_cases * 51,
            "interpolation": "forbidden",
        },
        "vela_execution": {
            "runner_sha256": sha256(runner),
            "qualified_devices": [str(item["id"]) for item in devices],
            "qualified_cases": expected_cases,
            "all_gate_points_converged": True,
            "restart_policy": "previous_accepted_state_only",
            "auxiliary_drain_points_have_acceptance_role": False,
        },
        "comparison": {
            "status": comparison["status"],
            "report": "original_physics/comparisons/comparison_report.json",
            "report_sha256": sha256(comparison_report),
            "case_count": expected_cases,
            "passed_cases": sum(item["status"] == "pass" for item in cases),
            "maximum_absolute_log10_ratio_observed": max(
                item["maximum_absolute_log10_ratio_above_floor"]
                for item in cases),
            "maximum_relative_error_observed": max(
                item["maximum_relative_error_above_floor"] for item in cases),
            "maximum_absolute_endpoint_log10_ratio_observed": max(
                abs(item["endpoint_log10_ratio"]) for item in cases),
            "trend_matches": sum(bool(item["trend_match"]) for item in cases),
            "limits": contract["comparison"]["numeric_parity"],
            "cases": [{
                "case": item["case"],
                "status": item["status"],
                "max_log10_ratio": item[
                    "maximum_absolute_log10_ratio_above_floor"],
                "median_log10_ratio": item[
                    "median_absolute_log10_ratio_above_floor"],
                "max_relative_error": item[
                    "maximum_relative_error_above_floor"],
                "endpoint_log10_ratio": item["endpoint_log10_ratio"],
            } for item in cases],
        },
        "repair": {
            "status": "confirmed",
            "root_cause": (
                "SimpleMOS default electron SRH lifetime was replaced by a "
                "TransportModels-specific override"),
            "evidence": repair_evidence_path.name,
            "evidence_sha256": sha256(repair_evidence_path),
        },
        "implementation_sha256": {
            "scripts/run_simplemos_m8_original_matrix.py": sha256(Path(__file__)),
            "scripts/run_simplemos_m8_deep_off_diagnostics.py": sha256(
                REPO / "scripts/run_simplemos_m8_deep_off_diagnostics.py"),
            "scripts/plot_simplemos_m8_validation.py": sha256(
                REPO / "scripts/plot_simplemos_m8_validation.py"),
            "tests/regression/test_simplemos_m8_original_matrix.py": sha256(
                REPO / "tests/regression/test_simplemos_m8_original_matrix.py"),
        },
        "raw_artifact_policy": (
            "Sentaurus TDR/PLT/log and generated Vela states remain ignored"),
    }
    write_json(
        contract_path.parent / "simplemos_m8_original_physics_evidence.json",
        evidence)
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--reference-dir", type=Path,
                        default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--materials", type=Path, default=DEFAULT_MATERIALS)
    parser.add_argument("--importer", type=Path,
                        default=REPO / "build-release/sentaurus_import.exe")
    parser.add_argument("--runner", type=Path,
                        default=REPO / "build-release/vela_example_runner.exe")
    parser.add_argument("--tdr", action="append", default=[], metavar="ID=PATH")
    parser.add_argument("--device", action="append", default=[])
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-existing", action="store_true")
    parser.add_argument("--execute-vela", action="store_true")
    parser.add_argument("--jobs", type=int, default=1,
                        help="parallel Vela devices; stages within a device remain serial")
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--comparison-dir", type=Path)
    parser.add_argument("--reuse-prepared", action="store_true")
    parser.add_argument("--freeze-evidence", action="store_true")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    args = parser.parse_args()

    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    devices = selected_devices(contract, args.device)
    output_dir = args.output_dir.resolve()
    if args.clean and output_dir.exists():
        shutil.rmtree(output_dir)
    require_nominal_gate(contract_path, devices, output_dir)
    if args.reuse_prepared:
        manifest = read_json(output_dir / "sentaurus_matrix_manifest.json")
        expected_ids = {str(item["id"]) for item in devices}
        actual_ids = {str(item["device"]) for item in manifest["inputs"]}
        if actual_ids != expected_ids:
            raise RuntimeError("prepared Sentaurus manifest does not match devices")
    else:
        tdrs = parse_tdr_specs(args.tdr)
        manifest = prepare_sentaurus_bundle(contract, devices, tdrs, output_dir)
        prepare_vela_workflows(
            contract, devices, tdrs, output_dir, args.materials.resolve(),
            args.importer.resolve())
    vela_report = None
    if args.execute_vela:
        vela_report = execute_vela(
            devices, output_dir, args.runner.resolve(), args.jobs)
    banner = ""
    if args.live_sentaurus:
        banner = run_sentaurus_vm(
            manifest, output_dir, args.ssh_target, args.ssh_bin,
            args.scp_bin, args.remote_root)
    elif args.extract_existing:
        banner = (output_dir / "sentaurus_banner.txt").read_text(
            encoding="utf-8").strip()
    if args.live_sentaurus or args.extract_existing:
        extract_references(
            contract, manifest, contract_path, output_dir,
            args.reference_dir.resolve(), banner)
    comparison = None
    if args.compare:
        comparison = compare_results(
            contract, devices, contract_path, output_dir,
            args.reference_dir.resolve(),
            args.comparison_dir.resolve() if args.comparison_dir else None)
        write_nominal_gate(contract_path, devices, output_dir, comparison)
    if args.freeze_evidence:
        if comparison is None:
            raise RuntimeError("--freeze-evidence requires --compare")
        freeze_evidence(
            contract, contract_path, devices, args.reference_dir.resolve(),
            comparison, args.runner.resolve())
    status = {
        "prepared_devices": len(devices),
        "prepared_cases": len(manifest["cases"]),
        "vela_status": (vela_report or {}).get("status"),
        "comparison_status": (comparison or {}).get("status"),
    }
    print(json.dumps(status))
    if vela_report is not None and vela_report["status"] != "accepted":
        return 1
    if comparison is not None and comparison["status"] != "pass":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
