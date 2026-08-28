#!/usr/bin/env python3
"""Run the SimpleMOS M8-A HFS confirmation matrix."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
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
import compare_simplemos_m4_controlled_matrix as compare_m4  # noqa: E402
import run_simplemos_m3_workflow as m3  # noqa: E402
import run_simplemos_m4_controlled_matrix as m4  # noqa: E402
import run_simplemos_m8_original_matrix as m8  # noqa: E402
import run_simplemos_m8a_model_ablation as m8a  # noqa: E402


DEFAULT_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8a_confirmation_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8a_confirmation"
)
DEFAULT_TDRS = {
    "n17": REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "source/n17_fps.tdr",
    "n21": REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8_upstream_tdrs/n21_fps.tdr",
}
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m8a_confirmation_20260828"
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8",
                    newline="\n")


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None)
    return completed.stdout or ""


def validate_contract(contract: dict[str, Any]) -> None:
    if [item["id"] for item in contract["devices"]] != ["n17", "n21"]:
        raise ValueError("confirmation devices must be n17 and n21")
    if [item["id"] for item in contract["variants"]] != [
            "full", "no_hfs", "phumob_only"]:
        raise ValueError("confirmation variants must be full/no_hfs/phumob_only")
    if [float(value) for value in contract["bias_matrix"][
            "drain_voltages_V"]] != [0.05, 1.0]:
        raise ValueError("confirmation drain biases must be 0.05 and 1 V")
    lattice = [float(value) for value in contract["bias_matrix"][
        "gate_lattice"]["values_V"]]
    if len(lattice) != 51 or any(not math.isclose(
            value, index * 0.05, rel_tol=0.0, abs_tol=1e-12)
            for index, value in enumerate(lattice)):
        raise ValueError("confirmation gate lattice must be 0:0.05:2.5 V")
    if contract["comparison"]["interpolation"] != "forbidden":
        raise ValueError("confirmation forbids interpolation")


def validate_tdrs(contract: dict[str, Any], tdrs: dict[str, Path]) -> None:
    expected = {item["id"]: item["tdr_sha256"]
                for item in contract["devices"]}
    if set(tdrs) != set(expected):
        raise ValueError("TDR set must match confirmation devices")
    for device, path in tdrs.items():
        if m8a.sha256(path) != expected[device]:
            raise ValueError(f"{device} TDR hash does not match the contract")


def parse_tdrs(values: list[str]) -> dict[str, Path]:
    if not values:
        return {key: value.resolve() for key, value in DEFAULT_TDRS.items()}
    result = {}
    for value in values:
        device, separator, raw_path = value.partition("=")
        if not separator or device in result:
            raise ValueError("TDR arguments must be unique DEVICE=PATH pairs")
        result[device] = Path(raw_path).resolve()
    return result


def prepare_sentaurus(contract: dict[str, Any], contract_path: Path,
                       tdrs: dict[str, Path], output: Path) -> dict[str, Any]:
    bundle = output / "sentaurus_bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    intervals = len(contract["bias_matrix"]["gate_lattice"]["values_V"]) - 1
    cases = []
    for device in contract["devices"]:
        device_id = device["id"]
        device_dir = bundle / device_id
        device_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(tdrs[device_id], device_dir / "input_fps.tdr")
        for variant in contract["variants"]:
            variant_id = variant["id"]
            case_dir = device_dir / variant_id
            case_dir.mkdir(parents=True, exist_ok=True)
            for vd in contract["bias_matrix"]["drain_voltages_V"]:
                tag = m8a.voltage_tag(float(vd))
                case = f"{device_id}_{variant_id}_vd_{tag}"
                deck = case_dir / f"{case}_des.cmd"
                text = m8a.sentaurus_deck(
                    case, float(vd), variant, intervals).replace(
                        'Grid="input_fps.tdr"', 'Grid="../input_fps.tdr"')
                deck.write_text(text, encoding="utf-8", newline="\n")
                cases.append({
                    "case": case,
                    "device": device_id,
                    "variant": variant_id,
                    "drain_voltage_V": float(vd),
                    "expected_plot": f"IdVg_{case}_des.plt",
                    "deck": m8a.portable_path(deck),
                    "deck_sha256": m8a.sha256(deck),
                })
    manifest = {
        "schema": "vela.simplemos.sdevice.m8a_confirmation_sentaurus.v1",
        "status": "prepared",
        "contract": m8a.portable_path(contract_path),
        "contract_sha256": m8a.sha256(contract_path),
        "tdr_sha256": {key: m8a.sha256(value) for key, value in tdrs.items()},
        "cases": cases,
    }
    write_json(output / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def run_sentaurus(manifest: dict[str, Any], output: Path, ssh_target: str,
                   ssh_bin: str, scp_bin: str, remote_root: str,
                   jobs: int) -> str:
    if jobs < 1:
        raise ValueError("Sentaurus jobs must be at least one")
    banner = run(
        [ssh_bin, ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
        capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, ssh_target,
         f"test ! -e {remote_root} && mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])

    def execute(index_case: tuple[int, dict[str, Any]]) -> None:
        index, case = index_case
        print(f"[{index}/{len(manifest['cases'])}] Sentaurus {case['case']}",
              flush=True)
        remote_dir = (f"{remote_root}/sentaurus_bundle/{case['device']}/"
                      f"{case['variant']}")
        command = (f"set -eu; cd {remote_dir}; "
                   f"sdevice {case['case']}_des.cmd > "
                   f"{case['case']}.console.log 2>&1")
        run([ssh_bin, ssh_target, command])

    indexed = list(enumerate(manifest["cases"], start=1))
    with ThreadPoolExecutor(max_workers=min(jobs, len(indexed))) as executor:
        list(executor.map(execute, indexed))
    archive_name = "simplemos_m8a_confirmation_results.tgz"
    run([ssh_bin, ssh_target,
         f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle"])
    raw = output / "sentaurus_raw"
    raw.mkdir(parents=True, exist_ok=True)
    archive = raw / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(raw, filter="data")
    (output / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def extract_references(contract: dict[str, Any], manifest: dict[str, Any],
                       output: Path, banner: str) -> dict[str, Any]:
    lattice = [float(value) for value in contract["bias_matrix"][
        "gate_lattice"]["values_V"]]
    tolerance = float(contract["comparison"]["exact_bias_tolerance_V"])
    reference_dir = output / "sentaurus_reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    artifacts = []
    for case in manifest["cases"]:
        source = (output / "sentaurus_raw/sentaurus_bundle" / case["device"]
                  / case["variant"] / case["expected_plot"])
        rows = m4.exact_plt_curve(source, lattice, tolerance)
        destination = reference_dir / f"{case['case']}_reference.csv"
        with destination.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]),
                                    lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        artifacts.append({
            "case": case["case"], "device": case["device"],
            "variant": case["variant"], "point_count": len(rows),
            "path": m8a.portable_path(destination),
            "sha256": m8a.sha256(destination),
        })
    report = {
        "schema": "vela.simplemos.sdevice.m8a_confirmation_reference.v1",
        "status": "qualified", "sentaurus_banner": banner,
        "contract_sha256": manifest["contract_sha256"],
        "tdr_sha256": manifest["tdr_sha256"],
        "interpolation": "forbidden", "artifacts": artifacts,
    }
    write_json(reference_dir / "reference_manifest.json", report)
    return report


def patch_workflow(contract: dict[str, Any], variant_id: str,
                   workflow: dict[str, Any]) -> None:
    lattice = [float(value) for value in contract["bias_matrix"][
        "gate_lattice"]["values_V"]]
    for stage in workflow["stages"]:
        config_path = Path(stage["config"])
        config = read_json(config_path)
        if stage["phase"] == "equilibrium":
            config["solver"]["quasi_fermi_update_limit_V"] = 0.1
        elif stage["phase"] == "drain_ramp":
            config["solver"].update({
                "line_search_mode": "block_filter",
                "residual_filter_gamma": 1e-4,
                "residual_filter_envelope_factor": 2.0,
            })
            config["sweep"]["step"] = 0.005
            config["sweep"]["max_step"] = 0.01
        else:
            config["sweep"]["bias_points"] = lattice
            config["sweep"]["step"] = 0.05
            config["sweep"]["min_step"] = min(
                float(config["sweep"]["min_step"]), 0.05)
            config["sweep"]["max_step"] = max(
                float(config["sweep"]["max_step"]), 0.05)
            stage["exact_gate_lattice_V"] = lattice
        config["simplemos_m8a_confirmation"] = {
            "variant": variant_id,
            "paired_sentaurus_configuration": True,
        }
        write_json(config_path, config)
        stage["config_sha256"] = m8a.sha256(config_path)


def prepare_vela(contract: dict[str, Any], tdrs: dict[str, Path], output: Path,
                 materials: Path, importer: Path) -> dict[str, Any]:
    workflows = []
    for device in contract["devices"]:
        device_id = device["id"]
        neutral = output / "neutral" / device_id
        m8.import_tdr(tdrs[device_id], neutral, importer)
        mesh = output / "vela" / device_id / "mesh.json"
        write_json(mesh, m4.neutral_mesh(neutral))
        for variant in contract["variants"]:
            variant_id = variant["id"]
            variant_dir = output / "vela" / device_id / variant_id
            base = variant_dir / "base.json"
            config = m4.base_config(
                mesh, neutral / "doping.csv", materials,
                m8a.mobility_config(variant))
            m8a.apply_variant_physics(config, variant, contract)
            write_json(base, config)
            workflow_dir = variant_dir / "workflow"
            workflow = m3.materialize(
                base, m8a.DEFAULT_VALIDATION_CONTRACT, workflow_dir,
                [float(value) for value in contract["bias_matrix"][
                    "drain_voltages_V"]])
            patch_workflow(contract, variant_id, workflow)
            workflow["simplemos_m8a_confirmation"] = {
                "device": device_id, "variant": variant_id,
            }
            manifest = workflow_dir / "workflow_manifest.json"
            write_json(manifest, workflow)
            workflows.append({
                "device": device_id, "variant": variant_id,
                "manifest": m8a.portable_path(manifest),
            })
    report = {
        "schema": "vela.simplemos.sdevice.m8a_confirmation_vela.v1",
        "status": "materialized", "workflows": workflows,
    }
    write_json(output / "vela/vela_matrix_manifest.json", report)
    return report


def execute_vela(output: Path, runner: Path, jobs: int) -> dict[str, Any]:
    if jobs < 1:
        raise ValueError("Vela jobs must be at least one")
    matrix = read_json(output / "vela/vela_matrix_manifest.json")

    def execute(index_item: tuple[int, dict[str, Any]]) -> dict[str, Any]:
        index, item = index_item
        print(f"[{index}/{len(matrix['workflows'])}] Vela "
              f"{item['device']} {item['variant']}", flush=True)
        manifest = REPO / item["manifest"]
        workflow_dir = manifest.parent
        workflow = m3.execute(read_json(manifest), runner.resolve(), workflow_dir)
        return {**item, "status": workflow["status"],
                "manifest_sha256": m8a.sha256(manifest)}

    indexed = list(enumerate(matrix["workflows"], start=1))
    with ThreadPoolExecutor(max_workers=min(jobs, len(indexed))) as executor:
        results = list(executor.map(execute, indexed))
    report = {
        "schema": "vela.simplemos.sdevice.m8a_confirmation_execution.v1",
        "status": "accepted" if all(item["status"] == "accepted"
                                    for item in results) else "fail",
        "runner": m8a.portable_path(runner),
        "runner_sha256": m8a.sha256(runner), "workflows": results,
    }
    write_json(output / "vela/vela_execution_report.json", report)
    return report


def compare(contract: dict[str, Any], output: Path) -> dict[str, Any]:
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
                result = compare_m4.compare_case(reference, candidate, contract)
                rows = result.pop("rows")
                result["diagnostic_metrics"] = m8a.diagnostic_metrics(
                    rows, contract)
                comparison_csv = destination / f"{case_id}_comparison.csv"
                compare_m4.write_case_csv(comparison_csv, rows)
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
                        and case["variant"] == "full")
        item["effect_vs_full"] = {
            "maximum_log10_ratio_change": (
                item["maximum_absolute_log10_ratio_above_floor"]
                - baseline["maximum_absolute_log10_ratio_above_floor"]),
            "deep_off_vg_0p05_log10_ratio_change": (
                item["diagnostic_metrics"]["deep_off"][1]["absolute_log10_ratio"]
                - baseline["diagnostic_metrics"]["deep_off"][1][
                    "absolute_log10_ratio"]),
            "weak_inversion_max_log10_ratio_change": (
                item["diagnostic_metrics"]["weak_inversion_max"][
                    "absolute_log10_ratio"]
                - baseline["diagnostic_metrics"]["weak_inversion_max"][
                    "absolute_log10_ratio"]),
        }
    report = {
        "schema": "vela.simplemos.sdevice.m8a_confirmation_comparison.v1",
        "status": "complete", "interpolation": "forbidden", "cases": cases,
    }
    write_json(destination / "comparison_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--tdr", action="append", default=[])
    parser.add_argument("--materials", type=Path, default=m8a.DEFAULT_MATERIALS)
    parser.add_argument("--importer", type=Path,
                        default=REPO / "build-release/sentaurus_import.exe")
    parser.add_argument("--runner", type=Path,
                        default=REPO / "build-release/vela_example_runner.exe")
    parser.add_argument("--reuse-prepared", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--extract-existing", action="store_true")
    parser.add_argument("--execute-vela", action="store_true")
    parser.add_argument("--compare", action="store_true")
    parser.add_argument("--vela-jobs", type=int, default=1)
    parser.add_argument("--sentaurus-jobs", type=int, default=2)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=m8a.executable("ssh"))
    parser.add_argument("--scp-bin", default=m8a.executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    args = parser.parse_args()

    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    tdrs = parse_tdrs(args.tdr)
    validate_tdrs(contract, tdrs)
    output = args.output_dir.resolve()
    if args.reuse_prepared:
        manifest = read_json(output / "sentaurus_matrix_manifest.json")
    else:
        manifest = prepare_sentaurus(contract, contract_path, tdrs, output)
        prepare_vela(contract, tdrs, output, args.materials.resolve(),
                     args.importer.resolve())
    vela_report = execute_vela(
        output, args.runner.resolve(), args.vela_jobs) if args.execute_vela else None
    if args.live_sentaurus:
        banner = run_sentaurus(
            manifest, output, args.ssh_target, args.ssh_bin, args.scp_bin,
            args.remote_root, args.sentaurus_jobs)
        extract_references(contract, manifest, output, banner)
    elif args.extract_existing:
        banner = (output / "sentaurus_banner.txt").read_text(
            encoding="utf-8").strip()
        extract_references(contract, manifest, output, banner)
    comparison = compare(contract, output) if args.compare else None
    print(json.dumps({
        "cases": len(manifest["cases"]),
        "vela_status": (vela_report or {}).get("status"),
        "comparison_status": (comparison or {}).get("status"),
    }))
    return 1 if vela_report is not None and vela_report["status"] != "accepted" else 0


if __name__ == "__main__":
    raise SystemExit(main())
