#!/usr/bin/env python3
"""Probe damped Newton settings for the two M8-B Eparallel failures."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m8a_model_ablation as m8a  # noqa: E402
import run_simplemos_m8b_hfs_controls as m8b  # noqa: E402


DEFAULT_CONTRACT = (
    REPO / "reference_tcad/simplemos_sentaurus2022"
    / "simplemos_m8b_eparallel_convergence_contract_v1.json"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/simplemos_sentaurus2022"
    / "m8b_eparallel_convergence"
)
DEFAULT_RAW = REPO / "build-release/m8b_ec_raw"
DEFAULT_REMOTE_ROOT = (
    "~/sentaurus_runs/vela_oracle/simplemos_m8b_eparallel_convergence_20260828"
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


def extract_archive(archive: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as stream:
        stream.extractall(destination, filter="data")


def validate_contract(contract: dict[str, Any]) -> None:
    if contract["devices"] != ["n17", "n21"]:
        raise ValueError("convergence probe devices must be n17 and n21")
    if not math.isclose(float(contract["drain_voltage_V"]), 1.0):
        raise ValueError("convergence probe drain voltage must be 1 V")
    damping = [float(item["line_search_damping"])
               for item in contract["probes"]]
    if damping != [0.5, 0.1, 0.01]:
        raise ValueError("convergence probe damping order is frozen")
    lattice = contract["gate_lattice"]
    if lattice["point_count"] != 51 or lattice["interpolation"] != "forbidden":
        raise ValueError("convergence probe must retain the exact 51-point lattice")


def damped_deck(case: str, variant: dict[str, Any], damping: float,
                numerics: dict[str, Any]) -> str:
    retry = {
        "coupled_iterations": int(numerics["coupled_iterations"]),
        "gate_min_step": float(numerics["gate_min_step"]),
    }
    text = m8b.sentaurus_deck(case, 1.0, variant, 50, retry)
    target = "      Coupled { Poisson Electron Hole }\n      CurrentPlot"
    replacement = (
        "      Coupled(Iterations="
        f"{retry['coupled_iterations']} LineSearchDamping={damping:g}) "
        "{ Poisson Electron Hole }\n      CurrentPlot"
    )
    if target not in text:
        raise ValueError("gate-sweep Coupled block was not found")
    return text.replace(target, replacement, 1)


def prepare(contract: dict[str, Any], contract_path: Path,
            output: Path) -> dict[str, Any]:
    bundle = output / "sentaurus_bundle"
    variant = {
        "id": "eparallel",
        "hfs_option": "HighFieldSaturation(Eparallel)",
        "math_options": [],
    }
    cases = []
    for device in contract["devices"]:
        device_dir = bundle / device
        device_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(m8b.DEFAULT_TDRS[device].resolve(),
                     device_dir / "input_fps.tdr")
        for probe in contract["probes"]:
            probe_id = str(probe["id"])
            case = f"{device}_eparallel_vd_1_{probe_id}"
            case_dir = device_dir / probe_id / case
            case_dir.mkdir(parents=True, exist_ok=True)
            deck = case_dir / f"{case}_des.cmd"
            text = damped_deck(
                case, variant, float(probe["line_search_damping"]),
                contract["common_numerics"]).replace(
                    'Grid="input_fps.tdr"', 'Grid="../../input_fps.tdr"')
            deck.write_text(text, encoding="utf-8", newline="\n")
            cases.append({
                "case": case, "device": device, "probe": probe_id,
                "line_search_damping": float(probe["line_search_damping"]),
                "deck": m8a.portable_path(deck),
                "deck_sha256": m8a.sha256(deck),
            })
    manifest = {
        "schema": "vela.simplemos.sdevice.m8b_eparallel_convergence_manifest.v1",
        "status": "prepared",
        "contract": m8a.portable_path(contract_path),
        "contract_sha256": m8a.sha256(contract_path),
        "cases": cases,
    }
    write_json(output / "sentaurus_matrix_manifest.json", manifest)
    return manifest


def execute_remote(manifest: dict[str, Any], output: Path, ssh_target: str,
                   ssh_bin: str, scp_bin: str, remote_root: str,
                   jobs: int, raw_extract: Path) -> str:
    if jobs < 1:
        raise ValueError("Sentaurus jobs must be at least one")
    banner = run([ssh_bin, "-n", ssh_target,
                  "sdevice -h 2>&1 | sed -n '1,5p'"], capture=True).strip()
    if "T-2022.03-SP2" not in banner:
        raise RuntimeError(f"unexpected Sentaurus release:\n{banner}")
    run([ssh_bin, "-n", ssh_target, f"mkdir -p {remote_root}"])
    run([scp_bin, "-r", str(output / "sentaurus_bundle"),
         f"{ssh_target}:{remote_root}/"])

    def execute(index_case: tuple[int, dict[str, Any]]) -> None:
        index, item = index_case
        case = str(item["case"])
        remote_dir = (f"{remote_root}/sentaurus_bundle/{item['device']}/"
                      f"{item['probe']}/{case}")
        print(f"[{index}/{len(manifest['cases'])}] Sentaurus {case}", flush=True)
        command = (
            f"set -u; cd {remote_dir}; "
            f"sdevice {case}_des.cmd > {case}.console.log 2>&1; "
            f"code=$?; printf '%s\\n' $code > {case}.exit_code.txt; exit 0"
        )
        run([ssh_bin, "-n", ssh_target, command])

    indexed = list(enumerate(manifest["cases"], start=1))
    with ThreadPoolExecutor(max_workers=min(jobs, len(indexed))) as executor:
        list(executor.map(execute, indexed))
    archive_name = "simplemos_m8b_eparallel_convergence_results.tgz"
    run([ssh_bin, "-n", ssh_target,
         f"cd {remote_root} && tar -czf {archive_name} sentaurus_bundle"])
    download = output / "sentaurus_raw"
    download.mkdir(parents=True, exist_ok=True)
    archive = download / archive_name
    run([scp_bin, f"{ssh_target}:{remote_root}/{archive_name}", str(archive)])
    extract_archive(archive, raw_extract)
    (output / "sentaurus_banner.txt").write_text(
        banner + "\n", encoding="utf-8", newline="\n")
    return banner


def summarize(contract: dict[str, Any], manifest: dict[str, Any],
              output: Path, raw_extract: Path) -> dict[str, Any]:
    cases = []
    curves: dict[str, list[dict[str, float]]] = {}
    lattice = [round(index * 0.05, 12) for index in range(51)]
    step_pattern = re.compile(r"Step-size less than MinStep \(step-size = ([^)]+)\)")
    t_pattern = re.compile(r"Computing step from t=([0-9.eE+-]+) to t=([0-9.eE+-]+)")
    for item in manifest["cases"]:
        case = str(item["case"])
        case_dir = (raw_extract / "sentaurus_bundle"
                    / item["device"] / item["probe"] / case)
        console_path = case_dir / f"{case}.console.log"
        exit_path = case_dir / f"{case}.exit_code.txt"
        console = console_path.read_text(encoding="utf-8", errors="replace")
        exit_code = int(exit_path.read_text(encoding="utf-8").strip())
        steps = t_pattern.findall(console)
        minimum = step_pattern.findall(console)
        success = (exit_code == 0 and "Exit due to failure" not in console
                   and "Good Bye" in console)
        exact_rows: list[dict[str, float]] = []
        if success:
            plot = case_dir / f"IdVg_{case}_des.plt"
            try:
                exact_rows = m8b.m4.exact_plt_curve(plot, lattice, 1e-10)
            except (ValueError, OSError):
                success = False
        if success:
            curves[case] = exact_rows
        cases.append({
            **item, "exit_code": exit_code,
            "success": success,
            "exact_gate_point_count": len(exact_rows),
            "last_continuation_t": (float(steps[-1][0]) if steps else None),
            "failed_step_size": (minimum[-1] if minimum else None),
            "console_log": m8a.portable_path(console_path),
            "console_log_sha256": m8a.sha256(console_path),
        })
    probe_summaries = []
    for probe in contract["probes"]:
        selected = [item for item in cases if item["probe"] == probe["id"]]
        probe_summaries.append({
            "probe": probe["id"],
            "line_search_damping": float(probe["line_search_damping"]),
            "successful_devices": [item["device"] for item in selected
                                   if item["success"]],
            "success_count": sum(item["success"] for item in selected),
        })
    passing = [item for item in probe_summaries
               if item["success_count"] == len(contract["devices"])]
    selected_probe = (max(passing, key=lambda item: item["line_search_damping"])
                      if passing else None)
    selected_by_device = {}
    for device in contract["devices"]:
        successful = [item for item in cases
                      if item["device"] == device and item["success"]]
        selected_by_device[device] = (
            max(successful, key=lambda item: item["line_search_damping"])["probe"]
            if successful else None)
    robustness_pairs = []
    for device in contract["devices"]:
        successful = sorted(
            [item for item in cases
             if item["device"] == device and item["success"]],
            key=lambda item: item["line_search_damping"], reverse=True)
        for left_index, left in enumerate(successful):
            for right in successful[left_index + 1:]:
                left_rows = curves[left["case"]]
                right_rows = curves[right["case"]]
                log_differences = []
                relative_differences = []
                for left_row, right_row in zip(left_rows, right_rows):
                    left_current = abs(left_row[
                        "drain_total_current_A_per_um"])
                    right_current = abs(right_row[
                        "drain_total_current_A_per_um"])
                    log_differences.append(abs(math.log10(
                        max(left_current, 1e-300) / max(right_current, 1e-300))))
                    relative_differences.append(
                        abs(left_current - right_current)
                        / max(left_current, right_current, 1e-300))
                robustness_pairs.append({
                    "device": device,
                    "left_probe": left["probe"],
                    "right_probe": right["probe"],
                    "maximum_absolute_log10_current_difference_dex": max(
                        log_differences),
                    "maximum_relative_current_difference": max(
                        relative_differences),
                })
    report = {
        "schema": "vela.simplemos.sdevice.m8b_eparallel_convergence_report.v1",
        "status": "complete",
        "selection_policy": contract["selection_policy"],
        "selected_probe": selected_probe,
        "selected_by_device": selected_by_device,
        "robustness_pairs": robustness_pairs,
        "probe_summaries": probe_summaries,
        "cases": cases,
    }
    write_json(output / "convergence_report.json", report)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--raw-dir", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--reuse-prepared", action="store_true")
    parser.add_argument("--live-sentaurus", action="store_true")
    parser.add_argument("--summarize-existing", action="store_true")
    parser.add_argument("--sentaurus-jobs", type=int, default=4)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=m8a.executable("ssh"))
    parser.add_argument("--scp-bin", default=m8a.executable("scp"))
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    args = parser.parse_args()
    contract_path = args.contract.resolve()
    contract = read_json(contract_path)
    validate_contract(contract)
    output = args.output_dir.resolve()
    raw_extract = args.raw_dir.resolve()
    manifest = (read_json(output / "sentaurus_matrix_manifest.json")
                if args.reuse_prepared else
                prepare(contract, contract_path, output))
    if args.live_sentaurus:
        execute_remote(manifest, output, args.ssh_target, args.ssh_bin,
                       args.scp_bin, args.remote_root, args.sentaurus_jobs,
                       raw_extract)
    elif args.summarize_existing:
        archive = (output / "sentaurus_raw"
                   / "simplemos_m8b_eparallel_convergence_results.tgz")
        if archive.is_file():
            extract_archive(archive, raw_extract)
    report = (summarize(contract, manifest, output, raw_extract)
              if args.live_sentaurus or args.summarize_existing else None)
    print(json.dumps({
        "cases": len(manifest["cases"]),
        "selected_probe": (report or {}).get("selected_probe"),
        "selected_by_device": (report or {}).get("selected_by_device"),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
