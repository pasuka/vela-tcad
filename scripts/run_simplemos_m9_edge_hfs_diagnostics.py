#!/usr/bin/env python3
"""Generate production-consistent per-edge HFS diagnostics for SimpleMOS."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess
import sys
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
import run_simplemos_m8a_model_ablation as m8a  # noqa: E402
import run_simplemos_m9_refdens_scan as m9  # noqa: E402


DEFAULT_CONTRACT = m9.DEFAULT_CONTRACT
DEFAULT_BASELINE = m9.DEFAULT_BASELINE
DEFAULT_OUTPUT = m9.DEFAULT_OUTPUT / "vela_edge_probes"
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner.exe"


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8", newline="\n")


def run(argv: Sequence[str]) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return completed.stdout


def gate_tag(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def exact_gate_points(stop: float) -> list[float]:
    count = int(round(stop / 0.05))
    points = [round(0.05 * index, 12) for index in range(count + 1)]
    if not math.isclose(points[-1], stop, abs_tol=1e-12):
        raise ValueError(f"gate target {stop:g} is not on the 0.05 V lattice")
    return points


def prepare_cases(contract: dict[str, Any], baseline: Path,
                  output: Path) -> list[dict[str, Any]]:
    cases = []
    for device in contract["edge_probe"]["devices"]:
        for drain_voltage in contract["edge_probe"]["drain_voltages_V"]:
            vd = float(drain_voltage)
            vd_tag = m8a.voltage_tag(vd)
            source = baseline / "vela" / device / "full" / "workflow" / f"vd_{vd_tag}"
            source_config = m9.read_json(source / "20_gate_sweep.json")
            drain_state = source / "10_drain_ramp_accepted_state.csv"
            final_state = source / "20_gate_sweep_accepted_state.csv"
            for gate_voltage in contract["edge_probe"]["gate_voltages_V"]:
                vg = float(gate_voltage)
                case = f"{device}_vd_{vd_tag}_vg_{gate_tag(vg)}"
                root = output / case
                root.mkdir(parents=True, exist_ok=True)
                needs_solve = vg not in (0.0, 2.5)
                state = drain_state if vg == 0.0 else final_state if vg == 2.5 else root / "accepted_state.csv"
                solve_config = None
                if needs_solve:
                    config = json.loads(json.dumps(source_config))
                    config["output_csv"] = str((root / "iv.csv").resolve())
                    config["log_file"] = str((root / "solve.log").resolve())
                    config["sweep"].update({
                        "start": 0.0,
                        "stop": vg,
                        "step": 0.05,
                        "bias_points": exact_gate_points(vg),
                        "initial_state_file": str(drain_state.resolve()),
                        "write_state_file": str(state.resolve()),
                        "write_vtk": False,
                    })
                    solve_config = root / "solve.json"
                    write_json(solve_config, config)
                probe = json.loads(json.dumps(source_config))
                probe.update({
                    "simulation_type": "edge_mobility_probe",
                    "state_file": str(state.resolve()),
                    "output_csv": str((root / "edge_hfs.csv").resolve()),
                })
                probe_config = root / "probe.json"
                write_json(probe_config, probe)
                cases.append({
                    "case": case,
                    "device": device,
                    "drain_voltage_V": vd,
                    "gate_voltage_V": vg,
                    "needs_solve": needs_solve,
                    "solve_config": str(solve_config.resolve()) if solve_config else None,
                    "state_file": str(state.resolve()),
                    "probe_config": str(probe_config.resolve()),
                    "edge_csv": str((root / "edge_hfs.csv").resolve()),
                })
    write_json(output / "edge_probe_manifest.json", {
        "schema": "vela.simplemos.sdevice.m9_edge_probe_manifest.v1",
        "status": "prepared",
        "cases": cases,
    })
    return cases


def execute_case(case: dict[str, Any], runner: Path) -> None:
    if case["needs_solve"]:
        run([str(runner), "--config", case["solve_config"]])
    if not Path(case["state_file"]).is_file():
        raise FileNotFoundError(f"missing accepted state for {case['case']}")
    status = json.loads(run([str(runner), "--config", case["probe_config"]]))
    if status.get("diagnostic_schema") != "vela.edge_mobility_probe.v2":
        raise RuntimeError(f"unexpected edge probe schema for {case['case']}")


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def summarize_case(case: dict[str, Any], required: list[str]) -> dict[str, Any]:
    with Path(case["edge_csv"]).open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    missing = [name for name in required if name not in rows[0]]
    if missing:
        raise ValueError(f"{case['case']} edge probe misses columns: {missing}")
    transport = [row for row in rows if float(
        row["electron_low_field_mobility_m2_V_s"]) > 0.0]
    e_drive = [float(row["electron_mobility_field_V_m"]) / 100.0
               for row in transport]
    h_drive = [float(row["hole_mobility_field_V_m"]) / 100.0
               for row in transport]
    e_ratio = [float(row["electron_high_field_ratio"]) for row in transport]
    h_ratio = [float(row["hole_high_field_ratio"]) for row in transport]
    e_limit = [float(row["electron_mobility_limiter"]) for row in transport]
    h_limit = [float(row["hole_mobility_limiter"]) for row in transport]
    e_error = [float(row["electron_limiter_reconstruction_error"])
               for row in transport]
    h_error = [float(row["hole_limiter_reconstruction_error"])
               for row in transport]
    e_aggregation = [float(row["electron_mean_mobility_aggregation_error"])
                     for row in transport]
    h_aggregation = [float(row["hole_mean_mobility_aggregation_error"])
                     for row in transport]
    return {
        "case": case["case"],
        "device": case["device"],
        "drain_voltage_V": case["drain_voltage_V"],
        "gate_voltage_V": case["gate_voltage_V"],
        "edge_count": len(rows),
        "transport_edge_count": len(transport),
        "maximum_electron_drive_V_per_cm": max(e_drive),
        "maximum_hole_drive_V_per_cm": max(h_drive),
        "median_electron_high_field_ratio": statistics.median(e_ratio),
        "p95_electron_high_field_ratio": percentile(e_ratio, 0.95),
        "maximum_electron_high_field_ratio": max(e_ratio),
        "median_hole_high_field_ratio": statistics.median(h_ratio),
        "p95_hole_high_field_ratio": percentile(h_ratio, 0.95),
        "maximum_hole_high_field_ratio": max(h_ratio),
        "minimum_electron_mobility_limiter": min(e_limit),
        "minimum_hole_mobility_limiter": min(h_limit),
        "maximum_electron_reconstruction_error": max(e_error),
        "maximum_hole_reconstruction_error": max(h_error),
        "maximum_electron_mean_mobility_aggregation_error": max(e_aggregation),
        "maximum_hole_mean_mobility_aggregation_error": max(h_aggregation),
        "edge_csv": m8a.portable_path(Path(case["edge_csv"])),
        "edge_csv_sha256": m8a.sha256(Path(case["edge_csv"])),
    }


def write_summary(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        raise ValueError("jobs must be at least one")
    contract = m9.read_json(args.contract.resolve())
    m9.validate_contract(contract)
    m9.validate_baseline_guard(contract)
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases = prepare_cases(contract, args.baseline_dir.resolve(), output)
    if args.prepare_only:
        print(json.dumps({"status": "prepared", "cases": len(cases)}))
        return 0
    runner = args.runner.resolve()
    if not runner.is_file():
        raise FileNotFoundError(f"missing Vela runner: {runner}")
    with ThreadPoolExecutor(max_workers=min(args.jobs, len(cases))) as executor:
        list(executor.map(lambda case: execute_case(case, runner), cases))
    required = list(contract["edge_probe"]["required_columns"])
    summaries = [summarize_case(case, required) for case in cases]
    write_summary(output / "edge_hfs_summary.csv", summaries)
    report = {
        "schema": "vela.simplemos.sdevice.m9_edge_hfs_diagnostics.v1",
        "status": "complete",
        "case_count": len(summaries),
        "intermediate_solve_count": sum(item["needs_solve"] for item in cases),
        "default_model_changed": False,
        "maximum_limiter_reconstruction_error": max(
            max(item["maximum_electron_reconstruction_error"],
                item["maximum_hole_reconstruction_error"])
            for item in summaries),
        "maximum_mean_mobility_aggregation_error": max(
            max(item["maximum_electron_mean_mobility_aggregation_error"],
                item["maximum_hole_mean_mobility_aggregation_error"])
            for item in summaries),
        "cases": summaries,
    }
    write_json(output / "edge_hfs_report.json", report)
    print(json.dumps({
        "status": report["status"],
        "cases": report["case_count"],
        "maximum_limiter_reconstruction_error": report[
            "maximum_limiter_reconstruction_error"],
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
