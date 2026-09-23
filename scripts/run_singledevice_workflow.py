#!/usr/bin/env python3
"""Materialize and optionally run the original SingleDevice Save/Load workflow.

The workflow preserves the dependency graph of the Sentaurus deck:

  equilibrium -> save
      save -> linear drain ramp -> linear Id-Vg
      save -> saturation drain ramp -> saturation Id-Vg

Extra gate-bias points may be inserted for convergence, but every committed
Sentaurus reference bias is solved directly and marked in the manifest.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
from typing import Any


INPUT_PATH_KEYS = ("mesh_file", "node_doping_file", "materials_file")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def reference_biases(path: Path) -> list[float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path}: reference curve is empty")
    key = next(
        (candidate for candidate in ("gate_voltage_V", "bias_V", "bias")
         if candidate in rows[0]), "")
    if not key:
        key = next((name for name in rows[0] if "voltage" in name.lower()), "")
    if not key:
        raise ValueError(f"{path}: cannot identify gate-bias column")
    values = [float(row[key]) for row in rows]
    if values != sorted(values) or len(set(values)) != len(values):
        raise ValueError(f"{path}: reference biases must be unique and increasing")
    return values


def expanded_biases(reference: list[float], subdivisions: int) -> list[float]:
    if subdivisions < 1:
        raise ValueError("subdivisions must be positive")
    expanded = [reference[0]]
    for left, right in zip(reference, reference[1:]):
        for index in range(1, subdivisions + 1):
            expanded.append(left + (right - left) * index / subdivisions)
    return [round(value, 14) for value in expanded]


def absolutize_inputs(config: dict[str, Any], base_dir: Path) -> None:
    for key in INPUT_PATH_KEYS:
        value = config.get(key)
        if value and not Path(value).is_absolute():
            config[key] = str((base_dir / value).resolve())


def contact_bias(config: dict[str, Any], name: str, value: float) -> None:
    for contact in config["contacts"]:
        if contact["name"] == name:
            contact["bias"] = value
            return
    raise ValueError(f"base configuration has no {name!r} contact")


def diagnostics(sweep: dict[str, Any], run_dir: Path, stem: str) -> None:
    sweep["diagnostics"] = {
        "transport": {"enabled": True},
        "contact_edge": {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": str((run_dir / f"{stem}_contact_edges.csv").resolve()),
        },
        "terminal_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": str((run_dir / f"{stem}_terminal_balance.csv").resolve()),
        },
    }
    sweep["write_state_every_point_prefix"] = str(
        (run_dir / f"{stem}_state").resolve())


def stage_config(base: dict[str, Any], base_dir: Path, run_dir: Path,
                 stem: str) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    absolutize_inputs(config, base_dir)
    config["output_csv"] = str((run_dir / f"{stem}.csv").resolve())
    config["log_file"] = str((run_dir / f"{stem}.log").resolve())
    config["sweep"].pop("initial_state_file", None)
    config["sweep"]["write_state_file"] = str(
        (run_dir / f"{stem}_final_state.h5").resolve())
    config["sweep"]["write_vtk"] = False
    diagnostics(config["sweep"], run_dir, stem)
    return config


def materialize(linear_base_path: Path, saturation_base_path: Path,
                linear_reference_path: Path, saturation_reference_path: Path,
                run_dir: Path) -> dict[str, Any]:
    linear_base = load_json(linear_base_path)
    saturation_base = load_json(saturation_base_path)
    for key in INPUT_PATH_KEYS:
        if linear_base.get(key) != saturation_base.get(key):
            raise ValueError(f"linear and saturation bases disagree on {key}")
    linear_biases = reference_biases(linear_reference_path)
    saturation_reference_biases = reference_biases(saturation_reference_path)
    if linear_biases != saturation_reference_biases:
        raise ValueError("linear and saturation reference bias lattices differ")
    saturation_biases = expanded_biases(saturation_reference_biases, 2)

    run_dir.mkdir(parents=True, exist_ok=True)
    common_state = (run_dir / "common_saved_equilibrium_state.h5").resolve()

    equilibrium = stage_config(
        linear_base, linear_base_path.parent, run_dir, "00_equilibrium")
    contact_bias(equilibrium, "gate", linear_biases[0])
    contact_bias(equilibrium, "drain", 0.0)
    equilibrium["sweep"].update({
        "contact": "gate", "current_contact": "drain",
        "start": linear_biases[0], "stop": linear_biases[0], "step": 0.01,
        "bias_points": [linear_biases[0]],
        "initialization": {"mode": "poisson_block"},
        "write_state_file": str(common_state),
    })
    # The source deck gives its first Coupled stage Iterations=100.  Curve
    # configurations are often tightened to 20 after a restart, so restore
    # the source-stage budget for cold equilibrium and drain pre-bias ramps.
    equilibrium["solver"]["max_iter"] = 100
    # The cold Poisson-block state can be several thermal voltages away from
    # coupled DD equilibrium.  The 0.1 V limiter also applies to the linear
    # gate sweep below; the proven saturation and drain-ramp settings remain
    # at their imported 0.025 V value.
    equilibrium["solver"]["quasi_fermi_update_limit_V"] = 0.1

    linear_ramp = stage_config(
        linear_base, linear_base_path.parent, run_dir, "10_linear_drain_ramp")
    contact_bias(linear_ramp, "gate", linear_biases[0])
    contact_bias(linear_ramp, "drain", 0.0)
    linear_ramp["sweep"].update({
        "contact": "drain", "current_contact": "drain",
        "start": 0.0, "stop": 0.1, "step": 0.01,
        "bias_points": [round(index * 0.01, 14) for index in range(11)],
        "initial_state_file": str(common_state),
    })
    linear_ramp["solver"]["max_iter"] = 100

    linear_curve = stage_config(
        linear_base, linear_base_path.parent, run_dir, "11_linear_idvg")
    contact_bias(linear_curve, "gate", linear_biases[0])
    contact_bias(linear_curve, "drain", 0.1)
    linear_curve["sweep"].update({
        "contact": "gate", "current_contact": "drain",
        "start": linear_biases[0], "stop": linear_biases[-1],
        "step": linear_biases[1] - linear_biases[0],
        "bias_points": linear_biases,
        "initial_state_file": linear_ramp["sweep"]["write_state_file"],
    })
    linear_curve["solver"]["quasi_fermi_update_limit_V"] = 0.1

    saturation_ramp = stage_config(
        saturation_base, saturation_base_path.parent, run_dir,
        "20_saturation_drain_ramp")
    contact_bias(saturation_ramp, "gate", saturation_reference_biases[0])
    contact_bias(saturation_ramp, "drain", 0.0)
    saturation_ramp_biases = [round(index * 0.05, 14) for index in range(23)]
    saturation_ramp["sweep"].update({
        "contact": "drain", "current_contact": "drain",
        "start": 0.0, "stop": 1.1, "step": 0.05,
        "bias_points": saturation_ramp_biases,
        "initial_state_file": str(common_state),
    })
    saturation_ramp["solver"]["max_iter"] = 100

    saturation_curve = stage_config(
        saturation_base, saturation_base_path.parent, run_dir,
        "21_saturation_idvg")
    contact_bias(saturation_curve, "gate", saturation_reference_biases[0])
    contact_bias(saturation_curve, "drain", 1.1)
    saturation_curve["sweep"].update({
        "contact": "gate", "current_contact": "drain",
        "start": saturation_biases[0], "stop": saturation_biases[-1],
        "step": saturation_biases[1] - saturation_biases[0],
        "bias_points": saturation_biases,
        "initial_state_file": saturation_ramp["sweep"]["write_state_file"],
    })

    configs = {
        "equilibrium": equilibrium,
        "linear_drain_ramp": linear_ramp,
        "linear_idvg": linear_curve,
        "saturation_drain_ramp": saturation_ramp,
        "saturation_idvg": saturation_curve,
    }
    stages: list[dict[str, Any]] = []
    dependencies = {
        "equilibrium": [],
        "linear_drain_ramp": ["equilibrium"],
        "linear_idvg": ["linear_drain_ramp"],
        "saturation_drain_ramp": ["equilibrium"],
        "saturation_idvg": ["saturation_drain_ramp"],
    }
    for index, (name, config) in enumerate(configs.items()):
        config_path = (run_dir / f"{index:02d}_{name}.json").resolve()
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        stages.append({
            "name": name,
            "config": str(config_path),
            "depends_on": dependencies[name],
            "initial_state_file": config["sweep"].get("initial_state_file"),
            "final_state_file": config["sweep"]["write_state_file"],
            "bias_points": config["sweep"]["bias_points"],
        })
    manifest = {
        "schema": "vela.singledevice.workflow.v1",
        "status": "materialized",
        "dependency_contract": "sentaurus_save_load_two_branch",
        "common_saved_state": str(common_state),
        "reference_gate_biases_V": linear_biases,
        "saturation_auxiliary_biases_V": [
            value for value in saturation_biases if value not in linear_biases],
        "stages": stages,
    }
    (run_dir / "workflow_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def comparison_csv(reference_path: Path, curve_path: Path, output: Path) -> None:
    with reference_path.open(newline="", encoding="utf-8") as handle:
        reference_rows = list(csv.DictReader(handle))
    with curve_path.open(newline="", encoding="utf-8") as handle:
        curve_rows = list(csv.DictReader(handle))
    if not reference_rows or not curve_rows:
        raise ValueError("reference and Vela curves must be non-empty")
    reference_bias_key = next(
        key for key in ("gate_voltage_V", "bias_V", "bias")
        if key in reference_rows[0])
    reference_current_key = next(
        key for key in ("current_total_A_per_um", "current_total")
        if key in reference_rows[0])
    curve_bias_key = "bias_V"
    curve_current_key = "current_total_A_per_um"
    output_rows = []
    for reference in reference_rows:
        bias = float(reference[reference_bias_key])
        matches = [row for row in curve_rows
                   if abs(float(row[curve_bias_key]) - bias) <= 1.0e-10
                   and row.get("converged") == "1"]
        if len(matches) != 1:
            raise ValueError(f"expected one converged Vela point at {bias} V")
        output_rows.append({
            "gate_voltage_V": bias,
            "sentaurus_current_A_per_um": abs(float(reference[reference_current_key])),
            "vela_current_A_per_um": abs(float(matches[0][curve_current_key])),
        })
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)


def execute(manifest: dict[str, Any], runner: Path, run_dir: Path,
            linear_reference: Path, saturation_reference: Path) -> dict[str, Any]:
    for stage in manifest["stages"]:
        if stage.pop("resume_skip", False):
            stage["status"] = "pass"
            stage["execution"] = "reused_verified"
            state = Path(stage["final_state_file"])
            stage["final_state_sha256"] = sha256(state)
            continue
        process = subprocess.run(
            [str(runner.resolve()), "--config", stage["config"]],
            cwd=run_dir, text=True, capture_output=True)
        (run_dir / f"{stage['name']}.console.log").write_text(
            process.stdout + process.stderr, encoding="utf-8")
        stage["returncode"] = process.returncode
        stage["status"] = "pass" if process.returncode == 0 else "fail"
        stage["execution"] = "executed"
        state = Path(stage["final_state_file"])
        if state.exists():
            stage["final_state_sha256"] = sha256(state)
        if process.returncode != 0:
            manifest["status"] = "fail"
            break
    else:
        manifest["status"] = "pass"
        comparisons = {
            "linear": run_dir / "linear_comparison.csv",
            "saturation": run_dir / "saturation_comparison.csv",
        }
        comparison_csv(
            linear_reference, run_dir / "11_linear_idvg.csv", comparisons["linear"])
        comparison_csv(
            saturation_reference, run_dir / "21_saturation_idvg.csv",
            comparisons["saturation"])
        manifest["comparisons"] = {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in comparisons.items()
        }
    common = Path(manifest["common_saved_state"])
    if common.exists():
        manifest["common_saved_state_sha256"] = sha256(common)
    for stage in manifest["stages"]:
        stage.pop("resume_skip", None)
    (run_dir / "workflow_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linear-base", type=Path, required=True)
    parser.add_argument("--saturation-base", type=Path, required=True)
    parser.add_argument("--linear-reference", type=Path, required=True)
    parser.add_argument("--saturation-reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument(
        "--resume", action="store_true",
        help="reuse passed stages only when their config hash and final state match")
    args = parser.parse_args()
    if args.execute and args.runner is None:
        parser.error("--execute requires --runner")
    if args.clean and args.resume:
        parser.error("--clean and --resume are mutually exclusive")
    previous = None
    previous_config_hashes: dict[str, str] = {}
    previous_manifest_path = args.output_dir.resolve() / "workflow_manifest.json"
    if args.resume:
        if not previous_manifest_path.is_file():
            parser.error("--resume requires an existing workflow_manifest.json")
        previous = load_json(previous_manifest_path)
        for stage in previous.get("stages", []):
            config_path = Path(stage["config"])
            if config_path.is_file():
                previous_config_hashes[stage["name"]] = sha256(config_path)
    if args.clean and args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    manifest = materialize(
        args.linear_base.resolve(), args.saturation_base.resolve(),
        args.linear_reference.resolve(), args.saturation_reference.resolve(),
        args.output_dir.resolve())
    previous_by_name = {
        stage["name"]: stage for stage in (previous or {}).get("stages", [])}
    for stage in manifest["stages"]:
        config_hash = sha256(Path(stage["config"]))
        stage["config_sha256"] = config_hash
        prior = previous_by_name.get(stage["name"], {})
        state = Path(stage["final_state_file"])
        stage["resume_skip"] = bool(
            args.resume and prior.get("status") == "pass" and state.is_file()
            and previous_config_hashes.get(stage["name"]) == config_hash)
    if args.resume:
        manifest["resume"] = {
            "enabled": True,
            "reused_stages": [stage["name"] for stage in manifest["stages"]
                              if stage["resume_skip"]],
        }
    if args.execute:
        manifest = execute(
            manifest, args.runner, args.output_dir.resolve(),
            args.linear_reference.resolve(), args.saturation_reference.resolve())
    else:
        previous_manifest_path.write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "stages": [
        {"name": stage["name"], "status": stage.get("status", "materialized")}
        for stage in manifest["stages"]]}))
    return 0 if manifest["status"] in {"materialized", "pass"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
