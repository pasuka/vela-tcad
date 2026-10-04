#!/usr/bin/env python3
"""Materialize and optionally execute the SimpleMOS M3 restart workflow.

This script is intentionally limited to SDevice orchestration.  It consumes a
Vela device-simulation base deck and the frozen SimpleMOS validation contract;
it neither runs nor models SProcess.  Every drain and gate stage can start only
from the immediately preceding stage's accepted, hash-identified state.
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


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = (
    REPO / "reference_tcad" / "simplemos_sentaurus2022"
    / "simplemos_sdevice_validation_contract_v1.json"
)
INPUT_PATH_KEYS = ("mesh_file", "node_doping_file", "materials_file")
REQUIRED_CONTACTS = ("source", "drain", "gate", "substrate")
REQUIRED_PHASES = ("equilibrium", "drain_ramp", "gate_sweep")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8", newline="\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def voltage_tag(voltage: float) -> str:
    text = format(voltage, ".12g").replace("-", "m").replace(".", "p")
    return f"vd_{text}"


def physical_step_control(normalized: dict[str, Any], start: float,
                          stop: float) -> dict[str, float]:
    """Convert Sentaurus normalized sweep-time controls to physical volts."""
    span = abs(stop - start)
    if not math.isfinite(span) or span <= 0.0:
        raise ValueError("quasistationary voltage span must be finite and positive")

    converted: dict[str, float] = {"voltage_span_V": span}
    for source, target in (
        ("InitialStep", "initial_step_V"),
        ("MinStep", "min_step_V"),
        ("MaxStep", "max_step_V"),
    ):
        value = float(normalized[source])
        if not math.isfinite(value) or value <= 0.0 or value > 1.0:
            raise ValueError(f"{source} must be in the normalized interval (0, 1]")
        converted[target] = span * value
    increment = float(normalized["Increment"])
    if not math.isfinite(increment) or increment < 1.0:
        raise ValueError("Increment must be finite and at least 1")
    converted["growth_factor"] = increment
    return converted


def set_contact_bias(config: dict[str, Any], name: str, value: float) -> None:
    for contact in config["contacts"]:
        if str(contact.get("name", "")).lower() == name:
            contact["bias"] = value
            return
    raise ValueError(f"base configuration has no {name!r} contact")


def validate_base(base: dict[str, Any]) -> None:
    if base.get("simulation_type") != "dc_sweep":
        raise ValueError("SimpleMOS M3 base configuration must be a dc_sweep deck")
    contacts = [str(item.get("name", "")).lower()
                for item in base.get("contacts", [])]
    if len(contacts) != len(set(contacts)):
        raise ValueError("base configuration contact names must be unique")
    if set(contacts) != set(REQUIRED_CONTACTS):
        raise ValueError(
            f"base configuration contacts {sorted(contacts)} do not match "
            f"required SimpleMOS contacts {sorted(REQUIRED_CONTACTS)}")
    if not isinstance(base.get("solver"), dict):
        raise ValueError("base configuration must define solver settings")
    if not isinstance(base.get("sweep"), dict):
        raise ValueError("base configuration must define sweep settings")


def prepare_config(base: dict[str, Any], base_dir: Path, run_dir: Path,
                   stem: str) -> dict[str, Any]:
    config = json.loads(json.dumps(base))
    for key in INPUT_PATH_KEYS:
        value = config.get(key)
        if value and not Path(value).is_absolute():
            config[key] = str((base_dir / value).resolve())
    for name in REQUIRED_CONTACTS:
        set_contact_bias(config, name, 0.0)
    config["output_csv"] = str((run_dir / f"{stem}.csv").resolve())
    config["log_file"] = str((run_dir / f"{stem}.log").resolve())
    config["sweep"] = {
        "mode": "iv",
        "contact": "gate",
        "current_contact": "drain",
        "start": 0.0,
        "stop": 0.0,
        "step": 0.01,
        "bias_points": [0.0],
        "warm_start": True,
        "write_vtk": False,
        "write_state_file": str((run_dir / f"{stem}_accepted_state.csv").resolve()),
    }
    config["simplemos_m3"] = {
        "previous_accepted_state_only": True,
        "pointwise_reclosure_forbidden": True,
    }
    return config


def configure_sweep(config: dict[str, Any], contact: str, start: float,
                    stop: float, normalized: dict[str, Any]) -> dict[str, float]:
    physical = physical_step_control(normalized, start, stop)
    sweep = config["sweep"]
    sweep.update({
        "contact": contact,
        "current_contact": "drain",
        "start": start,
        "stop": stop,
        "step": math.copysign(physical["initial_step_V"], stop - start),
        "min_step": physical["min_step_V"],
        "max_step": physical["max_step_V"],
        "growth_factor": physical["growth_factor"],
    })
    sweep.pop("bias_points", None)
    return physical


def validate_manifest_chain(manifest: dict[str, Any]) -> None:
    stages = manifest.get("stages", [])
    by_id: dict[str, dict[str, Any]] = {}
    previous_by_branch: dict[str, str] = {}
    phases_by_branch: dict[str, list[str]] = {}
    for stage in stages:
        stage_id = str(stage["id"])
        if stage_id in by_id:
            raise ValueError(f"duplicate workflow stage id: {stage_id}")
        branch = str(stage["branch"])
        phases_by_branch.setdefault(branch, []).append(str(stage["phase"]))
        expected = previous_by_branch.get(branch)
        actual = stage.get("predecessor")
        if actual != expected:
            raise ValueError(
                f"{stage_id}: predecessor {actual!r} is not the immediately "
                f"previous branch stage {expected!r}")
        initial = stage.get("initial_state_file")
        if expected is None:
            if initial is not None:
                raise ValueError(f"{stage_id}: first stage cannot use a restart state")
        else:
            expected_state = by_id[expected]["final_state_file"]
            if initial != expected_state:
                raise ValueError(
                    f"{stage_id}: restart state is not predecessor {expected}'s "
                    "accepted state")
        by_id[stage_id] = stage
        previous_by_branch[branch] = stage_id
    for branch, phases in phases_by_branch.items():
        if phases != list(REQUIRED_PHASES):
            raise ValueError(
                f"{branch}: workflow phases {phases} do not match "
                f"required sequence {list(REQUIRED_PHASES)}")


def materialize(base_path: Path, contract_path: Path, output_dir: Path,
                drain_voltages: list[float]) -> dict[str, Any]:
    base = read_json(base_path)
    contract = read_json(contract_path)
    validate_base(base)
    if not drain_voltages:
        raise ValueError("at least one drain voltage is required")
    if len(set(drain_voltages)) != len(drain_voltages):
        raise ValueError("drain voltages must be unique")
    if any(not math.isfinite(value) or value <= 0.0 for value in drain_voltages):
        raise ValueError("drain voltages must be finite and positive")

    solve = contract["solve"]
    drain_normalized = solve["drain_ramp"]["normalized_step_control"]
    gate_spec = solve["gate_sweep"]
    gate_start = float(gate_spec["start_voltage_V"])
    gate_stop = float(gate_spec["stop_voltage_V"])
    gate_normalized = gate_spec["normalized_step_control"]
    output_dir.mkdir(parents=True, exist_ok=True)

    stages: list[dict[str, Any]] = []
    for vd in drain_voltages:
        branch = voltage_tag(vd)
        branch_dir = (output_dir / branch).resolve()
        branch_dir.mkdir(parents=True, exist_ok=True)

        equilibrium = prepare_config(
            base, base_path.parent, branch_dir, "00_equilibrium")
        equilibrium["sweep"]["initialization"] = {"mode": "poisson_block"}
        equilibrium["solver"]["max_iter"] = max(
            int(equilibrium["solver"].get("max_iter", 0)),
            int(solve["initial"][0]["iterations"]),
        )

        drain = prepare_config(base, base_path.parent, branch_dir, "10_drain_ramp")
        drain_physical = configure_sweep(drain, "drain", 0.0, vd, drain_normalized)
        drain["solver"]["max_iter"] = max(
            int(drain["solver"].get("max_iter", 0)),
            int(solve["initial"][0]["iterations"]),
        )

        gate = prepare_config(base, base_path.parent, branch_dir, "20_gate_sweep")
        set_contact_bias(gate, "drain", vd)
        gate_physical = configure_sweep(
            gate, "gate", gate_start, gate_stop, gate_normalized)

        branch_configs = [
            ("equilibrium", equilibrium, None, None),
            ("drain_ramp", drain, "equilibrium", drain_physical),
            ("gate_sweep", gate, "drain_ramp", gate_physical),
        ]
        ids: dict[str, str] = {}
        for index, (phase, config, predecessor_phase, conversion) in enumerate(
                branch_configs):
            stage_id = f"{branch}.{phase}"
            ids[phase] = stage_id
            predecessor = ids.get(predecessor_phase) if predecessor_phase else None
            if predecessor is not None:
                previous = stages[-1]
                config["sweep"]["initial_state_file"] = previous["final_state_file"]
            config["simplemos_m3"].update({
                "branch": branch,
                "phase": phase,
                "predecessor": predecessor,
            })
            config_path = (branch_dir / f"{index * 10:02d}_{phase}.json").resolve()
            write_json(config_path, config)
            stage: dict[str, Any] = {
                "id": stage_id,
                "branch": branch,
                "phase": phase,
                "predecessor": predecessor,
                "config": str(config_path),
                "config_sha256": sha256(config_path),
                "output_csv": config["output_csv"],
                "initial_state_file": config["sweep"].get("initial_state_file"),
                "final_state_file": config["sweep"]["write_state_file"],
                "status": "materialized",
            }
            if conversion is not None:
                normalized = (
                    drain_normalized if phase == "drain_ramp" else gate_normalized)
                stage["vela_physical_step_control"] = conversion
                # Keep the source object immutable in the generated manifest.
                stage["sentaurus_normalized_step_control"] = json.loads(
                    json.dumps(normalized))
            stages.append(stage)

    manifest = {
        "schema": "vela.simplemos.sdevice.m3_workflow.v1",
        "status": "materialized",
        "scope": "sdevice_only",
        "base_config": str(base_path.resolve()),
        "base_config_sha256": sha256(base_path),
        "validation_contract": str(contract_path.resolve()),
        "validation_contract_sha256": sha256(contract_path),
        "drain_voltages_V": drain_voltages,
        "rules": {
            "previous_accepted_state_only": True,
            "pointwise_reclosure_forbidden": True,
            "sentaurus_steps_are_normalized_sweep_time": True,
            "vela_steps_are_physical_voltage": True,
        },
        "physics_readiness": {
            "inherited_from_base_config": True,
            "original_deck_phumob_status": contract["frontend_contract"]["phumob_status"],
            "workflow_acceptance_is_not_physics_parity": True,
        },
        "stages": stages,
    }
    validate_manifest_chain(manifest)
    write_json(output_dir / "workflow_manifest.json", manifest)
    return manifest


def accepted_terminal_output(stage: dict[str, Any]) -> tuple[bool, str]:
    output = Path(stage["output_csv"])
    if not output.is_file():
        return False, "missing stage output CSV"
    with output.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return False, "stage output CSV has no data rows"
    if "converged" not in rows[0] or "bias_V" not in rows[0]:
        return False, "stage output CSV lacks bias_V or converged"
    if any(row["converged"] != "1" for row in rows):
        return False, "stage contains a rejected bias point"
    config = read_json(Path(stage["config"]))
    expected = float(config["sweep"]["stop"])
    actual = float(rows[-1]["bias_V"])
    tolerance = 1.0e-10 * max(1.0, abs(expected))
    if abs(actual - expected) > tolerance:
        return False, f"terminal bias {actual} V does not reach {expected} V"
    return True, "accepted_terminal_state"


def execute(manifest: dict[str, Any], runner: Path, output_dir: Path) -> dict[str, Any]:
    validate_manifest_chain(manifest)
    by_id: dict[str, dict[str, Any]] = {}
    for stage in manifest["stages"]:
        predecessor_id = stage["predecessor"]
        if predecessor_id is not None:
            predecessor = by_id[predecessor_id]
            state = Path(stage["initial_state_file"])
            if predecessor.get("status") != "accepted":
                stage.update(status="blocked", reason="predecessor_not_accepted")
                manifest["status"] = "fail"
                break
            if not state.is_file():
                stage.update(status="blocked", reason="predecessor_state_missing")
                manifest["status"] = "fail"
                break
            if sha256(state) != predecessor.get("final_state_sha256"):
                stage.update(status="blocked", reason="predecessor_state_hash_mismatch")
                manifest["status"] = "fail"
                break

        for artifact_key in ("output_csv", "final_state_file"):
            Path(stage[artifact_key]).unlink(missing_ok=True)
        process = subprocess.run(
            [str(runner.resolve()), "--config", stage["config"]],
            cwd=output_dir, text=True, capture_output=True)
        console = output_dir / f"{stage['id'].replace('.', '_')}.console.log"
        console.write_text(
            process.stdout + process.stderr, encoding="utf-8", newline="\n")
        stage["returncode"] = process.returncode
        stage["console_log"] = str(console.resolve())
        if process.returncode != 0:
            stage.update(status="rejected", reason="runner_failed")
            manifest["status"] = "fail"
            break
        accepted, reason = accepted_terminal_output(stage)
        state = Path(stage["final_state_file"])
        if not accepted:
            stage.update(status="rejected", reason=reason)
            manifest["status"] = "fail"
            break
        if not state.is_file():
            stage.update(status="rejected", reason="accepted_state_missing")
            manifest["status"] = "fail"
            break
        stage.update(
            status="accepted",
            reason=reason,
            final_state_sha256=sha256(state),
        )
        by_id[stage["id"]] = stage
        write_json(output_dir / "workflow_manifest.json", manifest)
    else:
        manifest["status"] = "accepted"
    write_json(output_dir / "workflow_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--vd", type=float, action="append", dest="drain_voltages")
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--clean", action="store_true")
    args = parser.parse_args()
    if args.execute and args.runner is None:
        parser.error("--execute requires --runner")
    output_dir = args.output_dir.resolve()
    if args.clean and output_dir.exists():
        shutil.rmtree(output_dir)
    manifest = materialize(
        args.base_config.resolve(), args.contract.resolve(), output_dir,
        args.drain_voltages
        or [float(value) for value in read_json(args.contract.resolve())
            ["m3_workflow"]["drain_voltages_V"]],
    )
    if args.execute:
        manifest = execute(manifest, args.runner, output_dir)
    print(json.dumps({
        "status": manifest["status"],
        "manifest": str((output_dir / "workflow_manifest.json").resolve()),
        "stages": [{"id": item["id"], "status": item["status"]}
                   for item in manifest["stages"]],
    }))
    return 0 if manifest["status"] in {"materialized", "accepted"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
