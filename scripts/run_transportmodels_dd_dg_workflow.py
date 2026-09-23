#!/usr/bin/env python3
"""Materialize and run the Sentaurus TransportModels DD/DG MOS workflow.

The official project contains two independent SDevice solves per transport
branch.  Id-Vg starts at Vg=-1 V, solves equilibrium, ramps Vd to 1.1 V,
relaxes once more at the final bias, then sweeps the gate.  Id-Vd starts from
a separate equilibrium at Vg=1 V and sweeps the drain.  This driver preserves
that dependency graph and solves the exact 21-point reference bias lattices.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any


INPUT_PATH_KEYS = ("mesh_file", "node_doping_file", "materials_file")
BRANCHES = ("dd", "dg")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reference_biases(path: Path) -> list[float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path}: reference curve is empty")
    key = next((name for name in ("bias_V", "gate_voltage_V", "bias")
                if name in rows[0]), "")
    if not key:
        raise ValueError(f"{path}: cannot identify bias column")
    values = [float(row[key]) for row in rows]
    if values != sorted(values) or len(set(values)) != len(values):
        raise ValueError(f"{path}: biases must be unique and increasing")
    return values


def comparison_prefix_biases(path: Path) -> list[float]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"bias_V", "current_total_A_per_um"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: comparison prefix lacks {sorted(required)}")
    values = [float(row["bias_V"]) for row in rows]
    if values != sorted(values) or len(set(values)) != len(values):
        raise ValueError(f"{path}: prefix biases must be unique and increasing")
    return values


def absolutize_inputs(config: dict[str, Any], base_dir: Path) -> None:
    for key in INPUT_PATH_KEYS:
        value = config.get(key)
        if value and not Path(value).is_absolute():
            config[key] = str((base_dir / value).resolve())


def set_contact(config: dict[str, Any], name: str, bias: float) -> None:
    for contact in config.get("contacts", []):
        if contact.get("name") == name:
            contact["bias"] = bias
            return
    raise ValueError(f"base deck has no contact {name!r}")


def solver_signature(config: dict[str, Any]) -> dict[str, Any]:
    signature = json.loads(json.dumps(config.get("solver", {})))
    signature.pop("electron_quantum_potential", None)
    return signature


def validate_controlled_delta(generated_dir: Path) -> None:
    decks = {
        name: load_json(generated_dir / "vela" / f"simulation_{name}.json")
        for name in ("dd_idvg", "dd_idvd", "dg_idvg", "dg_idvd")
    }
    for branch in BRANCHES:
        if solver_signature(decks[f"{branch}_idvg"]) != solver_signature(
                decks[f"{branch}_idvd"]):
            raise ValueError(f"{branch} Id-Vg and Id-Vd solver physics differ")
    if solver_signature(decks["dd_idvg"]) != solver_signature(decks["dg_idvg"]):
        raise ValueError("DD/DG solver physics differ beyond electron quantum potential")
    if "electron_quantum_potential" in decks["dd_idvg"].get("solver", {}):
        raise ValueError("DD deck unexpectedly enables electron quantum potential")
    if not decks["dg_idvg"].get("solver", {}).get(
            "electron_quantum_potential", {}).get("enabled", False):
        raise ValueError("DG deck does not enable electron quantum potential")
    for key in INPUT_PATH_KEYS:
        values = {deck.get(key) for deck in decks.values()}
        if len(values) != 1:
            raise ValueError(f"DD/DG decks do not share {key}")


def diagnostics(run_dir: Path, stem: str) -> dict[str, Any]:
    return {
        "transport": {"enabled": True},
        "terminal_balance": {
            "enabled": True,
            "contacts": ["source", "drain", "gate", "substrate"],
            "csv_file": str((run_dir / f"{stem}_terminal_balance.csv").resolve()),
        },
    }


def stage_config(base: dict[str, Any], base_dir: Path, run_dir: Path,
                 stem: str, *, gate_bias: float, drain_bias: float,
                 sweep_contact: str, bias_points: list[float],
                 initial_state: Path | None, cold_start: bool,
                 quantum_outer_max_iterations: int = 40,
                 quantum_outer_acceleration: str = "aitken",
                 quantum_outer_relaxation: float = 1.0) -> dict[str, Any]:
    if not bias_points:
        raise ValueError(f"{stem}: bias_points must not be empty")
    config = json.loads(json.dumps(base))
    absolutize_inputs(config, base_dir)
    set_contact(config, "source", 0.0)
    set_contact(config, "substrate", 0.0)
    set_contact(config, "gate", gate_bias)
    set_contact(config, "drain", drain_bias)
    output = (run_dir / f"{stem}.csv").resolve()
    final_state = (run_dir / f"{stem}_final_state.h5").resolve()
    config["output_csv"] = str(output)
    config["log_file"] = str((run_dir / f"{stem}.log").resolve())
    step = (bias_points[1] - bias_points[0]) if len(bias_points) > 1 else 0.01
    sweep = config.setdefault("sweep", {})
    sweep.clear()
    sweep.update({
        "mode": "iv",
        "contact": sweep_contact,
        "current_contact": "drain",
        "start": bias_points[0],
        "stop": bias_points[-1],
        "step": step,
        "bias_points": bias_points,
        "write_vtk": False,
        "write_state_file": str(final_state),
        "write_state_every_point_prefix": str(
            (run_dir / f"{stem}_state").resolve()),
        "diagnostics": diagnostics(run_dir, stem),
    })
    if cold_start:
        sweep["initialization"] = {"mode": "poisson_block"}
    elif initial_state is not None:
        sweep["initial_state_file"] = str(initial_state.resolve())
    else:
        raise ValueError(f"{stem}: restarted stage requires initial_state")
    solver = config.setdefault("solver", {})
    solver["max_iter"] = max(int(solver.get("max_iter", 0)), 100)
    solver["quasi_fermi_update_limit_V"] = 0.1 if cold_start else 0.025
    mobility = solver.get("mobility", {})
    if mobility.get("high_field_driving_force") == "quasi_fermi_gradient":
        # Sentaurus GradQuasiFermi is a two-dimensional cell-vector
        # magnitude.  The edge-projection control is orientation dependent
        # and loses the DD branch near Vd=0.84 V on this imported MOS mesh.
        mobility["high_field_gradient_discretization"] = \
            "transport_cell_vector"
    quantum = solver.get("electron_quantum_potential", {})
    if quantum.get("enabled", False):
        # The TransportModels MOS drain ramp has a steadily contracting
        # quantum outer fixed-point iteration, but needs 24 iterations near
        # Vd=0.35 V.  At 0.50 V, vector Aitken acceleration converges in 29
        # iterations where the unaccelerated residual is still 1.1e-4 V at
        # iteration 20.  Keep explicit headroom without changing tolerances.
        quantum["outer_max_iterations"] = max(
            int(quantum.get("outer_max_iterations", 0)),
            quantum_outer_max_iterations)
        quantum["outer_acceleration"] = quantum_outer_acceleration
        quantum["outer_relaxation"] = quantum_outer_relaxation
        quantum["outer_relaxation_min"] = 0.1
        quantum["outer_relaxation_max"] = 1.5
    return config


def drain_ramp_biases() -> list[float]:
    # The saved equilibrium already represents Vd=0.  Starting at the next
    # bias avoids needlessly re-solving the identical state after restart.
    return [round(index * 0.05, 14) for index in range(1, 23)]


def materialize(generated_dir: Path, run_dir: Path,
                branches: list[str],
                idvg_ramp_states: dict[str, Path] | None = None,
                idvg_curve_restarts: dict[str, dict[str, Any]] | None = None,
                idvd_curve_restarts: dict[str, dict[str, Any]] | None = None,
                idvd_bridge_biases: dict[str, list[float]] | None = None,
                quantum_outer_max_iterations: int = 40,
                quantum_outer_acceleration: str = "aitken",
                quantum_outer_relaxation: float = 1.0) -> dict[str, Any]:
    validate_controlled_delta(generated_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    idvg_ramp_states = idvg_ramp_states or {}
    idvg_curve_restarts = idvg_curve_restarts or {}
    idvd_curve_restarts = idvd_curve_restarts or {}
    idvd_bridge_biases = idvd_bridge_biases or {}
    stages: list[dict[str, Any]] = []
    physics_hashes: dict[str, str] = {}
    stage_index = 0
    for branch in branches:
        idvg_base_path = (generated_dir / "vela" /
                          f"simulation_{branch}_idvg.json").resolve()
        idvd_base_path = (generated_dir / "vela" /
                          f"simulation_{branch}_idvd.json").resolve()
        idvg_reference = (generated_dir / "reference_curves" /
                          f"transportmodels_sentaurus2022_{branch}_idvg_reference.csv").resolve()
        idvd_reference = (generated_dir / "reference_curves" /
                          f"transportmodels_sentaurus2022_{branch}_idvd_reference.csv").resolve()
        idvg_biases = reference_biases(idvg_reference)
        idvd_biases = reference_biases(idvd_reference)
        idvg_base = load_json(idvg_base_path)
        idvd_base = load_json(idvd_base_path)
        physics_hashes[branch] = hashlib.sha256(json.dumps(
            idvg_base.get("solver", {}), sort_keys=True).encode()).hexdigest()

        eq_vg_name = f"{branch}_idvg_equilibrium"
        eq_vg_state = (run_dir / f"{eq_vg_name}_final_state.h5").resolve()
        ramp_name = f"{branch}_idvg_drain_ramp"
        ramp_state = (run_dir / f"{ramp_name}_final_state.h5").resolve()
        relax_name = f"{branch}_idvg_final_bias_relax"
        relax_state = (run_dir / f"{relax_name}_final_state.h5").resolve()
        eq_vd_name = f"{branch}_idvd_equilibrium"
        eq_vd_state = (run_dir / f"{eq_vd_name}_final_state.h5").resolve()

        external_ramp_state = idvg_ramp_states.get(branch)
        curve_restart = idvg_curve_restarts.get(branch)
        idvd_curve_restart = idvd_curve_restarts.get(branch)
        idvd_bridges = idvd_bridge_biases.get(branch, [])
        if idvd_curve_restart is not None and (
                external_ramp_state is not None or curve_restart is not None):
            raise ValueError(
                f"{branch}: Id-Vd curve restart cannot be combined with "
                "an Id-Vg restart")
        if external_ramp_state is not None and curve_restart is not None:
            raise ValueError(
                f"{branch}: Id-Vg ramp state and curve restart are mutually exclusive")
        if curve_restart is not None:
            restart_state = Path(curve_restart["state"]).resolve()
            restart_prefix = Path(curve_restart["prefix"]).resolve()
            restart_bias = float(curve_restart["bias_V"])
            if not restart_state.is_file():
                raise ValueError(
                    f"{branch}: external Id-Vg curve state does not exist: "
                    f"{restart_state}")
            if not restart_prefix.is_file():
                raise ValueError(
                    f"{branch}: external Id-Vg curve prefix does not exist: "
                    f"{restart_prefix}")
            matching_bias = next(
                (value for value in idvg_biases
                 if abs(value - restart_bias) <= 1.0e-10), None)
            if matching_bias is None or restart_bias >= idvg_biases[-1]:
                raise ValueError(
                    f"{branch}: curve restart bias {restart_bias} is not a "
                    "nonterminal reference bias")
            remaining_idvg_biases = [
                value for value in idvg_biases
                if value > restart_bias + 1.0e-10]
            idvg_prefix = []
            relax_initial = None
            relax_dependencies = []
        elif external_ramp_state is not None:
            external_ramp_state = external_ramp_state.resolve()
            if not external_ramp_state.is_file():
                raise ValueError(
                    f"{branch}: external Id-Vg ramp state does not exist: "
                    f"{external_ramp_state}")
            idvg_prefix = []
            relax_initial = external_ramp_state
            relax_dependencies: list[str] = []
        else:
            idvg_prefix = [
                (eq_vg_name, idvg_base, idvg_base_path.parent, -1.0, 0.0,
                 "gate", [idvg_biases[0]], None, True, [], None),
                (ramp_name, idvg_base, idvg_base_path.parent, -1.0, 0.0,
                 "drain", drain_ramp_biases(), eq_vg_state, False,
                 [eq_vg_name], None),
            ]
            relax_initial = ramp_state
            relax_dependencies = [ramp_name]

        if curve_restart is not None:
            idvg_specs = [
                (f"{branch}_idvg_curve", idvg_base, idvg_base_path.parent,
                 restart_bias, 1.1, "gate", remaining_idvg_biases,
                 restart_state, False, [],
                 {"external_prefix": str(restart_prefix),
                  "reference_biases": idvg_biases}),
            ]
        else:
            idvg_specs = idvg_prefix + [
                (relax_name, idvg_base, idvg_base_path.parent, -1.0, 1.1,
                 "gate", [idvg_biases[0]], relax_initial, False,
                 relax_dependencies, None),
                (f"{branch}_idvg_curve", idvg_base, idvg_base_path.parent,
                 -1.0, 1.1, "gate", idvg_biases[1:], relax_state, False,
                 [relax_name],
                 {"stage": relax_name, "bias_V": idvg_biases[0],
                  "reference_biases": idvg_biases}),
            ]

        if idvd_curve_restart is not None:
            restart_state = Path(idvd_curve_restart["state"]).resolve()
            restart_prefix = Path(idvd_curve_restart["prefix"]).resolve()
            restart_bias = float(idvd_curve_restart["bias_V"])
            if not restart_state.is_file():
                raise ValueError(
                    f"{branch}: external Id-Vd curve state does not exist: "
                    f"{restart_state}")
            if not restart_prefix.is_file():
                raise ValueError(
                    f"{branch}: external Id-Vd curve prefix does not exist: "
                    f"{restart_prefix}")
            matching_bias = next(
                (value for value in idvd_biases
                 if abs(value - restart_bias) <= 1.0e-10), None)
            if matching_bias is None or restart_bias >= idvd_biases[-1]:
                raise ValueError(
                    f"{branch}: Id-Vd curve restart bias {restart_bias} is "
                    "not a nonterminal reference bias")
            prefix_biases = comparison_prefix_biases(restart_prefix)
            expected_prefix = [
                value for value in idvd_biases
                if value <= restart_bias + 1.0e-10]
            if len(prefix_biases) != len(expected_prefix) or any(
                    abs(left - right) > 1.0e-10
                    for left, right in zip(prefix_biases, expected_prefix)):
                raise ValueError(
                    f"{branch}: Id-Vd prefix lattice {prefix_biases} != "
                    f"completed reference prefix {expected_prefix}")
            if idvd_bridges != sorted(idvd_bridges) or len(
                    set(idvd_bridges)) != len(idvd_bridges):
                raise ValueError(
                    f"{branch}: Id-Vd bridge biases must be unique and increasing")
            for bridge in idvd_bridges:
                if not restart_bias + 1.0e-10 < bridge < idvd_biases[-1] - 1.0e-10:
                    raise ValueError(
                        f"{branch}: Id-Vd bridge bias {bridge} must lie between "
                        f"restart {restart_bias} and terminal {idvd_biases[-1]}")
                if any(abs(bridge - value) <= 1.0e-10 for value in idvd_biases):
                    raise ValueError(
                        f"{branch}: Id-Vd bridge bias {bridge} is already a "
                        "reference bias")
            remaining_idvd_biases = [
                value for value in idvd_biases
                if value > restart_bias + 1.0e-10]
            stage_idvd_biases = idvd_bridges + remaining_idvd_biases
            # An Id-Vd restart is a terminal workflow slice: all Id-Vg work
            # and the Id-Vd equilibrium/prefix are supplied as completed
            # external evidence, so only the missing exact biases are emitted.
            idvg_specs = []
            idvd_specs = [
                (f"{branch}_idvd_curve", idvd_base, idvd_base_path.parent,
                 1.0, restart_bias, "drain", stage_idvd_biases,
                 restart_state, False, [],
                 {"external_prefix": str(restart_prefix),
                  "reference_biases": idvd_biases}),
            ]
        else:
            if idvd_bridges:
                raise ValueError(
                    f"{branch}: Id-Vd bridge biases require a curve restart")
            idvd_specs = [
                (eq_vd_name, idvd_base, idvd_base_path.parent, 1.0, 0.0,
                 "gate", [1.0], None, True, [], None),
                (f"{branch}_idvd_curve", idvd_base, idvd_base_path.parent,
                 1.0, 0.0, "drain", idvd_biases[1:], eq_vd_state, False,
                 [eq_vd_name],
                 {"stage": eq_vd_name, "bias_V": idvd_biases[0],
                  "reference_biases": idvd_biases}),
            ]
        specs = idvg_specs + idvd_specs
        for (name, base, base_dir, gate, drain, contact, biases, initial,
             cold, dependencies, comparison_seed) in specs:
            config = stage_config(
                base, base_dir, run_dir, name, gate_bias=gate,
                drain_bias=drain, sweep_contact=contact,
                bias_points=biases, initial_state=initial, cold_start=cold,
                quantum_outer_max_iterations=quantum_outer_max_iterations,
                quantum_outer_acceleration=quantum_outer_acceleration,
                quantum_outer_relaxation=quantum_outer_relaxation)
            config_path = (run_dir / f"{stage_index:02d}_{name}.json").resolve()
            write_json(config_path, config)
            stages.append({
                "name": name,
                "branch": branch,
                "config": str(config_path),
                "config_sha256": sha256(config_path),
                "depends_on": dependencies,
                "bias_points": biases,
                "reference": str(idvg_reference if "idvg_curve" in name else
                                 idvd_reference if "idvd_curve" in name else ""),
                "output_csv": config["output_csv"],
                "initial_state_file": config["sweep"].get("initial_state_file"),
                "final_state_file": config["sweep"]["write_state_file"],
            })
            if comparison_seed is not None:
                stages[-1]["comparison_seed"] = comparison_seed
            if name == relax_name and external_ramp_state is not None:
                stages[-1]["external_initial_state"] = {
                    "role": "completed_idvg_drain_ramp",
                    "path": str(external_ramp_state),
                    "sha256": sha256(external_ramp_state),
                    "terminal_bias_V": 1.1,
                }
            if (name == f"{branch}_idvg_curve" and
                    curve_restart is not None):
                stages[-1]["external_initial_state"] = {
                    "role": "completed_idvg_curve_prefix",
                    "path": str(restart_state),
                    "sha256": sha256(restart_state),
                    "terminal_bias_V": restart_bias,
                    "prefix_csv": str(restart_prefix),
                    "prefix_sha256": sha256(restart_prefix),
                }
            if (name == f"{branch}_idvd_curve" and
                    idvd_curve_restart is not None):
                stages[-1]["external_initial_state"] = {
                    "role": "completed_idvd_curve_prefix",
                    "path": str(restart_state),
                    "sha256": sha256(restart_state),
                    "terminal_bias_V": restart_bias,
                    "prefix_csv": str(restart_prefix),
                    "prefix_sha256": sha256(restart_prefix),
                }
            stage_index += 1

    manifest = {
        "schema": "vela.transportmodels.dd_dg_workflow.v1",
        "status": "materialized",
        "generated_dir": str(generated_dir.resolve()),
        "branches": branches,
        "controlled_delta": {
            "shared_mesh_doping_materials": True,
            "dd_to_dg_only_electron_quantum_potential": True,
            "solver_sha256": physics_hashes,
        },
        "stages": stages,
    }
    write_json(run_dir / "workflow_manifest.json", manifest)
    return manifest


def prepare_comparison_candidate(stage: dict[str, Any],
                                 stages_by_name: dict[str, dict[str, Any]],
                                 run_dir: Path) -> Path:
    seed = stage["comparison_seed"]
    if "external_prefix" in seed:
        with Path(seed["external_prefix"]).open(
                newline="", encoding="utf-8") as handle:
            prefix_rows = list(csv.DictReader(handle))
        required = {"bias_V", "current_total_A_per_um"}
        if not prefix_rows or not required.issubset(prefix_rows[0]):
            raise ValueError(
                f"{stage['name']}: external prefix lacks {sorted(required)}")
        rows = [{
            "bias_V": float(row["bias_V"]),
            "current_total_A_per_um": float(row["current_total_A_per_um"]),
        } for row in prefix_rows]
    else:
        seed_stage = stages_by_name[seed["stage"]]
        with Path(seed_stage["output_csv"]).open(
                newline="", encoding="utf-8") as handle:
            seed_rows = [row for row in csv.DictReader(handle)
                         if row.get("converged") == "1"]
        if not seed_rows:
            raise ValueError(f"{stage['name']}: missing converged seed row")
        rows = [{
            "bias_V": float(seed["bias_V"]),
            "current_total_A_per_um": float(
                seed_rows[-1]["current_total_A_per_um"]),
        }]
    with Path(stage["output_csv"]).open(
            newline="", encoding="utf-8") as handle:
        curve_rows = [row for row in csv.DictReader(handle)
                      if row.get("converged") == "1"]
    if not curve_rows:
        raise ValueError(f"{stage['name']}: missing converged curve rows")
    expected = [float(value) for value in seed["reference_biases"]]
    prefix_biases = [row["bias_V"] for row in rows]
    rows.extend({
        "bias_V": float(row["bias_V"]),
        "current_total_A_per_um": float(row["current_total_A_per_um"]),
    } for row in curve_rows if any(
        abs(float(row["bias_V"]) - value) <= 1.0e-10
        for value in expected) and not any(
            abs(float(row["bias_V"]) - value) <= 1.0e-10
            for value in prefix_biases))
    actual = [row["bias_V"] for row in rows]
    if len(actual) != len(expected) or any(
            abs(left - right) > 1.0e-10 for left, right in zip(actual, expected)):
        raise ValueError(
            f"{stage['name']}: comparison lattice {actual} != reference {expected}")
    output = run_dir / f"{stage['name']}_comparison_candidate.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output


def compare_curve(stage: dict[str, Any], stages_by_name: dict[str, dict[str, Any]],
                  run_dir: Path) -> dict[str, Any]:
    report_dir = run_dir / "reports"
    report_dir.mkdir(parents=True, exist_ok=True)
    output_json = report_dir / f"{stage['name']}_comparison.json"
    output_md = report_dir / f"{stage['name']}_comparison.md"
    candidate = prepare_comparison_candidate(stage, stages_by_name, run_dir)
    command = [
        sys.executable,
        str(Path(__file__).with_name("compare_reference_curves.py")),
        "--reference", stage["reference"],
        "--candidate", str(candidate),
        "--output-json", str(output_json),
        "--output-md", str(output_md),
        "--kind", "iv",
        "--min-points", str(len(stage["comparison_seed"]["reference_biases"])),
        "--reference-column", "current_total",
        "--candidate-column", "current_total_A_per_um",
        "--interpolation", "log_current",
        "--require-trend-match",
    ]
    process = subprocess.run(
        command, text=True, capture_output=True, cwd=run_dir)
    return {
        "status": "pass" if process.returncode == 0 else "fail",
        "returncode": process.returncode,
        "json": str(output_json.resolve()),
        "markdown": str(output_md.resolve()),
        "candidate": str(candidate.resolve()),
        "stdout": process.stdout,
        "stderr": process.stderr,
    }


def execute(manifest: dict[str, Any], runner: Path, run_dir: Path,
            previous: dict[str, Any] | None = None) -> dict[str, Any]:
    prior = {stage["name"]: stage for stage in (previous or {}).get("stages", [])}
    for stage in manifest["stages"]:
        old = prior.get(stage["name"], {})
        state = Path(stage["final_state_file"])
        can_reuse = bool(
            old.get("status") == "pass" and state.is_file()
            and old.get("config_sha256") == stage["config_sha256"])
        if can_reuse:
            stage["status"] = "pass"
            stage["execution"] = "reused_verified"
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
        if state.is_file():
            stage["final_state_sha256"] = sha256(state)
        write_json(run_dir / "workflow_manifest.json", manifest)
        if process.returncode != 0:
            manifest["status"] = "fail"
            break
    else:
        manifest["status"] = "pass"
    comparisons: dict[str, Any] = {}
    stages_by_name = {stage["name"]: stage for stage in manifest["stages"]}
    for stage in manifest["stages"]:
        if stage.get("status") == "pass" and stage.get("reference"):
            comparisons[stage["name"]] = compare_curve(
                stage, stages_by_name, run_dir)
    manifest["comparisons"] = comparisons
    if any(item["status"] == "fail" for item in comparisons.values()):
        manifest["comparison_status"] = "fail"
    elif comparisons:
        manifest["comparison_status"] = "pass"
    write_json(run_dir / "workflow_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--branches", default="dd,dg")
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--clean", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--idvg-ramp-state", action="append", default=[], metavar="BRANCH=PATH",
        help=("restart a branch after its completed Id-Vg drain ramp; may be "
              "specified more than once"),
    )
    parser.add_argument(
        "--idvg-curve-state", action="append", default=[],
        metavar="BRANCH=BIAS=PATH",
        help=("restart a branch after a completed exact Id-Vg point; requires "
              "the matching --idvg-curve-prefix"),
    )
    parser.add_argument(
        "--idvg-curve-prefix", action="append", default=[],
        metavar="BRANCH=PATH",
        help="completed comparison-candidate prefix for an Id-Vg curve restart",
    )
    parser.add_argument(
        "--idvd-curve-state", action="append", default=[],
        metavar="BRANCH=BIAS=PATH",
        help=("restart a branch after a completed exact Id-Vd point; requires "
              "the matching --idvd-curve-prefix"),
    )
    parser.add_argument(
        "--idvd-curve-prefix", action="append", default=[],
        metavar="BRANCH=PATH",
        help="completed comparison-candidate prefix for an Id-Vd curve restart",
    )
    parser.add_argument(
        "--idvd-bridge-bias", action="append", default=[],
        metavar="BRANCH=BIAS",
        help=("non-reference Id-Vd continuation bias to solve and checkpoint "
              "before the remaining exact comparison biases; repeatable"),
    )
    parser.add_argument(
        "--quantum-outer-max-iterations", type=int, default=40,
        help="minimum DG outer fixed-point iteration budget (default: 40)",
    )
    parser.add_argument(
        "--quantum-outer-acceleration", choices=("aitken", "none"),
        default="aitken",
        help="DG outer fixed-point acceleration policy (default: aitken)",
    )
    parser.add_argument(
        "--quantum-outer-relaxation", type=float, default=1.0,
        help="base DG outer relaxation in [0.1, 1.5] (default: 1.0)",
    )
    args = parser.parse_args()
    if args.execute and args.runner is None:
        parser.error("--execute requires --runner")
    if args.clean and args.resume:
        parser.error("--clean and --resume are mutually exclusive")
    branches = [item.strip().lower() for item in args.branches.split(",")
                if item.strip()]
    unknown = sorted(set(branches) - set(BRANCHES))
    if unknown or not branches:
        parser.error("--branches must be a non-empty comma-separated subset of dd,dg")
    if args.quantum_outer_max_iterations < 1:
        parser.error("--quantum-outer-max-iterations must be positive")
    if not 0.1 <= args.quantum_outer_relaxation <= 1.5:
        parser.error("--quantum-outer-relaxation must lie in [0.1, 1.5]")
    idvg_ramp_states: dict[str, Path] = {}
    for value in args.idvg_ramp_state:
        branch, separator, raw_path = value.partition("=")
        branch = branch.strip().lower()
        if not separator or not branch or not raw_path.strip():
            parser.error("--idvg-ramp-state must use BRANCH=PATH")
        if branch not in branches:
            parser.error(
                f"--idvg-ramp-state branch {branch!r} is not selected by --branches")
        if branch in idvg_ramp_states:
            parser.error(f"duplicate --idvg-ramp-state for branch {branch!r}")
        idvg_ramp_states[branch] = Path(raw_path.strip())
    curve_states: dict[str, tuple[float, Path]] = {}
    for value in args.idvg_curve_state:
        parts = value.split("=", 2)
        if len(parts) != 3:
            parser.error("--idvg-curve-state must use BRANCH=BIAS=PATH")
        branch, raw_bias, raw_path = (item.strip() for item in parts)
        branch = branch.lower()
        if branch not in branches:
            parser.error(
                f"--idvg-curve-state branch {branch!r} is not selected")
        if branch in curve_states:
            parser.error(f"duplicate --idvg-curve-state for branch {branch!r}")
        try:
            bias = float(raw_bias)
        except ValueError:
            parser.error(f"invalid Id-Vg curve restart bias {raw_bias!r}")
        curve_states[branch] = (bias, Path(raw_path))
    curve_prefixes: dict[str, Path] = {}
    for value in args.idvg_curve_prefix:
        branch, separator, raw_path = value.partition("=")
        branch = branch.strip().lower()
        if not separator or not branch or not raw_path.strip():
            parser.error("--idvg-curve-prefix must use BRANCH=PATH")
        if branch in curve_prefixes:
            parser.error(f"duplicate --idvg-curve-prefix for branch {branch!r}")
        curve_prefixes[branch] = Path(raw_path.strip())
    if set(curve_states) != set(curve_prefixes):
        parser.error(
            "--idvg-curve-state and --idvg-curve-prefix branches must match")
    idvg_curve_restarts = {
        branch: {"bias_V": bias, "state": path,
                 "prefix": curve_prefixes[branch]}
        for branch, (bias, path) in curve_states.items()
    }
    idvd_curve_states: dict[str, tuple[float, Path]] = {}
    for value in args.idvd_curve_state:
        parts = value.split("=", 2)
        if len(parts) != 3:
            parser.error("--idvd-curve-state must use BRANCH=BIAS=PATH")
        branch, raw_bias, raw_path = (item.strip() for item in parts)
        branch = branch.lower()
        if branch not in branches:
            parser.error(
                f"--idvd-curve-state branch {branch!r} is not selected")
        if branch in idvd_curve_states:
            parser.error(f"duplicate --idvd-curve-state for branch {branch!r}")
        try:
            bias = float(raw_bias)
        except ValueError:
            parser.error(f"invalid Id-Vd curve restart bias {raw_bias!r}")
        idvd_curve_states[branch] = (bias, Path(raw_path))
    idvd_curve_prefixes: dict[str, Path] = {}
    for value in args.idvd_curve_prefix:
        branch, separator, raw_path = value.partition("=")
        branch = branch.strip().lower()
        if not separator or not branch or not raw_path.strip():
            parser.error("--idvd-curve-prefix must use BRANCH=PATH")
        if branch in idvd_curve_prefixes:
            parser.error(f"duplicate --idvd-curve-prefix for branch {branch!r}")
        idvd_curve_prefixes[branch] = Path(raw_path.strip())
    if set(idvd_curve_states) != set(idvd_curve_prefixes):
        parser.error(
            "--idvd-curve-state and --idvd-curve-prefix branches must match")
    idvd_curve_restarts = {
        branch: {"bias_V": bias, "state": path,
                 "prefix": idvd_curve_prefixes[branch]}
        for branch, (bias, path) in idvd_curve_states.items()
    }
    idvd_bridge_biases: dict[str, list[float]] = {}
    for value in args.idvd_bridge_bias:
        branch, separator, raw_bias = value.partition("=")
        branch = branch.strip().lower()
        if not separator or branch not in branches or not raw_bias.strip():
            parser.error("--idvd-bridge-bias must use selected BRANCH=BIAS")
        try:
            bias = float(raw_bias)
        except ValueError:
            parser.error(f"invalid Id-Vd bridge bias {raw_bias!r}")
        idvd_bridge_biases.setdefault(branch, []).append(bias)
    output_dir = args.output_dir.resolve()
    manifest_path = output_dir / "workflow_manifest.json"
    previous = load_json(manifest_path) if args.resume and manifest_path.is_file() else None
    if args.resume and previous is None:
        parser.error("--resume requires an existing workflow_manifest.json")
    if args.clean and output_dir.exists():
        shutil.rmtree(output_dir)
    manifest = materialize(
        args.generated_dir.resolve(), output_dir, branches, idvg_ramp_states,
        idvg_curve_restarts, idvd_curve_restarts, idvd_bridge_biases,
        args.quantum_outer_max_iterations, args.quantum_outer_acceleration,
        args.quantum_outer_relaxation)
    if args.execute:
        manifest = execute(
            manifest, args.runner.resolve(), output_dir, previous if args.resume else None)
    print(json.dumps({
        "status": manifest["status"],
        "comparison_status": manifest.get("comparison_status"),
        "stages": [
            {"name": stage["name"], "status": stage.get("status", "materialized")}
            for stage in manifest["stages"]
        ],
    }))
    return 0 if manifest["status"] in {"materialized", "pass"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
