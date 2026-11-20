#!/usr/bin/env python3
"""Run the device reference numerical regressions through sd-sim."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

from check_sdevice_gate_contact_mapping import flatband_from_samples, read_vtu_potential
from close_singledevice_validation import branch_metrics, load_comparison
from compare_genius_bjt_sentaurus_vela import compare_model
from run_schottky_reference_workflow import audit_kcl, audit_stage_b, compare_curves, read_csv
from transportmodels_idvd_fixture import stage_transportmodels_inputs
from vela_sdevice.plt import read_plt_document


REFERENCE_ROOT = Path(__file__).resolve().parents[1] / "reference_tcad"
CASES = {
    "schottky_iv": "schottky_charon_sentaurus2018",
    "bjt_m1": "genius_bjt_sentaurus2022",
    "gate_contacts": "gate_reference",
    "transportmodels_idvd": "transportmodels_sentaurus2022",
    "transportmodels_dg_idvd": "transportmodels_sentaurus2022",
    "transportmodels_idvg": "transportmodels_sentaurus2022",
    "singledevice_idvg_linear": "singledevice_sentaurus2018",
}


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def verify_inputs(fixture):
    manifest = read_json(fixture / "manifest.json")
    for name, expected in manifest["files_sha256_lf"].items():
        data = (fixture / name).read_bytes().replace(b"\r\n", b"\n")
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f"Reference input checksum mismatch: {fixture / name}")


def stage_inputs(fixture, output, stem, grid):
    output.mkdir()
    for name in (stem + ".cmd", stem + ".vela.json", "models.par"):
        shutil.copyfile(fixture / name, output / name)
    for path in grid:
        shutil.copyfile(path, output / path.name)


def run_sdevice(launcher, output, stem, env, report, *, final_only=False):
    command = [str(launcher), "--physical-constants-profile", "sentaurus_u2022_12",
               "--threads", "1", "--exit-on-failure", "-q"]
    if final_only:
        command.append("--plot-final-only")
    log_path = output / "sd-sim.log"
    started = time.monotonic()
    with log_path.open("wb") as log:
        result = subprocess.run([*command, stem + ".cmd"], cwd=output, env=env,
                                stdout=log, stderr=subprocess.STDOUT)
    run = output / (stem + "_vela")
    status_path = run / "run_manifest.json"
    status = read_json(status_path) if status_path.is_file() else {}
    report["runs"].append({"deck": stem, "exit_code": result.returncode,
                           "seconds": time.monotonic() - started,
                           "stages": status.get("stages", [])})
    if result.returncode != 0:
        tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-12:])
        raise RuntimeError(f"{stem}: sd-sim failed. See {log_path}\n{tail}")
    plan = read_json(run / "stage_plan.json")["stages"]
    if (status["exit_code"] != 0 or len(status["stages"]) != len(plan)
            or not all(stage["goal_reached"] and stage["runner_exit_code"] == 0
                       and stage["accepted_points"] > 0 for stage in status["stages"])
            or not all((run / stage["output_state"]).is_file() for stage in plan)):
        raise RuntimeError(f"{stem}: incomplete sd-sim stage sequence")
    return run, plan


def gate_cases(launcher, fixture, output, env, report):
    reference = read_json(fixture.parent / "source/wp14_gate_inputs.json")
    grid = [fixture / f"gate_msh.{ext}" for ext in ("grd", "dat")]
    results = {}
    report["flatband"] = results
    for case in ("barrier_zero", "barrier_positive", "barrier_negative",
                 "workfunction", "workfunction_zero"):
        stem = case + "_sdevice"
        work = output / case
        stage_inputs(fixture, work, stem, grid)
        run, plan = run_sdevice(launcher, work, stem, env, report)
        nodes = read_json(run / "mesh.json")["nodes"]
        silicon_top = reference["geometry_um"]["silicon_thickness"]
        interface = [row["id"] for row in nodes if math.isclose(row["y"], silicon_top, abs_tol=1e-12)]
        substrate = [row["id"] for row in nodes if math.isclose(row["y"], 0.0, abs_tol=1e-12)]
        if not interface or not substrate:
            raise ValueError("Gate mesh has no silicon interface or substrate probe")
        [series] = plan[-1]["vtu_write"]["series"]
        spec = reference["variants"][case]
        if len(series["t"]) != spec["samples"]:
            raise ValueError(f"{case}: incorrect number of requested field snapshots")
        lo, hi = spec["scan_biases_V"]
        samples = {}
        for index, t in enumerate(series["t"], start=series["global_index_start"]):
            path = work / series["pattern"].replace("NNNN", f"{index:04d}")
            bias, potential = read_vtu_potential(path)
            if (not math.isclose(bias, lo + (hi - lo) * t, abs_tol=1e-9)
                    or len(potential) != len(nodes) or not all(map(math.isfinite, potential))):
                raise ValueError(f"{path.name}: incorrect bias or invalid potential")
            samples[bias] = (math.fsum(potential[i] for i in interface) / len(interface)
                             - math.fsum(potential[i] for i in substrate) / len(substrate))
        if len(samples) != spec["samples"] or not samples[min(samples)] < 0 < samples[max(samples)]:
            raise ValueError(f"{case}: scan does not bracket flat band")
        measured = flatband_from_samples(samples)
        expected = spec["expected_flatband_V"]
        results[case] = {"measured_V": measured, "expected_V": expected,
                         "error_V": measured - expected, "samples": len(samples),
                         "surface_minus_bulk_V": samples}
    zero_barrier = results["barrier_zero"]["measured_V"]
    expected_shifts = (
        ("barrier_shift_V", "barrier_positive", 0.55),
        ("negative_barrier_shift_V", "barrier_negative", -0.30),
        ("zero_workfunction_shift_V", "workfunction_zero", 0.0),
    )
    for field, case, _ in expected_shifts:
        report[field] = results[case]["measured_V"] - zero_barrier
    return (all(abs(report[field] - target) <= 1e-3 for field, _, target in expected_shifts)
            and all(abs(row["error_V"]) <= 1e-3 for row in results.values()))


def schottky_gates(fixture, run, stage, report):
    curve_path = run / stage["csv"]
    rows = [row for row in read_csv(curve_path) if row["converged"] == "1"]
    points = [
        {key: float(row[key]) for key in ("bias_V", "current_total_A_per_um")}
        for row in rows
    ]
    curve = compare_curves(points, fixture / "sentaurus_reference.csv")
    branch = audit_stage_b(rows)
    components = audit_kcl(rows)
    terminals = {}
    for row in read_csv(curve_path.with_name(curve_path.stem + "_terminal_balance.csv")):
        if row["converged"] == "1":
            bias = float(row["bias_V"])
            terminals.setdefault(bias, {})[row["contact"]] = float(row["current_total_A_per_um"])
    port_kcl = []
    for target in (0.0, 0.4, 1.0):
        bias = min(terminals, key=lambda value: abs(value - target))
        currents = terminals[bias]
        if currents.keys() != {"anode", "cathode"} or not all(map(math.isfinite, currents.values())):
            raise ValueError(f"Missing or invalid terminal currents at {bias} V")
        absolute = abs(math.fsum(currents.values()))
        relative = absolute / max(max(map(abs, currents.values())), 1e-300)
        port_kcl.append({
            "target_V": target, "bias_V": bias, "currents_A_per_um": currents,
            "absolute_A_per_um": absolute, "relative": relative,
            "passed": absolute <= 1e-14 or relative <= 0.01,
        })
    report.update(curve=curve, monotonic=branch, components=components, port_kcl=port_kcl)
    return all(gate["status"] == "pass" for gate in (curve, branch, components)) and all(
        point["passed"] for point in port_kcl
    )


def bjt_gates(fixture, run, stage, report):
    curve_path = run / stage["csv"]
    contract = read_json(fixture.parent / "contracts" / "comparison_thresholds.json")
    parity = contract["wp3_wp5_vela_comparison"]["numerical_parity"]
    points, summary = compare_model(
        "M1", fixture / "sentaurus_reference.csv", curve_path.with_name(curve_path.stem + "_terminal_balance.csv"),
        tuple(parity["active_region_VCE_range_V"]), parity["models"]["M1"],
    )
    with (run.parent / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)
    report["comparison"] = summary
    return summary["operational_pass"] and summary["numerical_parity_gate"]["pass"]


def transportmodels_idvd_gates(launcher, fixture, output, env, report, *, density_gradient=False):
    work = output / "inputs"
    stem = "pp13_des" if density_gradient else "pp12_des"
    stage_transportmodels_inputs(fixture, work, stem)
    run, plan = run_sdevice(launcher, work, stem, env, report, final_only=True)
    contract = read_json(fixture.parent / "vela/contracts/transportmodels_dd_dg_sentaurus2022_v1.json")
    equilibrium, sweep_stage = plan
    biases = {"gate": 1.0, "drain": 0.0, "source": 0.0, "substrate": 0.0}
    if (equilibrium["kind"] != "equilibrium"
            or equilibrium["electrode_voltages_start_V"] != biases
            or sweep_stage["electrode_voltages_goal_V"] != {**biases, "drain": 2.0}
            or sweep_stage["sweep"]["contact"] != "drain"
            or sweep_stage["input_state"] != equilibrium["output_state"]):
        raise ValueError("Unexpected TransportModels equilibrium or sweep")
    coordinates = ({"state_variables": "potential_quasi_fermi", "quasi_fermi_reference": "contact_basin"}
                   if density_gradient else {"state_variables": "potential_density", "primary_variable_coordinates": "unit"})
    max_iter = read_json(fixture / f"{stem}.vela.json")["solver"].get("max_iter", 20)
    for stage in plan:
        solver = read_json(run / stage["deck"])["solver"]
        expected = {**coordinates, "max_iter": max_iter, "temperature_K": 300.0}
        if any(solver.get(key) != value for key, value in expected.items()):
            raise ValueError(f"Unexpected TransportModels solver settings in {stage['deck']}")
        if density_gradient:
            quantum = solver.get("electron_quantum_potential", {})
            required = ("enabled", "include_insulators", "global_discretization", "coefficient_mass_ratio",
                        "outer_absolute_tolerance_V", "outer_max_iterations", "max_iterations")
            if any(quantum.get(key) != contract["dg_quantum_contract"][key] for key in required):
                raise ValueError(f"Unexpected density-gradient settings in {stage['deck']}")
    poisson = read_json(run / "stage00_poisson_initialization.json")
    if poisson["max_iter"] != 100 or not poisson["converged"]:
        raise ValueError("TransportModels Poisson initialization did not match the CMD budget")
    sweep = read_json(run / sweep_stage["deck"])["sweep"]
    controls = {
        "initial_step": 0.002, "min_step": 0.00002, "max_step": 2.0,
        "growth_factor": 2.0, "shrink_factor": 0.5,
        "step_growth_mode": "newton_iterations", "stop_below_minimum_step": True,
        "preserve_step_after_output_clip": True,
    }
    if (any(sweep.get(key) != value for key, value in controls.items())
            or sweep["continuation"]["predictor"]["mode"] != "linear"):
        raise ValueError("Unexpected TransportModels Quasistationary controls")
    report["sweep_controls"] = {**controls, "predictor": "linear"}

    contacts, samples = read_plt_document(
        (work / sweep_stage["current_file"]).read_text(encoding="utf-8"))
    reference_name = "sentaurus_dg_idvd_reference.csv" if density_gradient else "sentaurus_reference.csv"
    reference = read_csv(fixture / reference_name)
    targets = contract["bias_contract"]["idvd"]["drain_bias_V"]
    if (set(contacts) != biases.keys() or len(samples) != len(targets)
            or [float(row["bias_V"]) for row in reference] != targets):
        raise ValueError("Expected 21 measured TransportModels CurrentPlot points and all four contacts")

    comparison = []
    for bias, sample, oracle in zip(targets, samples, reference):
        expected_voltages = {**biases, "drain": bias}
        if any(not math.isclose(sample[name]["outer_voltage_V"], voltage, rel_tol=0, abs_tol=1e-12)
               for name, voltage in expected_voltages.items()):
            raise ValueError(f"Incorrect electrode voltage at target Vd={bias:g} V")
        currents = {name: sample[name]["total_current_A_per_um"] for name in contacts}
        native = float(oracle["drain_A_per_um"])
        measured = currents["drain"]
        comparison.append({
            "bias_V": bias, "vela_A_per_um": measured, "sentaurus_A_per_um": native,
            "relative_error": abs(measured - native) / abs(native) if bias > 0 else None,
            "kcl_residual_A_per_um": math.fsum(currents.values()),
        })
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)
    limits = contract["acceptance"]
    max_error = max(row["relative_error"] for row in comparison[1:])
    end_error = comparison[-1]["relative_error"]
    report["comparison"] = {
        "points": len(comparison),
        "max_relative_error": max_error, "endpoint_relative_error": end_error,
        "zero_bias_drain_A_per_um": comparison[0]["vela_A_per_um"],
        "max_abs_kcl_residual_A_per_um": max(abs(row["kcl_residual_A_per_um"]) for row in comparison),
    }
    if density_gradient:
        qualify_density_gradient(run, plan, report)
    branch = "DG" if density_gradient else "DD"
    print(f"{branch} Id-Vd max error: {max_error:.3%}, at 2 V: {end_error:.3%}", flush=True)
    return (max_error <= limits["idvd_max_absolute_relative_error"]
            and end_error <= limits["idvd_2V_max_absolute_relative_error"])


def qualify_density_gradient(run, plan, report):
    from state_archive import inspect

    node_count = len(read_json(run / "mesh.json")["nodes"])
    potentials = []
    report["electron_quantum_potential"] = []
    for stage in plan:
        points = [row for row in read_csv(run / stage["csv"]) if row["converged"] == "1"]
        if not points or any(not row["newton_convergence_reason"].endswith("+electron_density_gradient")
                            for row in points):
            raise ValueError("Expected self-consistent density-gradient convergence at every accepted point")
        fields, metadata = inspect(run / stage["output_state"])
        values = fields.get("electron_quantum_potential_V")
        contact = stage.get("sweep", {}).get("contact", "drain")
        bias = (stage["electrode_voltages_goal_V"][contact] if "sweep" in stage
                else float(points[-1]["bias_V"]))
        if (values is None or len(values) != node_count or not all(math.isfinite(v) for v in values)
                or not any(values) or not math.isclose(metadata["bias_V"], bias, abs_tol=1e-12)):
            raise ValueError(f"Missing or inconsistent quantum state at {contact}={bias:g} V")
        potentials.append(values)
        report["electron_quantum_potential"].append({
            f"{contact}_bias_V": bias, "nodes": node_count, "accepted_dg_points": len(points),
            "minimum_V": float(min(values)), "maximum_V": float(max(values)),
        })
    change = max(abs(a - b) for a, b in zip(potentials[0], potentials[-1], strict=True))
    if change <= 0:
        raise ValueError("The electron quantum potential did not respond to the applied bias")
    report["max_quantum_potential_change_V"] = float(change)


def transportmodels_idvg_samples(launcher, fixture, work, env, report, bias_contract, *, deep_off=False):
    stem = "pp6_des"
    stage_transportmodels_inputs(fixture, work, stem)
    targets = bias_contract["gate_bias_V"]
    gate_span = targets[-1] - targets[0]
    if deep_off:
        targets = targets[:3]
        # Preserve the absolute initial/minimum step on the shortened ramp.
        factor = gate_span / (targets[-1] - targets[0])
        deck = (work / f"{stem}.cmd").read_text(encoding="utf-8")
        replacements = {
            "InitialStep=1e-3 MinStep=1e-5 MaxStep=1":
                f"InitialStep={1e-3 * factor:g} MinStep={1e-5 * factor:g} MaxStep=1",
            'Goal{ Name="gate" Voltage= 2.2 }': f'Goal{{ Name="gate" Voltage= {targets[-1]:g} }}',
            "Intervals=20": f"Intervals={len(targets)-1}",
        }
        for old, new in replacements.items():
            if deck.count(old) != 1:
                raise ValueError("The staged Id-Vg CMD no longer matches the deep-off sweep controls")
            deck = deck.replace(old, new)
        stem = "pp6_deep_off"
        (work / f"{stem}.cmd").write_text(deck, encoding="utf-8")
        shutil.copyfile(fixture / f"{stem}.vela.json", work / f"{stem}.vela.json")

    run, plan = run_sdevice(launcher, work, stem, env, report, final_only=True)
    equilibrium, drain_stage, gate_stage = plan
    initial = {"gate": targets[0], "drain": 0.0, "source": 0.0, "substrate": 0.0}
    drain_goal = {**initial, "drain": bias_contract["drain_bias_V"]}
    if equilibrium["kind"] != "equilibrium" or equilibrium["electrode_voltages_start_V"] != initial:
        raise ValueError("Expected cold equilibrium at Vg=-1 V and Vd=0 V")

    coordinates = ({"state_variables": "potential_quasi_fermi", "quasi_fermi_reference": "contact_basin"}
                   if deep_off else {"state_variables": "potential_density", "primary_variable_coordinates": "unit"})
    for index, stage in enumerate(plan):
        solver = read_json(run / stage["deck"])["solver"]
        expected = {**coordinates, "max_iter": 20, "temperature_K": 300.0}
        if any(solver.get(key) != value for key, value in expected.items()):
            raise ValueError(f"Unexpected TransportModels solver settings in {stage['deck']}")
        if index and stage["input_state"] != plan[index-1]["output_state"]:
            raise ValueError("TransportModels stages must consume their preceding state")
    poisson = read_json(run / "stage00_poisson_initialization.json")
    if poisson["max_iter"] != 100 or not poisson["converged"]:
        raise ValueError("TransportModels Poisson initialization did not match the CMD budget")

    controls = report.setdefault("sweep_controls", {})[stem] = {}
    for stage, contact, start, goal, initial_step in (
            (drain_stage, "drain", initial, drain_goal, 0.01),
            (gate_stage, "gate", drain_goal, {**drain_goal, "gate": targets[-1]}, 0.001)):
        if (stage["kind"] != "sweep" or stage["sweep"]["contact"] != contact
                or stage["electrode_voltages_start_V"] != start or stage["electrode_voltages_goal_V"] != goal):
            raise ValueError(f"Unexpected TransportModels {contact} sweep")
        span = gate_span if contact == "gate" else goal[contact] - start[contact]
        numeric = {"initial_step": initial_step * span, "min_step": 1e-5 * span,
                   "max_step": goal[contact] - start[contact],
                   "growth_factor": 2.0, "shrink_factor": 0.5}
        sweep = read_json(run / stage["deck"])["sweep"]
        discrete = {"step_growth_mode": "newton_iterations", "stop_below_minimum_step": True,
                    "preserve_step_after_output_clip": True}
        predictor = sweep["continuation"]["predictor"]["mode"]
        if (any(not math.isclose(sweep[key], value, rel_tol=1e-12) for key, value in numeric.items())
                or any(sweep.get(key) != value for key, value in discrete.items()) or predictor != "linear"):
            raise ValueError(f"Unexpected Quasistationary controls in {stage['deck']}")
        controls[contact] = {**numeric, **discrete, "predictor": predictor}

    contacts, samples = read_plt_document((work / gate_stage["current_file"]).read_text(encoding="utf-8"))
    if set(contacts) != initial.keys() or len(samples) != len(targets):
        raise ValueError(f"Expected {len(targets)} measured CurrentPlot points and all four contacts in {stem}")
    for bias, sample in zip(targets, samples):
        if any(not math.isclose(sample[name]["outer_voltage_V"], voltage, rel_tol=0, abs_tol=1e-12)
               for name, voltage in {**drain_goal, "gate": bias}.items()):
            raise ValueError(f"Incorrect electrode voltage at target Vg={bias:g} V in {stem}")
    return contacts, samples


def transportmodels_idvg_gates(launcher, fixture, output, env, report):
    contract = read_json(fixture.parent / "vela/contracts/transportmodels_dd_dg_sentaurus2022_v1.json")
    bias_contract = contract["bias_contract"]["idvg"]
    targets = bias_contract["gate_bias_V"]
    reference = read_csv(fixture / "sentaurus_idvg_reference.csv")
    if [float(row["bias_V"]) for row in reference] != targets:
        raise ValueError("Expected the 21 native TransportModels Id-Vg biases")
    contacts, population = transportmodels_idvg_samples(
        launcher, fixture, output / "inputs", env, report, bias_contract)
    _, deep = transportmodels_idvg_samples(
        launcher, fixture, output / "deep_off_inputs", env, report, bias_contract, deep_off=True)
    samples = deep + population[len(deep):]
    limits = contract["acceptance"]
    comparison = []
    for index, (bias, sample, oracle) in enumerate(zip(targets, samples, reference)):
        currents = {name: sample[name]["total_current_A_per_um"] for name in contacts}
        if not all(math.isfinite(value) for value in (*currents.values(), population[index]["drain"]["total_current_A_per_um"])):
            raise ValueError(f"Nonfinite current at target Vg={bias:g} V")
        native = float(oracle["drain_A_per_um"])
        measured = currents["drain"]
        kcl = math.fsum(currents.values())
        ratio = measured / abs(kcl) if kcl else None
        log_error = abs(math.log10(measured / native)) if measured > 0 else None
        relative_error = abs(measured / native - 1)
        regime = "deep_off" if index < 3 else "transition" if index < 8 else "on"
        if regime == "on":
            passed = relative_error <= limits["idvg_on_max_absolute_relative_error"]
        else:
            passed = log_error is not None and log_error <= limits["idvg_transition_max_absolute_log_error_dex"]
            if regime == "deep_off":
                passed = passed and (ratio is None or ratio >= limits["deep_off_min_id_to_kcl_ratio"])
        comparison.append({
            "bias_V": bias, "regime": regime, "coordinate_mode": "contact_basin" if index < len(deep) else "population",
            "vela_A_per_um": measured, "sentaurus_A_per_um": native,
            "population_drain_A_per_um": population[index]["drain"]["total_current_A_per_um"],
            "relative_error": relative_error, "log_error_dex": log_error, "kcl_residual_A_per_um": kcl,
            "drain_to_abs_kcl_ratio": ratio, "pass": passed,
        })
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)
    report["comparison"] = {
        "points": len(comparison), "failed_biases_V": [row["bias_V"] for row in comparison if not row["pass"]],
        "max_on_relative_error": max(row["relative_error"] for row in comparison if row["regime"] == "on"),
        "max_abs_kcl_residual_A_per_um": max(abs(row["kcl_residual_A_per_um"]) for row in comparison),
    }
    return not report["comparison"]["failed_biases_V"]


def singledevice_idvg_gates(launcher, fixture, output, env, report):
    stem = "singledevice_idvg"
    inputs = output / "inputs"
    stage_inputs(fixture, inputs, stem, [fixture / f"device.{ext}" for ext in ("grd", "dat")])
    shutil.copyfile(fixture.parent.parent / "transportmodels_sentaurus2022/ci_u2022_12/polysilicon.par",
                    inputs / "polysilicon.par")
    run, plan = run_sdevice(launcher, inputs, stem, env, report, final_only=True)

    manifest = read_json(fixture / "manifest.json")
    mesh = read_json(run / "mesh.json")
    inventory = {"nodes": len(mesh["nodes"]), "triangles": len(mesh["triangles"])}
    if inventory != manifest["mesh"]["inventory"]:
        raise ValueError(f"Unexpected SingleDevice mesh inventory: {inventory}")
    report["mesh"] = manifest["mesh"]

    equilibrium, drain, gate = plan
    initial = dict(source=0.0, drain=0.0, gate=-0.5, substrate=0.0)
    prebiased = {**initial, "drain": 0.1}
    endpoint = {**prebiased, "gate": 2.2}
    if (equilibrium["kind"] != "equilibrium" or equilibrium["input_state"] is not None
            or equilibrium["electrode_voltages_goal_V"] != initial
            or equilibrium["output_state"] != "singledevice_init_des.h5"):
        raise ValueError("Expected a cold SingleDevice equilibrium saved at Vgs=-0.5 V and Vds=0 V")
    for stage, preceding, contact, start, goal, initial_step, min_step, max_step, growth in (
        (drain, equilibrium, "drain", initial, prebiased, 0.001, 1e-6, 0.01, 2.0),
        (gate, drain, "gate", prebiased, endpoint, 0.0027, 2.7e-5, 0.27, 1.3),
    ):
        if (stage["kind"] != "sweep" or stage["sweep"]["contact"] != contact
                or stage["input_state"] != preceding["output_state"]
                or stage["electrode_voltages_start_V"] != start
                or stage["electrode_voltages_goal_V"] != goal):
            raise ValueError(f"Unexpected SingleDevice {contact} sweep or restart state")
        sweep = read_json(run / stage["deck"])["sweep"]
        numeric = dict(initial_step=initial_step, min_step=min_step, max_step=max_step,
                       growth_factor=growth, shrink_factor=0.5)
        controls = dict(step_growth_mode="newton_iterations", stop_below_minimum_step=True,
                        preserve_step_after_output_clip=True)
        if (any(not math.isclose(sweep[name], value, rel_tol=1e-12) for name, value in numeric.items())
                or any(sweep.get(name) != value for name, value in controls.items())
                or sweep["continuation"]["predictor"]["mode"] != "linear"):
            raise ValueError(f"Unexpected SingleDevice Quasistationary controls in {stage['deck']}")
    poisson = read_json(run / "stage00_poisson_initialization.json")
    if poisson["max_iter"] != 100 or not poisson["converged"]:
        raise ValueError("SingleDevice Poisson initialization did not match the CMD budget")
    for stage in plan:
        solver = read_json(run / stage["deck"])["solver"]
        if (solver["state_variables"] != "potential_quasi_fermi"
                or solver["temperature_K"] != 300
                or not solver["electron_quantum_potential"]["enabled"]):
            raise ValueError(f"Unexpected SingleDevice solver settings in {stage['deck']}")
    qualify_density_gradient(run, plan, report)

    reference = read_json(fixture.parent / "singledevice_sentaurus2018_reference.json")
    linear = next(branch for branch in reference["simulations"] if branch["name"] == "idvg_lin")
    oracle = read_csv(fixture.parent / linear["reference_csv"])
    contacts, samples = read_plt_document((inputs / gate["current_file"]).read_text(encoding="utf-8"))
    if set(contacts) != initial.keys() or len(samples) != len(oracle) or len(oracle) != linear["point_count"]:
        raise ValueError("Expected 21 SingleDevice CurrentPlot biases and all four contacts")

    comparison = []
    for sample, expected in zip(samples, oracle, strict=True):
        bias = float(expected["bias_V"])
        if any(not math.isclose(sample[name]["outer_voltage_V"], voltage, rel_tol=0, abs_tol=1e-12)
               for name, voltage in {**prebiased, "gate": bias}.items()):
            raise ValueError(f"Incorrect electrode voltage at target Vgs={bias:g} V")
        currents = [sample[name]["total_current_A_per_um"] for name in contacts]
        current = sample["drain"]["total_current_A_per_um"]
        if not all(map(math.isfinite, currents)) or current <= 0.0:
            raise ValueError(f"Invalid terminal currents at Vgs={bias:g} V")
        residual = abs(math.fsum(currents))
        comparison.append({
            "gate_voltage_V": bias, "sentaurus_current_A_per_um": float(expected["current_total"]),
            "vela_current_A_per_um": current, "abs_kcl_residual_A_per_um": residual,
            "relative_kcl_residual": residual / current,
        })
    path = output / "comparison.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(comparison[0]))
        writer.writeheader()
        writer.writerows(comparison)

    curve = branch_metrics(load_comparison(path), reference["low_current_acceptance"])
    limits = reference["acceptance"]
    checkpoints = (comparison[0], comparison[6], comparison[-1])
    checks = {
        "pointwise": curve["status"] == "pass",
        "trend": curve["trend_match"],
        "median_log_error": curve["median_abs_log_error_dex"] <= limits["median_log10_current_error_max_dex"],
        "p95_log_error": curve["p95_abs_log_error_dex"] <= limits["p95_log10_current_error_max_dex"],
        "ion": curve["ion"]["relative_error"] <= limits["ion_relative_error_max"],
        "threshold": abs(curve["threshold"]["delta_V"]) <= limits["threshold_voltage_error_max_V"],
        "subthreshold_swing": curve["subthreshold_swing"]["relative_error"]
        <= limits["subthreshold_swing_relative_error_max"],
        "terminal_balance": all(row["abs_kcl_residual_A_per_um"] <= 1e-14
                                or row["relative_kcl_residual"] <= 0.01 for row in checkpoints),
    }
    report.update(curve=curve, acceptance=checks, terminal_balance_checkpoints=checkpoints,
                  reference_version=reference["sentaurus_version"])
    print(f"Linear DG Id-Vg at Vds=0.1 V: {len(samples)} output points, "
          f"Ion error {curve['ion']['relative_error']:.3%}", flush=True)
    return all(checks.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--launcher", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--linear-solver", required=True, choices=("pardiso_mkl", "umfpack", "sparselu"))
    args = parser.parse_args()
    launcher = args.launcher.resolve(strict=True)
    fixture = REFERENCE_ROOT / CASES[args.case] / "ci_u2022_12"
    args.output.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix=f"{args.case}-", dir=args.output.resolve()))
    env = dict(
        os.environ, VELA_PYTHON=sys.executable, VELA_BIN=str(launcher.parent),
        VELA_LINEAR_SOLVER=args.linear_solver, VELA_LINEAR_THREADS="1",
        MKL_THREADING_LAYER="SEQUENTIAL", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
    )
    env.pop("MKL_CBWR", None)
    report = {"case": args.case, "linear_solver": args.linear_solver, "status": "incomplete", "runs": []}
    print(f"Reports: {output}", flush=True)
    try:
        verify_inputs(fixture)
        if args.case == "gate_contacts":
            passed = gate_cases(launcher, fixture, output, env, report)
        elif args.case in ("transportmodels_idvd", "transportmodels_dg_idvd"):
            passed = transportmodels_idvd_gates(
                launcher, fixture, output, env, report, density_gradient=args.case == "transportmodels_dg_idvd")
        elif args.case == "transportmodels_idvg":
            passed = transportmodels_idvg_gates(launcher, fixture, output, env, report)
        elif args.case == "singledevice_idvg_linear":
            passed = singledevice_idvg_gates(launcher, fixture, output, env, report)
        else:
            schottky = args.case == "schottky_iv"
            stem = "schottky_n_des" if schottky else "bjt_m1_des"
            mesh = "schottky_n_msh" if schottky else "bjt_msh"
            grid_dir = REFERENCE_ROOT / "sdevice_check_only" if schottky else fixture
            work = output / "inputs"
            stage_inputs(fixture, work, stem, [grid_dir / f"{mesh}.{ext}" for ext in ("grd", "dat")])
            run, plan = run_sdevice(launcher, work, stem, env, report, final_only=True)
            expected = [("anode", 1.0)] if schottky else [("base", 0.7), ("collector", 3.0)]
            actual = [(stage["sweep"]["contact"], stage["sweep"]["stop_V"]) for stage in plan[1:]]
            if plan[0]["kind"] != "equilibrium" or actual != expected:
                raise ValueError(f"Unexpected {args.case} stage sequence: {actual}")
            gates = schottky_gates if schottky else bjt_gates
            passed = gates(fixture, run, plan[-1], report)
        report["status"] = "passed" if passed else "gate_failed"
    finally:
        write_json(output / "result.json", report)
    print(f"{args.case}: {report['status']}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
