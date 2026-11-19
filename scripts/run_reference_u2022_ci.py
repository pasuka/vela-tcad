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
from compare_genius_bjt_sentaurus_vela import compare_model
from run_schottky_reference_workflow import audit_kcl, audit_stage_b, compare_curves, read_csv
from transportmodels_idvd_fixture import stage_idvd_inputs
from vela_sdevice.plt import read_plt_document


REFERENCE_ROOT = Path(__file__).resolve().parents[1] / "reference_tcad"
CASES = {
    "schottky_iv": "schottky_charon_sentaurus2018",
    "bjt_m1": "genius_bjt_sentaurus2022",
    "gate_contacts": "gate_reference",
    "transportmodels_idvd": "transportmodels_sentaurus2022",
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


def transportmodels_idvd_gates(launcher, fixture, output, env, report):
    work = output / "inputs"
    stage_idvd_inputs(fixture, work)
    run, plan = run_sdevice(launcher, work, "pp12_des", env, report, final_only=True)
    equilibrium, sweep_stage = plan
    biases = {"gate": 1.0, "drain": 0.0, "source": 0.0, "substrate": 0.0}
    if (equilibrium["kind"] != "equilibrium"
            or equilibrium["electrode_voltages_start_V"] != biases
            or sweep_stage["electrode_voltages_goal_V"] != {**biases, "drain": 2.0}
            or sweep_stage["sweep"]["contact"] != "drain"
            or sweep_stage["input_state"] != equilibrium["output_state"]):
        raise ValueError("Unexpected TransportModels equilibrium or sweep")
    for stage in plan:
        solver = read_json(run / stage["deck"])["solver"]
        expected = {"state_variables": "potential_density", "primary_variable_coordinates": "unit",
                    "max_iter": 20, "temperature_K": 300.0}
        if any(solver.get(key) != value for key, value in expected.items()):
            raise ValueError(f"Unexpected TransportModels solver settings in {stage['deck']}")
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
    reference = read_csv(fixture / "sentaurus_reference.csv")
    targets = [point / 10 for point in range(21)]
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
    contract = read_json(fixture.parent / "vela/contracts/transportmodels_dd_dg_sentaurus2022_v1.json")
    limits = contract["acceptance"]
    max_error = max(row["relative_error"] for row in comparison[1:])
    end_error = comparison[-1]["relative_error"]
    report["comparison"] = {
        "points": len(comparison),
        "max_relative_error": max_error, "endpoint_relative_error": end_error,
        "zero_bias_drain_A_per_um": comparison[0]["vela_A_per_um"],
        "max_abs_kcl_residual_A_per_um": max(abs(row["kcl_residual_A_per_um"]) for row in comparison),
    }
    print(f"DD Id–Vd max error: {max_error:.3%}, at 2 V: {end_error:.3%}", flush=True)
    return (max_error <= limits["idvd_max_absolute_relative_error"]
            and end_error <= limits["idvd_2V_max_absolute_relative_error"])


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
        elif args.case == "transportmodels_idvd":
            passed = transportmodels_idvd_gates(launcher, fixture, output, env, report)
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
