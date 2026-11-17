#!/usr/bin/env python3
"""Run the three U-2022.12 BVMethods method regressions with sd-sim."""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

from analyze_sentaurus_bvmethods import first_upward_crossing
from run_bvmethods_ci import FAMILY, check_inputs, reference_grid
from sdevice_acceptance import csv_rows, finite, read_json, sha
from vela_sdevice.orchestrator.support import write_json

REFERENCE = FAMILY.with_name("u2022_12_20261117_methods")
METHODS = ("ABA_poisson", "ABA_coupled", "voltage2current")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def reference(method, inputs):
    check_inputs(inputs)
    receipt = read_json(REFERENCE / "reference.json")
    require(receipt["version"] == "U-2022.12" and
            receipt["mesh_sha256"] == reference_grid()["source_tdr"]["sha256"] and
            receipt["par_sha256"] == sha(inputs / "pp5_des.par"),
            "BVMethods method reference and shared inputs differ")
    base = REFERENCE / method
    contract = receipt["cases"][method]
    for name in ("run.cmd", "native.csv"):
        require(sha(base / name) == contract["files"][name], f"{method} {name} SHA-256 differs")
    return base, contract


def prepare(method, inputs, target):
    base, contract = reference(method, inputs)
    for name in ("n1_msh.grd", "n1_msh.dat", "pp5_des.par"):
        shutil.copyfile(inputs / name, target / name)
    shutil.copyfile(inputs / "pp5_des.vela.json", target / "run.vela.json")
    cmd = (base / "run.cmd").read_text(encoding="utf-8")
    deck, count = re.subn(r'(\bGrid\s*=\s*")n1_msh\.tdr(")', r'\g<1>n1_msh.grd\2', cmd)
    require(count == 1, "Expected one original grid declaration")
    (target / "run.cmd").write_text(deck, encoding="utf-8")
    return contract, csv_rows(base / "native.csv")


def qualify_aba(method, contract, native, run, plan, statuses, rows):
    require(statuses[-1].get("stop_criterion_reached") == "ionization_integral",
            "Expected a physical ionization-integral stop")
    last_stage = plan[-1]
    cfg = read_json(run / last_stage["deck"])
    solver = cfg["solver"]
    expected = "poisson_only" if method == "ABA_poisson" else "newton"
    require(solver["method"] == expected, f"Expected {expected} for ABA scan")
    require(solver["impact_ionization"]["coupling_mode"] == "postprocess_only",
            "ABA must exclude the avalanche source from the continuity solve")
    ionization = cfg["sweep"]["diagnostics"]["path_ionization_integrals"]
    require(ionization["break_rank"] == contract["rank"] and
            ionization["break_value"] == contract["threshold"],
            "ABA criterion differs from U-2022.12 CMD")
    ranks = [finite(row["path_ionization_break_rank_value"]) for row in rows]
    threshold = contract["threshold"]
    require(min(ranks) >= 0 and max(ranks[:-1]) < threshold <= ranks[-1] and
            not any(row["ionization_stop_criterion"] for row in rows[:-1]),
            "ABA did not end at the first accepted threshold crossing")
    for row in rows:
        require(min(finite(row[k]) for k in ("path_ionization_max_electron",
                                           "path_ionization_max_hole")) >= 0,
                "Negative leading path integral")
    candidate = [finite(row["bias_V"]) for row in rows[-2:]]
    expected_bracket = contract["crossing_bracket_V"]
    require(max(candidate[0], expected_bracket[0]) <= min(candidate[1], expected_bracket[1]),
            f"ABA crossing intervals do not overlap: {candidate} versus {expected_bracket}")
    result = dict(crossing_bracket_V=candidate, native_crossing_bracket_V=expected_bracket,
                  final_rank_value=ranks[-1], ranking=ionization.get("break_ordering", "path_mean"),
                  final_electron=finite(rows[-1]["path_ionization_max_electron"]),
                  final_hole=finite(rows[-1]["path_ionization_max_hole"]))
    if method == "ABA_coupled":
        reference = [r for r in native if r["stage"] == "1" and finite(r["time"]) > 0]
        q = 1.602192e-19
        native_crossing = first_upward_crossing(
            [finite(r["inner_V"]) for r in reference],
            [q * 1e-12 * finite(r["ionization_per_s"]) - finite(r["drain_A_per_um"])
             for r in reference], 0.)
        candidate = [r for r in rows if finite(r["bias_V"]) > 0]
        vela_crossing = first_upward_crossing(
            [finite(r["bias_V"]) for r in candidate],
            [finite(r["path_avalanche_current_A_per_um"]) - finite(r["current_total_A_per_um"])
             for r in candidate], 0.)
        require(native_crossing is not None and vela_crossing is not None,
                "Required avalanche-current intersection is absent")
        result.update(native_current_intersection_V=native_crossing, current_intersection_V=vela_crossing,
                      current_intersection_relative_error=abs(vela_crossing / native_crossing - 1.0))
        require(result["current_intersection_relative_error"] <= .03,
                f"Avalanche-current intersection differs: {result}")
    return result


def qualify_switch(contract, native, run, plan, statuses, rows):
    require(statuses[-2]["goal_reached"] and
            (statuses[-1]["goal_reached"] or statuses[-1].get("stop_criterion_reached") == "current"),
            "Voltage-to-current stage chain did not reach its CMD goals")
    before = csv_rows(run / plan[-2]["csv"])[-1]
    require(finite(before.get("inner_voltage_V") or before["bias_V"]) == 6.0,
            "Voltage stage did not reach 6 V")
    final = read_json(run / plan[-1]["deck"])
    sweep = final["sweep"]
    circuit = sweep["external_circuit"]
    require(circuit["mode"] == "prescribed_current" and
            sweep["initial_state_file"] == plan[-2]["output_state"],
            "Current stage must use the preceding accepted state")
    currents = [finite(r["target_current_A_per_um"]) for r in rows]
    require(len(rows) > 2 and all(r["boundary_control_mode"] == "current" for r in rows) and
            all(a < b for a, b in zip(currents, currents[1:])), "Invalid adaptive current sequence")
    require(math.isclose(currents[0], finite(before["current_total_A_per_um"]), abs_tol=1e-16),
            "Current ramp did not start at the measured switch current")
    require(math.isclose(currents[1] - currents[0], sweep["initial_step"] / circuit["resistance_ohm_um"],
                         rel_tol=1e-10, abs_tol=1e-16), "First current step differs from the configured step")
    stop = contract["final_current_A_per_um"]
    native_ramp = [r for r in native if r["stage"] == "2" and finite(r["time"]) > 0.0]
    first_fraction = (currents[1] - currents[0]) / (stop - currents[0])
    require(native_ramp and math.isclose(first_fraction, finite(native_ramp[0]["time"]),
                                         rel_tol=1e-10, abs_tol=1e-14),
            "First current step has a different normalized size from the native CMD ramp")
    require(math.isclose(currents[-1], stop, abs_tol=1e-13), "Current sweep did not reach the native goal")
    filmed = [finite(r["current_total_A_per_um"]) for r in rows]
    tolerance = finite(circuit.get("residual_tolerance_V", 1e-6)) / circuit["resistance_ohm_um"]
    max_error = max(abs(goal - measured) for goal, measured in zip(currents, filmed))
    require(max_error <= tolerance, "Current tracking exceeds the active circuit tolerance")
    for row, measured in zip(rows, filmed):
        quality = json.loads(row["validation_diagnostics"])
        require(quality["terminal_kcl_satisfied"] and quality["physical_residuals_finite"] and
                quality["actual_load_line_satisfied"] and
                abs(finite(row["current_boundary_residual_A_per_um"])) <= tolerance,
                f"Current point {measured} did not pass the independent quality checks")
    bv = first_upward_crossing([finite(r["inner_voltage_V"]) for r in rows], filmed,
                               contract["breakdown_current_A_per_um"])
    require(bv is not None, "No measured current brackets the BV criterion")
    error = abs(bv / contract["BV_inner_V"] - 1)
    require(error <= contract["BV_relative_tolerance"], f"BV differs by {error:g}")
    return dict(switch_voltage_V=6.0, initial_current_A_per_um=currents[0],
                native_initial_current_A_per_um=contract["initial_current_A_per_um"],
                first_current_step_fraction=first_fraction,
                final_inner_V=finite(rows[-1]["inner_voltage_V"]),
                native_final_inner_V=contract["final_inner_V"],
                final_current_A_per_um=filmed[-1], max_current_error_A_per_um=max_error,
                BV_inner_V=bv, native_BV_inner_V=contract["BV_inner_V"], BV_relative_error=error)


def qualify(method, contract, native, run):
    manifest = read_json(run / "run_manifest.json")
    plan = read_json(run / "stage_plan.json")["stages"]
    statuses = manifest["stages"]
    require(manifest["exit_code"] == 0 and len(statuses) == len(plan) and len(plan) >= 2 and
            all(s["runner_exit_code"] == 0 and s["accepted_points"] for s in statuses) and
            all(s["goal_reached"] for s in statuses[:-1]), "Incomplete sd-sim stage chain")
    counts = []
    for stage, status in zip(plan, statuses):
        cfg = read_json(run / stage["deck"])["solver"]
        require(cfg["state_variables"] == "potential_density" and
                cfg["primary_variable_coordinates"] == "unit", "Expected Unit population coordinates")
        rows = csv_rows(run / stage["csv"])
        count = [r for r in rows if r["converged"] in ("1", "true")]
        require(len(count) == status["accepted_points"], "CSV and manifest accepted counts differ")
        counts.append(dict(points=len(count), newton=sum(int(r["newton_iterations"]) for r in count)))
    rows = count
    require(len(rows) >= 2, "Expected at least two accepted points in the final stage")
    if method == "voltage2current":
        result = qualify_switch(contract, native, run, plan, statuses, rows)
    else:
        result = qualify_aba(method, contract, native, run, plan, statuses, rows)
    return dict(passed=True, stages=counts, **result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("method", choices=METHODS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--launcher", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    launcher = args.launcher.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="run-", dir=out) as folder:
        work = Path(folder)
        run = work / "run_vela"
        try:
            contract, native = prepare(args.method, out.parent / "inputs", work)
            env = dict(os.environ, VELA_PYTHON=sys.executable, VELA_BIN=str(launcher.parent),
                       VELA_LINEAR_SOLVER="pardiso_mkl", VELA_LINEAR_THREADS="1", MKL_NUM_THREADS="1",
                       MKL_THREADING_LAYER="SEQUENTIAL", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
            with (out / "launch.log").open("w", encoding="utf-8") as log:
                subprocess.run([str(launcher), "--physical-constants-profile", "sentaurus_u2022_12",
                                "--threads", "1", "--allow-lossy", "--exit-on-failure", "run.cmd"],
                               cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
            report = qualify(args.method, contract, native, run)
        except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
            write_json(out / "qualification.json", {"passed": False, "error": str(error)})
            raise
        finally:
            paths = [work / "run_des.log", *(run / name for name in ("run_manifest.json", "stage_plan.json",
                                                                    "execution_ir.json", "fail_report.json"))]
            paths.extend(run.glob("simulation_stage*.json"))
            paths.extend(run.glob("stage[0-9][0-9].csv"))
            paths.extend(run.glob("stage[0-9][0-9]_path_ionization_integrals*.csv"))
            paths.extend(run.glob("stage[0-9][0-9]_state.h5"))
            for path in paths:
                if path.is_file():
                    shutil.copyfile(path, out / path.name)
        write_json(out / "qualification.json", report)
    print(f"{args.method} passed")


if __name__ == "__main__":
    main()
