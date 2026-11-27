"""Run an opt-in SLOT LDMOS CMD package and qualify the sd-sim sweep."""

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from vela_sdevice.cmd_parser import parse_cmd_file


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def finite_number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite sweep result: {value}")
    return result


def check_precision(precision, manifest):
    audit = manifest["arithmetic"]
    if audit["requested"] != precision:
        raise ValueError("sd-sim precision request differs from the CMD Math section")
    if precision != "extended80":
        return
    capabilities = audit["reported_capabilities"]
    if (audit["capability_check"] != "passed"
            or audit["frontend_input_check"]["status"] != "passed"
            or (capabilities["arithmetic_precision"], capabilities["radix"],
                capabilities["significand_bits"], capabilities["max_exponent"])
            != ("extended80", 2, 64, 16384)):
        raise ValueError("The runtime did not qualify the requested float80 solve")


def qualify(case, cmd):
    work = cmd.with_name(cmd.stem + "_vela")
    manifest = read_json(work / "run_manifest.json")
    plan = read_json(work / "stage_plan.json")
    source = parse_cmd_file(cmd)
    precision = source["math"].get("arithmetic_precision", "double")
    if manifest["check_only"] or manifest["exit_code"] != 0:
        raise ValueError("sd-sim did not complete a simulation")
    if plan["cmd"]["sha256"] != source["source"]["sha256"]:
        raise ValueError("The executed stage plan does not match the input CMD")
    check_precision(precision, manifest)
    stages = {item["index"]: item for item in manifest["stages"]}
    if len(stages) != len(manifest["stages"]) or list(stages) != [
            item["index"] for item in plan["stages"]]:
        raise ValueError("The stage plan and executed stages differ")

    contacts = ("drain", "gate") if case == "idsat" else ("gate" if case == "idvg" else "drain",)
    source_goals = {stage["index"]: stage["goal"]
                    for stage in source["stages"]
                    if stage.get("goal", {}).get("contact") in contacts
                    and "voltage_V" in stage["goal"]}
    sweeps = []
    idvd_gate = None
    for stage in plan["stages"]:
        run = stages[stage["index"]]
        if run["runner_exit_code"] != 0 or not run["accepted_points"]:
            raise ValueError(f"Stage {stage['index']} did not produce an accepted solution")
        deck = read_json(work / stage["deck"])
        solver = deck["solver"]
        if solver.get("state_variables") != "potential_density":
            raise ValueError(f"Stage {stage['index']} did not use potential and density unknowns")
        if deck.get("format_version") == 2:
            # V2 uses UnitScaling implicitly; the legacy 'scaling' key is forbidden.
            scaled = "scaling" not in deck
        else:
            scaled = deck.get("scaling", {}).get("mode") == "unit_scaling"
        if not scaled:
            raise ValueError(f"Stage {stage['index']} did not use unit scaling")
        if (solver.get("arithmetic_precision", "double") != precision
                or run["requested_arithmetic_precision"] != precision
                or run["reported_arithmetic_precision"] != precision):
            raise ValueError(f"Stage {stage['index']} did not use the CMD arithmetic precision")
        if precision == "extended80":
            state = run.get("accepted_extended80_state", {})
            if (solver.get("require_end_to_end_precision") is not True
                    or state.get("schema") != "vela.extended80-stage-state.v1"
                    or state.get("role") != "accepted"
                    or state.get("archive") != stage["output_state"]
                    or state.get("accepted_points") != run["accepted_points"]
                    or not (work / stage["output_state"]).is_file()):
                raise ValueError(f"Stage {stage['index']} did not produce a verified float80 accepted state")
        sweep = stage.get("sweep", {})
        if sweep.get("contact") not in contacts or sweep.get("mode") != "voltage":
            if not run["goal_reached"]:
                raise ValueError(f"Preparation stage {stage['index']} did not reach its goal")
            if case == "idvd" and sweep.get("contact") == "gate" and sweep.get("mode") == "voltage":
                vstart = stage["electrode_voltages_start_V"]
                vend = run["reached_voltages_V"]
                if vstart["drain"] == vend["drain"] == 0 and vend["gate"] != vstart["gate"]:
                    idvd_gate = vend["gate"]
            continue
        if case == "idvd" and not sweeps:
            start = stage["electrode_voltages_start_V"]
            if idvd_gate is None or start["gate"] != idvd_gate or start["drain"] != 0.0:
                raise ValueError("IdVd must start from the accepted zero-drain gate sweep")
        series_resistor = deck["sweep"].get("external_circuit", {}).get("mode") == "series_resistor"
        if (case == "bvds") != series_resistor:
            raise ValueError(f"Stage {stage['index']} has the wrong drain circuit for {case}")
        source_goal = source_goals.get(stage["ir_stage_index"])
        if (source_goal is None or sweep["contact"] != source_goal["contact"]
                or finite_number(sweep["stop_V"]) != source_goal["voltage_V"]):
            raise ValueError(f"Stage {stage['index']} voltage goal differs from the CMD")
        sweeps.append((stage, run, sweep))

    if not sweeps:
        raise ValueError(f"No {case} voltage sweep in {cmd.name}")
    if [stage["ir_stage_index"] for stage, _, _ in sweeps] != list(source_goals):
        raise ValueError(f"The executed {case} voltage sweeps differ from the CMD")
    if case == "idsat":
        if tuple(sweep["contact"] for _, _, sweep in sweeps) != contacts:
            raise ValueError("IDSAT requires a drain sweep followed by a gate sweep")
        (drain_stage, drain_run, _), (gate_stage, gate_run, _) = sweeps
        load = plan["stages"][drain_stage["index"] + 1]
        drain_bias = drain_run["reached_voltages_V"]["drain"]
        gate_start = gate_stage["electrode_voltages_start_V"]
        if (drain_stage["electrode_voltages_start_V"]["drain"] != 0
                or drain_stage["electrode_voltages_start_V"]["gate"] != 0
                or drain_run["reached_voltages_V"]["gate"] != 0
                or load["kind"] != "load" or load["index"] + 1 != gate_stage["index"]
                or load["input_state"] != drain_stage["output_state"]
                or gate_stage["input_state"] != load["output_state"]
                or gate_start["gate"] != 0 or gate_start["drain"] != drain_bias
                or gate_run["reached_voltages_V"]["drain"] != drain_bias):
            raise ValueError("IDSAT must save and reload the zero-gate drain endpoint before sweeping the gate")
    for stage, run, sweep in sweeps:
        contact = sweep["contact"]
        with (work / stage["csv"]).open(newline="", encoding="utf-8") as stream:
            rows = [row for row in csv.DictReader(stream) if row["converged"] == "1"]
        if len(rows) != run["accepted_points"]:
            raise ValueError(f"Stage {stage['index']} accepted-point count differs from the manifest")
        start = finite_number(sweep["start_V"])
        stop = finite_number(sweep["stop_V"])
        if start == stop:
            raise ValueError(f"Stage {stage['index']} has a zero voltage span")
        for row in rows:
            finite_number(row["inner_voltage_V"])
            finite_number(row["current_total_A_per_um"])
        end = finite_number(rows[-1]["outer_voltage_V"])
        if (end - start) * (stop - start) <= 0:
            raise ValueError(f"Stage {stage['index']} made no voltage progress")
        reached_goal = run["goal_reached"] and math.isclose(end, stop, rel_tol=1e-10, abs_tol=1e-9)
        if not reached_goal:
            current = abs(finite_number(rows[-1]["current_total_A_per_um"]))
            reached_break = (case == "bvds" and run.get("stop_criterion_reached") == "current"
                             and any(criterion["contact"] == contact
                                     and current >= finite_number(criterion["abs_val_A_per_um"])
                                     for criterion in stage.get("break", {}).get("current", [])))
            if not reached_break:
                raise ValueError(f"Stage {stage['index']} did not reach its voltage or breakdown target")
        print(f"{case} stage {stage['index']} ({precision}): "
              f"{len(rows)} accepted points, final {contact}={end:.12g} V")

    if not any(Path(name).suffix.lower() == ".plt" and (cmd.parent / name).is_file()
               for name in manifest["outputs"]):
        raise ValueError("sd-sim did not publish a current curve")
    if not (cmd.parent / plan["outputs"]["log"]).is_file():
        raise ValueError("sd-sim did not publish its runtime log")
    if not (cmd.parent / plan["outputs"]["final_vtu"]).is_file():
        raise ValueError("sd-sim did not publish its final physical fields")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("bvds", "idvg", "idvd", "idsat"))
    parser.add_argument("--cmd", type=Path, required=True)
    parser.add_argument("--launcher", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--default-par", type=Path)
    args = parser.parse_args()
    cmd = args.cmd.resolve(strict=True)
    launcher = args.launcher.resolve(strict=True)
    output = args.output.resolve()
    if output.is_relative_to(cmd.parent):
        parser.error("The output directory must be outside the input package")
    output.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix=args.case + "-", dir=output))
    print(f"SLOT LDMOS {args.case} artifacts: {work}", flush=True)
    shutil.copytree(cmd.parent, work, dirs_exist_ok=True, ignore=shutil.ignore_patterns(cmd.stem + "_vela"))
    env = dict(os.environ, VELA_PYTHON=sys.executable, VELA_BIN=str(launcher.parent))
    command = [str(launcher), "--exit-on-failure", "--keep-intermediate", "--threads", "1"]
    if args.default_par:
        defaults = args.default_par.resolve(strict=True)
        if defaults.is_relative_to(cmd.parent):
            defaults = work / defaults.relative_to(cmd.parent)
        command.extend(["--default-par", str(defaults)])
    with (work / "ctest-sd-sim.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(
            [*command, cmd.name],
            cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        print((work / "ctest-sd-sim.log").read_text(encoding="utf-8", errors="replace")[-6000:],
              file=sys.stderr)
        return 1
    qualify(args.case, work / cmd.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
