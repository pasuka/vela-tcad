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


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def finite_number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite sweep result: {value}")
    return result


def qualify(case, cmd):
    work = cmd.with_name(cmd.stem + "_vela")
    manifest = read_json(work / "run_manifest.json")
    plan = read_json(work / "stage_plan.json")
    if manifest["check_only"] or manifest["exit_code"] != 0:
        raise ValueError("sd-sim did not complete a simulation")
    stages = {item["index"]: item for item in manifest["stages"]}
    if set(stages) != {item["index"] for item in plan["stages"]}:
        raise ValueError("The stage plan and executed stages differ")

    contact = "drain" if case == "bvds" else "gate"
    sweeps = []
    for stage in plan["stages"]:
        run = stages[stage["index"]]
        if run["runner_exit_code"] != 0 or not run["accepted_points"]:
            raise ValueError(f"Stage {stage['index']} did not produce an accepted solution")
        deck = read_json(work / stage["deck"])
        if deck["solver"].get("state_variables") != "potential_density":
            raise ValueError(f"Stage {stage['index']} did not use potential and density unknowns")
        if deck["scaling"]["mode"] != "unit_scaling":
            raise ValueError(f"Stage {stage['index']} did not use unit scaling")
        sweep = stage.get("sweep", {})
        if sweep.get("contact") != contact or sweep.get("mode") != "voltage":
            continue
        if case == "bvds" and deck["sweep"].get("external_circuit", {}).get("mode") != "series_resistor":
            raise ValueError("The BVDS drain sweep must use the external series resistor")
        sweeps.append((stage, run, sweep))

    if not sweeps:
        raise ValueError(f"No {contact} voltage sweep in {cmd.name}")
    for stage, run, sweep in sweeps:
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
        print(f"{case} stage {stage['index']}: {len(rows)} accepted points, final {contact}={end:.12g} V")

    if not any(Path(name).suffix.lower() == ".plt" and (cmd.parent / name).is_file()
               for name in manifest["outputs"]):
        raise ValueError("sd-sim did not publish a current curve")
    if not (cmd.parent / plan["outputs"]["log"]).is_file():
        raise ValueError("sd-sim did not publish its runtime log")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("bvds", "idvg"))
    parser.add_argument("--cmd", type=Path, required=True)
    parser.add_argument("--launcher", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
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
    with (work / "ctest-sd-sim.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(
            [str(launcher), "--exit-on-failure", "--threads", "1", cmd.name],
            cwd=work, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        print((work / "ctest-sd-sim.log").read_text(encoding="utf-8", errors="replace")[-6000:],
              file=sys.stderr)
        return 1
    qualify(args.case, work / cmd.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
