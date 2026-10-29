"""Prepare/run isolated stopping-policy controls from an existing BVMethods deck.

Requires an already qualified mesh/material/physics input and a Release runner.
Generated inputs/results belong in an ignored build directory. No native state
is imported. Each arm starts with its own Poisson initialization.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time


def prepare(baseline: Path, output: Path, stop: float) -> dict[str, Path]:
    original = json.loads(baseline.read_text(encoding="utf-8"))
    arms = {}
    for name, rows, closure, update in (
        ("closure_report", "enforce", "report", "report"),
        ("conservation_report", "report", "report", "report"),
        ("density_update", "report", "report", "enforce"),
    ):
        folder = output / name
        if (folder / "config.json").exists():
            raise FileExistsError(f"Refusing to overwrite an existing control: {folder}")
        folder.mkdir(parents=True, exist_ok=True)
        config = copy.deepcopy(original)
        for key in ("mesh_file", "node_doping_file", "materials_file"):
            source = Path(config[key])
            if not source.is_absolute():
                source = baseline.parent / source
            config[key] = str(source.resolve(strict=True))
        solver, sweep = config["solver"], config["sweep"]
        solver["carrier_row_convergence"]["mode"] = rows
        solver["global_continuity_closure"]["mode"] = closure
        solver["carrier_row_convergence"]["recovery"] = {"mode": "off"}
        sweep.pop("initial_state_file", None)
        sweep.pop("bias_points", None)
        sweep.update(start=0., stop=stop, write_vtk=False)
        sweep["initialization"] = {
            "mode": "poisson_block", "poisson_max_iter": 100,
            "diagnostic_csv": "poisson_init.csv", "write_state_file": "poisson_init.h5",
        }
        circuit = sweep["external_circuit"]
        if circuit["solver"] != "coupled_newton":
            raise ValueError("Baseline must select coupled_newton")
        circuit["coupled_update_convergence"] = {
            "mode": update, "digits": 5, "potential_reference_V": .025852,
            "electron_reference_cm3": 1e10, "hole_reference_cm3": 1e10,
            "contact_reference_V": .025852, "diagnostic_csv": "density_update.csv",
        }
        config["output_csv"] = "points.csv"
        (folder/"config.json").write_text(json.dumps(config, indent=2)+"\n", encoding="utf-8")
        arms[name] = folder
    return arms


def run(runner: Path, folder: Path, timeout: float) -> int:
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               VELA_LINEAR_SOLVER="umfpack")
    receipt = {
        "status": "running", "started_unix": time.time(),
        "runner_sha256": hashlib.sha256(runner.read_bytes()).hexdigest(),
        "config_sha256": hashlib.sha256((folder/"config.json").read_bytes()).hexdigest(),
    }
    with (folder/"run.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen([str(runner), "--config", "config.json"],
                                   cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT)
        receipt["pid"] = process.pid
        (folder/"execution.json").write_text(json.dumps(receipt, indent=2)+"\n")
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait()
            code = 124
    receipt.update(status="finished", exit_code=code,
                   seconds=time.time()-receipt["started_unix"])
    (folder/"execution.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"arm": folder.name, "exit_code": code, "seconds": receipt["seconds"]}), flush=True)
    return code


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-config", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--stop-V", type=float, default=100000.)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--arm", choices=("closure_report", "conservation_report", "density_update"))
    parser.add_argument("--timeout", type=float, default=3600.)
    args = parser.parse_args()
    output = args.output.resolve()
    if args.baseline_config:
        if not (args.stop_V > 0.) or not math.isfinite(args.stop_V):
            parser.error("stop-V must be positive and finite")
        prepare(args.baseline_config.resolve(strict=True), output, args.stop_V)
    if args.arm:
        if args.runner is None:
            parser.error("--arm requires --runner")
        return run(args.runner.resolve(strict=True), output/args.arm, args.timeout)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
