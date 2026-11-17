#!/usr/bin/env python3
"""Run the U-2022.12 Schottky and BJT M1 numerical regressions."""

from __future__ import annotations

import argparse
from copy import deepcopy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time

from compare_genius_bjt_sentaurus_vela import compare_model
from run_schottky_reference_workflow import audit_kcl, audit_stage_b, compare_curves, read_csv


REFERENCE_ROOT = Path(__file__).resolve().parents[1] / "reference_tcad"
CASES = {
    "schottky_iv": "schottky_charon_sentaurus2018",
    "bjt_m1": "genius_bjt_sentaurus2022",
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


def run_stages(runner, fixture, output, env, report):
    config = read_json(fixture / "simulation.json")
    simulation = config["simulation"]
    for field in ("mesh_file", "node_doping_file", "materials_file"):
        simulation[field] = str((fixture / simulation[field]).resolve(strict=True))
    contacts = {entry["name"]: entry for entry in simulation["contacts"]}
    state = None
    for stage in config["stages"]:
        name = stage["name"]
        deck = deepcopy(simulation)
        deck["solver"]["max_iter"] = stage.get("max_iter", deck["solver"]["max_iter"])
        sweep = deck["sweep"]
        sweep.update(stage["sweep"])
        contact = contacts[sweep["contact"]]
        sweep["start"] = contact["bias"]
        sweep["diagnostics"] = {
            "terminal_balance": {
                "enabled": True,
                "contacts": list(contacts),
                "csv_file": str(output / f"{name}-terminal.csv"),
            }
        }
        if state is not None:
            sweep["initial_state_file"] = str(state)
        state = output / f"{name}-final.h5"
        sweep["write_state_file"] = str(state)
        curve = output / f"{name}.csv"
        deck["output_csv"] = str(curve)
        config_path = output / f"{name}.json"
        write_json(config_path, deck)
        log_path = output / f"{name}-runner.log"
        started = time.monotonic()
        with log_path.open("wb") as log:
            result = subprocess.run(
                [str(runner), "--config", str(config_path)],
                cwd=output, env=env, stdout=log, stderr=subprocess.STDOUT,
            )
        rows = read_csv(curve) if curve.is_file() else []
        accepted = [row for row in rows if row["converged"] == "1"]
        last_bias = float(accepted[-1]["bias_V"]) if accepted else None
        reached = last_bias is not None and math.isclose(last_bias, sweep["stop"], rel_tol=0.0, abs_tol=1e-9)
        report["stages"].append({
            "stage": name, "exit_code": result.returncode, "seconds": time.monotonic() - started,
            "accepted": len(accepted), "failed_rows": len(rows) - len(accepted),
            "goal_V": sweep["stop"], "last_accepted_V": last_bias,
        })
        if result.returncode or not reached or not state.is_file():
            tail = "\n".join(log_path.read_text(errors="replace").splitlines()[-12:])
            raise RuntimeError(f"{name} failed to finish at {sweep['stop']} V. See {log_path}\n{tail}")
        contact["bias"] = sweep["stop"]
        write_json(output / "result.json", report)


def schottky_gates(fixture, output, report):
    rows = [row for row in read_csv(output / "forward.csv") if row["converged"] == "1"]
    points = [
        {key: float(row[key]) for key in ("bias_V", "current_total_A_per_um")}
        for row in rows
    ]
    curve = compare_curves(points, fixture / "sentaurus_reference.csv")
    branch = audit_stage_b(rows)
    components = audit_kcl(rows)
    terminals = {}
    for row in read_csv(output / "forward-terminal.csv"):
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


def bjt_gates(fixture, output, report):
    contract = read_json(fixture.parent / "contracts" / "comparison_thresholds.json")
    parity = contract["wp3_wp5_vela_comparison"]["numerical_parity"]
    points, summary = compare_model(
        "M1", fixture / "sentaurus_reference.csv", output / "collector-terminal.csv",
        tuple(parity["active_region_VCE_range_V"]), parity["models"]["M1"],
    )
    with (output / "comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(points[0]))
        writer.writeheader()
        writer.writerows(points)
    report["comparison"] = summary
    return summary["operational_pass"] and summary["numerical_parity_gate"]["pass"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=CASES)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--linear-solver", required=True, choices=("pardiso_mkl", "umfpack", "sparselu"))
    args = parser.parse_args()
    runner = args.runner.resolve(strict=True)
    fixture = REFERENCE_ROOT / CASES[args.case] / "ci_u2022_12"
    args.output.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix=f"{args.case}-", dir=args.output.resolve()))
    env = dict(
        os.environ, VELA_LINEAR_SOLVER=args.linear_solver, VELA_LINEAR_THREADS="1",
        MKL_THREADING_LAYER="SEQUENTIAL", MKL_NUM_THREADS="1", OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
    )
    env.pop("MKL_CBWR", None)
    report = {"case": args.case, "linear_solver": args.linear_solver, "status": "incomplete", "stages": []}
    print(f"Reports: {output}", flush=True)
    try:
        verify_inputs(fixture)
        run_stages(runner, fixture, output, env, report)
        gates = schottky_gates if args.case == "schottky_iv" else bjt_gates
        report["status"] = "passed" if gates(fixture, output, report) else "gate_failed"
    finally:
        write_json(output / "result.json", report)
    print(f"{args.case}: {report['status']}")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
