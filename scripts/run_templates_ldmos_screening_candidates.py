"""Serial, same-program screening-root controls using explicitly staged D0 inputs.

The input directory contains vg4.json/vg8.json, deck_vg4.json/deck_vg8.json,
mesh.json and ialmob_geometry.json. Neutral initialization is required: archived
seed references are removed because this protocol computes initialization anew.
No external/native qualification is inferred from baseline agreement.
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time
import traceback

from electrothermal_state import read_bound_record
from run_templates_ldmos_electrothermal_curve import state_gate


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--points", type=int, choices=(8, 31), default=8)
    args = parser.parse_args()
    runner, inputs, output = (p.resolve() for p in (args.runner, args.inputs, args.output))
    output.mkdir(parents=True, exist_ok=False)
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
                      VELA_LINEAR_THREADS="1", VELA_BLAS_THREADS="1")
    frozen = {str(p): digest(p) for p in [runner, Path(__file__), *sorted(inputs.glob("*.json"))]}
    status = dict(status="running", pid=os.getpid(), started=time.time(), points=args.points,
                  frozen_sha256=frozen, runs=[], comparisons=[], methods_rejected=[])
    publish = lambda: save(output/"summary.json", status)
    publish()
    try:
        fields = ("state_interleaved", "referenced_state_interleaved", "electron_qf_reference_V",
                  "hole_qf_reference_V", "temperature_K")
        for gate in (4, 8):
            baseline = None
            methods = ("legacy", "newton", "halley", "toms748") if gate == 4 else ("legacy", "toms748", "halley", "newton")
            for method in methods:
                name = f"vg{gate}_{method}"
                folder = output/name
                folder.mkdir()
                deck = read(inputs/f"deck_vg{gate}.json")
                if deck["initialization"]["mode"] != "neutral_300K":
                    raise ValueError("This protocol requires neutral_300K initialization")
                cfg = copy.deepcopy(read(inputs/f"vg{gate}.json"))
                for key in (*fields, "state_archive"):
                    cfg.pop(key, None)
                cfg["mesh_file"] = str(inputs/"mesh.json")
                cfg["mobility_SI"]["ialmob"]["geometry_file"] = str(inputs/"ialmob_geometry.json")
                cfg["diagnostic_ialmob_screening_method"] = method
                cfg["diagnostic_electrothermal_cost"] = "off"
                save(folder/"input.json", cfg)
                deck.update(input_file=str(folder/"input.json"), output_directory=str(folder/"results"),
                            resume=False, pause_after_attempts=0)
                deck["sweep"]["bias_points_V"] = deck["sweep"]["bias_points_V"][:args.points]
                save(folder/"deck.json", deck)
                status.update(current=name)
                start = time.perf_counter()
                with (folder/"run.log").open("w") as log:
                    child = subprocess.Popen([str(runner), "--config", str(folder/"deck.json")],
                                             cwd=folder, stdout=log, stderr=subprocess.STDOUT)
                    status["child_pid"] = child.pid
                    publish()
                    rc = child.wait()
                elapsed = time.perf_counter()-start
                status.pop("child_pid", None)
                if rc:
                    raise RuntimeError(f"{name}: return code {rc}")
                ledger = read(folder/"results/ledger.json")
                if ledger["status"] != "complete" or len(ledger["exact_points"]) != args.points:
                    raise ValueError(f"{name}: incomplete curve")
                totals = {}
                for section, records in (("initialization", ledger.get("initialization_runs", [])), ("drain", ledger["runs"])):
                    sums = dict(updates=0, attempts=len(records), rejected_attempts=0, performance={})
                    for entry in records:
                        result = read_bound_record(entry["result"] if section == "initialization" else Path(entry["directory"])/"output.json")
                        sums["updates"] += result["newton_updates"]
                        sums["rejected_attempts"] += not entry.get("gate", {}).get("pass_gate", True)
                        perf = result["performance"]
                        if perf["linear_solver"] != "umfpack" or perf["ialmob_screening_method"] != method:
                            raise ValueError("Unexpected backend or screening method")
                        for key in ("assembly_seconds", "assembly_calls", "factorization_seconds", "factorizations",
                                    "ialmob_screening_candidate_calls", "ialmob_screening_function_evaluations", "ialmob_screening_fallbacks"):
                            sums["performance"][key] = sums["performance"].get(key, 0)+perf.get(key, 0)
                        if section == "drain" and entry["gate"]["pass_gate"]:
                            if not state_gate(result, entry["bias_V"])["pass_gate"]:
                                raise ValueError("Accepted continuation state failed original gate")
                    totals[section] = sums
                states = []
                for point in ledger["exact_points"]:
                    result = read_bound_record(point["result"])
                    if not state_gate(result, point["bias_V"])["pass_gate"]:
                        raise ValueError("Exact point failed original gate")
                    states.append(result)
                row = dict(name=name, gate=gate, method=method, wall_seconds=elapsed, totals=totals,
                           original_point_gates_pass=True)
                status["runs"].append(row)
                if baseline is None:
                    baseline = states
                else:
                    for index, (old, new) in enumerate(zip(baseline, states)):
                        errors = {}
                        for key in fields:
                            if len(old[key]) != len(new[key]):
                                raise ValueError("State size changed")
                            errors[key] = max(abs(x-y) for x, y in zip(old[key], new[key]))
                        passed = all(v <= 1e-8 for v in errors.values())
                        status["comparisons"].append(dict(gate=gate, method=method, point=index,
                                                         maximum_absolute_errors=errors, passed=passed))
                        if not passed and method not in status["methods_rejected"]:
                            status["methods_rejected"].append(method)
                publish()
        if any(digest(Path(p)) != h for p, h in frozen.items()):
            raise ValueError("Frozen executable, inputs or driver changed")
        status.update(status="completed", finished=time.time(),
                      qualification="pass" if not status["methods_rejected"] else "candidate_rejected")
        publish()
    except BaseException as error:
        status.update(status="failed", error=repr(error), traceback=traceback.format_exc(), finished=time.time())
        publish()
        raise


if __name__ == "__main__":
    main()
