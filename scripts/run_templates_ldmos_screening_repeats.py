"""Qualify four full screening-root curves, then repeat the best candidate.

Consumes the frozen first-stage batch without altering it. Native temperature
fields must already have portable paths and verified hashes. Three fresh rounds
are timed serially; first-stage timings select a candidate but are not repeats.
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import time
import traceback

from analyze_templates_ldmos_d0 import analyze
from electrothermal_state import read_bound_record
from run_templates_ldmos_electrothermal_curve import state_gate
from run_templates_ldmos_screening_candidates import read, save, digest

FIELDS = ("state_interleaved", "referenced_state_interleaved", "electron_qf_reference_V",
          "hole_qf_reference_V", "temperature_K")
METHODS = ("legacy", "newton", "halley", "toms748")


def select_candidate(batch, native_status):
    if batch.get("status") != "completed" or batch.get("points") != 31:
        raise ValueError("Completed full first-stage batch required")
    keys = [(r["method"], r["gate"]) for r in batch["runs"]]
    if len(keys) != 8 or set(keys) != {(m, g) for m in METHODS for g in (4, 8)}:
        raise ValueError("Exactly eight distinct first-stage curves required")
    if native_status.get("legacy") != "pass":
        raise ValueError("Legacy native qualification failed")
    if not all(r["original_point_gates_pass"] for r in batch["runs"]):
        raise ValueError("Original point gates failed")
    eligible = []
    for method in METHODS[1:]:
        comparisons = [c for c in batch["comparisons"] if c["method"] == method]
        keys = [(c["gate"], c["point"]) for c in comparisons]
        complete = len(keys) == 62 and set(keys) == {(g, p) for g in (4, 8) for p in range(31)}
        if complete and all(c["passed"] for c in comparisons) and native_status.get(method) == "pass":
            wall = sum(r["wall_seconds"] for r in batch["runs"] if r["method"] == method)
            if not math.isfinite(wall) or wall <= 0:
                raise ValueError("Invalid first-stage time")
            eligible.append((wall, method))
    if not eligible:
        raise ValueError("No qualified candidate")
    return min(eligible)[1]


def paired_order(round_index, gate, candidate):
    return (candidate, "legacy") if (round_index+gate//4) % 2 else ("legacy", candidate)


def compare_states(left, right):
    errors = {}
    for key in FIELDS:
        a, b = left[key], right[key]
        if not a or len(a) != len(b) or not all(math.isfinite(v) for v in (*a, *b)):
            raise ValueError("Invalid state vector: "+key)
        errors[key] = max(abs(x-y) for x, y in zip(a, b))
    if any(v > 1e-8 for v in errors.values()):
        raise ValueError("Original state agreement gate failed: "+str(errors))
    return errors


def inspect_curve(folder, method):
    ledger = read(folder/"results/ledger.json")
    if ledger["status"] != "complete" or len(ledger["exact_points"]) != 31:
        raise ValueError("Incomplete repeated curve")
    totals = {}
    for section, records in (("initialization", ledger.get("initialization_runs", [])), ("drain", ledger["runs"])):
        sums = dict(updates=0, attempts=len(records), rejected_attempts=0, trials=0, performance={})
        for entry in records:
            path = entry["result"] if section == "initialization" else Path(entry["directory"])/"output.json"
            result = read_bound_record(path)
            sums["updates"] += result["newton_updates"]
            sums["trials"] += sum(h.get("line_search_trials", 0) for h in result.get("history", []))
            sums["rejected_attempts"] += not entry.get("gate", {}).get("pass_gate", True)
            perf = result["performance"]
            if perf["linear_solver"] != "umfpack" or perf["ialmob_screening_method"] != method:
                raise ValueError("Unexpected backend or root method")
            for key in ("assembly_seconds", "assembly_calls", "factorization_seconds", "factorizations",
                        "ialmob_screening_candidate_calls", "ialmob_screening_function_evaluations", "ialmob_screening_fallbacks"):
                sums["performance"][key] = sums["performance"].get(key, 0)+perf[key]
            if section == "drain" and entry["gate"]["pass_gate"] and not state_gate(result, entry["bias_V"])["pass_gate"]:
                raise ValueError("Accepted continuation failed original gates")
        totals[section] = sums
    states = []
    for point in ledger["exact_points"]:
        result = read_bound_record(point["result"])
        if not state_gate(result, point["bias_V"])["pass_gate"]:
            raise ValueError("Exact point failed original gates")
        states.append(result)
    return ledger, states, totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("runner", "qualification", "native", "contract", "output"):
        parser.add_argument("--"+key, type=Path, required=True)
    args = parser.parse_args()
    runner, qualification, native, contract, output = (getattr(args, k).resolve() for k in
                                                       ("runner", "qualification", "native", "contract", "output"))
    output.mkdir(parents=True, exist_ok=False)
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", OMP_DYNAMIC="FALSE",
                      VELA_LINEAR_THREADS="1", VELA_BLAS_THREADS="1")
    batch = read(qualification/"summary.json")
    frozen = {p: h for p, h in batch["frozen_sha256"].items()}
    if str(runner) not in frozen:
        raise ValueError("Runner must be the qualified executable")
    for p in [Path(__file__), qualification/"summary.json", contract, *native.rglob("*.json"), *native.rglob("*.csv"),
              *qualification.glob("vg*/input.json"), *qualification.glob("vg*/deck.json"),
              *qualification.glob("vg*/results/ledger.json"),
              *Path(__file__).parent.glob("analyze_templates_ldmos*.py"),
              Path(__file__).parent/"electrothermal_state.py", Path(__file__).parent/"state_archive.py"]:
        frozen[str(p)] = digest(p)
    def check_frozen():
        for p, h in frozen.items():
            if digest(Path(p)) != h:
                raise ValueError("Frozen file changed: "+p)
    check_frozen()
    status = dict(status="running", stage="native_qualification", pid=os.getpid(), started=time.time(),
                  runs=[], comparisons=[], joint={}, frozen_sha256=frozen)
    publish = lambda: save(output/"summary.json", status)
    publish()
    try:
        native_status = {}
        for method in METHODS:
            status["current"] = "qualify_"+method
            publish()
            result = analyze(native, {g: qualification/f"vg{g}_{method}/results" for g in (4, 8)}, read(contract), 15)
            save(output/f"qualification_{method}.json", result)
            native_status[method] = result["status"]
            status["joint"]["qualification_"+method] = result["status"]
        candidate = select_candidate(batch, native_status)
        status.update(candidate=candidate, stage="repeats", selection="minimum qualified first-stage dual-gate total wall")
        publish()
        for repetition in range(3):
            for gate in (4, 8):
                pair = {}
                for method in paired_order(repetition, gate, candidate):
                    check_frozen()
                    if shutil.disk_usage(output).free < 2*1024**3:
                        raise RuntimeError("Less than 2 GiB disk space remains")
                    name = f"r{repetition}_vg{gate}_{method}"
                    folder = output/name
                    folder.mkdir()
                    original = qualification/f"vg{gate}_{method}"
                    cfg = read(original/"input.json")
                    if cfg["diagnostic_ialmob_screening_method"] != method:
                        raise ValueError("Input root-method mismatch")
                    save(folder/"input.json", cfg)
                    deck = read(original/"deck.json")
                    deck.update(input_file=str(folder/"input.json"), output_directory=str(folder/"results"), resume=False)
                    save(folder/"deck.json", deck)
                    status["current"] = name
                    cpu_before = None
                    if os.name == "posix":
                        import resource
                        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
                        cpu_before = usage.ru_utime+usage.ru_stime
                    load_before = os.getloadavg() if hasattr(os, "getloadavg") else None
                    start = time.perf_counter()
                    with (folder/"run.log").open("w") as log:
                        child = subprocess.Popen([str(runner), "--config", str(folder/"deck.json")], cwd=folder,
                                                 stdout=log, stderr=subprocess.STDOUT)
                        status["child_pid"] = child.pid
                        publish()
                        rc = child.wait()
                    elapsed = time.perf_counter()-start
                    cpu_seconds = None
                    if cpu_before is not None:
                        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
                        cpu_seconds = usage.ru_utime+usage.ru_stime-cpu_before
                    status.pop("child_pid", None)
                    if rc:
                        raise RuntimeError(f"{name}: exit {rc}")
                    ledger, states, totals = inspect_curve(folder, method)
                    old = read(original/"results/ledger.json")
                    for p, point, state in zip(old["exact_points"], ledger["exact_points"], states):
                        if p["bias_V"] != point["bias_V"]:
                            raise ValueError("Repeated bias mismatch")
                        status["comparisons"].append(dict(scope="repeat_vs_qualification", name=name, bias_V=p["bias_V"],
                                                         errors=compare_states(read_bound_record(p["result"]), state)))
                    trajectory = lambda x: [(v["bias_V"], v["gate"]["pass_gate"]) for v in x["runs"]]
                    status["runs"].append(dict(name=name, repetition=repetition, gate=gate, method=method,
                                               wall_seconds=elapsed, cpu_seconds=cpu_seconds, load_before=load_before,
                                               totals=totals, trajectory_equal=trajectory(old) == trajectory(ledger)))
                    pair[method] = states
                    publish()
                for index, (left, right) in enumerate(zip(pair["legacy"], pair[candidate])):
                    status["comparisons"].append(dict(scope="paired", repetition=repetition, gate=gate, point=index,
                                                     errors=compare_states(left, right)))
                publish()
            for method in ("legacy", candidate):
                result = analyze(native, {g: output/f"r{repetition}_vg{g}_{method}/results" for g in (4, 8)}, read(contract), 15)
                name = f"r{repetition}_{method}"
                save(output/f"joint_{name}.json", result)
                status["joint"][name] = result["status"]
                publish()
                if result["status"] != "pass":
                    raise ValueError("Repeated native joint qualification failed: "+name)
            check_frozen()
        status.update(status="completed", stage="qualified", finished=time.time())
        publish()
    except BaseException as error:
        status.update(status="failed", error=repr(error), traceback=traceback.format_exc(), finished=time.time())
        publish()
        raise


if __name__ == "__main__":
    main()
