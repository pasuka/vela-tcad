"""Separate convergence acceptance from the global-closure line-search merit."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import math
from pathlib import Path
import subprocess

import run_simplemos_convergence_audit as audit
import run_simplemos_strict_convergence_contrast as first

REPO, ROOT = audit.REPO, audit.ROOT
LOCAL = REPO / "build-release/simplemos_convergence_acceptance_isolation"
OUT = ROOT / "convergence_acceptance_isolation"
CONTRACT = ROOT / "simplemos_convergence_acceptance_isolation_contract_v1.json"
FREEZE = ROOT / "simplemos_convergence_acceptance_isolation_freeze_v1.json"


def prepare():
    audit.verify(first.FREEZE)
    contract = audit.read(first.CONTRACT)
    paths = [Path(__file__), first.CONTRACT, first.FREEZE, first.OUT/"case_ledger.csv", first.OUT/"summary.json",
             REPO/"src/solver/NewtonSolver.cpp", first.PRODUCTION]
    for case in contract["cases"]:
        source = first.LOCAL/case["case"]/"strict/config.json"
        deck = audit.read(source)
        dest = LOCAL/case["case"]
        deck["solver"]["global_continuity_closure"]["mode"] = "off"
        deck["solver"]["diagnostics"] = True
        deck["output_state_file"] = str(dest/"state.csv")
        audit.write(dest/"config.json", deck)
        paths.extend((source, dest/"config.json"))
    audit.write(CONTRACT, {"status": "frozen_before_execution", "cases": contract["cases"],
        "additional_nonlinear_reclosures": 4, "new_bias_sweeps": 0, "sentaurus_runs": 0,
        "trigger": "All four first strict solves failed. Inspection of NewtonSolver.cpp globalClosureLineSearchNorm shows both report and enforce alter line-search merit, even when sources are unqualified. Thus first contrast is not an acceptance-only intervention.",
        "change": "Only global_continuity_closure.mode -> off within Newton, plus diagnostics logging. Keep local enforce at 1e-6, identical state, tolerances, row scaling, damping and iteration limit. Independently recompute the same global closure and KCL after each solve; original release thresholds unchanged.",
        "gates": contract["gates"], "global_profile": audit.read(first.LOCAL/contract["cases"][0]["case"]/"strict/config.json")["solver"]["global_continuity_closure"],
        "no_posthoc_relaxation": True, "m82_released": False})
    paths.append(CONTRACT)
    audit.write(FREEZE, {"input_hashes": {audit.rel(p): audit.sha(p) for p in sorted(set(paths))}})


def run_case(case):
    dest = LOCAL/case["case"]
    process = subprocess.run([str(first.PRODUCTION), "--config", str(dest/"config.json"), "--log", "off"], capture_output=True, text=True)
    (dest/"stdout.txt").write_text(process.stdout); (dest/"stderr.txt").write_text(process.stderr)
    status = json.loads(process.stdout.strip().splitlines()[-1]); status["exit_code"] = process.returncode
    audit.write(dest/"status.json", status)
    deck = audit.read(dest/"config.json"); deck.pop("output_state_file", None)
    deck["state_file"] = str(dest/"state.csv")
    deck["solver"]["global_continuity_closure"]["mode"] = "report"
    deck["solver"]["carrier_row_convergence"]["mode"] = "report"
    checks = audit.execute(dict(deck, simulation_type="newton_carrier_term_probe", output_csv=str(dest/"terms.csv"),
                               carrier_term_probe={"solved_equation_terms": True}), dest/"audit.json")
    currents = {c: audit.execute(dict(deck, simulation_type="terminal_current_functional_probe", contact=c), dest/(c+".json"))
                for c in ("source", "drain", "gate", "substrate")}
    current = currents["drain"]["current_A_per_um"]
    kcl = abs(math.fsum(r["current_A_per_um"] for r in currents.values()))/max(abs(current),1e-20)
    local, glob = checks["carrier_row_convergence"], checks["global_continuity_closure"]
    result = {"case": case["case"], "device": case["device"], "vd": case["vd"], "vg": case["vg"],
        "native_converged": status["converged"], "exit_code": status["exit_code"], "iterations": status["iterations"],
        "reason": status["convergence_reason"], "failure": status["failure_reason"], "current_A_per_um": current,
        "delta_log10_Id_dex": math.log10(abs(current/case["reference_current_A_per_um"])),
        "local_satisfied": local["satisfied"], "local_violations": local["violation_count"], "local_max_ratio": local["max_ratio"],
        "global_satisfied": glob["satisfied"], "electron_global_qualified": glob["electron"]["qualified"],
        "hole_global_qualified": glob["hole"]["qualified"], "kcl_over_Id": kcl, "final_residual": status["final_residual"]}
    result["qualified"] = bool(status["converged"] and local["satisfied"] and glob["satisfied"] and kcl <= 1e-8)
    audit.write(dest/"result.json", result)
    print(case["case"], "qualified", result["qualified"], "local", result["local_violations"], "delta", result["delta_log10_Id_dex"], flush=True)
    return result


def run():
    audit.verify(FREEZE)
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run_case, audit.read(CONTRACT)["cases"]))
    audit.write_csv(OUT/"case_ledger.csv", results)
    audit.write(OUT/"summary.json", {"status": "qualified" if all(r["qualified"] for r in results) else "strict_qualification_incomplete",
        "nonlinear_reclosures": 4, "native_converged": sum(r["native_converged"] for r in results),
        "qualified": sum(r["qualified"] for r in results), "m82_released": False})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("action", choices=("prepare", "run", "verify"))
    a = p.parse_args(); {"prepare": prepare, "run": run, "verify": lambda: audit.verify(FREEZE)}[a.action]()
