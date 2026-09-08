"""Four-workpoint original-versus-qualified-convergence reclosure; no physics change."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import json
import math
from pathlib import Path

import run_simplemos_convergence_audit as audit

REPO, ROOT = audit.REPO, audit.ROOT
LOCAL = REPO / "build-release/simplemos_strict_convergence"
OUT = ROOT / "strict_convergence_contrast"
CONTRACT = ROOT / "simplemos_strict_convergence_contract_v1.json"
FREEZE = ROOT / "simplemos_strict_convergence_freeze_v1.json"
PRODUCTION = REPO / "build-release/vela_example_runner.exe"


def prepare():
    audit.verify()
    selected = [r for r in audit.rows(audit.OUT/"state_audit.csv")
                if r["branch"] == "m65" and r["stage"] == "endpoint" and r["device"] in ("n19", "n23")]
    if len(selected) != 4 or any(r["identity_pass"] != "True" for r in selected):
        raise ValueError("Four endpoint identities required")
    cases, paths = [], [Path(__file__), PRODUCTION, audit.OUT/"state_audit.csv", audit.FREEZE]
    for source in selected:
        case = {k: source[k] for k in ("case", "device", "source_config", "state")}
        case.update({k: float(source[k]) for k in ("vd", "vg", "reference_current_A_per_um")})
        folder = audit.LOCAL/"cases"/case["case"]
        scale = audit.read(folder/"drain.status.json")["current_scale"]
        case["current_scale"] = scale
        cases.append(case)
        for arm in ("original", "strict"):
            dest = LOCAL/case["case"]/arm
            deck = audit.deck_for(REPO/case["source_config"], REPO/case["state"], case["vd"], case["vg"])
            deck.update(simulation_type="newton_solve_from_state", output_state_file=str(dest/"state.csv"))
            if arm == "strict":
                solver = deck["solver"]
                solver["max_iter"] = 200
                solver["carrier_row_convergence"] = dict(audit.PROFILE, mode="enforce", eps_row=1e-6,
                    min_source_scale=1e-20/scale, min_flux_scale=1e-20/scale, min_newton_max_iter=200)
                solver["global_continuity_closure"] = dict(mode="enforce", tolerance=1e-6, source_floor=1e-18/scale)
            audit.write(dest/"config.json", deck); paths.append(dest/"config.json")
        paths.extend((REPO/case["state"], REPO/case["source_config"], folder/"audit.status.json", folder/"drain.status.json"))
    audit.write(CONTRACT, {"status": "frozen_before_nonlinear_execution", "cases": cases,
        "arms": ["original", "strict"], "counts": {"nonlinear_reclosures": 8, "new_bias_sweeps": 0, "sentaurus_runs": 0},
        "changes": "Same original saved state and physical model in both arms. Strict arm enforces local source/flux row residual <=1e-6 and enables global source-relative gate above a declared absolute relevance floor. max_iter 100->200; original block tolerances and stall floors unchanged. No recovery algorithm, row-scaling, damping, mesh or physical parameter changes.",
        "reviewed_resolution": {"local_flux_or_source_floor_A_per_um": 1e-20,
            "global_source_floor_A_per_um": 1e-18,
            "rationale": "Audit observed net sources 3.8e-21..4.8e-20 A/um and contact-source mismatch 3.4e-20..5.6e-19 A/um. Source-relative global percentages here resolve cancellation noise, not Id accuracy. Global qualification flags must be reported; unqualified is not evidence of source-relative closure. Independent terminal KCL/Id remains mandatory. Local absolute floor suppresses sub-1e-20 A/um flow rows, while retaining observed electron residual violations and the resolvable high-Vd hole rows.",
            "density_screening": "disabled; use current/flux normalization instead of ambiguous internal density units"},
        "gates": {"native_converged": True, "local_satisfied": True, "global_satisfied": True,
                  "terminal_kcl_over_Id": 1e-8, "current_materiality_dex": 1e-4, "pair_materiality_dex": 1e-4},
        "failure_policy": "Keep last accepted failed state for diagnostics only. A failed strict solve cannot establish a converged Id correction or validate M80. Do not relax gates after results.",
        "m82_released": False})
    paths.append(CONTRACT)
    audit.write(FREEZE, {"input_hashes": {audit.rel(p): audit.sha(p) for p in sorted(set(paths))}})


def run_one(item):
    case, arm = item
    dest = LOCAL/case["case"]/arm
    config = dest/"config.json"
    # Frozen config is executed directly, never rewritten by audit.execute.
    import subprocess
    process = subprocess.run([str(PRODUCTION), "--config", str(config), "--log", "off"], capture_output=True, text=True)
    (dest/"stdout.txt").write_text(process.stdout); (dest/"stderr.txt").write_text(process.stderr)
    status = json.loads(process.stdout.strip().splitlines()[-1]) if process.stdout.strip() else {}
    status["exit_code"] = process.returncode
    audit.write(dest/"status.json", status)
    if not (dest/"state.csv").exists():
        raise ValueError("No accepted state returned: " + str(dest))
    deck = audit.read(config)
    deck.pop("output_state_file", None); deck["state_file"] = str(dest/"state.csv")
    deck["solver"]["carrier_row_convergence"] = dict(audit.read(LOCAL/case["case"]/"strict/config.json")["solver"]["carrier_row_convergence"], mode="report")
    deck["solver"]["global_continuity_closure"] = dict(audit.read(LOCAL/case["case"]/"strict/config.json")["solver"]["global_continuity_closure"], mode="report")
    checks = audit.execute(dict(deck, simulation_type="newton_carrier_term_probe", output_csv=str(dest/"terms.csv"),
                               carrier_term_probe={"solved_equation_terms": True}), dest/"audit.json")
    currents = {}
    for contact in ("source", "drain", "gate", "substrate"):
        currents[contact] = audit.execute(dict(deck, simulation_type="terminal_current_functional_probe", contact=contact), dest/(contact+".json"))
    current = currents["drain"]["current_A_per_um"]
    current_extractor = currents["drain"]["contact_current_extractor_A_per_um"]
    kcl = math.fsum(v["current_A_per_um"] for v in currents.values())
    local, glob = checks["carrier_row_convergence"], checks["global_continuity_closure"]
    result = {"case": case["case"], "device": case["device"], "vd": case["vd"], "vg": case["vg"], "arm": arm,
        "exit_code": status["exit_code"], "native_converged": status.get("converged", False),
        "iterations": status.get("iterations"), "convergence_reason": status.get("convergence_reason"),
        "failure_reason": status.get("failure_reason"), "current_A_per_um": current,
        "extractor_A_per_um": current_extractor, "original_current_A_per_um": case["reference_current_A_per_um"],
        "delta_log10_Id_dex": math.log10(abs(current/case["reference_current_A_per_um"])),
        "kcl_over_Id": abs(kcl)/max(abs(current),1e-20), "local_satisfied": local["satisfied"],
        "local_max_ratio": local["max_ratio"], "local_violations": local["violation_count"],
        "qualified_rows": local["qualified_row_count"], "global_satisfied": glob["satisfied"],
        "electron_global_qualified": glob["electron"]["qualified"], "hole_global_qualified": glob["hole"]["qualified"],
        "final_residual": status.get("final_residual")}
    result["qualified"] = bool(status.get("converged") and local["satisfied"] and glob["satisfied"] and result["kcl_over_Id"] <= 1e-8)
    audit.write(dest/"result.json", result)
    print(case["case"], arm, "converged", result["native_converged"], "qualified", result["qualified"], "delta", result["delta_log10_Id_dex"], flush=True)
    return result


def run():
    audit.verify(FREEZE)
    work = [(c, arm) for c in audit.read(CONTRACT)["cases"] for arm in ("original", "strict")]
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run_one, work))
    audit.write_csv(OUT/"case_ledger.csv", results)
    strict = [r for r in results if r["arm"] == "strict"]
    audit.write(OUT/"summary.json", {"status": "qualified" if all(r["qualified"] for r in strict) else "strict_qualification_incomplete",
        "nonlinear_reclosures": 8, "strict_native_converged": sum(r["native_converged"] for r in strict),
        "strict_qualified": sum(r["qualified"] for r in strict),
        "max_observed_strict_current_shift_dex": max(abs(r["delta_log10_Id_dex"]) for r in strict),
        "failed_states_are_not_converged_corrections": True, "m82_released": False})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("action", choices=("prepare", "run", "verify"))
    a = p.parse_args()
    {"prepare": prepare, "run": run, "verify": lambda: audit.verify(FREEZE)}[a.action]()
