"""Reconcile audit identity, strict failure and paired current shifts without overclaiming."""
import argparse
import math
from pathlib import Path

import run_simplemos_convergence_audit as audit
import run_simplemos_strict_convergence_contrast as first
import run_simplemos_convergence_acceptance_isolation as isolation

REPO, ROOT = audit.REPO, audit.ROOT
OUT = ROOT / "convergence_acceptance_review"
EVIDENCE = ROOT / "simplemos_convergence_acceptance_evidence.json"


def paired_change(low_shift, high_shift, low_qualified, high_qualified):
    return {"paired_error_change_dex": high_shift-low_shift,
            "qualified_pair": bool(low_qualified and high_qualified)}


def summarize():
    audit.verify(); audit.verify(first.FREEZE); audit.verify(isolation.FREEZE)
    states = audit.rows(audit.OUT/"state_audit.csv")
    groups = []
    for branch in ("m65", "m77"):
        for stage in ("endpoint", "stall"):
            subset = [r for r in states if r["branch"] == branch and r["stage"] == stage]
            valid = [r for r in subset if r["identity_pass"] == "True"]
            groups.append({"branch": branch, "stage": stage, "states": len(subset), "identity_pass": len(valid),
                "local_pass_on_valid_states": sum(r["local_pass"] == "True" for r in valid),
                "max_kcl_over_Id_on_valid_states": max(float(r["kcl_over_Id"]) for r in valid) if valid else None})
    audit.write_csv(OUT/"audit_identity_groups.csv", groups)
    sent = {(r["device"], float(r["drain_voltage_V"])): r
            for r in audit.rows(ROOT/"electron_poisson_charge_volume/m74_case_ledger.csv")}
    allcases = audit.rows(first.OUT/"case_ledger.csv")
    allcases += [dict(r, arm="acceptance_only") for r in audit.rows(isolation.OUT/"case_ledger.csv")]
    paired = []
    for arm in ("original", "strict", "acceptance_only"):
        for vd in (.05, 1.):
            case = {r["device"]: r for r in allcases if r["arm"] == arm and float(r["vd"]) == vd}
            errors, original = {}, {}
            for device in ("n19", "n23"):
                sr = sent[(device, vd)]
                if abs(float(sr["diagnostic_gate_voltage_V"])-float(case[device]["vg"])) > 1e-10:
                    raise ValueError("Sentaurus comparison bias changed")
                errors[device] = math.log10(abs(float(case[device]["current_A_per_um"])/float(sr["sentaurus_current_A_per_um"])))
                original[device] = float(sr["baseline_error_dex"])
            paired.append({"arm": arm, "vd": vd,
                "baseline_low_error_dex": original["n19"], "baseline_high_error_dex": original["n23"],
                "observed_low_error_dex": errors["n19"], "observed_high_error_dex": errors["n23"],
                "baseline_pair_error_dex": original["n23"]-original["n19"],
                "observed_pair_error_dex": errors["n23"]-errors["n19"],
                **paired_change(float(case["n19"]["delta_log10_Id_dex"]), float(case["n23"]["delta_log10_Id_dex"]),
                                case["n19"]["qualified"] == "True", case["n23"]["qualified"] == "True")})
    audit.write_csv(OUT/"paired_current_observations.csv", paired)
    state_identities = []
    for c in audit.read(first.CONTRACT)["cases"]:
        original = first.LOCAL/c["case"]/"original/state.csv"
        isolated = isolation.LOCAL/c["case"]/"state.csv"
        state_identities.append({"case": c["case"], "original_sha256": audit.sha(original),
                                 "acceptance_only_sha256": audit.sha(isolated),
                                 "identical": audit.sha(original) == audit.sha(isolated)})
    audit.write(OUT/"control_state_identity.json", state_identities)
    results = audit.rows(isolation.OUT/"case_ledger.csv")
    audit.write(OUT/"review_summary.json", {
        "status": "requested_audit_and_minimal_contrast_executed_strict_qualification_incomplete",
        "read_only_states": 64, "endpoint_identities_passed": 32, "endpoint_count": 32,
        "stall_identities_failed": 19, "stall_count": 32,
        "strict_initial_passed": 0, "acceptance_only_passed": sum(r["qualified"] == "True" for r in results),
        "acceptance_only_count": 4, "qualified_high_low_pairs": sum(r["qualified_pair"] for r in paired if r["arm"] == "acceptance_only"),
        "new_nonlinear_reclosures": 12, "new_sentaurus_runs": 0, "new_bias_sweeps": 0,
        "max_observed_original_control_shift_dex": max(abs(float(r["delta_log10_Id_dex"])) for r in allcases if r["arm"] == "original"),
        "control_and_acceptance_only_state_hashes_identical": all(r["identical"] for r in state_identities),
        "limitations": ["VTK-reconstructed states failing current identity are excluded from quantitative conclusions about the original saved solve.",
            "Global source-relative gates are unqualified below the explicit source floor; satisfied=true is not source-relative closure evidence.",
            "Failed strict states and small observed Id movement cannot establish an upper bound on the fully converged correction.",
            "global_continuity_closure report/enforce both alter Newton line-search merit in the current implementation."],
        "m79_m81_historical_evidence_preserved": True, "m82_released": False, "m83_released": False})


def freeze():
    audit.verify(); audit.verify(first.FREEZE); audit.verify(isolation.FREEZE)
    paths = list((REPO/"scripts").glob("*simplemos*convergence*.py"))
    paths += list((REPO/"tests/regression").glob("test_simplemos_convergence*.py"))
    paths += list(ROOT.glob("simplemos*convergence*contract*.json")) + list(ROOT.glob("simplemos*convergence*freeze*.json"))
    paths += [REPO/"docs/validation/simplemos_convergence_audit_and_strict_contrast_2026-09-05.md"]
    for folder in (audit.OUT, first.OUT, isolation.OUT, OUT):
        paths += [p for p in folder.rglob("*") if p.is_file()]
    for folder in (audit.LOCAL, first.LOCAL, isolation.LOCAL):
        paths += [p for p in folder.rglob("*") if p.is_file()]
    audit.write(EVIDENCE, {"status": "frozen_completed_execution_with_failed_strict_gates",
        "input_hashes": {audit.rel(p): audit.sha(p) for p in sorted(set(paths))}})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__); p.add_argument("action", choices=("summarize", "freeze", "verify"))
    a = p.parse_args(); {"summarize": summarize, "freeze": freeze, "verify": lambda: audit.verify(EVIDENCE)}[a.action]()
