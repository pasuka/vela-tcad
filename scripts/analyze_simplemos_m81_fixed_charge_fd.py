"""Independently compare native current IFM against signed bulk-charge DC solves."""
import argparse
import csv
import json
import math
from pathlib import Path

from run_simplemos_m81_output_amendment import REPO, ROOT, LOCAL, OUT, rows, sha, write
from analyze_simplemos_m81_native_pilot import verify_hashes


def compare(plus, minus, zero, native, amplitude, pilot_current):
    central = (plus - minus) / 2
    predicted = native * amplitude / 1e13
    if predicted == 0 or central == 0:
        raise ValueError("Zero calibration signal")
    drift = abs(zero - pilot_current)
    return {"amplitude_cm_3": amplitude, "plus_current_A_per_um": plus,
            "minus_current_A_per_um": minus, "zero_current_A_per_um": zero,
            "central_delta_A_per_um": central, "ifm_delta_A_per_um": predicted,
            "relative_error": abs(central - predicted) / abs(predicted),
            "same_sign": central * predicted > 0,
            "positive_perturbation_sign_pass": (plus-zero)*predicted > 0,
            "negative_perturbation_sign_pass": (minus-zero)*predicted < 0,
            "even_nonlinear_fraction": abs((plus + minus) / 2 - zero) / abs(central),
            "zero_control_drift_A_per_um": drift,
            "signal_to_zero_control_drift": abs(central) / max(drift, 1e-300)}


def analyze():
    verify_hashes(ROOT / "simplemos_m81_fixed_charge_fd_freeze_v1.json")
    contract = json.loads((ROOT / "simplemos_m81_fixed_charge_fd_contract_v1.json").read_text())
    pilot = json.loads((OUT / "m81_native_pilot_result.json").read_text())
    currents = {}
    executions = []
    for case in contract["cases"]:
        device, name = case["device"], case["case"]
        source = LOCAL / "fixed_charge_bundle" / device / name
        raw = LOCAL / "fixed_charge_raw" / device / name
        for p in source.iterdir():
            if sha(p) != sha(raw / p.name):
                raise ValueError("FD remote input changed")
        exit_code = int((raw / "exit_code.txt").read_text())
        if exit_code != 0:
            raise ValueError(f"FD failed: {device}/{name}")
        dc = rows(raw / "m81_des.plt")
        if len(dc) != 1:
            raise ValueError("Expected exactly one FD DC point")
        r = dc[0]
        if abs(r["gate OuterVoltage"] - .9) > 1e-12 or abs(r["drain OuterVoltage"] - .05) > 1e-12:
            raise ValueError("FD bias changed")
        log = (raw / "m81.console.log").read_text()
        if "T-2022.03-SP2" not in log or "Good Bye" not in log:
            raise ValueError("Wrong runtime or incomplete log")
        current = r["drain TotalCurrent"]
        currents[(device, name)] = current
        executions.append({**case, "exit_code": exit_code, "current_A_per_um": current})
    results = []
    for p in pilot["cases"]:
        device = p["device"]
        zero = currents[(device, "zero")]
        block = []
        for suffix, amplitude in (("1e13", 1e13), ("5e12", 5e12)):
            result = compare(currents[(device, "plus_" + suffix)], currents[(device, "minus_" + suffix)],
                             zero, p["native_dI_A"], amplitude, p["dc_current_A_per_um"])
            result.update({"device": device,
                           "zero_control_reclosure_dex": abs(math.log10(zero / p["dc_current_A_per_um"]))})
            block.append(result)
        convergence = abs((block[1]["central_delta_A_per_um"] / 5e12) /
                          (block[0]["central_delta_A_per_um"] / 1e13) - 1)
        gates = contract["gates"]
        for result in block:
            result["two_amplitude_derivative_relative_change"] = convergence
            result["pass"] = (result["relative_error"] <= gates["fd_vs_ifm_relative_error"] and
                result["same_sign"] and result["positive_perturbation_sign_pass"] and result["negative_perturbation_sign_pass"] and
                result["zero_control_reclosure_dex"] <= gates["zero_control_reclosure_dex"] and
                result["signal_to_zero_control_drift"] >= gates["minimum_signal_to_zero_control_drift"] and
                convergence <= gates["two_amplitude_derivative_relative_change"])
        results.extend(block)
    passed = all(r["pass"] for r in results)
    for filename, data in (("m81_fd_case_ledger.csv", executions), ("m81_fd_calibration_ledger.csv", results)):
        path = OUT / filename
        if path.exists():
            raise FileExistsError(path)
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]))
            writer.writeheader(); writer.writerows(data)
    report = {"status": "qualified_native_uniform_charge_response" if passed else "failed_fd_qualification",
        "passed": passed, "calibration_rows": len(results), "passed_rows": sum(r["pass"] for r in results),
        "max_relative_error": max(r["relative_error"] for r in results),
        "max_two_amplitude_derivative_change": max(r["two_amplitude_derivative_relative_change"] for r in results),
        "sign_convention": "Compare native dI directly to device drain TotalCurrent difference; AC i(,N_drain) has the opposite sign and is not substituted for dI.",
        "units": "Native AC dI in A for the default 1 um 2D device depth; numerically equals the device current difference in A/um. Fixed charge density in cm^-3.",
        "scope": "Qualifies uniform Silicon fixed-charge integrated response at n19/n23 Vd=.05 V Vg=.9 V only. This FD comparison alone does not qualify raw Green volume weighting (see the separate n23 native Green integral audit), local source shapes, high Vd, or a Vela physics correction.",
        "counts_this_execution": {"initial_dc_reclosures": 2, "amendment_dc_reclosures": 2,
            "initial_ac_ifm_points": 2, "amendment_ac_ifm_points": 2,
            "perturbed_dc_solves": 8, "zero_charge_control_dc_solves": 2, "new_bias_sweeps": 0},
        "m82_released": False, "m83_released": False,
        "supersedes_status": "sentaurus_ifm_pilot/m81_execution_status.json: previous upload block resolved by explicit user authorization; historical file retained."}
    write(OUT / "m81_fd_calibration_result.json", report)
    evidence_files = [Path(__file__), ROOT / "simplemos_m81_fixed_charge_fd_contract_v1.json",
                      ROOT / "simplemos_m81_fixed_charge_fd_freeze_v1.json"]
    evidence_files += [OUT / n for n in ("m81_fd_case_ledger.csv", "m81_fd_calibration_ledger.csv", "m81_fd_calibration_result.json")]
    evidence_files += [p for p in (LOCAL / "fixed_charge_raw").rglob("*") if p.is_file()]
    write(OUT / "m81_fd_execution_evidence.json", {"status": report["status"], "input_hashes": {
        str(p.relative_to(REPO)).replace("\\", "/"): sha(p) for p in evidence_files}})
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        verify_hashes(ROOT / "simplemos_m81_output_amendment_freeze_v1.json")
        verify_hashes(ROOT / "simplemos_m81_fixed_charge_fd_freeze_v1.json")
        verify_hashes(OUT / "m81_fd_execution_evidence.json")
        print((OUT / "m81_fd_calibration_result.json").read_text())
    else:
        analyze()
