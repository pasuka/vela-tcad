"""Seal or verify M81 completion without rewriting prior blocked/failed records."""
import argparse
import json
from pathlib import Path

import run_simplemos_m81_sentaurus_ifm_pilot as pilot
from run_simplemos_m81_output_amendment import REPO, ROOT, LOCAL, OUT, sha, write
from analyze_simplemos_m81_native_pilot import verify_hashes

EVIDENCE = ROOT / "simplemos_m81_native_execution_evidence.json"


def audit():
    pilot.validate()
    for path in (ROOT / "simplemos_m81_output_amendment_freeze_v1.json",
                 ROOT / "simplemos_m81_fixed_charge_fd_freeze_v1.json",
                 ROOT / "simplemos_m81_green_integral_freeze_v1.json",
                 OUT / "m81_fd_execution_evidence.json", OUT / "m81_state_closure_audit.json"):
        verify_hashes(path)
    for name in ("m81_initial_native_result.json", "m81_native_pilot_result.json"):
        data = json.loads((OUT / name).read_text())
        for path, digest in data["raw_hashes"].items():
            if sha(REPO / path) != digest:
                raise ValueError("Native evidence changed")
    for path, digest in json.loads((OUT / "m81_preparation_evidence.json").read_text())["artifacts"].items():
        if sha(REPO / path) != digest:
            raise ValueError("Historical preparation evidence changed")
    for attempt in ("raw", "output_amendment_raw"):
        for device in ("n19", "n23"):
            raw = LOCAL / attempt / device
            if int((raw / "exit_code.txt").read_text()) != 0 or "T-2022.03-SP2" not in (raw / "m81.console.log").read_text():
                raise ValueError("Wrong native runtime or unsuccessful execution")
    final = json.loads((OUT / "m81_fd_calibration_result.json").read_text())
    integral = json.loads((OUT / "m81_native_green_integral_result.json").read_text())
    if not final["passed"] or final["passed_rows"] != 4 or not integral["pass"]:
        raise ValueError("M81 incomplete")
    return final


def freeze():
    final = audit()
    paths = list(ROOT.glob("simplemos_m81*.json")) + list(OUT.glob("m81*"))
    paths += list((REPO / "scripts").glob("*simplemos_m81*.py"))
    paths += list((REPO / "tests/regression").glob("test_simplemos_m81*.py"))
    paths += [REPO / "docs/validation/simplemos_m81_native_ifm_execution_2026-09-05.md"]
    write(EVIDENCE, {"schema": "vela.simplemos.m81.native_execution_evidence.v1",
        "status": "m81_two_case_calibration_complete", "fd_max_relative_error": final["max_relative_error"],
        "green_integral_relative_error": json.loads((OUT / "m81_native_green_integral_result.json").read_text())["relative_error"],
        "m82_released": False, "m83_released": False,
        "input_hashes": {str(p.relative_to(REPO)).replace("\\", "/"): sha(p)
                         for p in sorted(set(paths)) if p.is_file() and p != EVIDENCE}})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.verify:
        audit(); verify_hashes(EVIDENCE)
        print("M81 complete: all native, FD, geometry, state and historical evidence hashes verified.")
    else:
        freeze()
