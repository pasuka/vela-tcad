"""Freeze two-amplitude fixed-charge FD after native M81 pilot passes."""
import json
from pathlib import Path
import shutil

from run_simplemos_m81_output_amendment import REPO, ROOT, LOCAL, OUT, sha, write
from analyze_simplemos_m81_native_pilot import verify_hashes


def prepare():
    verify_hashes(ROOT / "simplemos_m81_output_amendment_freeze_v1.json")
    report_path = OUT / "m81_native_pilot_result.json"
    report = json.loads(report_path.read_text())
    if report["status"] != "passed" or not report["fixed_charge_fd_released"]:
        raise ValueError("Native pilot has not passed")
    for name, digest in report["raw_hashes"].items():
        if sha(REPO / name) != digest:
            raise ValueError("Pilot evidence changed")
    bundle = LOCAL / "fixed_charge_bundle"
    cases = []
    for device in ("n19", "n23"):
        for name, charge in (("zero", 0.), ("plus_1e13", 1e13), ("minus_1e13", -1e13),
                             ("plus_5e12", 5e12), ("minus_5e12", -5e12)):
            target = bundle / device / name
            shutil.copytree(LOCAL / "bundle" / device, target)
            deck = target / "m81_des.cmd"
            text = deck.read_text()
            start = text.index("  DeterministicVariation(")
            end = text.index('Physics(Material="Silicon")')
            text = text[:start] + "}\n" + text[end:]
            text = text.replace('Physics(Material="Silicon") {',
                                'Physics(Material="Silicon") {\n  Traps(FixedCharge Conc=' + format(charge, ".16g") + ")")
            text = text.split("NoisePlot {", 1)[0]
            text += '''Solve {
  Load(FilePrefix="state")
  Coupled { Poisson Electron Hole }
  Plot(FilePrefix="reclosed")
}
'''
            deck.write_text(text, encoding="utf-8", newline="\n")
            cases.append({"device": device, "case": name, "fixed_charge_cm_3": charge,
                          "deck": str(deck.relative_to(REPO)).replace("\\", "/")})
    contract = ROOT / "simplemos_m81_fixed_charge_fd_contract_v1.json"
    write(contract, {"schema": "vela.simplemos.m81.fixed_charge_fd.v1", "status": "frozen_before_execution",
        "cases": cases, "remote_root": "/tmp/vela_simplemos_m81_ifm_20260905/fixed_charge_bundle",
        "perturbation": "Uniform bulk Silicon FixedCharge; Conc is signed charge number density in cm^-3, UG pp543-545. No dopant profile, mobility, BGN or SRH parameter change.",
        "baseline": "Keep ImplicitACSystem and exact M81 electrode/Load/DC configuration; remove diagnostic noise and deterministic variations; no AC solve.",
        "amplitude_rationale": "Pilot predicts 5.17e-4 and 7.01e-4 relative Id at 1e13 cm^-3; half-amplitude plus central differences check linearity while staying above measured closure drift.",
        "additional_counts": {"perturbed_dc_solves": 8, "zero_charge_control_dc_solves": 2,
                              "ac_ifm_points": 0, "bias_sweeps": 0},
        "comparison": "central_delta=(Id(+a)-Id(-a))/2; compare directly to native dI*(a/1e13), with no fitted sign or scaling; device depth is the default 1 um.",
        "gates": {"all_exit_zero": True, "zero_control_reclosure_dex": 1e-5,
                  "fd_vs_ifm_relative_error": .05, "two_amplitude_derivative_relative_change": .05,
                  "same_sign": True, "minimum_signal_to_zero_control_drift": 100},
        "m82_released": False})
    inputs = [contract, Path(__file__), report_path, REPO / "scripts/analyze_simplemos_m81_native_pilot.py"]
    inputs += [p for p in bundle.rglob("*") if p.is_file()]
    write(ROOT / "simplemos_m81_fixed_charge_fd_freeze_v1.json", {
        "status": "frozen_before_execution", "input_hashes": {
            str(p.relative_to(REPO)).replace("\\", "/"): sha(p) for p in inputs}})


if __name__ == "__main__":
    prepare()
