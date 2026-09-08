"""Freeze the minimal M81 output amendment; preserve the first native run."""
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil

import sentaurus_import as si
import run_simplemos_m81_sentaurus_ifm_pilot as pilot

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "reference_tcad/simplemos_sentaurus2022"
LOCAL = REPO / "build-release/m81_sentaurus_ifm_pilot"
OUT = ROOT / "sentaurus_ifm_pilot"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, obj):
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def rows(path):
    text = path.read_text()
    names = si.parse_quoted_list(text, "datasets")
    if len(names) != len(set(names)):
        raise ValueError("Duplicate dataset")
    result = [dict(zip(names, r)) for r in si.parse_values_block(text, len(names))]
    if not result or any(not math.isfinite(v) for r in result for v in r.values()):
        raise ValueError("Empty or non-finite data")
    return result


def initial_result():
    pilot.validate()
    refs = list(csv.DictReader((ROOT / "electron_poisson_charge_volume/m74_case_ledger.csv").open()))
    results = []
    for device in ("n19", "n23"):
        folder = LOCAL / "raw" / device
        for source in (LOCAL / "bundle" / device).iterdir():
            if sha(source) != sha(folder / source.name):
                raise ValueError("Remote input identity failed")
        dc = rows(folder / "m81_des.plt")
        ac = rows(folder / "m81_ac_ac_des.plt")
        ref = next(r for r in refs if r["device"] == device and float(r["drain_voltage_V"]) == .05)
        reference = float(ref["sentaurus_current_A_per_um"])
        delta = max(abs(math.log10(r["drain TotalCurrent"] / reference)) for r in dc)
        if len(ac) != 1 or ac[0]["frequency"] != 0:
            raise ValueError("Expected one zero-frequency point")
        results.append({"device": device, "exit_code": int((folder / "exit_code.txt").read_text()),
                        "dc_current_A_per_um": dc[-1]["drain TotalCurrent"],
                        "reference_current_A_per_um": reference, "dc_reclosure_dex": delta,
                        "dc_pass": delta <= 1e-5,
                        "native_dI_A": ac[0]["dIuniform_charge(N_drain)"],
                        "green_files": [p.name for p in folder.glob("*acgf*.tdr")]})
    write(OUT / "m81_initial_native_result.json", {
        "status": "failed_missing_green_output", "upload_authorized_by_user": True,
        "counts": {"dc_reclosures": 2, "ac_ifm_points": 2}, "cases": results,
        "raw_hashes": {str(p.relative_to(REPO)).replace("\\", "/"): sha(p)
                       for p in (LOCAL / "raw").rglob("*") if p.is_file()}})


def prepare():
    initial_result()
    bundle = LOCAL / "output_amendment_bundle"
    cases = []
    for device in ("n19", "n23"):
        target = bundle / device
        shutil.copytree(LOCAL / "bundle" / device, target)
        deck = target / "m81_des.cmd"
        text = deck.read_text()
        needle = "  EffectiveIntrinsicDensity(NoBandGapNarrowing)"
        assert text.count(needle) == 1
        deck.write_text(text.replace(needle, needle + "\n  Noise ( DiffusionNoise ( LatticeTemperature ) )"),
                        encoding="utf-8", newline="\n")
        cases.append(device)
    contract = ROOT / "simplemos_m81_output_amendment_contract_v1.json"
    write(contract, {
        "schema": "vela.simplemos.m81.output_amendment.v1", "status": "frozen_before_execution",
        "reason": "Initial DIFM completed and emitted native dI but no spatial Green TDR; unnamed Noise enables plot output (UG p815).",
        "change": "Add unnamed Noise(DiffusionNoise(LatticeTemperature)); retain the exact DC physics, mesh, state, deterministic variation and solver.",
        "manual_pages": [787, 788, 815], "cases": cases,
        "remote_root": "/tmp/vela_simplemos_m81_ifm_20260905/output_amendment_bundle",
        "additional_counts": {"dc_reclosures": 2, "ac_ifm_points": 2, "bias_sweeps": 0},
        "gates": {"dc_reclosure_dex": 1e-5, "native_response_relative_change": 1e-8,
                  "three_finite_nonzero_green_fields": True},
        "m82_released": False})
    inputs = [contract, Path(__file__), OUT / "m81_initial_native_result.json"]
    inputs += [p for p in bundle.rglob("*") if p.is_file()]
    write(ROOT / "simplemos_m81_output_amendment_freeze_v1.json", {
        "status": "frozen_before_execution", "input_hashes": {
            str(p.relative_to(REPO)).replace("\\", "/"): sha(p) for p in inputs}})


if __name__ == "__main__":
    prepare()
