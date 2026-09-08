"""Calibrate n23 Green volume weighting using frozen M34 native measures."""
import argparse
import csv
import json
import math
from pathlib import Path

import run_simplemos_m34_sentaurus_interface_box_probe as m34
from run_simplemos_m81_output_amendment import REPO, ROOT, LOCAL, OUT, sha, write
from analyze_simplemos_m81_native_pilot import verify_hashes

CONTRACT = ROOT / "simplemos_m81_green_integral_contract_v1.json"
FREEZE = ROOT / "simplemos_m81_green_integral_freeze_v1.json"
M34 = REPO / "build-release/reference_tcad/simplemos_sentaurus2022/m34_interface_box_probe"
DEBUG = M34 / "sentaurus_raw/sentaurus_bundle/n23/MeasureCoefficients.debug"
M34REPORT = ROOT / "sentaurus_interface_box_probe/m34_sentaurus_interface_box_probe_report.json"


def csvrows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def prepare():
    report = json.loads(M34REPORT.read_text())
    if sha(DEBUG) != report["debug_sha256"] or sha(LOCAL / "bundle/n23/input_fps.tdr") != report["input_tdr_sha256"]:
        raise ValueError("M34/native grid identity changed")
    write(CONTRACT, {"status": "frozen_before_read_only_evaluation", "device": "n23",
        "source": "Reuse existing M34 native region-local MeasureCoefficients; no new geometric export or solve.",
        "formula": "delta_I = q * delta_N_cm^-3 * sum_Silicon_nodes(Gpot_s^-1 * native_Si_measure_um2) * 1e-12; 1 um depth; signed measures; multiply volume exactly once; no fitted sign/scale.",
        "q_C": 1.602176634e-19, "delta_N_cm_3": 1e13, "um3_to_cm3": 1e-12,
        "measure_permutation": report["probe"]["input_to_debug_measure_local_permutation"],
        "gates": {"relative_error_vs_native_dI": .05, "same_sign": True, "coordinate_tolerance_um": 1e-10},
        "m82_released": False})
    paths = [CONTRACT, Path(__file__), DEBUG, M34REPORT, OUT / "m81_native_pilot_result.json"]
    paths += [M34 / "n23_import" / n for n in ("nodes.csv", "elements.csv")]
    paths += [LOCAL / "green_exports/n23" / n for n in ("nodes.csv", "elements.csv", "fields/CurPotReACGreenFunction_region0.csv")]
    write(FREEZE, {"status": "frozen_before_evaluation", "input_hashes": {
        str(p.relative_to(REPO)).replace("\\", "/"): sha(p) for p in paths}})


def audit():
    verify_hashes(FREEZE)
    contract = json.loads(CONTRACT.read_text())
    grid = LOCAL / "green_exports/n23"
    original = csvrows(M34 / "n23_import/nodes.csv")
    final = csvrows(grid / "nodes.csv")
    if len(original) != len(final):
        raise ValueError("Node count changed")
    for a, b in zip(original, final):
        if a["id"] != b["id"] or max(abs(float(a[k])-float(b[k])) for k in ("x_um", "y_um")) > 1e-10:
            raise ValueError("Node identity changed")
    elements = csvrows(M34 / "n23_import/elements.csv")
    if elements != csvrows(grid / "elements.csv"):
        raise ValueError("Element ordering changed")
    measures = m34.parse_debug_block(DEBUG.read_text(), "Measure")
    green = {int(r["node_id"]): float(r["component0"]) for r in csvrows(grid / "fields/CurPotReACGreenFunction_region0.csv")}
    volumes = {node: 0. for node in green}
    perm = contract["measure_permutation"]
    for element in elements:
        if element["material"] != "Si":
            continue
        measure = measures[int(element["id"])]["values"]
        for i in range(3):
            volumes[int(element[f"node{i}"])] += measure[perm[i]]
    integrated = contract["q_C"] * contract["delta_N_cm_3"] * contract["um3_to_cm3"] * math.fsum(green[n]*v for n, v in volumes.items())
    native = next(r["native_dI_A"] for r in json.loads((OUT / "m81_native_pilot_result.json").read_text())["cases"] if r["device"] == "n23")
    error = abs(integrated-native)/abs(native)
    report = {"device": "n23", "integrated_delta_A_per_um": integrated, "native_dI_A": native,
              "relative_error": error, "same_sign": integrated*native > 0,
              "silicon_nodes": len(green), "silicon_measure_sum_um2": sum(volumes.values()),
              "pass": error <= .05 and integrated*native > 0, "new_solves": 0,
              "scope": "n23 uniform Silicon charge only; existing M34 signed native volumes, one volume factor. Does not qualify arbitrary source shapes or n19 native geometry.",
              "m82_released": False}
    write(OUT / "m81_native_green_integral_result.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare", action="store_true")
    args = parser.parse_args()
    prepare() if args.prepare else audit()
