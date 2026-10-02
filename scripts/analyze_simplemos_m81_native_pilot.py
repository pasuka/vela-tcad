"""Audit native M81 output and release only the fixed-charge calibration."""
import csv
import json
import math
import subprocess

from run_simplemos_m81_output_amendment import REPO, ROOT, LOCAL, OUT, rows, sha, write

GREEN = ("CurPotReACGreenFunction", "CurECReACGreenFunction", "CurHCReACGreenFunction")


def verify_hashes(path):
    for name, digest in json.loads(path.read_text())["input_hashes"].items():
        if sha(REPO / name) != digest:
            raise ValueError(f"Frozen input changed: {name}")


def analyze():
    verify_hashes(ROOT / "simplemos_m81_output_amendment_freeze_v1.json")
    initial = json.loads((OUT / "m81_initial_native_result.json").read_text())
    results = []
    for base in initial["cases"]:
        device = base["device"]
        raw = LOCAL / "output_amendment_raw" / device
        for source in (LOCAL / "output_amendment_bundle" / device).iterdir():
            if sha(source) != sha(raw / source.name):
                raise ValueError("Remote amended input mismatch")
        dc = rows(raw / "m81_des.plt")
        ac = rows(raw / "m81_ac_ac_des.plt")
        if len(ac) != 1 or ac[0]["frequency"] != 0:
            raise ValueError("Expected one zero-frequency point")
        dc_error = max(abs(math.log10(r["drain TotalCurrent"] / base["reference_current_A_per_um"])) for r in dc)
        dI = ac[0]["dIuniform_charge(N_drain)"]
        response_change = abs(dI / base["native_dI_A"] - 1)
        # AC circuit current is opposite to device current. Do not transfer this
        # convention to dI: its sign is independently qualified by fixed-charge FD.
        circuit_current_error = abs(ac[0]["i(,N_drain)"] / dc[-1]["drain TotalCurrent"] + 1)
        green_files = list(raw.glob("*acgf*.tdr"))
        if len(green_files) != 1:
            raise ValueError(f"Expected one Green TDR: {green_files}")
        export = LOCAL / "green_exports" / device
        export.mkdir(parents=True, exist_ok=True)
        subprocess.run([str(REPO / "build-release/sentaurus_import.exe"), "--tdr", str(green_files[0]),
                        "--export-dir", str(export), "--inventory-json", str(export / "inventory.json")], check=True)
        manifest = json.loads((export / "field_manifest.json").read_text())["fields"]
        stats = []
        for name in GREEN:
            fields = [f for f in manifest if f["name"] == name]
            if not fields:
                raise ValueError(f"Missing field {name}")
            for field in fields:
                values = list(csv.DictReader((export / "fields" / field["csv_file"]).open()))
                if field["mapping_status"] != "complete" or len(values) != field["region_node_count"]:
                    raise ValueError("Incomplete Green mapping")
                key = "component0"
                data = [float(v[key]) for v in values]
                if not all(math.isfinite(v) for v in data):
                    raise ValueError("Non-finite Green value")
                stats.append({"field": name, "region": field["region"], "region_name": field["region_name"],
                              "unit": field["unit"], "nodes": len(data), "min": min(data), "max": max(data),
                              "max_abs": max(abs(v) for v in data)})
            if not any(r["max_abs"] > 0 for r in stats if r["field"] == name):
                raise ValueError("Identically zero Green field")
        passed = (int((raw / "exit_code.txt").read_text()) == 0 and dc_error <= 1e-5
                  and response_change <= 1e-8 and circuit_current_error <= 1e-8)
        results.append({"device": device, "dc_current_A_per_um": dc[-1]["drain TotalCurrent"],
                        "dc_reclosure_dex": dc_error, "native_dI_A": dI,
                        "response_change_relative": response_change,
                        "ac_circuit_vs_device_current_opposite_sign_error": circuit_current_error,
                        "green": stats, "pass": passed})
    report = {"status": "passed" if all(r["pass"] for r in results) else "failed",
              "cases": results, "fixed_charge_fd_released": all(r["pass"] for r in results),
              "m82_released": False,
              "raw_hashes": {str(p.relative_to(REPO)).replace("\\", "/"): sha(p)
                             for p in (LOCAL / "output_amendment_raw").rglob("*") if p.is_file()}}
    write(OUT / "m81_native_pilot_result.json", report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    analyze()
