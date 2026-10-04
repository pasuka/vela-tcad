"""Read-only state identity ledger for M69 Load -> M81 native DC reclosure."""
import csv
import json
import math
from pathlib import Path
import subprocess

from run_simplemos_m81_output_amendment import REPO, ROOT, LOCAL, OUT, sha, write, rows


def export(tdr, target):
    target.mkdir(parents=True, exist_ok=True)
    subprocess.run([str(REPO / "build-release/sentaurus_import.exe"), "--tdr", str(tdr),
                    "--export-dir", str(target), "--inventory-json", str(target / "inventory.json")], check=True)
    manifest = json.loads((target / "field_manifest.json").read_text())["fields"]
    result = {}
    for field in manifest:
        if field["name"] not in ("ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential", "eDensity", "hDensity"):
            continue
        if field["mapping_status"] != "complete":
            raise ValueError("State mapping incomplete")
        values = list(csv.DictReader((target / "fields" / field["csv_file"]).open()))
        data = {int(r["node_id"]): float(r["component0"]) for r in values}
        if len(data) != field["region_node_count"] or not all(math.isfinite(v) for v in data.values()):
            raise ValueError("Invalid state data")
        result[(field["name"], field["region"])] = data
    return result


def audit():
    contract = json.loads((ROOT / "simplemos_m81_sentaurus_ifm_pilot_contract_v1.json").read_text())
    entries, currents, inputs = [], [], []
    for case in contract["cases"]:
        device = case["device"]
        source = REPO / case["source_state"]
        reclosed = LOCAL / "output_amendment_raw" / device / "reclosed_des.tdr"
        states = [export(path, LOCAL / "state_closure_exports" / device / name)
                  for path, name in ((source, "source"), (reclosed, "reclosed"))]
        if states[0].keys() != states[1].keys():
            raise ValueError("State field sets differ")
        for key, original in states[0].items():
            final = states[1][key]
            if original.keys() != final.keys():
                raise ValueError("State node mapping differs")
            density = key[0] in ("eDensity", "hDensity")
            deltas = [abs(math.log10(final[n] / v)) if density else abs(final[n] - v) for n, v in original.items()]
            entries.append({"device": device, "field": key[0], "region": key[1], "nodes": len(original),
                            "metric": "absolute_log10_ratio_dex" if density else "absolute_difference_V",
                            "maximum": max(deltas)})
        current_file = next(source.parent.glob("IdVg*des.plt"))
        current_rows = [r for r in rows(current_file) if abs(r["gate OuterVoltage"]-.9) < 1e-10]
        if len(current_rows) != 1:
            raise ValueError("M69 source current ambiguous")
        final_current = rows(LOCAL / "output_amendment_raw" / device / "m81_des.plt")[-1]["drain TotalCurrent"]
        original_current = current_rows[0]["drain TotalCurrent"]
        currents.append({"device": device, "m69_source_current_A_per_um": original_current,
                         "m81_reclosed_current_A_per_um": final_current,
                         "absolute_log10_ratio_dex": abs(math.log10(final_current/original_current))})
        inputs += [source, reclosed, current_file]
    write(OUT / "m81_state_closure_audit.json", {"status": "read_only_supplemental_audit",
        "note": "State differences are supplemental diagnostics; no post-hoc change to the frozen Id release gate.",
        "fields": entries, "currents": currents,
        "input_hashes": {str(p.relative_to(REPO)).replace("\\", "/"): sha(p) for p in inputs + [Path(__file__)]}})
    print(json.dumps({"currents": currents, "max_potential_difference_V": max(r["maximum"] for r in entries if r["metric"]=="absolute_difference_V"),
                      "max_density_difference_dex": max(r["maximum"] for r in entries if r["metric"]=="absolute_log10_ratio_dex")}, indent=2))


if __name__ == "__main__":
    audit()
