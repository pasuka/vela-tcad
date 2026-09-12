#!/usr/bin/env python3
"""Apply the approved D0 thermal gates to aligned, per-unit-width SI evidence."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path


def assess(reference_K: list[float], candidate_K: list[float], volumes: list[float],
           source_W_per_m: float, outward_W_per_m: float, contract: dict) -> dict:
    if not reference_K or len(reference_K) != len(candidate_K) or len(reference_K) != len(volumes):
        raise ValueError("Temperature and volume vectors must be nonempty and exactly aligned")
    if not all(math.isfinite(x) and x > 0 for x in reference_K + candidate_K + volumes):
        raise ValueError("Temperatures and control volumes must be finite and positive")
    if not all(math.isfinite(x) for x in (source_W_per_m, outward_W_per_m)):
        raise ValueError("Heat integrals must be finite")
    ambient = contract["ambient_temperature_K"]
    weight = math.fsum(volumes)
    rms_error = math.sqrt(math.fsum(v * (c-r)**2 for r,c,v in zip(reference_K,candidate_K,volumes))/weight)
    reference_rms = math.sqrt(math.fsum(v * (r-ambient)**2 for r,v in zip(reference_K,volumes))/weight)
    peak_error = abs(max(candidate_K)-max(reference_K))
    peak = contract["peak_temperature_rise"]
    field = contract["temperature_rise_field"]
    peak_limit = max(peak["absolute_error_floor_K"],peak["relative_error_max"]*(max(reference_K)-ambient))
    rms_limit = max(field["absolute_rms_error_floor_K"],field["relative_rms_error_max"]*reference_rms)
    scale = max(abs(source_W_per_m),abs(outward_W_per_m))
    balance = abs(source_W_per_m-outward_W_per_m)/scale if scale else None
    gates = {
        "heat_balance": {"observed":balance,"limit":contract["heat_balance"]["relative_error_max"],
                         "status":"not_evaluable" if balance is None else "pass" if balance<=contract["heat_balance"]["relative_error_max"] else "fail"},
        "peak_temperature_rise": {"error_K":peak_error,"limit_K":peak_limit,"status":"pass" if peak_error<=peak_limit else "fail"},
        "temperature_rise_field": {"rms_error_K":rms_error,"limit_K":rms_limit,"reference_rms_rise_K":reference_rms,"status":"pass" if rms_error<=rms_limit else "fail"},
    }
    statuses={x["status"] for x in gates.values()}
    return {"schema":"vela.templates_ldmos.thermal_gate_results.v1",
            "status":"fail" if "fail" in statuses else "not_evaluable" if "not_evaluable" in statuses else "pass",
            "scope":"Thermal gates only; electrical qualification must be supplied separately",
            "nodes":len(volumes),"gates":gates,
            "reference_peak_K":max(reference_K),"candidate_peak_K":max(candidate_K),
            "absolute_heat_imbalance_W_per_m":abs(source_W_per_m-outward_W_per_m)}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference",type=Path,required=True,help="JSON with node_id and temperature_K arrays")
    parser.add_argument("--candidate",type=Path,required=True,help="Thermal probe JSON; requires node_id array")
    parser.add_argument("--contract",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    reference=json.loads(args.reference.read_text());candidate=json.loads(args.candidate.read_text());contract=json.loads(args.contract.read_text())
    ids=reference["node_id"]
    if ids!=candidate["node_id"] or len(set(ids))!=len(ids) or len(ids)!=len(reference["temperature_K"]):
        raise ValueError("Node mapping must be identical and unique; no interpolation or truncation")
    result=assess(reference["temperature_K"],candidate["temperature_K"],candidate["nodal_area_m2"],
                  candidate["integrated_source_W_per_m"],candidate["outward_boundary_heat_W_per_m"],contract)
    result["provenance_sha256"]={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.reference,args.candidate,args.contract)}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(result["status"])
    raise SystemExit(0 if result["status"]=="pass" else 2)


if __name__=="__main__":
    main()
