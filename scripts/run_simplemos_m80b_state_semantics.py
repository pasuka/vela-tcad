"""Qualify source thermal convention, then use a coherent Vela mapped state.

No material parameter is changed. Raw Sentaurus carrier densities remain in
the identity audit; residual attribution uses densities reconstructed from
the mapped potentials with Vela's unchanged ni and thermal voltage.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import numpy as np

import run_simplemos_m80_absolute_current_attribution as original
old=original.old
REPO,ROOT=original.REPO,original.ROOT
LOCAL=REPO/"build-release/m80b_state_semantics"
OUT=ROOT/"absolute_current_state_semantics"
CONTRACT=ROOT/"simplemos_m80b_state_semantics_contract_v1.json"
FREEZE=ROOT/"simplemos_m80b_state_semantics_contract_freeze_v1.json"
EVIDENCE=ROOT/"simplemos_m80b_state_semantics_evidence.json"
DOC=REPO/"docs/validation/simplemos_m80b_state_semantics_2026-09-05.md"
SCRIPT=Path(__file__).resolve()
RAW_READER=original.m73.sentaurus_state
VT=8.617333262145e-5*300


def intrinsic():
    return next(m["ni"] for m in old.read_json(original.m73.M65_MATERIALS)["materials"] if m["name"]=="Si")*1e6


def audit_rows():
    rows=[];ni=intrinsic()
    for (device,vd,stage),export in original.m73.exports().items():
        if stage!="gate":continue
        geo=original.m73.Geometry(device);psi,n,p,spread=RAW_READER(export,geo)
        en=original.m73.scalar(export/"fields/eQuasiFermiPotential_region0.csv")
        hp=original.m73.scalar(export/"fields/hQuasiFermiPotential_region0.csv")
        ids=np.array(sorted(en));logn=np.log(n[ids]/ni);logp=np.log(p[ids]/ni)
        x=np.concatenate((logn,logp));y=np.concatenate(([psi[i]-en[i] for i in ids],[hp[i]-psi[i] for i in ids]))
        # A source-state identity inference, never a device parameter fit.
        vt=float(np.dot(x,y)/np.dot(x,x))
        rows.append({"device":device,"drain_voltage_V":vd,"source_inferred_thermal_voltage_V":vt,
                     "vela_thermal_voltage_V":VT,"thermal_relative_difference":vt/VT-1,
                     "source_identity_max_error_V":float(np.max(np.abs(y-vt*x))),
                     "raw_density_identity_max_dex":float(np.max(np.abs(x-y/VT))/math.log(10)),
                     "source_coordinate_error_um":original.m73.coordinate_error(export,geo),
                     "shared_potential_spread_V":spread})
    return rows


def record_preflight():
    original.validate()
    if original.EVIDENCE.exists():raise ValueError("M80 preflight already recorded")
    rows=audit_rows()
    old.write_csv(original.OUT/"m80_state_identity_failure.csv",rows)
    report={"status":"stopped_at_state_identity_gate","cases_checked":len(rows),"new_nonlinear_solves":0,
            "new_sentaurus_runs":0,"runner_probes":0,"m81_released":False,"m82_released":False,
            "maximum_raw_boltzmann_error_dex":max(r["raw_density_identity_max_dex"] for r in rows),
            "original_threshold_dex":old.read_json(original.CONTRACT)["thresholds"]["boltzmann_dex"],
            "maximum_source_identity_error_V":max(r["source_identity_max_error_V"] for r in rows),
            "source_inferred_thermal_voltage_V":float(np.median([r["source_inferred_thermal_voltage_V"] for r in rows])),
            "meaning":"Raw Sentaurus densities and Vela-recomputed densities are not identical. Stop M80 v1 before residual replay."}
    old.write_json(original.OUT/"m80_report.json",report)
    original.DOC.write_text("# M80 原始映射状态身份检查\n\n"
        "M80 v1 在任何残差回放前停止：16 个原生状态的坐标一致，按 Vela 常数检查原生 n/p 的 Boltzmann 关系超过冻结容差。\n\n"
        f"最大差异 {report['maximum_raw_boltzmann_error_dex']:.9g} dex；原阈值 1e-6 dex。"
        f"由电子和空穴关系共同反推的 Sentaurus 热电压为 {report['source_inferred_thermal_voltage_V']:.17g} V，"
        f"Vela 为 {VT:.17g} V；源数据关系的最大残差 {report['maximum_source_identity_error_V']:.3e} V。\n\n"
        "这是导出状态的常数约定差异证据，不能据此宣称找到了 Id–Vg 主误差。原始参数、状态及阈值保留。\n",encoding="utf-8")
    artifacts=[original.OUT/"m80_state_identity_failure.csv",original.OUT/"m80_report.json",original.DOC]
    old.write_json(original.EVIDENCE,{"status":"frozen","contract_sha256":old.sha256(original.CONTRACT),
        "artifacts":{old.portable(p):old.sha256(p) for p in artifacts},"preflight_recorder_sha256":old.sha256(SCRIPT)})
    return report


def coherent_state(export,geometry):
    psi,n,p,spread=RAW_READER(export,geometry)
    en=original.m73.scalar(export/"fields/eQuasiFermiPotential_region0.csv")
    hp=original.m73.scalar(export/"fields/hQuasiFermiPotential_region0.csv")
    ni=intrinsic()
    for i in en:n[i]=ni*math.exp((psi[i]-en[i])/VT)
    for i in hp:p[i]=ni*math.exp((hp[i]-psi[i])/VT)
    return psi,n,p,spread


def freeze():
    original.verify()
    if CONTRACT.exists() or FREEZE.exists():raise ValueError("M80b already frozen")
    c=old.read_json(original.CONTRACT)
    c.update(schema="vela.simplemos.m80b_state_semantics.v1")
    c["state_mapping"]="Map raw Sentaurus psi/phin/phip unchanged, but reconstruct n/p with Vela's unchanged ni/Vt for a coherent Vela residual state. Keep raw carrier identity differences in separate ledger; never describe raw and recomputed densities as identical."
    c["source_identity_gate"]={"common_inferred_thermal_voltage_spread_V":1e-14,
        "source_relation_error_V":1e-12,"raw_to_reconstructed_max_dex":1e-5,
        "note":"Source Vt is inferred solely to verify exported state identity; never passed into a solver or material file."}
    paths=[SCRIPT,original.SCRIPT,original.CONTRACT,original.FREEZE,original.EVIDENCE,
           original.OUT/"m80_state_identity_failure.csv",original.m73.M65_MATERIALS,
           REPO/"tests/regression/test_simplemos_m80b_state_semantics.py"]
    old.write_json(CONTRACT,c)
    old.write_json(FREEZE,{"status":"frozen_before_execution","contract_sha256":old.sha256(CONTRACT),
                           "input_hashes":{old.portable(p):old.sha256(p) for p in paths}})
    return {"status":"frozen_before_execution","physical_parameters_changed":False}


def validate():
    original.verify();f=old.read_json(FREEZE);old.m78.check_hash(CONTRACT,f["contract_sha256"])
    for rel,d in f["input_hashes"].items():old.m78.check_hash(REPO/rel,d)


def run(workers):
    validate();rows=audit_rows();gate=old.read_json(CONTRACT)["source_identity_gate"]
    values=[r["source_inferred_thermal_voltage_V"] for r in rows]
    if (max(values)-min(values)>gate["common_inferred_thermal_voltage_spread_V"] or
        max(r["source_identity_max_error_V"] for r in rows)>gate["source_relation_error_V"] or
        max(r["raw_density_identity_max_dex"] for r in rows)>gate["raw_to_reconstructed_max_dex"]):
        raise ValueError("M80b source identity failed")
    old.write_csv(OUT/"m80b_source_identity.csv",rows)
    # Reuse the frozen M80 implementation; only the explicitly declared
    # carrier representation and output namespace differ.
    original.LOCAL,original.OUT,original.DOC=LOCAL,OUT,DOC
    original.CONTRACT,original.FREEZE,original.EVIDENCE=CONTRACT,FREEZE,EVIDENCE
    original.m73.sentaurus_state=coherent_state
    return original.run(workers)


def verify():
    validate();e=old.read_json(EVIDENCE);old.m78.check_hash(CONTRACT,e["contract_sha256"])
    for rel,d in e["artifacts"].items():old.m78.check_hash(REPO/rel,d)
    return old.read_json(OUT/"m80_report.json")


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--record-preflight",action="store_true");g.add_argument("--freeze-contract",action="store_true")
    g.add_argument("--run",action="store_true");g.add_argument("--verify",action="store_true")
    p.add_argument("--workers",type=int,default=2);a=p.parse_args()
    result=record_preflight() if a.record_preflight else (freeze() if a.freeze_contract else (run(a.workers) if a.run else verify()))
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
