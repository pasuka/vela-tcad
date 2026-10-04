"""Frozen full-DD adjoint ledger for the actual matched-ni cross-code gap."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import numpy as np

import run_simplemos_m79c_current_difference_calibration as upstream
import run_simplemos_m15_operator_adjoint as m15
old = upstream.old
m73 = old.m78.m73
REPO, ROOT = old.REPO, old.ROOT
LOCAL = REPO / "build-release/m80_absolute_current_attribution"
OUT = ROOT / "absolute_current_adjoint_attribution"
CONTRACT = ROOT / "simplemos_m80_absolute_current_attribution_contract_v1.json"
FREEZE = ROOT / "simplemos_m80_absolute_current_attribution_contract_freeze_v1.json"
EVIDENCE = ROOT / "simplemos_m80_absolute_current_attribution_evidence.json"
DOC = REPO / "docs/validation/simplemos_m80_absolute_current_attribution_2026-09-05.md"
SCRIPT = Path(__file__).resolve()


def freeze():
    if not upstream.verify()["m80_released"]:
        raise ValueError("M79 has not qualified; stop before M80")
    if CONTRACT.exists() or FREEZE.exists():
        raise ValueError("M80 already frozen")
    paths = {SCRIPT, upstream.SCRIPT, upstream.CONTRACT, upstream.FREEZE, upstream.EVIDENCE,
             m73.SCRIPT, Path(m15.__file__), upstream.previous.RUNNER, old.m74.CASES,
             REPO / "tests/regression/test_simplemos_m80_absolute_current_attribution.py"}
    for (device,vd,stage),root in m73.exports().items():
        if stage == "gate":
            paths.add(root / "nodes.csv")
            for pattern in ("ElectrostaticPotential_region*.csv", "eQuasiFermiPotential_region0.csv",
                            "hQuasiFermiPotential_region0.csv", "eDensity_region0.csv", "hDensity_region0.csv"):
                matches = list((root / "fields").glob(pattern))
                if not matches: raise ValueError(f"missing fields: {root}/{pattern}")
                paths.update(matches)
    contract = {
        "schema":"vela.simplemos.m80_absolute_current_attribution.v1",
        "matrix":old.read_json(old.m74.CONTRACT)["matrix"],
        "execution":{"new_nonlinear_solves":0,"new_sentaurus_runs":0,"new_physical_parameters":0},
        "state_mapping":"Use shared-node Sentaurus psi with region0 priority; Si phin/phip/n/p from same export; retain baseline qf only on insulating nodes. Never overwrite historical states.",
        "definitions":{
            "exact_identity":"gV(uV)-gS(uS) = (gV(uV)-gV(uS)) + (gV(uS)-gS(uS))",
            "prediction":"-lambda_at_uS dot (FV(uS)-FV(uV)), plus separately measured target difference; baseline adjoint prediction is a nonlinear sensitivity control.",
            "scope":"FV(uS) is Vela reaction to the mapped state, not native Sentaurus residual and not proof of a defective operator.",
            "components":"Poisson dielectric/electron/hole/dopant/boundary/reconstruction remainder; solved electron/hole transport/SRH/impact/gauge/boundary terms; remaining equation mismatch separately retained.",
            "space":"Exclusive SiSiO2/SiNitride/triple/otherSi/nonSi tags; independent shallow source/channel/drain/deep-body partition at .05/.10/.20 um. Do not double-count triple nodes.",
            "thresholds_are_screening":"Residual representation and local linear approximation qualify numerical attribution, not causal physical correction. Independent native evidence remains mandatory."
        },
        "thresholds":{"coordinate_um":1e-10,"boltzmann_dex":1e-6,"current_replay_relative":1e-8,
                      "adjoint_relative_residual":1e-10,"poisson_scale_relative":1e-6,
                      "component_reconstruction_relative":1e-5,
                      "linear_remainder_fraction":0.2,"linear_remainder_floor_dex":1e-4,
                      "stable_pairs":6,"max_cancellation_ratio_for_candidate":10},
        "release":"No candidate A/B from state substitution alone. M81 localization pilot only if qualified linear remainder and at least6/8 pair stable dominant signs with cancellation ratio <=10; otherwise stop and report unresolved representation/nonlinearity."
    }
    old.write_json(CONTRACT,contract)
    old.write_json(FREEZE,{"status":"frozen_before_execution","contract_sha256":old.sha256(CONTRACT),
                           "input_hashes":{old.portable(p):old.sha256(p) for p in sorted(paths)}})
    return {"status":"frozen_before_execution","inputs":len(paths)}


def validate():
    if not upstream.verify()["m80_released"]:raise ValueError("M79 not qualified")
    f=old.read_json(FREEZE);old.m78.check_hash(CONTRACT,f["contract_sha256"])
    for rel,digest in f["input_hashes"].items():old.m78.check_hash(REPO/rel,digest)


def blocks(rows, names):
    return np.array([[float(r[k]) for r in rows] for k in names])


def exact_split(vela, mapped, sentaurus):
    return {"actual":vela-sentaurus,"state":vela-mapped,"target":mapped-sentaurus}


def case(workflow):
    device,vd=workflow["device"],float(workflow["drain_voltage_V"])
    root=LOCAL/device/f"vd_{old.m74.voltage_tag(vd)}"
    source=REPO/workflow["phases"][-1]["source_config"]
    base=old.read_csv(source.parent/"state.csv")
    geo=m73.Geometry(device);export=m73.exports()[(device,vd,"gate")]
    coord=m73.coordinate_error(export,geo)
    psi,n,p,spread=m73.sentaurus_state(export,geo)
    en=m73.scalar(export/"fields/eQuasiFermiPotential_region0.csv")
    hp=m73.scalar(export/"fields/hQuasiFermiPotential_region0.csv")
    ni=next(m["ni"] for m in old.read_json(m73.M65_MATERIALS)["materials"] if m["name"]=="Si")*1e6
    vt=8.617333262145e-5*300
    boltzmann=max(max(abs(math.log10(n[i]/ni)-(psi[i]-en[i])/(vt*math.log(10))),
                      abs(math.log10(p[i]/ni)-(hp[i]-psi[i])/(vt*math.log(10)))) for i in en if n[i]>0 and p[i]>0)
    contract=old.read_json(CONTRACT)
    if coord>contract["thresholds"]["coordinate_um"] or boltzmann>contract["thresholds"]["boltzmann_dex"]:
        raise ValueError(f"M80 mapped-state identity failed {device} {vd}: coordinate={coord}, Boltzmann={boltzmann}")
    mapped=[]
    for row in base:
        i=int(row["node_id"])
        mapped.append({"node_id":i,"psi":psi[i],"phin":en.get(i,float(row["phin"])),
                       "phip":hp.get(i,float(row["phip"])),"electrons_m3":n[i],"holes_m3":p[i]})
    state_path=root/"sentaurus_mapped_state.csv";old.write_csv(state_path,mapped)
    statuses={};residuals={};terms={};adjoints={}
    for role,state in (("baseline",source.parent/"state.csv"),("sentaurus",state_path)):
        deck=old.base_deck(source,state,workflow);dest=root/role
        statuses[role]=old.run_probe(dict(deck,simulation_type="terminal_current_functional_probe",
            residual_output_csv=str(dest/"residual.csv")),dest/"functional.json")
        residuals[role]=old.read_csv(dest/"residual.csv")
        old.run_probe(dict(deck,simulation_type="newton_carrier_term_probe",output_csv=str(dest/"carrier.csv"),
            carrier_term_probe={"solved_equation_terms":True}),dest/"carrier.json")
        terms[role]=old.read_csv(dest/"carrier.csv")
        ast=old.run_probe(dict(deck,simulation_type="terminal_current_adjoint_probe",
            output_csv=str(dest/"adjoint.csv")),dest/"adjoint.json")
        statuses[role+"_adjoint"]=ast
        adjoints[role]=old.read_csv(dest/"adjoint.csv")
    fv=blocks(residuals["baseline"],("psi_residual","phin_residual","phip_residual"))
    fs=blocks(residuals["sentaurus"],("psi_residual","phin_residual","phip_residual"))
    delta=fs-fv
    ls=blocks(adjoints["sentaurus"],("lambda_poisson","lambda_electron","lambda_hole"))
    lv=blocks(adjoints["baseline"],("lambda_poisson","lambda_electron","lambda_hole"))
    # Calibrate the physical-to-scaled Poisson source using the same M79 known
    # parameter source, avoiding assumptions about the solver's unit scaling.
    adjroot=(upstream.previous.LOCAL if device in old.PILOT else upstream.LOCAL)/device/f"vd_{old.m74.voltage_tag(vd)}"/"baseline"
    known=blocks(old.read_csv(adjroot/"response.csv"),("source_poisson",))[0]
    bv=blocks(base,("psi","electrons_m3","holes_m3"))
    physical=m73.Q*bv[1]*(geo.volumes["barycentric_si"]-geo.volumes["all_cell"])
    free=geo.free
    active=free[np.abs(physical[free])>np.max(np.abs(physical[free]))*1e-6]
    factor=float(np.median(known[active]/physical[active]))
    scale_error=float(np.max(np.abs(known[active]/physical[active]/factor-1)))
    comp={}
    poisson={"dielectric":factor*(geo.matrices["legacy"]@(psi-bv[0])),
             "electron":factor*m73.Q*(n-bv[1])*geo.volumes["all_cell"],
             "hole":-factor*m73.Q*(p-bv[2])*geo.volumes["all_cell"],
             "dopant":np.zeros(geo.count)}
    for name,values in poisson.items():
        a=np.zeros_like(delta);a[0,free]=values[free];comp["poisson_"+name]=a
    boundary=np.zeros_like(delta);boundary[0,geo.contact_nodes]=delta[0,geo.contact_nodes]
    comp["poisson_boundary"]=boundary
    rest=np.zeros_like(delta);rest[0]=delta[0]-sum(a[0] for a in comp.values())
    comp["poisson_reconstruction_remainder"]=rest
    for name,(column,_lambda) in m15.COMPONENT_COLUMNS.items():
        a=np.zeros_like(delta);block=1 if name.startswith("electron") else 2
        a[block]=blocks(terms["sentaurus"],(column,))[0]-blocks(terms["baseline"],(column,))[0]
        comp[name]=a
    comp["continuity_reconstruction_remainder"]=delta-sum(comp.values())
    reconstruction=float(np.linalg.norm(comp["poisson_reconstruction_remainder"]+comp["continuity_reconstruction_remainder"])/max(np.linalg.norm(delta),1e-300))
    sent=next(r for r in old.read_csv(old.m74.CASES) if r["device"]==device and float(r["drain_voltage_V"])==vd)
    iv=float(sent["baseline_vela_current_A_per_um"]);is_=float(sent["sentaurus_current_A_per_um"])
    gv=statuses["baseline"]["current_A_per_um"];gs=statuses["sentaurus"]["current_A_per_um"]
    split=exact_split(gv,gs,is_)
    predicted=-float(np.sum(ls*delta));predicted_v=-float(np.sum(lv*delta))
    remainder=split["state"]-predicted
    tolerance=max(.2*abs(split["actual"]),1e-4*abs(is_)*math.log(10))
    sums={name:-float(np.sum(ls*a)) for name,a in comp.items()}
    cancellation=sum(abs(x) for x in sums.values())/max(abs(predicted),1e-300)
    ledger=[];spatial=[]
    masks,_,xy=old.m78.supports(device,geo)
    triple=masks["interface"] & masks["si_nitride_interface"]
    zones={"SiSiO2":masks["interface"] & ~triple,"SiNitride":masks["si_nitride_interface"] & ~triple,
           "triple":triple,"otherSi":masks["all_si"] & ~(masks["interface"]|masks["si_nitride_interface"]),
           "nonSi":~masks["all_si"]}
    for name,a in comp.items():
        influence=-np.sum(ls*a,axis=0)
        ledger.append({"device":device,"drain_voltage_V":vd,"component":name,
                       "signed_current_A_per_um":sums[name],"absolute_node_sum_A_per_um":float(np.sum(np.abs(influence)))})
        for zone,mask in zones.items():
            spatial.append({"device":device,"drain_voltage_V":vd,"component":name,"partition":"interface_tags",
                            "depth_um":0,"region":zone,"signed_current_A_per_um":float(np.sum(influence[mask]))})
        for depth in (.05,.1,.2):
            spatial_masks,_,_=old.m78.supports(device,geo,depth)
            for zone in ("source","channel","drain","substrate"):
                spatial.append({"device":device,"drain_voltage_V":vd,"component":name,"partition":"physical_regions",
                                "depth_um":depth,"region":zone,"signed_current_A_per_um":float(np.sum(influence[spatial_masks[zone]]))})
    th=contract["thresholds"]
    max_adj=max(statuses[role+"_adjoint"]["adjoint_relative_residual"] for role in ("baseline","sentaurus"))
    row={"device":device,"drain_voltage_V":vd,"gate_voltage_V":workflow["gate_voltage_V"],
         "vela_current_A_per_um":gv,"sentaurus_current_A_per_um":is_,"mapped_current_A_per_um":gs,
         "actual_error_A_per_um":split["actual"],"state_difference_A_per_um":split["state"],"target_difference_A_per_um":split["target"],
         "predicted_state_difference_A_per_um":predicted,"baseline_adjoint_prediction_A_per_um":predicted_v,
         "linear_remainder_A_per_um":remainder,"linear_remainder_fraction":abs(remainder)/max(abs(split["actual"]),1e-300),
         "component_cancellation_ratio":cancellation,"component_reconstruction_relative":reconstruction,
         "poisson_scale_relative_error":scale_error,"boltzmann_identity_dex":boltzmann,
         "current_replay_relative":abs(gv/iv-1),"adjoint_relative_residual":max_adj,
         "representation_pass":scale_error<=th["poisson_scale_relative"] and reconstruction<=th["component_reconstruction_relative"],
         "linear_pass":abs(remainder)<=tolerance,
         "candidate_cancellation_pass":cancellation<=th["max_cancellation_ratio_for_candidate"]}
    row["numerical_pass"]=row["representation_pass"] and row["current_replay_relative"]<=th["current_replay_relative"] and max_adj<=th["adjoint_relative_residual"]
    old.write_csv(root/"components.csv",ledger);old.write_csv(root/"spatial.csv",spatial)
    old.write_json(root/"case_summary.json",row)
    print(f"M80 {device} Vd={vd}: numerical={row['numerical_pass']} linear={row['linear_pass']} remainder={row['linear_remainder_fraction']:.3g}",flush=True)
    return row,ledger,spatial


def run(workers):
    validate()
    if (OUT/"m80_report.json").exists():raise ValueError("M80 results exist")
    old.m74.RUNNER=upstream.previous.RUNNER
    workflows=old.read_json(old.m74.PORTABLE_MANIFEST)["workflows"]
    rows=[];components=[];spatial=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for r,c,s in pool.map(case,workflows):rows.append(r);components.extend(c);spatial.extend(s)
    old.write_csv(OUT/"m80_case_ledger.csv",rows);old.write_csv(OUT/"m80_component_ledger.csv",components);old.write_csv(OUT/"m80_spatial_ledger.csv",spatial)
    report={"status":"attribution_requires_review","cases":len(rows),"new_nonlinear_solves":0,"new_sentaurus_runs":0,
            "numerical_pass_count":sum(r["numerical_pass"] for r in rows),"linear_pass_count":sum(r["linear_pass"] for r in rows),
            "candidate_cancellation_pass_count":sum(r["candidate_cancellation_pass"] for r in rows),
            "maximum_linear_remainder_fraction":max(r["linear_remainder_fraction"] for r in rows),
            "m81_released":False,"m82_released":False,
            "scope":"Mapped-state reaction; not native Sentaurus residual. Pair/spatial stability must be reviewed before any release."}
    old.write_json(OUT/"m80_report.json",report)
    lines=["# M80 当前绝对 Id 差异的完整 DD 伴随账本","",f"已完成 {len(rows)} 工况，无新非线性求解。",
           f"数值账本通过 {report['numerical_pass_count']}/16，局部线性余项通过 {report['linear_pass_count']}/16。",
           "M81/M82 尚未放行；映射状态残差不是 Sentaurus 原生残差，不能直接宣称物理根因。","",
           "| 器件 | Vd | 实际误差 A/um | 状态项 | 端口目标项 | 线性余项/实际误差 | 分项抵消比 |",
           "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:lines.append(f"| {r['device']} | {r['drain_voltage_V']} | {r['actual_error_A_per_um']:.4e} | {r['state_difference_A_per_um']:.4e} | {r['target_difference_A_per_um']:.4e} | {r['linear_remainder_fraction']:.3g} | {r['component_cancellation_ratio']:.3g} |")
    DOC.write_text("\n".join(lines)+"\n",encoding="utf-8")
    artifacts=list(OUT.glob('*.csv'))+[OUT/"m80_report.json",DOC]+list(LOCAL.rglob('*.json'))+list(LOCAL.rglob('*.csv'))
    old.write_json(EVIDENCE,{"status":"frozen","contract_sha256":old.sha256(CONTRACT),"artifacts":{old.portable(p):old.sha256(p) for p in sorted(artifacts)}})
    return report


def verify():
    validate();e=old.read_json(EVIDENCE);old.m78.check_hash(CONTRACT,e["contract_sha256"])
    for rel,digest in e["artifacts"].items():old.m78.check_hash(REPO/rel,digest)
    return old.read_json(OUT/"m80_report.json")


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--freeze-contract",action="store_true");g.add_argument("--run",action="store_true");g.add_argument("--verify",action="store_true")
    p.add_argument("--workers",type=int,default=2);args=p.parse_args()
    print(json.dumps(freeze() if args.freeze_contract else (run(args.workers) if args.run else verify()),indent=2))


if __name__=="__main__":main()
