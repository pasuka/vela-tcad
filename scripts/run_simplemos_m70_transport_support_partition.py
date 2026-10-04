"""Freeze and run the conditional read-only SimpleMOS M70 transport support partition."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

REPO=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(REPO/"scripts"))
import run_simplemos_m47_default_bgn_self_consistent_attribution as m47  # noqa: E402
import run_simplemos_m65_nobgn_intrinsic_density_attribution as m65  # noqa: E402
ROOT=REPO/"reference_tcad/simplemos_sentaurus2022"; CONTRACT=ROOT/"simplemos_m70_transport_support_partition_contract_v1.json"; FREEZE=ROOT/"simplemos_m70_transport_support_partition_contract_freeze.json"
M67_EVIDENCE=ROOT/"simplemos_m67_residual_barrier_partition_evidence.json"; M67_REPORT=ROOT/"residual_barrier_partition/m67_residual_barrier_partition_report.json"; M69_EVIDENCE=ROOT/"simplemos_m69_stage_residual_localization_evidence.json"
M65_VELA=REPO/"build-release/m65_ni/vela_manifest.json"; M65_SENT=REPO/"build-release/m65_ni/sentaurus_export_manifest.json"; M65_PAIRS=ROOT/"nobgn_intrinsic_density_attribution/m65_nwell_pair_closure_ledger.csv"
PORTABLE=ROOT/"transport_support_partition"; REPORT=PORTABLE/"m70_transport_support_partition_report.json"; STATES=PORTABLE/"m70_transport_support_state_ledger.csv"; PAIRS=PORTABLE/"m70_transport_support_pair_ledger.csv"; NODES=PORTABLE/"m70_transport_support_node_ledger.csv"
DOC=REPO/"docs/validation/simplemos_m70_transport_support_partition_2026-09-03.md"; ARTIFACT=REPO/"docs/validation/reports/simplemos_m70/artifact.json"; EVIDENCE=ROOT/"simplemos_m70_transport_support_partition_evidence.json"; SCRIPT=Path(__file__).resolve(); Q=1.602176634e-19

def read_json(path:Path)->dict[str,Any]: return json.loads(path.read_text(encoding="utf-8-sig"))
def write_json(path:Path,value:Any)->None: path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8",newline="\n")
def read_csv(path:Path)->list[dict[str,str]]:
    with path.open(newline="",encoding="utf-8-sig") as stream:return list(csv.DictReader(stream))
def write_csv(path:Path,rows:list[dict[str,Any]])->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as stream: writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
def sha256(path:Path)->str:return hashlib.sha256(path.read_bytes()).hexdigest()
def portable(path:Path)->str:return path.resolve().relative_to(REPO.resolve()).as_posix()
def median(values:list[float])->float:
    values=sorted(values);middle=len(values)//2;return values[middle] if len(values)%2 else .5*(values[middle-1]+values[middle])
def vector_csv(path:Path)->dict[int,tuple[float,float]]:return {int(row["node_id"]):(float(row["component0"]),float(row["component1"])) for row in read_csv(path)}


def vtk_vector(path:Path,name:str,count:int)->list[tuple[float,float]]:
    lines=path.read_text(encoding="utf-8",errors="replace").splitlines();anchor=f"VECTORS {name} double"
    try:index=lines.index(anchor)+1
    except ValueError as exc:raise RuntimeError(f"missing VTK vector {name}: {path}") from exc
    values=[]
    for line in lines[index:index+count]:
        parts=line.split()
        if len(parts)<2:raise RuntimeError(f"short VTK vector {name}")
        values.append((float(parts[0]),float(parts[1])))
    if len(values)!=count:raise RuntimeError(f"VTK vector count {name}")
    return values


def sources()->list[Path]:
    paths=[M67_EVIDENCE,M67_REPORT,M69_EVIDENCE,M65_VELA,M65_SENT,M65_PAIRS]
    for row in read_json(M65_VELA)["workflows"]:paths.extend([REPO/row["state"],REPO/row["vtk"]])
    for row in read_json(M65_SENT)["states"]:
        root=REPO/row["export_dir"];paths.extend([root/"nodes.csv",root/"elements.csv",root/"fields/eDensity_region0.csv",root/"fields/eMobility_region0.csv",root/"fields/eGradQuasiFermi_region0.csv",root/"fields/eCurrentDensity_region0.csv"])
    return paths


def freeze_contract()->None:
    contract=read_json(CONTRACT);m67=read_json(M67_REPORT);m69=read_json(M69_EVIDENCE)
    if contract.get("schema")!="vela.simplemos.sdevice.m70_transport_support_partition_contract.v1":raise ValueError("unexpected M70 contract")
    if read_json(M67_EVIDENCE).get("status")!="frozen" or float(m67["summary"]["median_pair_closure_fraction"])>=float(contract["activation"]["required_m67_median_pair_closure_below"]):raise ValueError("M70 activation rule not met by M67")
    if m69.get("status")!=contract["activation"]["required_m69_status"]:raise ValueError("M69 is not frozen")
    paths=sources();missing=[path for path in paths if not path.is_file()]
    if missing:raise FileNotFoundError(missing[0])
    write_json(FREEZE,{"schema":"vela.simplemos.sdevice.m70_transport_support_partition_contract_freeze.v1","status":"frozen_before_execution","contract":portable(CONTRACT),"contract_sha256":sha256(CONTRACT),"activation":{"m67_median_pair_closure_fraction":m67["summary"]["median_pair_closure_fraction"],"m69_classification":m69["classification"]},"upstream_hashes":{portable(path):sha256(path) for path in paths}})


def validate_contract()->dict[str,Any]:
    contract,freeze=read_json(CONTRACT),read_json(FREEZE)
    if freeze.get("status")!="frozen_before_execution" or freeze.get("contract_sha256")!=sha256(CONTRACT):raise ValueError("M70 contract not frozen")
    for relative,expected in freeze["upstream_hashes"].items():
        if sha256(REPO/relative)!=expected:raise ValueError(f"M70 input changed: {relative}")
    return contract


def analyze(contract:dict[str,Any])->tuple[dict[str,Any],dict[str,Any]]:
    vela_manifest,sent_manifest=read_json(M65_VELA),read_json(M65_SENT);vela_by={(row["device"],float(row["drain_voltage_V"])):row for row in vela_manifest["workflows"]};state_rows=[];node_rows=[]
    for item in sent_manifest["states"]:
        key=(item["device"],float(item["drain_voltage_V"]));workflow=vela_by[key];export=REPO/item["export_dir"];state=m65.vela_state(REPO/workflow["state"]);coords=m47.coordinates(export);_,support=m47.interface_nodes(export);sent_n=m47.scalar_field(export,"eDensity");sent_mu=m47.scalar_field(export,"eMobility");sent_grad=vector_csv(export/"fields/eGradQuasiFermi_region0.csv");sent_j=vector_csv(export/"fields/eCurrentDensity_region0.csv")
        count=max(state["psi"])+1;vela_mu=m47.vtk_scalar(REPO/workflow["vtk"],"ElectronMobilityCm2PerVs",count);vela_grad=vtk_vector(REPO/workflow["vtk"],"ElectronGradQuasiFermiVector",count)
        common=sorted(set(support)&set(sent_n)&set(state["psi"]));jmag={node:math.hypot(*sent_j[node]) for node in common};threshold=max(jmag.values())*1e-3
        active=[node for node in common if jmag[node]>=threshold and sent_n[node]>0 and sent_mu[node]>0 and math.hypot(*sent_grad[node])>0 and state["electrons_m3"][node]>0 and vela_mu[node]>0 and math.hypot(*vela_grad[node])>0]
        components={"density":[],"mobility":[],"gradient":[],"product":[]};identities=[]
        for node in active:
            dn=math.log10((state["electrons_m3"][node]/1e6)/sent_n[node]);dmu=math.log10(vela_mu[node]/sent_mu[node]);dg=math.log10((math.hypot(*vela_grad[node])/100.0)/math.hypot(*sent_grad[node]));dp=math.log10((Q*(state["electrons_m3"][node]/1e6)*vela_mu[node]*(math.hypot(*vela_grad[node])/100.0))/(Q*sent_n[node]*sent_mu[node]*math.hypot(*sent_grad[node])));identity=dp-dn-dmu-dg
            components["density"].append(dn);components["mobility"].append(dmu);components["gradient"].append(dg);components["product"].append(dp);identities.append(identity);node_rows.append({"device":item["device"],"drain_voltage_V":item["drain_voltage_V"],"gate_voltage_V":item["gate_voltage_V"],"node_id":node,"x_um":coords[node][0],"y_um":coords[node][1],"sentaurus_current_density_A_per_cm2":jmag[node],"density_log_ratio_dex":dn,"mobility_log_ratio_dex":dmu,"qf_gradient_log_ratio_dex":dg,"transport_product_log_ratio_dex":dp,"component_identity_residual_dex":identity})
        state_rows.append({"device":item["device"],"drain_voltage_V":item["drain_voltage_V"],"gate_voltage_V":item["gate_voltage_V"],"active_node_count":len(active),"median_density_log_ratio_dex":median(components["density"]),"median_mobility_log_ratio_dex":median(components["mobility"]),"median_qf_gradient_log_ratio_dex":median(components["gradient"]),"median_transport_product_log_ratio_dex":median(components["product"]),"maximum_component_identity_residual_dex":max(abs(value) for value in identities)})
    states={(row["device"],float(row["drain_voltage_V"])):row for row in state_rows};pair_rows=[];closures=[]
    for pair in read_csv(M65_PAIRS):
        low,high,drain=pair["low_device"],pair["high_device"],float(pair["drain_voltage_V"]);lo,hi=states[(low,drain)],states[(high,drain)];current=float(pair["matched_ni_pair_growth_dex"])
        growth={name:float(hi[f"median_{name}_log_ratio_dex"])-float(lo[f"median_{name}_log_ratio_dex"]) for name in ("density","mobility","qf_gradient","transport_product")};raw=1-abs(growth["transport_product"]-current)/max(abs(current),1e-300);closure=min(1.0,max(0.0,raw));closures.append(closure)
        pair_rows.append({"low_device":low,"high_device":high,"drain_voltage_V":drain,"diagnostic_gate_voltage_V":pair["diagnostic_gate_voltage_V"],"current_pair_growth_dex":current,"density_support_pair_growth_dex":growth["density"],"mobility_support_pair_growth_dex":growth["mobility"],"qf_gradient_support_pair_growth_dex":growth["qf_gradient"],"transport_product_pair_growth_dex":growth["transport_product"],"component_sum_identity_residual_dex":growth["transport_product"]-growth["density"]-growth["mobility"]-growth["qf_gradient"],"proxy_minus_current_growth_dex":growth["transport_product"]-current,"raw_pair_closure_fraction":raw,"classification_pair_closure_fraction":closure})
    metric=median(closures);analysis=contract["analysis"]
    classification="transport_support_proxy_dominant" if metric>=float(analysis["dominant_minimum_median_pair_closure_fraction"]) else "transport_support_proxy_material_but_not_complete" if metric>=float(analysis["material_minimum_median_pair_closure_fraction"]) else "transport_support_proxy_not_material"
    checks={"contract_frozen":read_json(FREEZE)["contract_sha256"]==sha256(CONTRACT),"state_count":len(state_rows)==int(contract["acceptance"]["required_state_count"]),"pair_count":len(pair_rows)==int(contract["acceptance"]["required_pair_count"]),"active_node_count":min(int(row["active_node_count"]) for row in state_rows)>=int(contract["acceptance"]["minimum_active_node_count_per_state"]),"component_identity":max(float(row["maximum_component_identity_residual_dex"]) for row in state_rows)<=float(contract["acceptance"]["maximum_component_identity_residual_dex"]),"classification_declared":classification in analysis["classifications"],"production_reference_not_replaced":True};checks["all_checks_pass"]=all(checks.values())
    if not checks["all_checks_pass"]:classification="execution_or_identity_failure"
    report={"schema":"vela.simplemos.sdevice.m70_transport_support_partition_report.v1","status":"accepted" if checks["all_checks_pass"] else "failed","classification":classification,"contract":{"path":portable(CONTRACT),"sha256":sha256(CONTRACT)},"execution":{"new_sentaurus_solves":0,"new_vela_solves":0,"production_default_changed":False},"summary":{"median_pair_closure_fraction":metric,"minimum_pair_closure_fraction":min(closures),"maximum_pair_closure_fraction":max(closures),"median_current_pair_growth_dex":median([float(row["current_pair_growth_dex"]) for row in pair_rows]),"median_density_support_pair_growth_dex":median([float(row["density_support_pair_growth_dex"]) for row in pair_rows]),"median_mobility_support_pair_growth_dex":median([float(row["mobility_support_pair_growth_dex"]) for row in pair_rows]),"median_qf_gradient_support_pair_growth_dex":median([float(row["qf_gradient_support_pair_growth_dex"]) for row in pair_rows]),"median_transport_product_pair_growth_dex":median([float(row["transport_product_pair_growth_dex"]) for row in pair_rows]),"minimum_active_node_count":min(int(row["active_node_count"]) for row in state_rows)},"causal_scope":{"closed":"The native state support is partitioned into density, mobility, and qf-gradient factors.","not_closed":"This diagnostic does not alter or revalidate SG/contact extraction; a remaining mismatch is spatial integration/discretization support, not authorization to reopen closed topics."},"acceptance":checks}
    return report,{"states":state_rows,"pairs":pair_rows,"nodes":node_rows}


def freeze_results(report:dict[str,Any],rows:dict[str,Any])->None:
    write_csv(STATES,rows["states"]);write_csv(PAIRS,rows["pairs"]);write_csv(NODES,rows["nodes"]);write_json(REPORT,report);s=report["summary"]
    DOC.parent.mkdir(parents=True,exist_ok=True);DOC.write_text(f"""# SimpleMOS M70 输运支撑分解

M70 分类为 `{report['classification']}`。在源侧有效电流支撑上，8 个 NWell 配对的中位电流增幅为 `{float(s['median_current_pair_growth_dex']):.6f} dex`；`n*mu*|grad(phin)|` 代理为 `{float(s['median_transport_product_pair_growth_dex']):.6f} dex`，中位闭合率 `{float(s['median_pair_closure_fraction']):.2%}`。代理分量中，载流子浓度 `{float(s['median_density_support_pair_growth_dex']):.6f} dex`、迁移率 `{float(s['median_mobility_support_pair_growth_dex']):.6f} dex`、准费米梯度 `{float(s['median_qf_gradient_support_pair_growth_dex']):.6f} dex`。

这是现有状态的输运支撑审计，不是 SG 或接触电流提取的重新排查。机器报告：`{portable(REPORT)}`。
""",encoding="utf-8",newline="\n")
    write_json(ARTIFACT,{"schema":"vela.validation.artifact.v1","title":"SimpleMOS M70 transport support partition","status":report["status"],"classification":report["classification"],"report":portable(REPORT),"ledgers":[portable(STATES),portable(PAIRS),portable(NODES)]});artifacts=[REPORT,STATES,PAIRS,NODES,DOC,ARTIFACT]
    write_json(EVIDENCE,{"schema":"vela.simplemos.sdevice.m70_transport_support_partition_evidence.v1","status":"frozen" if report["acceptance"]["all_checks_pass"] else "failed","classification":report["classification"],"contract_sha256":sha256(CONTRACT),"implementation_hashes":{portable(SCRIPT):sha256(SCRIPT)},"artifacts":{portable(path):sha256(path) for path in artifacts},"new_sentaurus_execution":False,"new_vela_execution":False,"production_reference_replaced":False,"acceptance":report["acceptance"]})


def verify()->dict[str,Any]:
    validate_contract();e,r=read_json(EVIDENCE),read_json(REPORT)
    if e["status"]!="frozen" or r["status"]!="accepted":raise ValueError("M70 not frozen")
    if e["implementation_hashes"][portable(SCRIPT)]!=sha256(SCRIPT):raise ValueError("M70 implementation changed")
    for relative,expected in e["artifacts"].items():
        if sha256(REPO/relative)!=expected:raise ValueError(f"M70 artifact changed: {relative}")
    return r


def main()->None:
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--freeze-contract",action="store_true");parser.add_argument("--analyze",action="store_true");parser.add_argument("--verify",action="store_true");args=parser.parse_args()
    if args.freeze_contract:freeze_contract();print(json.dumps({"status":"frozen_before_execution","contract_sha256":sha256(CONTRACT)}));return
    contract=validate_contract()
    if args.verify:print(json.dumps(verify()["summary"],indent=2));return
    if args.analyze:
        report,rows=analyze(contract);freeze_results(report,rows);print(json.dumps({"status":report["status"],"classification":report["classification"],"summary":report["summary"],"acceptance":report["acceptance"]},indent=2))


if __name__=="__main__":main()
