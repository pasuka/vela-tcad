"""Prepare/freeze a two-state native Sentaurus zero-frequency current IFM pilot.

Remote transport/execution uses explicit shell tools; this script only prepares
and audits local files. No full-curve sweep or parameter fitting is introduced.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import shutil

import run_simplemos_m80b_state_semantics as upstream
old=upstream.old
REPO,ROOT=upstream.REPO,upstream.ROOT
LOCAL=REPO/"build-release/m81_sentaurus_ifm_pilot"
BUNDLE=LOCAL/"bundle"
OUT=ROOT/"sentaurus_ifm_pilot"
CONTRACT=ROOT/"simplemos_m81_sentaurus_ifm_pilot_contract_v1.json"
FREEZE=ROOT/"simplemos_m81_sentaurus_ifm_pilot_contract_freeze_v1.json"
SCRIPT=Path(__file__).resolve()
REMOTE="/tmp/vela_simplemos_m81_ifm_20260905"


def review_m80():
    report=upstream.verify()
    if report["numerical_pass_count"]!=16 or report["linear_pass_count"]!=16:
        raise ValueError("M80 numerical/linear qualification incomplete")
    rows=old.read_csv(upstream.OUT/"m80_case_ledger.csv")
    comp=old.read_csv(upstream.OUT/"m80_component_ledger.csv")
    spatial=old.read_csv(upstream.OUT/"m80_spatial_ledger.csv")
    cases={(r["device"],float(r["drain_voltage_V"])):r for r in rows}
    terms={(r["device"],float(r["drain_voltage_V"]),r["component"]):float(r["signed_current_A_per_um"]) for r in comp}
    names=sorted({r["component"] for r in comp})
    def weight(r):return math.log10(float(r["vela_current_A_per_um"])/float(r["sentaurus_current_A_per_um"]))/float(r["actual_error_A_per_um"])
    pairs=[]
    for pair in old.read_json(upstream.CONTRACT)["matrix"]["pairs"]:
        lo,hi,vd=pair["low_device"],pair["high_device"],pair["drain_voltage_V"]
        values={name:terms[(hi,vd,name)]*weight(cases[(hi,vd)])-terms[(lo,vd,name)]*weight(cases[(lo,vd)]) for name in names}
        dominant=max(values,key=lambda name:abs(values[name]))
        stable=True
        for device in (lo,hi):
            for depth in (.05,.1,.2):
                subset=[r for r in spatial if r["device"]==device and float(r["drain_voltage_V"])==vd and
                        r["component"]==dominant and r["partition"]=="physical_regions" and float(r["depth_um"])==depth]
                winner=max(subset,key=lambda r:abs(float(r["signed_current_A_per_um"])))
                stable &= winner["region"]=="channel" and float(winner["signed_current_A_per_um"])>0
        pairs.append({**pair,"dominant_component":dominant,"dominant_pair_contribution_dex":values[dominant],
                      "channel_dominance_stable":stable,"pass":dominant=="electron_transport" and values[dominant]>0 and stable})
    if sum(r["pass"] for r in pairs)<6:raise ValueError("M80 paired spatial stability failed")
    return pairs


def prepare():
    if CONTRACT.exists() or FREEZE.exists():raise ValueError("M81 already frozen")
    pairs=review_m80();old.write_csv(OUT/"m80_pair_release_review.csv",pairs)
    manifests=old.read_json(upstream.original.m73.M69_EXPORTS)["states"]
    cases=[]
    for device in ("n19","n23"):
        row=next(r for r in manifests if r["device"]==device and float(r["drain_voltage_V"])==.05 and r["stage"]=="gate")
        root=BUNDLE/device;root.mkdir(parents=True,exist_ok=False)
        state=REPO/row["tdr"];saved=state.with_name(state.name.replace("_des.tdr","_circuit_des.sav"))
        grid=state.parent/"input_fps.tdr"
        for source,name in ((grid,"input_fps.tdr"),(state,"state_des.tdr"),(saved,"state_circuit_des.sav")):
            shutil.copy2(source,root/name)
        original=(REPO/row["deck"]).read_text(encoding="utf-8")
        prefix=original.split("Solve {",1)[0]
        case=row["case"]
        prefix=prefix.replace(case,"m81")
        prefix=prefix.replace('{ Name="drain" Voltage=0.0 }','{ Name="drain" Voltage=0.05 }')
        prefix=prefix.replace('{ Name="gate" Voltage=0.0 }','{ Name="gate" Voltage=0.9 }')
        prefix=prefix.replace("Physics { EffectiveIntrinsicDensity(NoBandGapNarrowing) }",'''Physics {
  EffectiveIntrinsicDensity(NoBandGapNarrowing)
  DeterministicVariation(
    DopingVariation "uniform_charge" (
      Conc=1e13 Type=Donor -Mobility -BandgapNarrowing
    )
  )
}''')
        prefix=prefix.replace("Math { Extrapolate", "Math { ImplicitACSystem Extrapolate")
        prefix+='''
NoisePlot {
  Potential eDensity hDensity
  CurPotReACGreenFunction CurECReACGreenFunction CurHCReACGreenFunction
}

Solve {
  Load(FilePrefix="state")
  Coupled { Poisson Electron Hole }
  Plot(FilePrefix="reclosed")
  ACCoupled(
    StartFrequency=0. EndFrequency=0. NumberOfPoints=1
    Node("source" "drain" "gate" "substrate")
    ObservationNode("drain")
    ACExtract="m81_ac" NoisePlot="m81_green"
  ) { Poisson Electron Hole Contact Circuit }
}
'''
        (root/"m81_des.cmd").write_text(prefix,encoding="utf-8",newline="\n")
        cases.append({"device":device,"drain_voltage_V":.05,"gate_voltage_V":.9,
                      "source_state":old.portable(state),"source_deck":row["deck"],
                      "deck":old.portable(root/"m81_des.cmd")})
    contract={"schema":"vela.simplemos.m81_sentaurus_ifm_pilot.v1","cases":cases,"remote_root":REMOTE,
              "physics":"M65/M69 no-BGN original mobility/SRH. ImplicitACSystem and deterministic donor variation disable only mobility/BGN variation pathways.",
              "execution":{"initial_dc_reclosures":2,"initial_ac_ifm_points":2,"new_bias_sweeps":0,
                           "explicit_fd_solves":0,"full_matrix_expansion":False},
              "gates":{"runtime":"T-2022.03-SP2","dc_reclosure_dex":1e-5,"finite_green_fields":True,
                       "native_current_variation_present":True,"m82_release":False},
              "next_step":"After successful native pilot, separately freeze explicit fixed-charge finite-difference calibration at two amplitudes; no candidate A/B or production change before that comparison.",
              "independence":"M80 mapped residual is not native residual. This pilot obtains native current Green functions and native uniform doping response.",
              "manual_pages":[152,153,783,784,809,810,815,817,1589]}
    old.write_json(CONTRACT,contract)
    inputs=[SCRIPT,upstream.SCRIPT,upstream.CONTRACT,upstream.FREEZE,upstream.EVIDENCE,
            REPO/"docs/validation/simplemos_post_m78_research_debug_plan_2026-09-05.md"]+list(BUNDLE.rglob('*'))
    inputs=[p for p in inputs if p.is_file()]
    old.write_json(FREEZE,{"status":"frozen_before_execution","contract_sha256":old.sha256(CONTRACT),
                           "input_hashes":{old.portable(p):old.sha256(p) for p in inputs}})
    return {"status":"frozen_before_execution","cases":2,"m80_stable_pairs":sum(r["pass"] for r in pairs),"remote_root":REMOTE}


def validate():
    upstream.verify();f=old.read_json(FREEZE);old.m78.check_hash(CONTRACT,f["contract_sha256"])
    for rel,d in f["input_hashes"].items():old.m78.check_hash(REPO/rel,d)
    return old.read_json(CONTRACT)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument("--prepare",action="store_true");g.add_argument("--verify-inputs",action="store_true")
    a=p.parse_args();print(json.dumps(prepare() if a.prepare else validate(),indent=2))
