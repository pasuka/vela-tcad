"""Four PhuMob curves, gated by the unified sixteen-point production cohort."""
import argparse
import copy
import math
import shutil
import tarfile
from pathlib import Path
import simplemos_masetti_curve_native_20260908 as n
import simplemos_masetti_curve_vela_20260908 as v
import simplemos_phumob_supported_export_20260908 as native_anchor
import simplemos_bgn_restore_vela_20260908 as bgn
import validate_simplemos_split_production_20260912 as qualification

REPO=Path(__file__).resolve().parents[1]
LOCAL=REPO/'build-release/phumob_curves_20260912'
OUT=REPO/'reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912'
REMOTE='/tmp/vela_simplemos_phumob_curves_20260912'
ANCHORS=native_anchor.LOCAL/'native_raw/bundle/phumob'
a,d=n.a,n.d

def configure():
    n.LOCAL,n.OUT,n.REMOTE=LOCAL,OUT,REMOTE
    v.LOCAL,v.OUT=LOCAL,OUT
    v.V.post_config=qualification.previous.post_config
    v.native_seed=native_seed

def gate():
    evidence=qualification.OUT/'comparison_evidence.json'
    a.verify(evidence)
    summary=a.read(qualification.OUT/'comparison_summary.json')
    assert summary['points']==16 and summary['qualified']==16 and summary['qualified_selected']==32,summary
    assert not summary['acceptance_changed']
    a.verify(qualification.OUT/'jvp_evidence.json')
    derivatives=a.read(qualification.OUT/'jvp_summary.json')
    assert derivatives['states']==16 and derivatives['failures']==0,derivatives
    return evidence

def deck(header):
    assert 'Mobility(PhuMob)' in header and 'OldSlotboom' in header
    assert 'Enormal' not in header and 'HighFieldSaturation' not in header
    assert n.prior.MATH in header
    # The old Tcl callback exported calibration geometry at every solver point.
    # Plot retains the physical fields; curve qualification uses TDR and currents.
    header=header.replace('CurrentPlot { Tcl(tcl="source runtime.tcl") }\n','')
    text=header+'''Solve {
 Load(FilePrefix="final")
 Coupled { Poisson Electron Hole }
 Quasistationary(InitialStep=0.05 Increment=1.4 MinStep=1e-6 MaxStep=0.1 Goal { Name="gate" Voltage=0 }) {
  Coupled { Poisson Electron Hole }
 }
'''
    for i,vg in enumerate(n.GRID):
        if i:
            text+=f''' NewCurrentPrefix="ramp_{i:03d}_"
 Quasistationary(InitialStep=0.5 Increment=1.4 MinStep=1e-6 MaxStep=1 Goal {{ Name="gate" Voltage={vg:.12g} }}) {{
  Coupled {{ Poisson Electron Hole }}
 }}
'''
        text+=f''' NewCurrentPrefix="point_{i:03d}_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="vg_{i:03d}")
 Save(FilePrefix="vg_{i:03d}")
'''
    return text+'}\n'

def prepare_native():
    configure();qualified=gate();assert not (OUT/'native_freeze.json').exists()
    cases=[];files=[qualified,qualification.OUT/'jvp_evidence.json',Path(__file__).resolve()]
    reference=a.rows(native_anchor.OUT/'native_points.csv')
    for device in ('n19','n23'):
        for vd in (.05,1.):
            case=f'm65_{device}_vd_{vd:.6f}_endpoint'.replace('.','p')
            point=next(r for r in reference if r['case']==case and int(r['index'])==40 and r['model']=='phumob')
            assert point['native_qualified']=='True'
            src=ANCHORS/(case+'_vg_040');dest=LOCAL/'bundle'/case;dest.mkdir(parents=True,exist_ok=False)
            header=(src/'native_des.cmd').read_text().split('Solve {',1)[0]
            (dest/'native_des.cmd').write_text(deck(header),newline='\n')
            for name in ('input_fps.tdr','final_des.sav','final_circuit_des.sav'):
                shutil.copyfile(src/name,dest/name);files.extend([src/name,dest/name])
            files.extend([src/'native_des.cmd',dest/'native_des.cmd'])
            cases.append(dict(case=case,device=device,vd=vd,vg=n.GRID))
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
pids=()
for path in bundle/*; do
 (
  cd "$path"
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1 &
  solver_pid=$!
  printf '%s\\n' "$solver_pid" > solver.pid
  wait "$solver_pid"
  code=$?
  printf '%s\\n' "$code" > exit_code.txt
 ) &
 pids+=("$!")
 if [ "${#pids[@]}" -eq 2 ]; then
  for pid in "${pids[@]}"; do wait "$pid"; done
  pids=()
 fi
done
wait
tar -czf results.tgz bundle
printf 'complete\\n' > complete.txt
'''
    (LOCAL/'run.sh').write_text(shell,newline='\n');files.append(LOCAL/'run.sh')
    a.write(OUT/'native_contract.json',dict(cases=cases,remote_root=REMOTE,curves=4,target_states=204,
        scope='300 K plain PhuMob, OldSlotboom, doping SRH; same native .8 V anchors and physical parameters as qualified controls.',
        io_change='Omit repeated Tcl calibration geometry exports; retain native Plot fields and every saved target.',
        math=n.prior.MATH,gate_step_V=.02,workers=2,
        path='Reclose native .8 V state, descend to zero, ascend in .02 V steps to one.',
        gates=dict(native_exit_zero=True,runtime='T-2022.03-SP2',bias_error_V=1e-10,kcl_over_Id=1e-8),
        acceptance_changed=False,production_changed=False))
    files.extend([OUT/'native_contract.json',native_anchor.OUT/'native_points.csv'])
    d.matrix.freeze(OUT/'native_freeze.json',files)
    with tarfile.open(LOCAL/'input.tgz','w:gz') as tar:
        tar.add(LOCAL/'bundle',arcname='bundle');tar.add(LOCAL/'run.sh',arcname='run.sh')
    print(LOCAL/'input.tgz',flush=True)

def prepare_vela():
    configure();qualified=gate();a.verify(OUT/'native_freeze.json')
    cases=[];selected=a.rows(qualification.OUT/'selected_states.csv')
    files=[qualified,OUT/'native_freeze.json',v.RUNNER,REPO/'build-release/libvela_core.a']
    for c in a.read(OUT/'native_contract.json')['cases']:
        chosen=next(r for r in selected if r['case']==c['case'] and int(r['index'])==40 and r['arm']=='vela')
        assert chosen['qualified']=='True';root=Path(chosen['dest']);cfg=a.read(root/'config.json')
        assert cfg['solver']['split_dd_state'] and cfg['solver']['mobility']['model']=='phumob'
        cfg.pop('state_file');cfg.pop('output_state_file')
        dest=LOCAL/'vela'/c['case'];a.write(dest/'template.json',cfg)
        cases.append(dict(**c,template=str(dest/'template.json'),anchor=str(root/'state.csv')))
        files.extend([dest/'template.json',root/'config.json',root/'state.csv',root/'independent_acceptance.json'])
        files.extend(Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file'))
    files.extend((REPO/'scripts').glob('*simplemos*.py'))
    files.extend(REPO/p for p in qualification.SOURCES)
    files.extend(p for directory in ('src','include') for p in (REPO/directory).rglob('*')
                 if p.suffix in ('.cpp','.h'))
    a.write(OUT/'vela_contract.json',dict(cases=cases,gates=v.GATES,target_states=204,paths=2,
        current_metric='100*(Id_Vela/Id_Sentaurus-1); preserve signed/absolute currents and do not invent an additional percentage gate.',
        qualification='Native plus both strict Vela paths; every 1814 free carrier rows, source closure, KCL, port and dual fields/current.',
        continuation='Reclose .8 anchor, descend to zero and independently ascend to one; propagate qualified states only.',
        independent='Each native target supplies psi and both QF fields; recompute densities with unchanged Vela matched ni, OldSlotboom and Vt.',
        recovery='One same-bias reload maximum. Preserve every attempt. No change to caps, tolerances or the explicit numerical bundle.',
        coverage='Current plain PhuMob four curves. Enormal, saturation and original 16 cases remain later stages.',
        gates_changed=False,defaults_changed=False,source_changes=False))
    d.matrix.freeze(OUT/'vela_freeze.json',files+[OUT/'vela_contract.json'])

def native_seed(c,index):
    dest=LOCAL/'initial'/c['case']/f'vg_{index:03d}';path=dest/'state.csv'
    if path.exists():a.verify(dest/'freeze.json');return path
    src=LOCAL/'native_exports'/c['case']/f'vg_{index:03d}'
    geo,_=v.V.m.previous.prior.support(c)
    psi,en,hp,spread=bgn.mapping.RAW_READER(src,geo);assert spread<1e-12
    scalar=bgn.mapping.original.m73.scalar
    fn=scalar(src/'fields/eQuasiFermiPotential_region0.csv');fp=scalar(src/'fields/hQuasiFermiPotential_region0.csv')
    cfg=a.read(Path(c['template']))
    doping={int(r['node_id']):float(r['donors_cm3'])+float(r['acceptors_cm3']) for r in a.rows(Path(cfg['node_doping_file']))}
    ni=next(m['ni'] for m in a.read(Path(cfg['materials_file']))['materials'] if m['name']=='Si')*1e6
    for i in fn:
        effective=ni*math.exp(bgn.delta_eg(doping[i],'old_slotboom')/(2*bgn.VT))
        en[i]=effective*math.exp((psi[i]-fn[i])/bgn.VT);hp[i]=effective*math.exp((fp[i]-psi[i])/bgn.VT)
    a.write_csv(path,[dict(node_id=i,psi=psi[i],phin=fn.get(i,0.),phip=fp.get(i,0.),electrons_m3=en[i],holes_m3=hp[i]) for i in range(geo.count)])
    d.matrix.freeze(dest/'freeze.json',[path,Path(cfg['materials_file']),Path(cfg['node_doping_file'])]+list((src/'fields').glob('*.csv')))
    return path

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare_native','prepare_vela','unpack','export','continuation','native','analyze'));action=parser.parse_args().action
    configure()
    if action in ('prepare_native','prepare_vela'):globals()[action]()
    elif action in ('unpack','export'):getattr(n,action)()
    elif action in ('continuation','native'):v.run(action)
    else:v.analyze()
