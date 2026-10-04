"""Enormal curves after response gates, with low-Vg controls before Vela sweeps."""
import argparse,math,shutil,tarfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_enormal_candidate_20260912 as q
import validate_simplemos_enormal_same_response_20260912 as response
a,d,c,v,n=q.a,q.d,q.c,q.v,q.c.n
R=q.R;L=R/'build-release/enormal_curves_20260912';O=R/'reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912'
REMOTE='/tmp/vela_simplemos_enormal_curves_20260912'
PHUMOB=R/'reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912'

def configure():
    c.LOCAL,c.OUT,c.REMOTE=L,O,REMOTE
    n.LOCAL,n.OUT,n.REMOTE=L,O,REMOTE
    v.LOCAL,v.OUT,v.RUNNER=L,O,q.RUNNER
    v.V.post_config=c.qualification.previous.post_config;v.native_seed=c.native_seed

def gate():
    for p in (q.O/'comparison_evidence.json',q.O/'real_column_evidence.json',response.O/'comparison_evidence.json'):a.verify(p)
    assert a.read(q.O/'summary.json')['all_qualified']
    assert a.read(q.O/'real_column_summary.json')['qualified']
    assert a.read(response.O/'summary.json')['all_qualified']
    manifest=a.read(q.O/'candidate_freeze.json')
    assert a.sha(q.RUNNER)==manifest['input_hashes'][a.rel(q.RUNNER)]

def prepare_native():
    gate();assert not (O/'native_freeze.json').exists()
    cases=[];files=[Path(__file__).resolve(),q.O/'comparison_evidence.json',response.O/'comparison_evidence.json']
    for cc in a.read(q.O/'inputs.json')['cases']:
        if cc['index']!=40:continue
        src=q.p.L/'native_raw/bundle'/(cc['key']+'_baseline')
        dest=L/'bundle'/cc['case'];dest.mkdir(parents=True,exist_ok=False)
        header=(src/'native_des.cmd').read_text().split('Solve {',1)[0]
        assert 'Mobility(PhuMob Enormal)' in header and 'HighFieldSaturation' not in header and 'ExtendedPrecision(128)' in header
        text=header+'''Solve {
 Load(FilePrefix="final")
 Coupled { Poisson Electron Hole }
 Quasistationary(InitialStep=0.05 Increment=1.4 MinStep=1e-6 MaxStep=0.1 Goal { Name="gate" Voltage=0 }) { Coupled { Poisson Electron Hole } }
'''
        for i,vg in enumerate(n.GRID):
            if i:text+=f''' NewCurrentPrefix="ramp_{i:03d}_"
 Quasistationary(InitialStep=0.5 Increment=1.4 MinStep=1e-6 MaxStep=1 Goal {{ Name="gate" Voltage={vg:.12g} }}) {{ Coupled {{ Poisson Electron Hole }} }}
'''
            text+=f''' NewCurrentPrefix="point_{i:03d}_"
 Coupled {{ Poisson Electron Hole }}
 Plot(FilePrefix="vg_{i:03d}")
 Save(FilePrefix="vg_{i:03d}")
'''
        (dest/'native_des.cmd').write_text(text+'}\n',newline='\n')
        for name in ('input_fps.tdr','final_des.sav','final_circuit_des.sav'):
            shutil.copyfile(src/name,dest/name);files += [src/name,dest/name]
        files += [src/'native_des.cmd',dest/'native_des.cmd']
        cases.append(dict(case=cc['case'],device=cc['device'],vd=cc['vd'],vg=n.GRID))
    shell='''#!/bin/bash
set -u
cd "$(dirname "$0")"
printf '%s\\n' "$$" > launcher.pid
running=0
for path in bundle/*; do
 (cd "$path" || exit 91
  test ! -f exit_code.txt || exit 92
  /atctools/Synopsys/tcad/T-2022.03/bin/sdevice native_des.cmd > console.log 2>&1 &
  pid=$!; printf '%s\\n' "$pid" > solver.pid
  wait "$pid"; code=$?; printf '%s\\n' "$code" > exit_code.txt
 ) &
 running=$((running+1))
 if [ "$running" -eq 2 ]; then wait; running=0; fi
done
wait
tar czf results.tgz bundle
printf 'complete; inspect per-job exits\\n' > complete.txt
'''
    (L/'run.sh').write_text(shell,newline='\n');files.append(L/'run.sh')
    a.write(O/'native_contract.json',dict(cases=cases,remote_root=REMOTE,target_states=204,workers=2,physics='300 K PhuMob + Enormal + OldSlotboom + doping-dependent SRH; no HFS.',path='Same .8 V Enormal anchors, reclose, descend to zero and ascend in .02 V increments to one.',gates=dict(native_exit_zero=True,bias_V=1e-10,kcl_over_Id=1e-8),acceptance_changed=False))
    d.matrix.freeze(O/'native_freeze.json',files+[O/'native_contract.json'])
    with tarfile.open(L/'input.tgz','w:gz') as t:
        t.add(L/'bundle',arcname='bundle');t.add(L/'run.sh',arcname='run.sh')
    print(L/'input.tgz',a.sha(L/'input.tgz'),flush=True)

def prepare_vela():
    gate();a.verify(O/'export_evidence.json');assert not (O/'vela_freeze.json').exists()
    native=a.rows(O/'native_points.csv');assert len(native)==204 and all(r['native_qualified']=='True' for r in native)
    selected=a.rows(q.O/'selected_states.csv');cases=[]
    files=[Path(__file__).resolve(),q.O/'comparison_evidence.json',response.O/'comparison_evidence.json',O/'export_evidence.json',q.RUNNER]
    for cc in a.read(O/'native_contract.json')['cases']:
        picked=next(r for r in selected if r['case']==cc['case'] and int(r['index'])==40 and r['arm']=='phumob_seed');assert picked['qualified']=='True'
        src=Path(picked['dest']);cfg=a.read(src/'config.json');cfg.pop('state_file');cfg.pop('output_state_file')
        path=L/'vela'/cc['case']/'template.json';a.write(path,cfg)
        cases.append(dict(**cc,template=str(path),anchor=str(src/'state.csv')))
        files += [path,src/'config.json',src/'state.csv']+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    a.write(O/'vela_contract.json',dict(cases=cases,gates=v.GATES,target_states=204,paths=2,gate_step_V=.02,low_control_indices=[0,10],low_control_required_before_curves=True,acceptance_changed=False,current_metric='100*(Id_Vela/Id_Sentaurus-1), without a new absolute-current acceptance tolerance.',recovery='At most one same-bias saved-state reload. Retain failures; only qualified states seed later biases.',native_initialization='Native potentials and coherent densities using unchanged Vela ni, OldSlotboom and thermal voltage.'))
    d.matrix.freeze(O/'vela_freeze.json',files+[O/'vela_contract.json'])

def controls():
    a.verify(O/'vela_freeze.json');a.verify(PHUMOB/'completion_evidence.json')
    v.LOCAL=L/'controls';v.OUT=O/'controls'
    phumob=a.rows(PHUMOB/'continuation_attempts.csv');jobs=[]
    for cc in a.read(O/'vela_contract.json')['cases']:
        for index in (0,10):
            row=next(r for r in phumob if r['case']==cc['case'] and int(r['index'])==index and r['qualified']=='True')
            jobs += [(cc,index,'phumob_seed',Path(row['dest'])/'state.csv'),(cc,index,'native_seed',c.native_seed(cc,index))]
    with ThreadPoolExecutor(max_workers=2) as pool:attempts=[r for group in pool.map(lambda j:v.solve_target(*j),jobs) for r in group]
    v.csv_union(O/'low_control_attempts.csv',attempts)
    d.matrix.freeze(O/'low_control_evidence.json',[O/'vela_freeze.json',PHUMOB/'completion_evidence.json',O/'low_control_attempts.csv']+[f for f in (L/'controls').rglob('*') if f.is_file()])
    comparison=[];native=a.rows(O/'native_points.csv')
    for cc in a.read(O/'vela_contract.json')['cases']:
        geo,mask=v.V.m.previous.prior.support(cc)
        for index in (0,10):
            chosen=[]
            for arm in ('phumob_seed','native_seed'):
                group=[r for r in attempts if r['case']==cc['case'] and r['index']==index and r['arm']==arm]
                chosen.append(next((r for r in group if r['qualified']),group[-1]))
            delta=dict(psi_max_V=math.inf,phin_max_V=math.inf,phip_max_V=math.inf,density_max_relative=math.inf)
            if all((Path(r['dest'])/'state.csv').exists() for r in chosen):
                delta=v.V.m.previous.prior.delta_states(*[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen],mask)
            currents=[float(r.get('current_A_per_um',math.nan)) for r in chosen]
            err=abs(currents[0]/currents[1]-1) if currents[1] else math.inf
            ref=next(r for r in native if r['case']==cc['case'] and int(r['index'])==index)
            comparison.append(dict(case=cc['case'],device=cc['device'],vd=cc['vd'],vg=n.GRID[index],index=index,**delta,dual_Id_relative=err,error_percent=100*(currents[0]/float(ref['Id_A_per_um'])-1),qualified=v.dual_qualified(chosen[0]['qualified'],chosen[1]['qualified'],delta,err) and ref['native_qualified']=='True'))
    a.write_csv(O/'low_control_comparison.csv',comparison)
    summary=dict(low_points=8,low_qualified=sum(r['qualified'] for r in comparison),high_points=8,high_qualified=a.read(q.O/'summary.json')['comparison_qualified'],failed_attempts=sum(not r['qualified'] for r in attempts),all_qualified=all(r['qualified'] for r in comparison) and a.read(q.O/'summary.json')['all_qualified'],acceptance_changed=False)
    a.write(O/'control_summary.json',summary)
    d.matrix.freeze(O/'control_evidence.json',[O/'low_control_evidence.json',q.O/'comparison_evidence.json',O/'low_control_comparison.csv',O/'control_summary.json'])
    print(summary,comparison,flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=('prepare_native','unpack','export','prepare_vela','controls','continuation','native','analyze'));ap.add_argument('--sha');args=ap.parse_args();configure()
    if args.action=='unpack':
        assert args.sha and a.sha(L/'results.tgz')==args.sha;n.unpack()
    elif args.action=='export':n.export()
    elif args.action in ('continuation','native'):
        a.verify(O/'control_evidence.json');assert a.read(O/'control_summary.json')['all_qualified'];v.run(args.action)
    elif args.action=='analyze':v.analyze()
    else:globals()[args.action]()
