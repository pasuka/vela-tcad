"""Reclose qualified D0 points with explicit finite source/drain hole exchange.

This is point reclosure from a frozen ideal-contact baseline, not a sweep from
zero. Original electrical and approved thermal gates remain unchanged.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from prepare_templates_ldmos_electrothermal import contact_boundary_lengths
from run_templates_ldmos_electrothermal_curve import state_gate
from analyze_templates_ldmos_thermal import assess


def read(p):
    return json.loads(Path(p).read_text(encoding='utf-8'))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','probe','output','contract'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--indices',nargs='+',type=int,default=[0,1,10,30])
    p.add_argument('--gates',nargs='+',type=int,choices=[4,8],default=[4,8])
    p.add_argument('--velocity-cm-s',type=float,default=1.93e6)
    p.add_argument('--seed-points',type=Path,help='Previously qualified representative-point summary for incremental reclosure')
    p.add_argument('--auger-with-generation',action='store_true',help='Explicit legacy signed Auger control; original D0 leaves this off')
    a=p.parse_args()
    if not 0<=a.velocity_cm_s<float('inf'):raise ValueError('Invalid velocity')
    if any(not 0<=i<=30 for i in a.indices) or len(set(a.indices))!=len(a.indices):raise ValueError('Invalid point indices')
    a.output.mkdir(parents=True,exist_ok=False)
    contract=read(a.contract);baseline=read(a.baseline/'full_r4_final.json')
    assert baseline['status']=='pass'
    native=read(a.baseline/'native_fields_r1/manifest.json')
    seeds=read(a.seed_points) if a.seed_points else None
    if seeds and seeds['status']!='pass':raise ValueError('Qualified seed points required')
    result=dict(scope=__doc__,status='running',velocity_cm_s=a.velocity_cm_s,points=[],
                source_sha256={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in (a.probe,a.contract,a.baseline/'full_r4_final.json',Path(__file__))})
    result['auger_with_generation']=a.auger_with_generation
    if a.seed_points:result['source_sha256'][str(a.seed_points)]=hashlib.sha256(a.seed_points.read_bytes()).hexdigest()
    def save():
        tmp=a.output/'summary.tmp';tmp.write_text(json.dumps(result,indent=2),encoding='utf-8');tmp.replace(a.output/'summary.json')
    env=dict(os.environ)
    if os.name=='nt':env['PATH']='D:/msys64/ucrt64/bin;D:/msys64/usr/bin;'+env['PATH']
    for gate in a.gates:
        ledger=read(a.baseline/f'full_vg{gate}_r4/ledger.json')
        for index in a.indices:
            point=ledger['exact_points'][index];state_path=Path(point['result'])
            if seeds:state_path=Path(next(row['result'] for row in seeds['points'] if row['gate']==gate and row['index']==index))
            state=read(state_path)
            cfg=read(state_path.parent/'input.json');mesh=read(cfg['mesh_file'])
            contact_nodes={}
            for contact in mesh['contacts']:
                if contact['name'] in ('source','drain'):
                    measures=contact_boundary_lengths(mesh,contact,cfg['coordinate_to_metres'])
                    for node,length in measures.items():
                        if node in contact_nodes:raise ValueError('Overlapping finite contacts')
                        contact_nodes[node]=length
            for b in cfg['boundaries']:
                if b['node'] in contact_nodes and b['kind']=='neutral_contact':
                    b.update(hole_recombination_velocity_m_per_s=a.velocity_cm_s*.01,boundary_length_m=contact_nodes[b['node']])
            for key in ('state_interleaved','referenced_state_interleaved','electron_qf_reference_V','hole_qf_reference_V'):
                cfg[key]=state[key]
            cfg['diagnostic_newton_max_iterations']=60
            cfg['auger_with_generation']=a.auger_with_generation
            target=a.output/f'vg{gate}_{index:02d}';target.mkdir()
            input_path=target/'input.json';input_path.write_text(json.dumps(cfg),encoding='utf-8')
            start=time.perf_counter()
            with (target/'run.log').open('w',encoding='utf-8') as f:
                process=subprocess.Popen([str(a.probe.resolve()),str(input_path.resolve()),str((target/'output.json').resolve())],env=env,stdout=f,stderr=subprocess.STDOUT)
                result['active_child']=dict(pid=process.pid,gate=gate,index=index);save()
                code=process.wait()
            result['active_child']=None
            if code:result['status']='failed';save();raise RuntimeError(f'Probe exit {code}')
            out=read(target/'output.json');numeric=state_gate(out,point['bias_V'])
            field=next(f for f in native['fields'] if f['gate_V']==gate and f['point_index']==index)
            native_path=Path(field['temperature_file'])
            assert hashlib.sha256(native_path.read_bytes()).hexdigest()==field['temperature_sha256']
            ref=read(native_path)
            thermal=assess(ref['temperature_K'],out['temperature_K'],out['nodal_area_m2'],out['lattice_source_W_per_m'],out['boundary_heat_W_per_m'],contract)
            thermal_pass=all(v['status']=='pass' for k,v in thermal['gates'].items() if not(index==0 and k=='heat_balance'))
            drain=lambda s:next(c['total_outflow_A_per_m'] for c in s['contacts'] if c['contact']=='drain')
            old_current=drain(state);current=drain(out)
            row=dict(gate=gate,index=index,bias_V=point['bias_V'],numeric=numeric,thermal=thermal,
                     pass_gate=numeric['pass_gate'] and thermal_pass,updates=out['newton_updates'],wall_seconds=time.perf_counter()-start,
                     current_A_per_m=current,baseline_current_A_per_m=old_current,
                     relative_current_change=current/old_current-1 if old_current else None,
                     max_temperature_change_K=max(abs(x-y) for x,y in zip(out['temperature_K'],state['temperature_K'])),
                     contacts=out['contacts'],baseline_contacts=state['contacts'],result=str((target/'output.json').resolve()))
            result['points'].append(row)
            for path in (state_path,state_path.parent/'input.json',input_path,target/'output.json',native_path):
                result['source_sha256'][str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
            save();print(json.dumps({k:row[k] for k in ('gate','index','pass_gate','updates','relative_current_change','max_temperature_change_K')}),flush=True)
            if not row['pass_gate']:result['status']='failed';save();raise RuntimeError('Original point gates failed')
    result['status']='pass';save()


if __name__=='__main__':
    main()
