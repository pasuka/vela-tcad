"""Unfitted native cell mobility and conservative edge/terminal calibration."""
import argparse,math,os,subprocess
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.spatial import cKDTree
import simplemos_phumob_supported_export_20260908 as s
import audit_simplemos_subset_mobility_export_20260907 as masetti

a,d=s.a,s.d
LOCAL=s.LOCAL/'analysis_20260909';OUT=s.OUT/'analysis_20260909'
MODES=('vertex_local','cell_arithmetic_carriers','cell_geometric_carriers','all_cell_arithmetic')
Q=1.602192e-19;VT=1.380662e-23*300/Q


def scalar(path,index='node_id',value='component0'):
    return {int(x[index]):float(x[value]) for x in a.rows(path)}


def load(job):
    raw=s.LOCAL/'native_raw/bundle'/job['model']/job['key']
    export=s.LOCAL/'exports'/job['model']/job['key'];prefix=job['runtime_prefix']
    manifest=a.read(export/'field_manifest.json')['fields']
    for name,unit,support in [('eDensity','cm^-3','node'),('hDensity','cm^-3','node'),('eMobility','cm^2*V^-1*s^-1','cell'),('hMobility','cm^2*V^-1*s^-1','cell')]:
        f=next(x for x in manifest if x['name']==name and x['region']==0)
        assert f['unit']==unit and f['support_kind']==support and f['mapping_status']=='complete',f
    coords={int(x['id']):(float(x['x_um']),float(x['y_um'])) for x in a.rows(export/'nodes.csv')}
    ids=sorted(coords);tree=cKDTree([coords[i] for i in ids])
    vertices=a.rows(raw/(prefix+'_vertices.csv'))
    dist,idx=tree.query([[float(x['x_um']),float(x['y_um'])] for x in vertices]);assert max(dist)<=1e-12
    mapping={int(x['index']):ids[int(j)] for x,j in zip(vertices,idx)}
    lookup={tuple(sorted(int(x['node'+str(k)]) for k in range(3))):int(x['id']) for x in a.rows(export/'elements.csv') if x['material']=='Si'}
    cells=defaultdict(list)
    for x in a.rows(raw/(prefix+'_element_vertices.csv')):
        if x['region']=='Silicon_1':cells[int(x['element'])].append(x)
    # Runtime exports include coincident auxiliary vertices not used by cells.
    # Qualify the actual silicon topology, not the unused global vertex pool.
    silicon_vertices={int(x['vertex']) for vs in cells.values() for x in vs}
    assert len({mapping[i] for i in silicon_vertices})==len(silicon_vertices)
    assert all(len(vs)==3 for vs in cells.values())
    cellmap={cid:lookup[tuple(sorted(mapping[int(x['vertex'])] for x in vs))] for cid,vs in cells.items()}
    assert len(cellmap)==len(lookup)==len(set(cellmap.values()))
    native={car:scalar(raw/(prefix+'_element_'+car+'Mobility.csv'),'index','value') for car in ('e','h')}
    plot={car:scalar(export/f'fields/{car}Mobility_region0_cells.csv','cell_id') for car in ('e','h')}
    plot_error=max(abs(native[car][cid]/plot[car][cellmap[cid]]-1) for car in ('e','h') for cid in cells)
    fields={key:scalar(export/'fields'/f'{name}_region0.csv') for key,name in (('nd','DonorConcentration'),('na','AcceptorConcentration'),('n','eDensity'),('p','hDensity'),('psi','ElectrostaticPotential'),('fn','eQuasiFermiPotential'),('fp','hQuasiFermiPotential'),('ni','EffectiveIntrinsicDensity'))}
    edges={int(x['index']):tuple(sorted((mapping[int(x['start'])],mapping[int(x['end'])]))) for x in a.rows(raw/(prefix+'_edges.csv'))}
    parts=defaultdict(list)
    for x in a.rows(raw/(prefix+'_element_edges.csv')):
        if int(x['element']) in cells:parts[edges[int(x['edge'])]].append((int(x['element']),float(x['coefficient'])))
    contacts={x['name']:set(map(int,x['node_ids'].split(';'))) for x in a.rows(export/'contacts.csv')}
    return dict(cells=cells,mapping=mapping,native=native,fields=fields,parts=dict(parts),contacts=contacts,plot_error=plot_error)


def prepare():
    a.verify(s.OUT/'export_evidence.json')
    jobs=a.rows(s.OUT/'native_points.csv');assert len(jobs)==16 and all(x['native_qualified']=='True' for x in jobs)
    inputs=[];mappings=[]
    for job in jobs:
        if job['model']!='phumob':continue
        data=load(job);f=data['fields'];node_ids={}
        def add(nd,na,n,p):
            key=str(len(inputs));inputs.append(dict(id=key,nd=nd,na=na,n=n,p=p));return key
        for node in sorted(f['n']):node_ids[node]=add(f['nd'][node],f['na'][node],f['n'][node],f['p'][node])
        for cid,vs in data['cells'].items():
            nodes=[data['mapping'][int(x['vertex'])] for x in vs]
            measures=[float(x['measure_um2']) for x in vs];volume=math.fsum(measures);assert volume>0
            w=[x/volume for x in measures]
            avg=lambda key:math.fsum(weight*f[key][node] for node,weight in zip(nodes,w))
            ng=math.exp(math.fsum(weight*math.log(f['n'][node]) for node,weight in zip(nodes,w)))
            pg=math.exp(math.fsum(weight*math.log(f['p'][node]) for node,weight in zip(nodes,w)))
            records=[(MODES[0],[node_ids[n] for n in nodes],w),
                     (MODES[1],[add(f['nd'][n],f['na'][n],avg('n'),avg('p')) for n in nodes],w),
                     (MODES[2],[add(f['nd'][n],f['na'][n],ng,pg) for n in nodes],w),
                     (MODES[3],[add(avg('nd'),avg('na'),avg('n'),avg('p'))],[1.])]
            for mode,ids,weights in records:mappings.append(dict(key=job['key'],cell=cid,mode=mode,ids=ids,weights=weights))
    a.write_csv(LOCAL/'scalar_input.csv',inputs)
    a.write(OUT/'mapping.json',dict(records=mappings))
    a.write(OUT/'contract.json',dict(jobs=jobs,modes=MODES,parameters_fitted=False,
        scalar_probe='Directly compiled current production PhuMob scalar kernel, default arsenic and 300 K.',
        geometry='Native runtime element-edge coefficients and vertex measures; validated topology mapped independently to the TDR export.',
        q_C=Q,thermal_voltage_V=VT,port_convention='Conventional total current, SG particle transport times -q*1e-4 for A/um.',
        gates=a.read(s.OUT/'native_contract.json')['gates']))
    d.matrix.freeze(OUT/'freeze.json',[Path(__file__).resolve(),s.OUT/'export_evidence.json',OUT/'contract.json',OUT/'mapping.json',LOCAL/'scalar_input.csv',s.original.LOCAL/'phumob_scalar_probe.exe',s.original.REPO/'src/physics/MobilityModel.cpp'])


def B(x):
    if abs(x)<1e-6:return 1-x/2+x*x/12
    return x/math.expm1(x)


def analyze():
    a.verify(OUT/'freeze.json')
    env=os.environ.copy();env['PATH']='D:/msys64/ucrt64/bin'+os.pathsep+env['PATH']
    done=subprocess.run([str(s.original.LOCAL/'phumob_scalar_probe.exe'),str(LOCAL/'scalar_input.csv'),str(LOCAL/'scalar_output.csv')],env=env,capture_output=True,text=True)
    a.write(OUT/'process.json',dict(exit_code=done.returncode,stdout=done.stdout,stderr=done.stderr));assert done.returncode==0
    scalar_mu={x['id']:{car:float(x['mu_'+car]) for car in ('e','h')} for x in a.rows(LOCAL/'scalar_output.csv')}
    predictions={}
    for x in a.read(OUT/'mapping.json')['records']:
        predictions[x['key'],x['cell'],x['mode']]={car:math.fsum(w*scalar_mu[i][car] for i,w in zip(x['ids'],x['weights'])) for car in ('e','h')}
    cellrows=[];summary=[];ports=[];edgechecks=[]
    before=a.rows(s.original.b.OUT/'native_points.csv')
    for job in a.read(OUT/'contract.json')['jobs']:
        data=load(job);f=data['fields'];calculated={}
        if job['model']=='masetti_control':
            for cid,vs in data['cells'].items():
                weights=[float(x['measure_um2']) for x in vs];nodes=[data['mapping'][int(x['vertex'])] for x in vs]
                calculated[cid]={car:math.fsum(w*masetti.formula(f['nd'][n]+f['na'][n],masetti.PARAMETERS[car]) for n,w in zip(nodes,weights))/math.fsum(weights) for car in ('e','h')}
        else:calculated={cid:predictions[job['key'],cid,MODES[0]] for cid in data['cells']}
        for mode in (MODES if job['model']=='phumob' else ('masetti_vertex_local',)):
            for car in ('e','h'):
                errors=[]
                for cid in data['cells']:
                    value=predictions[job['key'],cid,mode][car] if job['model']=='phumob' else calculated[cid][car]
                    native=data['native'][car][cid];error=value/native-1;errors.append(abs(error))
                    cellrows.append(dict(key=job['key'],model=job['model'],mode=mode,carrier=car,cell=cid,native_mu=native,formula_mu=value,signed_relative=error))
                summary.append(dict(key=job['key'],model=job['model'],mode=mode,carrier=car,cells=len(errors),max_relative=max(errors),p95_relative=float(np.percentile(errors,95)),runtime_plot_relative=data['plot_error'],qualified=max(errors)<=1e-7 and data['plot_error']<=1e-7))
        current={name:[] for name in data['contacts']};edge_error=0
        for (i,j),parts in data['parts'].items():
            factors={car:math.fsum(g*data['native'][car][cid] for cid,g in parts) for car in ('e','h')}
            trial={car:math.fsum(g*calculated[cid][car] for cid,g in parts) for car in ('e','h')}
            for car in ('e','h'):
                if factors[car]:edge_error=max(edge_error,abs(trial[car]/factors[car]-1))
            eta=(f['psi'][j]-f['psi'][i])/VT
            en=-VT*f['n'][i]*B(-eta-math.log(f['ni'][j]/f['ni'][i]))*math.expm1((f['fn'][i]-f['fn'][j])/VT)*factors['e']
            hp=-VT*f['p'][i]*B(eta+math.log(f['ni'][i]/f['ni'][j]))*math.expm1((f['fp'][j]-f['fp'][i])/VT)*factors['h']
            for name,ids in data['contacts'].items():current[name].append(-Q*1e-4*(int(i in ids)-int(j in ids))*(en-hp))
        currents={name:math.fsum(values) for name,values in current.items()}
        ref=float(job['Id_A_per_um']);err=currents['drain']/ref-1
        old=next(r for r in before if r['model']=='old_slotboom' and r['case']==job['case'] and int(r['index'])==int(job['index']))
        control_drift=ref/float(old['Id_A_per_um'])-1
        ports.append(dict(key=job['key'],model=job['model'],device=job['device'],vd=job['vd'],vg=job['vg'],native_Id_A_per_um=ref,replayed_Id_A_per_um=currents['drain'],signed_replay_relative=err,kcl_over_Id=abs(math.fsum(currents.values()))/abs(currents['drain']),replay_qualified=abs(err)<=1e-6,versus_Masetti_current_relative=control_drift,control_reclosure_qualified=abs(control_drift)<=1e-8 if job['model']=='masetti_control' else 'not_control'))
        edgechecks.append(dict(key=job['key'],model=job['model'],edges=len(data['parts']),max_relative=edge_error,qualified=edge_error<=1e-7))
    a.write_csv(OUT/'cell_comparison.csv',cellrows);a.write_csv(OUT/'cell_summary.csv',summary);a.write_csv(OUT/'port_replay.csv',ports);a.write_csv(OUT/'edge_comparison.csv',edgechecks)
    chosen=[r for r in summary if r['mode'] in (MODES[0],'masetti_vertex_local')]
    result=dict(native_states=len(ports),chosen_cell_checks=len(chosen),chosen_cell_failures=sum(not r['qualified'] for r in chosen),chosen_cell_max_relative=max(r['max_relative'] for r in chosen),edge_failures=sum(not r['qualified'] for r in edgechecks),port_replay_failures=sum(not r['replay_qualified'] for r in ports),port_replay_max_relative=max(abs(r['signed_replay_relative']) for r in ports),production_element_box_phumob_enabled=False)
    a.write(OUT/'summary.json',result)
    d.matrix.freeze(OUT/'evidence.json',[OUT/'freeze.json',LOCAL/'scalar_output.csv',OUT/'process.json']+[OUT/name for name in ('cell_comparison.csv','cell_summary.csv','port_replay.csv','edge_comparison.csv','summary.json')])
    print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('prepare','analyze'))
    globals()[parser.parse_args().action]()
