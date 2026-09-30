"""Read-only reduction of the frozen 816-point run; never invokes a solver."""
import argparse,csv,gzip,hashlib,json,math,shutil,statistics,tarfile
from pathlib import Path
import h5py
from simplemos_srh_cloud_20260926 import signed_si

def read(p):return json.loads(Path(p).read_text())
def rows(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def write(p,x):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def table(p,x):
    p.parent.mkdir(parents=True,exist_ok=True)
    if not x:return
    keys=list(dict.fromkeys(k for r in x for k in r))
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(x)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def collect(base,out):
    run=base/'continuation_validation_20260927-v2';matrix=run/'original_matrix';inputs=base/'original_inputs'
    out.mkdir(parents=True,exist_ok=True);provenance=[]
    def copy(p,name):
        target=out/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        provenance.append(dict(source=str(p),package_path=name,sha256=sha(p)))
    for name in ('summary.json','comparison.csv','seal.json'):
        copy(matrix/name,'frozen/'+name)
    for name in ('failure.json','regression.json','source_hashes.json'):
        if (run/name).exists():copy(run/name,'frozen/'+name)
    for name in ('contract.json','hashes.json','materials.json'):copy(inputs/name,'inputs/'+name)
    cases=read(inputs/'contract.json')['cases'];allpoints=[];stats=[];selected=set();outliers=[]
    for c in cases:
        gate=matrix/c['case']/'gate';rr=[read(p) for p in sorted(gate.glob('audit_*/result.json'))]
        assert len(rr)==51 and len({r['index'] for r in rr})==51
        for r in rr:
            r['absolute_delta_A_per_um']=r['current_A_per_um']-r['native_Id_A_per_um']
            allpoints.append({k:v for k,v in r.items() if not isinstance(v,(dict,list))})
            if not r['current_accepted']:
                outliers.append(r)
                for i in range(max(0,r['index']-2),min(50,r['index']+2)+1):selected.add((c['case'],i))
                selected.add((c['device']+'_vd_1',r['index']))
        stats.append(dict(case=c['case'],points=len(rr),qualified=sum(r['qualified'] for r in rr),within_2percent=sum(r['current_accepted'] for r in rr),max_abs_Id_error_percent=max(abs(r['Id_error_percent']) for r in rr),median_abs_Id_error_percent=statistics.median(abs(r['Id_error_percent']) for r in rr)))
        for name in ('curve.csv','config.json','attempts.csv','iterations.csv','config.status.json'):
            copy(gate/name,'sweeps/'+c['case']+'/'+name)
        copy(inputs/c['reference'],'references/'+Path(c['reference']).name)
    assert len(allpoints)==816 and len(outliers)==6
    table(out/'all_points.csv',allpoints);table(out/'case_summary.csv',stats)
    lookup={(r['case'],r['index']):r for r in allpoints}
    table(out/'neighborhoods.csv',[dict(focus_case=o['case'],focus_vg=o['vg'],**lookup[o['case'],i]) for o in outliers for i in range(max(0,o['index']-2),min(50,o['index']+2)+1)])
    table(out/'high_vd_controls.csv',[dict(focus_case=o['case'],**lookup[o['device']+'_vd_1',o['index']]) for o in outliers])
    terminals=[];fieldstats=[];toprows=[];topedges=[];state_rows=[];source_sums=[];attempts=[]
    q=1.602176634e-19
    for case,index in sorted(selected):
        c=next(c for c in cases if c['case']==case);gate=matrix/case/'gate';audit=gate/f'audit_{index:03d}'
        result=read(audit/'result.json');vg=result['vg'];label=dict(case=case,vg=vg,index=index)
        mesh=read(inputs/c['device']/'mesh.json');geo=read(inputs/c['device']/'geometry.json');free=set(geo['free_si']);areas=signed_si(mesh)
        terms=rows(audit/'all_row.csv');edges=rows(audit/'acceptance_edges.csv');copy(audit/'result.json',f'points/{case}/vg_{index:03d}/result.json')
        physical_edge=max(edges,key=lambda e:abs(float(e['electron_flux'])))
        factor=q*1e-6*float(physical_edge['electron_particle_line_flux_per_m_s'])/float(physical_edge['electron_flux'])
        for contact in mesh['contacts']:
            ids=set(contact['node_ids']);curr={}
            for car in ('electron','hole'):
                curr[car]=(-1 if car=='electron' else 1)*q*1e-6*math.fsum((int(int(e['node0']) in ids)-int(int(e['node1']) in ids))*float(e[car+'_particle_line_flux_per_m_s']) for e in edges)
            terminals.append(dict(**label,contact=contact['name'],electron_A_per_um=curr['electron'],hole_A_per_um=curr['hole'],total_A_per_um=curr['electron']+curr['hole']))
        source_sums.append(dict(**label,free_Si_SRH_A_per_um=math.fsum(float(t['electron_recombination'])*factor for t in terms if int(t['node_id']) in free),raw_residual_to_A_per_um=factor))
        for car in ('electron','hole'):
            ranked=[]
            for t in terms:
                if int(t['node_id']) not in free:continue
                scale=max(float(t[car+'_flux_abs_sum']),abs(float(t[car+'_recombination'])),abs(float(t[car+'_impact'])))
                ratio=abs(float(t[car+'_residual']))/scale if scale else 0.
                ranked.append(dict(**label,carrier=car,node_id=int(t['node_id']),x_um=float(t['x']),y_um=float(t['y']),row_ratio=ratio,residual_equivalent_A_per_um=float(t[car+'_residual'])*factor,SRH_equivalent_A_per_um=float(t[car+'_recombination'])*factor,flux_abs_sum_equivalent_A_per_um=float(t[car+'_flux_abs_sum'])*factor))
            toprows.extend(sorted(ranked,key=lambda x:x['row_ratio'],reverse=True)[:10])
        wanted=['edge_id','node0','node1','x0','y0','x1','y1','psi0_V','psi1_V','phin0_V','phin1_V','phip0_V','phip1_V','electron_mobility_m2_V_s','hole_mobility_m2_V_s','electron_sg_cancellation_condition','electron_sg_phin0_relative_V','electron_sg_phin1_relative_V','electron_qf_reference0_V','electron_qf_reference1_V']
        for e in sorted(edges,key=lambda e:abs(float(e['electron_particle_line_flux_per_m_s'])),reverse=True)[:15]:
            topedges.append(dict(**label,**{k:e[k] for k in wanted},electron_oriented_A_per_um=-q*1e-6*float(e['electron_particle_line_flux_per_m_s']),hole_oriented_A_per_um=q*1e-6*float(e['hole_particle_line_flux_per_m_s'])))
        path=gate/('state_bias_'+format(vg,'.6f').replace('.','p')+'.h5')
        with h5py.File(path) as f:
            fields={k:v[()] for k,v in f['fields'].items()};meta=json.loads(f.attrs['metadata_json'])
        assert sha(path)==read(audit/'config.status.json')['state_sha256']
        for name in ('psi','phin','phip','electrons_m3','holes_m3'):
            vals=[float(fields[name][i]) for i in free]
            fieldstats.append(dict(**label,field=name,min=min(vals),max=max(vals),median=statistics.median(vals),sample='free_si',count=len(vals)))
        if (case,index) in {(x['case'],x['index']) for x in outliers}:
            copy(path,f'detailed/{case}_vg_{index:03d}.h5')
            copy(audit/'all_row.json',f'points/{case}/vg_{index:03d}/all_row.json')
            for i,node in enumerate(mesh['nodes']):
                row=dict(**label,node_id=i,x_um=node['x'],y_um=node['y'],free_si=i in free,potential_scale_V=meta['packed_potential_scale_V'])
                row.update({k:float(v[i]) for k,v in fields.items()})
                row['SRH_integrated_A_per_um']=float(terms[i]['electron_recombination'])*factor if i in free else ''
                row['SRH_rate_per_m3_s']=float(terms[i]['electron_recombination'])*factor/(q*areas[i]*1e-6) if i in free and areas[i]>0 else ''
                state_rows.append(row)
        for attempt in rows(gate/'attempts.csv'):
            if abs(float(attempt['actual_target_bias_V'])-vg)<1e-10:attempts.append(dict(case=case,**attempt))
    for name,data in [('vela_terminals',terminals),('field_statistics',fieldstats),('top_residual_rows',toprows),('top_transport_edges',topedges),('srh_source_integrals',source_sums),('selected_attempts',attempts)]:table(out/(name+'.csv'),data)
    table(out/'detailed/outlier_node_fields.csv',state_rows)
    table(out/'outliers.csv',[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in outliers])
    for dev in sorted({r['device'] for r in outliers}):
        for name in ('mesh.json','doping.csv','geometry.json','template.json'):copy(inputs/dev/name,'inputs/'+dev+'/'+name)
    for name in ('simplemos_original_matrix_cloud_20260926.py','simplemos_continuation_cloud_20260927.py','simplemos_hfs_cloud_20260926.py'):
        copy(base/'source-continuation-v2/scripts'/name,'source/'+name)
    for name in ('src/equation/SplitDDRuntime.cpp','include/vela/numerics/SplitDDState.h','src/equation/SplitDDOperator.cpp','src/physics/RecombinationModel.cpp'):
        p=base/'source-continuation-v2'/name
        if p.exists():copy(p,'source/'+name)
    write(out/'provenance.json',provenance)
    write(out/'collection_summary.json',dict(points=816,outliers=6,selected_points=len(selected),node_rows=len(state_rows),solver_invoked=False,native_spatial_fields_available_at_outliers=False))
    write(out/'cloud_hashes.json',{p.relative_to(out).as_posix():sha(p) for p in out.rglob('*') if p.is_file() and p.name!='cloud_hashes.json'})
    archive=out.with_suffix('.tgz')
    with tarfile.open(archive,'w:gz') as t:t.add(out,arcname=out.name)
    print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,sha256=sha(archive),summary=read(out/'collection_summary.json'))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--base',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();collect(a.base,a.output)
