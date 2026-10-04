"""Explicit HFS production candidate, immutable inputs and unchanged DC gates.

Reuse completed native controls; never mutate their evidence or fitted floors.
"""
import argparse,copy,math,shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import validate_simplemos_enormal_candidate_20260912 as q

a,d,c,v,V,R=q.a,q.d,q.c,q.v,q.V,q.R
L=R/'build-release/hfs_candidate_20260914/v1'
O=R/'reference_tcad/simplemos_sentaurus2022/hfs_candidate_20260914/v1'
N=R/'reference_tcad/simplemos_sentaurus2022/hfs_restore_20260912'
NL=R/'build-release/hfs_restore_20260912'
RUNNER=L/'bin/vela_example_runner.exe'

def prepare():
    assert not (O/'inputs.json').exists()
    a.verify(N/'completion_20260914.json')
    selected=a.rows(q.O/'selected_states.csv');cases=[];files=[]
    for stage in ('pilot','rest'):
        for j in a.rows(N/f'{stage}_points.csv'):
            if j['arm']!='baseline':continue
            key=j['key'];dest=L/'inputs'/key
            old=next(r for r in selected if r['key']==key and r['arm']=='phumob_seed' and r['qualified']=='True')
            seed=Path(old['dest'])/'state.csv';cfg=a.read(Path(old['dest'])/'config.json')
            cfg['solver']['mobility'].update(model='phumob_field_lombardi',high_field_driving_force='quasi_fermi_gradient',high_field_gradient_discretization='element_vertex_partial_layer',jacobian_field_derivatives=True)
            cfg.pop('state_file');cfg.pop('output_state_file');a.write(dest/'template.json',cfg)
            cc=dict(case=j['case'],key=key,device=j['device'],vd=float(j['vd']),vg=float(j['vg']),index=int(j['index']),template=str(dest/'template.json'),enormal_seed=str(seed),native_Id_A_per_um=float(j['Id_A_per_um']))
            src=NL/f'{stage}_exports'/j['name'];geo,_=V.m.previous.prior.support(cc)
            psi,en,hp,spread=c.bgn.mapping.RAW_READER(src,geo);assert spread<=1e-12
            scalar=c.bgn.mapping.original.m73.scalar
            fn=scalar(src/'fields/eQuasiFermiPotential_region0.csv');fp=scalar(src/'fields/hQuasiFermiPotential_region0.csv')
            raw=[dict(node_id=i,psi=psi[i],phin=fn.get(i,0.),phip=fp.get(i,0.),electrons_m3=en[i],holes_m3=hp[i]) for i in range(geo.count)]
            a.write_csv(dest/'native_raw.csv',raw)
            doping={int(r['node_id']):float(r['donors_cm3'])+float(r['acceptors_cm3']) for r in a.rows(Path(cfg['node_doping_file']))}
            ni=next(m['ni'] for m in a.read(Path(cfg['materials_file']))['materials'] if m['name']=='Si')*1e6
            for i in fn:
                effective=ni*math.exp(c.bgn.delta_eg(doping[i],'old_slotboom')/(2*c.bgn.VT))
                raw[i]['electrons_m3']=effective*math.exp((psi[i]-fn[i])/c.bgn.VT);raw[i]['holes_m3']=effective*math.exp((fp[i]-psi[i])/c.bgn.VT)
            a.write_csv(dest/'native_coherent.csv',raw)
            cc.update(native_seed=str(dest/'native_coherent.csv'),raw_native=str(dest/'native_raw.csv'));cases.append(cc)
            files += [seed,Path(old['dest'])/'config.json',dest/'template.json',dest/'native_raw.csv',dest/'native_coherent.csv']+[Path(cfg[k]) for k in ('mesh_file','materials_file','node_doping_file')]
    a.write(O/'inputs.json',dict(cases=cases,gates=v.GATES,acceptance_changed=False,scope='300 K PhuMob -> alpha-zero element-distance Lombardi -> vertex Canali -> box cell/edge integration; native PartialLayer projection, any-contact-vertex electric override, 1 V/cm cutoff. Independent Enormal and coherent-native starts. Existing mathematical PhuMob floor and constants unchanged.',derivatives=dict(row_relative_gate=1e-4,synthetic_entry_gate=2e-6,weak_entries_excluded=0,cutoff='Branch-local derivative; no smoothness claim across the cutoff.')))
    d.matrix.freeze(O/'input_evidence.json',files+[O/'inputs.json',N/'completion_20260914.json',Path(__file__).resolve()])
    print('Prepared eight targets and two starts.',flush=True)

def freeze():
    a.verify(O/'input_evidence.json');assert not (O/'candidate_freeze.json').exists()
    RUNNER.parent.mkdir(parents=True,exist_ok=False);shutil.copyfile(R/'build-release/vela_example_runner.exe',RUNNER)
    files=[RUNNER,O/'input_evidence.json',Path(__file__).resolve(),R/'build-release/libvela_core.a',R/'CMakeLists.txt',R/'docs/config_schema.md']
    files += [f for folder in ('src','include','tests') for f in (R/folder).rglob('*') if f.suffix in ('.h','.cpp')]
    files += list((R/'scripts').glob('*simplemos*.py'))
    d.matrix.freeze(O/'candidate_freeze.json',files)

def config_case(cc):return copy.deepcopy(a.read(Path(cc['template'])))

def fixed():
    a.verify(O/'candidate_freeze.json');rows=[]
    expected={(r['key'],r['carrier'],int(r['cell'])):float(r['mobility']) for stage in ('pilot','rest') for r in a.rows(N/'chain_20260914'/f'{stage}_selection.csv')}
    for cc in a.read(O/'inputs.json')['cases']:
        dest=L/'fixed'/cc['key'];cfg=config_case(cc);cfg.update(simulation_type='element_box_probe',state_file=cc['raw_native'],output_csv=str(dest/'cells.csv'))
        a.write(dest/'config.json',cfg);s=V.execute(dest/'config.json',RUNNER,V.environment());assert s['exit_code']==0,s
        for r in a.rows(dest/'cells.csv'):
            if int(r['local_vertex'])!=0:continue
            for car,name in [('e','electron'),('h','hole')]:
                key=(cc['key'],car,int(r['cell_id']))
                if key not in expected:continue
                actual=float(r[name+'_mobility_m2_V_s'])*1e4;ref=expected[key];err=abs(actual/ref-1)
                rows.append(dict(key=key[0],carrier=car,cell=key[2],production_mu_cm2_V_s=actual,independent_mu_cm2_V_s=ref,relative=err,qualified=err<=1e-9))
        print(cc['key'],'fixed replay complete',flush=True)
    a.write_csv(O/'fixed_cells.csv',rows);d.matrix.freeze(O/'fixed_evidence.json',[O/'candidate_freeze.json',O/'fixed_cells.csv']+[f for f in (L/'fixed').rglob('*') if f.is_file()]);assert len(rows)==len(expected) and all(r['qualified'] for r in rows),(len(rows),len(expected),max(r['relative'] for r in rows))

def jvp():
    a.verify(O/'fixed_evidence.json');a.verify(O/'candidate_freeze.json')
    def run(cc):
        cfg=config_case(cc);count=len(a.read(Path(cfg['mesh_file']))['nodes']);dest=L/'jvp'/cc['key']
        cfg.update(simulation_type='newton_jvp_probe',state_file=cc['native_seed'],output_csv=str(dest/'jvp.csv'),row_output_csv=str(dest/'rows.csv'),
            directions=[dict(name=f'{label}_{mode}_{h:.0e}',mode=mode,amplitude_V=h,node_ids=nodes,exclude_contacts=False) for label,nodes in [('single',[1000]),('hotspots',[320,792,1000,1009,1057,1089,1091])] for mode in ('psi','phin','phip') for h in (1e-6,1e-20)],
            sample_rows=[dict(block=b,node_id=i) for b in ('psi','phin','phip') for i in range(count)])
        a.write(dest/'config.json',cfg);s=V.execute(dest/'config.json',RUNNER,V.environment());assert s['exit_code']==0,s
        groups={};failures=[]
        for r in a.rows(dest/'rows.csv'):
            an=float(r['analytic_derivative']);fd=float(r['finite_difference_derivative']);scale=max(abs(an),abs(fd));rel=abs(an-fd)/scale if scale else 0.
            key=(r['direction'],r['row_block']);g=groups.setdefault(key,dict(key=cc['key'],direction=key[0],block=key[1],rows=0,nonzero=0,max_relative=0.,qualified=True))
            good=math.isfinite(rel) and rel<=1e-4;g['rows']+=1;g['nonzero']+=scale>0;g['max_relative']=max(g['max_relative'],rel);g['qualified'] &= good
            if not good:failures.append(dict(key=cc['key'],relative=rel,**r))
        print(cc['key'],'Jv failures',len(failures),flush=True);return list(groups.values()),failures
    with ThreadPoolExecutor(max_workers=2) as pool:result=list(pool.map(run,a.read(O/'inputs.json')['cases']))
    groups=[r for g,_ in result for r in g];bad=[r for _,b in result for r in b]
    a.write_csv(O/'jvp_groups.csv',groups);a.write(O/'jvp_failures.json',bad)
    a.write(O/'jvp_summary.json',dict(groups=len(groups),nonzero_entries=sum(r['nonzero'] for r in groups),max_relative=max(r['max_relative'] for r in groups),failures=len(bad),qualified=not bad,weak_entries_excluded=0))
    d.matrix.freeze(O/'jvp_evidence.json',[O/'candidate_freeze.json',O/'jvp_groups.csv',O/'jvp_failures.json',O/'jvp_summary.json']+[f for f in (L/'jvp').rglob('*') if f.is_file()]);assert not bad,len(bad)

def dc():
    a.verify(O/'jvp_evidence.json');assert a.read(O/'jvp_summary.json')['qualified'];a.verify(O/'candidate_freeze.json')
    v.LOCAL=L;v.OUT=O;v.RUNNER=RUNNER;V.post_config=c.qualification.previous.post_config
    jobs=[(cc,arm,cc[arm]) for cc in a.read(O/'inputs.json')['cases'] for arm in ('enormal_seed','native_seed')]
    def run(job):cc,arm,seed=job;return v.solve_target(cc,cc['index'],arm,Path(seed))
    with ThreadPoolExecutor(max_workers=2) as pool:attempts=[r for result in pool.map(run,jobs) for r in result]
    v.csv_union(O/'attempts.csv',attempts);d.matrix.freeze(O/'dc_evidence.json',[O/'candidate_freeze.json',O/'attempts.csv']+[f for f in (L/'vela').rglob('*') if f.is_file()])

def summarize():
    a.verify(O/'dc_evidence.json');attempts=a.rows(O/'attempts.csv');selected=[];comparison=[]
    for cc in a.read(O/'inputs.json')['cases']:
        chosen=[]
        for arm in ('enormal_seed','native_seed'):
            group=[r for r in attempts if r['case']==cc['case'] and int(r['index'])==cc['index'] and r['arm']==arm]
            chosen.append(next((r for r in group if r['qualified']=='True'),group[-1]));selected.append(dict(key=cc['key'],**chosen[-1]))
        geo,mask=V.m.previous.prior.support(cc);delta=V.m.previous.prior.delta_states(*[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen],mask)
        currents=[float(r['current_A_per_um']) for r in chosen];dualI=abs(currents[0]/currents[1]-1)
        dual=v.dual_qualified(chosen[0]['qualified']=='True',chosen[1]['qualified']=='True',delta,dualI)
        comparison.append(dict(key=cc['key'],device=cc['device'],vd=cc['vd'],vg=cc['vg'],native_Id_A_per_um=cc['native_Id_A_per_um'],vela_Id_A_per_um=currents[0],error_percent=100*(currents[0]/cc['native_Id_A_per_um']-1),dual_Id_relative=dualI,**delta,dual_qualified=dual))
    a.write_csv(O/'comparison.csv',comparison);v.csv_union(O/'selected_states.csv',selected)
    summary=dict(states=len(selected),qualified_states=sum(r['qualified']=='True' for r in selected),attempts=len(attempts),failed_attempts=sum(r['qualified']!='True' for r in attempts),dual_qualified=sum(r['dual_qualified'] for r in comparison),max_row_ratio=max(float(r['max_row_ratio']) for r in selected),max_port_relative=max(float(r['port_relative']) for r in selected),max_abs_Id_error_percent=max(abs(r['error_percent']) for r in comparison),max_dual_Id_relative=max(r['dual_Id_relative'] for r in comparison),max_dual_phi_V=max(r[k] for r in comparison for k in ('psi_max_V','phin_max_V','phip_max_V')),max_dual_density_relative=max(r['density_max_relative'] for r in comparison),acceptance_changed=False,all_qualified=all(r['dual_qualified'] for r in comparison))
    a.write(O/'summary.json',summary);d.matrix.freeze(O/'comparison_evidence.json',[O/'dc_evidence.json',O/'comparison.csv',O/'selected_states.csv',O/'summary.json']);print(summary,comparison,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=('prepare','freeze','fixed','jvp','dc','summarize'));globals()[p.parse_args().action]()
