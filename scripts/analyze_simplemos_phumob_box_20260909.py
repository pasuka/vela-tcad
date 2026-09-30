"""Probe consistency, dual-initialization comparison and final-state Jv for PhuMob box."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import validate_simplemos_phumob_box_20260909 as v

a,d,q,OUT,LOCAL=v.a,v.d,v.q,v.OUT,v.LOCAL


def probes():
    a.verify(OUT/'fixed_raw_evidence.json')
    rows=[]
    for job in a.read(OUT/'contract.json')['fixed']:
        path=Path(job['dir'])
        mobility={int(r['edge_id']):r for r in a.rows(path/'mobility.csv')}
        for edge in a.rows(path/'edges.csv'):
            for carrier in ('electron','hole'):
                actual=float(edge[carrier+'_mobility_m2_V_s'])
                if int(edge['edge_id']) not in mobility:
                    assert actual==0.;continue
                expected=float(mobility[int(edge['edge_id'])][carrier+'_final_mobility_m2_V_s'])
                relative=abs(actual-expected)/max(abs(actual),abs(expected),1e-300)
                rows.append(dict(key=job['key'],edge=edge['edge_id'],carrier=carrier,relative=relative,qualified=relative<=1e-12))
    a.write_csv(OUT/'probe_consistency.csv',rows)
    summary=dict(checks=len(rows),failed=sum(not r['qualified'] for r in rows),max_relative=max(r['relative'] for r in rows))
    a.write(OUT/'probe_summary.json',summary)
    d.matrix.freeze(OUT/'probe_evidence.json',[OUT/'fixed_raw_evidence.json',OUT/'probe_consistency.csv',OUT/'probe_summary.json',Path(__file__).resolve()])
    print(summary,flush=True);assert summary['failed']==0


def compare():
    for arm in ('legacy','candidate'):a.verify(OUT/arm/'dc_evidence.json')
    sources={arm:a.rows(OUT/arm/'attempts.csv') for arm in ('legacy','candidate')}
    native_path=v.previous.p.OUT/'supported_export/native_points.csv'
    native=a.rows(native_path);rows=[];selection=[];activation=[]
    cases={(j['case'],j['index']):j for j in a.read(OUT/'contract.json')['jobs']}
    for (case,index),job in cases.items():
        key=f'{case}_vg_{index:03d}';geo,mask=q.run.V.m.previous.prior.support(job)
        reference=next(r for r in native if r['case']==case and r['model']=='phumob' and int(r['index'])==index)
        point=dict(key=key,case=case,device=job['device'],vd=job['vd'],vg=job['vg'],index=index,native_Id_A_per_um=float(reference['Id_A_per_um']))
        for stage in ('legacy','candidate'):
            chosen=[]
            for init in ('vela','native'):
                group=[r for r in sources[stage] if r['case']==case and int(r['index'])==index and r['arm']==init]
                selected=next((r for r in group if r['qualified']=='True'),group[-1]);chosen.append(selected)
                selection.append(dict(stage=stage,key=key,**selected))
                for carrier,value in a.read(Path(selected['dest'])/'independent_acceptance.json')['closure'].items():
                    activation.append(dict(stage=stage,key=key,initialization=init,carrier=carrier,**value))
            states=[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen]
            delta=q.run.V.m.previous.prior.delta_states(*states,mask)
            ids=[float(r['current_A_per_um']) for r in chosen];id_relative=abs(ids[0]/ids[1]-1)
            dual=q.run.v.dual_qualified(*(r['qualified']=='True' for r in chosen),delta,id_relative)
            point.update({stage+'_'+k:value for k,value in delta.items()})
            point.update({stage+'_dual_qualified':dual,stage+'_dual_Id_relative':id_relative,
                stage+'_Id_A_per_um':ids[0],stage+'_versus_native_percent':100*(ids[0]/point['native_Id_A_per_um']-1)})
        point['candidate_over_legacy_relative']=point['candidate_Id_A_per_um']/point['legacy_Id_A_per_um']-1
        rows.append(point)
    a.write_csv(OUT/'dc_comparison.csv',rows);a.write_csv(OUT/'selected_states.csv',selection);a.write_csv(OUT/'source_activation.csv',activation)
    summary=dict(points=len(rows),attempts={stage:len(rr) for stage,rr in sources.items()},
        failed_attempts={stage:sum(r['qualified']!='True' for r in rr) for stage,rr in sources.items()},
        dual_qualified={stage:sum(r[stage+'_dual_qualified'] for r in rows) for stage in sources},
        max_abs_native_error_percent={stage:max(abs(r[stage+'_versus_native_percent']) for r in rows) for stage in sources},
        source_components=len(activation),source_relative_active=sum(r['qualified'] for r in activation),
        native_G_floor_qualification=False,acceptance_changed=False)
    a.write(OUT/'dc_summary.json',summary)
    d.matrix.freeze(OUT/'dc_comparison_evidence.json',[OUT/'legacy/dc_evidence.json',OUT/'candidate/dc_evidence.json',native_path,
        OUT/'contract.json',Path(__file__).resolve()]+[OUT/n for n in ('dc_comparison.csv','selected_states.csv','source_activation.csv','dc_summary.json')])
    print(summary,rows,flush=True)


def post_jvp():
    a.verify(OUT/'dc_comparison_evidence.json')
    summary=a.read(OUT/'dc_summary.json');assert summary['dual_qualified']['candidate']==8
    jobs=[];files=[]
    for row in a.rows(OUT/'selected_states.csv'):
        if row['stage']!='candidate' or row['arm']!='vela':continue
        assert row['qualified']=='True'
        src=Path(row['dest']);cfg=a.read(src/'config.json');dest=LOCAL/'post_jvp'/row['key']
        _,mask=q.run.V.m.previous.prior.support(row)
        cfg.update(simulation_type='newton_jvp_probe',state_file=str(src/'state.csv'),output_csv=str(dest/'jvp.csv'))
        cfg.pop('output_state_file',None)
        cfg['directions']=[dict(name=f'{mode}_{step:.0e}',mode=mode,amplitude_V=step,node_ids=np.where(mask)[0][::3].tolist(),exclude_contacts=True)
            for mode in ('psi','phin','phip') for step in (1e-4,3e-5,1e-5,3e-6)]
        a.write(dest/'config.json',cfg);jobs.append(dict(key=row['key'],path=str(dest/'config.json')));files += [dest/'config.json',src/'state.csv']
    a.write(OUT/'post_contract.json',dict(jobs=jobs,gate=1e-4,weak_SRH_cross='Excluded from whole-residual FD qualification, retained in output.'))
    d.matrix.freeze(OUT/'post_freeze.json',files+[OUT/'post_contract.json',OUT/'dc_comparison_evidence.json',q.run.RUNNER,Path(__file__).resolve()])
    def one(job):
        path=Path(job['path']);status=q.run.V.execute(path,q.run.RUNNER,q.run.V.environment());assert status['exit_code']==0,status
        rows=[]
        for r in a.rows(path.parent/'jvp.csv'):
            for block in ('psi','phin','phip'):
                ana=float(r[f'analytic_{block}_norm']);fd=float(r[f'finite_difference_{block}_norm'])
                error=float(r[f'{block}_relative_error'])*max(1.,fd)/max(ana,fd,1e-300)
                weak=(r['mode'],block) in (('phin','phip'),('phip','phin'))
                rows.append(dict(key=job['key'],input=r['mode'],output=block,step_V=r['amplitude_V'],relative=error,
                    gated=not weak and float(r['amplitude_V'])<1e-4,qualified=error<=1e-4))
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for group in pool.map(one,jobs) for r in group]
    a.write_csv(OUT/'post_jvp.csv',rows);gated=[r for r in rows if r['gated']]
    summary=dict(states=len(jobs),gated_checks=len(gated),gated_failures=sum(not r['qualified'] for r in gated),max_relative=max(r['relative'] for r in gated))
    a.write(OUT/'post_summary.json',summary)
    d.matrix.freeze(OUT/'post_evidence.json',[OUT/'post_freeze.json',OUT/'post_jvp.csv',OUT/'post_summary.json']+[p for p in (LOCAL/'post_jvp').rglob('*') if p.is_file()])
    print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('probes','compare','post_jvp'))
    globals()[parser.parse_args().action]()
