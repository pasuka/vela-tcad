"""Qualify the unchanged-physics DC comparison and derivatives at final states."""
import argparse,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import validate_simplemos_phumob_chain_fix_20260909 as run

a,d,q=run.a,run.d,run.q;LOCAL,OUT=run.LOCAL,run.OUT


def compare():
    for arm in ('before','after'):a.verify(OUT/arm/'dc_evidence.json')
    sources={arm:a.rows(OUT/arm/'attempts.csv') for arm in ('before','after')}
    native=a.rows(run.p.LOCAL/'supported_export/native_raw/nonexistent.csv') if False else a.rows(run.p.OUT/'supported_export/native_points.csv')
    rows=[];selection=[];source_activation=[]
    cases={j['case']:j for j in a.read(OUT/'contract.json')['jobs']}
    for case,job in cases.items():
        geo,mask=q.run.V.m.previous.prior.support(job)
        reference=next(r for r in native if r['case']==case and r['model']=='phumob' and int(r['index'])==40)
        point=dict(case=case,device=job['device'],vd=job['vd'],native_Id_A_per_um=float(reference['Id_A_per_um']))
        for stage in ('before','after'):
            chosen=[]
            for init in ('vela','native'):
                group=[r for r in sources[stage] if r['case']==case and r['arm']==init]
                selected=next((r for r in group if r['qualified']=='True'),group[-1]);chosen.append(selected)
                selection.append(dict(stage=stage,**selected))
                for component,value in a.read(Path(selected['dest'])/'independent_acceptance.json')['closure'].items():
                    source_activation.append(dict(stage=stage,case=case,initialization=init,carrier=component,**value))
            states=[d.ordered(Path(r['dest'])/'state.csv',geo.count) for r in chosen]
            delta=q.run.V.m.previous.prior.delta_states(*states,mask)
            ids=[float(r['current_A_per_um']) for r in chosen];id_relative=abs(ids[0]/ids[1]-1)
            dual=q.run.v.dual_qualified(*(r['qualified']=='True' for r in chosen),delta,id_relative)
            point.update({stage+'_'+k:v for k,v in delta.items()})
            point.update({stage+'_dual_qualified':dual,stage+'_dual_Id_relative':id_relative,stage+'_Id_A_per_um':ids[0],
                stage+'_versus_native_percent':100*(ids[0]/point['native_Id_A_per_um']-1)})
        point['Id_change_relative']=point['after_Id_A_per_um']/point['before_Id_A_per_um']-1
        rows.append(point)
    a.write_csv(OUT/'dc_comparison.csv',rows);a.write_csv(OUT/'selected_states.csv',selection);a.write_csv(OUT/'source_gate_activation.csv',source_activation)
    before_cross=a.read(run.OLD_CROSS/'active_transport_analysis/summary.json')
    after_cross=a.read(OUT/'cross_analysis/summary.json')
    summary=dict(points=len(rows),initial_attempts_before=sum(r['attempt']=='0' for r in sources['before']),initial_attempts_after=sum(r['attempt']=='0' for r in sources['after']),
        failed_attempts_before=sum(r['qualified']!='True' for r in sources['before']),failed_attempts_after=sum(r['qualified']!='True' for r in sources['after']),
        dual_points_before=sum(r['before_dual_qualified'] for r in rows),dual_points_after=sum(r['after_dual_qualified'] for r in rows),
        max_Id_change_relative=max(abs(r['Id_change_relative']) for r in rows),
        cross_checks_before=before_cross['nonzero_entry_checks'],cross_failed_before=before_cross['derivative_failures'],cross_failed_after=after_cross['derivative_failures'],
        source_relative_active=sum(r['qualified'] for r in source_activation),source_components=len(source_activation),
        production_changed=True,physics_or_acceptance_changed=False,native_element_box_phumob_restored=False,
        metadata_note='Reused cross analyzer writes production_changed=false to mean the diagnostic itself does not edit source; this experiment explicitly changes the production Jacobian, with baseline/candidate identities frozen separately.')
    a.write(OUT/'repair_summary.json',summary)
    d.matrix.freeze(OUT/'comparison_evidence.json',[Path(__file__).resolve(),OUT/'contract.json',OUT/'before/dc_evidence.json',OUT/'after/dc_evidence.json',OUT/'cross_analysis/evidence.json']+
        [OUT/p for p in ('dc_comparison.csv','selected_states.csv','source_gate_activation.csv','repair_summary.json')])
    print(summary,flush=True);print(rows,flush=True)


def post_jvp():
    a.verify(OUT/'comparison_evidence.json');jobs=[];files=[]
    for row in a.rows(OUT/'selected_states.csv'):
        if row['stage']!='after' or row['arm']!='vela':continue
        assert row['qualified']=='True'
        src=Path(row['dest']);cfg=a.read(src/'config.json');dest=LOCAL/'post_jvp'/row['case']
        _,mask=q.run.V.m.previous.prior.support(row)
        cfg.update(simulation_type='newton_jvp_probe',state_file=str(src/'state.csv'),output_csv=str(dest/'jvp.csv'))
        cfg.pop('output_state_file',None)
        cfg['directions']=[dict(name=f'{mode}_{step:.0e}',mode=mode,amplitude_V=step,node_ids=np.where(mask)[0][::3].tolist(),exclude_contacts=True)
            for mode in ('psi','phin','phip') for step in (1e-4,3e-5,1e-5,3e-6)]
        a.write(dest/'config.json',cfg);jobs.append(dict(case=row['case'],path=str(dest/'config.json')));files += [dest/'config.json',src/'state.csv']
    a.write(OUT/'post_jvp_contract.json',dict(jobs=jobs,gate=1e-4,source='Qualified final Vela-initialized state, all SRH settings restored.',weak_source_blocks='Retained but separately labeled; no weak SRH calibration inferred.'))
    d.matrix.freeze(OUT/'post_jvp_freeze.json',files+[OUT/'post_jvp_contract.json',OUT/'comparison_evidence.json',run.RUNNER,Path(__file__).resolve()])
    def one(job):
        path=Path(job['path']);status=q.run.V.execute(path,run.RUNNER,q.run.V.environment());assert status['exit_code']==0,status
        rows=[]
        for r in a.rows(path.parent/'jvp.csv'):
            for block in ('psi','phin','phip'):
                ana=float(r[f'analytic_{block}_norm']);fd=float(r[f'finite_difference_{block}_norm']);error=float(r[f'{block}_relative_error'])*max(1.,fd)
                relative=error/max(ana,fd,1e-300);weak=(r['mode'],block) in (('phin','phip'),('phip','phin'))
                rows.append(dict(case=job['case'],input=r['mode'],output=block,step_V=r['amplitude_V'],true_relative=relative,gated=not weak and float(r['amplitude_V'])<1e-4,qualified=relative<=1e-4))
        return rows
    with ThreadPoolExecutor(max_workers=2) as pool:rows=[r for part in pool.map(one,jobs) for r in part]
    a.write_csv(OUT/'post_jvp.csv',rows);gated=[r for r in rows if r['gated']]
    summary=dict(states=len(jobs),checks=len(rows),gated_checks=len(gated),gated_failures=sum(not r['qualified'] for r in gated),gated_max_relative=max(r['true_relative'] for r in gated))
    a.write(OUT/'post_jvp_summary.json',summary)
    d.matrix.freeze(OUT/'post_jvp_evidence.json',[OUT/'post_jvp_freeze.json',OUT/'post_jvp.csv',OUT/'post_jvp_summary.json']+[p for p in (LOCAL/'post_jvp').rglob('*') if p.is_file()])
    print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=('compare','post_jvp'))
    globals()[parser.parse_args().action]()
